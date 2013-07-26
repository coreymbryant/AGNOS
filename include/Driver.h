
#ifndef DRIVER_H
#define DRIVER_H

#include "agnosDefines.h"


namespace AGNOS
{

  /********************************************//**
   * \brief Base class for driving the construction of SurrogateModel objects
   ***********************************************/
  class Driver
  { 

    public:

      Driver( );
      Driver( const Communicator& comm, const GetPot& input );

      virtual ~Driver( );

      void run( ) ;

      void printSolutionData( std::ostream& out ) ;
      void printSurrogateSettings( std::ostream& out ) ;
      void printParameterSettings( std::ostream& out ) ;
      void printDriverSettings( std::ostream& out  ) ;
      void printSettings( ) ;
      void printSolution( unsigned int iteration=1 ) ;


    protected:
      virtual void _buildPhysics( const GetPot& input );
      virtual void _buildSurrogate( const GetPot& input );

      const Communicator* m_comm;

      // DRIVER VARIABLES
      unsigned int m_maxIter;

      // PARAMETERS VARIABLES
      unsigned int              m_paramDim;
      std::vector<Parameter*>   m_parameters;
      
      // ADAPTIVE SETTINGS
      // TODO adaptive settings
      bool m_refinePhysical;
      bool m_refineSurrogate;

      // SURROGATE VARIABLES
      int                       m_surrogateType;
      std::vector<unsigned int> m_order;
      std::vector<unsigned int> m_errorOrder;
      SurrogateModel<T_S,T_P>*  m_surrogate;
      SurrogateModel<T_S,T_P>*  m_errorSurrogate;
      
      // PHYSICS VARIABLES
      PhysicsModel<T_S,T_P>*                              m_physics;
      std::map< std::string, PhysicsFunction<T_S,T_P>* >  m_physicsFunctions ;
      std::map< std::string, PhysicsFunction<T_S,T_P>* >  m_errorFunctions ;

      // OUTPUT VARIABLES
      std::string               m_outputFilename; 
      std::vector<std::string>  m_solutionsToPrint ;
      bool                      m_outputIterations  ;
      bool                      m_outputCoefficients  ;
      bool                      m_outputErrorCoefficients  ;
      bool                      m_outputWeights       ;
      bool                      m_outputPoints        ;
      bool                      m_outputIndexSet      ;
  };


/********************************************//**
 * \brief 
 * 
 ***********************************************/
  Driver::Driver( 
      const Communicator& comm, 
      const GetPot& input ) 
    : m_comm(&comm)
  {
    
    // DRIVER SETTINGS
    m_maxIter = input("driver/maxIter",1);

    
    // ADAPTIVE SETTINGS
    // TODO adaptive settings
    m_refinePhysical = input("adaptive/refinePhysical",false);
    m_refineSurrogate = input("adaptive/refineSurrogate",true);
    
    
    // PARAMETER SETTINGS
    m_paramDim = input("parameters/dimension", 1);

    m_parameters.resize(m_paramDim);
    for (unsigned int i=0; i < m_paramDim; i++)
      m_parameters[i] = new Parameter( 
          input("parameters/types",0,i),
          input("parameters/mins",-1.0,i),
          input("parameters/maxs", 1.0,i)
          );
    
    // BUILD PHYSICS
    _buildPhysics( input );
    

    // SURROGATE MODEL SETTINGS
    _buildSurrogate( input );


    // OUTPUT DATA SETTINGS
    m_outputFilename      = input("output/filename","cout");

    m_solutionsToPrint.resize( input.vector_variable_size("output/solutions") );
    for (unsigned int i=0; i < m_solutionsToPrint.size(); i++)
      m_solutionsToPrint[i] = input("output/solutions", " ",i);

    m_outputIterations    = input("output/iterations",false);

    m_outputCoefficients  = input("output/coefficients",true);
    m_outputErrorCoefficients  = input("output/errorCoefficients",true);
    m_outputWeights       = input("output/weights",true);
    m_outputPoints        = input("output/points",true);
    m_outputIndexSet      = input("output/index_set",true);

  }

/********************************************//**
 * \brief 
 * 
 ***********************************************/
  Driver::~Driver( )
  {

    delete m_physics;
    delete m_comm;
    /* delete m_surrogate; */
    /* delete m_errorSurrogate; */
  }


/********************************************//**
 * \brief build physics model from default classes.
 * Can be overidden to include new physicsModel classes.
 * 
 ***********************************************/
  void Driver::_buildPhysics( const GetPot& input)
  {
    //TODO make more general 
    // deal with parallel issue
    std::string physicsName = input("physics/type","");
    if ( physicsName == "viscousBurgers" )
      m_physics = new AGNOS::PhysicsViscousBurgers<T_S,T_P>( *m_comm, input );
    else if ( physicsName == "catenaryLibmesh" )
      m_physics = new AGNOS::PhysicsCatenaryLibmesh<T_S,T_P>( *m_comm, input );
    else if ( physicsName == "catenary" )
      m_physics = new AGNOS::PhysicsCatenary<T_S,T_P>(
          input("physics/forcing",-10.0) );
    else
    {
      std::cerr << " ERROR: unrecognized physics type: " << physicsName
        << "\n"
        << "please choose an appropriate physics type "
        << "(or implememnt it yourself) \n"
        << std::endl;
      exit(1);
    }

    return ;
  }

/********************************************//**
 * \brief build surrogateModel
 * Can be overidden to include new surrogateModel classes.
 * 
 ***********************************************/
  void Driver::_buildSurrogate( const GetPot& input)
  {

    for (unsigned int i=0; i < m_paramDim; i++)
      m_order.push_back( input("surrogateModel/order", 0, i) ) ;
    // TODO make this a functino like we have in matlab code
    // i.e. incremental order increase
    for (unsigned int i=0; i < m_paramDim; i++)
      m_errorOrder.push_back( input("surrogateModel/errorOrder", m_order[i]+1, i) ) ;

    std::string surrType  = 
      input("surrogateModel/type","PseudoSpectralTensorProduct");

    if (surrType == "PseudoSpectralTensorProduct")
      m_surrogateType =  PSEUDO_SPECTRAL_TENSOR_PRODUCT ;
    else if (surrType == "PseudoSpectralSparseGrid")
      m_surrogateType =  PSEUDO_SPECTRAL_SPARSE_GRID ;
    else if (surrType == "PseudoSpectralMonteCarlo")
      m_surrogateType =  PSEUDO_SPECTRAL_MONTE_CARLO ;
    else if (surrType == "Collocation")
      m_surrogateType =  COLLOCATION ;
    else
    {
      std::cerr << " ERROR: unrecognized SurrogateModelType " 
        << surrType << std::endl;
      exit(1);
    }
    
    // primal solution
    m_physicsFunctions.insert( 
        std::pair< std::string, PhysicsFunction<T_S,T_P>* >(
          "primal", 
          new PhysicsFunctionPrimal<T_S,T_P>( m_physics ) ) 
        );

    // adjoint solution
    m_physicsFunctions.insert( 
        std::pair< std::string, PhysicsFunction<T_S,T_P>* >(
          "adjoint", 
          new PhysicsFunctionAdjoint<T_S,T_P>( m_physics ) ) 
        ); 

    // qoi evaluation
    m_physicsFunctions.insert( 
        std::pair< std::string, PhysicsFunction<T_S,T_P>* >(
          "qoi", 
          new PhysicsFunctionQoi<T_S,T_P>( m_physics ) ) 
        ); 
    
    // error estimate
    m_physicsFunctions.insert( 
        std::pair< std::string, PhysicsFunction<T_S,T_P>* >(
          "error", 
          new PhysicsFunctionError<T_S,T_P>( m_physics ) ) 
        ); 

    // error indicators if needed
    if ( ! m_physics->useUniformRefinement() )
    {
      m_physicsFunctions.insert( 
          std::pair< std::string, PhysicsFunction<T_S,T_P>* >(
            "indicators", 
            new PhysicsFunctionIndicators<T_S,T_P>( m_physics ) ) 
          ); 
    }


    // type specific setup
    switch( m_surrogateType )
    {
      case(PSEUDO_SPECTRAL_TENSOR_PRODUCT):
        {

          m_surrogate = new PseudoSpectralTensorProduct<T_S,T_P>(
              m_comm,
              m_physicsFunctions, 
              m_parameters, 
              m_order  );

          break;
        }
      case(PSEUDO_SPECTRAL_SPARSE_GRID):
        {
          std::cerr 
            << " this SurrogateModelType is not yet implemented\n" ;
          exit(1);
          break;
        }
      case(PSEUDO_SPECTRAL_MONTE_CARLO):
        {
          std::cerr 
            << " this SurrogateModelType is not yet implemented\n" ;
          exit(1);
          break;
        }
      case(COLLOCATION):
        {
          std::cerr 
            << " this SurrogateModelType is not yet implemented\n" ;
          exit(1);
          break;
        }
    }
    

    // ----- ERROR SURROGATE
    //
    m_errorFunctions.insert( 
        std::pair< std::string, PhysicsFunction<T_S,T_P>* >(
          "error", 
          new PhysicsFunctionTotalError<T_S,T_P>( 
            m_physics, m_surrogate ) 
          ) 
        ); 

    // type specific setup
    switch( m_surrogateType )
    {
      case(PSEUDO_SPECTRAL_TENSOR_PRODUCT):
        {

          m_errorSurrogate = new PseudoSpectralTensorProduct<T_S,T_P>(
              m_comm,
              m_errorFunctions, 
              m_parameters, 
              m_errorOrder  );

          break;
        }
      case(PSEUDO_SPECTRAL_SPARSE_GRID):
        {
          std::cerr 
            << " this SurrogateModelType is not yet implemented\n" ;
          exit(1);
          break;
        }
      case(PSEUDO_SPECTRAL_MONTE_CARLO):
        {
          std::cerr 
            << " this SurrogateModelType is not yet implemented\n" ;
          exit(1);
          break;
        }
      case(COLLOCATION):
        {
          std::cerr 
            << " this SurrogateModelType is not yet implemented\n" ;
          exit(1);
          break;
        }
    }

    return;
  }

/********************************************//**
 * \brief an initial driver run routine for testing
 * 
 ***********************************************/
  void Driver::run( )
  {

    // print out settings
    printSettings();
    
    // build initial approximation
    m_surrogate->build();
    
    // build error surrogate
    m_errorSurrogate->build();

    
    // print out first iteration if requested
    if (this->m_outputIterations && (this->m_comm->rank() == 0) )
    {
      std::cout << "\n writing results to: " << this->m_outputFilename
        << " (iter = " << 1 << " )"
       << std::endl;
      std::cout << std::endl;
      printSolution(1);
    }
    
      // evaluate QoI
      T_S evalPoint(2);
      evalPoint(0) = 0.5;
      evalPoint(1) = 0.5;
      /* T_S evalPoint(1); */
      /* evalPoint(0) = 1.5; */
      T_P solutionVec = m_surrogate->evaluate("primal", evalPoint );
      T_P qoiValue = m_physics->evaluateQoi( evalPoint, solutionVec ) ;

      T_P l2normofphyerror = m_surrogate->l2Norm("error");
      T_P l2normoftotalerror = m_errorSurrogate->l2Norm("error");
      double normDiff = m_surrogate->l2NormDifference( *m_errorSurrogate, "error");

    if (this->m_comm->rank() == 0)
    {

      std::cout << "phyErrorNorm = " << l2normofphyerror << std::endl;
      std::cout << "totalErrorNorm = " << l2normoftotalerror << std::endl;
      std::cout << "| phyErrorNorm-totalErrorNorm | = " <<
        std::abs(l2normofphyerror(0)-l2normoftotalerror(0)) << std::endl;
      std::cout << "normDiff = " << normDiff << std::endl;

      std::cout << "\n Qoi = " << qoiValue(0) << std::endl;
    }

    // refine approximation
    for (unsigned int iter=2; iter <= this->m_maxIter; iter++)
    {
      if (this->m_comm->rank() == 0) 
        std::cout << "\n-------------  ITER " << iter << "  -------------\n " ;
      

      // if physical error dominates and we are allowed to refine physical
      // solution
      if ( 
          ( m_refinePhysical && (l2normofphyerror(0) >= normDiff) ) 
          ||
          ( m_refinePhysical && !m_refineSurrogate )
          )
      {
        if (this->m_comm->rank() == 0) 
          std::cout << "    refining physical solution " << std::endl;
        
        // if using uniform refinement we won't have error indicators in the
        // surrogate model (by design)
        if ( m_physicsFunctions.count("indicators") == 0 )
          m_physics->refine( );
        else
        {
          // retrieve first coefficient (i.e. the mean) of error inidcators
          T_P errorIndicators = (m_surrogate->mean())["indicators"];
          
          m_physics->refine( errorIndicators );
        }
        if (this->m_comm->rank() == 0) 
          std::cout << "-------------------------------------\n " ;

        m_surrogate->build();
        m_errorSurrogate->build();
      }


      // if surrogate error dominates and we are allowed to refine surrogate
      // model
      if ( 
          ( m_refineSurrogate && (l2normofphyerror(0) < normDiff) ) 
          ||
          ( m_refineSurrogate && !m_refinePhysical ) 
          )
      {
        if (this->m_comm->rank() == 0) 
        {
          std::cout << "    refining surrogate model " << std::endl;
          std::cout << "-------------------------------------\n " ;
        }
        m_surrogate->refine( );
        m_errorSurrogate->refine( );
      }





      if (this->m_outputIterations && (this->m_comm->rank() == 0) )
      {
        std::cout << "\n writing results to: " << this->m_outputFilename
          << " (iter = " << iter << " )"
          << std::endl;
        std::cout << std::endl;
        printSolution(iter);
      }

      this->m_comm->barrier();

      solutionVec = m_surrogate->evaluate("primal", evalPoint );
      qoiValue = m_physics->evaluateQoi( evalPoint, solutionVec ) ;

      l2normofphyerror = m_surrogate->l2Norm("error");
      l2normoftotalerror = m_errorSurrogate->l2Norm("error");
      normDiff = m_surrogate->l2NormDifference( *m_errorSurrogate, "error");

      if (this->m_comm->rank() == 0)
      {

        std::cout << "phyErrorNorm = " << l2normofphyerror << std::endl;
        std::cout << "totalErrorNorm = " << l2normoftotalerror << std::endl;
        /* std::cout << "| phyErrorNorm-totalErrorNorm | = " << */
        /* std::abs(l2normofphyerror(0)-l2normoftotalerror(0)) << std::endl; */
        std::cout << "normDiff = " << normDiff << std::endl;

        std::cout << "\n Qoi = " << qoiValue(0) << std::endl;
      }

      
      /* // evaluate QoI */
      /* T_S evalPoint(1); */
      /* evalPoint(0) = 1.5; */
      /* T_P solutionVec = m_surrogate->evaluate("primal", evalPoint ); */
      /* T_P qoiValue = m_physics->evaluateQoi( evalPoint, solutionVec ) ; */
      /* if (this->m_comm->rank() == 0) */
      /* { */
      /*   std::cout << "\n Qoi = " << qoiValue(0) << std::endl; */
      /* } */
    }
    
    // output whatever user asks for
    if (this->m_comm->rank() == 0)
    {
      std::cout << "\n writing final results to: " << this->m_outputFilename
        << std::endl;
      std::cout << std::endl;
    }
      printSolution(m_maxIter);






    return;
  }
    

/********************************************//**
 * \brief 
 ***********************************************/
  void Driver::printDriverSettings( std::ostream& out  ) 
  {
    out << std::endl;
    out << "#====================================================" <<
      std::endl;
    out << "# Driver settings: " << std::endl;
    out << "#    maxIter = " << m_maxIter << std::endl;

    out << std::endl;
    return;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  void Driver::printParameterSettings( std::ostream& out ) 
  {
    out << std::endl;
    out << "#====================================================" <<
      std::endl;
    out << "# Parameter settings: " << std::endl;
    out << "#     dimension = " << m_paramDim << std::endl;
    out << "#     order = " ;
    for(unsigned int i=0; i < m_paramDim; i++)
      out << m_order[i] << " " ;
    out << std::endl;
    out << "#     errorOrder = " ;
    for(unsigned int i=0; i < m_paramDim; i++)
      out << m_errorOrder[i] << " " ;
    out << std::endl;

    out << "#     mins = " ;
    for (unsigned int i=0; i < m_paramDim; i++)
      out << m_parameters[i]->min() << " " ;
    out << std::endl;

    out << "#     maxs = " ;
    for (unsigned int i=0; i < m_paramDim; i++)
      out << m_parameters[i]->max() << " " ;
    out << std::endl;

    out << std::endl;
    return;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  void Driver::printSurrogateSettings( std::ostream& out ) 
  {
    out << std::endl;
    return;
  }


/********************************************//**
 * \brief 
 ***********************************************/
  void Driver::printSolutionData( std::ostream& out ) 
  {
      if (m_outputCoefficients)
        m_surrogate->printCoefficients( m_solutionsToPrint, out );
      if (m_outputErrorCoefficients)
      {
        std::vector<std::string> errorSols;
        errorSols.push_back("error");
        m_errorSurrogate->printCoefficients( errorSols, out );
      }
      if (m_outputWeights)
        m_surrogate->printIntegrationWeights( out );
      if (m_outputPoints)
        m_surrogate->printIntegrationPoints( out );
      if (m_outputIndexSet)
        m_surrogate->printIndexSet( out );

    out << std::endl;
    return;
  }



/********************************************//**
 * \brief 
 ***********************************************/
  void Driver::printSettings( ) 
  {

    // set output steam
    if (m_outputFilename == "cout" )
    {
      printDriverSettings( std::cout  ) ;
      printParameterSettings( std::cout ) ;
      printSurrogateSettings( std::cout ) ;
    }
    else
    {
      std::ofstream out;
      out.open( m_outputFilename.c_str() );

      printDriverSettings( out  ) ;
      printParameterSettings( out ) ;
      printSurrogateSettings( out ) ;

      out.close( );
    }

    return;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  void Driver::printSolution( unsigned int iteration ) 
  {
    // set output steam
    if (m_outputFilename == "cout" )
    {
      std::cout << std::endl;
      std::cout << "#====================================================\n" 
         << "#      Solution data for ITER " << iteration << std::endl;
      printSolutionData( std::cout ) ;
    }
    else
    {
      std::ofstream out;
      out.open( m_outputFilename.c_str(), std::ofstream::out | std::ofstream::app );

      out << std::endl;
      out << "#====================================================\n" 
         << "#       Solution data for ITER " << iteration << std::endl;
      printSolutionData( out ) ;

      out.close( );
    }
    return;
  }


}
#endif // DRIVER_H


