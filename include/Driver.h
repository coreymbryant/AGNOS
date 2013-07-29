
#ifndef DRIVER_H
#define DRIVER_H

#include "agnosDefines.h"
#include "PseudoSpectralTensorProduct.h"
/* #include "PseudoSpectralMonteCarlo.h" */
/* #include "PseudoSpectralSparseGrid.h" */




namespace AGNOS
{

  /********************************************//**
   * \brief Base class for driving the construction of SurrogateModel objects
   ***********************************************/
  class Driver
  { 

    public:

      Driver( );
      Driver( Communicator& comm, Communicator& physicsComm, const GetPot& input );

      virtual ~Driver( );

      virtual void run( ) ;

      void printSolutionData( std::ostream& out ) ;
      void printSurrogateSettings( std::ostream& out ) ;
      void printParameterSettings( std::ostream& out ) ;
      void printDriverSettings( std::ostream& out  ) ;
      void printSettings( ) ;
      void printSolution( unsigned int iteration=1 ) ;

    protected:
      virtual void _buildPhysics( const GetPot& input );
      virtual void _buildSurrogate( const GetPot& input );

      Communicator& _comm;
      Communicator& _physicsComm;


      // DRIVER VARIABLES
      unsigned int _maxIter;

      // PARAMETERS VARIABLES
      unsigned int              _paramDim;
      std::vector<Parameter*>   _parameters;
      
      // ADAPTIVE SETTINGS
      // TODO adaptive settings
      bool _refinePhysical;
      bool _refineSurrogate;

      // SURROGATE VARIABLES
      int                       _surrogateType;
      std::vector<unsigned int> _order;
      std::vector<unsigned int> _errorOrder;
      SurrogateModel<T_S,T_P>*  _surrogate;
      SurrogateModel<T_S,T_P>*  _errorSurrogate;
      //
      // PHYSICS VARIABLES
      std::vector< PhysicsModel<T_S,T_P>* >                _physics;

      // OUTPUT VARIABLES
      std::string               _outputFilename; 
      std::vector<std::string>  _solutionsToPrint ;
      bool                      _outputIterations  ;
      bool                      _outputCoefficients  ;
      bool                      _outputErrorCoefficients  ;
      bool                      _outputWeights       ;
      bool                      _outputPoints        ;
      bool                      _outputIndexSet      ;
  };


/********************************************//**
 * \brief 
 * 
 ***********************************************/
  Driver::Driver( 
      Communicator& comm, 
      Communicator& physicsComm,
      const GetPot& input 
      ) : 
    _comm(comm), _physicsComm(physicsComm) 
  {
    if(DEBUG)
      std::cout << "rank: " << _comm.rank() <<  std::endl;
    
    // DRIVER SETTINGS
    _maxIter = input("driver/maxIter",1);
    
    // ADAPTIVE SETTINGS
    // TODO adaptive settings
    _refinePhysical = input("adaptive/refinePhysical",false);
    _refineSurrogate = input("adaptive/refineSurrogate",true);
    
    
    // PARAMETER SETTINGS
    _paramDim = input("parameters/dimension", 1);

    _parameters.resize(_paramDim);
    for (unsigned int i=0; i < _paramDim; i++)
      _parameters[i] = new Parameter( 
          input("parameters/types",0,i),
          input("parameters/mins",-1.0,i),
          input("parameters/maxs", 1.0,i)
          );
    
    // BUILD PHYSICS
    _buildPhysics( input );
    

    // SURROGATE MODEL SETTINGS
    _buildSurrogate( input );


    

    // OUTPUT DATA SETTINGS
    _outputFilename      = input("output/filename","cout");

    _solutionsToPrint.resize( input.vector_variable_size("output/solutions") );
    for (unsigned int i=0; i < _solutionsToPrint.size(); i++)
      _solutionsToPrint[i] = input("output/solutions", " ",i);

    _outputIterations    = input("output/iterations",false);

    _outputCoefficients  = input("output/coefficients",true);
    _outputErrorCoefficients  = input("output/errorCoefficients",true);
    _outputWeights       = input("output/weights",true);
    _outputPoints        = input("output/points",true);
    _outputIndexSet      = input("output/index_set",true);

  }

/********************************************//**
 * \brief 
 * 
 ***********************************************/
  Driver::~Driver( )
  {
    for(unsigned int i=0; i < _physics.size(); i++)
      delete _physics[i];
    _physics.clear();

    /* delete _comm; */
    /* delete _surrogate; */
    /* delete _errorSurrogate; */
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
    {
      if(DEBUG)
        std::cout << "_buildPhysics() rank: " << _comm.rank() << std::endl;
      _physics.push_back( 
          new AGNOS::PhysicsViscousBurgers<T_S,T_P>(
            _physicsComm, input )
          );
    }
    else if ( physicsName == "catenaryLibmesh" )
    {
      /* _physics.push_back( */
      /*     new AGNOS::PhysicsCatenaryLibmesh<T_S>( input ) */
      /*     ); */
    }
    else if ( physicsName == "catenary" )
    {
      /* _physics.push_back( */
      /*     new AGNOS::PhysicsCatenary<T_S>( input("physics/forcing",-10.0) ) */
      /*     ); */
    }
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

    int orderDim = input.vector_variable_size("surrogateModel/order");
    int errorOrderDim = input.vector_variable_size("surrogateModel/errorOrder");

    if (orderDim == 1)
      for (unsigned int i=0; i < _paramDim; i++)
        _order.push_back( input("surrogateModel/order", 0) ) ;
    else
      for (unsigned int i=0; i < _paramDim; i++)
        _order.push_back( input("surrogateModel/order", 0, i) ) ;

    // TODO make this a functino like we have in matlab code
    // i.e. incremental order increase
    if (errorOrderDim == 1)
      for (unsigned int i=0; i < _paramDim; i++)
        _errorOrder.push_back( input("surrogateModel/errorOrder", _order[i]+1) ) ;
    else
      for (unsigned int i=0; i < _paramDim; i++)
        _errorOrder.push_back( input("surrogateModel/errorOrder", _order[i]+1, i) ) ;

    std::string surrType  = 
      input("surrogateModel/type","PseudoSpectralTensorProduct");

    if (surrType == "PseudoSpectralTensorProduct")
      _surrogateType =  PSEUDO_SPECTRAL_TENSOR_PRODUCT ;
    else if (surrType == "PseudoSpectralSparseGrid")
      _surrogateType =  PSEUDO_SPECTRAL_SPARSE_GRID ;
    else if (surrType == "PseudoSpectralMonteCarlo")
      _surrogateType =  PSEUDO_SPECTRAL_MONTE_CARLO ;
    else if (surrType == "Collocation")
      _surrogateType =  COLLOCATION ;
    else
    {
      std::cerr << " ERROR: unrecognized SurrogateModelType " 
        << surrType << std::endl;
      exit(1);
    }
    


    // type specific setup
    switch( _surrogateType )
    {
      case(PSEUDO_SPECTRAL_TENSOR_PRODUCT):
        {

          _surrogate = new PseudoSpectralTensorProduct<T_S,T_P>(
              _comm,
              _physics[0], 
              _parameters, 
              _order  );

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

    // type specific setup
    switch( _surrogateType )
    {
      case(PSEUDO_SPECTRAL_TENSOR_PRODUCT):
        {

          _errorSurrogate = new PseudoSpectralTensorProduct<T_S,T_P>(
              _comm,
              _physics[0], 
              _parameters, 
              _errorOrder  );

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
    _surrogate->build();

    // build error surrogate
    /* _errorSurrogate->build(); */

    
    /* // print out first iteration if requested */
    /* if (this->_outputIterations && (_comm.rank() == 0) ) */
    /* { */
    /*   std::cout << "\n writing results to: " << this->_outputFilename */
    /*     << " (iter = " << 1 << " )" */
    /*    << std::endl; */
    /*   std::cout << std::endl; */
    /*   printSolution(1); */
    /* } */
    
    /*   // evaluate QoI */
    /*   T_S evalPoint(2); */
    /*   evalPoint(0) = 0.5; */
    /*   evalPoint(1) = 0.5; */
    /*   /1* T_S evalPoint(1); *1/ */
    /*   /1* evalPoint(0) = 1.5; *1/ */
    /*   T_P solutionVec = _surrogate->evaluate("primal", evalPoint ); */
    /*   T_P qoiValue = _physics[0]->evaluateQoi( evalPoint, solutionVec ) ; */

    /*   T_P l2normofphyerror = _surrogate->l2Norm("error"); */
    /*   T_P l2normoftotalerror = _errorSurrogate->l2Norm("error"); */
    /*   double normDiff = _surrogate->l2NormDifference( *_errorSurrogate, "error"); */

    /* if (_comm.rank() == 0) */
    /* { */

    /*   std::cout << "phyErrorNorm = " << l2normofphyerror << std::endl; */
    /*   std::cout << "totalErrorNorm = " << l2normoftotalerror << std::endl; */
    /*   std::cout << "| phyErrorNorm-totalErrorNorm | = " << */
    /*     std::abs(l2normofphyerror(0)-l2normoftotalerror(0)) << std::endl; */
    /*   std::cout << "normDiff = " << normDiff << std::endl; */

    /*   std::cout << "\n Qoi = " << qoiValue(0) << std::endl; */
    /* } */

    /* // refine approximation */
    /* for (unsigned int iter=2; iter <= this->_maxIter; iter++) */
    /* { */
    /*   if (_comm.rank() == 0) */ 
    /*     std::cout << "\n-------------  ITER " << iter << "  -------------\n " ; */
      

    /*   // if physical error dominates and we are allowed to refine physical */
    /*   // solution */
    /*   if ( */ 
    /*       ( _refinePhysical && (l2normofphyerror(0) >= normDiff) ) */ 
    /*       || */
    /*       ( _refinePhysical && !m_refineSurrogate ) */
    /*       ) */
    /*   { */
    /*     if (_comm.rank() == 0) */ 
    /*       std::cout << "    refining physical solution " << std::endl; */
        
    /*     // if using uniform refinement we won't have error indicators in the */
    /*     // surrogate model (by design) */
    /*     if ( _physicsFunctions.count("indicators") == 0 ) */
    /*       _physics[0]->refine( ); */
    /*     else */
    /*     { */
    /*       // retrieve first coefficient (i.e. the mean) of error inidcators */
    /*       T_P errorIndicators = (_surrogate->mean())["indicators"]; */
          
    /*       _physics[0]->refine( errorIndicators ); */
    /*     } */
    /*     if (_comm.rank() == 0) */ 
    /*       std::cout << "-------------------------------------\n " ; */

    /*     _surrogate->build(); */


    /*     _errorSurrogate->build(); */
    /*   } */


    /*   // if surrogate error dominates and we are allowed to refine surrogate */
    /*   // model */
    /*   if ( */ 
    /*       ( _refineSurrogate && (l2normofphyerror(0) < normDiff) ) */ 
    /*       || */
    /*       ( _refineSurrogate && !m_refinePhysical ) */ 
    /*       ) */
    /*   { */
    /*     if (_comm.rank() == 0) */ 
    /*     { */
    /*       std::cout << "    refining surrogate model " << std::endl; */
    /*       std::cout << "-------------------------------------\n " ; */
    /*     } */
    /*     _surrogate->refine( ); */
    /*     _errorSurrogate->refine( ); */
    /*   } */





    /*   if (this->_outputIterations && (_comm.rank() == 0) ) */
    /*   { */
    /*     std::cout << "\n writing results to: " << this->_outputFilename */
    /*       << " (iter = " << iter << " )" */
    /*       << std::endl; */
    /*     std::cout << std::endl; */
    /*     printSolution(iter); */
    /*   } */

    /*   solutionVec = _surrogate->evaluate("primal", evalPoint ); */
    /*   qoiValue = _physics[0]->evaluateQoi( evalPoint, solutionVec ) ; */

    /*   l2normofphyerror = _surrogate->l2Norm("error"); */
    /*   l2normoftotalerror = _errorSurrogate->l2Norm("error"); */
    /*   normDiff = _surrogate->l2NormDifference( *_errorSurrogate, "error"); */

    /*   if (_comm.rank() == 0) */
    /*   { */

    /*     std::cout << "phyErrorNorm = " << l2normofphyerror << std::endl; */
    /*     std::cout << "totalErrorNorm = " << l2normoftotalerror << std::endl; */
    /*     /1* std::cout << "| phyErrorNorm-totalErrorNorm | = " << *1/ */
    /*     /1* std::abs(l2normofphyerror(0)-l2normoftotalerror(0)) << std::endl; *1/ */
    /*     std::cout << "normDiff = " << normDiff << std::endl; */

    /*     std::cout << "\n Qoi = " << qoiValue(0) << std::endl; */
    /*   } */

      
      /* // evaluate QoI */
      /* T_S evalPoint(1); */
      /* evalPoint(0) = 1.5; */
      /* T_P solutionVec = _surrogate->evaluate("primal", evalPoint ); */
      /* T_P qoiValue = _physics->evaluateQoi( evalPoint, solutionVec ) ; */
      /* if (_comm.rank() == 0) */
      /* { */
      /*   std::cout << "\n Qoi = " << qoiValue(0) << std::endl; */
      /* } */
    /* } */
    
    /* // output whatever user asks for */
    /* if (_comm.rank() == 0) */
    /* { */
    /*   std::cout << "\n writing final results to: " << this->_outputFilename */
    /*     << std::endl; */
    /*   std::cout << std::endl; */
    /* } */
    /*   printSolution(_maxIter); */






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
    out << "#    maxIter = " << _maxIter << std::endl;

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
    out << "#     dimension = " << _paramDim << std::endl;
    out << "#     order = " ;
    for(unsigned int i=0; i < _paramDim; i++)
      out << _order[i] << " " ;
    out << std::endl;
    out << "#     errorOrder = " ;
    for(unsigned int i=0; i < _paramDim; i++)
      out << _errorOrder[i] << " " ;
    out << std::endl;

    out << "#     mins = " ;
    for (unsigned int i=0; i < _paramDim; i++)
      out << _parameters[i]->min() << " " ;
    out << std::endl;

    out << "#     maxs = " ;
    for (unsigned int i=0; i < _paramDim; i++)
      out << _parameters[i]->max() << " " ;
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
      if (_outputCoefficients)
        _surrogate->printCoefficients( _solutionsToPrint, out );
      if (_outputErrorCoefficients)
      {
        std::vector<std::string> errorSols;
        errorSols.push_back("error");
        _errorSurrogate->printCoefficients( errorSols, out );
      }
      if (_outputWeights)
        _surrogate->printIntegrationWeights( out );
      if (_outputPoints)
        _surrogate->printIntegrationPoints( out );
      if (_outputIndexSet)
        _surrogate->printIndexSet( out );

    out << std::endl;
    return;
  }



/********************************************//**
 * \brief 
 ***********************************************/
  void Driver::printSettings( ) 
  {

    // set output steam
    if (_outputFilename == "cout" )
    {
      printDriverSettings( std::cout  ) ;
      printParameterSettings( std::cout ) ;
      printSurrogateSettings( std::cout ) ;
    }
    else
    {
      std::ofstream out;
      out.open( _outputFilename.c_str() );

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
    if (_outputFilename == "cout" )
    {
      std::cout << std::endl;
      std::cout << "#====================================================\n" 
         << "#      Solution data for ITER " << iteration << std::endl;
      printSolutionData( std::cout ) ;
    }
    else
    {
      std::ofstream out;
      out.open( _outputFilename.c_str(), std::ofstream::out | std::ofstream::app );

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


