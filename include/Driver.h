
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

      // TODO default inputs settings
      Driver( );
      Driver( const Communicator& comm, GetPot& input );

      virtual ~Driver( );

      virtual void run( ) = 0 ;

      // TODO clean up output headings
      void printSolutionData( std::ostream& out ) ;
      void printSurrogateSettings( std::ostream& out ) ;
      void printParameterSettings( std::ostream& out ) ;
      void printDriverSettings( std::ostream& out  ) ;
      void printSettings( ) ;
      void printSolution( unsigned int iteration=1 ) ;

    protected:
      const Communicator* m_comm;
      GetPot m_input;


      // DRIVER VARIABLES
      unsigned int m_maxIter;
      // TODO adaptivity settings

      // PARAMETERS VARIABLES
      unsigned int              m_paramDim;
      std::vector<Parameter*>   m_parameters;


      // SURROGATE VARIABLES
      int                       m_surrogateType;
      std::vector<unsigned int> m_order;
      SurrogateModel<T_S,T_P>*  m_surrogate;

      // OUTPUT VARIABLES
      std::string               m_outputFilename; 
      std::vector<std::string>  m_solutionsToPrint ;
      // TODO some sort of update output used for each iteration
      bool                      m_outputIterations  ;
      bool                      m_outputCoefficients  ;
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
      GetPot& input ) 
    : m_comm(&comm), m_input(input)
  {
    
    // DRIVER SETTINGS
    m_input.set_prefix("driver/");
    m_maxIter = m_input("maxIter",1);
    
    // TODO adaptive settings
    
    
    // PARAMETER SETTINGS
    m_input.set_prefix("parameters/");
    m_paramDim = m_input("dimension", 1);

    m_parameters.resize(m_paramDim);
    for (unsigned int i=0; i < m_paramDim; i++)
      m_parameters[i] = new Parameter( 
          m_input("types",0,i),
          m_input("mins",-1.0,i),
          m_input("maxs", 1.0,i)
          );

    // SURROGATE MODEL SETTINGS
    m_input.set_prefix("surrogateModel/");
    for (unsigned int i=0; i < m_paramDim; i++)
      m_order.push_back( m_input("order", 0, i) ) ;

    std::string surrType  = 
      m_input("type","PseudoSpectralTensorProduct");

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
      std::cerr << "unrecognized SurrogateModelType " 
        << surrType << std::endl;
      exit(1);
    }

    

    // OUTPUT DATA SETTINGS
    m_input.set_prefix("output/");
    m_outputFilename      = m_input("filename","cout");

    m_solutionsToPrint.resize( m_input.vector_variable_size("solutions") );
    for (unsigned int i=0; i < m_solutionsToPrint.size(); i++)
      m_solutionsToPrint[i] = m_input("solutions", " ",i);

    m_outputIterations    = m_input("iterations",false);

    m_outputCoefficients  = m_input("coefficients",true);
    m_outputWeights       = m_input("weights",true);
    m_outputPoints        = m_input("points",true);
    m_outputIndexSet      = m_input("index_set",true);

  }

/********************************************//**
 * \brief 
 * 
 ***********************************************/
  Driver::~Driver( )
  {
    for (unsigned int i=0; i < m_paramDim; i++)
      delete m_parameters[i];
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


