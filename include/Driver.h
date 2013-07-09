
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
      Driver( const Communicator& comm, const GetPot& input );

      virtual ~Driver( );

      virtual void run( ) = 0 ;

      void printSolutionData( std::ostream& out ) ;
      void printSurrogateSettings( std::ostream& out ) ;
      void printParameterSettings( std::ostream& out ) ;
      void printDriverSettings( std::ostream& out  ) ;
      void printSettings( ) ;
      void printSolution( unsigned int iteration=1 ) ;

    protected:
      const Communicator* m_comm;


      // DRIVER VARIABLES
      unsigned int m_maxIter;

      // PARAMETERS VARIABLES
      unsigned int              m_paramDim;
      std::vector<Parameter*>   m_parameters;
      
      // ADAPTIVE SETTINGS
      // TODO adaptive settings

      // SURROGATE VARIABLES
      int                       m_surrogateType;
      std::vector<unsigned int> m_order;
      std::vector<unsigned int> m_errorOrder;
      SurrogateModel<T_S,T_P>*  m_surrogate;
      SurrogateModel<T_S,T_P>*  m_errorSurrogate;

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
    
    
    // PARAMETER SETTINGS
    m_paramDim = input("parameters/dimension", 1);

    m_parameters.resize(m_paramDim);
    for (unsigned int i=0; i < m_paramDim; i++)
      m_parameters[i] = new Parameter( 
          input("parameters/types",0,i),
          input("parameters/mins",-1.0,i),
          input("parameters/maxs", 1.0,i)
          );

    // SURROGATE MODEL SETTINGS
    for (unsigned int i=0; i < m_paramDim; i++)
      m_order.push_back( input("surrogateModel/order", 0, i) ) ;
    // TODO make this a functino like we have in matlab code
    for (unsigned int i=0; i < m_paramDim; i++)
      m_errorOrder.push_back( input("surrogateModel/errorOrder", m_order[i], i) ) ;

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
      std::cerr << "unrecognized SurrogateModelType " 
        << surrType << std::endl;
      exit(1);
    }

    

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

    /* delete m_comm; */
    /* delete m_surrogate; */
    /* delete m_errorSurrogate; */
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


