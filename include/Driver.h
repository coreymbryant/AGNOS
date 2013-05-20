
#ifndef DRIVER_H
#define DRIVER_H


#include <iostream>
#include <fstream>
#include <stdio.h>
#include <assert.h>
#include <cstring>
#include <GetPot>


// TODO find a better way to include all necessary submodels?
#include "Parameter.h"
#include "PseudoSpectralTensorProduct.h"
/* #include "PseudoSpectralSparseGrid.h" */
/* #include "PseudoSpectralMonteCarlo.h" */
/* #include "SurrogateCollocation.h" */

#include "libmesh/dense_vector.h"
// move these into same file as above includes ??
typedef libMesh::DenseVector<double> T_P ;
typedef libMesh::DenseVector<double> T_S ;


namespace AGNOS
{

  /********************************************//**
   * \brief Base class for driving the construction of SurrogateModel objects
   ***********************************************/
  class Driver
  { 

    public:

      // TODO default inputs settings
      Driver(  );
      Driver( GetPot& input );

      virtual ~Driver( );

      virtual void run( ) = 0 ;

      void printSolution( std::ostream& out ) ;
      void printSurrogateSettings( std::ostream& out ) ;
      void printParameterSettings( std::ostream& out ) ;
      void printDriverSettings( std::ostream& out  ) ;
      void printOutput( ) ;

    protected:
      GetPot m_input;
      void initializeSurrogateModel( int surrogateType, GetPot& input );


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
  Driver::Driver( GetPot& input ) : m_input(input)
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
    //TODO header
    out << "#maxIter = " << m_maxIter << std::endl;

    out << std::endl;
    return;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  void Driver::printParameterSettings( std::ostream& out ) 
  {
    //TODO header
    out << "#dimension = " << m_paramDim << std::endl;
    out << "#order = " ;
    for(unsigned int i=0; i < m_paramDim; i++)
      out << m_order[i] << " " ;
    out << std::endl;

    out << "#mins = " ;
    for (unsigned int i=0; i < m_paramDim; i++)
      out << "m_parameters[i]->min() << " ;
    out << std::endl;

    out << "#maxs = " ;
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
  void Driver::printSolution( std::ostream& out ) 
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
  void Driver::printOutput( ) 
  {

    // set output steam
    if (m_outputFilename == "cout" )
    {
      printDriverSettings( std::cout  ) ;
      printParameterSettings( std::cout ) ;
      printSurrogateSettings( std::cout ) ;
      printSolution( std::cout ) ;
    }
    else
    {
      std::ofstream out;
      out.open( m_outputFilename.c_str() );

      printDriverSettings( out  ) ;
      printParameterSettings( out ) ;
      printSurrogateSettings( out ) ;
      printSolution( out ) ;

      out.close( );
    }

    return;
  }

    
}





#endif // DRIVER_H


