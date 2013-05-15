
#ifndef DRIVER_H
#define DRIVER_H


#include <iostream>
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
    // TODO much later

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




}


#endif // DRIVER_H


