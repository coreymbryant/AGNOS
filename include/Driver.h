
#ifndef DRIVER_H
#define DRIVER_H


#include <iostream>
#include <stdio.h>
#include <assert.h>
#include <GetPot>

// TODO find a better way to include all necessary submodels?
#include "Parameter.h"
#include "PseudoSpectralTensorProduct.h"
#include "PhysicsCatenary.h"
#include "PhysicsFunction.h"

#include "libmesh/dense_vector.h"


namespace AGNOS
{
  // move these into same file as above includes ??
  typedef libMesh::DenseVector<double> T_P ;
  typedef libMesh::DenseVector<double> T_S ;

  /********************************************//**
   * \brief Base class for driving the construction of SurrogateModel objects
   *
   * This Base class is designed to be a non-adaptive (i.e. single construction)
   * method for surrogate construction.  Several different adaptive strategies
   * can be implemented as derived classes of this routine. 
   *
   * TODO this should be defined based on user inputs. Or at the very least
   * require a user provided derived class to be defined for the particular
   * application.
   * 
   *
   * 
   ***********************************************/
  class Driver
  { 

    public:

      Driver( );
      Driver( GetPot& input );
      virtual ~Driver( );

      void run( ) ;

    private:
      GetPot m_input;
      void initializeSurrogateModel( int surrogateType, GetPot& input );
      void initializePhysicsModel( int physicsType, GetPot& input );


      // DRIVER VARIABLES
      unsigned int m_maxIter;

      // PARAMETERS VARIABLES
      unsigned int              m_paramDim;
      std::vector<Parameter*>   m_parameters;


      // PHYSICS VARIABLES
      enum PhysicsModelType{
        PHYSICS_CATENARY=0,
        PHYSICS_VISCOUS_BURGERS,
        PHYSICS_GRINS };
      int                       m_physicsType;
      PhysicsModel<T_S,T_P>*    m_physicsModel;
      PhysicsFunction<T_S,T_P>* m_primalPhysics;
      PhysicsFunction<T_S,T_P>* m_adjointPhysics;
      PhysicsFunction<T_S,T_P>* m_qoiFunction;


      // SURROGATE VARIABLES
      enum SurrogateModelType{
        PSEUDO_SPECTRAL_TENSOR_PRODUCT=0,
        PSEUDO_SPECTRAL_SPARSE_GRID,
        PSEUDO_SPECTRAL_MONTE_CARLO };

      int                       m_surrogateType;
      std::vector<unsigned int> m_order;
      SurrogateModel<T_S,T_P>*  m_primalSurrogate;
      SurrogateModel<T_S,T_P>*  m_adjointSurrogate;
      SurrogateModel<T_S,T_P>*  m_qoiSurrogate;

      // OUTPUT VARIABLES
  };

/********************************************//**
 * \brief 
 * 
 ***********************************************/
  Driver::Driver(
      GetPot& input
      ) : m_input(input)
  {
    
    // DRIVER SETTINGS
    m_input.set_prefix("driver/");
    m_maxIter = m_input("maxIter",1);
    
    
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

    
    
    // PHYSICS MODEL SETTINGS
    std::string physicsType  = 
      m_input("physics/type","catenary");

    if (physicsType == "catenary")
      m_physicsType =  PHYSICS_CATENARY ;
    else if (physicsType == "burgers")
      m_physicsType =  PHYSICS_VISCOUS_BURGERS ;
    else if (physicsType == "grins")
      m_physicsType =  PHYSICS_GRINS ;
    else
    {
      std::cerr << "unrecognized PhysicsModelType " 
        << physicsType << std::endl;
      exit(1);
    }

    initializePhysicsModel( m_physicsType, m_input );
    

    // SURROGATE MODEL SETTINGS
    std::string surrType  = 
      m_input("surrogateModels/type","PseudoSpectralTensorProduct");

    if (surrType == "PseudoSpectralTensorProduct")
      m_surrogateType =  PSEUDO_SPECTRAL_TENSOR_PRODUCT ;
    else if (surrType == "PseudoSpectralSparseGrid")
      m_surrogateType =  PSEUDO_SPECTRAL_SPARSE_GRID ;
    else if (surrType == "PseudoSpectralMonteCarlo")
      m_surrogateType =  PSEUDO_SPECTRAL_MONTE_CARLO ;
    else
    {
      std::cerr << "unrecognized SurrogateModelType " 
        << surrType << std::endl;
      exit(1);
    }

    initializeSurrogateModel( m_surrogateType, m_input );


    

    // OUTPUT DATA SETTINGS
    // TODO much later

  }


/********************************************//**
 * \brief 
 * 
 ***********************************************/
  Driver::~Driver( )
  {
  }

/********************************************//**
 * \brief an initial driver run routine for testing
 * 
 ***********************************************/
  void Driver::run( )
  {
    //TODO
    std::cout << "maxIter = " << m_maxIter << std::endl;
    std::cout << "dimension = " << m_paramDim << std::endl;
    std::cout << "order = " ;
      for(unsigned int i=0; i < m_paramDim; i++)
        std::cout << m_order[i] << " " ;
    std::cout << std::endl;
    return;
  }


/********************************************//**
 * \brief 
 ***********************************************/
  void Driver::initializeSurrogateModel(
      int surrogateType,
      GetPot& input)
  {
    // TODO multiple surrogate models for each sol (primal,adjoint,qoi,etc)
    input.set_prefix("surrogateModels/");
    // Setup common to all surrogate models

    // type specific setup
    switch( m_surrogateType )
    {
      case(PSEUDO_SPECTRAL_TENSOR_PRODUCT):
        {
          // TODO do we really want to define a sub section in the input file ??
          input.set_prefix("surrogateModels/PseudoSpectralTensorProduct/");

          for (unsigned int i=0; i < m_paramDim; i++)
            m_order.push_back( m_input("order", 0, i) ) ;

          m_primalSurrogate = new PseudoSpectralTensorProduct<T_S,T_P>(
                *m_primalPhysics, 
                m_parameters, 
                m_order  );
          m_adjointSurrogate = new PseudoSpectralTensorProduct<T_S,T_P>(
                *m_adjointPhysics, 
                m_parameters, 
                m_order  );
          m_qoiSurrogate = new PseudoSpectralTensorProduct<T_S,T_P>(
                *m_qoiFunction, 
                m_parameters, 
                m_order  );

          break;
        }
      case(PSEUDO_SPECTRAL_SPARSE_GRID):
        {
          std::cerr << " this SurrogateModelType is not yet implemented\n" ;
          exit(1);
          break;
        }
      case(PSEUDO_SPECTRAL_MONTE_CARLO):
        {
          std::cerr << " this SurrogateModelType is not yet implemented\n" ;
          exit(1);
          break;
        }
    }


  }

/********************************************//**
 * \brief 
 ***********************************************/
  void Driver::initializePhysicsModel(
      int physicsType,
      GetPot& input)
  {
    input.set_prefix("physics/");
    // setup common to all physics models
    
    // type-specific physics setup
    switch( m_physicsType )
    {
      input.set_prefix("physics/catenary/");
      case(PHYSICS_CATENARY):
        {

          m_physicsModel = new PhysicsCatenary<T_S,T_P>( input("forcing",-10.0) );
          m_primalPhysics = new PhysicsFunctionPrimal<T_S,T_P>( 
              *m_physicsModel) ;
          m_adjointPhysics = new PhysicsFunctionAdjoint<T_S,T_P>( 
              *m_physicsModel) ;
          m_qoiFunction = new PhysicsFunctionQoi<T_S,T_P>( 
              *m_physicsModel) ;

          break;
        }
      case(PHYSICS_VISCOUS_BURGERS):
        {
          std::cerr << " this PhysicsModelType is not yet implemented\n" ;
          exit(1);
          break;
        }
      case(PHYSICS_GRINS):
        {
          std::cerr << " this PhysicsModelType is not yet implemented\n" ;
          exit(1);
          break;
        }
    }

  }

}


#endif // DRIVER_H


