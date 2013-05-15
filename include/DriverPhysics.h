
#ifndef DRIVER_PHYSICS_H
#define DRIVER_PHYSICS_H


#include "Driver.h"


namespace AGNOS
{

/********************************************//**
 * \brief Driver based on user defined PhysicsModel class
 ***********************************************/
  class DriverPhysics : public Driver
  { 

    public:

      // either physics model or physics function must be defined
      DriverPhysics( 
          PhysicsModel<T_S,T_P>*    myPhysics,          
          GetPot&                   input );

      virtual ~DriverPhysics( );

      void run( ) ;

    protected:


  };


/********************************************//**
 * \brief 
 * 
 ***********************************************/
  DriverPhysics::DriverPhysics( 
      PhysicsModel<T_S,T_P>*    myPhysics,          
      GetPot&                   input 
      ) : Driver( input )
  {
    
    // TODO multiple surrogate models for each sol (primal,adjoint,qoi,etc)
    input.set_prefix("surrogateModel/");
    // Setup common to all surrogate models
    
    std::vector< PhysicsFunction<T_S,T_P>* > physicsFunctions;

    // primal solution
    physicsFunctions.push_back( 
      new PhysicsFunctionPrimal<T_S,T_P>( *myPhysics ) ); 

    // adjoint solution
    physicsFunctions.push_back( 
      new PhysicsFunctionAdjoint<T_S,T_P>( *myPhysics ) ); 

    // qoi evaluation
    physicsFunctions.push_back( 
      new PhysicsFunctionQoi<T_S,T_P>( *myPhysics ) ); 
    


    // type specific setup
    switch( m_surrogateType )
    {
      case(PSEUDO_SPECTRAL_TENSOR_PRODUCT):
        {

          m_surrogate = new PseudoSpectralTensorProduct<T_S,T_P>(
                physicsFunctions, 
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

  }

/********************************************//**
 * \brief 
 * 
 ***********************************************/
  DriverPhysics::~DriverPhysics( )
  {
    for (unsigned int i=0; i < physicsFunctions.size(); i++)
      delete physicsFunctions[i];
    delete m_surrogate;
  }

/********************************************//**
 * \brief an initial driver run routine for testing
 * 
 ***********************************************/
  void DriverPhysics::run( )
  {
    //TODO
    std::cout << "maxIter = " << m_maxIter << std::endl;
    std::cout << "dimension = " << m_paramDim << std::endl;
    std::cout << "order = " ;
    for(unsigned int i=0; i < m_paramDim; i++)
      std::cout << m_order[i] << " " ;
    std::cout << std::endl;

    std::cout << "mins " ;
    for (unsigned int i=0; i < m_paramDim; i++)
      std::cout << m_parameters[i]->min() << " " ;
    std::cout << std::endl;

    std::cout << "maxs " ;
    for (unsigned int i=0; i < m_paramDim; i++)
      std::cout << m_parameters[i]->max() << " " ;
    std::cout << std::endl;


    m_surrogate->build();




    return;
  }





}


#endif // DRIVER_PHYSICS_H


