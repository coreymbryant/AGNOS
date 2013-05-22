
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
          const Communicator&       comm,
          PhysicsModel<T_S,T_P>*    myPhysics,          
          GetPot&                   input );

      virtual ~DriverPhysics( );

      void run( ) ;

    protected:
      std::map< std::string, PhysicsFunction<T_S,T_P>* > m_physicsFunctions ;


  };


/********************************************//**
 * \brief 
 * 
 ***********************************************/
  DriverPhysics::DriverPhysics( 
      const Communicator& comm,
      PhysicsModel<T_S,T_P>*    myPhysics,          
      GetPot&                   input 
      ) : Driver( comm, input )
  {
    
    // TODO multiple surrogate models for each sol (primal,adjoint,qoi,etc)
    input.set_prefix("surrogateModel/");
    // Setup common to all surrogate models
    
    std::map< std::string, PhysicsFunction<T_S,T_P>* > m_physicsFunctions;

    // primal solution
    m_physicsFunctions.insert( 
        std::pair< std::string, PhysicsFunction<T_S,T_P>* >(
          "primal", 
          new PhysicsFunctionPrimal<T_S,T_P>( *myPhysics ) ) 
        );

    // adjoint solution
    m_physicsFunctions.insert( 
        std::pair< std::string, PhysicsFunction<T_S,T_P>* >(
          "adjoint", 
          new PhysicsFunctionAdjoint<T_S,T_P>( *myPhysics ) ) 
        ); 

    // qoi evaluation
    m_physicsFunctions.insert( 
        std::pair< std::string, PhysicsFunction<T_S,T_P>* >(
          "qoi", 
          new PhysicsFunctionQoi<T_S,T_P>( *myPhysics ) ) 
        ); 
    


    // type specific setup
    switch( m_surrogateType )
    {
      case(PSEUDO_SPECTRAL_TENSOR_PRODUCT):
        {

          this->m_surrogate = new PseudoSpectralTensorProduct<T_S,T_P>(
              this->m_comm,
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

  }

/********************************************//**
 * \brief 
 * 
 ***********************************************/
  DriverPhysics::~DriverPhysics( )
  {
    std::map< std::string, PhysicsFunction<T_S,T_P>* >::iterator id;
    for (id=m_physicsFunctions.begin(); id !=m_physicsFunctions.end(); id++)
      delete id->second;
    delete m_surrogate;
  }

/********************************************//**
 * \brief an initial driver run routine for testing
 * 
 ***********************************************/
  void DriverPhysics::run( )
  {

    // TODO add some cout statmenents for progress 
    
    // build initial approximation
    m_surrogate->build();
    
    // output whatever user asks for
    printOutput();

    // refine approximation
    for (unsigned int iter=2; iter <= this->m_maxIter; iter++)
    {
      // TODO may need to pass some info to refine
      m_surrogate->refine();

      if (this->m_outputIterations)
      {
        printOutput();
      }
    }
    



    return;
  }





}


#endif // DRIVER_PHYSICS_H


