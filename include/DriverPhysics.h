
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
          const GetPot&             input );

      virtual ~DriverPhysics( );

      void run( ) ;

    protected:
      PhysicsModel<T_S,T_P>*   m_myPhysics;
      std::map< std::string, PhysicsFunction<T_S,T_P>* > m_physicsFunctions ;
      std::map< std::string, PhysicsFunction<T_S,T_P>* > m_errorFunctions ;


  };


/********************************************//**
 * \brief 
 * 
 ***********************************************/
  DriverPhysics::DriverPhysics( 
      const Communicator& comm,
      PhysicsModel<T_S,T_P>*    myPhysics,          
      const GetPot&                   input 
      ) : Driver( comm, input ), m_myPhysics(myPhysics)
  {
    
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
    


    std::map< std::string, PhysicsFunction<T_S,T_P>* > m_errorFunctions;
    // error estimate
    m_errorFunctions.insert( 
        std::pair< std::string, PhysicsFunction<T_S,T_P>* >(
          "error", 
          new PhysicsFunctionError<T_S,T_P>( *myPhysics ) ) 
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

          this->m_errorSurrogate = new PseudoSpectralTensorProduct<T_S,T_P>(
              this->m_comm,
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
    delete m_errorSurrogate;
  }

/********************************************//**
 * \brief an initial driver run routine for testing
 * 
 ***********************************************/
  void DriverPhysics::run( )
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

    // refine approximation
    for (unsigned int iter=2; iter <= this->m_maxIter; iter++)
    {
      std::cout << "\n-------------  ITER "
        << iter << "  -------------\n " ;
      
      // TODO control what type of refinement to perform
      /* m_surrogate->refine(); */
      /* m_errorSurrogate->refine(); */


      // use mean of errorSurrogate as error indicators
      std::map< std::string, std::vector<T_P> >   errorCoeff = m_errorSurrogate->getCoefficients() ;
      libMesh::ErrorVector errorIndicators(errorCoeff["error"][0].size()) ;
      for (int i=0; i<errorCoeff["error"][0].size(); i++)
      {
        errorIndicators[i] = errorCoeff["error"][0](i);
        std::cout << "error(" << i << ") = " << errorIndicators[i] << std::endl;
      }

      m_myPhysics->refine( errorIndicators );
      m_surrogate->build();
      m_errorSurrogate->build();



      /* double totalError = errorIndicators->l2_norm(); */
      /* std::cout << "totalError = " << totalError << std::endl; */

      if (this->m_outputIterations)
      {
        std::cout << "\n writing results to: " << this->m_outputFilename
          << " (iter = " << iter << " )"
          << std::endl;
        std::cout << std::endl;
        printSolution(iter);
      }
    }
    
    // output whatever user asks for
    if (this->m_comm->rank() == 0)
    {
      std::cout << "\n writing final results to: " << this->m_outputFilename
        << std::endl;
      std::cout << std::endl;
    }
      printSolution(m_maxIter);



    // evaluate QoI
      /* T_S evalPoint(1); */
      /* evalPoint(0) = 1.5; */
      /* T_P qoiValue = m_surrogate->evaluate("qoi", evalPoint ); */
    /* if (this->m_comm->rank() == 0) */
    /* { */
      /* std::cout << "\n Qoi = " << qoiValue(0) << std::endl; */
    /* } */



    return;
  }





}


#endif // DRIVER_PHYSICS_H


