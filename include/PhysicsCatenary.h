

#ifndef PHYSICS_CATENARY_H
#define PHYSICS_CATENARY_H

#include "PhysicsModel.h"

namespace AGNOS
{

  /********************************************//**
   * \brief Example PhysicsModel class - catenary chain
   *
   * A simple, single dof, system useful for testing purposes
   *
   * 1 dof corresponds to middle node in two (linear) element approximation on
   * domain [0,1]
   *
   * 
   ***********************************************/
  template<class T_S, class T_P>
  class PhysicsCatenary : public PhysicsModel<T_S,T_P>
  {

  public:
    /** Default constructor Default */
    PhysicsCatenary( 
        const Communicator& comm_in,
        const GetPot& input
        );

    /** Destructor */
    ~PhysicsCatenary( );

    /** Function called by SurrogateModel to solve for requested solution
     * vectors at each evaluation point in parameter space  */
    virtual void compute( 
        const T_S& paramVector, 
        std::map<std::string, T_P >& solutionVectors 
        ) ;

    /** Refinement methods for the physics model */
    virtual void refine( );
  
  protected:
    /** communicator reference */
    const Communicator &_communicator;
    /** reference to input file just incase its needed after initialization */
    const GetPot&         _input;

    /** forcing magnitude */
    double m_forcing;
  };




/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsCatenary<T_S,T_P>::PhysicsCatenary(
        const Communicator& comm_in,
        const GetPot& input
      ) :
    PhysicsModel<T_S,T_P>(comm_in,input),
    _communicator(comm_in),_input(input)
  {
    m_forcing = input("forcing",-10.0) ;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsCatenary<T_S,T_P>::~PhysicsCatenary( )
  {
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    void PhysicsCatenary<T_S,T_P>::compute(
        const T_S& paramVector,
        std::map<std::string, T_P >& solutionVectors 
        ) 
    {
      if(AGNOS_DEBUG)
        std::cout << "entering compute routine" << std::endl;
      solutionVectors.clear();

      T_P primalSolution(1);
      primalSolution(0) =  m_forcing / (8. *  paramVector(0) ) ;
      solutionVectors.insert( 
          std::pair<std::string,T_P>( "primal", primalSolution )
            );
      if(AGNOS_DEBUG)
        std::cout << "compute -- primal solution size:" << primalSolution.size()
          << std::endl;

      T_P adjointSolution(1);
      adjointSolution(0) =  1.0 / (4. * paramVector(0))  ;
      solutionVectors.insert( 
          std::pair<std::string,T_P>( "adjoint", adjointSolution )
            );
      
      T_P qoiValue(1);
      qoiValue(0) = primalSolution(0);
      solutionVectors.insert( 
                std::pair<std::string,T_P>( "qoi", qoiValue )
                  );

      // in this case the FE solution interpolates at x=1/2 so the QoI is
      // evaluated exactly
      T_P error(1);
      error.zero();
      solutionVectors.insert( 
                std::pair<std::string,T_P>( "errorEstimate", error )
                  );

      if(AGNOS_DEBUG)
      {
        std::cout << "solutionVectors['primal'].size():" 
          << solutionVectors["primal"].size() << std::endl;

        std::cout << "leaving compute routine" << std::endl;
      }
    }


  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S,class T_P> 
    void PhysicsCatenary<T_S,T_P>::refine( ) 
    {
      return;
    }


}

#endif // PHYSICS_CATENARY_H
