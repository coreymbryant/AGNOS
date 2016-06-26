

#ifndef PHYSICS_CATENARY_H
#define PHYSICS_CATENARY_H

#include "PhysicsUser.h"

namespace AGNOS
{

  /********************************************//**
   * \brief Example PhysicsUser class - catenary chain
   *
   * A simple, single dof, system useful for testing purposes
   *
   * 1 dof corresponds to middle node in two (linear) element approximation on
   * domain [0,1]
   *
   * 
   ***********************************************/
  template<class T_S, class T_P>
  class PhysicsCatenary : public PhysicsUser<T_S,T_P>
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
        std::set<std::string>& computeSolutions,
        const T_S& paramVector, std::map<std::string, T_P >& solutionVectors 
        ) ;

    /** Refinement methods for the physics model */
    virtual void refine( );

    /** get number of degrees of freedom */
    virtual unsigned int getNDofs( ) { return 1; }

  protected:
    /** communicator reference */
    const Communicator &_communicator;
    /** reference to input file just incase its needed after initialization */
    const GetPot&         _input;

    /** forcing magnitude */
    double m_forcing;

  };


}

#endif // PHYSICS_CATENARY_H
