

#ifndef PHYSICS_USER_H
#define PHYSICS_USER_H

#include "PhysicsModel.h"

namespace AGNOS
{

  /********************************************//**
   * \brief Derived PhysicsModel class that takes care of pure virtuals and
   * allows user defined compute function to be set. Note: this class should
   * only be used for simple test problems, any sophisticated model should be
   * derived directly from the base PhysicsModel class. See PhysicsLibmesh as an
   * example. 
   ***********************************************/
  template<class T_S, class T_P>
  class PhysicsUser : public PhysicsModel<T_S,T_P>
  {

  public:
    /** Default constructor Default */
    PhysicsUser( 
        const Communicator& comm_in,
        const GetPot& input
        );

    /** Destructor */
    ~PhysicsUser( );

    /** Function called by SurrogateModel to solve for requested solution
     * vectors at each evaluation point in parameter space  */
    virtual void compute( 
        std::set<std::string>& computeSolutions,
        const T_S& paramVector, std::map<std::string, T_P >& solutionVectors ); 

    /** attach a user defined compute function */
    void attach_compute_function(
        void fptr(
          std::set<std::string>& computeSolutions,
          const T_S& paramVector,std::map<std::string,T_P>& solutionVectors)
        )
    {
      libmesh_assert(fptr);
      _compute_function = fptr;
    }
  
  protected:
    /** _setParameterValues must be redefined because it is a pure virtual in
     * PhysicsModel. But it does nothing in this class. */
    void _setParameterValues( const T_S& parameterValues ){} ;

    /** Function pointer to compute function */
    void(* _compute_function)(
        std::set<std::string>& computeSolutions,
        const T_S& paramVector, 
        std::map<std::string, T_P >& solutionVectors ) ;
  };




}

#endif // PHYSICS_USER_H
