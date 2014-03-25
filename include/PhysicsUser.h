

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
        const T_S& paramVector, std::map<std::string, T_P >& solutionVectors ); 

    /** attach a user defined compute function */
    void attach_compute_function(
        void fptr(const T_S& paramVector,std::map<std::string,T_P>& solutionVectors)
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
        const T_S& paramVector, 
        std::map<std::string, T_P >& solutionVectors ) ;
  };




/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsUser<T_S,T_P>::PhysicsUser(
        const Communicator& comm_in,
        const GetPot& input
      ) :
    PhysicsModel<T_S,T_P>(comm_in,input),
    _compute_function(NULL)
  {
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsUser<T_S,T_P>::~PhysicsUser( )
  {
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  void PhysicsUser<T_S,T_P>::compute( 
      const T_S& paramVector, 
      std::map<std::string, T_P >& solutionVectors 
      )
  { 
    if(AGNOS_DEBUG)
      std::cout << "DEBUG: calling provided compute function\n" ;

    // set parameter values
    this->_setParameterValues( paramVector );

    if(_compute_function != NULL)
      this->_compute_function(paramVector, solutionVectors);
    else
    {
      std::cerr << std::endl;
      std::cerr 
        << "ERROR: No compute function is defined for the PhysicsUser class."
        << "       Either redefine compute or use attach_compute_function(..)" 
        << "       to attach a user defined compute routine. "  ;
      std::cerr << std::endl;
      exit(1);
    }

    if(AGNOS_DEBUG)
      std::cout << "DEBUG: ending call to provided compute function\n" ;

  }

}

#endif // PHYSICS_USER_H
