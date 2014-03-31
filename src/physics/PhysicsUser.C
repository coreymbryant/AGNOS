


#include "PhysicsUser.h"

namespace AGNOS
{
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
      std::set<std::string>& computeSolutions,
      const T_S& paramVector, 
      std::map<std::string, T_P >& solutionVectors 
      )
  { 
    if(AGNOS_DEBUG)
      std::cout << "DEBUG: calling provided compute function\n" ;

    // set parameter values
    this->_setParameterValues( paramVector );

    if(_compute_function != NULL)
      this->_compute_function(computeSolutions, paramVector, solutionVectors);
    else
    {
      std::cerr << std::endl;
      std::cerr 
        << "ERROR: No compute function is defined for the PhysicsUser class.\n"
        << "       Either redefine compute or use attach_compute_function(..)\n" 
        << "       to attach a user defined compute routine. \n"  ;
      std::cerr << std::endl;
      exit(1);
    }

    if(AGNOS_DEBUG)
      std::cout << "DEBUG: ending call to provided compute function\n" ;

  }

  template class
    PhysicsUser<libMesh::DenseVector<double>, libMesh::DenseVector<double> >;

}

