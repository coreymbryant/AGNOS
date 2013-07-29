
#ifndef PHYSICS_VISCOUS_BURGERS_H
#define PHYSICS_VISCOUS_BURGERS_H


#include "agnosDefines.h"
#include "PhysicsLibmesh.h"
#include "ResidualViscousBurgers.h"
#include "JacobianViscousBurgers.h"
#include "QoiViscousBurgers.h"
#include "QoiDerivativeViscousBurgers.h"

// libmesh includes
#include "libmesh/mesh_generation.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/petsc_nonlinear_solver.h"

namespace AGNOS
{


  /********************************************//**
   * \brief Basic 1D Burger's PhysicsModel class
   *
   * This example is given in the book
   * "Spectral Methods for Uncertainty Quantification" by Le Maitre and Knio
   ***********************************************/
  template<class T_S, class T_P>
  class PhysicsViscousBurgers : public PhysicsLibmesh<T_S,T_P>
  {

    public:
      PhysicsViscousBurgers( 
        Parallel::Communicator& comm_in, 
        const GetPot& input );

      ~PhysicsViscousBurgers( );

      void setParameterValues( const T_S& parameterValues ) ;

    protected:
      double          _L;
      double          _uMinus;
      double          _uPlus;
      int _n;

      PhysicsJacobian<T_S>*           _physicsJacobian;
      PhysicsResidual<T_S>*           _physicsResidual;

      // Initialization
      void _init( );
      void _setParameterValues( const T_S& parameterValues ) ;

      // solver settings
      unsigned int    _nonlinearTolerance;
      unsigned int    _nonlinearSteps;


  };




/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsViscousBurgers<T_S,T_P>::PhysicsViscousBurgers(
        Parallel::Communicator& comm_in, 
        const GetPot& input 
        ) :
    PhysicsLibmesh<T_S,T_P>(comm_in,input)
  {
    _L      = input("physics/L",10.);
    _uMinus = input("physics/uMinus",(0.5 * ( 1 + std::tanh( -1.*_L / 4. / 1.0) ) ));
    _uPlus  = input("physics/uPlus",(0.5 * ( 1 + std::tanh( _L / 4. / 1.0) ) ) );

    _nonlinearSteps      = input("physics/nNonlinearSteps",15);
    _nonlinearTolerance  = input("physics/nonlinearTolerance",1.e-9);

    _init( );
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsViscousBurgers<T_S,T_P>::~PhysicsViscousBurgers( )
  {
    delete _physicsResidual;
    delete _physicsJacobian;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  void PhysicsViscousBurgers<T_S,T_P>::_init( )
  {
    //----------------------------------------------
    this->_mesh = new Mesh( this->_communicator );

    libMesh::MeshTools::Generation::build_line(
        *this->_mesh, _n, -1.*_L, _L, EDGE3);

    // define equation system
    this->_equationSystems 
      = new libMesh::EquationSystems(*this->_mesh);

    // add system and set parent pointer 
    this->_system = 
      &( this->_equationSystems->add_system(
            "NonlinearImplicit","Burgers") );

    // add variable (first order)
    this->_system->add_variable("u",FIRST);
    //----------------------------------------------

    //----------------------------------------------
    //---- set up nonlinear solver
    // TODO do we need this?
    // initalize the nonlinear solver
    /* this->_equationSystems->template */
    /*   get_system<NonlinearImplicitSystem>("Burgers").nonlinear_solver->init(); */

    //----------------------------------------------
    // TODO set from input file ?
    // set solver settings
    this->_equationSystems->parameters.template set<Real>
      ("nonlinear solver relative residual tolerance") = _nonlinearTolerance;
    this->_equationSystems->parameters.template set<Real>
      ("nonlinear solver absolute residual tolerance") = 1.e-35;
    this->_equationSystems->parameters.template set<Real>
      ("nonlinear solver absolute step tolerance") = 1.e-12;
    this->_equationSystems->parameters.template set<Real>
      ("nonlinear solver relative step tolerance") = 1.e-12;
    this->_equationSystems->parameters.template set<unsigned int>
      ("nonlinear solver maximum iterations") = 50;
    //----------------------------------------------
    

    //----------------------------------------------
    // provide pointer to residual object
    _physicsResidual = new ResidualViscousBurgers<T_S>( );
    _physicsResidual->setSystemData( this->_input );
    this->_equationSystems->template
      get_system<NonlinearImplicitSystem>("Burgers").nonlinear_solver->residual_object 
      = _physicsResidual;

    // provide pointer to jacobian object
    _physicsJacobian = new JacobianViscousBurgers<T_S>( );
    _physicsJacobian->setSystemData( this->_input );
    this->_equationSystems->template
      get_system<NonlinearImplicitSystem>("Burgers").nonlinear_solver->jacobian_object
      = _physicsJacobian ;
    //----------------------------------------------


    //----------------------------------------------
    // QoISet
    this->_qois = new libMesh::QoISet;
    std::vector<unsigned int> qoi_indices;
    qoi_indices.push_back(0);
    this->_qois->add_indices(qoi_indices);
    this->_qois->set_weight(0, 1.0);
                            
    //pointer to qoi
    this->_qoi = new QoiViscousBurgers<T_S>(
        *this->_equationSystems, "Burgers");

    // pointer to qoi derivative assembly
    this->_system->qoi.resize(1);
    this->_qoiDerivative = new QoiDerivativeViscousBurgers<T_S>(
        *this->_equationSystems, "Burgers");
    //----------------------------------------------

    /* std::cout << "test: system initialized" << std::endl; */
  }

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S, class T_P>
    void PhysicsViscousBurgers<T_S,T_P>::setParameterValues( 
        const T_S& parameterValues ) 
    {
      _physicsResidual->setParameterValues( parameterValues );
      _physicsJacobian->setParameterValues( parameterValues );
      return;
    }

}

#endif // PHYSICS_VISCOUS_BURGERS_H
