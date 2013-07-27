
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
          const GetPot& input
          );
      ~PhysicsViscousBurgers( );

      void setParameterValues( const T_S& parameterValues ) ;

    protected:
      PhysicsJacobian<T_S>*           m_physicsJacobian;
      PhysicsResidual<T_S>*           m_physicsResidual;


      double          m_L;
      double          m_uMinus;
      double          m_uPlus;

      // solver settings
      unsigned int    m_nonlinearTolerance;
      unsigned int    m_nonlinearSteps;

      void _initializeSystem( );
      void _constructMesh( );


  };




/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsViscousBurgers<T_S,T_P>::PhysicsViscousBurgers(
      const GetPot&             physicsInput
      ) : PhysicsLibmesh<T_S,T_P>(physicsInput)
  {
    m_L                 = physicsInput("physics/L",10.);
    m_uMinus                 = physicsInput("physics/uMinus",(0.5 * ( 1 + std::tanh( -1.*m_L / 4. / 1.0) ) ));
    m_uPlus                 = physicsInput("physics/uPlus",(0.5 * ( 1 + std::tanh( m_L / 4. / 1.0) ) ) );
    m_nonlinearSteps      = physicsInput("physics/nNonlinearSteps",15);
    m_nonlinearTolerance  = physicsInput("physics/nonlinearTolerance",1.e-9);


    this->init( );
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsViscousBurgers<T_S,T_P>::~PhysicsViscousBurgers( )
  {
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  void PhysicsViscousBurgers<T_S,T_P>::_constructMesh( )
  {
    libMesh::MeshTools::Generation::build_line(
        this->m_mesh,this->m_n,-1.*m_L,m_L,EDGE3);

    return;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  void PhysicsViscousBurgers<T_S,T_P>::_initializeSystem( )
  {
    // define equation system
    this->m_equation_systems 
      = new libMesh::EquationSystems(this->m_mesh);

    // add system and set parent pointer 
    this->m_system = 
      &( this->m_equation_systems->add_system(
            "NonlinearImplicit","Burgers") );

    // add variable (first order)
    this->m_system->add_variable("u",FIRST);

    //---- set up nonlinear solver
    // TODO do we need this?
    // initalize the nonlinear solver
    /* this->m_equation_systems->template */
    /*   get_system<NonlinearImplicitSystem>("Burgers").nonlinear_solver->init(); */

    // TODO set from input file ?
    // set solver settings
    this->m_equation_systems->parameters.template set<Real>
      ("nonlinear solver relative residual tolerance") = m_nonlinearTolerance;
    this->m_equation_systems->parameters.template set<Real>
      ("nonlinear solver absolute residual tolerance") = 1.e-35;
    this->m_equation_systems->parameters.template set<Real>
      ("nonlinear solver absolute step tolerance") = 1.e-12;
    this->m_equation_systems->parameters.template set<Real>
      ("nonlinear solver relative step tolerance") = 1.e-12;
    this->m_equation_systems->parameters.template set<unsigned int>
      ("nonlinear solver maximum iterations") = 50;
    
    // provide pointer to residual object
    m_physicsResidual = new ResidualViscousBurgers<T_S>( );
    m_physicsResidual->setSystemData( this->m_input );
    this->m_equation_systems->template
      get_system<NonlinearImplicitSystem>("Burgers").nonlinear_solver->residual_object 
      = m_physicsResidual;

    // provide pointer to jacobian object
    m_physicsJacobian = new JacobianViscousBurgers<T_S>( );
    m_physicsJacobian->setSystemData( this->m_input );
    this->m_equation_systems->template
      get_system<NonlinearImplicitSystem>("Burgers").nonlinear_solver->jacobian_object
      = m_physicsJacobian ;

    // QoISet
    this->m_qois = new libMesh::QoISet;
    std::vector<unsigned int> qoi_indices;
    qoi_indices.push_back(0);
    this->m_qois->add_indices(qoi_indices);
    this->m_qois->set_weight(0, 1.0);
                            
    //pointer to qoi
    this->m_qoi = new QoiViscousBurgers<T_S>(
        *this->m_equation_systems, "Burgers");

    // pointer to qoi derivative assembly
    this->m_system->qoi.resize(1);
    this->m_qoiDerivative = new QoiDerivativeViscousBurgers<T_S>(
        *this->m_equation_systems, "Burgers");
    /* std::cout << "test: system initialized" << std::endl; */
  }


  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S, class T_P>
    void PhysicsViscousBurgers<T_S,T_P>::setParameterValues( 
        const T_S& parameterValues ) 
    {
      m_physicsResidual->setParameterValues( parameterValues );
      m_physicsJacobian->setParameterValues( parameterValues );
      return;
    }

}

#endif // PHYSICS_VISCOUS_BURGERS_H
