
#ifndef PHYSICS_VISCOUS_BURGERS_H
#define PHYSICS_VISCOUS_BURGERS_H


#include "agnosDefines.h"
#include "PhysicsLibmesh.h"
#include "PhysicsAssembly.h"
#include "PhysicsQoi.h"
#include "PhysicsQoiDerivative.h"

// libmesh includes
#include "libmesh/mesh_generation.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/linear_implicit_system.h"

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
          const Communicator&       comm,
          const GetPot& input
          );
      ~PhysicsViscousBurgers( );

    protected:
      double          m_min;
      double          m_max;

      void _initializeSystem( );
      void _constructMesh( );


  };




/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsViscousBurgers<T_S,T_P>::PhysicsViscousBurgers(
      const Communicator&       comm,
      const GetPot&             physicsInput
      ) : PhysicsLibmesh(comm,physicsInput)
  {
    m_min                 = physicsInput("physics/min",-10.);
    m_max                 = physicsInput("physics/max",10.);
    m_nonlinearSteps      = physicsInput("physics/nNonlinearSteps",15);
    m_nonlineartolerance  = physicsInput("physics/nonlinearTolerance",1.e-3);

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
        this->m_mesh,this->m_n,m_min,m_max,EDGE3);

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
    this->m_system = &( 
        this->m_equation_systems->add_system("NonlinearImplicit", "1D")
        ) ;

    // add variable (first order)
    this->m_system->add_variable("u",FIRST);

    // provide pointer to assemly routine
    this->m_physicsAssembly = new PhysicsAssembly<T_S>( 
        *this->m_equation_systems, "1D");
    const T_S tempParamValues(1);
    this->m_physicsAssembly->setParameterValues( tempParamValues );


    // QoISet
    this->m_qois = new libMesh::QoISet;
    std::vector<unsigned int> qoi_indices;
    qoi_indices.push_back(0);
    this->m_qois->add_indices(qoi_indices);
    this->m_qois->set_weight(0, 1.0);
                            
    //pointer to qoi
    this->m_qoi = new PhysicsQoi<T_S>(
        *this->m_equation_systems, "1D");

    // pointer to qoi derivative assembly
    this->m_system->qoi.resize(1);
    this->m_qoiDerivative = new PhysicsQoiDerivative<T_S>(
        *this->m_equation_systems, "1D");
  }





}

#endif // PHYSICS_VISCOUS_BURGERS_H
