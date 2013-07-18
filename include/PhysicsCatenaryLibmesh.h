

#ifndef PHYSICS_CATENARY_LIBMESH_H
#define PHYSICS_CATENARY_LIBMESH_H

#include "agnosDefines.h"
#include "PhysicsLibmesh.h"
#include "AssemblyCatenary.h"
#include "QoiCatenary.h"
#include "QoiDerivativeCatenary.h"

// libmesh includes
#include "libmesh/mesh_generation.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/linear_implicit_system.h"
        

namespace AGNOS
{

  /********************************************//**
   * \brief Example PhysicsLibmesh class - catenary chain (solution computed
   * using libmesh)
   *
   * A simple 1D example useful for testing purposes
   *
   * 
   ***********************************************/
  template<class T_S, class T_P>
  class PhysicsCatenaryLibmesh : public PhysicsLibmesh<T_S,T_P>
  {

    public:
      PhysicsCatenaryLibmesh( 
          const Communicator&       comm,
          const GetPot& input
          );
      ~PhysicsCatenaryLibmesh( );

      void setParameterValues( const T_S& parameterValues ) ;

    protected:
      PhysicsAssembly<T_S>*           m_physicsAssembly;

      double          m_forcing;
      double          m_min;
      double          m_max;

      void _initializeSystem( );
      void _constructMesh( );


  };




/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsCatenaryLibmesh<T_S,T_P>::PhysicsCatenaryLibmesh(
      const Communicator&       comm,
      const GetPot&             physicsInput
      ) : PhysicsLibmesh<T_S,T_P>(comm,physicsInput)
  {
    m_min             = physicsInput("physics/min",0.);
    m_max             = physicsInput("physics/max",1.);
    m_forcing         = physicsInput("physics/forcing",-10.);

    this->init( );
    

  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsCatenaryLibmesh<T_S,T_P>::~PhysicsCatenaryLibmesh( )
  {
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  void PhysicsCatenaryLibmesh<T_S,T_P>::_constructMesh( )
  {
    libMesh::MeshTools::Generation::build_line(
        this->m_mesh,this->m_n,m_min,m_max,EDGE3);

    return;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  void PhysicsCatenaryLibmesh<T_S,T_P>::_initializeSystem( )
  {
    // define equation system
    this->m_equation_systems 
      = new libMesh::EquationSystems(this->m_mesh);
    this->m_system = &( 
        this->m_equation_systems->add_system("LinearImplicit", "1D")
        ) ;

    // add variable (first order)
    this->m_system->add_variable("u",FIRST);

    // provide pointer to assemly routine
    m_physicsAssembly = new AssemblyCatenary<T_S>( 
        *this->m_equation_systems, "1D");
    m_physicsAssembly->setSystemData( this->m_input );

    // provide pointer to assemly routine
    this->m_system->attach_assemble_object( *m_physicsAssembly );

    // QoISet
    this->m_qois = new libMesh::QoISet;
    std::vector<unsigned int> qoi_indices;
    qoi_indices.push_back(0);
    this->m_qois->add_indices(qoi_indices);
    this->m_qois->set_weight(0, 1.0);
                            
    //pointer to qoi
    this->m_qoi = new QoiCatenary<T_S>(
        *this->m_equation_systems, "1D");

    // pointer to qoi derivative assembly
    this->m_system->qoi.resize(1);
    this->m_qoiDerivative = new QoiDerivativeCatenary<T_S>(
        *this->m_equation_systems, "1D");
  }

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S, class T_P>
    void PhysicsCatenaryLibmesh<T_S,T_P>::setParameterValues( 
        const T_S& parameterValues ) 
    {
      m_physicsAssembly->setParameterValues( parameterValues );
      return;
    }


}

#endif // PHYSICS_CATENARY_LIBMESH_H
