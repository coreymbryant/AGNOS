

#ifndef PHYSICS_CATENARY_LIBMESH_H
#define PHYSICS_CATENARY_LIBMESH_H

#include "agnosDefines.h"
#include "PhysicsModel.h"
#include "PhysicsAssembly.h"

// libmesh includes
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
/* #include "libmesh/error_vector.h" */
/* #include "libmesh/kelly_error_estimator.h" */
#include "libmesh/mesh_refinement.h"
        

namespace AGNOS
{

  /********************************************//**
   * \brief Example PhysicsModel class - catenary chain (solution computed using
   * libmesh)
   *
   * A simple 1D example useful for testing purposes
   *
   * 
   ***********************************************/
  template<class T_S, class T_P>
  class PhysicsCatenaryLibmesh : public PhysicsModel<T_S,T_P>
  {

    public:
      PhysicsCatenaryLibmesh( 
          const Communicator&       comm,
          const GetPot& input
          );
      ~PhysicsCatenaryLibmesh( );


      T_P solvePrimal( 
          const T_S& parameterValue  
          );

      T_P solveAdjoint( 
          const T_S& parameterValue,  
          const T_P& primalSolution    
          );

      T_P evaluateQoi( 
          const T_S& parameterValue,
          const T_P& primalSolution    
          );

      T_P estimateError( 
          const T_S& parameterValue,  
          const T_P& primalSolution,   
          const T_P& adjointSolution  
          );
      
    protected:
      double          m_forcing;
      double          m_min;
      double          m_max;
      unsigned int    m_n;    // number of elements
      unsigned int    m_maxRefineSteps;

      libMesh::Mesh                   m_mesh; // mesh
      libMesh::EquationSystems*       m_equation_systems;
      libMesh::LinearImplicitSystem*  m_system;
      libMesh::MeshRefinement*        m_mesh_refinement;

      PhysicsAssembly<T_S>*           m_physicsAssembly;

  };




/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsCatenaryLibmesh<T_S,T_P>::PhysicsCatenaryLibmesh(
      const Communicator&       comm,
      const GetPot&             physicsInput
      ) 
  {
    //-------- get physics model settings
    if (comm.rank() == 0)
      std::cout << "\n-- Reading Physics model data\n";

    m_min             = physicsInput("physics/min",-1.);
    m_max             = physicsInput("physics/max",1.);
    m_forcing         = physicsInput("physics/forcing",-10.);
    m_n               = physicsInput("physics/n",2);
    m_maxRefineSteps  = physicsInput("physics/maxRefineSteps",1);


    //-------- crate libmesh mesh object
    if (comm.rank() == 0)
      std::cout << "\n-- Creating mesh\n";

    libMesh::MeshTools::Generation::build_line(m_mesh,m_n,m_min,m_max,EDGE2);

    
    //-------- create equation system
    if (comm.rank() == 0)
      std::cout << "\n-- Setting up equation system and refinement strategy.\n";

    // define equation system
    m_equation_systems = new libMesh::EquationSystems(m_mesh);
    m_system = &( m_equation_systems->add_system<LinearImplicitSystem>("1D") ) ;

    // add variable (first order)
    m_system->add_variable("u",FIRST);

    // provide pointer to assemly routine
    // TODO define this
    m_physicsAssembly = new PhysicsAssembly<T_S>( 
        *m_equation_systems, "1D", m_forcing);
    const T_S tempParamValues(1);
    m_physicsAssembly->setParameterValues( tempParamValues );
    m_system->attach_assemble_object( *m_physicsAssembly );


    //TODO mesh refinement stuff
    // define mesh refinement object
    m_mesh_refinement = new libMesh::MeshRefinement(m_mesh);
    m_mesh_refinement->refine_fraction()  = 0.7;
    m_mesh_refinement->coarsen_fraction() = 0.3;
    m_mesh_refinement->max_h_level()      = 5;

    //------ initialize data structures
    m_equation_systems->init();

  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsCatenaryLibmesh<T_S,T_P>::~PhysicsCatenaryLibmesh( )
  {
    delete m_equation_systems;
    delete m_mesh_refinement;
    delete m_physicsAssembly;
  }




/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    T_P PhysicsCatenaryLibmesh<T_S,T_P>::solvePrimal( 
        const T_S& parameterValue  
        )
    {
      T_P imageValue;
      /* m_currentCoeff = parameterValue(0); */
      m_physicsAssembly->setParameterValues( parameterValue );

      m_equation_systems->reinit();
      m_equation_systems->get_system("1D").solve();
      
      return imageValue;
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    T_P PhysicsCatenaryLibmesh<T_S,T_P>::solveAdjoint( 
        const T_S& parameterValue,  
        const T_P& primalSolution    
        )
    {

      T_P imageValue(1);

      imageValue(0) =  1.0 / (4. * parameterValue(0))  ;
      return imageValue;
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    T_P PhysicsCatenaryLibmesh<T_S,T_P>::evaluateQoi( 
        const T_S& parameterValue,
        const T_P& primalSolution    
        )
    {
      libMesh::Point evalPoint(0.5);
      libMesh::Number qoiValue = m_system->point_value(
          0, evalPoint);
      
      T_P returnVec;
      returnVec.resize(1);
      returnVec(0) = qoiValue;
      return returnVec;
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    T_P PhysicsCatenaryLibmesh<T_S,T_P>::estimateError( 
        const T_S& parameterValue,  
        const T_P& primalSolution,   
        const T_P& adjointSolution  
        )
    {
      // in this case the FE solution interpolates at x=1/2 so the QoI is
      // evaluated exactly
      T_P imageValue(1);
      imageValue.zero();

      return imageValue;
    }



}

#endif // PHYSICS_CATENARY_LIBMESH_H
