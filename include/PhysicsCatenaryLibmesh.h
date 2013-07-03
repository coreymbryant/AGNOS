

#ifndef PHYSICS_CATENARY_LIBMESH_H
#define PHYSICS_CATENARY_LIBMESH_H

#include "agnosDefines.h"
#include "PhysicsModel.h"
#include "PhysicsAssembly.h"
#include "PhysicsQoi.h"
#include "PhysicsQoiDerivative.h"

// libmesh includes
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/error_vector.h"
#include "libmesh/adjoint_refinement_estimator.h"
#include "libmesh/mesh_refinement.h"
        

namespace AGNOS
{

  libMesh::AutoPtr<libMesh::AdjointRefinementEstimator>
    build_adjoint_refinement_error_estimator(libMesh::QoISet &qois)
    {
      libMesh::AutoPtr<libMesh::AdjointRefinementEstimator> error_estimator;

      /* std::cout */
      /*   <<"Computing the error estimate using the Adjoint Refinement Error Estimator" */ 
      /*   <<std::endl<<std::endl; */

      libMesh::AdjointRefinementEstimator *adjoint_refinement_estimator = new
        libMesh::AdjointRefinementEstimator;

      error_estimator.reset(adjoint_refinement_estimator);

      adjoint_refinement_estimator->qoi_set() = qois;
      adjoint_refinement_estimator->number_h_refinements = 0;
      adjoint_refinement_estimator->number_p_refinements = 1;
                    
      return error_estimator;
    }

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

      void refine ( libMesh::ErrorVector errorIndicators );
      
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
      PhysicsQoi<T_S>*                m_qoi;
      PhysicsQoiDerivative<T_S>*      m_qoiDerivative;
      libMesh::QoISet*                m_qois;

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

    libMesh::MeshTools::Generation::build_line(m_mesh,m_n,m_min,m_max,EDGE3);

    
    //-------- create equation system
    if (comm.rank() == 0)
      std::cout << "\n-- Setting up equation system and refinement strategy.\n";

    // define equation system
    m_equation_systems = new libMesh::EquationSystems(m_mesh);
    m_system = &( m_equation_systems->add_system<LinearImplicitSystem>("1D") ) ;

    // add variable (first order)
    m_system->add_variable("u",FIRST);

    // provide pointer to assemly routine
    m_physicsAssembly = new PhysicsAssembly<T_S>( 
        *m_equation_systems, "1D", m_forcing);
    const T_S tempParamValues(1);
    m_physicsAssembly->setParameterValues( tempParamValues );
    m_system->attach_assemble_object( *m_physicsAssembly );

    // QoISet
    m_qois = new libMesh::QoISet;
    std::vector<unsigned int> qoi_indices;
    qoi_indices.push_back(0);
    m_qois->add_indices(qoi_indices);
    m_qois->set_weight(0, 1.0);
                            
    //pointer to qoi
    m_qoi = new PhysicsQoi<T_S>(
        *m_equation_systems, "1D");
    m_system->attach_QOI_object( *m_qoi );

    // pointer to qoi derivative assembly
    m_system->qoi.resize(1);
    m_qoiDerivative = new PhysicsQoiDerivative<T_S>(
        *m_equation_systems, "1D");
    m_system->attach_QOI_derivative_object( *m_qoiDerivative );

    // define mesh refinement object
    m_mesh_refinement = new libMesh::MeshRefinement(m_mesh);
    m_mesh_refinement->refine_fraction()  = 0.7;
    m_mesh_refinement->coarsen_fraction() = 0.3;
    m_mesh_refinement->max_h_level()      = 15;
    m_mesh_refinement->absolute_global_tolerance() = 1e-3;

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
    delete m_qoi;
    delete m_qoiDerivative;
    delete m_qois;
  }




/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    T_P PhysicsCatenaryLibmesh<T_S,T_P>::solvePrimal( 
        const T_S& parameterValue  
        )
    {

      // reference to system 
      /* libMesh::LinearImplicitSystem& system = */
      /*   m_equation_systems->get_system<LinearImplicitSystem>("1D"); */

      // solve system
      m_physicsAssembly->setParameterValues( parameterValue );
      m_system->solve();

      // convert solution to T_P framework
      std::set<libMesh::dof_id_type> dofIndices;
      m_system->local_dof_indices( 
          m_system->variable_number("u"), dofIndices);

      T_P imageValue(dofIndices.size());
      std::set<libMesh::dof_id_type>::iterator dofIt = dofIndices.begin();
      for (unsigned int i=0; dofIt != dofIndices.end(); ++dofIt, i++)
      {
        imageValue(i) =  m_system->current_solution(*dofIt) ;
      }

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
      // reference to system 
      /* libMesh::LinearImplicitSystem& system = */
      /*   m_equation_systems->get_system<LinearImplicitSystem>("1D"); */

      // set solution to provided value
      /* m_physicsAssembly->setParameterValues( parameterValue ); */
      for (unsigned int i=0; i<m_system->solution->size(); i++)
        m_system->solution->set(i, primalSolution(i) );
      m_system->solution->close();

      // solve adjoint
      m_system->adjoint_solve( );
      
      // convert solution to T_P
      libMesh::NumericVector<double>& libmeshSol 
        = m_system->get_adjoint_solution( ) ;

      T_P imageValue( libmeshSol.size() );
      for (unsigned int i=0; i< libmeshSol.size(); i++)
        imageValue(i) = libmeshSol(i);

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
      // reference to system 
      /* libMesh::LinearImplicitSystem& system = */
      /*   m_equation_systems->get_system<LinearImplicitSystem>("1D"); */
      
      // set solution to provided value
      /* m_physicsAssembly->setParameterValues( parameterValue ); */
      for (unsigned int i=0; i<m_system->solution->size(); i++)
        m_system->solution->set(i, primalSolution(i) );
      m_system->solution->close();

      // evaluate QoI and get value
      m_system->assemble_qoi();
      std::vector< libMesh::Number > qoiValue = m_system->qoi;

      // covert to output format
      T_P returnVec( qoiValue.size() );
      for (unsigned int i=0; i<qoiValue.size(); i++)
        returnVec(i) = qoiValue[i];

      /* std::cout << "qoi = " << returnVec(0) << std::endl; */

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
      // set solution to provided value
      for (unsigned int i=0; i<m_system->solution->size(); i++)
      {
        m_system->solution->set(i, primalSolution(i) );
        /* std::cout << "primal(" << i << ") = " << (*m_system->solution)(i) << std::endl; */
      }
      m_system->solution->close();
      
      // set adjoint solutions
      libMesh::NumericVector<libMesh::Number>* adjSolution 
        = &( m_system->get_adjoint_solution(0) ) ;
      for (unsigned int i=0; i<adjointSolution.size(); i++)
      {
        adjSolution->set(i, adjointSolution(i) );
        /* std::cout << "adjoint(" << i << ") = " << (*adjSolution)(i) << std::endl; */
      }
      adjSolution->close();

      // error indicators
      libMesh::AutoPtr<libMesh::AdjointRefinementEstimator> 
        adjoint_refinement_error_estimator =
        build_adjoint_refinement_error_estimator(*m_qois);

      libMesh::ErrorVector qoiErrorIndicators;

      adjoint_refinement_error_estimator->estimate_error(*m_system,
          qoiErrorIndicators);

      /* std::cout << "errorEst = " << qoiErrorIndicators->l2_norm() */
      /*   << std::endl; */


      // copy to solution vector
      T_P imageValue(qoiErrorIndicators.size());
      for (unsigned int i=0; i<imageValue.size(); i++)
      {
        imageValue(i) = qoiErrorIndicators[i];
        std::cout << "error(" << i << ") = " << qoiErrorIndicators[i] <<
          std::endl;
      }

      return imageValue;
    }


  template<class T_S, class T_P>
    void PhysicsCatenaryLibmesh<T_S,T_P>::refine(
        libMesh::ErrorVector errorIndicators
        ) 
    {


      //TODO need to use some sort of average as error indicators
      /* m_mesh_refinement->flag_elements_by_error_tolerance */
      /*    (this->m_meanErrorIndicator); */              
      /* m_mesh_refinement->flag_elements_by_error_fraction */
      /*    (this->m_meanErrorIndicator); */              
                   
      /* m_mesh_refinement->refine_and_coarsen_elements(); */
      /* m_mesh_refinement->refine_elements(); */

      m_mesh_refinement->uniformly_refine(1);

      m_equation_systems->reinit();

      std::cout << "nDofs = " << m_equation_systems->n_active_dofs() << std::endl;


      return;
    }

}

#endif // PHYSICS_CATENARY_LIBMESH_H
