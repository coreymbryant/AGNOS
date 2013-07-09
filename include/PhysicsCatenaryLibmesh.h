

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
      adjoint_refinement_estimator->number_h_refinements = 2;
      adjoint_refinement_estimator->number_p_refinements = 0;

                    

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

      void refine (  );
      void refine ( libMesh::ErrorVector errorIndicators );
      
    protected:
      const Communicator*   m_comm;
      double          m_forcing;
      double          m_min;
      double          m_max;
      unsigned int    m_n;    // number of elements
      unsigned int    m_maxRefineSteps;

      libMesh::Mesh                   m_mesh; // mesh
      libMesh::EquationSystems*       m_equation_systems;
      libMesh::LinearImplicitSystem*  m_system;
      libMesh::MeshRefinement*        m_mesh_refinement;
      libMesh::AutoPtr<libMesh::AdjointRefinementEstimator> 
        m_adjointRefinementErrorEstimator ;

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
      ) : m_comm(&comm)
  {
    //-------- get physics model settings
    if (comm.rank() == 0)
      std::cout << "\n-- Reading Physics model data\n";

    m_min             = physicsInput("physics/min",0.);
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

    std::cout << "test" << std::endl;
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

    m_mesh_refinement->coarsen_by_parents()         = true;
    m_mesh_refinement->absolute_global_tolerance()  = 0.0;
    m_mesh_refinement->nelem_target()               = 21;  
    m_mesh_refinement->refine_fraction()            = 0.7;
    m_mesh_refinement->coarsen_fraction()           = 0.3;  
    m_mesh_refinement->coarsen_threshold()          = 5;
    m_mesh_refinement->max_h_level()                = 15;
    
    // error indicators
    m_adjointRefinementErrorEstimator =
      build_adjoint_refinement_error_estimator(*m_qois);

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
      m_system->reinit();
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
      //TODO may need to refine and project solution (like what is done in the
      //adjoint refinement estimator)
      // reference to system 
      /* libMesh::LinearImplicitSystem& system = */
      /*   m_equation_systems->get_system<LinearImplicitSystem>("1D"); */

      // set solution to provided value
      /* m_physicsAssembly->setParameterValues( parameterValue ); */
      for (unsigned int i=0; i<m_system->solution->size(); i++)
        m_system->solution->set(i, primalSolution(i) );
      m_system->solution->close();

      // We have to break the rules here, because we can't refine a const System
      /* LinearImplicitSystem& system = const_cast<LinearImplicitSystem&>(*m_system); */

      // An EquationSystems reference will be convenient.
      EquationSystems& es = m_system->get_equation_systems();

      // The current mesh
      MeshBase& mesh = es.get_mesh();

      // We'll want to back up all coarse grid vectors
      std::map<std::string, NumericVector<Number> *> coarse_vectors;
      for (System::vectors_iterator vec = m_system->vectors_begin(); vec !=
           m_system->vectors_end(); ++vec)
        {
          // The (string) name of this vector
          const std::string& var_name = vec->first;

          coarse_vectors[var_name] = vec->second->clone().release();
        }
      // Back up the coarse solution and coarse local solution
      NumericVector<Number> * coarse_solution =
        m_system->solution->clone().release();
      NumericVector<Number> * coarse_local_solution =
        m_system->current_local_solution->clone().release();
      // And make copies of the projected solution
      NumericVector<Number> * projected_solution;

      // And we'll need to temporarily change solution projection settings
      bool old_projection_setting;
      old_projection_setting = m_system->project_solution_on_reinit();

      // Make sure the solution is projected when we refine the mesh
      m_system->project_solution_on_reinit() = true;

      // And it'll be best to avoid any repartitioning
      AutoPtr<Partitioner> old_partitioner = mesh.partitioner();
      mesh.partitioner().reset(NULL);

      // And we can't allow any renumbering
      const bool old_renumbering_setting = mesh.allow_renumbering();
      mesh.allow_renumbering(false);


#ifndef NDEBUG
      // n_coarse_elem is only used in an assertion later so
      // avoid declaring it unless asserts are active.
      const dof_id_type n_coarse_elem = mesh.n_elem();
#endif

      // Uniformly refine the mesh
      MeshRefinement mesh_refinement(mesh);



      libmesh_assert (m_adjointRefinementErrorEstimator->number_h_refinements > 0
          || m_adjointRefinementErrorEstimator->number_p_refinements > 0);

      // FIXME: this may break if there is more than one System
      // on this mesh but estimate_error was still called instead of
      // estimate_errors
      for (unsigned int i = 0; i !=
          m_adjointRefinementErrorEstimator->number_h_refinements; ++i)
        {
          mesh_refinement.uniformly_refine(1);
          es.reinit();
        }

      for (unsigned int i = 0; i !=
          m_adjointRefinementErrorEstimator->number_p_refinements; ++i)
        {
          mesh_refinement.uniformly_p_refine(1);
          es.reinit();
        }

      // Copy the projected coarse grid solutions, which will be
      // overwritten by solve()
      projected_solution = NumericVector<Number>::build(mesh.comm()).release();
      projected_solution->init(m_system->solution->size(), true, SERIAL);
      m_system->solution->localize(*projected_solution,
              m_system->get_dof_map().get_send_list());

      // Rebuild the rhs with the projected primal solution
      (dynamic_cast<ImplicitSystem&>(*m_system)).assembly(true, false);
      NumericVector<Number> & projected_residual = (dynamic_cast<ExplicitSystem&>(*m_system)).get_vector("RHS Vector");
      projected_residual.close();

      // solve adjoint
      m_system->adjoint_solve( );
      m_system->set_adjoint_already_solved(true);

      
      // Don't bother projecting the solution; we'll restore from backup
      // after coarsening
      m_system->project_solution_on_reinit() = false;

      // Uniformly coarsen the mesh, without projecting the solution
      libmesh_assert (m_adjointRefinementErrorEstimator->number_h_refinements > 0
          || m_adjointRefinementErrorEstimator->number_p_refinements > 0);

      for (unsigned int i = 0; i !=
          m_adjointRefinementErrorEstimator->number_h_refinements; ++i)
        {
          mesh_refinement.uniformly_coarsen(1);
          // FIXME - should the reinits here be necessary? - RHS
          es.reinit();
        }

      for (unsigned int i = 0; i !=
          m_adjointRefinementErrorEstimator->number_p_refinements; ++i)
        {
          mesh_refinement.uniformly_p_coarsen(1);
          es.reinit();
        }

      /* std::cout << "primalSize: " << m_system->solution->size() << std::endl; */
      /* for(unsigned int i=0; i < m_system->solution->size() ; i++) */
      /* { */
      /*   std::cout << "solution(" << i << ") = " << (*m_system->solution)(i) << */
      /*     std::endl; */
      /* } */
      /* std::cout << "adjointSize: " << m_system->get_adjoint_solution().size() << std::endl; */
      /* for(unsigned int i=0; i < m_system->get_adjoint_solution().size() ; i++) */
      /* { */
      /*   std::cout << "adjoint(" << i << ") = " << */
      /*     (m_system->get_adjoint_solution())(i) << std::endl; */
      /* } */

      // We should be back where we started
      libmesh_assert_equal_to (n_coarse_elem, mesh.n_elem());

      // Restore old solutions and clean up the heap
      m_system->project_solution_on_reinit() = old_projection_setting;

      // Restore the coarse solution vectors and delete their copies
      *m_system->solution = *coarse_solution;
      delete coarse_solution;
      *m_system->current_local_solution = *coarse_local_solution;
      delete coarse_local_solution;
      delete projected_solution;

      std::map<std::string, NumericVector<Number> *>::iterator vec;
      for (vec = coarse_vectors.begin(); vec != coarse_vectors.end(); ++vec)
        {
          /* std::cout << "var_naem: <" << vec->first << std::endl; */
          // The (string) name of this vector
          const std::string& var_name = vec->first;

          m_system->get_vector(var_name) = *coarse_vectors[var_name];

          coarse_vectors[var_name]->clear();
          delete coarse_vectors[var_name];
        }

      // Restore old partitioner and renumbering settings
      mesh.partitioner() = old_partitioner;
      mesh.allow_renumbering(old_renumbering_setting);


      
      // convert solution to T_P
      libMesh::NumericVector<double>& libmeshSol 
        = m_system->get_adjoint_solution( ) ;

      /* std::cout << "primalSize: " << m_system->solution->size() << std::endl; */
      /* for(unsigned int i=0; i < m_system->solution->size() ; i++) */
      /* { */
      /*   std::cout << "solution(" << i << ") = " << (*m_system->solution)(i) << */
      /*     std::endl; */
      /* } */
      /* std::cout << "adjointSize: " << m_system->get_adjoint_solution().size() << std::endl; */
      /* for(unsigned int i=0; i < m_system->get_adjoint_solution().size() ; i++) */
      /* { */
      /*   std::cout << "adjoint(" << i << ") = " << */
      /*     (m_system->get_adjoint_solution())(i) << std::endl; */
      /* } */

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
      m_system->update();

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
      m_system->set_adjoint_already_solved(true);


      libMesh::ErrorVector qoiErrorIndicators(m_mesh.n_elem());

      //  modified estimate_error routine and rebuilding libmesh
      m_adjointRefinementErrorEstimator->estimate_error(
          *m_system,
          qoiErrorIndicators);

      /* std::cout << "errorEst = " << qoiErrorIndicators->l2_norm() */
      /*   << std::endl; */


      // copy to solution vector
      T_P imageValue(qoiErrorIndicators.size());
      for (unsigned int i=0; i<imageValue.size(); i++)
      {
        imageValue(i) = qoiErrorIndicators[i];
        /* std::cout << "error(" << i << ") = " << qoiErrorIndicators[i] << */
        /*   std::endl; */
      }

      return imageValue;
    }

  /********************************************//**
   * \brief 
   *
   * 
   ***********************************************/
  template<class T_S, class T_P>
    void PhysicsCatenaryLibmesh<T_S,T_P>::refine() 
    {
      m_mesh_refinement->uniformly_refine(1);

      m_equation_systems->reinit();
    }

  /********************************************//**
   * \brief 
   *
   * 
   ***********************************************/
  template<class T_S, class T_P>
    void PhysicsCatenaryLibmesh<T_S,T_P>::refine(
        libMesh::ErrorVector errorIndicators
        ) 
    {
      if (m_comm->rank() == 0)
      {
        std::cout << " error.size = " << errorIndicators.size() << std::endl;
        std::cout << "nElem = " << m_mesh.n_active_elem() << std::endl;
      }

      /* for(unsigned int i=0; i<errorIndicators.size(); i++) */
      /*   std::cout << "error(" << i << ") = " << errorIndicators[i] << std::endl; */
      /* for(unsigned int j=0; j < this->m_primalSolution->size(); j++) */
      /*   std::cout << "primal(" << j << ") = " << (*this->m_primalSolution)(j) << */
      /*     std::endl; */
      /* for(unsigned int k=0; k < this->m_adjointSolution->size(); k++) */
      /*   std::cout << "adjoint(" << k << ") = " << (*this->m_adjointSolution)(k) << */
      /*     std::endl; */

      //TODO need to use some sort of average as error indicators
      /* m_mesh_refinement->clean_refinement_flags(); */
      /* m_mesh_refinement->flag_elements_by_error_tolerance (errorIndicators); */              
      /* m_mesh_refinement->flag_elements_by_error_fraction (errorIndicators); */              
      /* m_mesh_refinement->flag_elements_by_elem_fraction (errorIndicators); */              
                   
      /* std::cout << "unflagged: " << m_mesh_refinement->test_unflagged( ) << std::endl; */
      /* m_mesh_refinement->refine_and_coarsen_elements(); */
      /* m_mesh_refinement->refine_elements(); */

      m_mesh_refinement->uniformly_refine(1);

      m_equation_systems->reinit();





      return;
    }

}

#endif // PHYSICS_CATENARY_LIBMESH_H
