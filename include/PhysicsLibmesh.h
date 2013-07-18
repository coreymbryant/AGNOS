#ifndef PHYSICS_LIBMESH_H
#define PHYSICS_LIBMESH_H

#include "agnosDefines.h"
#include "PhysicsModel.h"
#include "PhysicsAssembly.h"
#include "PhysicsJacobian.h"
#include "PhysicsResidual.h"
#include "PhysicsAssembly.h"
#include "PhysicsQoi.h"
#include "PhysicsQoiDerivative.h"

// libmesh includes
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/error_vector.h"
#include "libmesh/adjoint_refinement_estimator.h"
#include "libmesh/mesh_refinement.h"
        

namespace AGNOS
{

  /********************************************//**
   * \brief Base libmesh physics model class. All libmesh physics, and maybe
   * GRINS, PhysicsModel classes should be derived from here.
   *
   ***********************************************/
  template<class T_S, class T_P>
  class PhysicsLibmesh : public PhysicsModel<T_S,T_P>
  {

    public:
      PhysicsLibmesh( 
          const Communicator& comm, //< global comm (not libmesh object specific)
          const GetPot&       input
          );
      ~PhysicsLibmesh( );


      T_P solvePrimal( 
          const T_S& parameterValue  
          ) ;

      T_P solveAdjoint( 
          const T_S& parameterValue,  
          const T_P& primalSolution    
          ) ;

      T_P evaluateQoi( 
          const T_S& parameterValue,
          const T_P& primalSolution    
          ) ;

      T_P estimateError( 
          const T_S& parameterValue,  
          const T_P& primalSolution,   
          const T_P& adjointSolution  
          ) ;

      void refine (  ) ;
      void refine ( T_P& errorIndicators ) ;
      
      libMesh::Mesh getMesh( ) ;

      void init( ) ;

      virtual void setParameterValues( const T_S& parameterValues ) = 0;


    protected:
      const Communicator*   m_comm; //< global comm (not libmesh object specific)
      const GetPot&               m_input;

      // mesh and equation variables
      unsigned int    m_n;              // number of elements
      libMesh::Mesh                   m_mesh; // mesh
      libMesh::EquationSystems*       m_equation_systems;
      libMesh::System*                m_system;
      libMesh::MeshRefinement*        m_mesh_refinement;

      // refinement options
      unsigned int m_numberHRefinements;
      unsigned int m_numberPRefinements;
      unsigned int m_maxRefineSteps;

      // problem specific routines
      PhysicsQoi<T_S>*                m_qoi;
      PhysicsQoiDerivative<T_S>*      m_qoiDerivative;
      libMesh::QoISet*                m_qois;


      virtual void _constructMesh( ) = 0;
      virtual void _initializeSystem( ) = 0;

  };




/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsLibmesh<T_S,T_P>::PhysicsLibmesh(
      const Communicator&       comm,
      const GetPot&             physicsInput
      ) : m_comm(&comm), m_input(physicsInput)
  {
    //-------- get physics model settings
    if (comm.rank() == 0)
      std::cout << "\n-- Reading Physics model data\n";

    // read in mesh settings
    m_n               = physicsInput("physics/n",1);

    // read refinement options
    this->m_useUniformRefinement =
      physicsInput("physics/useUniformRefinement",true);
    m_numberHRefinements = physicsInput("physics/numberHRefinements",0);
    m_numberPRefinements = physicsInput("physics/numberPRefinements",1);
    m_maxRefineSteps  = physicsInput("physics/maxRefineSteps",1);

    // other options
    this->m_resolveAdjoint =
      physicsInput("physics/resolveAdjoint",false);
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  void PhysicsLibmesh<T_S,T_P>::init( )
  {

    //-------- create libmesh mesh object
    if (m_comm->rank() == 0)
      std::cout << "\n-- Creating mesh\n";
    _constructMesh();

    
    //-------- create equation system
    if (m_comm->rank() == 0)
      std::cout << "\n-- Setting up equation system and refinement strategy.\n";
    _initializeSystem();


    //-------- set up specific routines

    //pointer to qoi
    m_system->attach_QOI_object( *m_qoi );

    // pointer to qoi derivative assembly
    m_system->attach_QOI_derivative_object( *m_qoiDerivative );

    // define mesh refinement object
    m_mesh_refinement = new libMesh::MeshRefinement(m_mesh);

    // TODO read these from input file
    m_mesh_refinement->coarsen_by_parents()         = true;
    m_mesh_refinement->absolute_global_tolerance()  = 1e-6;
    /* m_mesh_refinement->nelem_target()               = 64; */  
    m_mesh_refinement->refine_fraction()            = 0.7;
    m_mesh_refinement->coarsen_fraction()           = 0.3;  
    m_mesh_refinement->coarsen_threshold()          = 1e-5;
    m_mesh_refinement->max_h_level()                = 15;
    
    /* std::cout << "test: mesh info " << std::endl; */
    m_mesh.print_info();


    // warning for nonlinear problems
    bool hasNonlinearSystem = false;
    for(unsigned int i=0; i<m_equation_systems->n_systems(); i++)
    {
      if (m_equation_systems->get_system(i).system_type() ==
          "NonlinearImplicit")
      {
        hasNonlinearSystem = true;
        break;
      }
    }

    if (hasNonlinearSystem && !this->m_resolveAdjoint)
    {
      if(this->m_comm->rank() == 0 )
        std::cout << std::endl << " WARNING: "
          << "You may want to set resolveAdjoint = true for nonlinear "
          << "problems.  \n" << std::endl;
    }

    //------ initialize data structures
    m_equation_systems->init();
    /* std::cout << "test: ES info " << std::endl; */
    m_equation_systems->print_info();
    /* std::cout << "test: solver info " << std::endl; */
    /* std::cout << "test: post init " << std::endl; */
    /* std::cout << "test: default solver package: " */ 
      /* << libMesh::default_solver_package() << std::endl; */

  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsLibmesh<T_S,T_P>::~PhysicsLibmesh( )
  {
    delete m_equation_systems;
    delete m_mesh_refinement;
    delete m_qoi;
    delete m_qoiDerivative;
    delete m_qois;
  }




/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    T_P PhysicsLibmesh<T_S,T_P>::solvePrimal( 
        const T_S& parameterValue  
        )
    {
      /* std::cout << "test: solvePrimal begin" << std::endl; */
      /* std::cout << "test:       solution.size(): " << */
      /*   m_system->solution->size() << std::endl; */
      /* std::cout << "test:               n_dofs(): " << */
      /*   m_system->n_dofs() << std::endl; */
      /* std::cout << "test:        n_active_dofs(): " << */
      /*   m_system->n_active_dofs() << std::endl; */

      // solve system
      this->setParameterValues( parameterValue );
      /* std::cout << "test: pre reinit" << std::endl; */
      m_system->init();
      /* std::cout << "test: pre solve" << std::endl; */
      m_system->solve();

      /* std::cout << "test:       solution.size(): " << */
      /*   m_system->solution->size() << std::endl; */
      /* std::cout << "test:               n_dofs(): " << */
      /*   m_system->n_dofs() << std::endl; */
      /* std::cout << "test:        n_active_dofs(): " << */
      /*   m_system->n_active_dofs() << std::endl; */

      /* std::cout << "test: final residual:" << std::endl; */
      /*   m_equation_systems->template get_system<NonlinearImplicitSystem>("Burgers").nonlinear_solver->print_converged_reason() ; */

      // convert solution to T_P framework
      std::set<libMesh::dof_id_type> dofIndices;
      m_system->local_dof_indices( 0, dofIndices);

      /* std::cout << "test:      dofIndices.size(): " << */
      /*   dofIndices.size() << std::endl; */
      /* std::cout << "test:               n_dofs(): " << */
      /*   m_system->n_dofs() << std::endl; */
      /* std::cout << "test:        n_active_dofs(): " << */
      /*   m_system->n_active_dofs() << std::endl; */
      /* std::cout << "test:        n_local _dofs(): " << */
      /*   m_system->n_local_dofs() << std::endl; */

      /* for (unsigned int i=0; i < m_system->n_active_dofs(); i++) */
      /* { */
      /*   std::cout << "solution(" << i << ")=" */ 
      /*     << m_system->current_solution(i) */
      /*     << std::endl; */
      /* } */

      T_P imageValue(dofIndices.size());
      std::set<libMesh::dof_id_type>::iterator dofIt 
        = dofIndices.begin();
      for (unsigned int i=0; dofIt != dofIndices.end(); ++dofIt, i++)
      {
        /* std::cout << "solution(dof)(" << i << ")=" */ 
        /*   << m_system->current_solution(*dofIt) */
        /*   << std::endl; */
        imageValue(i) =  m_system->current_solution(*dofIt) ;
      }

      return imageValue;
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    T_P PhysicsLibmesh<T_S,T_P>::solveAdjoint( 
        const T_S& parameterValue,  
        const T_P& primalSolution    
        )
    {
      /* std::cout << "test: solveAdjoint begin" << std::endl; */

      /* std::cout << " test: pre set parameter values " << std::endl; */
      this->setParameterValues(parameterValue);

      m_system->set_adjoint_already_solved(false);
      
      m_system->update();

      /* std::cout << "test:       solution.size(): " << */
      /*   m_system->solution->size() << std::endl; */
      /* std::cout << "test:               n_dofs(): " << */
      /*   m_system->n_dofs() << std::endl; */
      /* std::cout << "test:        n_active_dofs(): " << */
      /*   m_system->n_active_dofs() << std::endl; */
      /* std::cout << "test: primalSolution.size(): " << */
      /*   primalSolution.size() << std::endl; */

      // set solution to provided value
      std::set<libMesh::dof_id_type> dofIndices;
      m_system->local_dof_indices( 0, dofIndices);
      std::set<libMesh::dof_id_type>::iterator dofIt 
        = dofIndices.begin();
      for (unsigned int i=0; dofIt != dofIndices.end(); ++dofIt, i++)
        m_system->solution->set(*dofIt, primalSolution(i) );
      m_system->solution->close();
      /* std::cout << "test: post setPrimal" << std::endl; */

      // An EquationSystems reference will be convenient.
      EquationSystems& es = m_system->get_equation_systems();

      // The current mesh
      MeshBase& mesh = es.get_mesh();

      // We'll want to back up all coarse grid vectors
      std::map<std::string, NumericVector<Number> *> coarse_vectors;
      for (System::vectors_iterator vec = m_system->vectors_begin(); vec !=
           m_system->vectors_end(); ++vec)
        {
          /* std::cout << " var_name: " << vec->first << std::endl; */
          // The (string) name of this vector
          const std::string& var_name = vec->first;

          if (var_name != "adjoint_solution0")
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



      libmesh_assert (m_numberHRefinements > 0 || m_numberPRefinements > 0);

      for (unsigned int i = 0; i != m_numberHRefinements; ++i)
        {
          mesh_refinement.uniformly_refine(1);
          es.reinit();
        }

      for (unsigned int i = 0; i != m_numberPRefinements; ++i)
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


      /* std::cout << " test: pre adjoint_solve " << std::endl; */
      // solve adjoint
      m_system->adjoint_solve( );
      m_system->set_adjoint_already_solved(true);
      /* std::cout << " 336:post adjoint_solve" << std::endl; */
      /* std::cout << " test: post adjoint_solve " << std::endl; */


      // convert solution to T_P
      libMesh::NumericVector<double>& computedAdjointSolution 
        = m_system->get_adjoint_solution( ) ;

      /* std::cout << "test: computedAdjointSolution.size(): " << */
      /*   computedAdjointSolution.size() << std::endl; */
      /* for (unsigned int i=0; i < m_system->n_active_dofs(); i++) */
      /* { */
      /*   std::cout << "adjoint(" << i << ")=" */ 
      /*     << computedAdjointSolution(i) */
      /*     << std::endl; */
      /* } */

      // get dof indices
      m_system->local_dof_indices( 0, dofIndices);
      dofIt = dofIndices.begin();

      /* std::cout << "test:      dofIndices.size(): " << */
      /*   dofIndices.size() << std::endl; */

      T_P imageValue(dofIndices.size());
      for (unsigned int i=0; dofIt != dofIndices.end(); ++dofIt, i++)
      {
        /* std::cout << "adjoint(dof)(" << i << ")=" */ 
        /*   << computedAdjointSolution(*dofIt) */
        /*   << std::endl; */
        imageValue(i) = computedAdjointSolution(*dofIt);
      }


      // Don't bother projecting the solution; we'll restore from backup
      // after coarsening
      m_system->project_solution_on_reinit() = false;

      // Uniformly coarsen the mesh, without projecting the solution
      libmesh_assert (m_numberHRefinements > 0 || m_numberPRefinements > 0);

      for (unsigned int i = 0; i != m_numberHRefinements; ++i)
        {
          mesh_refinement.uniformly_coarsen(1);
          es.reinit();
        }

      for (unsigned int i = 0; i != m_numberPRefinements; ++i)
        {
          mesh_refinement.uniformly_p_coarsen(1);
          es.reinit();
        }

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
          // The (string) name of this vector
          const std::string& var_name = vec->first;

          m_system->get_vector(var_name) = *coarse_vectors[var_name];

          coarse_vectors[var_name]->clear();
          delete coarse_vectors[var_name];
        }

      // Restore old partitioner and renumbering settings
      mesh.partitioner() = old_partitioner;
      mesh.allow_renumbering(old_renumbering_setting);



      /* std::cout << " test: solveAdjoint end " << std::endl; */
      return imageValue;
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    T_P PhysicsLibmesh<T_S,T_P>::evaluateQoi( 
        const T_S& parameterValue,
        const T_P& primalSolution    
        )
    {
      this->setParameterValues(parameterValue);
      
      // set solution to provided value
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


      return returnVec;
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    T_P PhysicsLibmesh<T_S,T_P>::estimateError( 
        const T_S& parameterValue,  
        const T_P& primalSolution,   
        const T_P& adjointSolution  
        )
    {
      /* std::cout << "test: estimateError() beginning" << std::endl; */
      
      m_system->set_adjoint_already_solved(false);
      this->setParameterValues(parameterValue);

      /* std::cout << "test: estimateError() pre set solution" << std::endl; */
      /* std::cout << "   primalSoution.size(): " << primalSolution.size() << */
      /*                                             std::endl; */
      /* std::cout << "   m_system->solution->size(): " */ 
      /*   << m_system->solution->size() << std::endl; */
      
      // set solution to provided value
      std::set<libMesh::dof_id_type> dofIndices;
      m_system->local_dof_indices( 0, dofIndices);
      std::set<libMesh::dof_id_type>::iterator dofIt 
        = dofIndices.begin();
      for (unsigned int i=0; dofIt != dofIndices.end(); ++dofIt, i++)
        m_system->solution->set(*dofIt, primalSolution(i) );
      m_system->solution->close();
      m_system->update();

      /* std::cout << "test: estimateError() post set solution" << std::endl; */
      


      // An EquationSystems reference will be convenient.
      EquationSystems& es = m_system->get_equation_systems();

      // The current mesh
      MeshBase& mesh = es.get_mesh();

      // Resize the error_per_cell vector to be
      // the number of elements, initialized to 0.
      libMesh::ErrorVector error_per_cell;
      error_per_cell.clear();
      error_per_cell.resize (mesh.max_elem_id(), 0.);

      // We'll want to back up all coarse grid vectors
      std::map<std::string, NumericVector<Number> *> coarse_vectors;
      for (System::vectors_iterator vec = m_system->vectors_begin(); vec !=
           m_system->vectors_end(); ++vec)
        {
          // The (string) name of this vector
          const std::string& var_name = vec->first;

          if (var_name != "adjoint_solution0")
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



      libmesh_assert (m_numberHRefinements > 0 || m_numberPRefinements > 0);

      for (unsigned int i = 0; i != m_numberHRefinements; ++i)
        {
          mesh_refinement.uniformly_refine(1);
          es.reinit();
        }

      for (unsigned int i = 0; i != m_numberPRefinements; ++i)
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

      /* std::cout << "test: estimateError() pre set adjiont" << std::endl; */


      // check if system has adjoint_solution vector
      //    this will fail if adjoint hasn't beend solved on this process yet.
      //    For example if only needed for higher order error estimate.
      if (!m_system->have_vector("adjoint_solution0"))
      {
        // if not present initialize it
        m_system->add_adjoint_solution(0);
      }
    
      // set solution to provided value
      libMesh::NumericVector<libMesh::Number>* adjSolution 
        = &( m_system->get_adjoint_solution(0) ) ;
      // shouldn't need to do this
      /* adjSolution->clear(); */
      /* adjSolution->init(adjointSolution.size()); */

      m_system->local_dof_indices( 0, dofIndices);
      dofIt = dofIndices.begin();
      for (unsigned int i=0; dofIt != dofIndices.end(); ++dofIt, i++)
        adjSolution->set(*dofIt, adjointSolution(i) );
      m_system->solution->close();
      adjSolution->close();
      m_system->set_adjoint_already_solved(true);
      m_system->update();


      /* std::cout << "test: estimateError() post set adjiont" << std::endl; */


      // Now that we have the refined adjoint solution and the projected primal solution,
      // we first compute the global QoI error estimate

      // Resize the computed_global_QoI_errors vector to hold the error estimates for each QoI
      std::vector<Number> computed_global_QoI_errors;
      computed_global_QoI_errors.resize(m_system->qoi.size());

      // Loop over all the adjoint solutions and get the QoI error
      // contributions from all of them
      
      
      for (unsigned int j=0; j != m_system->qoi.size(); j++)
        {
          computed_global_QoI_errors[j] = projected_residual.dot(m_system->get_adjoint_solution(j));
          /* std::cout << "computed_global_QoI_errors = " << */
          /*   computed_global_QoI_errors[j] << std::endl; */
        }
      /* std::cout << "test: estimateError() post compute error " << std::endl; */

      // Done with the global error estimates, now construct the element wise error indicators

      // We ought to account for 'spill-over' effects while computing the
      // element error indicators This happens because the same dof is
      // shared by multiple elements, one way of mitigating this is to
      // scale the contribution from each dof by the number of elements it
      // belongs to We first obtain the number of elements each node
      // belongs to

      bool split_shared_dofs = false;

      if (split_shared_dofs) {

      // A map that relates a node id to an int that will tell us how many elements it is a node of
      LIBMESH_BEST_UNORDERED_MAP<dof_id_type, unsigned int>shared_element_count;

      // To fill this map, we will loop over elements, and then in each element, we will loop
      // over the nodes each element contains, and then query it for the number of coarse
      // grid elements it was a node of

      // We will be iterating over all the active elements in the fine mesh that live on
      // this processor
      MeshBase::const_element_iterator elem_it = mesh.active_local_elements_begin();
      const MeshBase::const_element_iterator elem_end = mesh.active_local_elements_end();

      // Keep track of which nodes we have already dealt with
      std::vector<dof_id_type> processed_node_ids;

      // Start loop over elems
      for(; elem_it != elem_end; ++elem_it)
        {
          // Pointer to this element
          const Elem* elem = *elem_it;

          // Loop over the nodes in the element
          for(unsigned int n=0; n != elem->n_nodes(); ++n)
      {
        // Get a pointer to the current node
        Node* node = elem->get_node(n);

        // Get the id of this node
        dof_id_type node_id = node->id();

        // A processed node flag
        bool processed_node = false;

        // Loop over existing processed nodes and see if
        // we have already done this one
        for(std::size_t i = 0; i != processed_node_ids.size(); i++)
          {
            if(node_id == processed_node_ids[i])
        {
          processed_node = true;
        }
          }

        // If we havent already processed this node, do so now
        if(!processed_node)
          {
            // Declare a neighbor_set to be filled by the find_point_neighbors
            std::set<const Elem *> fine_grid_neighbor_set;

            // Call find_point_neighbors to fill the neighbor_set
            elem->find_point_neighbors(*node, fine_grid_neighbor_set);

            // A vector to hold the coarse grid parents neighbors
            std::vector<dof_id_type> coarse_grid_neighbors;

            // Iterators over the fine grid neighbors set
            std::set<const Elem*>::iterator fine_neighbor_it = fine_grid_neighbor_set.begin();
            const std::set<const Elem*>::iterator fine_neighbor_end = fine_grid_neighbor_set.end();

            // Loop over all the fine neighbors of this node
            for(; fine_neighbor_it != fine_neighbor_end ; ++fine_neighbor_it)
        {
          // Pointer to the current fine neighbor element
          const Elem* fine_elem = *fine_neighbor_it;

          // Find the element id for the corresponding coarse grid element
          const Elem* coarse_elem = fine_elem;
          for (unsigned int j = 0; j != m_numberHRefinements; ++j)
            {
              libmesh_assert (coarse_elem->parent());

              coarse_elem = coarse_elem->parent();
            }

          // Loop over the existing coarse neighbors and check if this one is
          // already in there
                      const dof_id_type coarse_id = coarse_elem->id();
          std::size_t j = 0;
          for (; j != coarse_grid_neighbors.size(); j++)
            {
              // If the set already contains this element break out of the loop
              if(coarse_grid_neighbors[j] == coarse_id)
          {
            break;
          }
            }

          // If we didn't leave the loop even at the last element,
          // this is a new neighbour, put in the coarse_grid_neighbor_set
          if(j == coarse_grid_neighbors.size())
            {
              coarse_grid_neighbors.push_back(coarse_id);
            }

        } // End loop over fine neighbors

            // Set the shared_neighbour index for this node to the
            // size of the coarse grid neighbor set
            shared_element_count[node_id] =
              libmesh_cast_int<unsigned int>(coarse_grid_neighbors.size());

            // Add this node to processed_node_ids vector
            processed_node_ids.push_back(node_id);

          } // End if processed node

      } // End loop over nodes

        }  // End loop over elems

      } // if (split_shared_dofs)

      // Get a DoF map, will be used to get the nodal dof_indices for each element
      DofMap &dof_map = m_system->get_dof_map();

      // The global DOF indices, we will use these later on when we compute the element wise indicators
      std::vector<dof_id_type> dof_indices;

      // Localize the global rhs and adjoint solution vectors (which might be shared on multiple processsors) onto a
      // local ghosted vector, this ensures each processor has all the dof_indices to compute an error indicator for
      // an element it owns
      AutoPtr<NumericVector<Number> > localized_projected_residual = NumericVector<Number>::build(m_system->comm());
      localized_projected_residual->init(m_system->n_dofs(), m_system->n_local_dofs(), m_system->get_dof_map().get_send_list(), false, GHOSTED);
      projected_residual.localize(*localized_projected_residual, m_system->get_dof_map().get_send_list());

      // Each adjoint solution will also require ghosting; for efficiency we'll reuse the same memory
      AutoPtr<NumericVector<Number> > localized_adjoint_solution = NumericVector<Number>::build(m_system->comm());
      localized_adjoint_solution->init(m_system->n_dofs(), m_system->n_local_dofs(), m_system->get_dof_map().get_send_list(), false, GHOSTED);

      // We will loop over each adjoint solution, localize that adjoint
      // solution and then loop over local elements
      for (unsigned int i=0; i != m_system->qoi.size(); i++)
        {
          // Skip this QoI if not in the QoI Set
          if (m_qois->has_index(i))
      {
        // Get the weight for the current QoI
        Real error_weight = m_qois->weight(i);

        (m_system->get_adjoint_solution(i)).localize(*localized_adjoint_solution, m_system->get_dof_map().get_send_list());

        // Loop over elements
              MeshBase::const_element_iterator elem_it = mesh.active_local_elements_begin();
              const MeshBase::const_element_iterator elem_end = mesh.active_local_elements_end();

        for(; elem_it != elem_end; ++elem_it)
          {
            // Pointer to the element
            const Elem* elem = *elem_it;

            // Go up number_h_refinements levels up to find the coarse parent
            const Elem* coarse = elem;

            for (unsigned int j = 0; j != m_numberHRefinements; ++j)
              {
          libmesh_assert (coarse->parent());

          coarse = coarse->parent();
        }

                  const dof_id_type e_id = coarse->id();

            // Get the local to global degree of freedom maps for this element
            dof_map.dof_indices (elem, dof_indices);

            // We will have to manually do the dot products.
                  Number local_contribution = 0.;

            for (unsigned int j=0; j != dof_indices.size(); j++)
        {
          // The contribution to the error indicator for this element from the current QoI
          local_contribution += (*localized_projected_residual)(dof_indices[j]) * (*localized_adjoint_solution)(dof_indices[j]);
        }

            // Multiply by the error weight for this QoI
                  local_contribution *= error_weight;

            error_per_cell[e_id] += static_cast<ErrorVectorReal>
              (libmesh_real(local_contribution));

          } // End loop over elements

      } // End if belong to QoI set

        } // End loop over QoIs

      // Don't bother projecting the solution; we'll restore from backup
      // after coarsening
      m_system->project_solution_on_reinit() = false;

      // Uniformly coarsen the mesh, without projecting the solution
      libmesh_assert (m_numberHRefinements > 0 || m_numberPRefinements > 0);

      for (unsigned int i = 0; i != m_numberHRefinements; ++i)
        {
          mesh_refinement.uniformly_coarsen(1);
          es.reinit();
        }

      for (unsigned int i = 0; i !=
          m_numberPRefinements; ++i)
        {
          mesh_refinement.uniformly_p_coarsen(1);
          es.reinit();
        }

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
          // The (string) name of this vector
          const std::string& var_name = vec->first;

          m_system->get_vector(var_name) = *coarse_vectors[var_name];

          coarse_vectors[var_name]->clear();
          delete coarse_vectors[var_name];
        }

      // Restore old partitioner and renumbering settings
      mesh.partitioner() = old_partitioner;
      mesh.allow_renumbering(old_renumbering_setting);


      // computed_global_QoI_errors holds error estimate for each qoi
      // error_per_cell holds error indicators
      
      // set errorIndicators if we need them for adaptive refinement
      if (!this->m_useUniformRefinement)
      {
        this->m_errorIndicators = new T_P(error_per_cell.size());
        for(unsigned int i=0; i<error_per_cell.size(); i++)
          (*this->m_errorIndicators)(i) =  error_per_cell[i] ;

        /* std::cout << "test: rank" << this->m_comm->rank() */ 
        /*   << ": error_per_cell.size(): " << error_per_cell.size() << std::endl; */

      }

      // convert to return type
      T_P errorEstimate(computed_global_QoI_errors.size());
      for (unsigned int i=0; i<errorEstimate.size(); i++)
        errorEstimate(i) = computed_global_QoI_errors[i];


      /* std::cout << "test: estimateError() end" << std::endl; */
      return errorEstimate;
    }

  /********************************************//**
   * \brief 
   *
   * 
   ***********************************************/
  template<class T_S, class T_P>
    void PhysicsLibmesh<T_S,T_P>::refine() 
    {
      if ( this->m_comm->rank() == 0 )
        std::cout << "  previous n_active_elem(): " 
          << m_mesh.n_active_elem() << std::endl;

      m_mesh_refinement->uniformly_refine(1);

      m_equation_systems->reinit();

      if ( this->m_comm->rank() == 0 )
        std::cout << "   refined n_active_elem(): " 
          << m_mesh.n_active_elem() << std::endl;

    }

  /********************************************//**
   * \brief 
   *
   * 
   ***********************************************/
  template<class T_S, class T_P>
    void PhysicsLibmesh<T_S,T_P>::refine(
        T_P& errorIndicators
        ) 
    {

      if ( this->m_comm->rank() == 0 )
        std::cout << "  previous n_active_elem(): " 
          << m_mesh.n_active_elem() << std::endl;

      libMesh::ErrorVector error_per_cell(errorIndicators.size());
      for(unsigned int i=0; i<errorIndicators.size();i++)
        error_per_cell[i] = errorIndicators(i);
      m_mesh_refinement->clean_refinement_flags();

      
      //TODO input options to control type of refinement
      m_mesh_refinement->flag_elements_by_error_tolerance (error_per_cell);            
      /* m_mesh_refinement->flag_elements_by_error_fraction (error_per_cell); */             
      /* m_mesh_refinement->flag_elements_by_elem_fraction (error_per_cell); */              
                   
      //TODO input options to control whether we coarsen as well
      /* m_mesh_refinement->refine_and_coarsen_elements(); */
      m_mesh_refinement->refine_elements();


      m_equation_systems->reinit();

      if ( this->m_comm->rank() == 0 )
        std::cout << "   refined n_active_elem(): " 
          << m_mesh.n_active_elem() << std::endl;

      return;
    }

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S,class T_P> 
    libMesh::Mesh PhysicsLibmesh<T_S,T_P>::getMesh( )
    {
      return m_mesh;
    }


}

#endif // PHYSICS_LIBMESH_H
