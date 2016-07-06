
#include "PhysicsLibmesh.h"

// libmesh includes
#include "libmesh/dof_map.h"
#include "libmesh/fem_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/mesh.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/adjoint_refinement_estimator.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/error_vector.h"

#include LIBMESH_INCLUDE_UNORDERED_MAP
#include LIBMESH_INCLUDE_UNORDERED_SET

namespace AGNOS
{

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsLibmesh<T_S,T_P>::PhysicsLibmesh(
      const Communicator& comm_in,
      const GetPot&           input
      ) : 
    _system(NULL),
    _equationSystems(NULL),
    _mesh(NULL), // mesh
    _meshRefinement(NULL),
    _errorEstimator(NULL),
    _qois(NULL),
    _timeStep(1),
    _primalExio(NULL),
    _adjointExio(NULL),
    PhysicsModel<T_S,T_P>(comm_in,input)
  {

    // default to only primal solution
    if(this->_availableSolutions.size() == 0)
      this->_availableSolutions.insert("primal");

    if(AGNOS_DEBUG)
      std::cout << "DEBUG: libMesh initialized?:" << libMesh::initialized() 
        << std::endl;


    // -----------------------------------------------------------
    //              get physics model settings
    // -----------------------------------------------------------
    std::cout << "\n-- Reading Physics model data\n";
    if(AGNOS_DEBUG)
    {
      int worldRank;
      MPI_Comm_rank(this->_communicator.get(),&worldRank);
      std::cout << "DEBUG: worldRank:" << worldRank << std::endl;
    }

    // read refinement options
    _useUniformRefinement = input("useUniformRefinement",true);
    _numberHRefinements   = input("numberHRefinements",1);
    _numberPRefinements   = input("numberPRefinements",0);
    _maxRefineSteps       = input("maxRefineSteps",1);

    // other options
    _writePrimalViz        = input("writePrimalViz",false);
    _writeAdjointViz       = input("writeAdjointViz",false);
    _resolveAdjoint       = input("resolveAdjoint",false);
    // -----------------------------------------------------------



  }
  
  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S, class T_P>
  void PhysicsLibmesh<T_S,T_P>::_buildMeshRefinement( )
  {
    // set up mesh refinement object
    _meshRefinement = new libMesh::MeshRefinement(*_mesh);
    
    // TODO read these from input file?
    _meshRefinement->coarsen_by_parents()         = true;
    _meshRefinement->absolute_global_tolerance()  = 1e-6;
    /* _meshRefinement->nele_target()               = 64; */  
    _meshRefinement->refine_fraction()            = 0.7;
    _meshRefinement->coarsen_fraction()           = 0.3;
    _meshRefinement->coarsen_threshold()          = 1e-5;
    _meshRefinement->max_h_level()                = 15;
  }

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S, class T_P>
  void PhysicsLibmesh<T_S,T_P>::_buildErrorEstimator( )
  {
    // set up error estimator object
    _errorEstimator = new AdjointRefinementEstimator; //TODO make this an option
    _errorEstimator->qoi_set() = *_qois; 
    _errorEstimator->number_h_refinements = _numberHRefinements;
    _errorEstimator->number_p_refinements = _numberPRefinements;

  }
    

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsLibmesh<T_S,T_P>::~PhysicsLibmesh( )
  {
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    void PhysicsLibmesh<T_S,T_P>::_setPrimalSolution(
        T_P& solutionVector ) 
    {
      // set primal solution with value from solutionVectors
      NumericVector<Number>& solution = *(_system->solution) ;
      solution.close();

      // make sure sizes agree
      agnos_assert( (solution.size() == solutionVector.size())) ;

      for (unsigned int i=0; i<solution.size(); i++)
        solution.set(i, solutionVector(i) ) ;
      solution.close();

      this->_system->update();

      return;
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    void PhysicsLibmesh<T_S,T_P>::_setAdjointSolution( 
        T_P& solutionVector, unsigned int j )
    {
      // set adjoint solution from solutionVectors
      NumericVector<Number>& adjoint_solution =
        _system->get_adjoint_solution(j) ;

      // make sure sizes agree
      agnos_assert( (adjoint_solution.size() == solutionVector.size())) ;

      for (unsigned int i=0; i<adjoint_solution.size(); i++)
        adjoint_solution.set(i, solutionVector(i)) ;
      adjoint_solution.close();

      return ;
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    void PhysicsLibmesh<T_S,T_P>::_insertSolVec( 
        libMesh::NumericVector<Number>& vec,
        std::string vecName,
        std::map<std::string, T_P >& mapOfVecs 
        ) 
      {
        std::vector<Number> localVec;
        vec.localize(localVec);
        mapOfVecs.insert( 
            std::pair<std::string,T_P>( vecName, T_P(localVec) )
              );

        return;
      }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    void PhysicsLibmesh<T_S,T_P>::_solve( ) 
    {
      // solve system
      _system->solve();
      return;
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    void PhysicsLibmesh<T_S,T_P>::_adjointSolve(
        libMesh::QoISet& qois ) 
    {
      // solve adjoint equation
      _system->adjoint_solve( qois );
      return;
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    void PhysicsLibmesh<T_S,T_P>::compute(
        std::set<std::string>& computeSolutions,
        const T_S& paramVector,
        std::map<std::string, T_P >& solutionVectors 
        ) 
    {
      if (AGNOS_DEBUG)
      {
        std::cout << "DEBUG:  begin: PhysicsModel::compute(...)" << std::endl;
        std::cout << "DEBUG:         rank: " << this->_communicator.rank() 
          << std::endl;
      }


      //TODO check that every computSolution requested is present for this model
      


      // An EquationSystems reference will be convenient.
      System& system = *_system;
      EquationSystems& es = system.get_equation_systems();

      // The current mesh
      MeshBase& mesh = es.get_mesh();

      
      // set parameter values
      this->_setParameterValues( paramVector );
      // reinitialize system
      es.reinit();

      if(AGNOS_DEBUG)
      {
        std::cout << "DEBUG: in compute routine:" << std::endl;
        mesh.print_info();
        es.print_info();
      }
      
      // broadcast solutionVectors to all physics procs
      int groupRank;
      MPI_Comm_rank(this->_communicator.get(),&groupRank);
      typename std::map<std::string,T_P>::iterator cid ;
      if(groupRank==0)
      {
        int nVectors = solutionVectors.size();
        this->_communicator.broadcast( nVectors );
        /* std::cout << "send nVectors: " << nVectors << std::endl; */
        for(cid=solutionVectors.begin();cid!=solutionVectors.end();cid++)
        {
          std::string solName = cid->first;
          this->_communicator.broadcast(solName);
          int solSize = cid->second.size();
          this->_communicator.broadcast(solSize);
          std::vector<double> dbuf(solSize,0.);
          /* std::cout << "solSize: " << solSize << std::endl; */
          for(unsigned int i=0; i<cid->second.size();i++)
          {
            dbuf[i]=cid->second(i) ;
            /* std::cout << "dbuf[" << i << "] = " << dbuf[i] << std::endl; */
          }
          this->_communicator.broadcast(dbuf);
        }
      }
      else
      {
        solutionVectors.clear();
        int nVectors;
        this->_communicator.broadcast( nVectors );
        /* std::cout << "recv nVectors: " << nVectors << std::endl; */
        for(unsigned int i=0; i<nVectors; i++)
        {
          std::string sbuf;
          this->_communicator.broadcast(sbuf);
          int solSize;
          this->_communicator.broadcast(solSize);
          std::vector<double> dbuf(solSize,0.);
          this->_communicator.broadcast(dbuf);

          /* std::cout << "solSize: " << solSize << std::endl; */
          for(unsigned int i=0; i<dbuf.size();i++)
          {
            /* std::cout << "dbuf[" << i << "] = " << dbuf[i] << std::endl; */
          }

          solutionVectors.insert( std::pair<std::string,T_P>(
                sbuf, T_P(dbuf) ) ) ;
        }
      }
      

      // PRIMAL SOLUTION
      // Primal solution must always be computed unless it was provided
      if (!solutionVectors.count("primal"))
      {
        
        // solve system
        _solve( ) ;

        //save the vizualizaton if requested
        if (_writePrimalViz)
          _outputPrimalViz( paramVector ) ;
        
        // save solution in solutionVectors
        _insertSolVec( *(system.solution), "primal", solutionVectors );
        
        if (AGNOS_DEBUG)
        {
          int globalRank;
          MPI_Comm_rank(MPI_COMM_WORLD,&globalRank);
          std::cout << "DEBUG: primal solution:\n " ;
          std::cout << "DEBUG: rank " << globalRank << std::endl;
          std::cout << "DEBUG: point " << paramVector(0) << std::endl;
          system.solution->print_global();
        }

      }


      

      // ADJOINT SOLUTION AND ERROR ESTIMATE/INDICATORS
      // copied from libmesh adjoint_error_refinement_estimator estimate_errors
      // routine
      /** if adjoint is requested, and not provided, or resolve flag is set */
      if( computeSolutions.count("adjoint" ) 
          || computeSolutions.count("errorEstimate") 
          || computeSolutions.count("errorIndicators") 
          )
      {
        if(AGNOS_DEBUG)
          std::cout << "Refining for adjoint|error solve" << std::endl;

        // set primal solution with value from solutionVectors
        _setPrimalSolution( solutionVectors["primal"] ) ;


        // Initialize and resize the error_per_cell in case we need it later
        // wont have access to original number of elements
        std::vector<float> error_per_cell;
        error_per_cell.resize (mesh.max_elem_id(), 0.);
        
        // We'll want to back up all coarse grid vectors
        std::map<std::string, NumericVector<Number> *> coarse_vectors;
        for (System::vectors_iterator vec = system.vectors_begin(); vec !=
             system.vectors_end(); ++vec)
          {
            // The (string) name of this vector
            const std::string& var_name = vec->first;

            coarse_vectors[var_name] = vec->second->clone().release();
          }
        // Back up the coarse solution and coarse local solution
        NumericVector<Number> * coarse_solution =
          system.solution->clone().release();
        NumericVector<Number> * coarse_local_solution =
          system.current_local_solution->clone().release();
        // And make copies of the projected solution
        NumericVector<Number> * projected_solution;

        // And we'll need to temporarily change solution projection settings
        bool old_projection_setting;
        old_projection_setting = system.project_solution_on_reinit();

        // Make sure the solution is projected when we refine the mesh
        system.project_solution_on_reinit() = true;

        // And it'll be best to avoid any repartitioning
        UniquePtr<Partitioner> old_partitioner(mesh.partitioner().release());

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

        libmesh_assert (_numberHRefinements > 0 || _numberPRefinements > 0);

        for (unsigned int i = 0; i != _numberHRefinements; ++i)
          {
            mesh_refinement.uniformly_refine(1);
            es.reinit();
          }

        for (unsigned int i = 0; i != _numberPRefinements; ++i)
          {
            mesh_refinement.uniformly_p_refine(1);
            es.reinit();
          }

        // Copy the projected coarse grid solutions, which will be
        // overwritten by solve()
        projected_solution = NumericVector<Number>::build(mesh.comm()).release();
        projected_solution->init(system.solution->size(), true, SERIAL);
        system.solution->localize(*projected_solution,
                system.get_dof_map().get_send_list());


        // Rebuild the rhs with the projected primal solution
        (dynamic_cast<ImplicitSystem&>(system)).assembly(true, false);
        NumericVector<Number> & projected_residual 
          = (dynamic_cast<ExplicitSystem&>(system)).get_vector("RHS Vector");
        projected_residual.close();
        /* projected_residual.print_global(); */

        if ( solutionVectors.count("adjoint") && !_resolveAdjoint )
        {
          // Make sure an adjoint solution exists
          // if solving in parallel we may not have an adjoint solution in this
          // system reference yet
          bool has_adjoint = system.have_vector("adjoint_solution0");
          if (!has_adjoint)
            system.add_adjoint_solution();
          
          // set adjoint solution from solutionVectors
          _setAdjointSolution( solutionVectors["adjoint"], 0 );
        }
        else // solve adjoint
        {
          if(AGNOS_DEBUG)
            std::cout << "Solving adjoint solution" << std::endl;
          _adjointSolve( *_qois ) ;
          
          //save the vizualizaton if requested
          if (_writeAdjointViz)
            _outputAdjointViz( paramVector );
          
          // save adjoint solution if its in the set of requested vectors
          if ( computeSolutions.count("adjoint") )
          {
            // save adjoint in solutionVectors
            _insertSolVec( system.get_adjoint_solution(0), "adjoint", solutionVectors );
          }
        }
        
        

        if(AGNOS_DEBUG)
        {
          std::cout << "DEBUG: adjoint solution\n:" ;
          system.get_adjoint_solution(0).print_global();
        }


        // ----------------------------------
        // ERROR ESTIMATE AND ERROR INDICATORS
        // Compute global error estimate if requested
        if ( computeSolutions.count("errorEstimate") ) 
        {
          if(AGNOS_DEBUG)
            std::cout << "DEBUG: Computing errorEstimate" << std::endl;
          // Now that we have the refined adjoint solution and the projected
          // primal solution, we first compute the global QoI error estimate

          // Resize the computed_global_QoI_errors vector to hold the error
          // estimates for each QoI
          T_P errorEstimate(system.qoi.size());

          // Loop over all the adjoint solutions and get the QoI error
          // contributions from all of them
          for (unsigned int j=0; j != system.qoi.size(); j++) 
          {
              errorEstimate(j) = projected_residual.dot(system.get_adjoint_solution(j));
          }

          // save error estimate in solutionVectors
          _insertSolVec( errorEstimate, "errorEstimate", solutionVectors );

          if(AGNOS_DEBUG)
            std::cout << "DEBUG: error[0]:" << errorEstimate(0) << std::endl;
        } // end if errorEstimate


        // compute indicators if they were requested
        if (  computeSolutions.count("errorIndicators") )
        {
          if(AGNOS_DEBUG)
            std::cout << "DEBUG: computing error indicators" << std::endl;

          // A map that relates a node id to an int that will tell us how many
          // elements it is a node of
          LIBMESH_BEST_UNORDERED_MAP<dof_id_type, unsigned int>shared_element_count;

          // To fill this map, we will loop over elements, and then in each
          // element, we will loop over the nodes each element contains, and
          // then query it for the number of coarse grid elements it was a node
          // of

          // We will be iterating over all the active elements in the fine mesh
          // that live on this processor
          MeshBase::const_element_iterator elem_it = mesh.active_local_elements_begin();
          const MeshBase::const_element_iterator elem_end 
            = mesh.active_local_elements_end();

          // Keep track of which nodes we have already dealt with
          LIBMESH_BEST_UNORDERED_SET<dof_id_type> processed_node_ids;

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

              // If we havent already processed this node, do so now
              if(processed_node_ids.find(node_id) == processed_node_ids.end())
              {
                // Declare a neighbor_set to be filled by the find_point_neighbors
                std::set<const Elem *> fine_grid_neighbor_set;

                // Call find_point_neighbors to fill the neighbor_set
                elem->find_point_neighbors(*node, fine_grid_neighbor_set);

                // A vector to hold the coarse grid parents neighbors
                std::vector<dof_id_type> coarse_grid_neighbors;

                // Iterators over the fine grid neighbors set
                std::set<const Elem*>::iterator fine_neighbor_it 
                  = fine_grid_neighbor_set.begin();
                const std::set<const Elem*>::iterator fine_neighbor_end 
                  = fine_grid_neighbor_set.end();

                // Loop over all the fine neighbors of this node
                for(; fine_neighbor_it != fine_neighbor_end ;
                    ++fine_neighbor_it)
                {
                  // Pointer to the current fine neighbor element
                  const Elem* fine_elem = *fine_neighbor_it;

                  // Find the element id for the corresponding coarse grid
                  // element
                  const Elem* coarse_elem = fine_elem;
                  for (unsigned int j = 0; j != _numberHRefinements; ++j)
                  {
                    libmesh_assert (coarse_elem->parent());

                    coarse_elem = coarse_elem->parent();
                  }

                  // Loop over the existing coarse neighbors and check if this
                  // one is already in there
                  const dof_id_type coarse_id = coarse_elem->id();
                  std::size_t j = 0;
                  for (; j != coarse_grid_neighbors.size(); j++)
                  {
                    // If the set already contains this element break out of the
                    // loop
                    if(coarse_grid_neighbors[j] == coarse_id)
                    {
                      break;
                    }
                  }

                  // If we didn't leave the loop even at the last element, this
                  // is a new neighbour, put in the coarse_grid_neighbor_set
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
                processed_node_ids.insert(node_id);

              } // End if not processed node

            } // End loop over nodes

          }  // End loop over elems

          // Get a DoF map, will be used to get the nodal dof_indices for each
          // element
          DofMap &dof_map = system.get_dof_map();

          // The global DOF indices, we will use these later on when we compute
          // the element wise indicators
          std::vector<dof_id_type> dof_indices;

          // Localize the global rhs and adjoint solution vectors (which might
          // be shared on multiple processsors) onto a local ghosted vector,
          // this ensures each processor has all the dof_indices to compute an
          // error indicator for an element it owns
          UniquePtr<NumericVector<Number> > localized_projected_residual 
            = NumericVector<Number>::build(system.comm());
          localized_projected_residual->init(system.n_dofs(),
              system.n_local_dofs(), system.get_dof_map().get_send_list(),
              false, GHOSTED);
          projected_residual.localize(*localized_projected_residual,
              system.get_dof_map().get_send_list());

          // Each adjoint solution will also require ghosting; for efficiency
          // we'll reuse the same memory
          UniquePtr<NumericVector<Number> > localized_adjoint_solution =
            NumericVector<Number>::build(system.comm());
          localized_adjoint_solution->init(system.n_dofs(),
              system.n_local_dofs(), system.get_dof_map().get_send_list(),
              false, GHOSTED);

          // We will loop over each adjoint solution, localize that adjoint
          // solution and then loop over local elements
          for (unsigned int i=0; i != system.qoi.size(); i++)
          {
            // Skip this QoI if not in the QoI Set
            if (_qois->has_index(i))
            {
              // Get the weight for the current QoI
              Real error_weight = _qois->weight(i);

              (system.get_adjoint_solution(i)).localize(*localized_adjoint_solution,
                  system.get_dof_map().get_send_list());

              // Loop over elements
              MeshBase::const_element_iterator elem_it =
                mesh.active_local_elements_begin();
              const MeshBase::const_element_iterator elem_end =
                mesh.active_local_elements_end();

              for(; elem_it != elem_end; ++elem_it)
              {
                // Pointer to the element
                const Elem* elem = *elem_it;

                // Go up number_h_refinements levels up to find the coarse parent
                const Elem* coarse = elem;

                for (unsigned int j = 0; j != _numberHRefinements; ++j)
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
                  // The contribution to the error indicator for this element
                  // from the current QoI
                  local_contribution +=
                    (*localized_projected_residual)(dof_indices[j]) *
                    (*localized_adjoint_solution)(dof_indices[j]);
                }

                // Multiply by the error weight for this QoI
                local_contribution *= error_weight;

                // FIXME: we're throwing away information in the
                // --enable-complex case
                error_per_cell[e_id] += static_cast<ErrorVectorReal>
                  (libmesh_real(local_contribution));

              } // End loop over elements

            } // End if belong to QoI set

          } // End loop over QoIs
          
          // Fiinally sum the vector of estimated error values.
          /* this->reduce_error(error_per_cell, system.comm()); */
          this->_communicator.sum( error_per_cell );

          
          // save errorIndicators in solutionVectors
          T_P errorIndicators( error_per_cell.size() );
          for (unsigned int i=0; i<errorIndicators.size(); i++)
          {
            errorIndicators(i) = std::abs(error_per_cell[i]);
            if(AGNOS_DEBUG)
              std::cout << "DEBUG: error_per_cell[" << i << "]:" <<
                error_per_cell[i] << std::endl;
          }

          // save errorIndicators in solutionVectors
          _insertSolVec( errorIndicators, "errorIndicators", solutionVectors );

        } // end if error indicators

        // Don't bother projecting the solution; we'll restore from backup
        // after coarsening
        system.project_solution_on_reinit() = false;

        // Uniformly coarsen the mesh, without projecting the solution
        libmesh_assert (_numberHRefinements > 0 || _numberPRefinements > 0);

        for (unsigned int i = 0; i != _numberHRefinements; ++i)
          {
            mesh_refinement.uniformly_coarsen(1);
            // FIXME - should the reinits here be necessary? - RHS
            es.reinit();
          }

        for (unsigned int i = 0; i != _numberPRefinements; ++i)
          {
            mesh_refinement.uniformly_p_coarsen(1);
            es.reinit();
          }

        // We should be back where we started
        libmesh_assert_equal_to (n_coarse_elem, mesh.n_elem());

        // Restore old solutions and clean up the heap
        system.project_solution_on_reinit() = old_projection_setting;

        // Restore the coarse solution vectors and delete their copies
        *system.solution = *coarse_solution;
        delete coarse_solution;
        *system.current_local_solution = *coarse_local_solution;
        delete coarse_local_solution;
        delete projected_solution;

        for (System::vectors_iterator vec = system.vectors_begin(); vec !=
         system.vectors_end(); ++vec)
          {
            // The (string) name of this vector
            const std::string& var_name = vec->first;

            // If it's a vector we already had (and not a newly created
            // vector like an adjoint rhs), we need to restore it.
            std::map<std::string, NumericVector<Number> *>::iterator it =
              coarse_vectors.find(var_name);
            if (it != coarse_vectors.end())
              {
                NumericVector<Number> *coarsevec = it->second;
                system.get_vector(var_name) = *coarsevec;

                coarsevec->clear();
                delete coarsevec;
              }
          }

        // Restore old partitioner and renumbering settings
        mesh.partitioner().reset(old_partitioner.release());
        mesh.allow_renumbering(old_renumbering_setting);

        
      } // end adjoint|errorEstimate|errorIndicators
        

      

      // QUANTITY OF INTEREST
      if ( computeSolutions.count("qoi") )
      {
        if(AGNOS_DEBUG)
          std::cout << "DEBUG computing qoi value" << std::endl;

        // set primal solution with value from solutionVectors
        NumericVector<Number>& solution = *(system.solution) ;
        solution.close();
        for (unsigned int i=0; i<solution.size(); i++)
          solution.set(i,solutionVectors["primal"](i)) ;
        solution.close();

        // Compute QoI with this value
        system.assemble_qoi();

        // save QoI value to solutionVectors
        T_P qoiValue(1);
        qoiValue(0) = _getQoIValue( );
        solutionVectors.insert( 
            std::pair<std::string,T_P >(
              "qoi", qoiValue)
              );
        if(AGNOS_DEBUG)
          std::cout << "DEBUG: qoi[0]:" << qoiValue(0) << std::endl;
      }

      // EXACT Qoi
      /** exact solution is optional, but it also depends on if it has been
       * defined in the derived class
       */
      if( computeSolutions.count("exactQoi" ) )
      {
        // save solution in solutionVectors
        solutionVectors.insert( 
            std::pair<std::string,T_P>( "exactQoi", this->exactQoi() )
              );
      }

      

      
      // incremenent dummy timestep
      _timeStep++;


      return;
    }


  /********************************************//**
   * \brief 
   *
   * 
   ***********************************************/
  template<class T_S, class T_P>
    void PhysicsLibmesh<T_S,T_P>::refine() 
    {
        std::cout << "physics refining" << std::endl;
      /* if ( this->_communicator.rank() == 0 ) */
        std::cout << "  previous n_active_elem(): " 
          << _mesh->n_active_elem() << std::endl;

      _meshRefinement->uniformly_refine(1);

      _equationSystems->reinit();

      /* if ( this->_communicator.rank() == 0 ) */
        std::cout << "un_refined n_active_elem(): " 
          << _mesh->n_active_elem() << std::endl;

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

      if ( this->_communicator.rank() == 0 )
        std::cout << "  previous n_active_elem(): " 
          << _mesh->n_active_elem() << std::endl;

      libMesh::ErrorVector error_per_cell(errorIndicators.size());
      for(unsigned int i=0; i<errorIndicators.size();i++)
        error_per_cell[i] = errorIndicators(i);
      _meshRefinement->clean_refinement_flags();

      
      //TODO input options to control type of refinement
      /* _meshRefinement->flag_elements_by_error_tolerance (error_per_cell); */            
      _meshRefinement->flag_elements_by_error_fraction (error_per_cell);             
      /* _meshRefinement->flag_elements_by_elem_fraction (error_per_cell); */              
                   
      //TODO input options to control whether we coarsen as well
      _meshRefinement->refine_and_coarsen_elements();
      /* _meshRefinement->refine_elements(); */


      _equationSystems->reinit();

      if ( this->_communicator.rank() == 0 )
        std::cout << "   refined n_active_elem(): " 
          << _mesh->n_active_elem() << std::endl;

      return;
    }

  /********************************************//**
   * \brief Output visulization of primal sol
   ***********************************************/
  template<class T_S, class T_P>
  void PhysicsLibmesh<T_S,T_P>::_outputPrimalViz( const T_S& paramVector )
    {
      if ( this->_mesh->mesh_dimension() == 1)
        _writeGnuplot(*this->_equationSystems,paramVector,"primal");
      else
        _writeExodus(*this->_equationSystems, paramVector, _primalExio,
            "primal.exo"); 
      return;
    }

  /********************************************//**
   * \brief Output visulization of adjoint sol
   ***********************************************/
  template<class T_S, class T_P>
  void PhysicsLibmesh<T_S,T_P>::_outputAdjointViz( const T_S& paramVector )
    {
      NumericVector<Number> &primal_solution = *this->_system->solution;
      NumericVector<Number> &dual_solution 
        = this->_system->get_adjoint_solution(0);
      primal_solution.swap(dual_solution);

      if ( this->_mesh->mesh_dimension() == 1)
        _writeGnuplot(*this->_equationSystems,paramVector,"adjoint");
      else
        _writeExodus( *this->_equationSystems, paramVector, _adjointExio,
            "adjoint.exo");

      primal_solution.swap(dual_solution);
      return ;
    }

  /********************************************//**
   * \brief Save libmesh solution as gnuplot format
   *
   * 
   ***********************************************/
  template<class T_S, class T_P>
  void PhysicsLibmesh<T_S,T_P>::_writeGnuplot(
      EquationSystems &es,
      const T_S& paramVector,       // The parameter vector
      std::string solution_type ) // primal or adjoint solve
  {
    MeshBase &mesh = es.get_mesh();

    std::ostringstream file_name_gp;
    /* file_name_gp << solution_type */
    /*               << ".out.gp." */
    /*               << std::setw(2) */
    /*               << std::setfill('0') */
    /*               << std::right */
    /*               << index; */
    file_name_gp << solution_type
                  << ".out.gp" ;
    for (unsigned int i=0; i<paramVector.size(); i++)
    {
      file_name_gp 
        << "_"
        << std::setiosflags(std::ios::scientific) 
        << std::setprecision(1) 
        << std::setw(3) 
        << paramVector(i) ;
    }

    GnuPlotIO(mesh).write_equation_systems
      (file_name_gp.str(), es);
  }

  /********************************************//**
   * \brief Output data in exodus format
   *
   * 
   ***********************************************/
  template<class T_S, class T_P>
  void PhysicsLibmesh<T_S,T_P>::_writeExodus(
      EquationSystems &es,
      const T_S& paramVector,       // The parameter vector
      libMesh::ExodusII_IO* exio,   // primal or adjoint solve
      std::string fileName) 
  {
    
    // output timestep info
    exio->write_timestep(
        fileName,
        *this->_equationSystems,
        this->_timeStep, this->_timeStep );

    // output solution vectors
    std::vector<std::string> solNames;
    std::vector<libMesh::Number> solVec;
    this->_equationSystems->build_variable_names(solNames,NULL,NULL);
    this->_equationSystems->build_solution_vector(solVec,NULL);
    exio->write_nodal_data(
        fileName,
        solVec,
        solNames
        );

    // output parameters as global vars
    std::vector<std::string> paramNames;
    for(unsigned int i=0; i<paramVector.size();i++)
      paramNames.push_back("parameter"+std::to_string(i));

    std::vector<Number> soln = paramVector.get_values();

    exio->write_global_data( soln, paramNames );

    return;
  }


  template class
    PhysicsLibmesh<libMesh::DenseVector<double>, libMesh::DenseVector<double> >;

}

