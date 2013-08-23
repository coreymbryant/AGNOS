#ifndef PHYSICS_LIBMESH_H
#define PHYSICS_LIBMESH_H

#include "agnosDefines.h"
#include "PhysicsModel.h"

// libmesh includes
#include "libmesh/fem_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/adjoint_refinement_estimator.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/error_vector.h"
        

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
      /** TODO initialize variables specific to this physics in constructor */

      /** Constructor. Pass input file to provide setting to physics class */
      PhysicsLibmesh( const Communicator& comm_in, const GetPot& input );

      /** Destructor */
      virtual ~PhysicsLibmesh( );

      /** Function called by SurrogateModel to solve for requested solution
       * vectors at each evaluation point in parameter space  */
      void compute( 
          const T_S& paramVector, 
          std::map<std::string, T_P >& solutionVectors 
          ) ;

      /** Refinement methods for the physics model */
      void refine (  ) ;
      void refine ( T_P& errorIndicators ) ;

      
      /** Return libMesh mesh object */
      const libMesh::Mesh getMesh( ) const { return *_mesh; }

    protected:
      /** communicator reference */
      const Communicator&   _communicator;
      /** reference to input file just incase its needed after initialization */
      const GetPot&         _input;

      /** mesh and equation pointers */
      libMesh::FEMSystem*                      _system;
      libMesh::EquationSystems*             _equationSystems;
      libMesh::Mesh*                        _mesh; // mesh
      libMesh::MeshRefinement*              _meshRefinement;
      libMesh::AdjointRefinementEstimator*  _errorEstimator;
      libMesh::QoISet*                      _qois;

      // problem specific routines
      // TODO add residual and jacobian?
      /* PhysicsQoi<T_S>*                _qoi; */
      /* PhysicsQoiDerivative<T_S>*      _qoiDerivative; */
      
      /** refinement options */
      bool         _useUniformRefinement;
      bool         _resolveAdjoint;
      unsigned int _numberHRefinements;
      unsigned int _numberPRefinements;
      unsigned int _maxRefineSteps;

      /** derived PhysicsModel classes need to handle settig parameter values
       * themselves */
      virtual void _setParameterValues( const T_S& parameterValues ) = 0;

      /** build mesh refinement object (can be overidden in derived class) */
      virtual void _buildMeshRefinement();

      /** build error estimator object. Defaults to adjoint_residual_estimator
       * (can be overridden in derived class) */
      virtual void _buildErrorEstimator();


  };




/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsLibmesh<T_S,T_P>::PhysicsLibmesh(
      const Communicator& comm_in,
      const GetPot&           input
      ) : 
    PhysicsModel<T_S,T_P>(comm_in,input),
    _communicator(comm_in),_input(input)
  {


    if(AGNOS_DEBUG)
      std::cout << "DEBUG: libMesh initialized?:" << libMesh::initialized() << std::endl;


    // -----------------------------------------------------------
    //              get physics model settings
    // -----------------------------------------------------------
    std::cout << "\n-- Reading Physics model data\n";
    if(AGNOS_DEBUG)
    {
      int worldRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&worldRank);
      std::cout << "DEBUG: worldRank:" << worldRank << std::endl;
    }

    // read refinement options
    _useUniformRefinement = input("physics/useUniformRefinement",true);
    _numberHRefinements   = input("physics/numberHRefinements",1);
    _numberPRefinements   = input("physics/numberPRefinements",0);
    _maxRefineSteps       = input("physics/maxRefineSteps",1);

    // other options
    _resolveAdjoint       = input("physics/resolveAdjoint",false);
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
    /* delete _system; */
    delete _mesh;
    delete _meshRefinement;
    delete _errorEstimator;
    delete _qois;
    delete _equationSystems;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    void PhysicsLibmesh<T_S,T_P>::compute(
        const T_S& paramVector,
        std::map<std::string, T_P >& solutionVectors 
        ) 
    {
      if (AGNOS_DEBUG)
        std::cout << "DEBUG:  begin: PhysicsModel::compute(...)" << std::endl;
      
      // clear out any old data
      solutionVectors.clear();


      // An EquationSystems reference will be convenient.
      FEMSystem& system = *_system;
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

      if(AGNOS_DEBUG)
        std::cout << "DEBUG: pre solve: " 
          << " n_vars: " << system.n_vars()
          << " n_dofs: " << system.n_dofs()
          << std::endl;

      std::cout << " number of systems: " << es.n_systems() << std::endl;
      std::cout << " number of vars: " << system.n_vars() << std::endl;
      std::cout << " number of comp: " << system.n_components() << std::endl;
      libMesh::MeshBase::element_iterator elemIt = mesh.elements_begin();
      for(;elemIt!=mesh.elements_end();elemIt++)
      {
        Elem* elem = *elemIt;
        std::cout << "elem.n_vars(s): " 
          << elem->n_vars(0) << std::endl;
        std::cout << "elem.n_dofs(s,var): " 
          << elem->n_dofs(0,0) << std::endl;
        std::cout << "elem.n_comp(s,vars): " 
          << elem->n_comp(0,0) << std::endl;
      }

      // solve system
      system.solve();

      // save solution in solutionVectors
      std::vector<Number> primalSolution;
      system.solution->localize(primalSolution);
      solutionVectors.insert( 
          std::pair<std::string,T_P>( "primal", T_P(primalSolution) )
            );
      
      if (AGNOS_DEBUG)
      {
        std::cout << "DEBUG: primal solution\n:" ;
        system.solution->print_global();
      }

      // TODO if adjoint requested
      // TODO if error requested
      // TODO if errorIndicators requested
      //
      
      system.adjoint_solve();
      system.set_adjoint_already_solved(true);

      if(AGNOS_DEBUG)
      {
        std::cout << "DEBUG: adjoint solution\n:" ;
        system.get_adjoint_solution(0).print_global();
      }


      
      // We will declare an error vector for passing to the adjoint refinement error estimator
      this->_buildErrorEstimator();
      ErrorVector QoI_elementwise_error;
      
      // Estimate the error in each element using the Adjoint Refinement estimator
      _errorEstimator->estimate_error(system, QoI_elementwise_error);

      if(AGNOS_DEBUG)
        std::cout << "DEBUG: post estimate_error" << std::endl;

      // save adjoint in solutionVectors
      std::vector<Number> adjointSolution ;
      system.get_adjoint_solution(0).localize(adjointSolution);
      if(AGNOS_DEBUG)
        std::cout << "DEBUG: post estimate_error:adjoint solution size" << adjointSolution.size() << std::endl;
      
      solutionVectors.insert( 
          std::pair<std::string,T_P >( "adjoint", T_P(adjointSolution) )
            );

      // save errorEstimate in solutionVectors
      T_P errorEstimate;
      errorEstimate.resize(1);
      errorEstimate(0) = _errorEstimator->get_global_QoI_error_estimate(0);
      solutionVectors.insert( 
          std::pair<std::string,T_P >(
            "errorEstimate", errorEstimate)
            );
      if(AGNOS_DEBUG)
        std::cout << "DEBUG: error[0]:" << errorEstimate(0) << std::endl;

      // save errorIndicators in solutionVectors
      T_P errorIndicators( QoI_elementwise_error.size() );
      for (unsigned int i=0; i<errorIndicators.size(); i++)
      {
        if(AGNOS_DEBUG)
          std::cout << "DEBUG: error_per_cell[" << i << "]:" <<
            QoI_elementwise_error[i] << std::endl;
        errorIndicators(i) = std::abs(QoI_elementwise_error[i]);
      }
      solutionVectors.insert( 
          std::pair<std::string,T_P >(
            "errorIndicators", errorIndicators)
            );

      // save QoI value to solutionVectors
      system.assemble_qoi();
      T_P qoiValue;
      qoiValue.resize(1);
      qoiValue(0) = system.qoi[0] ;
      solutionVectors.insert( 
          std::pair<std::string,T_P >(
            "qoi", qoiValue)
            );
      if(AGNOS_DEBUG)
        std::cout << "DEBUG: qoi[0]:" << qoiValue(0) << std::endl;

      

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
      if ( this->_communicator.rank() == 0 )
        std::cout << "  previous n_active_elem(): " 
          << _mesh->n_active_elem() << std::endl;

      _meshRefinement->uniformly_refine(1);

      _equationSystems->reinit();

      if ( this->_communicator.rank() == 0 )
        std::cout << "   refined n_active_elem(): " 
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
      /* _meshRefinement->flag_elements_by_error_fraction (error_per_cell); */             
      _meshRefinement->flag_elements_by_elem_fraction (error_per_cell);              
                   
      //TODO input options to control whether we coarsen as well
      _meshRefinement->refine_and_coarsen_elements();
      /* _meshRefinement->refine_elements(); */


      _equationSystems->reinit();

      if ( this->_communicator.rank() == 0 )
        std::cout << "   refined n_active_elem(): " 
          << _mesh->n_active_elem() << std::endl;

      return;
    }


}

#endif // PHYSICS_LIBMESH_H
