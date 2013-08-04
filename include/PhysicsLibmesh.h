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
      PhysicsLibmesh( const Parallel::Communicator& comm_in, const GetPot& input );
      ~PhysicsLibmesh( );

      void compute( 
          const T_S& paramVector, 
          std::map<std::string, T_P > solutionVectors 
          ) ;

      void refine (  ) ;
      void refine ( T_P& errorIndicators ) ;

      libMesh::Mesh getMesh( ) const { return _mesh; }

    protected:
      const GetPot&  _input;

      // mesh and equation variables
      libMesh::Mesh*            _mesh; // mesh
      libMesh::EquationSystems* _equationSystems;
      libMesh::System*          _system;
      libMesh::MeshRefinement*  _meshRefinement;
      libMesh::AdjointRefinementEstimator*  _errorEstimator;

      // problem specific routines
      // TODO add residual and jacobian?
      PhysicsQoi<T_S>*                _qoi;
      PhysicsQoiDerivative<T_S>*      _qoiDerivative;
      libMesh::QoISet*                _qois;
      
      // refinement options
      bool         _useUniformRefinement;
      bool         _resolveAdjoint;
      unsigned int _numberHRefinements;
      unsigned int _numberPRefinements;
      unsigned int _maxRefineSteps;


      // derived PhysicsModel classes need to handle setting parameter values
      virtual void _setParameterValues( const T_S& parameterValues ) = 0;

      // build mesh refinement object (can be overidden in derived class)
      virtual void _build_mesh_refinement();

      // build error estimator object. Defaults to adjoint_residual_estimator 
      // (can be overidden in derived class)
      virtual void _build_error_estimator();
  };




/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsLibmesh<T_S,T_P>::PhysicsLibmesh(
      const Parallel::Communicator& comm_in,
      const GetPot&           input
      ) : 
    PhysicsModel<T_S,T_P>(comm_in),
    _input(input)
  {

    if(DEBUG)
      std::cout << "libMesh initialized?:" << libMesh::initialized() << std::endl;


    // -----------------------------------------------------------
    //              get physics model settings
    // -----------------------------------------------------------
    std::cout << "\n-- Reading Physics model data\n";
    if(DEBUG)
    {
      int worldRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&worldRank);
      std::cout << "worldRank:" << worldRank << std::endl;
    }

    // read refinement options
    _useUniformRefinement = input("physics/useUniformRefinement",true);
    _numberHRefinements   = input("physics/numberHRefinements",0);
    _numberPRefinements   = input("physics/numberPRefinements",1);
    _maxRefineSteps       = input("physics/maxRefineSteps",1);

    // other options
    _resolveAdjoint       = input("physics/resolveAdjoint",false);
    // -----------------------------------------------------------


    // -----------------------------------------------------------
    // which solutions do we want to compute a surrogate for
    this->_solutionNames.clear( );

    for(unsigned int i=0; i<input.vector_variable_size("physics/solutions") ; i++)
      this->_solutionNames.insert( input("physics/solutions"," ",i) ) ;

    // default to only primal solution
    if(this->_solutionNames.size() == 0)
      this->_solutionNames.insert("primal");

    if(DEBUG)
      std::cout << "solutionNames.size()" << this->_solutionNames.size() << std::endl;
    // -----------------------------------------------------------
    



  }
  
  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S, class T_P>
  void PhysicsLibmesh<T_S,T_P>::_build_mesh_refinement( )
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
  void PhysicsLibmesh<T_S,T_P>::_build_error_estimator( )
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
    /* delete _mesh; */
    /* delete _equationSystems; */
    /* delete _system; */
    /* delete _meshRefinement; */

    /* delete _qoi; */
    /* delete _qoiDerivative; */
    /* delete _qois; */
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    void PhysicsLibmesh<T_S,T_P>::compute(
        const T_S& paramVector,
        std::map<std::string, T_P > solutionVectors 
        ) 
    {
      if (DEBUG)
        std::cout << " begin: PhysicsModel::compute(...)" << std::endl;
      
      // clear out any old data
      solutionVectors.clear();


      // An EquationSystems reference will be convenient.
      System& system = const_cast<System&>(*_system);
      EquationSystems& es = system.get_equation_systems();

      // The current mesh
      MeshBase& mesh = es.get_mesh();
      if(DEBUG)
        mesh.print_info();

      // set parameter values
      this->_setParameterValues( paramVector );
      // reinitialize system
      es.reinit();

      // solve system
      system.solve();

      // save solution in solutionVectors
      std::vector<Number> primalSolution;
      system.solution->localize(primalSolution);
      solutionVectors.insert( 
          std::pair<std::string,T_P>( "primal", T_P(primalSolution) )
            );
      
      if (DEBUG)
        system.solution->print_global();

      // TODO if adjoint requested
      // TODO if error requested
      // TODO if errorIndicators requested
      
      // We will declare an error vector for passing to the adjoint refinement error estimator
      ErrorVector QoI_elementwise_error;
      
      // Estimate the error in each element using the Adjoint Refinement estimator
      _errorEstimator->estimate_error(system, QoI_elementwise_error);

      // save adjoint in solutionVectors
      std::vector<Number> adjointSolution ;
      system.get_adjoint_solution(0).localize(adjointSolution);
      solutionVectors.insert( 
          std::pair<std::string,T_P >( "adjoint", T_P(adjointSolution) )
            );
      if(DEBUG)
        system.get_adjoint_solution(0).print_global();

      // save errorEstimate in solutionVectors
      T_P errorEstimate;
      errorEstimate.resize(1);
      errorEstimate(0) = _errorEstimator->get_global_QoI_error_estimate(0);
      solutionVectors.insert( 
          std::pair<std::string,T_P >(
            "errorEstimate", errorEstimate)
            );
      if(DEBUG)
        std::cout << "error[0]:" << errorEstimate(0) << std::endl;

      // save errorIndicators in solutionVectors
      T_P errorIndicators( QoI_elementwise_error.size() );
      for (unsigned int i=0; i<errorIndicators.size(); i++)
      {
        if(DEBUG)
          std::cout << "error_per_cell[" << i << "]:" <<
            QoI_elementwise_error[i] << std::endl;
        errorIndicators(i) = QoI_elementwise_error[i];
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
            "qoiValue", qoiValue)
            );
      if(DEBUG)
        std::cout << "qoi[0]:" << qoiValue(0) << std::endl;

      

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
