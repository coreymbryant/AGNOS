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
  template<class T_S>
  class PhysicsLibmesh : public PhysicsModel<T_S>
  {

    public:

      /** constructor */
      PhysicsLibmesh( Parallel::Communicator& comm_in, const GetPot& input );
      ~PhysicsLibmesh( );

      void compute( 
          const T_S& paramVector, 
          std::map<std::string, NumericVector<Number>* > solutionVectors 
          ) ;


      void refine (  ) ;
      void refine ( std::vector<Number>& errorIndicators ) ;
      

      libMesh::Mesh getMesh( ) const ;


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


      // virtual functions for constructing and solving model
      virtual void _init( ) = 0;
      virtual void _setParameterValues( const T_S& parameterValues ) = 0;

      /* virtual NumericVector<Number>* _solvePrimal( )   = 0; */
      /* virtual NumericVector<Number>* _solveAdjoint( )  = 0; */
      /* virtual NumericVector<Number>* _estimateError( ) = 0; */
      /* virtual NumericVector<Number>* _evaluateQoi( )   = 0; */

  };




/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S>
  PhysicsLibmesh<T_S>::PhysicsLibmesh(
      Parallel::Communicator& comm_in,
      const GetPot&           input
      ) : 
    PhysicsModel<T_S>(comm_in),
    _input(input)
  {


    if(DEBUG)
      std::cout << "initialized:" << libMesh::initialized() << std::endl;


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
    _numberHRefinements = input("physics/numberHRefinements",0);
    _numberPRefinements = input("physics/numberPRefinements",1);
    _maxRefineSteps  = input("physics/maxRefineSteps",1);

    // other options
    _resolveAdjoint = input("physics/resolveAdjoint",false);
    // -----------------------------------------------------------


    // -----------------------------------------------------------
    // which solutions do we want to compute a surrogate for
    this->_solutionNames.clear( );
    for(unsigned int i; i< input.vector_variable_size("physics/solutions"); i++)
      this->_solutionNames.insert( input("physics/solutions"," ",i) ) ;
    if(this->_solutionNames.size() == 0)
      this->_solutionNames.insert("primal");

    if(DEBUG)
      std::cout << "solutionNames.size()" << this->_solutionNames.size() << std::endl;
    // -----------------------------------------------------------

    // -----------------------------------------------------------
    // initialize model
    /* _init( ); */
    // -----------------------------------------------------------


    // -----------------------------------------------------------
    // set up mesh refinement object
    _meshRefinement = new MeshRefinement(*_mesh);
    // TODO read these from input file?
    _meshRefinement->coarsen_by_parents()         = true;
    _meshRefinement->absolute_global_tolerance()  = 1e-6;
    /* _meshRefinement->nele_target()               = 64; */  
    _meshRefinement->refine_fraction()            = 0.7;
    _meshRefinement->coarsen_fraction()           = 0.3;  
    _meshRefinement->coarsen_threshold()          = 1e-5;
    _meshRefinement->max_h_level()                = 15;
    // -----------------------------------------------------------
    

    // -----------------------------------------------------------
    // set up error estimator object
    _errorEstimator = new AdjointRefinementEstimator; //TODO make this an option
    _errorEstimator->qoi_set() = *_qois; 
    _errorEstimator->number_h_refinements = _numberHRefinements;
    _errorEstimator->number_p_refinements = _numberPRefinements;
    // -----------------------------------------------------------

  }


/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S>
  PhysicsLibmesh<T_S>::~PhysicsLibmesh( )
  {
    delete _mesh;
    delete _equationSystems;
    delete _system;
    delete _meshRefinement;

    delete _qoi;
    delete _qoiDerivative;
    delete _qois;
  }


/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S>
    void PhysicsLibmesh<T_S>::compute(
        const T_S& paramVector,
        std::map<std::string, NumericVector<Number>* > solutionVectors 
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
      _system->solve();

      // save solution in solutionVectors
      NumericVector<Number>* primalSolution = system.solution->clone().release();
      solutionVectors.insert( 
          std::pair<std::string,NumericVector<Number>* >(
            "primal", primalSolution)
            );
      
      if (DEBUG)
        _system->solution->print_global();

      // TODO if adjoint requested
      // TODO if error requested
      // TODO if errorIndicators requested
      
      // We will declare an error vector for passing to the adjoint refinement error estimator
      ErrorVector QoI_elementwise_error;
      
      // Estimate the error in each element using the Adjoint Refinement estimator
      _errorEstimator->estimate_error(system, QoI_elementwise_error);

      // save adjoint in solutionVectors
      NumericVector<Number>* adjointSolution =
        system.get_adjoint_solution(0).clone().release();
      solutionVectors.insert( 
          std::pair<std::string,NumericVector<Number>* >(
            "adjoint", adjointSolution)
            );

      // save errorEstimate in solutionVectors
      NumericVector<Number>* errorEstimate =
        NumericVector<Number>::build(_mesh->comm()).release();
      errorEstimate->init( 1, 1, false, SERIAL);
      errorEstimate->set(0, _errorEstimator->get_global_QoI_error_estimate(0) );
      solutionVectors.insert( 
          std::pair<std::string,NumericVector<Number>* >(
            "errorEstimate", errorEstimate)
            );

      // save errorIndicators in solutionVectors
      NumericVector<Number>* errorIndicators =
        NumericVector<Number>::build(_mesh->comm()).release();
      errorIndicators->init( mesh.n_elem(), mesh.n_local_elem(),
          system.get_dof_map().get_send_list() );
      for(unsigned int i=0; i<QoI_elementwise_error.size();i++)
        errorIndicators->set(i, QoI_elementwise_error[i] );
      errorIndicators->close();
      solutionVectors.insert( 
          std::pair<std::string,NumericVector<Number>* >(
            "errorIndicators", errorIndicators)
            );

      // save QoI value to solutionVectors
      system.assemble_qoi();
      NumericVector<Number>* qoiValue =
        NumericVector<Number>::build(_mesh->comm()).release();
      qoiValue->init( 1, 1, false, SERIAL);
      qoiValue->set(0, system.qoi[0] );
      solutionVectors.insert( 
          std::pair<std::string,NumericVector<Number>* >(
            "qoiValue", qoiValue)
            );

      

      return;
    }




  /********************************************//**
   * \brief 
   *
   * 
   ***********************************************/
  template<class T_S>
    void PhysicsLibmesh<T_S>::refine() 
    {
      std::cout << "  previous n_active_elem(): " 
        << _mesh->n_active_elem() << std::endl;

      _meshRefinement->uniformly_refine(1);

      _equationSystems->reinit();

      std::cout << "   refined n_active_elem(): " 
        << _mesh->n_active_elem() << std::endl;

    }

  /********************************************//**
   * \brief 
   *
   * 
   ***********************************************/
  template<class T_S>
    void PhysicsLibmesh<T_S>::refine(
        std::vector<Number>& errorIndicators
        ) 
    {

      std::cout << "  previous n_active_elem(): " 
        << _mesh->n_active_elem() << std::endl;

      libMesh::ErrorVector error_per_cell(errorIndicators.size());
      for(unsigned int i=0; i<errorIndicators.size();i++)
        error_per_cell[i] = errorIndicators[i];
      _meshRefinement->clean_refinement_flags();

      
      //TODO input options to control type of refinement
      /* _meshRefinement->flag_elements_by_error_tolerance (error_per_cell); */            
      /* _meshRefinement->flag_elements_by_error_fraction (error_per_cell); */             
      _meshRefinement->flag_elements_by_elem_fraction (error_per_cell);              
                   
      //TODO input options to control whether we coarsen as well
      _meshRefinement->refine_and_coarsen_elements();
      /* _meshRefinement->refine_elements(); */


      _equationSystems->reinit();

      std::cout << "   refined n_active_elem(): " 
        << _mesh->n_active_elem() << std::endl;

      return;
    }

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S> 
    libMesh::Mesh PhysicsLibmesh<T_S>::getMesh( ) const
    {
      return *_mesh;
    }


}

#endif // PHYSICS_LIBMESH_H
