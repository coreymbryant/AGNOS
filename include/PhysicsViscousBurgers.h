
#ifndef PHYSICS_VISCOUS_BURGERS_H
#define PHYSICS_VISCOUS_BURGERS_H


#include "agnosDefines.h"
#include "PhysicsLibmesh.h"
#include "BurgersSystem.h"

// libmesh includes
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/newton_solver.h"
#include "libmesh/steady_solver.h"

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
      /** Constructor. Pass input file to provide setting to physics class */
      PhysicsViscousBurgers( const Communicator& comm_in, const GetPot& input );

      /** Destructor */
      /* virtual ~PhysicsViscousBurgers( ); */


    protected:
      /** Geometry and boundary data */
      double          _L;
      double          _uMinus;
      double          _uPlus;
      int             _nElem;

      Number          _mu;

      /** solver settings */ 
      unsigned int    _nonlinearTolerance;
      unsigned int    _nonlinearSteps;


      /** set parameter values */
      virtual void _setParameterValues( const T_S& parameterValues ) ;

      

  };




/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsViscousBurgers<T_S,T_P>::PhysicsViscousBurgers(
        const Communicator& comm_in, 
        const GetPot& input 
        ) 
  :
    PhysicsLibmesh<T_S,T_P>(comm_in,input)
  {

    // read in parameters unique to this model
    _L      = input("L",10.);
    _uMinus = input("uMinus",(0.5 * ( 1 + std::tanh( -1.*_L / 4. / 1.0) ) ));
    _uPlus  = input("uPlus",(0.5 * ( 1 + std::tanh( _L / 4. / 1.0) ) ) );
    _nElem  = input("nElem",4);

    // and some nonlinear solver parameters
    _nonlinearSteps      = input("nNonlinearSteps",15);
    _nonlinearTolerance  = input("nonlinearTolerance",1.e-9);

    //----------------------------------------------
    // construct a 1d mesh 
    this->_mesh = new libMesh::Mesh(this->_communicator,1);
    /* libMesh::Mesh* mesh = new libMesh::Mesh(this->_communicator); */
    libMesh::MeshTools::Generation::build_line(
        *this->_mesh, _nElem, -1.*_L, _L, EDGE3);
    this->_mesh->print_info();

    // define equation system
    this->_equationSystems 
      = new libMesh::EquationSystems(*this->_mesh);

    // add system and set parent pointer 
    this->_system = 
      &( this->_equationSystems->template add_system<BurgersSystem>("Burgers") );
    
    // No transient time solver
    this->_system->time_solver =
        AutoPtr<TimeSolver>(new SteadySolver(*this->_system));

    // Nonlinear solver options
    {
      NewtonSolver *solver = new NewtonSolver(*this->_system);
      this->_system->time_solver->diff_solver() = AutoPtr<DiffSolver>(solver);

      //TODO read in these setting?
    solver->quiet                       = true;
    solver->verbose                     = false;
    /* solver->max_nonlinear_iterations    = param.max_nonlinear_iterations; */
    /* solver->minsteplength               = param.min_step_length; */
    /* solver->relative_step_tolerance     = param.relative_step_tolerance; */
    /* solver->absolute_residual_tolerance = param.absolute_residual_tolerance; */
    /* solver->relative_residual_tolerance = param.relative_residual_tolerance; */
    /* solver->require_residual_reduction  = param.require_residual_reduction; */
    /* solver->linear_tolerance_multiplier = param.linear_tolerance_multiplier; */
    /* if (system.time_solver->reduce_deltat_on_diffsolver_failure) */
    /*   { */
	solver->continue_after_max_iterations = true;
	solver->continue_after_backtrack_failure = true;
      /* } */

    // And the linear solver options
    /* solver->max_linear_iterations       = param.max_linear_iterations; */
    /* solver->initial_linear_tolerance    = param.initial_linear_tolerance; */
    /* solver->minimum_linear_tolerance    = param.minimum_linear_tolerance; */
    }

    this->_equationSystems->init ();
    //----------------------------------------------

    //----------------------------------------------
    //---- set up nonlinear solver
    // TODO do we need to initialize ourselves
    // initalize the nonlinear solver
    /* this->_equationSystems->template */
    /*   get_system<NonlinearImplicitSystem>("Burgers").nonlinear_solver->init(); */

    /* //---------------------------------------------- */
    /* // TODO set from input file ? */
    /* // set solver settings */
    /* this->_equationSystems->parameters.template set<Real> */
    /*   ("nonlinear solver relative residual tolerance") = _nonlinearTolerance; */
    /* this->_equationSystems->parameters.template set<Real> */
    /*   ("nonlinear solver absolute residual tolerance") = 1.e-35; */
    /* this->_equationSystems->parameters.template set<Real> */
    /*   ("nonlinear solver absolute step tolerance") = 1.e-12; */
    /* this->_equationSystems->parameters.template set<Real> */
    /*   ("nonlinear solver relative step tolerance") = 1.e-12; */
    /* this->_equationSystems->parameters.template set<unsigned int> */
    /*   ("nonlinear solver maximum iterations") = 50; */
    /* //---------------------------------------------- */
    

    //----------------------------------------------
    // set up QoISet object 
    this->_qois = new libMesh::QoISet;
    std::vector<unsigned int> qoi_indices;
    qoi_indices.push_back(0);
    this->_qois->add_indices(qoi_indices);

    // weight the qois (in case we have more than 1)
    this->_qois->set_weight(0, 1.0);
    //----------------------------------------------

    if (AGNOS_DEBUG)
      std::cout << "test: pre mesh_refinement " << std::endl;
    
    //----------------------------------------------
    // build mesh refinement object 
    this->_buildMeshRefinement();

    // build error estimator object
    this->_buildErrorEstimator();
    //----------------------------------------------

    if (AGNOS_DEBUG)
      std::cout << "test: post error estimator " << std::endl;
    

    this->_equationSystems->init ();
    if (AGNOS_DEBUG)
      std::cout << "test: pre print info " << std::endl;

    // Print information about the mesh and system to the screen.

    this->_equationSystems->print_info();

    std::cout << "test: system initialized" << std::endl;
  }


  template<class T_S,class T_P>
  void PhysicsViscousBurgers<T_S,T_P>::_setParameterValues(
    const T_S& parameterValues )
  {
    //TODO 
  }


}

#endif // PHYSICS_VISCOUS_BURGERS_H
