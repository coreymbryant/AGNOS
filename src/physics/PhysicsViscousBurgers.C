


#include "PhysicsViscousBurgers.h"
#include "BurgersSystem.h"

// libmesh includes
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/petsc_diff_solver.h"
#include "libmesh/newton_solver.h"
#include "libmesh/steady_solver.h"

namespace AGNOS
{

  template<class T_S, class T_P>
  T_P PhysicsViscousBurgers<T_S,T_P>::exactQoi()
  {
    T_P resultVector(1);
    /* resultVector(0) = 10.; */
    Number mu = dynamic_cast<BurgersSystem*>(this->_system)->_mu ;
    resultVector(0) = 0.5 
      + 2.*mu * std::log( std::cosh( 1./4./mu));
    if(AGNOS_DEBUG)
    {
      std::cout << "this->mu: " << mu << std::endl;
      std::cout << "cosh( ) : " << std::cosh( 1./4./mu) << std::endl;
      std::cout << "log( ): : " << std::log( std::cosh( 1./4./mu)) << std::endl;
      std::cout << "exactQoi: " << resultVector(0) << std::endl;
    }
    return resultVector;
  }

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
    // define available solution names
    this->_availableSolutions.insert("primal");
    this->_availableSolutions.insert("adjoint");
    this->_availableSolutions.insert("qoi");
    this->_availableSolutions.insert("errorEstimate");
    this->_availableSolutions.insert("errorIndicators");
    this->_availableSolutions.insert("exactQoi");


    // read in parameters unique to this model
    _L      = input("L",10.);
    _nElem  = input("nElem",4);
    /* _uMinus = input("uMinus",(0.5 * ( 1 + std::tanh( -1.*_L / 4. / 1.0) ) )); */
    /* _uPlus  = input("uPlus",(0.5 * ( 1 + std::tanh( _L / 4. / 1.0) ) ) ); */

    // and some nonlinear solver parameters
    _nonlinearSteps      = input("nNonlinearSteps",15);
    _nonlinearTolerance  = input("nonlinearTolerance",1.e-9);
    //----------------------------------------------

    //----------------------------------------------
    // initialize mesh object
    // build temporary mesh 
    libMesh::Mesh mesh(this->_communicator);
    libMesh::MeshTools::Generation::build_line(
        mesh,this->_nElem,-1.*_L,_L,EDGE2);
    // deep copy to PhysicsLibmesh object pointer
    this->_mesh = new libMesh::Mesh(mesh);
    this->_mesh->print_info();
    //----------------------------------------------
    
    // build mesh refinement object 
    if (AGNOS_DEBUG)
      std::cout << "DEBUG: pre mesh_refinement " << std::endl;
    this->_buildMeshRefinement();
    //----------------------------------------------
    
    // define equation system
    this->_equationSystems 
      = new libMesh::EquationSystems(*this->_mesh);
    BurgersSystem& burgersSystem = 
      this->_equationSystems->template add_system<BurgersSystem>("Burgers") ;
    this->_system = &( burgersSystem );
    burgersSystem._L = _L ;
    if (AGNOS_DEBUG)
      std::cout << "DEBUG: post add system" << std::endl;
    //----------------------------------------------
    

    // No transient time solver
    burgersSystem.time_solver =
        UniquePtr<TimeSolver>(new SteadySolver(burgersSystem) );
    {
      NewtonSolver *solver = new NewtonSolver(burgersSystem);
      /* PetscDiffSolver *solver = new PetscDiffSolver(burgersSystem); */
      burgersSystem.time_solver->diff_solver() = UniquePtr<DiffSolver>(solver);
      
      //TODO read in these setting?
      solver->quiet                       = true;
      solver->verbose                     = false;
      solver->max_nonlinear_iterations    = 100;
      solver->relative_step_tolerance     = 1.e-99;
      solver->absolute_residual_tolerance = 1.e-16;
      solver->relative_residual_tolerance = 1.e-99;

      // continue if fail to converge 
      solver->continue_after_max_iterations = true;
      solver->continue_after_backtrack_failure = true;
      
      // And the linear solver options
      solver->max_linear_iterations       = 10000;
      solver->initial_linear_tolerance    = 1.e-16;
      solver->minimum_linear_tolerance    = TOLERANCE*TOLERANCE;
      
    }
    if (AGNOS_DEBUG)
      std::cout << "post solver set up" << std::endl;

    //---------------------------------------------
    /** initialize equation system */
    this->_equationSystems->init ();
    if (AGNOS_DEBUG)
      std::cout << "post init system" << std::endl;
    


    //----------------------------------------------
    // set up QoISet object 
    this->_qois = new libMesh::QoISet;
    std::vector<unsigned int> qoi_indices;
    qoi_indices.push_back(0);
    this->_qois->add_indices(qoi_indices);

    // weight the qois (in case we have more than 1)
    this->_qois->set_weight(0, 1.0);
    //----------------------------------------------

    
    // build error estimator object
    this->_buildErrorEstimator();
    if (AGNOS_DEBUG)
      std::cout << "debug: post error estimator " << std::endl;
    //----------------------------------------------
    
    // build mesh refinement object 
    this->_buildMeshRefinement();
    if (AGNOS_DEBUG)
      std::cout << "test: pre mesh_refinement " << std::endl;
    //----------------------------------------------



    // Print information about the mesh and system to the screen.
    this->_equationSystems->print_info();
    if (AGNOS_DEBUG)
      std::cout << "DEBUG: leaving model specific setup" << std::endl;
    //----------------------------------------------
  }


  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S,class T_P>
  void PhysicsViscousBurgers<T_S,T_P>::_setParameterValues(
    const T_S& parameterValues )
  {
    /* *(this->_system)._mu */ 
    /*   = 1.0 + 0.62 * parameterValues(0) + 0.36 * parameterValues(1) ; */
    dynamic_cast<BurgersSystem*>(this->_system)->_mu 
      = 1.0 + 0.62 * parameterValues(0) + 0.36 * parameterValues(1) ;

    if (AGNOS_DEBUG)
      std::cout << "DEBUG: system mu: " 
        << dynamic_cast<BurgersSystem*>(this->_system)->_mu  << std::endl;

    dynamic_cast<BurgersSystem*>(this->_system)->init_bcs(); 

    this->_system->reinit();
  }

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S,class T_P>
  PhysicsViscousBurgers<T_S,T_P>::~PhysicsViscousBurgers()
  {
    if( this->_mesh != NULL            ){ delete this->_mesh; }
    if( this->_meshRefinement != NULL  ){ delete this->_meshRefinement; }
    if( this->_errorEstimator != NULL  ){ delete this->_errorEstimator; }
    if( this->_qois != NULL            ){ delete this->_qois; }
    if( this->_equationSystems != NULL ){ delete this->_equationSystems; }
  }

  template class
    PhysicsViscousBurgers<libMesh::DenseVector<double>, libMesh::DenseVector<double> >;

}

