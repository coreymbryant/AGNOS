
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
    _nElem  = input("nElem",4);
    /* _uMinus = input("uMinus",(0.5 * ( 1 + std::tanh( -1.*_L / 4. / 1.0) ) )); */
    /* _uPlus  = input("uPlus",(0.5 * ( 1 + std::tanh( _L / 4. / 1.0) ) ) ); */

    // and some nonlinear solver parameters
    _nonlinearSteps      = input("nNonlinearSteps",15);
    _nonlinearTolerance  = input("nonlinearTolerance",1.e-9);
    //----------------------------------------------
    


    // initialize mesh object
    this->_mesh = new libMesh::Mesh(this->_communicator);
    //----------------------------------------------
    


    // build mesh refinement object 
    if (AGNOS_DEBUG)
      std::cout << "DEBUG: pre mesh_refinement " << std::endl;
    this->_buildMeshRefinement();
    //----------------------------------------------
    

    // build mesh 
    libMesh::MeshTools::Generation::build_line(
        *this->_mesh,this->_nElem,-1.*_L,_L,EDGE2);
    this->_mesh->print_info();

    //----------------------------------------------


    // define equation system
    this->_equationSystems 
      = new libMesh::EquationSystems(*this->_mesh);
    this->_system = 
      &( this->_equationSystems->template add_system<BurgersSystem>("Burgers") );
    static_cast<BurgersSystem*>(this->_system)->_L = _L ; 
    if (AGNOS_DEBUG)
      std::cout << "DEBUG: post add system" << std::endl;
    //----------------------------------------------
    

    // No transient time solver
    this->_system->time_solver =
        AutoPtr<TimeSolver>(new SteadySolver(*this->_system));
    {
      NewtonSolver *solver = new NewtonSolver(*this->_system);
      this->_system->time_solver->diff_solver() = AutoPtr<DiffSolver>(solver);
      
      //TODO read in these setting?
      solver->quiet                       = true;
      solver->verbose                     = false;
      solver->max_nonlinear_iterations    = 1;
      solver->continue_after_max_iterations = true;
      solver->continue_after_backtrack_failure = true;
      
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
    static_cast<BurgersSystem*>(this->_system)->_mu 
      = 1.0 + 0.62 * parameterValues(0) + 0.36 * parameterValues(1) ;
  }


}

#endif // PHYSICS_VISCOUS_BURGERS_H
