

#ifndef PHYSICS_CATENARY_LIBMESH_H
#define PHYSICS_CATENARY_LIBMESH_H


#include "agnosDefines.h"
#include "PhysicsLibmesh.h"
#include "CatenarySystem.h"

// libmesh includes
#include "libmesh/mesh_generation.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/steady_solver.h"
#include "libmesh/newton_solver.h"
        



namespace AGNOS
{

  /********************************************//**
   * \brief Example PhysicsLibmesh class - catenary chain (solution computed
   * using libmesh)
   *
   * A simple 1D example useful for testing purposes
   *
   * 
   ***********************************************/
  template<class T_S, class T_P>
  class PhysicsCatenaryLibmesh : public PhysicsLibmesh<T_S,T_P>
  {

    public:
      /** Constructor. Pass input file to provide setting to physics class */
      PhysicsCatenaryLibmesh( const Communicator& comm_in, const GetPot& input );

      /** destructor */
      virtual ~PhysicsCatenaryLibmesh( );

      /** Redefine exactQoi for this model */
      T_P exactQoi( )
      {
        T_P resultVector(1);
        resultVector(0) = _forcing / 
          (8. * static_cast<CatenarySystem*>(this->_system)->_coeff)  ;
        return resultVector;
      }

    protected:
      /** Geometry and boundary data */
      double          _min;
      double          _max;
      unsigned int    _nElem;

      double          _forcing;

      /** set parameter values */
      virtual void _setParameterValues( const T_S& parameterValues ) ;

  };

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsCatenaryLibmesh<T_S,T_P>::~PhysicsCatenaryLibmesh( )
  {
  }



/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsCatenaryLibmesh<T_S,T_P>::PhysicsCatenaryLibmesh(
        const Communicator& comm_in, 
        const GetPot& input 
      )
  :
    PhysicsLibmesh<T_S,T_P>(comm_in,input)
  {
    _min    = input("min",0.);
    _max    = input("max",1.);
    _nElem  = input("nElem",4);

    _forcing = input("forcing",-10.);

    
    //------------------------
    // initialize mesh object
    this->_mesh = new libMesh::Mesh(this->_communicator);


    // build mesh refinement object 
    if (AGNOS_DEBUG)
      std::cout << "DEBUG: pre mesh_refinement " << std::endl;
    this->_buildMeshRefinement();
    //----------------------------------------------


    // build mesh 
    libMesh::MeshTools::Generation::build_line(
        *this->_mesh,this->_nElem,_min,_max,EDGE2);
    this->_mesh->print_info();
    //----------------------------------------------


    // define equation system
    this->_equationSystems 
      = new libMesh::EquationSystems(*this->_mesh);
    this->_system = &( 
        this->_equationSystems->template add_system<CatenarySystem>("Catenary")
        ) ;
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

    


    // Print information about the mesh and system to the screen.
    this->_equationSystems->print_info();
    if (AGNOS_DEBUG)
      std::cout << "debug: leaving model specific setup" << std::endl;
  }

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S, class T_P>
    void PhysicsCatenaryLibmesh<T_S,T_P>::_setParameterValues( 
        const T_S& parameterValues ) 
    {
      static_cast<CatenarySystem*>(this->_system)->_coeff = parameterValues(0) ;
    }


}

#endif // PHYSICS_CATENARY_LIBMESH_H
