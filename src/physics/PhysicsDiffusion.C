

#include "PhysicsDiffusion.h"
#include "DiffusionSystem.h"

// libmesh includes
#include "libmesh/mesh_generation.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/steady_solver.h"
#include "libmesh/newton_solver.h"
        
#define _USE_MATH_DEFINES



namespace AGNOS
{

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsDiffusion<T_S,T_P>::~PhysicsDiffusion( )
  {
  }



/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsDiffusion<T_S,T_P>::PhysicsDiffusion(
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

    _nElem  = input("nElem",100);

    
    //------------------------
    // initialize mesh object
    this->_mesh = new libMesh::Mesh(this->_communicator,2);


    // build mesh refinement object 
    if (AGNOS_DEBUG)
      std::cout << "DEBUG: pre mesh_refinement " << std::endl;
    this->_buildMeshRefinement();
    //----------------------------------------------


    // build mesh 
    // MeshTools::Generation::build_square (
    // nElem per side
    libMesh::MeshTools::Generation::build_square(
        *static_cast<libMesh::Mesh*>(this->_mesh),
        this->_nElem,
        this->_nElem,
        0.,1.,
        0.,1.,
        TRI6);
    this->_mesh->print_info();

    //----------------------------------------------


    // define equation system
    this->_equationSystems 
      = new libMesh::EquationSystems(*this->_mesh);
    DiffusionSystem& diffusionSystem = 
        this->_equationSystems->template
        add_system<DiffusionSystem>("Diffusion") ;
    this->_system = &( diffusionSystem ) ;
    if (AGNOS_DEBUG)
      std::cout << "DEBUG: post add system" << std::endl;
    //----------------------------------------------
    

    // No transient time solver
    diffusionSystem.time_solver =
        AutoPtr<TimeSolver>(new SteadySolver(diffusionSystem));
    {
      NewtonSolver *solver = new NewtonSolver(diffusionSystem);
      diffusionSystem.time_solver->diff_solver() = AutoPtr<DiffSolver>(solver);
      
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
    void PhysicsDiffusion<T_S,T_P>::_setParameterValues( 
        const T_S& parameterValues ) 
    {

      if ( parameterValues.size() 
          == 
          static_cast<DiffusionSystem*>(this->_system)->_xik.size() )
      {
        for (unsigned int i=0;
            i< static_cast<DiffusionSystem*>(this->_system)->_xik.size(); 
            i++)
          static_cast<DiffusionSystem*>(this->_system)->_xik[i] 
            = parameterValues(i) ;
        std::cout << "xi = " ;
        for(unsigned int i=0; i<parameterValues.size();i++)
          std::cout << parameterValues(i) << "  " ;
        std::cout << std::endl;
      }
      else
      {
        std::cerr << " ERROR: incorrect dimension for parameter space " 
          << "\n"
          << " this model requires a parameter dimension of 10. " 
          << " given size = " << parameterValues.size() 
          << " model size = " <<
            static_cast<DiffusionSystem*>(this->_system)->_xik.size()        
          << std::endl;
        exit(1);
      }

    }


  template class
    PhysicsDiffusion<libMesh::DenseVector<double>, libMesh::DenseVector<double> >;

}
