

#ifndef PHYSICS_HIGHER_ORDER_DIFFUSION_H
#define PHYSICS_HIGHER_ORDER_DIFFUSION_H

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

  class DiffusionSystem : public FEMSystem
  {
  public:
    // Constructor
    DiffusionSystem(EquationSystems& es,
                 const std::string& name_in,
                 const unsigned int number_in)
    : FEMSystem(es, name_in, number_in)
      { 
        qoi.resize(1); 
        _coeff=1.0;
      }


    double _coeff;
    double exact_solution (const Point&);

    protected:
    // System initialization
    virtual void init_data ();

    // Context initialization
    virtual void init_context (DiffContext &context);

    // Element residual and jacobian calculations
    // Time dependent parts
    virtual bool element_time_derivative (bool request_jacobian,
            DiffContext &context);


    // Overloading the qoi function on elements
    virtual void element_qoi_derivative
      (DiffContext &context,
       const QoISet & qois);
    void element_qoi (DiffContext &context, const QoISet & /* qois */);


  };

  void DiffusionSystem::init_data ()
  {
    unsigned int u_var = this->add_variable ("u",FIRST) ;
    this->time_evolving(u_var);


    // Do the parent's initialization after variables are defined
    FEMSystem::init_data();

  }

  void DiffusionSystem::init_context(DiffContext &context)
  {
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

    // Now make sure we have requested all the data
    // we need to build the linear system.
    FEBase* elem_fe = NULL;
    c.get_element_fe( 0, elem_fe );
    elem_fe->get_JxW();
    elem_fe->get_phi();
    elem_fe->get_dphi();

  }

#define optassert(X) {if (!(X)) libmesh_error();}

  // Assemble the element contributions to the stiffness matrix
  bool DiffusionSystem::element_time_derivative (bool request_jacobian,
              DiffContext &context)
  {
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

    // First we get some references to cell-specific data that
    // will be used to assemble the linear system.
    FEBase* elem_fe = NULL;
    c.get_element_fe( 0, elem_fe );

    // Element Jacobian * quadrature weights for interior integration
    const std::vector<Real> &JxW = elem_fe->get_JxW();

    // Element basis functions
    const std::vector<std::vector<RealGradient> > &dphi = elem_fe->get_dphi();
    const std::vector<std::vector<Real> > &phi = elem_fe->get_phi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = c.get_dof_indices(0).size();

    // The subvectors and submatrices we need to fill:
    DenseSubMatrix<Number> &K = c.get_elem_jacobian(0,0);
    DenseSubVector<Number> &F = c.get_elem_residual(0);

    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    unsigned int n_qpoints = c.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        // Compute the solution gradient at the Newton iterate
        Gradient grad_u = c.interior_gradient(0, qp);
        Number u = c.interior_value(0, qp);

        // The residual contribution from this element
        for (unsigned int i=0; i != n_u_dofs; i++)
          F(i) += JxW[qp] * ( 
              // (f,v) or (-K gradu_true dot gradv) 
               10.0 * phi[i][qp]
              // - (u_x,v_x)
              + _coeff * grad_u * dphi[i][qp] 
              );
        if (request_jacobian)
          for (unsigned int i=0; i != n_u_dofs; i++)
            for (unsigned int j=0; j != n_u_dofs; ++j)
        // The analytic jacobian
              K(i,j) += JxW[qp]*(
                  // (du_x,v_x)
                  1.0 * _coeff * (dphi[j][qp] * dphi[i][qp] )       
                );
      } // end of the quadrature point qp-loop

    return request_jacobian;
  }


  // exact solution
  double DiffusionSystem::exact_solution(const Point& p)// xyz location
  {
    std::cerr << "\nERROR: not defined\n" << std::endl;
    return 0.;

  }


  void DiffusionSystem::element_qoi (DiffContext &context,
                                              const QoISet & /* qois */)
  {
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

    // First we get some references to cell-specific data that
    // will be used to assemble the linear system.
    FEBase* elem_fe = NULL;
    c.get_element_fe( 0, elem_fe );

    // Element Jacobian * quadrature weights for interior integration
    const std::vector<Real> &JxW = elem_fe->get_JxW();

    // The element quadrature points
    const std::vector<Point > &q_point = elem_fe->get_xyz();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = c.get_dof_indices(0).size();
    unsigned int n_qpoints = c.get_element_qrule().n_points();

    // Fill the QoI RHS corresponding to this QoI. Since this is the 0th QoI
    // we fill in the [0][i] subderivatives, i corresponding to the variable index.
    // Our system has only one variable, so we only have to fill the [0][0] subderivative
    std::vector<Number> &Q = c.get_qois();

    // Loop over the qps
    for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      const Real x = q_point[qp](0);
      const Real y = q_point[qp](1);
      Number u = c.interior_value(0, qp);

      Q[0] += JxW[qp] * u * 10./3.14159 * 
        std::exp( -10. * std::pow(x-0.5,2.) -10. * std::pow(y-0.5,2.) );

    } // end of the quadrature point qp-loop
  }

  void DiffusionSystem::element_qoi_derivative (DiffContext &context,
                                              const QoISet & /* qois */)
  {
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

    // First we get some references to cell-specific data that
    // will be used to assemble the linear system.
    FEBase* elem_fe = NULL;
    c.get_element_fe( 0, elem_fe );

    // Element Jacobian * quadrature weights for interior integration
    const std::vector<Real> &JxW = elem_fe->get_JxW();

    // The basis functions for the element
    const std::vector<std::vector<Real> > &phi = elem_fe->get_phi();

    // The element quadrature points
    const std::vector<Point > &q_point = elem_fe->get_xyz();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = c.get_dof_indices(0).size();
    unsigned int n_qpoints = c.get_element_qrule().n_points();

    // Fill the QoI RHS corresponding to this QoI. Since this is the 0th QoI
    // we fill in the [0][i] subderivatives, i corresponding to the variable index.
    // Our system has only one variable, so we only have to fill the [0][0] subderivative
    DenseSubVector<Number> &Q = c.get_qoi_derivatives(0,0);

    // Loop over the qps
    for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      const Real x = q_point[qp](0);
      const Real y = q_point[qp](1);

      for (unsigned int i=0; i != n_u_dofs; i++)
        Q(i) += JxW[qp] * phi[i][qp] * 10./3.14159 * 
          std::exp( -10. * std::pow(x-0.5,2.) -10. * std::pow(x-0.5,2.) );

    } // end of the quadrature point qp-loop
  }


  /********************************************//**
   * \brief Example PhysicsLibmesh class - higher order diffusion example from
   * paper
   ***********************************************/
  template<class T_S, class T_P>
  class PhysicsHigherOrderDiffusion : public PhysicsLibmesh<T_S,T_P>
  {

    public:
      /** Constructor. Pass input file to provide setting to physics class */
      PhysicsHigherOrderDiffusion( const Communicator& comm_in, const GetPot& input );

      /** destructor */
      virtual ~PhysicsHigherOrderDiffusion( );


    protected:
      /** Geometry and boundary data */
      unsigned int    _nElem;
      
      /** coefficients */
      std::vector<double> _ck ;

      /** set parameter values */
      virtual void _setParameterValues( const T_S& parameterValues ) ;

  };

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsHigherOrderDiffusion<T_S,T_P>::~PhysicsHigherOrderDiffusion( )
  {
  }



/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsHigherOrderDiffusion<T_S,T_P>::PhysicsHigherOrderDiffusion(
        const Communicator& comm_in, 
        const GetPot& input 
      )
  :
    PhysicsLibmesh<T_S,T_P>(comm_in,input)
  {
    _nElem  = input("nElem",4);

    _ck.resize(10);
    _ck[0] = 10.;
    for(unsigned int i=1; i<_ck.size(); i++)
      _ck[i] = _ck[i-1]/10.;
    
    //------------------------
    // initialize mesh object
    this->_mesh = new libMesh::Mesh(this->_communicator);


    // build mesh refinement object 
    if (AGNOS_DEBUG)
      std::cout << "DEBUG: pre mesh_refinement " << std::endl;
    this->_buildMeshRefinement();
    //----------------------------------------------


    // build mesh 
    // MeshTools::Generation::build_square (
    libMesh::MeshTools::Generation::build_square(
        *this->_mesh,
        this->_nElem,this->_nElem,
        0.,1.,
        0.,1.,
        QUAD9);
    this->_mesh->print_info();
    //----------------------------------------------


    // define equation system
    this->_equationSystems 
      = new libMesh::EquationSystems(*this->_mesh);
    this->_system = &( 
        this->_equationSystems->template
        add_system<CatenarySystem>("HigherOrderDiffusion")
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
    void PhysicsHigherOrderDiffusion<T_S,T_P>::_setParameterValues( 
        const T_S& parameterValues ) 
    {
      static_cast<CatenarySystem*>(this->_system)->_coeff = parameterValues(0) ;
    }


}
#endif // PHYSICS_HIGHER_ORDER_DIFFUSION_H
