
#ifndef DIFFUSION_SYSTEM_H
#define DIFFUSION_SYSTEM_H

#include "agnosDefines.h"

// libmesh includes
#include "libmesh/enum_fe_family.h"
#include "libmesh/fem_system.h"
#include "libmesh/qoi_set.h"
#include "libmesh/system.h"
#include "libmesh/fe_base.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/parallel.h"
#include "libmesh/fem_context.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/elem.h"
#include "libmesh/point.h"
#include "libmesh/gnuplot_io.h"
#include "libmesh/zero_function.h"
#include "libmesh/const_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"

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
      }


    std::vector<double> _ck,_lambdak,_xik;
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
    unsigned int u_var = this->add_variable ("u",FIRST,HIERARCHIC) ;
    this->time_evolving(u_var);

    // define coeffs for KLE
    _ck.resize(10);
    _ck[0] = 10.;
    for(unsigned int i=1; i<_ck.size(); i++)
      _ck[i] = _ck[i-1]/10.;

    // define initial param values
    _xik.resize(10);
    for(unsigned int i=0; i<_xik.size(); i++)
      _xik[i] = 1.;


    // define periods for KLE
    _lambdak.resize(10);
    _lambdak[0] = 2.3926505e-01 ; 
    _lambdak[1] = 3.1340351e-01 ;
    _lambdak[2] = 4.0664199e-01;
    _lambdak[3] = 6.4388576e-01;
    _lambdak[4] = 7.4053791e-01;
    _lambdak[5] = 2.0492245e+00;
    _lambdak[6] = 3.1054506e+00;
    _lambdak[7] = 3.3967834e+00;
    _lambdak[8] = 8.7311801e+00;
    _lambdak[9] = 9.9791010e+00;

    /** set up boundary conditions */
    std::set<boundary_id_type> allBoundaries;
    allBoundaries.insert(0);
    allBoundaries.insert(1);
    allBoundaries.insert(2);
    allBoundaries.insert(3);

    std::vector<unsigned int> variables(1,u_var);
    
    ZeroFunction<double> zf;

    this->get_dof_map().add_dirichlet_boundary(
        DirichletBoundary(allBoundaries,variables,&zf) );
    if (AGNOS_DEBUG)
      std::cout << "post BC set up" << std::endl;

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
    const std::vector<Point > &q_point = elem_fe->get_xyz();

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
    unsigned int n_qpoints = c.get_element_qrule()->n_points();



    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        // Compute the solution gradient at the Newton iterate
        Gradient grad_u = c.interior_gradient(0, qp);
        Number u = c.interior_value(0, qp);


        // calculate coeff value
        const Real x1 = q_point[qp](0);
        const Real x2 = q_point[qp](1);
        double logCoeff = 0;
        for (unsigned int i=0; i<_ck.size();i++)
          logCoeff += _ck[i] * _xik[i] * 
            std::sin(_lambdak[i]*M_PI*x1) * std::cos(_lambdak[i]*M_PI*x2)  ;


        // The residual contribution from this element
        for (unsigned int i=0; i != n_u_dofs; i++)
          F(i) += JxW[qp] * ( 
              // (f,v) or (-K gradu_true dot gradv) 
               -10.0 * phi[i][qp]
              // - (u_x,v_x)
              + std::exp(logCoeff) * grad_u * dphi[i][qp] 
              );
        if (request_jacobian)
          for (unsigned int i=0; i != n_u_dofs; i++)
            for (unsigned int j=0; j != n_u_dofs; ++j)
        // The analytic jacobian
              K(i,j) += JxW[qp]*(
                  // (du_x,v_x)
                  1.0 * std::exp(logCoeff) * (dphi[j][qp] * dphi[i][qp] )       
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
    unsigned int n_qpoints = c.get_element_qrule()->n_points();

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

      Q[0] += JxW[qp] * u * 10./M_PI * 
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
    unsigned int n_qpoints = c.get_element_qrule()->n_points();

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
        Q(i) += JxW[qp] * phi[i][qp] * 10./M_PI * 
          std::exp( -10. * std::pow(x-0.5,2.) -10. * std::pow(y-0.5,2.) );

    } // end of the quadrature point qp-loop
  }


#endif // DIFFUSION_SYSTEM_H
