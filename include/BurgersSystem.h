
#ifndef BURGERS_SYSTEM_H
#define BURGERS_SYSTEM_H

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

// Bring in everything from the libMesh namespace
using namespace libMesh;

namespace AGNOS{

void write_output(EquationSystems &es,
		  unsigned int a_step,       // The adaptive step count
		  std::string solution_type) // primal or adjoint solve
{
  MeshBase &mesh = es.get_mesh();

  std::ostringstream file_name_gp;
  file_name_gp << solution_type
                << ".out.gp."
                << std::setw(2)
                << std::setfill('0')
                << std::right
                << a_step;

  GnuPlotIO(mesh).write_equation_systems
    (file_name_gp.str(), es);
}



// FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
// but we must specify element residuals
class BurgersSystem : public FEMSystem
{
public:
  // Constructor
  BurgersSystem(EquationSystems& es,
               const std::string& name_in,
               const unsigned int number_in)
  : FEMSystem(es, name_in, number_in)
  { 
      qoi.resize(1); 
      _mu=1.0;
  }



  double _mu;
  double _L;
  double _uMinus;
  double _uPlus;

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


  Number exact_solution (const Point&, const Real);

};

void BurgersSystem::init_data ()
{
  unsigned int u_var = this->add_variable ("u", FIRST);
  this->time_evolving(u_var);

  /* this->extra_quadrature_order = 2; */


  //---------------------------------------------
  /** set up boundary conditions */
  if (AGNOS_DEBUG)
    std::cout << "DEBUG: pre BC set up" << std::endl;
  std::set<boundary_id_type> minusBoundaries;
  minusBoundaries.insert(0);
  std::set<boundary_id_type> plusBoundaries;
  plusBoundaries.insert(1);

  std::vector<unsigned int> variables(1,u_var);

  Point lp(-1.*_L);
  Point rp(1.*_L);
  _uPlus  = this->exact_solution(rp,0);
  _uMinus = this->exact_solution(lp,0);

  std::cout << std::setprecision(17) << "uPlus = " << _uPlus << std::endl;
  std::cout << std::setprecision(17) << "uMinus = " << _uMinus << std::endl;

  ConstFunction<double> uPlus(_uPlus);
  ConstFunction<double> uMinus(_uMinus);

  this->get_dof_map().add_dirichlet_boundary(
      DirichletBoundary(minusBoundaries,variables,&uMinus) );
  this->get_dof_map().add_dirichlet_boundary(
      DirichletBoundary(plusBoundaries,variables,&uPlus) );

  if (AGNOS_DEBUG)
    std::cout << "DEBUG: post BC set up" << std::endl;
  //---------------------------------------------

  // Do the parent's initialization after variables are defined
  FEMSystem::init_data();

}

void BurgersSystem::init_context(DiffContext &context)
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
bool BurgersSystem::element_time_derivative (bool request_jacobian,
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
        F(i) += JxW[qp]*(
            // mu(u_x,v_x)
              -1.0 * _mu * (grad_u * dphi[i][qp] )    
            // -1/2 (u (1-u) , v_x)
            + 0.5 * u * ( 1. - u ) * dphi[i][qp](0)
            );
      if (request_jacobian)
        for (unsigned int i=0; i != n_u_dofs; i++)
          for (unsigned int j=0; j != n_u_dofs; ++j)
	    // The analytic jacobian
            K(i,j) += JxW[qp]*(
                // mu(du_x,v_x)
                -1.0 * _mu * (dphi[j][qp] * dphi[i][qp] )       
                // -1/2 (du - 2 u du , v_x)
                + 0.5 * phi[j][qp] * (1. - 2. * u) * dphi[i][qp](0)
              );
    } // end of the quadrature point qp-loop

  return request_jacobian;
}


// exact solution
Number BurgersSystem::exact_solution(const Point& p, const Real t)// xyz location
{
  const Real x1 = p(0);

  return 0.5 * (1 + tanh(x1 / 4. / _mu) ) ;

}


void BurgersSystem::element_qoi (DiffContext &context,
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
      Number u = c.interior_value(0, qp);

      Q[0] += JxW[qp] * u ;

    } // end of the quadrature point qp-loop
}



void BurgersSystem::element_qoi_derivative (DiffContext &context,
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

      for (unsigned int i=0; i != n_u_dofs; i++)
        Q(i) += JxW[qp] *phi[i][qp] ;

    } // end of the quadrature point qp-loop
}


}
#endif // BURGERS_SYSTEM_H
