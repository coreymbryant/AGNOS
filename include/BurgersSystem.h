
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
  : FEMSystem(es, name_in, number_in),
    _fe_family("LAGRANGE"), _fe_order(1),
    _analytic_jacobians(true) { 
      
      qoi.resize(1); }

  std::string & fe_family() { return _fe_family;  }
  unsigned int & fe_order() { return _fe_order;  }
  bool & analytic_jacobians() { return _analytic_jacobians; }

  // Postprocessing function which we are going to override for this application

  virtual void postprocess(void);

  Number &get_QoI_value(std::string type, unsigned int QoI_index)
    {
      if(type == "exact")
	{
	  return exact_QoI[QoI_index];
	}
      else
	{
	  return computed_QoI[QoI_index];
	}
    }


  double _mu;

  protected:
  // System initialization
  virtual void init_data ();

  // Context initialization
  virtual void init_context (DiffContext &context);

  // Element residual and jacobian calculations
  // Time dependent parts
  virtual bool element_time_derivative (bool request_jacobian,
					DiffContext &context);


  // Overloading the postprocess function

  virtual void element_postprocess(DiffContext &context);


  // Overloading the qoi function on elements

  virtual void element_qoi_derivative
    (DiffContext &context,
     const QoISet & qois);
  void element_qoi (DiffContext &context, const QoISet & /* qois */);

  // Overloading the qoi function on sides


  Number exact_solution (const Point&);

  // Variables to hold the computed QoIs

  Number computed_QoI[1];

  // Variables to read in the exact QoIs from l-shaped.in

  Number exact_QoI[1];

  // The FE type to use
  std::string _fe_family;
  unsigned int _fe_order;

  // Calculate Jacobians analytically or not?
  bool _analytic_jacobians;
};

void BurgersSystem::init_data ()
{
  this->add_variable ("u", static_cast<Order>(_fe_order),
                      Utility::string_to_enum<FEFamily>(_fe_family));
  /* this->add_variable ("u", FIRST); */

  /* GetPot infile("burgers.in"); */
  /* exact_QoI[0] = infile("QoI_0", 0.0); */
  // TODO set from input file
  // TODO update using setParameters
  _mu = 1.0;
  exact_QoI[0] = 5. + 0.5 * 4. * 0.1 * (log(cosh(10./4./_mu))) ;

  // Do the parent's initialization after variables are defined
  FEMSystem::init_data();

  this->extra_quadrature_order = 2;

  this->time_evolving(0);
  std::cout << "test"<< std::endl; 

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

  FEBase* side_fe = NULL;
  c.get_side_fe( 0, side_fe );

  side_fe->get_JxW();
  side_fe->get_phi();
  side_fe->get_dphi();
}

#define optassert(X) {if (!(X)) libmesh_error();}

// Assemble the element contributions to the stiffness matrix
bool BurgersSystem::element_time_derivative (bool request_jacobian,
					  DiffContext &context)
{
  // Are the jacobians specified analytically ?
  bool compute_jacobian = request_jacobian && _analytic_jacobians;
  /* std::cout << "compute_jacobian?:" << compute_jacobian << std::endl; */

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
      if (compute_jacobian)
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

  return compute_jacobian;
}

// Override the default DiffSystem postprocess function to compute the
// approximations to the QoIs
void BurgersSystem::postprocess()
{
  // Reset the array holding the computed QoIs
  computed_QoI[0] = 0.0;

  FEMSystem::postprocess();

  this->comm().sum(computed_QoI[0]);

}

// The exact solution to the singular problem,
// u_exact = r^(2/3)*sin(2*theta/3). We use this to set the Dirichlet boundary conditions
Number BurgersSystem::exact_solution(const Point& p)// xyz location
{
  const Real x1 = p(0);

  return 0.5 * (1 + tanh(x1 / 4. / _mu) ) ;

}

void BurgersSystem::element_postprocess (DiffContext &context)

{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  FEBase* elem_fe = NULL;
  c.get_element_fe( 0, elem_fe );

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = elem_fe->get_JxW();

  const std::vector<Point> &xyz = elem_fe->get_xyz();

  // The number of local degrees of freedom in each variable

  unsigned int n_qpoints = c.get_element_qrule().n_points();

  // The function R = int_{omega} T dR
  // omega is a subset of Omega (the whole domain), omega = [0.75, 1.0] x [0.0, 0.25]

  Number dQoI_0 = 0.;

  // Loop over quadrature points

  for (unsigned int qp = 0; qp != n_qpoints; qp++)
    {
      // Get co-ordinate locations of the current quadrature point
      const Real x = xyz[qp](0);

  	  // Get the solution value at the quadrature point
  	  Number u = c.interior_value(0, qp);

  	  // Update the elemental increment dR for each qp
      if ( x >=0.0)
        dQoI_0 += JxW[qp] * u;
    }

  // Update the computed value of the global functional R, by adding the contribution from this element

  computed_QoI[0] = computed_QoI[0] + dQoI_0;

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
      const Real x = q_point[qp](0);
      Number u = c.interior_value(0, qp);

      if (x >=0.0)
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
      const Real x = q_point[qp](0);

      if ( x >=0.0 )
        for (unsigned int i=0; i != n_u_dofs; i++)
          Q(i) += JxW[qp] *phi[i][qp] ;

    } // end of the quadrature point qp-loop
}


}
#endif // BURGERS_SYSTEM_H
