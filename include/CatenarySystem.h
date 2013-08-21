

#ifndef CATENARY_SYSTEM_H
#define CATENARY_SYSTEM_H

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



// FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
// but we must specify element residuals
class CatenarySystem : public FEMSystem
{
public:
  // Constructor
  CatenarySystem(EquationSystems& es,
               const std::string& name_in,
               const unsigned int number_in)
  : FEMSystem(es, name_in, number_in)
    { 
      qoi.resize(1); 
      _coeff=1.0;
    }


  /* // Postprocessing function which we are going to override for this application */

  /* virtual void postprocess(void); */

  double _forcing; 
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

  // Overloading the postprocess function
  /* virtual void element_postprocess(DiffContext &context); */


  // Overloading the qoi function on elements
  virtual void element_qoi_derivative
    (DiffContext &context,
     const QoISet & qois);
  void element_qoi (DiffContext &context, const QoISet & /* qois */);


};

void CatenarySystem::init_data ()
{
  this->add_variable ("u",FIRST) ;
  // TODO update using setParameters
  /* _coeff = 1.0; */
  _forcing = -10.0;

  // Do the parent's initialization after variables are defined
  FEMSystem::init_data();

  this->time_evolving(0);
}

void CatenarySystem::init_context(DiffContext &context)
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
bool CatenarySystem::element_time_derivative (bool request_jacobian,
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
            // (f,v) 
            phi[i][qp] * _forcing 
            // - (u_x,v_x)
            -1.0 * _coeff * (grad_u * dphi[i][qp] )    
            );
      if (request_jacobian)
        for (unsigned int i=0; i != n_u_dofs; i++)
          for (unsigned int j=0; j != n_u_dofs; ++j)
	    // The analytic jacobian
            K(i,j) += JxW[qp]*(
                // (du_x,v_x)
                -1.0 * _coeff * (dphi[j][qp] * dphi[i][qp] )       
              );
    } // end of the quadrature point qp-loop

  return request_jacobian;
}

/* // Override the default DiffSystem postprocess function to compute the */
/* // approximations to the QoIs */
/* void CatenarySystem::postprocess() */
/* { */
/*   // Reset the array holding the computed QoIs */
/*   computed_QoI[0] = 0.0; */

/*   FEMSystem::postprocess(); */

/*   this->comm().sum(computed_QoI[0]); */

/* } */


// exact solution
double CatenarySystem::exact_solution(const Point& p)// xyz location
{
  const Real x1 = p(0);

  return _forcing / (2.*_coeff) * x1 * (1.-x1) ;

}

/* void CatenarySystem::element_postprocess (DiffContext &context) */

/* { */
/*   FEMContext &c = libmesh_cast_ref<FEMContext&>(context); */

/*   FEBase* elem_fe = NULL; */
/*   c.get_element_fe( 0, elem_fe ); */

/*   // Element Jacobian * quadrature weights for interior integration */
/*   const std::vector<Real> &JxW = elem_fe->get_JxW(); */

/*   const std::vector<Point> &xyz = elem_fe->get_xyz(); */

/*   // The number of local degrees of freedom in each variable */

/*   unsigned int n_qpoints = c.get_element_qrule().n_points(); */

/*   // The function R = int_{omega} T dR */
/*   // omega is a subset of Omega (the whole domain), omega = [0.75, 1.0] x [0.0, 0.25] */

/*   Number dQoI_0 = 0.; */

/*   // Loop over quadrature points */

/*   for (unsigned int qp = 0; qp != n_qpoints; qp++) */
/*     { */
/*       // Get co-ordinate locations of the current quadrature point */
/*       const Real x = xyz[qp](0); */

/*   	  // Get the solution value at the quadrature point */
/*   	  Number u = c.interior_value(0, qp); */

/*   	  // Update the elemental increment dR for each qp */
/*       if ( x >=0.0) */
/*         dQoI_0 += JxW[qp] * u; */
/*     } */

/*   // Update the computed value of the global functional R, by adding the contribution from this element */

/*   computed_QoI[0] = computed_QoI[0] + dQoI_0; */

/* } */

void CatenarySystem::element_qoi (DiffContext &context,
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

  for(unsigned int s=0; s<c.get_elem().n_sides(); s++)
  {
    // Get co-ordinate locations of the current quadrature point
    Point p = c.get_elem().point(s) ;
    const Real x = p(0);

    // Get the solution value at the quadrature point
    Number u = c.point_value(0, p);

    // Update the elemental increment dR for each qp
    if ( x == 0.5 )
      Q[0] += u * 0.5; /** 1/2 because node will get counted once for each
                           element it belongs to  */
  }
}

void CatenarySystem::element_qoi_derivative (DiffContext &context,
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

  for(unsigned int s=0; s<c.get_elem().n_sides(); s++)
  {
      Point p = c.get_elem().point(s) ;
      const Real x = p(0);
      Number u = c.point_value(0, p);

      if ( x == 0.5 )
      {
        Q(s) += u*0.5  ; /** 1/2 because node will get counted once for each
                           element it belongs to  */
      }
  }
}


}
#endif // CATENARY_SYSTEM_H
