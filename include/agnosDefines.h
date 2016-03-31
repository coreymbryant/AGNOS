

#ifndef AGNOS_DEFINES_H
#define AGNOS_DEFINES_H

#ifdef  DEBUG
#define AGNOS_DEBUG 1
#else
#define AGNOS_DEBUG 0
#endif

#include "agnos_config.h"

#include <iostream>
#include <memory>
#include <fstream>
#include <stdio.h>
#include <assert.h>
#include <cstring>
#include <list>
#include <queue>
#include <iterator>

// libmesh includes
#include "libmesh/getpot.h"
#include "libmesh/libmesh.h"
#include "libmesh/parallel.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
/* #include "libmesh/numeric_vector.h" */
/* #include "libmesh/mesh.h" */
/* #include "libmesh/equation_systems.h" */
/* #include "libmesh/nonlinear_solver.h" */
/* #include "libmesh/error_vector.h" */
/* #include "libmesh/adjoint_refinement_estimator.h" */
/* #include "libmesh/mesh_refinement.h" */

typedef libMesh::DenseVector<double> T_S ;
typedef libMesh::DenseVector<double> T_P ;
typedef libMesh::Parallel::Communicator Communicator;
typedef libMesh::PetscVector<double> Vector;
typedef libMesh::PetscMatrix<double> DistMatrix;
typedef libMesh::DenseMatrix<double> LocalMatrix;

/* #include "Parameter.h" */

/* #include "PhysicsModel.h" */
/* #include "PseudoSpectralTensorProduct.h" */
/* #include "PhysicsFunction.h" */
/* #include "PhysicsFunctionTotalError.h" */
/* #include "PhysicsViscousBurgers.h" */
/* #include "PhysicsCatenary.h" */
/* #include "PhysicsCatenaryLibmesh.h" */

#define agnos_assert(asserted)                                        \
  do {                                                                  \
    if (!(asserted)) {                                                  \
      std::cerr << "Assertion `" #asserted "' failed." << std::endl; \
      std::abort() ; \
    } } while(0)

namespace AGNOS{
  /* /** Comparison functor for index sets *1/ */
  /* bool indexSetCompare( */ 
  /*     const std::vector<unsigned int>& a, */ 
  /*     const std::vector<unsigned int>& b */ 
  /*     ) */
  /* { */
  /*   assert( a.size() == b.size() ) ; */
  /*   bool aLessThanb = false; */

  /*   for(unsigned int i=0; i<a.size(); i++) */
  /*     if ( a[i] < b[i] ) */
  /*     { */
  /*       aLessThanb = true ; */
  /*       break; */
  /*     } */
  /*     else if ( a[i] > b[i] ) */
  /*       break; */

  /*   return aLessThanb; */
  /* } */

}


#endif // AGNOS_DEFINES_H
