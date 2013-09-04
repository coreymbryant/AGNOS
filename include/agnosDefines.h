
#ifndef AGNOS_DEFINES_H
#define AGNOS_DEFINES_H

#ifndef  NDEBUG
#define AGNOS_DEBUG 0
#else
#define AGNOS_DEBUG 1
#endif

#include "agnos_config.h"

#include <iostream>
#include <memory>
#include <fstream>
#include <stdio.h>
#include <assert.h>
#include <cstring>

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
typedef libMesh::NumericVector<double> Vector;
typedef libMesh::PetscMatrix<double> DistMatrix;
typedef libMesh::DenseMatrix<double> LocalMatrix;

/* #include "Parameter.h" */



namespace AGNOS{
}


#endif // AGNOS_DEFINES_H
