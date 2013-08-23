
#ifndef AGNOS_DEFINES_H
#define AGNOS_DEFINES_H

#define AGNOS_DEBUG 1

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <assert.h>
#include <cstring>

// libmesh includes
#include "libmesh/getpot.h"
#include "libmesh/libmesh.h"
#include "libmesh/parallel.h"
#include "libmesh/dense_vector.h"
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

#include "Parameter.h"



namespace AGNOS{
}


#endif // AGNOS_DEFINES_H
