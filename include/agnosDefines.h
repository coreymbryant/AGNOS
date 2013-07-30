
#ifndef AGNOS_DEFINES_H
#define AGNOS_DEFINES_H

#define DEBUG 1

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <assert.h>
#include <cstring>
/* #include <GetPot> */

#include "libmesh/getpot.h"
#include "libmesh/libmesh.h"
#include "libmesh/parallel.h"
#include "libmesh/dense_vector.h"
#include "libmesh/numeric_vector.h"
// libmesh includes
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/error_vector.h"
#include "libmesh/adjoint_refinement_estimator.h"
#include "libmesh/mesh_refinement.h"
typedef libMesh::DenseVector<double> T_S ;
typedef libMesh::DenseVector<double> T_P ;
typedef libMesh::Parallel::Communicator Communicator;

#include "Parameter.h"

/* #include "PhysicsModel.h" */
/* #include "PseudoSpectralTensorProduct.h" */
/* /1* #include "PhysicsFunction.h" *1/ */
/* /1* #include "PhysicsFunctionTotalError.h" *1/ */
/* #include "PhysicsViscousBurgers.h" */
/* /1* #include "PhysicsCatenary.h" *1/ */
/* /1* #include "PhysicsCatenaryLibmesh.h" *1/ */



namespace AGNOS{
}


#endif // AGNOS_DEFINES_H
