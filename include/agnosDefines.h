

#ifndef AGNOS_DEFINES_H
#define AGNOS_DEFINES_H

#define DEBUG false

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
typedef libMesh::Parallel::Communicator Communicator;


#include "PhysicsModel.h"
#include "PseudoSpectralTensorProduct.h"
/* #include "PhysicsFunction.h" */
/* #include "PhysicsFunctionTotalError.h" */
#include "PhysicsViscousBurgers.h"
/* #include "PhysicsCatenary.h" */
/* #include "PhysicsCatenaryLibmesh.h" */


namespace AGNOS{
}


#endif // AGNOS_DEFINES_H
