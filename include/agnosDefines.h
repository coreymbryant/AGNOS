
#ifndef AGNOS_DEFINES_H
#define AGNOS_DEFINES_H

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
typedef libMesh::DenseVector<double> T_P ;
typedef libMesh::DenseVector<double> T_S ;
typedef libMesh::Parallel::Communicator Communicator;


#include "PseudoSpectralTensorProduct.h"
#include "PhysicsFunction.h"
#include "PhysicsFunctionTotalError.h"
#include "PhysicsViscousBurgers.h"
#include "PhysicsCatenary.h"
#include "PhysicsCatenaryLibmesh.h"


namespace AGNOS{
}


#endif // AGNOS_DEFINES_H
