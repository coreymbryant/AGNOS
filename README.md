AGNOS [![Build Status](https://travis-ci.org/coreymbryant/AGNOS.svg?branch=master)](https://travis-ci.org/coreymbryant/AGNOS)
https://travis-ci.org/coreymbryant/libmesh.svg?branch=build
=====

AGNOS supports non-intrusive construction of surrogate models (response surface
approximations) for uncertainty propagation. It is agnostic to the
discretization method being used for the physical problem. It requires
the user to define a PhysicsModel class with a compute function providing a
std::map of name-vector pairs. 

The adaptive framework utilizes adjoint information, therfore To utilize the
full functionality of adaptivity and error estimation the users physics model
class should compute:"
* primal
* adjoint
* qoi
as well as a refine function to adapt the physical discretization is ideal.
Alternatively, one can supply a primal solution, qoi,  and their own
errorEstimate to be used in the adaptive procedure. 

AGNOS has a PhysicsModelLibmesh class that supports the
[libMesh](https://github.com/libMesh/libmesh.git) finite element library.
Existing [libMesh](https://github.com/libMesh/libmesh.git) simulations can be
wrapped in this model class with limited modification.  Limited functionality is
also implemented for [GRINS](http://grinsfem.github.io): A C++ Multiphysics
Finite Element Package based on the libMesh Finite Element Library

AGNOS has the capability to write/read HDF5 files for the surrogate model. Once
a response surface has been created, it can be saved in HDF5 format and it can
be used in external (or subsequent) codes to construct a SurrogateEvaluator
provided by the library. 

A number of examples are provided of varying complexity. Development efforts of
the library have been put on hold, but a number of future improvements are
planned, or are already in some stage of implementation, including:
* the addition of Gaussian random variables as parameters
* additional SurrogateModel classes (MonteCarlo, SparseGrids, Collocation)

Dependicies: (recommend being built in this order)
* gcc
* mpi (open-mpi or mpich)
* mkl (optional)
* petsc
* trilinos
* hdf5
* vtk
* tbb
* boost 
* cpp-unit (for tests)
