AGNOS
=====

AGNOS supports non-intrusive construction of surrogate models for uncertainty
propagation. It is agnostic to the discretization method being used, but
requires the user to provide a physics class that supports a number of
different function. 

(Agnostos is the Greek word for unknown.)

Build order:
* gcc
* mpi (open-mpi or mpich)
* (mkl)
* petsc
* trilinos
* vtk
* tbb
* boost 
* cpp-unit (for tests only)