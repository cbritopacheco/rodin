#ifndef RODIN_PETSC_ASSEMBLY_H
#define RODIN_PETSC_ASSEMBLY_H

#include "Assembly/Sequential.h"

#ifdef RODIN_USE_MPI
#include "Assembly/MPI.h"
#endif

#ifdef RODIN_USE_OPENMP
#include "Assembly/OpenMP.h"
#endif

#endif
