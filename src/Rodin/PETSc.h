#ifndef RODIN_PETSC_H
#define RODIN_PETSC_H

#include "PETSc/ForwardDecls.h"

#include "PETSc/Math/Vector.h"
#include "PETSc/Math/Matrix.h"
#include "PETSc/Math/LinearSystem.h"

#include "PETSc/Assembly/Sequential.h"
#include "PETSc/Assembly/MPI.h"

#include "PETSc/Solver/CG.h"

#include "PETSc/FormLanguage/Traits.h"

#ifdef RODIN_USE_OPENMP
#include "PETSc/Assembly/OpenMP.h"
#endif

#endif
