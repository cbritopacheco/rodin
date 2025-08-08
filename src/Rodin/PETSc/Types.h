#ifndef RODIN_PETSC_TYPES_H
#define RODIN_PETSC_TYPES_H

#include <petsc.h>

namespace Rodin::PETSc
{
  using Scalar = PetscScalar;

  using Integer = PetscInt;

  using Real = PetscReal;

  using Complex = PetscComplex;
}

#endif
