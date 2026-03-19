#ifndef RODIN_PETSC_TYPES_H
#define RODIN_PETSC_TYPES_H

/**
 * @file
 * @brief PETSc scalar and index aliases used by Rodin.
 */

#include <petsc.h>

namespace Rodin::PETSc
{
  /** @brief Scalar type used by PETSc vectors and matrices. */
  using Scalar = PetscScalar;

  /** @brief Integer index type used by PETSc APIs. */
  using Integer = PetscInt;

  /** @brief Real-valued floating-point type used by PETSc tolerances and norms. */
  using Real = PetscReal;

  /** @brief Complex scalar type exposed by PETSc builds with complex arithmetic. */
  using Complex = PetscComplex;
}

#endif
