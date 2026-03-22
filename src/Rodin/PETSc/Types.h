#ifndef RODIN_PETSC_TYPES_H
#define RODIN_PETSC_TYPES_H

/**
 * @file
 * @brief PETSc scalar, integer, and floating-point type aliases.
 *
 * These aliases expose the PETSc-configured numeric types so that
 * Rodin code can be written independently of PETSc's scalar mode
 * (real or complex).
 */

#include <petsc.h>

namespace Rodin::PETSc
{
  /// @brief Scalar type used by PETSc vectors and matrices.
  using Scalar = PetscScalar;

  /// @brief Integer index type used by PETSc row/column indices.
  using Integer = PetscInt;

  /// @brief Real-valued floating-point type (tolerances, norms, etc.).
  using Real = PetscReal;

  /// @brief Complex scalar type (available when PETSc is built with complex support).
  using Complex = PetscComplex;
}

#endif
