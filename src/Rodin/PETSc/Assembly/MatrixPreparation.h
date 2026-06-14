/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_ASSEMBLY_MATRIXPREPARATION_H
#define RODIN_PETSC_ASSEMBLY_MATRIXPREPARATION_H

#include <petscmat.h>

namespace Rodin::PETSc::Assembly
{
  /**
   * @brief Decides whether a matrix needs its structure (sizes/type/layout)
   *        set up before assembly, or whether the existing structure can be
   *        re-used.
   *
   * ## Contract
   *
   * The caller guarantees that, between successive assemblies of the *same*
   * matrix object, the sparsity pattern is unchanged. Concretely this means the
   * trial and test finite element spaces and the mesh connectivity are fixed.
   * Vertex coordinates may move (as in an ALE formulation) because that changes
   * geometry, not the connectivity graph, hence not the pattern. If the
   * contract is ever violated (adaptive remeshing, a changed constraint set that
   * alters the eliminated structure, a resized problem) the global dimensions
   * must differ, and this function then reports that a full setup is needed.
   *
   * ## Why this matters
   *
   * The structural calls @c MatSetSizes / @c MatSetType / @c MatSetFromOptions /
   * @c MatSetUp re-establish the matrix internals and bump PETSc's nonzero-state
   * counter. Preconditioners compare that counter
   * (@c pc->matnonzerostate != Amat->nonzerostate) to decide whether to redo the
   * symbolic factorization. Running the structural calls on every assembly
   * forces external direct solvers (e.g. MUMPS) to recompute their
   * ordering/analysis each time — which, inside a SNES/Newton loop, repeats the
   * expensive symbolic factorization on every residual/Jacobian evaluation.
   *
   * When this function reports that no setup is needed, the caller should skip
   * the structural calls and only @c MatZeroEntries + refill. That preserves the
   * nonzero state, so PETSc reports @c SAME_NONZERO_PATTERN and the symbolic
   * factorization is computed once and reused across assemblies.
   *
   * @param[in] A            Matrix to inspect.
   * @param[in] globalRows   Expected global number of rows.
   * @param[in] globalCols   Expected global number of columns.
   * @param[out] needsSetup  Set to @c true when the structural setup must be
   *                         performed (first assembly or dimensions changed),
   *                         @c false when the existing structure can be re-used.
   * @returns The first non-zero PETSc error code encountered, else
   *          @c PETSC_SUCCESS.
   */
  inline PetscErrorCode matrixNeedsStructuralSetup(
      ::Mat A, PetscInt globalRows, PetscInt globalCols, bool& needsSetup)
  {
    PetscInt curRows = 0;
    PetscInt curCols = 0;
    PetscErrorCode ierr = MatGetSize(A, &curRows, &curCols);
    if (ierr)
      return ierr;

    PetscBool assembled = PETSC_FALSE;
    ierr = MatAssembled(A, &assembled);
    if (ierr)
      return ierr;

    needsSetup = !((assembled == PETSC_TRUE) && (curRows == globalRows) &&
                   (curCols == globalCols));
    return PETSC_SUCCESS;
  }
}

#endif
