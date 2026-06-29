/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_ASSEMBLY_MATRIXSETUP_H
#define RODIN_PETSC_ASSEMBLY_MATRIXSETUP_H

#include <petscmat.h>

namespace Rodin::PETSc::Assembly
{
  /**
   * @brief Sets up a PETSc matrix for assembly while reusing compatible
   *        existing structure.
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
   * must differ, and @c prepare then performs a full setup.
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
   * When the existing structure is reusable, @c prepare skips the structural
   * calls and only does @c MatZeroEntries + refill. That preserves the nonzero
   * state, so PETSc reports @c SAME_NONZERO_PATTERN and the symbolic
   * factorization is computed once and reused across assemblies.
   */
  class MatrixSetup
  {
    public:
      struct Options
      {
        PetscInt localRows;
        PetscInt localCols;
        PetscInt globalRows;
        PetscInt globalCols;
        MatType type = nullptr;
        bool setFromOptions = false;
        bool keepNonzeroPattern = false;
      };

      explicit MatrixSetup(::Mat matrix)
        : m_matrix(matrix)
      {}

      PetscErrorCode prepare(const Options& options) const
      {
        bool needsSetup = true;
        PetscErrorCode ierr = needsStructuralSetup(options, needsSetup);
        if (ierr)
          return ierr;

        if (needsSetup)
        {
          ierr = MatSetSizes(
              m_matrix,
              options.localRows,
              options.localCols,
              options.globalRows,
              options.globalCols);
          if (ierr)
            return ierr;

          if (options.type)
          {
            ierr = MatSetType(m_matrix, options.type);
            if (ierr)
              return ierr;
          }

          if (options.setFromOptions)
          {
            ierr = MatSetFromOptions(m_matrix);
            if (ierr)
              return ierr;
          }

          if (options.keepNonzeroPattern)
          {
            ierr = MatSetOption(
                m_matrix, MAT_IGNORE_ZERO_ENTRIES, PETSC_FALSE);
            if (ierr)
              return ierr;
          }

          ierr = MatSetUp(m_matrix);
          if (ierr)
            return ierr;

          if (options.keepNonzeroPattern)
          {
            ierr = MatSetOption(
                m_matrix, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
            if (ierr)
              return ierr;
          }
        }

        return MatZeroEntries(m_matrix);
      }

    private:
      PetscErrorCode needsStructuralSetup(
          const Options& options, bool& needsSetup) const
      {
        PetscInt curRows = 0;
        PetscInt curCols = 0;
        PetscErrorCode ierr = MatGetSize(m_matrix, &curRows, &curCols);
        if (ierr)
          return ierr;

        PetscBool assembled = PETSC_FALSE;
        ierr = MatAssembled(m_matrix, &assembled);
        if (ierr)
          return ierr;

        needsSetup = !((assembled == PETSC_TRUE) &&
                       (curRows == options.globalRows) &&
                       (curCols == options.globalCols));
        return PETSC_SUCCESS;
      }

      ::Mat m_matrix;
  };
}

#endif
