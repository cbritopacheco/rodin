/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_ASSEMBLY_MATRIXSETUP_H
#define RODIN_PETSC_ASSEMBLY_MATRIXSETUP_H

#include <petscmat.h>

#include "Rodin/PETSc/Check.h"

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
   * geometry, not the connectivity graph, hence not the pattern.
   *
   * Because a @c LinearSystem is bound to fixed finite element spaces, its
   * global dimensions never change over its lifetime: a different mesh or space
   * means a different @c LinearSystem. PETSc has no in-place resize for an
   * already-assembled matrix (@c MatSetSizes is only valid before the layout is
   * established), so @c prepare does not attempt one. If it ever observes an
   * already-assembled matrix whose dimensions differ from the requested ones —
   * which can only happen if a @c LinearSystem is illegally reused across
   * different spaces — it raises rather than silently corrupting state.
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
        RODIN_PETSC_CHECK_OK(needsStructuralSetup(options, needsSetup));

        if (needsSetup)
        {
          RODIN_PETSC_CHECK_OK(MatSetSizes(
              m_matrix,
              options.localRows,
              options.localCols,
              options.globalRows,
              options.globalCols));

          if (options.type)
            RODIN_PETSC_CHECK_OK(MatSetType(m_matrix, options.type));

          if (options.setFromOptions)
            RODIN_PETSC_CHECK_OK(MatSetFromOptions(m_matrix));

          if (options.keepNonzeroPattern)
            RODIN_PETSC_CHECK_OK(MatSetOption(
                m_matrix, MAT_IGNORE_ZERO_ENTRIES, PETSC_FALSE));

          RODIN_PETSC_CHECK_OK(MatSetUp(m_matrix));

          if (options.keepNonzeroPattern)
            RODIN_PETSC_CHECK_OK(MatSetOption(
                m_matrix, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE));
        }

        return MatZeroEntries(m_matrix);
      }

    private:
      // Decides whether the matrix needs a structural (re)setup. A virgin
      // matrix must be set up. An already-assembled matrix whose dimensions
      // match is reused (only zeroed). An already-assembled matrix whose
      // dimensions differ violates the constant-space contract and cannot be
      // resized in place, so we raise.
      PetscErrorCode needsStructuralSetup(
          const Options& options, bool& needsSetup) const
      {
        PetscInt curRows = 0;
        PetscInt curCols = 0;
        RODIN_PETSC_CHECK_OK(MatGetSize(m_matrix, &curRows, &curCols));

        PetscBool assembled = PETSC_FALSE;
        RODIN_PETSC_CHECK_OK(MatAssembled(m_matrix, &assembled));

        if (assembled == PETSC_TRUE &&
            (curRows != options.globalRows || curCols != options.globalCols))
        {
          Alert::Exception()
            << "MatrixSetup: cannot resize an assembled matrix from "
            << static_cast<long long>(curRows) << "x"
            << static_cast<long long>(curCols) << " to "
            << static_cast<long long>(options.globalRows) << "x"
            << static_cast<long long>(options.globalCols)
            << ". A LinearSystem is bound to fixed finite element spaces; "
            << "use a fresh LinearSystem for a different space or mesh."
            << Alert::Raise;
        }

        needsSetup = (assembled != PETSC_TRUE);
        return PETSC_SUCCESS;
      }

      ::Mat m_matrix;
  };
}

#endif
