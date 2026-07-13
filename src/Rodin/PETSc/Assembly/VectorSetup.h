/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_ASSEMBLY_VECTORSETUP_H
#define RODIN_PETSC_ASSEMBLY_VECTORSETUP_H

#include <petscvec.h>

#include <cassert>

namespace Rodin::PETSc::Assembly
{
  /**
   * @brief Sets up a PETSc vector for assembly while reusing compatible
   *        existing structure.
   *
   * The @c Vec counterpart of @ref MatrixSetup. A @c LinearSystem is bound to
   * fixed finite element spaces, so its right-hand side @f$ \mathbf{b} @f$ and
   * solution @f$ \mathbf{x} @f$ keep their sizes for the whole lifetime. The
   * first call lays the vector out (sizes, type, options); subsequent calls
   * reuse it. PETSc has no in-place resize for a vector whose layout is set, so
   * if @c prepare ever observes a typed vector whose size differs from the
   * requested one, which is only possible if a @c LinearSystem is illegally
   * reused across different spaces, it fails the debug assertion.
   *
   * ## Zeroing
   *
   * The right-hand side and residual vectors are rebuilt from scratch on every
   * assembly and must therefore be zeroed each time (@c zeroOnReuse = true, the
   * default). The solution vector must *not* be zeroed on reuse: it carries the
   * previous iterate, which iterative/Newton solvers consume as the initial
   * guess. It is still zeroed once, when first laid out, to provide a defined
   * starting point.
   */
  class VectorSetup
  {
    public:
      /**
       * @brief Requested PETSc vector layout and setup policy.
       *
       * Sizes and optional type/options are applied only when the vector has no
       * PETSc type yet. Reused vectors may still be zeroed, depending on
       * @ref zeroOnReuse.
       */
      struct Options
      {
        /// @brief Local vector size owned by this MPI rank, or @c PETSC_DECIDE.
        PetscInt localSize;
        /// @brief Global vector size.
        PetscInt globalSize;
        /// @brief Optional PETSc vector type to set during initial setup.
        VecType type = nullptr;
        /// @brief Whether to call @c VecSetFromOptions during initial setup.
        bool setFromOptions = false;
        /// @brief Whether compatible reused vectors are zeroed before assembly.
        bool zeroOnReuse = true;
      };

      /**
       * @brief Wraps an existing PETSc vector handle.
       * @param[in] vector Vector to prepare before assembly or solve reuse.
       */
      explicit VectorSetup(::Vec vector)
        : m_vector(vector)
      {}

      /**
       * @brief Prepares the vector for use by assembly or a solver.
       *
       * Virgin vectors receive their sizes, optional type, and optional command
       * line options. Reused vectors keep their structure and are zeroed only
       * when requested by @ref Options::zeroOnReuse.
       *
       * @param[in] options Requested layout and reuse policy.
       * @returns PETSc error code from zeroing, or @c PETSC_SUCCESS when no
       *          zeroing is requested.
       */
      PetscErrorCode prepare(const Options& options) const
      {
        bool needsSetup = true;
        PetscErrorCode ierr = needsStructuralSetup(options, needsSetup);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;

        if (needsSetup)
        {
          ierr = VecSetSizes(m_vector, options.localSize, options.globalSize);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;

          if (options.type)
          {
            ierr = VecSetType(m_vector, options.type);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
          }

          if (options.setFromOptions)
          {
            ierr = VecSetFromOptions(m_vector);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
          }
        }

        if (needsSetup || options.zeroOnReuse)
          return VecZeroEntries(m_vector);

        return PETSC_SUCCESS;
      }

    private:
      // A vector with no type yet is virgin and must be laid out. A typed
      // vector whose size matches is reused. A typed vector whose size differs
      // violates the constant-space contract and cannot be resized in place, so
      // it fails the debug assertion.
      PetscErrorCode needsStructuralSetup(
          const Options& options, bool& needsSetup) const
      {
        VecType curType = nullptr;
        PetscErrorCode ierr = VecGetType(m_vector, &curType);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;

        const bool hasType = (curType != nullptr);
        if (hasType)
        {
          PetscInt curSize = 0;
          ierr = VecGetSize(m_vector, &curSize);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
          assert(
              curSize == options.globalSize &&
              "VectorSetup cannot resize an assembled vector; use a fresh "
              "LinearSystem for a different space or mesh.");
        }

        needsSetup = !hasType;
        return PETSC_SUCCESS;
      }

      ::Vec m_vector;
  };
}

#endif
