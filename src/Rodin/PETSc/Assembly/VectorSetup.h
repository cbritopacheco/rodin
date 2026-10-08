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
   * A @c LinearSystem is associated with fixed finite element spaces. Its
   * right-hand side @f$ \mathbf{b} @f$ and solution @f$ \mathbf{x} @f$ therefore
   * retain their global sizes throughout its lifetime. The first call
   * establishes the vector layout; subsequent calls reuse that layout.
   *
   * @note PETSc does not support resizing a vector after its layout has been
   * established. Reuse requires the existing global size to match the requested
   * size; this precondition is checked by a debug assertion. A change of finite
   * element space or mesh requires a new @c LinearSystem.
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
       * Vectors without a PETSc type receive their initial sizes, optional type,
       * and optional command line options. Reused vectors retain their structure
       * and are zeroed only when requested by @ref Options::zeroOnReuse.
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
        (void)ierr;

        if (needsSetup)
        {
          ierr = VecSetSizes(m_vector, options.localSize, options.globalSize);
          assert(ierr == PETSC_SUCCESS);
          (void)ierr;

          if (options.type)
          {
            ierr = VecSetType(m_vector, options.type);
            assert(ierr == PETSC_SUCCESS);
            (void)ierr;
          }

          if (options.setFromOptions)
          {
            ierr = VecSetFromOptions(m_vector);
            assert(ierr == PETSC_SUCCESS);
            (void)ierr;
          }
        }

        if (needsSetup || options.zeroOnReuse)
          return VecZeroEntries(m_vector);

        return PETSC_SUCCESS;
      }

    private:
      /**
       * @brief Determines whether the vector needs structural setup.
       *
       * @note Structural setup is required only when the vector has no PETSc
       * type. Otherwise, the existing global size must equal the requested size.
       * A debug assertion checks this precondition; the vector is not resized.
       *
       * @param[in] options Requested vector layout.
       * @param[out] needsSetup True when the vector has no PETSc type yet.
       * @returns @c PETSC_SUCCESS after checking the existing layout.
       */
      PetscErrorCode needsStructuralSetup(const Options& options, bool& needsSetup) const
      {
        VecType curType = nullptr;
        PetscErrorCode ierr = VecGetType(m_vector, &curType);
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;

        const bool hasType = (curType != nullptr);
        if (hasType)
        {
          PetscInt curSize = 0;
          ierr = VecGetSize(m_vector, &curSize);
          assert(ierr == PETSC_SUCCESS);
          (void)ierr;
          assert(curSize == options.globalSize &&
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
