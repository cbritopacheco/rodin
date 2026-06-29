/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_ASSEMBLY_VECTORSETUP_H
#define RODIN_PETSC_ASSEMBLY_VECTORSETUP_H

#include <petscvec.h>

#include "Rodin/Alert/Exception.h"
#include "Rodin/Alert/Raise.h"
#include "Rodin/PETSc/Check.h"

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
   * requested one — only possible if a @c LinearSystem is illegally reused
   * across different spaces — it raises rather than silently corrupting state.
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
      struct Options
      {
        PetscInt localSize;
        PetscInt globalSize;
        VecType type = nullptr;
        bool setFromOptions = false;
        bool zeroOnReuse = true;
      };

      explicit VectorSetup(::Vec vector)
        : m_vector(vector)
      {}

      PetscErrorCode prepare(const Options& options) const
      {
        bool needsSetup = true;
        RODIN_PETSC_CHECK_OK(needsStructuralSetup(options, needsSetup));

        if (needsSetup)
        {
          RODIN_PETSC_CHECK_OK(VecSetSizes(
              m_vector, options.localSize, options.globalSize));

          if (options.type)
            RODIN_PETSC_CHECK_OK(VecSetType(m_vector, options.type));

          if (options.setFromOptions)
            RODIN_PETSC_CHECK_OK(VecSetFromOptions(m_vector));
        }

        if (needsSetup || options.zeroOnReuse)
          return VecZeroEntries(m_vector);

        return PETSC_SUCCESS;
      }

    private:
      // A vector with no type yet is virgin and must be laid out. A typed
      // vector whose size matches is reused. A typed vector whose size differs
      // violates the constant-space contract and cannot be resized in place, so
      // we raise.
      PetscErrorCode needsStructuralSetup(
          const Options& options, bool& needsSetup) const
      {
        VecType curType = nullptr;
        RODIN_PETSC_CHECK_OK(VecGetType(m_vector, &curType));

        const bool hasType = (curType != nullptr);
        if (hasType)
        {
          PetscInt curSize = 0;
          RODIN_PETSC_CHECK_OK(VecGetSize(m_vector, &curSize));
          if (curSize != options.globalSize)
          {
            Alert::Exception()
              << "VectorSetup: cannot resize a vector from "
              << static_cast<long long>(curSize) << " to "
              << static_cast<long long>(options.globalSize)
              << ". A LinearSystem is bound to fixed finite element spaces; "
              << "use a fresh LinearSystem for a different space or mesh."
              << Alert::Raise;
          }
        }

        needsSetup = !hasType;
        return PETSC_SUCCESS;
      }

      ::Vec m_vector;
  };
}

#endif
