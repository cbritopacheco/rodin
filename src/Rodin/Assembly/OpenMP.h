/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef RODIN_ASSEMBLY_OPENMP_H
#define RODIN_ASSEMBLY_OPENMP_H

#include <variant>
#include <ostream>

#include "AssemblyBase.h"

#ifdef RODIN_USE_OPENMP
namespace Rodin::Variational::Assembly
{
  template <>
  class OpenMP<BilinearFormBase<mfem::SparseMatrix>>
    : public AssemblyBase<BilinearFormBase<mfem::SparseMatrix>>
  {
    public:
      using Parent = AssemblyBase<BilinearFormBase<mfem::SparseMatrix>>;
      using OperatorType = mfem::SparseMatrix;

      OpenMP() = default;

      OpenMP(const OpenMP& other)
        : Parent(other)
      {}

      OpenMP(OpenMP&& other)
        : Parent(std::move(other))
      {}

      OperatorType execute(const Input& input) const override;

      OpenMP* copy() const noexcept override
      {
        return new OpenMP(*this);
      }
  };
}

#endif // #ifdef RODIN_USE_OPENMP
#endif
