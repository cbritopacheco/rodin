/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ASSEMBLY_NATIVE_H
#define RODIN_ASSEMBLY_NATIVE_H

#include <variant>
#include <ostream>

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "AssemblyBase.h"

namespace Rodin::Variational::Assembly
{
  template <>
  class Native<BilinearFormBase<Math::SparseMatrix>>
    : public AssemblyBase<BilinearFormBase<Math::SparseMatrix>>
  {
    public:
      using Parent = AssemblyBase<BilinearFormBase<Math::SparseMatrix>>;
      using OperatorType = Math::SparseMatrix;

      Native() = default;

      Native(const Native& other)
        : Parent(other)
      {}

      Native(Native&& other)
        : Parent(std::move(other))
      {}

      OperatorType execute(const Input& input) const override;

      Native* copy() const noexcept override
      {
        return new Native(*this);
      }
  };

  template <>
  class Native<LinearFormBase<Math::Vector>>
    : public AssemblyBase<LinearFormBase<Math::Vector>>
  {
    public:
      using Parent = AssemblyBase<LinearFormBase<Math::Vector>>;
      using VectorType = Math::Vector;

      Native() = default;

      Native(const Native& other)
        : Parent(other)
      {}

      Native(Native&& other)
        : Parent(std::move(other))
      {}

      VectorType execute(const Input& input) const override;

      Native* copy() const noexcept override
      {
        return new Native(*this);
      }
  };
}

#endif

