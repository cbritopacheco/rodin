/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ASSEMBLY_NATIVE_H
#define RODIN_ASSEMBLY_NATIVE_H

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "AssemblyBase.h"

namespace Rodin::Variational::Assembly
{
  /**
   * @brief Native assembly of the Math::SparseMatrix associated to a
   * BilinearFormBase object.
   */
  template <>
  class Native<BilinearFormBase<Math::SparseMatrix>>
    : public AssemblyBase<BilinearFormBase<Math::SparseMatrix>>
  {
    /**
     * @internal
     */
    static void add(Math::SparseMatrix& out, const Math::Matrix& in, const IndexArray& rows, const IndexArray& cols);

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

      /**
       * @brief Executes the assembly and returns the linear operator
       * associated to the bilinear form.
       */
      OperatorType execute(const Input& input) const override;

      Native* copy() const noexcept override
      {
        return new Native(*this);
      }
  };

  /**
   * @brief %Native assembly of the Math::Vector associated to a LinearFormBase
   * object.
   */
  template <>
  class Native<LinearFormBase<Math::Vector>>
    : public AssemblyBase<LinearFormBase<Math::Vector>>
  {

    static void add(Math::Vector& out, const Math::Vector& in, const IndexArray& s);

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

      /**
       * @brief Executes the assembly and returns the vector associated to the
       * linear form.
       */
      VectorType execute(const Input& input) const override;

      Native* copy() const noexcept override
      {
        return new Native(*this);
      }
  };
}

#endif

