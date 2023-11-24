/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ASSEMBLY_SERIAL_H
#define RODIN_ASSEMBLY_SERIAL_H

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "AssemblyBase.h"

namespace Rodin::Assembly
{
  template <>
  class Serial<Variational::BilinearFormBase<std::vector<Eigen::Triplet<Scalar>>>>
    : public AssemblyBase<Variational::BilinearFormBase<std::vector<Eigen::Triplet<Scalar>>>>
  {
    /**
     * @internal
     */
    static void add(
        std::vector<Eigen::Triplet<Scalar>>& out, const Math::Matrix& in,
        const IndexArray& rows, const IndexArray& cols);

    public:
      using Parent =
        AssemblyBase<Variational::BilinearFormBase<std::vector<Eigen::Triplet<Scalar>>>>;
      using OperatorType = std::vector<Eigen::Triplet<Scalar>>;

      Serial() = default;

      Serial(const Serial& other)
        : Parent(other)
      {}

      Serial(Serial&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Executes the assembly and returns the linear operator
       * associated to the bilinear form.
       */
      OperatorType execute(const BilinearAssemblyInput& input) const override;

      Serial* copy() const noexcept override
      {
        return new Serial(*this);
      }
  };

  /**
   * @brief Serial assembly of the Math::SparseMatrix associated to a
   * BilinearFormBase object.
   */
  template <>
  class Serial<Variational::BilinearFormBase<Math::SparseMatrix>>
    : public AssemblyBase<Variational::BilinearFormBase<Math::SparseMatrix>>
  {
    public:
      using Parent = AssemblyBase<Variational::BilinearFormBase<Math::SparseMatrix>>;
      using OperatorType = Math::SparseMatrix;

      Serial() = default;

      Serial(const Serial& other)
        : Parent(other)
      {}

      Serial(Serial&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Executes the assembly and returns the linear operator
       * associated to the bilinear form.
       */
      OperatorType execute(const BilinearAssemblyInput& input) const override;

      Serial* copy() const noexcept override
      {
        return new Serial(*this);
      }
  };

  /**
   * @brief Serial assembly of the Math::SparseMatrix associated to a
   * BilinearFormBase object.
   */
  template <>
  class Serial<Variational::BilinearFormBase<Math::Matrix>>
    : public AssemblyBase<Variational::BilinearFormBase<Math::Matrix>>
  {
    static void add(
        Math::Matrix& out, const Math::Matrix& in,
        const IndexArray& rows, const IndexArray& cols);

    public:
      using Parent = AssemblyBase<Variational::BilinearFormBase<Math::Matrix>>;
      using OperatorType = Math::Matrix;

      Serial() = default;

      Serial(const Serial& other)
        : Parent(other)
      {}

      Serial(Serial&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Executes the assembly and returns the linear operator
       * associated to the bilinear form.
       */
      OperatorType execute(const BilinearAssemblyInput& input) const override;

      Serial* copy() const noexcept override
      {
        return new Serial(*this);
      }
  };

  /**
   * @brief %Serial assembly of the Math::Vector associated to a LinearFormBase
   * object.
   */
  template <>
  class Serial<Variational::LinearFormBase<Math::Vector>>
    : public AssemblyBase<Variational::LinearFormBase<Math::Vector>>
  {

    static void add(Math::Vector& out, const Math::Vector& in, const IndexArray& s);

    public:
      using Parent = AssemblyBase<Variational::LinearFormBase<Math::Vector>>;
      using VectorType = Math::Vector;

      Serial() = default;

      Serial(const Serial& other)
        : Parent(other)
      {}

      Serial(Serial&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Executes the assembly and returns the vector associated to the
       * linear form.
       */
      VectorType execute(const Input& input) const override;

      Serial* copy() const noexcept override
      {
        return new Serial(*this);
      }
  };
}

#endif

