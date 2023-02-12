/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_TRANSPOSE_H
#define RODIN_VARIATIONAL_TRANSPOSE_H

#include "ShapeFunction.h"
#include "MatrixFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup TransposeSpecializations Transpose Template Specializations
   * @brief Template specializations of the Transpose class.
   * @see Transpose
   */

  /**
   * @brief Transpose of a FunctionBase object.
   * @ingroup TransposeSpecializations
   */
  template <>
  class Transpose<FunctionBase> : public FunctionBase
  {
    public:
      /**
       * @brief Constructs the Transpose matrix of the given matrix.
       */
      Transpose(const FunctionBase& m)
        : m_matrix(m.copy())
      {}

      Transpose(const Transpose& other)
        :  FunctionBase(other),
          m_matrix(other.m_matrix->copy())
      {}

      Transpose(Transpose&& other)
        : FunctionBase(std::move(other)),
          m_matrix(std::move(other.m_matrix))
      {}

      RangeShape getRangeShape() const override
      {
        return m_matrix->getRangeShape().transpose();
      }

      FunctionValue getValue(const Geometry::Point& p) const override
      {
        Math::Matrix res = m_matrix->getValue(p).matrix().transpose();
        return res;
      }

      Transpose* copy() const noexcept override
      {
        return new Transpose(*this);
      }

    private:
      std::unique_ptr<FunctionBase> m_matrix;
  };
  Transpose(const FunctionBase&) -> Transpose<MatrixFunctionBase>;

  /**
   * @brief Transpose of a ShapeFunctionBase object.
   * @ingroup TransposeSpecializations
   */
  template <ShapeFunctionSpaceType Space>
  class Transpose<ShapeFunctionBase<Space>> : public ShapeFunctionBase<Space>
  {
    public:
      Transpose(const ShapeFunctionBase<Space>& op)
        : m_shape(op.copy())
      {}

      Transpose(const Transpose& other)
        : m_shape(other.m_shape->copy())
      {}

      Transpose(Transpose&& other)
        : m_shape(std::move(other.m_shape))
      {}

      const ShapeFunctionBase<Space>& getLeaf() const override
      {
        return m_shape->getLeaf();
      }

      int getRows() const override
      {
        return m_shape->getColumns();
      }

      int getColumns() const override
      {
        return m_shape->getRows();
      }

      int getDOFs(const Geometry::Simplex& element) const override
      {
        return m_shape->getDOFs(element);
      }

      TensorBasis getOperator(
          ShapeComputator& compute, const Geometry::Point& p) const override
      {
        const auto& fe = getFiniteElementSpace().getFiniteElement(p.getSimplex());
        return m_shape->getOperator(compute, p).shuffle(std::array<int, 3>{0, 2, 1});
      }

      FiniteElementSpaceBase& getFiniteElementSpace() override
      {
        return m_shape->getFiniteElementSpace();
      }

      const FiniteElementSpaceBase& getFiniteElementSpace() const override
      {
        return m_shape->getFiniteElementSpace();
      }

      Transpose* copy() const noexcept override
      {
        return new Transpose(*this);
      }

    private:
      std::unique_ptr<ShapeFunctionBase<Space>> m_shape;
  };
}

#endif
