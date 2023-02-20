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
  template <class NestedDerived>
  class Transpose<FunctionBase<NestedDerived>> final
    : public FunctionBase<Transpose<FunctionBase<NestedDerived>>>
  {
    public:
      using Operand = FunctionBase<NestedDerived>;
      using Parent = FunctionBase<Operand>;

      /**
       * @brief Constructs the Transpose matrix of the given matrix.
       */
      constexpr
      Transpose(const Operand& m)
        : m_matrix(m)
      {}

      constexpr
      Transpose(const Transpose& other)
        : Parent(other),
          m_matrix(other.m_matrix)
      {}

      constexpr
      Transpose(Transpose&& other)
        : Parent(std::move(other)),
          m_matrix(std::move(other.m_matrix))
      {}

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        return m_matrix->getRangeShape().transpose();
      }

      inline
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        using OperandRange =
          FormLanguage::RangeOf<typename FormLanguage::Traits<Operand>::ResultType>;
        static_assert(std::is_same_v<OperandRange, Math::Matrix>);
        return m_matrix.getValue(p).transpose();
      }

    private:
      Operand m_matrix;
  };

  template <class NestedDerived>
  Transpose(const FunctionBase<NestedDerived>&) -> Transpose<FunctionBase<NestedDerived>>;

  /**
   * @brief Transpose of a ShapeFunctionBase object.
   * @ingroup TransposeSpecializations
   */
  template <class NestedDerived, ShapeFunctionSpaceType Space>
  class Transpose<ShapeFunctionBase<NestedDerived, Space>> final
    : public ShapeFunctionBase<Transpose<ShapeFunctionBase<NestedDerived, Space>>, Space>
  {
    public:
      using Operand = ShapeFunctionBase<NestedDerived, Space>;
      using Parent = ShapeFunctionBase<Transpose<Operand>, Space>;

      constexpr
      Transpose(const Operand& op)
        : m_shape(op)
      {}

      constexpr
      Transpose(const Transpose& other)
        : Parent(other),
          m_shape(other.m_shape)
      {}

      constexpr
      Transpose(Transpose&& other)
        : Parent(std::move(other)),
          m_shape(std::move(other.m_shape))
      {}

      inline
      constexpr
      const auto& getLeaf() const
      {
        return m_shape.getLeaf();
      }

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        return m_shape.getRangeShape().transpose();
      }

      inline
      constexpr
      size_t getDOFs(const Geometry::Simplex& simplex) const
      {
        return m_shape.getDOFs(simplex);
      }

      inline
      constexpr
      auto getOperator(ShapeComputator& compute, const Geometry::Point& p) const
      {
        return m_shape.getOperator(compute, p).shuffle(std::array<int, 3>{0, 2, 1});
      }

      inline
      constexpr
      auto& getFiniteElementSpace()
      {
        return m_shape.getFiniteElementSpace();
      }

      inline
      constexpr
      const auto& getFiniteElementSpace() const
      {
        return m_shape.getFiniteElementSpace();
      }

    private:
      Operand m_shape;
  };

  template <class NestedDerived, ShapeFunctionSpaceType Space>
  Transpose(const ShapeFunctionBase<NestedDerived, Space>&)
    -> Transpose<ShapeFunctionBase<NestedDerived, Space>>;
}

#endif
