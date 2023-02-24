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
        : m_operand(m.copy())
      {}

      constexpr
      Transpose(const Transpose& other)
        : Parent(other),
          m_operand(other.m_matrix->copy())
      {}

      constexpr
      Transpose(Transpose&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_matrix))
      {}

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        return m_operand->getRangeShape().transpose();
      }

      inline
      constexpr
      const auto& getOperand() const
      {
        assert(m_operand);
        return *m_operand;
      }

      inline
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        using OperandRange = typename FormLanguage::Traits<Operand>::RangeType;
        static_assert(std::is_same_v<OperandRange, Math::Vector> || std::is_same_v<OperandRange, Math::Matrix>);
        return this->object(getOperand().getValue(p)).transpose();
      }

      inline
      Transpose* copy() const noexcept override
      {
        return new Transpose(*this);
      }

    private:
      std::unique_ptr<Operand> m_operand;
  };

  template <class NestedDerived>
  Transpose(const FunctionBase<NestedDerived>&) -> Transpose<FunctionBase<NestedDerived>>;

  /**
   * @brief Transpose of a ShapeFunctionBase object.
   * @ingroup TransposeSpecializations
   */
  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  class Transpose<ShapeFunctionBase<NestedDerived, FES, Space>> final
    : public ShapeFunctionBase<Transpose<ShapeFunctionBase<NestedDerived, FES, Space>>, FES, Space>
  {
    public:
      using Operand = ShapeFunctionBase<NestedDerived, FES, Space>;
      using Parent = ShapeFunctionBase<Transpose<Operand>, FES, Space>;

      constexpr
      Transpose(const Operand& op)
        : m_operand(op.copy())
      {}

      constexpr
      Transpose(const Transpose& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      constexpr
      Transpose(Transpose&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}

      inline
      constexpr
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        return getOperand().getRangeShape().transpose();
      }

      inline
      constexpr
      size_t getDOFs(const Geometry::Simplex& simplex) const
      {
        return getOperand().getDOFs(simplex);
      }

      inline
      constexpr
      const Operand& getOperand() const
      {
        assert(m_operand);
        return *m_operand;
      }

      inline
      constexpr
      auto getTensorBasis(const Geometry::Point& p) const
      {
        using OperandRange = typename FormLanguage::Traits<Operand>::RangeType;
        static_assert(std::is_same_v<OperandRange, Math::Vector> || std::is_same_v<OperandRange, Math::Matrix>);
        const auto& op = this->object(getOperand().getTensorBasis(p));
        return TensorBasis(getDOFs(p.getSimplex()), [&](size_t i){ return op(i).transpose(); });
      }

      inline
      constexpr
      const FES& getFiniteElementSpace() const
      {
        return m_operand.getFiniteElementSpace();
      }

    private:
      std::unique_ptr<Operand> m_operand;
  };

  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  Transpose(const ShapeFunctionBase<NestedDerived, FES, Space>&)
    -> Transpose<ShapeFunctionBase<NestedDerived, FES, Space>>;
}

#endif
