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
      using OperandType = FunctionBase<NestedDerived>;

      using Parent = FunctionBase<Transpose<OperandType>>;

      using OperandRangeType = typename FormLanguage::Traits<OperandType>::RangeType;

      /**
       * @brief Constructs the Transpose matrix of the given matrix.
       */
      constexpr
      Transpose(const OperandType& m)
        : m_operand(m.copy())
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

      constexpr
      RangeShape getRangeShape() const
      {
        return m_operand->getRangeShape().transpose();
      }

      constexpr
      const OperandType& getOperand() const
      {
        assert(m_operand);
        return *m_operand;
      }

      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return getOperand().getValue(p).transpose();
      }

      template <class T>
      constexpr
      void getValue(T& out, const Geometry::Point& p) const
      {
        getOperand().getValue(out, p);
        out.transposeInPlace();
      }

      Transpose* copy() const noexcept override
      {
        return new Transpose(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
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
      using FESType = FES;

      using OperandType = ShapeFunctionBase<NestedDerived, FES, Space>;

      using Parent = ShapeFunctionBase<Transpose<OperandType>, FES, Space>;

      constexpr
      Transpose(const OperandType& op)
        : Parent(op.getFiniteElementSpace()),
          m_operand(op.copy())
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

      constexpr
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }

      constexpr
      RangeShape getRangeShape() const
      {
        return getOperand().getRangeShape().transpose();
      }

      constexpr
      size_t getDOFs(const Geometry::Polytope& simplex) const
      {
        return getOperand().getDOFs(simplex);
      }

      constexpr
      const OperandType& getOperand() const
      {
        assert(m_operand);
        return *m_operand;
      }

      const Geometry::Point& getPoint() const
      {
        return m_operand->getPoint();
      }

      Transpose& setPoint(const Geometry::Point& p)
      {
        m_operand->setPoint(p);
        return *this;
      }

      constexpr
      auto getBasis(size_t local) const
      {
        return getOperand().getBasis(local).transpose();
      }

      constexpr
      const FES& getFiniteElementSpace() const
      {
        return m_operand.getFiniteElementSpace();
      }

      Transpose* copy() const noexcept override
      {
        return new Transpose(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
  };

  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  Transpose(const ShapeFunctionBase<NestedDerived, FES, Space>&)
    -> Transpose<ShapeFunctionBase<NestedDerived, FES, Space>>;
}

#endif
