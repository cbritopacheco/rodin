/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_TRACE_H
#define RODIN_VARIATIONAL_TRACE_H

#include "ForwardDecls.h"

#include "RealFunction.h"
#include "ShapeFunction.h"

namespace Rodin::FormLanguage
{
  template <class NestedDerived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct Traits<Variational::Trace<Variational::ShapeFunctionBase<NestedDerived, FES, Space>>>
  {
    using FESType = FES;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup TraceSpecializations
   * @brief Template specializations of the Trace class.
   * @see Trace
   */

  /**
   * @ingroup TraceSpecializations
   * @brief Trace of a FunctionBase instance.
   */
  template <class NestedDerived>
  class Trace<FunctionBase<NestedDerived>> final
    : public RealFunctionBase<Trace<FunctionBase<NestedDerived>>>
  {
    public:
      using OperandType = FunctionBase<NestedDerived>;

      using Parent = RealFunctionBase<Trace<OperandType>>;

      /**
       * @brief Constructs the Trace of the given matrix
       * @param[in] m Square matrix
       */
      constexpr
      Trace(const OperandType& m)
        : m_operand(m.copy())
      {}

      constexpr
      Trace(const Trace& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      constexpr
      Trace(Trace&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}

      inline
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        using OperandRange = typename FormLanguage::Traits<OperandType>::RangeType;
        static_assert(std::is_same_v<OperandRange, Math::Matrix<Real>>);
        return getOperand().getValue(p).trace();
      }

      inline
      constexpr
      const OperandType& getOperand() const
      {
        assert(m_operand);
        return *m_operand;
      }

      inline
      Trace& traceOf(Geometry::Attribute attrs)
      {
        m_operand.traceOf(attrs);
        return *this;
      }

      inline Trace* copy() const noexcept override
      {
        return new Trace(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
  };

  template <class NestedDerived>
  Trace(const FunctionBase<NestedDerived>&) -> Trace<FunctionBase<NestedDerived>>;

  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  class Trace<ShapeFunctionBase<NestedDerived, FES, Space>> final
    : public ShapeFunctionBase<Trace<ShapeFunctionBase<NestedDerived, FES, Space>>>
  {
    public:
      using FESType = FES;
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      using OperandType = ShapeFunctionBase<NestedDerived, FES, Space>;

      using Parent = ShapeFunctionBase<Trace<OperandType>>;

      constexpr
      Trace(const OperandType& operand)
        : Parent(operand.getFiniteElementSpace()),
          m_operand(operand.copy())
      {}

      constexpr
      Trace(const Trace& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      constexpr
      Trace(Trace&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}

      inline
      constexpr
      const OperandType& getOperand() const
      {
        return *m_operand;
      }

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
        return { 1, 1 };
      }

      inline
      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getOperand().getDOFs(element);
      }

      inline
      const Geometry::Point& getPoint() const
      {
        return m_operand->getPoint();
      }

      inline
      Trace& setPoint(const Geometry::Point& p)
      {
        m_operand->setPoint(p);
        return *this;
      }

      inline
      constexpr
      auto getBasis(size_t local) const
      {
        return this->object(getOperand().getBasis(local)).transpose();
      }

      const FES& getFiniteElementSpace() const
      {
        return getOperand().getFiniteElementSpace();
      }

      inline Trace* copy() const noexcept override
      {
        return new Trace(*this);
      }
    private:
      std::unique_ptr<OperandType> m_operand;
  };

  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  Trace(const ShapeFunctionBase<NestedDerived, FES, Space>&)
    -> Trace<ShapeFunctionBase<NestedDerived, FES, Space>>;
}

#endif
