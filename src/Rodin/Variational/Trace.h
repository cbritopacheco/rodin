/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_TRACE_H
#define RODIN_VARIATIONAL_TRACE_H

#include "ForwardDecls.h"

#include "ScalarFunction.h"
#include "ShapeFunction.h"

namespace Rodin::FormLanguage
{
  template <class NestedDerived, class FESType, Variational::ShapeFunctionSpaceType SpaceType>
  struct Traits<Variational::Trace<Variational::ShapeFunctionBase<NestedDerived, FESType, SpaceType>>>
  {
    using FES = FESType;
    static constexpr Variational::ShapeFunctionSpaceType Space = SpaceType;
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
    : public ScalarFunctionBase<Trace<FunctionBase<NestedDerived>>>
  {
    public:
      using Operand = FunctionBase<NestedDerived>;
      using Parent = ScalarFunctionBase<Trace<Operand>>;

      /**
       * @brief Constructs the Trace of the given matrix
       * @param[in] m Square matrix
       */
      constexpr
      Trace(const Operand& m)
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
        using OperandRange = typename FormLanguage::Traits<Operand>::RangeType;
        static_assert(std::is_same_v<OperandRange, Math::Matrix<Scalar>>);
        return getOperand().getValue(p).trace();
      }

      inline
      constexpr
      const Operand& getOperand() const
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
      std::unique_ptr<Operand> m_operand;
  };

  template <class NestedDerived>
  Trace(const FunctionBase<NestedDerived>&) -> Trace<FunctionBase<NestedDerived>>;

  template <class NestedDerived, class FESType, ShapeFunctionSpaceType SpaceType>
  class Trace<ShapeFunctionBase<NestedDerived, FESType, SpaceType>> final
    : public ShapeFunctionBase<Trace<ShapeFunctionBase<NestedDerived, FESType, SpaceType>>>
  {
    public:
      using FES = FESType;
      static constexpr ShapeFunctionSpaceType Space = SpaceType;

      using Operand = ShapeFunctionBase<NestedDerived, FES, Space>;
      using Parent = ShapeFunctionBase<Trace<Operand>>;

      constexpr
      Trace(const Operand& op)
        : Parent(op.getFiniteElementSpace()),
          m_op(op.copy())
      {}

      constexpr
      Trace(const Trace& other)
        : Parent(other),
          m_op(other.m_op->copy())
      {}

      constexpr
      Trace(Trace&& other)
        : Parent(std::move(other)),
          m_op(std::move(other.m_op))
      {}

      inline
      constexpr
      const Operand& getOperand() const
      {
        return *m_op;
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
      constexpr
      auto getTensorBasis(const Geometry::Point& p) const
      {
        const size_t d = p.getPolytope().getDimension();
        const Index i = p.getPolytope().getIndex();
        const auto& fe = this->getFiniteElementSpace().getFiniteElement(d, i);
        const auto& tb = m_op->getTensorBasis(p);
        return TensorBasis(fe.getCount(),
            [&](size_t local) { return this->object(tb(local)).transpose(); } );
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
      std::unique_ptr<Operand> m_op;
  };

  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  Trace(const ShapeFunctionBase<NestedDerived, FES, Space>&)
    -> Trace<ShapeFunctionBase<NestedDerived, FES, Space>>;
}

#endif
