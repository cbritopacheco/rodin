/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_CONJUGATE_H
#define RODIN_VARIATIONAL_CONJUGATE_H

#include "Rodin/Math/Traits.h"

#include "ForwardDecls.h"

#include "Function.h"
#include "RealFunction.h"
#include "ShapeFunction.h"

namespace Rodin::FormLanguage
{
  template <class NestedDerived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Conjugate<Variational::ShapeFunctionBase<NestedDerived, FES, Space>>>
  {
    using FESType = FES;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup ConjugateSpecializations Conjugate Template Specializations
   * @brief Template specializations of the Conjugate class.
   * @see Conjugate
   */

  /**
   * @ingroup ConjugateSpecializations
   */
  template <class NestedDerived>
  class Conjugate<FunctionBase<NestedDerived>>
    : public FunctionBase<Conjugate<FunctionBase<NestedDerived>>>
  {
    public:
      using OperandType = FunctionBase<NestedDerived>;

      using Parent = FunctionBase<Conjugate<OperandType>>;

      Conjugate(const OperandType& v)
        : m_v(v.copy())
      {}

      Conjugate(const Conjugate& other)
        : Parent(other),
          m_v(other.m_v->copy())
      {}

      Conjugate(Conjugate&& other)
        : Parent(std::move(other)),
          m_v(std::move(other.m_v))
      {}

      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return Math::conj(this->object(getOperand().getValue(p)));
      }

      const OperandType& getOperand() const
      {
        assert(m_v);
        return *m_v;
      }

      Conjugate* copy() const noexcept override
      {
        return new Conjugate(*this);
      }

    private:
      std::unique_ptr<OperandType> m_v;
  };

  template <class NestedDerived>
  Conjugate(const FunctionBase<NestedDerived>&) -> Conjugate<FunctionBase<NestedDerived>>;

  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  class Conjugate<ShapeFunctionBase<NestedDerived, FES, Space>> final
    : public ShapeFunctionBase<Conjugate<ShapeFunctionBase<NestedDerived, FES, Space>>>
  {
    public:
      using FESType = FES;
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      using OperandType =
        ShapeFunctionBase<NestedDerived, FES, Space>;

      using RangeType =
        typename FormLanguage::Traits<OperandType>::RangeType;

      using Parent =
        ShapeFunctionBase<Conjugate<ShapeFunctionBase<NestedDerived, FES, Space>>>;

      constexpr
      Conjugate(const OperandType& op)
        : Parent(op.getFiniteElementSpace()),
          m_operand(op.copy())
      {}

      constexpr
      Conjugate(const Conjugate& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      constexpr
      Conjugate(Conjugate&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}

      constexpr
      const OperandType& getOperand() const
      {
        return *m_operand;
      }

      constexpr
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }

      constexpr
      RangeShape getRangeShape() const
      {
        return getOperand().getRangeShape();
      }

      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getOperand().getDOFs(element);
      }

      constexpr
      Conjugate& setPoint(const Geometry::Point& p)
      {
        m_operand->setPoint(p);
        return *this;
      }

      const Geometry::Point& getPoint() const
      {
        return m_operand->getPoint();
      }

      constexpr
      auto getBasis(size_t local) const
      {
        return Math::conj(this->object(getOperand().getBasis(local)));
      }

      const FES& getFiniteElementSpace() const
      {
        return getOperand().getFiniteElementSpace();
      }

      Conjugate* copy() const noexcept override
      {
        return new Conjugate(*this);
      }
    private:
      std::unique_ptr<OperandType> m_operand;
  };

  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  Conjugate(const ShapeFunctionBase<NestedDerived, FES, Space>&)
    -> Conjugate<ShapeFunctionBase<NestedDerived, FES, Space>>;
}

#endif



