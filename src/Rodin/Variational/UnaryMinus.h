/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_UNARYMINUS_H
#define RODIN_VARIATIONAL_UNARYMINUS_H

#include "Rodin/FormLanguage/Base.h"
#include "Rodin/FormLanguage/List.h"

#include "ForwardDecls.h"
#include "Function.h"
#include "ShapeFunction.h"
#include "LinearFormIntegrator.h"
#include "BilinearFormIntegrator.h"

namespace Rodin::FormLanguage
{
  template <class NestedDerived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct Traits<Variational::UnaryMinus<Variational::ShapeFunctionBase<NestedDerived, FES, Space>>>
  {
    using FESType = FES;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup UnaryMinusSpecializations UnaryMinus Template Specializations
   * @brief Template specializations of the UnaryMinus class.
   * @see UnaryMinus
   */

  /**
   * @ingroup UnaryMinusSpecializations
   */
  template <class NestedDerived>
  class UnaryMinus<FunctionBase<NestedDerived>> final
    : public FunctionBase<UnaryMinus<FunctionBase<NestedDerived>>>
  {
    public:
      using OperandType = FunctionBase<NestedDerived>;

      using OperandRangeType = typename FormLanguage::Traits<OperandType>::RangeType;

      using Parent = FunctionBase<UnaryMinus<OperandType>>;

      constexpr
      UnaryMinus(const OperandType& op)
        : m_op(op.copy())
      {}

      constexpr
      UnaryMinus(const UnaryMinus& other)
        : Parent(other),
          m_op(other.m_op->copy())
      {}

      constexpr
      UnaryMinus(UnaryMinus&& other)
        : Parent(std::move(other)),
          m_op(std::move(other.m_op))
      {}

      constexpr
      RangeShape getRangeShape() const
      {
        return getOperand().getRangeShape();
      }

      constexpr
      const OperandType& getOperand() const
      {
        assert(m_op);
        return *m_op;
      }

      auto getValue(const Geometry::Point& p) const
      {
        return -this->object(getOperand().getValue(p));
      }

      template <class T>
      constexpr
      void getValue(T& res, const Geometry::Point& p) const
      {
        getOperand().getValue(res, p);
        res *= -1;
      }

      constexpr
      UnaryMinus& traceOf(Geometry::Attribute attr)
      {
        Parent::traceOf(attr);
        m_op->traceOf(attr);
        return *this;
      }

      constexpr
      UnaryMinus& traceOf(const FlatSet<Geometry::Attribute>& attrs)
      {
        Parent::traceOf(attrs);
        m_op->traceOf(attrs);
        return *this;
      }

      UnaryMinus* copy() const noexcept override
      {
        return new UnaryMinus(*this);
      }

    private:
      std::unique_ptr<OperandType> m_op;
  };

  template <class NestedDerived>
  UnaryMinus(const FunctionBase<NestedDerived>&) -> UnaryMinus<FunctionBase<NestedDerived>>;

  template <class NestedDerived>
  constexpr
  auto operator-(const FunctionBase<NestedDerived>& op)
  {
    return UnaryMinus(op);
  }

  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  class UnaryMinus<ShapeFunctionBase<NestedDerived, FES, Space>> final
    : public ShapeFunctionBase<UnaryMinus<ShapeFunctionBase<NestedDerived, FES, Space>>>
  {
    public:
      using FESType = FES;
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      using OperandType = ShapeFunctionBase<NestedDerived, FES, Space>;

      using Parent = ShapeFunctionBase<UnaryMinus<ShapeFunctionBase<NestedDerived, FES, Space>>, FES, Space>;

      constexpr
      UnaryMinus(const OperandType& op)
        : Parent(op.getFiniteElementSpace()),
          m_operand(op.copy())
      {}

      constexpr
      UnaryMinus(const UnaryMinus& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      constexpr
      UnaryMinus(UnaryMinus&& other)
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
      UnaryMinus& setPoint(const Geometry::Point& p)
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
        return -this->object(getOperand().getBasis(local));
      }

      const FES& getFiniteElementSpace() const
      {
        return getOperand().getFiniteElementSpace();
      }

      UnaryMinus* copy() const noexcept override
      {
        return new UnaryMinus(*this);
      }
    private:
      std::unique_ptr<OperandType> m_operand;
  };

  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  UnaryMinus(const ShapeFunctionBase<NestedDerived, FES, Space>&)
    -> UnaryMinus<ShapeFunctionBase<NestedDerived, FES, Space>>;


  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  constexpr
  auto operator-(const ShapeFunctionBase<NestedDerived, FES, Space>& op)
  {
    return UnaryMinus(op);
  }

  template <class Number>
  class UnaryMinus<LinearFormIntegratorBase<Number>> final
    : public LinearFormIntegratorBase<Number>
  {
    public:
      using ScalarType = Number;

      using OperandType = LinearFormIntegratorBase<ScalarType>;

      using Parent = LinearFormIntegratorBase<ScalarType>;

      UnaryMinus(const OperandType& op)
        : Parent(op),
          m_op(op.copy())
      {}

      UnaryMinus(const UnaryMinus& other)
        : Parent(other),
          m_op(other.m_op->copy())
      {}

      UnaryMinus(UnaryMinus&& other)
        : Parent(std::move(other)),
          m_op(std::move(other.m_op))
      {}

      const OperandType& getOperand() const
      {
        assert(m_op);
        return *m_op;
      }

      Integrator::Region getRegion() const override
      {
        return getOperand().getRegion();
      }

      const Geometry::Polytope& getPolytope() const override
      {
        return getOperand().getPolytope();
      }

      UnaryMinus& setPolytope(const Geometry::Polytope& polytope) override
      {
        m_op->setPolytope(polytope);
        return *this;
      }

      ScalarType integrate(size_t local) override
      {
        return -m_op->integrate(local);
      }

      UnaryMinus* copy() const noexcept override
      {
        return new UnaryMinus(*this);
      }

    private:
      std::unique_ptr<OperandType> m_op;
  };

  template <class Number>
  UnaryMinus(const LinearFormIntegratorBase<Number>&)
    -> UnaryMinus<LinearFormIntegratorBase<Number>>;

  template <class Number>
  UnaryMinus<LinearFormIntegratorBase<Number>> operator-(const LinearFormIntegratorBase<Number>& lfi)
  {
    return UnaryMinus(lfi);
  }

  template <class Number>
  class UnaryMinus<FormLanguage::List<LinearFormIntegratorBase<Number>>>
    : public FormLanguage::List<LinearFormIntegratorBase<Number>>
  {
    public:
      using ScalarType = Number;

      using LinearFormIntegratorBaseType = LinearFormIntegratorBase<Number>;

      using OperandType = FormLanguage::List<LinearFormIntegratorBaseType>;

      using Parent = FormLanguage::List<LinearFormIntegratorBaseType>;

      UnaryMinus(const OperandType& op)
      {
        for (const auto& p : op)
          add(UnaryMinus<LinearFormIntegratorBaseType>(p));
      }

      UnaryMinus(const UnaryMinus& other)
        : Parent(other)
      {}

      UnaryMinus(UnaryMinus&& other)
        : Parent(std::move(other))
      {}

      UnaryMinus* copy() const noexcept override
      {
        return new UnaryMinus(*this);
      }
  };

  template <class Number>
  UnaryMinus(const FormLanguage::List<LinearFormIntegratorBase<Number>>&)
    -> UnaryMinus<FormLanguage::List<LinearFormIntegratorBase<Number>>>;

  template <class Number>
  UnaryMinus<FormLanguage::List<LinearFormIntegratorBase<Number>>>
  operator-(const FormLanguage::List<LinearFormIntegratorBase<Number>>& op)
  {
    return UnaryMinus(op);
  }

  template <class Number>
  class UnaryMinus<LocalBilinearFormIntegratorBase<Number>>
    : public LocalBilinearFormIntegratorBase<Number>
  {
    public:
      using ScalarType = Number;

      using OperandType = LocalBilinearFormIntegratorBase<ScalarType>;

      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      UnaryMinus(const OperandType& op)
        : Parent(op),
          m_op(op.copy())
      {}

      UnaryMinus(const UnaryMinus& other)
        : Parent(other),
          m_op(other.m_op->copy())
      {}

      UnaryMinus(UnaryMinus&& other)
        : Parent(std::move(other)),
          m_op(std::move(other.m_op))
      {}

      const OperandType& getOperand() const
      {
        assert(m_op);
        return *m_op;
      }

      Integrator::Region getRegion() const override
      {
        return getOperand().getRegion();
      }

      const Geometry::Polytope& getPolytope() const override
      {
        return getOperand().getPolytope();
      }

      UnaryMinus& setPolytope(const Geometry::Polytope& polytope) override
      {
        m_op->setPolytope(polytope);
        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) override
      {
        return -m_op->integrate(tr, te);
      }

      UnaryMinus* copy() const noexcept override
      {
        return new UnaryMinus(*this);
      }

    private:
      std::unique_ptr<OperandType> m_op;
  };

  template <class Number>
  UnaryMinus(const LocalBilinearFormIntegratorBase<Number>&)
    -> UnaryMinus<LocalBilinearFormIntegratorBase<Number>>;

  template <class Number>
  constexpr
  auto
  operator-(const LocalBilinearFormIntegratorBase<Number>& op)
  {
    return UnaryMinus(op);
  }

  template <class Number>
  class UnaryMinus<FormLanguage::List<LocalBilinearFormIntegratorBase<Number>>>
    : public FormLanguage::List<LocalBilinearFormIntegratorBase<Number>>
  {
    public:
      using ScalarType = Number;

      using LocalBilinearFormIntegratorBaseType =
        LocalBilinearFormIntegratorBase<ScalarType>;

      using OperandType =
        FormLanguage::List<LocalBilinearFormIntegratorBaseType>;

      using Parent =
        FormLanguage::List<LocalBilinearFormIntegratorBaseType>;

      UnaryMinus(const OperandType& op)
      {
        for (const auto& p : op)
          add(UnaryMinus<LocalBilinearFormIntegratorBaseType>(p));
      }

      UnaryMinus(const UnaryMinus& other)
        : Parent(other)
      {}

      UnaryMinus(UnaryMinus&& other)
        : Parent(std::move(other))
      {}

      UnaryMinus* copy() const noexcept override
      {
        return new UnaryMinus(*this);
      }
  };

  template <class Number>
  UnaryMinus(const FormLanguage::List<LocalBilinearFormIntegratorBase<Number>>&)
    -> UnaryMinus<FormLanguage::List<LocalBilinearFormIntegratorBase<Number>>>;

  template <class Number>
  constexpr
  auto
  operator-(const FormLanguage::List<LocalBilinearFormIntegratorBase<Number>>& op)
  {
    return UnaryMinus(op);
  }
}

#endif
