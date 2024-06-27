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
#include "ScalarFunction.h"
#include "LinearFormIntegrator.h"
#include "BilinearFormIntegrator.h"

namespace Rodin::FormLanguage
{
  template <class NestedDerived, class FESType, Variational::ShapeFunctionSpaceType SpaceType>
  struct Traits<Variational::UnaryMinus<Variational::ShapeFunctionBase<NestedDerived, FESType, SpaceType>>>
  {
    using FES = FESType;
    static constexpr Variational::ShapeFunctionSpaceType Space = SpaceType;
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
      using Operand = FunctionBase<NestedDerived>;
      using Parent = FunctionBase<UnaryMinus<Operand>>;
      using OperandRange = typename FormLanguage::Traits<Operand>::RangeType;

      constexpr
      UnaryMinus(const Operand& op)
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

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        return getOperand().getRangeShape();
      }

      inline
      constexpr
      const Operand& getOperand() const
      {
        assert(m_op);
        return *m_op;
      }

      inline
      auto getValue(const Geometry::Point& p) const
      {
        return -1 * this->object(getOperand().getValue(p));
      }

      inline
      constexpr
      void getValue(Math::Vector<Scalar>& res, const Geometry::Point& p) const
      {
        static_assert(FormLanguage::IsVectorRange<OperandRange>::Value);
        getOperand().getValue(res, p);
        res *= -1;
      }

      inline
      constexpr
      void getValue(Math::Matrix<Scalar>& res, const Geometry::Point& p) const
      {
        static_assert(FormLanguage::IsMatrixRange<OperandRange>::Value);
        getOperand().getValue(res, p);
        res *= -1;
      }

      inline
      constexpr
      UnaryMinus& traceOf(Geometry::Attribute attr)
      {
        Parent::traceOf(attr);
        m_op->traceOf(attr);
        return *this;
      }

      inline
      constexpr
      UnaryMinus& traceOf(const FlatSet<Geometry::Attribute>& attrs)
      {
        Parent::traceOf(attrs);
        m_op->traceOf(attrs);
        return *this;
      }

      inline UnaryMinus* copy() const noexcept override
      {
        return new UnaryMinus(*this);
      }

    private:
      std::unique_ptr<Operand> m_op;
  };

  template <class NestedDerived>
  UnaryMinus(const FunctionBase<NestedDerived>&) -> UnaryMinus<FunctionBase<NestedDerived>>;

  template <class NestedDerived>
  inline
  constexpr
  auto operator-(const FunctionBase<NestedDerived>& op)
  {
    return UnaryMinus(op);
  }

  template <class NestedDerived, class FESType, ShapeFunctionSpaceType SpaceType>
  class UnaryMinus<ShapeFunctionBase<NestedDerived, FESType, SpaceType>> final
    : public ShapeFunctionBase<UnaryMinus<ShapeFunctionBase<NestedDerived, FESType, SpaceType>>>
  {
    public:
      using FES = FESType;
      static constexpr ShapeFunctionSpaceType Space = SpaceType;

      using Operand = ShapeFunctionBase<NestedDerived, FES, Space>;
      using Parent = ShapeFunctionBase<UnaryMinus<ShapeFunctionBase<NestedDerived, FES, Space>>, FES, Space>;

      constexpr
      UnaryMinus(const Operand& op)
        : Parent(op.getFiniteElementSpace()),
          m_op(op.copy())
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
        return getOperand().getRangeShape();
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
        return -1 * this->object(getOperand().getTensorBasis(p));
      }

      const FES& getFiniteElementSpace() const
      {
        return getOperand().getFiniteElementSpace();
      }

      inline UnaryMinus* copy() const noexcept override
      {
        return new UnaryMinus(*this);
      }
    private:
      std::unique_ptr<Operand> m_op;
  };

  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  UnaryMinus(const ShapeFunctionBase<NestedDerived, FES, Space>&)
    -> UnaryMinus<ShapeFunctionBase<NestedDerived, FES, Space>>;


  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  inline
  constexpr
  auto operator-(const ShapeFunctionBase<NestedDerived, FES, Space>& op)
  {
    return UnaryMinus(op);
  }

  template <>
  class UnaryMinus<LinearFormIntegratorBase> final : public LinearFormIntegratorBase
  {
    public:
      using Parent = LinearFormIntegratorBase;

      UnaryMinus(const LinearFormIntegratorBase& op);

      UnaryMinus(const UnaryMinus& other);

      UnaryMinus(UnaryMinus&& other);

      Region getRegion() const override;

      void assemble(const Geometry::Polytope& element) override;

      UnaryMinus* copy() const noexcept override
      {
        return new UnaryMinus(*this);
      }

    private:
      std::unique_ptr<LinearFormIntegratorBase> m_op;
  };
  UnaryMinus(const LinearFormIntegratorBase&) -> UnaryMinus<LinearFormIntegratorBase>;

  UnaryMinus<LinearFormIntegratorBase> operator-(const LinearFormIntegratorBase& lfi);

  template <>
  class UnaryMinus<LocalBilinearFormIntegratorBase> : public LocalBilinearFormIntegratorBase
  {
    public:
      UnaryMinus(const LocalBilinearFormIntegratorBase& op);

      UnaryMinus(const UnaryMinus& other);

      UnaryMinus(UnaryMinus&& other);

      Region getRegion() const override;

      void assemble(const Geometry::Polytope& element) override;

      UnaryMinus* copy() const noexcept override
      {
        return new UnaryMinus(*this);
      }

    private:
      std::unique_ptr<LocalBilinearFormIntegratorBase> m_op;
  };
  UnaryMinus(const LocalBilinearFormIntegratorBase&) -> UnaryMinus<LocalBilinearFormIntegratorBase>;

  UnaryMinus<LocalBilinearFormIntegratorBase> operator-(const LocalBilinearFormIntegratorBase& op);

  template <>
  class UnaryMinus<FormLanguage::List<LinearFormIntegratorBase>>
    : public FormLanguage::List<LinearFormIntegratorBase>
  {
    public:
      UnaryMinus(const FormLanguage::List<LinearFormIntegratorBase>& op)
      {
        for (const auto& p : op)
          add(UnaryMinus<LinearFormIntegratorBase>(p));
      }

      UnaryMinus(const UnaryMinus& other)
        : FormLanguage::List<LinearFormIntegratorBase>(other)
      {}

      UnaryMinus(UnaryMinus&& other)
        : FormLanguage::List<LinearFormIntegratorBase>(std::move(other))
      {}

      UnaryMinus* copy() const noexcept override
      {
        return new UnaryMinus(*this);
      }
  };
  UnaryMinus(const FormLanguage::List<LinearFormIntegratorBase>&)
    -> UnaryMinus<FormLanguage::List<LinearFormIntegratorBase>>;

  UnaryMinus<FormLanguage::List<LinearFormIntegratorBase>>
  operator-(const FormLanguage::List<LinearFormIntegratorBase>& op);

  template <>
  class UnaryMinus<FormLanguage::List<LocalBilinearFormIntegratorBase>>
    : public FormLanguage::List<LocalBilinearFormIntegratorBase>
  {
    public:
      UnaryMinus(const FormLanguage::List<LocalBilinearFormIntegratorBase>& op)
      {
        for (const auto& p : op)
          add(UnaryMinus<LocalBilinearFormIntegratorBase>(p));
      }

      UnaryMinus(const UnaryMinus& other)
        : FormLanguage::List<LocalBilinearFormIntegratorBase>(other)
      {}

      UnaryMinus(UnaryMinus&& other)
        : FormLanguage::List<LocalBilinearFormIntegratorBase>(std::move(other))
      {}

      UnaryMinus* copy() const noexcept override
      {
        return new UnaryMinus(*this);
      }
  };

  UnaryMinus(const FormLanguage::List<LocalBilinearFormIntegratorBase>&)
    -> UnaryMinus<FormLanguage::List<LocalBilinearFormIntegratorBase>>;

  UnaryMinus<FormLanguage::List<LocalBilinearFormIntegratorBase>>
  operator-(const FormLanguage::List<LocalBilinearFormIntegratorBase>& op);
}

#endif
