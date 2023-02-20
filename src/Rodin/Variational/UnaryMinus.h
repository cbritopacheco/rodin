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

      constexpr
      UnaryMinus(const Operand& op)
        : m_op()
      {}

      constexpr
      UnaryMinus(const UnaryMinus& other)
        : Parent(other),
          m_op(other.m_op)
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
        return m_op.getRangeShape();
      }

      inline
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return -m_op.getValue(p);
      }

      inline
      UnaryMinus* copy() const noexcept
      override
      {
        return new UnaryMinus(*this);
      }

    private:
      Operand m_op;
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

  template <class NestedDerived, ShapeFunctionSpaceType Space>
  class UnaryMinus<ShapeFunctionBase<NestedDerived, Space>> final
    : public ShapeFunctionBase<NestedDerived, Space>
  {
    public:
      using Operand = ShapeFunctionBase<NestedDerived, Space>;
      using Parent = ShapeFunctionBase<NestedDerived, Space>;

      constexpr
      UnaryMinus(const Operand& op)
        : m_op(op)
      {}

      constexpr
      UnaryMinus(const UnaryMinus& other)
        : Parent(other),
          m_op(other.m_op)
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
        return m_op.getRangeShape();
      }

      inline
      constexpr
      size_t getDOFs(const Geometry::Simplex& element) const
      {
        return m_op.getDOFs(element);
      }

      inline
      constexpr
      auto getOperator(ShapeComputator& compute, const Geometry::Point& p) const
      {
        return m_op.getOperator(compute, p);
      }

      auto& getFiniteElementSpace()
      {
        return getOperand().getFiniteElementSpace();
      }

      const auto& getFiniteElementSpace() const
      {
        return getOperand().getFiniteElementSpace();
      }

      inline
      UnaryMinus* copy() const noexcept
      override
      {
        return new UnaryMinus(*this);
      }
    private:
      Operand m_op;
  };

  template <>
  class UnaryMinus<LinearFormIntegratorBase> : public LinearFormIntegratorBase
  {
    public:
      using Parent = LinearFormIntegratorBase;
      using Parent::Parent;

      UnaryMinus(const LinearFormIntegratorBase& op);

      UnaryMinus(const UnaryMinus& other);

      UnaryMinus(UnaryMinus&& other);

      Region getRegion() const override;

      Math::Vector getVector(const Geometry::Simplex& element) const override;

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
  class UnaryMinus<BilinearFormIntegratorBase> : public BilinearFormIntegratorBase
  {
    public:
      UnaryMinus(const BilinearFormIntegratorBase& op);

      UnaryMinus(const UnaryMinus& other);

      UnaryMinus(UnaryMinus&& other);

      Region getRegion() const override;

      Math::Matrix getMatrix(const Geometry::Simplex& element) const override;

      UnaryMinus* copy() const noexcept override
      {
        return new UnaryMinus(*this);
      }

    private:
      std::unique_ptr<BilinearFormIntegratorBase> m_op;
  };
  UnaryMinus(const BilinearFormIntegratorBase&) -> UnaryMinus<BilinearFormIntegratorBase>;

  UnaryMinus<BilinearFormIntegratorBase> operator-(const BilinearFormIntegratorBase& op);

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
  class UnaryMinus<FormLanguage::List<BilinearFormIntegratorBase>>
    : public FormLanguage::List<BilinearFormIntegratorBase>
  {
    public:
      UnaryMinus(const FormLanguage::List<BilinearFormIntegratorBase>& op)
      {
        for (const auto& p : op)
          add(UnaryMinus<BilinearFormIntegratorBase>(p));
      }

      UnaryMinus(const UnaryMinus& other)
        : FormLanguage::List<BilinearFormIntegratorBase>(other)
      {}

      UnaryMinus(UnaryMinus&& other)
        : FormLanguage::List<BilinearFormIntegratorBase>(std::move(other))
      {}

      UnaryMinus* copy() const noexcept override
      {
        return new UnaryMinus(*this);
      }
  };

  UnaryMinus(const FormLanguage::List<BilinearFormIntegratorBase>&)
    -> UnaryMinus<FormLanguage::List<BilinearFormIntegratorBase>>;

  UnaryMinus<FormLanguage::List<BilinearFormIntegratorBase>>
  operator-(const FormLanguage::List<BilinearFormIntegratorBase>& op);
}

#endif
