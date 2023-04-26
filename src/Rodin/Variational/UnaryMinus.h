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

  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  class UnaryMinus<ShapeFunctionBase<NestedDerived, FES, Space>> final
    : public ShapeFunctionBase<UnaryMinus<ShapeFunctionBase<NestedDerived, FES, Space>>, FES, Space>
  {
    public:
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
  class UnaryMinus<LinearFormIntegratorBase> : public LinearFormIntegratorBase
  {
    public:
      using Parent = LinearFormIntegratorBase;
      using Parent::Parent;

      UnaryMinus(const LinearFormIntegratorBase& op);

      UnaryMinus(const UnaryMinus& other);

      UnaryMinus(UnaryMinus&& other);

      Region getRegion() const override;

      Math::Vector getVector(const Geometry::Polytope& element) const override;

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

      Math::Matrix getMatrix(const Geometry::Polytope& element) const override;

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
