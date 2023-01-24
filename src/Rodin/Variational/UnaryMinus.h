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
  template <>
  class UnaryMinus<FunctionBase> : public FunctionBase
  {
    public:
      using Parent = FunctionBase;

      UnaryMinus(const FunctionBase& op);

      UnaryMinus(const UnaryMinus& other);

      UnaryMinus(UnaryMinus&& other);

      RangeShape getRangeShape() const override;

      FunctionValue getValue(const Geometry::Point& p) const override
      {
        return -m_op->getValue(p);
      }

      UnaryMinus* copy() const noexcept override
      {
        return new UnaryMinus(*this);
      }

    private:
      std::unique_ptr<FunctionBase> m_op;
  };
  UnaryMinus(const FunctionBase&) -> UnaryMinus<FunctionBase>;

  UnaryMinus<FunctionBase> operator-(const FunctionBase& op);

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

      mfem::Vector getVector(const Geometry::Simplex& element) const override;

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

      mfem::DenseMatrix getMatrix(const Geometry::Simplex& element) const override;

      UnaryMinus* copy() const noexcept override
      {
        return new UnaryMinus(*this);
      }

    private:
      std::unique_ptr<BilinearFormIntegratorBase> m_op;
  };
  UnaryMinus(const BilinearFormIntegratorBase&) -> UnaryMinus<BilinearFormIntegratorBase>;

  UnaryMinus<BilinearFormIntegratorBase> operator-(const BilinearFormIntegratorBase& op);

  template <ShapeFunctionSpaceType Space>
  class UnaryMinus<ShapeFunctionBase<Space>> : public ShapeFunctionBase<Space>
  {
    public:
      using Parent = ShapeFunctionBase<Space>;
      UnaryMinus(const ShapeFunctionBase<Space>& rhs)
        :  Parent(rhs),
          m_op(rhs.copy())
      {}

      UnaryMinus(const UnaryMinus& other)
        :  Parent(other),
          m_op(other.m_op->copy())
      {}

      UnaryMinus(UnaryMinus&& other)
        :  Parent(std::move(other)),
          m_op(std::move(other.m_op))
      {}

      ShapeFunctionBase<Space>& getOperand()
      {
        return *m_op;
      }

      const ShapeFunctionBase<Space>& getOperand() const
      {
        return *m_op;
      }

      ShapeFunctionBase<Space>& getLeaf() override
      {
        return getOperand().getLeaf();
      }

      const ShapeFunctionBase<Space>& getLeaf() const override
      {
        return getOperand().getLeaf();
      }

      int getRows(
          const mfem::FiniteElement& fe,
          const mfem::ElementTransformation& trans) const override
      {
        return getOperand().getRows(fe, trans);
      }

      int getColumns(
          const mfem::FiniteElement& fe,
          const mfem::ElementTransformation& trans) const override
      {
        return getOperand().getColumns(fe, trans);
      }

      int getDOFs(
          const mfem::FiniteElement& fe,
          const mfem::ElementTransformation& trans) const override
      {
        return getOperand().getDOFs(fe, trans);
      }

      std::unique_ptr<BasisOperator> getOperator(
          const mfem::FiniteElement& fe,
          mfem::ElementTransformation& trans) const override
      {
        auto result = getOperand().getOperator(fe, trans);
        (*result) *= -1.0;
        return result;
      }

      FiniteElementSpaceBase& getFiniteElementSpace() override
      {
        return getOperand().getFiniteElementSpace();
      }

      const FiniteElementSpaceBase& getFiniteElementSpace() const override
      {
        return getOperand().getFiniteElementSpace();
      }

      UnaryMinus* copy() const noexcept override
      {
        return new UnaryMinus(*this);
      }
    private:
      std::unique_ptr<ShapeFunctionBase<Space>> m_op;
  };

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
}

#endif
