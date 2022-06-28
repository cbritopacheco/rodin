/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_UNARYMINUS_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_UNARYMINUS_H

#include <memory>
#include <type_traits>

#include "Rodin/Variational/ScalarFunction.h"
#include "Rodin/Variational/VectorFunction.h"
#include "Rodin/Variational/LinearFormIntegrator.h"

#include "FormLanguage/LinearFormIntegratorSum.h"
#include "FormLanguage/BilinearFormIntegratorSum.h"
#include "FormLanguage/Base.h"

#include "LinearFormIntegrator.h"

#include "ForwardDecls.h"

#include "Sum.h"

namespace Rodin::Variational
{
   template <class Operand>
   class UnaryMinus : public FormLanguage::Base
   {
      static_assert(std::is_base_of_v<Base, Operand>,
            "Operand must be derived from FormLanguage::Base");

      public:
         UnaryMinus(const Operand& op)
            : m_op(op.copy())
         {}

         UnaryMinus(const UnaryMinus& other)
            :  Base(other),
               m_op(other.m_op->copy())
         {}

         UnaryMinus(UnaryMinus&& other)
            : Base(std::move(other)),
              m_op(std::move(other.m_op))
         {}

         Operand& getOperand()
         {
            return *m_op;
         }

         const Operand& getOperand() const
         {
            return *m_op;
         }

         UnaryMinus* copy() const noexcept override
         {
            return new UnaryMinus(*this);
         }

      private:
         std::unique_ptr<Operand> m_op;
   };

   template <>
   class UnaryMinus<ScalarFunctionBase>
      : public ScalarFunctionBase
   {
      public:
         UnaryMinus(const ScalarFunctionBase& op)
            : m_op(op.copy())
         {}

         UnaryMinus(const UnaryMinus& other)
            :  ScalarFunctionBase(other),
               m_op(other.m_op->copy())
         {}

         UnaryMinus(UnaryMinus&& other)
            : ScalarFunctionBase(std::move(other)),
              m_op(std::move(other.m_op))
         {}

         ScalarFunctionBase& getOperand()
         {
            return *m_op;
         }

         const ScalarFunctionBase& getOperand() const
         {
            return *m_op;
         }

         double getValue(
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const override;

         UnaryMinus* copy() const noexcept override
         {
            return new UnaryMinus(*this);
         }

      private:
         std::unique_ptr<ScalarFunctionBase> m_op;
   };
   UnaryMinus<ScalarFunctionBase> operator-(const ScalarFunctionBase& op);

   template <>
   class UnaryMinus<VectorFunctionBase>
      : public VectorFunctionBase
   {
      public:
         UnaryMinus(const VectorFunctionBase& op)
            : m_op(op.copy())
         {}

         UnaryMinus(const UnaryMinus& other)
            :  VectorFunctionBase(other),
               m_op(other.m_op->copy())
         {}

         UnaryMinus(UnaryMinus&& other)
            : VectorFunctionBase(std::move(other)),
              m_op(std::move(other.m_op))
         {}

         VectorFunctionBase& getOperand()
         {
            return *m_op;
         }

         const VectorFunctionBase& getOperand() const
         {
            return *m_op;
         }

         int getDimension() const override
         {
            return getOperand().getDimension();
         }

         void getValue(
               mfem::Vector& value,
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip) const override;

         UnaryMinus* copy() const noexcept override
         {
            return new UnaryMinus(*this);
         }

      private:
         std::unique_ptr<VectorFunctionBase> m_op;
   };
   UnaryMinus<VectorFunctionBase> operator-(const VectorFunctionBase& op);

   template <>
   class UnaryMinus<LinearFormIntegratorBase>
      : public LinearFormIntegratorBase
   {
      public:
         UnaryMinus(const LinearFormIntegratorBase& op)
            : m_op(op.copy())
         {}

         UnaryMinus(const UnaryMinus& other)
            :  LinearFormIntegratorBase(other),
               m_op(other.m_op->copy())
         {}

         UnaryMinus(UnaryMinus&& other)
            : LinearFormIntegratorBase(std::move(other)),
              m_op(std::move(other.m_op))
         {}

         LinearFormIntegratorBase& getOperand()
         {
            return *m_op;
         }

         const LinearFormIntegratorBase& getOperand() const
         {
            return *m_op;
         }

         const std::set<int>& getAttributes() const override
         {
            return getOperand().getAttributes();
         }

         IntegratorRegion getIntegratorRegion() const override
         {
            return getOperand().getIntegratorRegion();
         }

         const ShapeFunctionBase<Test>& getTestFunction() const override
         {
            return m_op->getTestFunction();
         }

         void getElementVector(const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans, mfem::Vector& vec) const override;

         UnaryMinus* copy() const noexcept override
         {
            return new UnaryMinus(*this);
         }

      private:
         std::unique_ptr<LinearFormIntegratorBase> m_op;
   };
   UnaryMinus<LinearFormIntegratorBase> operator-(const LinearFormIntegratorBase& lfi);

   template <>
   class UnaryMinus<FormLanguage::LinearFormIntegratorSum>
      : public FormLanguage::LinearFormIntegratorSum
   {
      public:
         UnaryMinus(const LinearFormIntegratorSum& op)
            : LinearFormIntegratorSum(op)
         {
            for (auto& p : getLinearFormDomainIntegratorList())
               p.reset(new UnaryMinus<LinearFormIntegratorBase>(*p));
            for (auto& p : getLinearFormBoundaryIntegratorList())
               p.reset(new UnaryMinus<LinearFormIntegratorBase>(*p));
         }

         UnaryMinus(const UnaryMinus& other)
            : LinearFormIntegratorSum(other)
         {}

         UnaryMinus(UnaryMinus&& other)
            : LinearFormIntegratorSum(std::move(other))
         {}

         UnaryMinus* copy() const noexcept override
         {
            return new UnaryMinus(*this);
         }
   };
   UnaryMinus<FormLanguage::LinearFormIntegratorSum>
   operator-(const FormLanguage::LinearFormIntegratorSum& lfi);

   template <>
   class UnaryMinus<BilinearFormIntegratorBase>
      : public BilinearFormIntegratorBase
   {
      public:
         UnaryMinus(const BilinearFormIntegratorBase& op)
            : m_op(op.copy())
         {}

         UnaryMinus(const UnaryMinus& other)
            :  BilinearFormIntegratorBase(other),
               m_op(other.m_op->copy())
         {}

         UnaryMinus(UnaryMinus&& other)
            : BilinearFormIntegratorBase(std::move(other)),
              m_op(std::move(other.m_op))
         {}

         BilinearFormIntegratorBase& getOperand()
         {
            return *m_op;
         }

         const BilinearFormIntegratorBase& getOperand() const
         {
            return *m_op;
         }

         const ShapeFunctionBase<Trial>& getTrialFunction() const override
         {
            return m_op->getTrialFunction();
         }

         const ShapeFunctionBase<Test>& getTestFunction() const override
         {
            return m_op->getTestFunction();
         }

         const std::set<int>& getAttributes() const override
         {
            return getOperand().getAttributes();
         }

         IntegratorRegion getIntegratorRegion() const override
         {
            return getOperand().getIntegratorRegion();
         }

         void getElementMatrix(const mfem::FiniteElement& trial, const mfem::FiniteElement& test,
               mfem::ElementTransformation& trans, mfem::DenseMatrix& vec) const override;

         UnaryMinus* copy() const noexcept override
         {
            return new UnaryMinus(*this);
         }

      private:
         std::unique_ptr<BilinearFormIntegratorBase> m_op;
   };
   UnaryMinus<BilinearFormIntegratorBase> operator-(const BilinearFormIntegratorBase& op);

   template <>
   class UnaryMinus<FormLanguage::BilinearFormIntegratorSum>
      : public FormLanguage::BilinearFormIntegratorSum
   {
      public:
         UnaryMinus(const BilinearFormIntegratorSum& op)
            : BilinearFormIntegratorSum(op)
         {
            for (auto& p : getBilinearFormDomainIntegratorList())
               p.reset(new UnaryMinus<BilinearFormIntegratorBase>(*p));
            for (auto& p : getBilinearFormBoundaryIntegratorList())
               p.reset(new UnaryMinus<BilinearFormIntegratorBase>(*p));
         }

         UnaryMinus(const UnaryMinus& other)
            : BilinearFormIntegratorSum(other)
         {}

         UnaryMinus(UnaryMinus&& other)
            : BilinearFormIntegratorSum(std::move(other))
         {}

         UnaryMinus* copy() const noexcept override
         {
            return new UnaryMinus(*this);
         }
   };
   UnaryMinus<FormLanguage::BilinearFormIntegratorSum> operator-(
         const FormLanguage::BilinearFormIntegratorSum& op);

   template <ShapeFunctionSpaceType Space>
   class UnaryMinus<ShapeFunctionBase<Space>> : public ShapeFunctionBase<Space>
   {
      public:
         UnaryMinus(const ShapeFunctionBase<Space>& rhs)
            : m_op(rhs.copy())
         {}

         UnaryMinus(const UnaryMinus& other)
            :  ShapeFunctionBase<Space>(other),
               m_op(other.m_op->copy())
         {}

         UnaryMinus(UnaryMinus&& other)
            :  ShapeFunctionBase<Space>(std::move(other)),
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

         std::unique_ptr<Rank3Operator> getOperator(
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
}

#endif
