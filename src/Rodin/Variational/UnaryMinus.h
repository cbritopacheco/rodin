/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_UNARYMINUS_H
#define RODIN_VARIATIONAL_UNARYMINUS_H

#include "Rodin/Variational/LinearFormIntegrator.h"

#include "Rodin/FormLanguage/Base.h"

#include "ForwardDecls.h"
#include "Function.h"
#include "ShapeFunction.h"
#include "ScalarFunction.h"
#include "LinearFormIntegrator.h"

#include "LinearFormIntegratorSum.h"
#include "BilinearFormIntegratorSum.h"

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
         UnaryMinus(const FunctionBase& op);

         UnaryMinus(const UnaryMinus& other);

         UnaryMinus(UnaryMinus&& other);

         RangeShape getRangeShape() const override;

         void getValue(
            mfem::DenseMatrix& value,
            mfem::ElementTransformation& trans,
            const mfem::IntegrationPoint& ip) const override;

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
         UnaryMinus(const LinearFormIntegratorBase& op);

         UnaryMinus(const UnaryMinus& other);

         UnaryMinus(UnaryMinus&& other);

         IntegratorRegion getIntegratorRegion() const override;

         bool isSupported(Linear::Assembly::Type t) const override;

         void getElementVector(const Linear::Assembly::Device& as) const override;

         void getElementVector(const Linear::Assembly::Common& as) const override;

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
   class UnaryMinus<LinearFormIntegratorSum> : public LinearFormIntegratorSum
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
   UnaryMinus<LinearFormIntegratorSum>
   operator-(const LinearFormIntegratorSum& lfi);

   template <>
   class UnaryMinus<BilinearFormIntegratorBase> : public BilinearFormIntegratorBase
   {
      public:
         UnaryMinus(const BilinearFormIntegratorBase& op);

         UnaryMinus(const UnaryMinus& other);

         UnaryMinus(UnaryMinus&& other);

         IntegratorRegion getIntegratorRegion() const override;

         void getElementMatrix(const Bilinear::Assembly::Common& as) const override;

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
   class UnaryMinus<BilinearFormIntegratorSum>
      : public BilinearFormIntegratorSum
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
   UnaryMinus<BilinearFormIntegratorSum> operator-(const BilinearFormIntegratorSum& op);

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
}

#endif
