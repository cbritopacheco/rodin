/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DOT_H
#define RODIN_VARIATIONAL_DOT_H

#include <memory>
#include <utility>
#include <optional>
#include <type_traits>

#include "Rodin/Alert.h"

#include "FormLanguage/Base.h"

#include "ForwardDecls.h"
#include "GridFunction.h"
#include "ScalarFunction.h"
#include "VectorFunction.h"
#include "MatrixFunction.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "Mult.h"

namespace Rodin::Variational
{
   /**
    * @brief Represents the dot product of two vectors
    *
    * For two given @f$ n @f$ dimensional vectors, their dot product @f$ u
    * \cdot v @f$ is defined as
    * @f[
    *    u \colon v = \sum_{i = 1} u_i v_i
    * @f]
    * where @f$ u_i @f$, @f$ v_i @f$ are the respective entries of each vector.
    */
   template <>
   class Dot<VectorFunctionBase, VectorFunctionBase> : public ScalarFunctionBase
   {
      public:
         /**
          * @brief Constructs the Dot product between two given matrices.
          * @param[in] a Derived instance of VectorFunctionBase
          * @param[in] b Derived instance of VectorFunctionBase
          */
         Dot(const VectorFunctionBase& a, const VectorFunctionBase& b)
            : m_a(a.copy()), m_b(b.copy())
         {}

         Dot(const Dot& other)
            :  ScalarFunctionBase(other),
               m_a(other.m_a->copy()), m_b(other.m_b->copy())
         {}

         Dot(Dot&& other)
            :  ScalarFunctionBase(std::move(other)),
               m_a(std::move(other.m_a)), m_b(std::move(other.m_b))
         {}

         double getValue(
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const override
         {
            mfem::Vector va, vb;
            m_a->getValue(va, trans, ip);
            m_b->getValue(vb, trans, ip);
            return va * vb;
         }

         Dot* copy() const noexcept override
         {
            return new Dot(*this);
         }
      private:
         std::unique_ptr<VectorFunctionBase> m_a, m_b;
   };
   Dot(const VectorFunctionBase&, const VectorFunctionBase&)
      -> Dot<VectorFunctionBase, VectorFunctionBase>;

   /**
    * @brief Represents the dot product between two matrices.
    *
    * For two given @f$ n \times m @f$ matrices @f$ A @f$ and @f$ B @f$, their
    * double-dot product @f$ A \colon B @f$ is defined as
    * @f[
    *    A \colon B = \sum_{i = 1}^n \sum_{j = 1}^m A_{ij} B_{ij} = \mathrm{tr}(B^T A)
    * @f]
    */
   template <>
   class Dot<MatrixFunctionBase, MatrixFunctionBase> : public ScalarFunctionBase
   {
      public:
         /**
          * @brief Constructs the Dot product between two given matrices.
          * @param[in] a Derived instance of MatrixFunctionBase
          * @param[in] b Derived instance of MatrixFunctionBase
          */
         Dot(const MatrixFunctionBase& a, const MatrixFunctionBase& b)
            : m_a(a.copy()), m_b(b.copy())
         {}

         Dot(const Dot& other)
            :  ScalarFunctionBase(other),
               m_a(other.m_a->copy()), m_b(other.m_b->copy())
         {}

         Dot(Dot&& other)
            :  ScalarFunctionBase(std::move(other)),
               m_a(std::move(other.m_a)), m_b(std::move(other.m_b))
         {}

         double getValue(
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const override
         {
            mfem::DenseMatrix ma, mb;
            m_a->getValue(ma, trans, ip);
            m_b->getValue(mb, trans, ip);
            return ma * mb;
         }

         Dot* copy() const noexcept override
         {
            return new Dot(*this);
         }
      private:
         std::unique_ptr<MatrixFunctionBase> m_a, m_b;
   };
   Dot(const MatrixFunctionBase&, const MatrixFunctionBase&)
      -> Dot<MatrixFunctionBase, MatrixFunctionBase>;

   /**
    * @brief Represents the dot product between insntances of
    * ScalarFunctionBase and ShapeFunctionBase.
    *
    * This corresponds to the usual multiplication.
    */
   template <ShapeFunctionSpaceType Space>
   class Dot<ScalarFunctionBase, ShapeFunctionBase<Space>>
      : public ShapeFunctionBase<Space>
   {
      public:
         Dot(const ShapeFunctionBase<Space>& lhs, const ScalarFunctionBase& rhs)
            : Dot(rhs, lhs)
         {}

         Dot(const ScalarFunctionBase& lhs, const ShapeFunctionBase<Space>& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         Dot(const Dot& other)
            :  ShapeFunctionBase<Space>(other),
               m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Dot(Dot&& other)
            :  ShapeFunctionBase<Space>(std::move(other)),
               m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
         {}

         ScalarFunctionBase& getLHS()
         {
            return *m_lhs;
         }

         ShapeFunctionBase<Space>& getRHS()
         {
            return *m_rhs;
         }

         const ScalarFunctionBase& getLHS() const
         {
            return *m_lhs;
         }

         const ShapeFunctionBase<Space>& getRHS() const
         {
            return *m_rhs;
         }

         ShapeFunctionBase<Space>& getLeaf() override
         {
            return getRHS().getLeaf();
         }

         const ShapeFunctionBase<Space>& getLeaf() const override
         {
            return getRHS().getLeaf();
         }

         int getRows(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const override
         {
            return getRHS().getRows(fe, trans);
         }

         int getColumns(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const override
         {
            return getRHS().getColumns(fe, trans);
         }

         int getDOFs(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const override
         {
            return getRHS().getDOFs(fe, trans);
         }

         std::unique_ptr<Rank3Operator> getOperator(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans) const override
         {
            assert(getRows(fe, trans) == 1);
            assert(getColumns(fe, trans) == 1);
            auto result = getRHS().getOperator(fe, trans);
            (*result) *= getLHS().getValue(trans, trans.GetIntPoint());
            return result;
         }

         FiniteElementSpaceBase& getFiniteElementSpace() override
         {
            return getRHS().getFiniteElementSpace();
         }

         const FiniteElementSpaceBase& getFiniteElementSpace() const override
         {
            return getRHS().getFiniteElementSpace();
         }

         Dot* copy() const noexcept override
         {
            return new Dot(*this);
         }
      private:
         std::unique_ptr<ScalarFunctionBase> m_lhs;
         std::unique_ptr<ShapeFunctionBase<Space>> m_rhs;
   };
   template <ShapeFunctionSpaceType Space>
   Dot(const ScalarFunctionBase&, const ShapeFunctionBase<Space>&)
      -> Dot<ScalarFunctionBase, ShapeFunctionBase<Space>>;
   template <ShapeFunctionSpaceType Space>
   Dot(const ShapeFunctionBase<Space>&, const ScalarFunctionBase&)
      -> Dot<ScalarFunctionBase, ShapeFunctionBase<Space>>;

   /**
    * @brief Dot product between VectorFunctionBase and ShapeFunctionBase.
    */
   template <ShapeFunctionSpaceType Space>
   class Dot<VectorFunctionBase, ShapeFunctionBase<Space>>
      : public ShapeFunctionBase<Space>
   {
      public:
         Dot(const ShapeFunctionBase<Space>& lhs, const VectorFunctionBase& rhs)
            : Dot(rhs, lhs)
         {}

         Dot(const VectorFunctionBase& lhs, const ShapeFunctionBase<Space>& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         Dot(const Dot& other)
            :  ShapeFunctionBase<Space>(other),
               m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Dot(Dot&& other)
            :  ShapeFunctionBase<Space>(std::move(other)),
               m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
         {}

         VectorFunctionBase& getLHS()
         {
            return *m_lhs;
         }

         ShapeFunctionBase<Space>& getRHS()
         {
            return *m_rhs;
         }

         const VectorFunctionBase& getLHS() const
         {
            return *m_lhs;
         }

         const ShapeFunctionBase<Space>& getRHS() const
         {
            return *m_rhs;
         }

         ShapeFunctionBase<Space>& getLeaf() override
         {
            return getRHS().getLeaf();
         }

         const ShapeFunctionBase<Space>& getLeaf() const override
         {
            return getRHS().getLeaf();
         }

         int getRows(
               const mfem::FiniteElement&,
               const mfem::ElementTransformation&) const override
         {
            return 1;
         }

         int getDOFs(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const override
         {
            return getRHS().getDOFs(fe, trans);
         }

         int getColumns(
               const mfem::FiniteElement&,
               const mfem::ElementTransformation&) const override
         {
            return 1;
         }

         std::unique_ptr<Rank3Operator> getOperator(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans) const override
         {
            assert(
               (getLHS().getDimension() == getRHS().getRows(fe, trans) && getRHS().getColumns(fe, trans) == 1)
                  || (getLHS().getDimension() == getRHS().getColumns(fe, trans) && getRHS().getRows(fe, trans) == 1));
            mfem::Vector v;
            getLHS().getValue(v, trans, trans.GetIntPoint());
            return getRHS().getOperator(fe, trans)->VectorDot(v);
         }

         FiniteElementSpaceBase& getFiniteElementSpace() override
         {
            return getRHS().getFiniteElementSpace();
         }

         const FiniteElementSpaceBase& getFiniteElementSpace() const override
         {
            return getRHS().getFiniteElementSpace();
         }

         Dot* copy() const noexcept override
         {
            return new Dot(*this);
         }
      private:
         std::unique_ptr<VectorFunctionBase> m_lhs;
         std::unique_ptr<ShapeFunctionBase<Space>> m_rhs;
   };
   template <ShapeFunctionSpaceType Space>
   Dot(const VectorFunctionBase&, const ShapeFunctionBase<Space>&)
      -> Dot<VectorFunctionBase, ShapeFunctionBase<Space>>;
   template <ShapeFunctionSpaceType Space>
   Dot(const ShapeFunctionBase<Space>&, const VectorFunctionBase&)
      -> Dot<VectorFunctionBase, ShapeFunctionBase<Space>>;

   /**
    * @brief Dot product between vector value ShapeFunctionBase instances.
    */
   template <>
   class Dot<
      ShapeFunctionBase<TrialSpace>,
      ShapeFunctionBase<TestSpace>>
      : public FormLanguage::Base
   {
      public:
         Dot(const ShapeFunctionBase<TrialSpace>& lhs,
               const ShapeFunctionBase<TestSpace>& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         Dot(const Dot& other)
            : m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Dot(Dot&& other)
            : m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
         {}

         ShapeFunctionBase<TrialSpace>& getLHS()
         {
            assert(m_lhs);
            return *m_lhs;
         }

         ShapeFunctionBase<TestSpace>& getRHS()
         {
            assert(m_rhs);
            return *m_rhs;
         }

         const ShapeFunctionBase<TrialSpace>& getLHS() const
         {
            assert(m_lhs);
            return *m_lhs;
         }

         const ShapeFunctionBase<TestSpace>& getRHS() const
         {
            assert(m_rhs);
            return *m_rhs;
         }

         mfem::DenseMatrix getElementMatrix(
               const mfem::FiniteElement& trialElement, const mfem::FiniteElement& testElement,
               mfem::ElementTransformation& trans) const
         {
            auto& trial = getLHS();
            auto& test  = getRHS();
            assert(trial.getRows(trialElement, trans) == test.getRows(testElement, trans));
            assert(trial.getColumns(trialElement, trans) == test.getColumns(testElement, trans));
            return test.getOperator(testElement, trans)->OperatorDot(
                     *trial.getOperator(trialElement, trans));
         }

         Dot* copy() const noexcept override
         {
            return new Dot(*this);
         }

      private:
         std::unique_ptr<ShapeFunctionBase<TrialSpace>> m_lhs;
         std::unique_ptr<ShapeFunctionBase<TestSpace>>  m_rhs;
   };
   Dot(const ShapeFunctionBase<TrialSpace>&,
         const ShapeFunctionBase<TestSpace>&)
      -> Dot<
         ShapeFunctionBase<TrialSpace>,
         ShapeFunctionBase<TestSpace>>;
}

#endif
