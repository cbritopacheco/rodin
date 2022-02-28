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
#include "Rodin/Variational/GridFunction.h"
#include "Rodin/Variational/ScalarCoefficient.h"
#include "Rodin/Variational/VectorCoefficient.h"
#include "Rodin/Variational/MatrixCoefficient.h"
#include "Rodin/Variational/TestFunction.h"
#include "Rodin/Variational/TrialFunction.h"
#include "FormLanguage/Base.h"
#include "ForwardDecls.h"

namespace Rodin::Variational
{
   // /**
   //  * @brief Dot product between instances of Lhs and Rhs
   //  * @tparam Lhs Left-hand side operand type
   //  * @tparam Rhs Right-hand side operand type
   //  */
   // template <class Lhs, class Rhs>
   // class Dot : public FormLanguage::Base
   // {
   //    static_assert(std::is_base_of_v<Base, Lhs>,
   //          "Lhs must be derived from FormLanguage::Base");
   //    static_assert(std::is_base_of_v<Base, Rhs>,
   //          "Rhs must be derived from FormLanguage::Base");

   //    public:
   //       Dot(const Lhs& lhs, const Rhs& rhs)
   //          : m_lhs(lhs.copy()), m_rhs(rhs.copy())
   //       {}

   //       Dot(const Dot& other)
   //          : m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
   //       {}

   //       Dot(Dot&&) = default;

   //       virtual ~Dot() = default;

   //       Lhs& getLHS()
   //       {
   //          assert(m_lhs);
   //          return *m_lhs;
   //       }

   //       Rhs& getRHS()
   //       {
   //          assert(m_rhs);
   //          return *m_rhs;
   //       }

   //       const Lhs& getLHS() const
   //       {
   //          assert(m_lhs);
   //          return *m_lhs;
   //       }

   //       const Rhs& getRHS() const
   //       {
   //          assert(m_lhs);
   //          return *m_rhs;
   //       }

   //       Dot* copy() const noexcept override
   //       {
   //          return new Dot(*this);
   //       }

   //    private:
   //       std::unique_ptr<Lhs> m_lhs;
   //       std::unique_ptr<Rhs> m_rhs;
   // };

   /**
    * @brief Represents the dot product between two matrices.
    *
    * For two given @f$ n \times m @f$ matrices @f$ A @f$ and @f$ B @f$, their
    * double-dot product @f$ A \colon B @f$ is defined as
    * @f[
    *    A \colon B = \sum_{i = 1}^n \sum_{j = 1}^m A_{ij} B_{ij} = \mathrm{tr}(B^T A)
    * @f]
    *
    * @tparam A Derived type from MatrixCoefficientBase
    * @tparam B Derived type from MatrixCoefficientBase
    */
   template <>
   class Dot<MatrixCoefficientBase, MatrixCoefficientBase> : public ScalarCoefficientBase
   {
      public:
         /**
          * @brief Constructs the Dot product between two given matrices.
          * @param[in] a Derived instance of MatrixCoefficientBase
          * @param[in] b Derived instance of MatrixCoefficientBase
          */
         Dot(const MatrixCoefficientBase& a, const MatrixCoefficientBase& b);

         Dot(const Dot& other);

         double getValue(
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const override;

         Dot* copy() const noexcept override
         {
            return new Dot(*this);
         }
      private:
         std::unique_ptr<MatrixCoefficientBase> m_a, m_b;
   };
   Dot(const MatrixCoefficientBase&, const MatrixCoefficientBase&)
      -> Dot<MatrixCoefficientBase, MatrixCoefficientBase>;

   template <ShapeFunctionSpaceType Space>
   class Dot<ScalarCoefficientBase, ShapeFunctionBase<Space>>
      : public ShapeFunctionBase<Space>
   {
      public:
         Dot(const ScalarCoefficientBase& lhs, const ShapeFunctionBase<Space>& rhs)
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

         ScalarCoefficientBase& getLHS()
         {
            return *m_lhs;
         }

         ShapeFunctionBase<Space>& getRHS()
         {
            return *m_rhs;
         }

         const ScalarCoefficientBase& getLHS() const
         {
            return *m_lhs;
         }

         const ShapeFunctionBase<Space>& getRHS() const
         {
            return *m_rhs;
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
               const mfem::ElementTransformation& trans) const
         {
            return getRHS().getDOFs(fe, trans);
         }

         std::unique_ptr<Internal::Rank3OperatorBase> getOperator(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans) const override
         {
            assert(getRows(fe, trans) == 1);
            assert(getColumns(fe, trans) == 1);

            auto result = getRHS().getOperator(fe, trans);
            (*result) *= getLHS().getValue(trans, trans.GetIntPoint());
            return std::move(result);
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
         std::unique_ptr<ScalarCoefficientBase> m_lhs;
         std::unique_ptr<ShapeFunctionBase<Space>> m_rhs;
   };
   template <ShapeFunctionSpaceType Space>
   Dot<ScalarCoefficientBase, ShapeFunctionBase<Space>>
   operator*(const ScalarCoefficientBase& lhs, const ShapeFunctionBase<Space>& rhs)
   {
      return Dot(lhs, rhs);
   }
   template <ShapeFunctionSpaceType Space>
   Dot(const ScalarCoefficientBase&, const ShapeFunctionBase<Space>&)
      -> Dot<ScalarCoefficientBase, ShapeFunctionBase<Space>>;

   /**
    * @brief Dot product between VectorCoefficientBase and ShapeFunctionBase.
    */
   template <ShapeFunctionSpaceType Space>
   class Dot<VectorCoefficientBase, ShapeFunctionBase<Space>>
      : public ShapeFunctionBase<Space>
   {
      public:
         Dot(const VectorCoefficientBase& lhs, const ShapeFunctionBase<Space>& rhs)
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

         VectorCoefficientBase& getLHS()
         {
            return *m_lhs;
         }

         ShapeFunctionBase<Space>& getRHS()
         {
            return *m_rhs;
         }

         const VectorCoefficientBase& getLHS() const
         {
            return *m_lhs;
         }

         const ShapeFunctionBase<Space>& getRHS() const
         {
            return *m_rhs;
         }

         int getRows(
               const mfem::FiniteElement&,
               const mfem::ElementTransformation&) const override
         {
            return 1;
         }

         int getDOFs(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const
         {
            return getRHS().getDOFs(fe, trans);
         }

         int getColumns(
               const mfem::FiniteElement&,
               const mfem::ElementTransformation&) const override
         {
            return 1;
         }

         std::unique_ptr<Internal::Rank3OperatorBase> getOperator(
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
         std::unique_ptr<VectorCoefficientBase> m_lhs;
         std::unique_ptr<ShapeFunctionBase<Space>> m_rhs;
   };
   template <ShapeFunctionSpaceType Space>
   Dot(const VectorCoefficientBase&, const ShapeFunctionBase<Space>&)
      -> Dot<VectorCoefficientBase, ShapeFunctionBase<Space>>;
   template <ShapeFunctionSpaceType Space>
   Dot<VectorCoefficientBase, ShapeFunctionBase<Space>>
   operator*(const VectorCoefficientBase& lhs, const ShapeFunctionBase<Space>& rhs)
   {
      return Dot(lhs, rhs);
   }

   /**
    * @brief Dot product between instances of Lhs and Rhs
    * @tparam Lhs Left-hand side operand type
    * @tparam Rhs Right-hand side operand type
    */
   template <>
   class Dot<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>
      : public FormLanguage::Base
   {
      public:
         Dot(const ShapeFunctionBase<Trial>& lhs, const ShapeFunctionBase<Test>& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         Dot(const Dot& other)
            : m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Dot(Dot&&) = default;

         ShapeFunctionBase<Trial>& getLHS()
         {
            assert(m_lhs);
            return *m_lhs;
         }

         ShapeFunctionBase<Test>& getRHS()
         {
            assert(m_rhs);
            return *m_rhs;
         }

         const ShapeFunctionBase<Trial>& getLHS() const
         {
            assert(m_lhs);
            return *m_lhs;
         }

         const ShapeFunctionBase<Test>& getRHS() const
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
            return test.getOperator(trialElement, trans)->OperatorDot(
                     *trial.getOperator(testElement, trans));
         }

         Dot* copy() const noexcept override
         {
            return new Dot(*this);
         }

      private:
         std::unique_ptr<ShapeFunctionBase<Trial>> m_lhs;
         std::unique_ptr<ShapeFunctionBase<Test>>  m_rhs;
   };
   Dot(const ShapeFunctionBase<Trial>&, const ShapeFunctionBase<Test>&)
      -> Dot<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>;
   Dot<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>
   operator*(const ShapeFunctionBase<Trial>& lhs, const ShapeFunctionBase<Test>& rhs);
}

#endif
