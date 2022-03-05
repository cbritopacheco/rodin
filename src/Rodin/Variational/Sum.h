/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_SUM_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_SUM_H

#include <memory>
#include <type_traits>

#include "Rodin/Variational/ScalarCoefficient.h"
#include "Rodin/Variational/MatrixCoefficient.h"
#include "FormLanguage/Base.h"
#include "ForwardDecls.h"

namespace Rodin::Variational
{
   template <class Lhs, class Rhs>
   class Sum : public FormLanguage::Base
   {
      static_assert(std::is_base_of_v<Base, Lhs>,
            "Lhs must be derived from FormLanguage::Base");
      static_assert(std::is_base_of_v<Base, Rhs>,
            "Rhs must be derived from FormLanguage::Base");

      public:
         Sum(const Lhs& lhs, const Rhs& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         Sum(const Sum& other)
            :  Base(other),
               m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Sum(Sum&& other)
            :  Base(std::move(other)),
               m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
         {}

         Lhs& getLHS()
         {
            return *m_lhs;
         }

         Rhs& getRHS()
         {
            return *m_rhs;
         }

         const Lhs& getLHS() const
         {
            return *m_lhs;
         }

         const Rhs& getRHS() const
         {
            return *m_rhs;
         }

         Sum* copy() const noexcept override
         {
            return new Sum(*this);
         }

      private:
         std::unique_ptr<Lhs> m_lhs;
         std::unique_ptr<Rhs> m_rhs;
   };

   template <>
   class Sum<ScalarCoefficientBase, ScalarCoefficientBase>
      : public ScalarCoefficientBase
   {
      public:
         Sum(const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         Sum(const Sum& other)
            :  ScalarCoefficientBase(other),
               m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Sum(Sum&& other)
            :  ScalarCoefficientBase(std::move(other)),
               m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
         {}

         ScalarCoefficientBase& getLHS()
         {
            return *m_lhs;
         }

         ScalarCoefficientBase& getRHS()
         {
            return *m_rhs;
         }

         const ScalarCoefficientBase& getLHS() const
         {
            return *m_lhs;
         }

         const ScalarCoefficientBase& getRHS() const
         {
            return *m_rhs;
         }

         double getValue(
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const override;

         Sum* copy() const noexcept override
         {
            return new Sum(*this);
         }

      private:
         std::unique_ptr<ScalarCoefficientBase> m_lhs;
         std::unique_ptr<ScalarCoefficientBase> m_rhs;
   };
   Sum<ScalarCoefficientBase, ScalarCoefficientBase>
   operator+(const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs);

   template <>
   class Sum<MatrixCoefficientBase, MatrixCoefficientBase>
      : public MatrixCoefficientBase
   {
      public:
         Sum(const MatrixCoefficientBase& lhs, const MatrixCoefficientBase& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         Sum(const Sum& other)
            :  MatrixCoefficientBase(other),
               m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Sum(Sum&& other)
            :  MatrixCoefficientBase(std::move(other)),
               m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
         {}

         MatrixCoefficientBase& getLHS()
         {
            return *m_lhs;
         }

         MatrixCoefficientBase& getRHS()
         {
            return *m_rhs;
         }

         const MatrixCoefficientBase& getLHS() const
         {
            return *m_lhs;
         }

         const MatrixCoefficientBase& getRHS() const
         {
            return *m_rhs;
         }

         int getRows() const override;

         int getColumns() const override;

         void getValue(
               mfem::DenseMatrix& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) override;

         Sum* copy() const noexcept override
         {
            return new Sum(*this);
         }

      private:
         std::unique_ptr<MatrixCoefficientBase> m_lhs;
         std::unique_ptr<MatrixCoefficientBase> m_rhs;
   };
   Sum<MatrixCoefficientBase, MatrixCoefficientBase>
   operator+(const MatrixCoefficientBase& lhs, const MatrixCoefficientBase& rhs);
}

#endif
