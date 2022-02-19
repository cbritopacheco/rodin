/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_MATRIXCOEFFICIENT_H
#define RODIN_VARIATIONAL_MATRIXCOEFFICIENT_H

#include <optional>
#include <mfem.hpp>

#include "ForwardDecls.h"
#include "FormLanguage/Base.h"

namespace Rodin::Variational
{
   namespace Internal
   {
      class MatrixCoefficient : public mfem::MatrixCoefficient
      {
         public:
            MatrixCoefficient(const MatrixCoefficientBase& mat);

            void Eval(
                  mfem::DenseMatrix& value,
                  mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) override;

         private:
            std::unique_ptr<MatrixCoefficientBase> m_mat;
      };
   }

   class MatrixCoefficientBase
      : public FormLanguage::Buildable<Internal::MatrixCoefficient>
   {
      public:
         constexpr
         MatrixCoefficientBase() = default;

         constexpr
         MatrixCoefficientBase(const MatrixCoefficientBase&) = default;

         constexpr
         MatrixCoefficientBase& setTraceDomain(int domain)
         {
            m_traceDomain = domain;
            return *this;
         }

         constexpr
         std::optional<int> getTraceDomain() const
         {
            return m_traceDomain;
         }

         virtual void getValueOnInteriorBoundary(
               mfem::DenseMatrix& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip);

         virtual void getValue(
               mfem::DenseMatrix& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) = 0;

         virtual ~MatrixCoefficientBase() = default;

         virtual Transpose T() const;

         std::unique_ptr<Internal::MatrixCoefficient> build() const override
         {
            return std::make_unique<Internal::MatrixCoefficient>(*this);
         }

         virtual int getRows() const = 0;

         virtual int getColumns() const = 0;

         virtual MatrixCoefficientBase* copy() const noexcept override = 0;

      private:
         std::optional<int> m_traceDomain;
   };
}

#endif
