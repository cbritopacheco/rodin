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
   class MatrixCoefficientBase
      : public FormLanguage::Buildable<mfem::MatrixCoefficient>
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

         std::unique_ptr<mfem::MatrixCoefficient> build() const override;

         virtual ~MatrixCoefficientBase() = default;

         virtual Transpose T() const;

         virtual int getRows() const = 0;

         virtual int getColumns() const = 0;

         virtual void getValue(
               mfem::DenseMatrix& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const = 0;

         virtual MatrixCoefficientBase* copy() const noexcept override = 0;

      private:
         std::optional<int> m_traceDomain;
   };
}

namespace Rodin::Variational::Internal
{
   class ProxyMatrixCoefficient : public mfem::MatrixCoefficient
   {
      public:
         ProxyMatrixCoefficient(const MatrixCoefficientBase& mat)
            :  mfem::MatrixCoefficient(mat.getRows(), mat.getColumns()),
               m_mat(mat)
         {}

         ProxyMatrixCoefficient(const ProxyMatrixCoefficient& other)
            :  mfem::MatrixCoefficient(other),
               m_mat(other.m_mat)
         {}

         ProxyMatrixCoefficient(ProxyMatrixCoefficient&& other)
            :  mfem::MatrixCoefficient(std::move(other)),
               m_mat(other.m_mat)
         {}

         void Eval(
               mfem::DenseMatrix& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) override
         {
            m_mat.getValue(value, trans, ip);
         }

      private:
         const MatrixCoefficientBase& m_mat;
   };
}

#endif
