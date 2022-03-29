/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_MATRIXCOEFFICIENT_H
#define RODIN_VARIATIONAL_MATRIXCOEFFICIENT_H

#include <set>
#include <optional>
#include <mfem.hpp>

#include "Rodin/Alert.h"
#include "ForwardDecls.h"
#include "FormLanguage/Base.h"

namespace Rodin::Variational
{
   /**
    * @brief Abstract base class for objects representing matrix coefficients.
    */
   class MatrixCoefficientBase
      : public FormLanguage::Buildable<mfem::MatrixCoefficient>
   {
      public:
         MatrixCoefficientBase() = default;

         MatrixCoefficientBase(const MatrixCoefficientBase& other)
            : m_traceDomain(other.m_traceDomain)
         {}

         MatrixCoefficientBase(MatrixCoefficientBase&& other)
            : m_traceDomain(std::move(other.m_traceDomain))
         {}

         /**
          * @brief Sets an attribute which will be interpreted as the domain to
          * trace.
          *
          * Convenience function to call traceOf(std::set<int>) with only one
          * attribute.
          *
          * @returns Reference to self (for method chaining)
          */
         MatrixCoefficientBase& traceOf(int attr)
         {
            return traceOf(std::set<int>{attr});
         }

         /**
          * @brief Sets which attributes will be interpreted as the domain to
          * trace.
          * @returns Reference to self (for method chaining)
          *
          * When integrating along interior boundaries sometimes it is
          * necessary to specify which attributes should be interpreted as the
          * respective "interior" domain. For example, coefficients which
          * involve the derivatives of a GridFunction need to know the element
          * to "trace".
          *
          * @note Setting the trace domain of a MatrixCoefficientBase instance
          * does not guarantee that it will taken into consideration when
          * computing its value. That said, it is up to the subclass to decide
          * how it will use this information which can be obtained via the
          * getTraceDomain() method.
          *
          * @see @ref MatrixCoefficientBase::getTraceDomain() "getTraceDomain()"
          *
          */
         MatrixCoefficientBase& traceOf(std::set<int> attrs)
         {
            m_traceDomain = attrs;
            return *this;
         }

         /**
          * @brief Gets the set of attributes which will be interpreted as the
          * domains to "trace".
          *
          * The domains to trace are interpreted as the domains where there
          * shall be a continuous extension from values to the interior
          * boundaries. If the trace domain is empty, then this has the
          * semantic value that it has not been specified yet.
          */
         const std::set<int>& getTraceDomain() const
         {
            return m_traceDomain;
         }

         std::unique_ptr<mfem::MatrixCoefficient> build() const override;

         virtual ~MatrixCoefficientBase() = default;

         virtual Transpose<MatrixCoefficientBase> T() const;

         virtual int getRows() const = 0;

         virtual int getColumns() const = 0;

         virtual void getValue(
               mfem::DenseMatrix& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const = 0;

         virtual MatrixCoefficientBase* copy() const noexcept override = 0;

      private:
         std::set<int> m_traceDomain;
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
