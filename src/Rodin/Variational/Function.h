/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FUNCTION_H
#define RODIN_VARIATIONAL_FUNCTION_H

#include <set>
#include <variant>
#include <type_traits>

#include <mfem.hpp>

#include "Rodin/FormLanguage/Base.h"

#include "ForwardDecls.h"
#include "Rodin/Mesh/Element.h"

namespace Rodin::Variational::Internal
{
   template <RangeType R>
   class ProxyFunction;

   using ScalarProxyFunction = ProxyFunction<RangeType::Scalar>;
   using VectorProxyFunction = ProxyFunction<RangeType::Vector>;
   using MatrixProxyFunction = ProxyFunction<RangeType::Matrix>;

   using MFEMFunction =
      std::variant<
         Internal::ScalarProxyFunction,
         Internal::VectorProxyFunction,
         Internal::MatrixProxyFunction>;
}

namespace Rodin::Variational
{
   class FunctionBase : public FormLanguage::Base
   {
      public:
         FunctionBase() = default;

         FunctionBase(const FunctionBase& other)
            : FormLanguage::Base(other),
              m_traceDomain(other.m_traceDomain)
         {}

         FunctionBase(FunctionBase&& other)
            : FormLanguage::Base(std::move(other)),
              m_traceDomain(std::move(other.m_traceDomain))
         {}

         FunctionBase& operator=(FunctionBase&& other)
         {
            m_traceDomain = std::move(other.m_traceDomain);
            return *this;
         }

         /**
          * @brief Sets an attribute which will be interpreted as the domain to
          * trace.
          *
          * Convenience function to call traceOf(std::set<int>) with only one
          * attribute.
          *
          * @returns Reference to self (for method chaining)
          */
         virtual FunctionBase& traceOf(int attr)
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
          * respective "interior" domain, since it is not clear which domain
          * attribute can be used to extend the value continuously up to the
          * boundary. To resolve this ambiguity the trace domain is interpreted
          * as the domain which shall be used to make this continuous
          * extension.
          *
          * @note Setting the trace domain of a FunctionBase instance
          * does not guarantee that it will taken into consideration when
          * computing its value. That said, it is up to the subclass to decide
          * how it will use this information which can be obtained via the
          * getTraceDomain() method.
          *
          * @see @ref getTraceDomain() "getTraceDomain()"
          */
         virtual FunctionBase& traceOf(const std::set<int>& attrs)
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

         virtual Transpose<FunctionBase> T() const;

         virtual RangeShape getRangeShape() const = 0;

         virtual RangeType getRangeType() const;

         /**
          * @note It is not necessary to set the size beforehand.
          */
         virtual void getValue(
               mfem::DenseMatrix& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const = 0;

         mfem::DenseMatrix operator()(const Vertex& v) const
         {
            mfem::DenseMatrix m;
            getValue(m, *v.getElementTransformation(), *v.getIntegrationPoint());
            return m;
         }

         virtual FunctionBase* copy() const noexcept override = 0;

         Internal::MFEMFunction build() const;

      protected:
         mfem::ElementTransformation& getTraceElementTrans(
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const;

      private:
         std::set<int> m_traceDomain;
   };
}

namespace Rodin::Variational::Internal
{
   template <>
   class ProxyFunction<RangeType::Scalar> : public mfem::Coefficient
   {
      public:
         ProxyFunction(const FunctionBase& s);

         ProxyFunction(const ProxyFunction& other);

         ProxyFunction(ProxyFunction&& other);

         double Eval(
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip) override;

      private:
         const FunctionBase& m_s;
   };

   template <>
   class ProxyFunction<RangeType::Vector> : public mfem::VectorCoefficient
   {
      public:
         ProxyFunction(const FunctionBase& s);

         ProxyFunction(const ProxyFunction& other);

         ProxyFunction(ProxyFunction&& other);

         void Eval(
               mfem::Vector& value,
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip) override;

      private:
         const FunctionBase& m_s;
   };

   template <>
   class ProxyFunction<RangeType::Matrix> : public mfem::MatrixCoefficient
   {
      public:
         ProxyFunction(const FunctionBase& s);

         ProxyFunction(const ProxyFunction& other);

         ProxyFunction(ProxyFunction&& other);

         void Eval(
               mfem::DenseMatrix& value,
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip) override;
      private:
         const FunctionBase& m_s;
   };
}

#endif
