/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_INTEGRAL_H
#define RODIN_VARIATIONAL_INTEGRAL_H

#include <set>
#include <utility>
#include <mfem.hpp>

#include "Dot.h"
#include "ForwardDecls.h"
#include "FormLanguage.h"
#include "GridFunction.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "ScalarCoefficient.h"
#include "LinearFormIntegrator.h"
#include "BilinearFormIntegrator.h"

namespace Rodin::Variational
{
   /**
    * @brief Integral of the dot product of a trial function and a test function
    *
    * @f[
    *    \int_\Omega A(u) : A(v) \ dx
    * @f]
    */
   template <>
   class Integral<Dot<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>>
      : public BilinearFormDomainIntegrator
   {
      public:
         using Integrand = Dot<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>;

         Integral(const ShapeFunctionBase<Trial>& lhs, const ShapeFunctionBase<Test>& rhs)
            : Integral(Dot(lhs, rhs))
         {}

         Integral(const Integrand& prod)
            : m_prod(prod),
              m_intOrder(
                    [](const mfem::FiniteElement& trial,
                       const mfem::FiniteElement& test,
                       mfem::ElementTransformation& trans)
                    { return trial.GetOrder() + test.GetOrder() + trans.OrderW(); })
         {}

         Integral(const Integral& other)
            : BilinearFormDomainIntegrator(other),
              m_prod(other.m_prod),
              m_intOrder(other.m_intOrder)
         {}

         Integral(Integral&& other)
            : BilinearFormDomainIntegrator(std::move(other)),
              m_prod(std::move(other.m_prod)),
              m_intOrder(std::move(other.m_intOrder))
         {}

         Integral& setIntegrationOrder(
            std::function<
               int(const mfem::FiniteElement&, const mfem::FiniteElement&, mfem::ElementTransformation&)> order)
         {
            m_intOrder = order;
            return *this;
         }

         int getIntegrationOrder(
               const mfem::FiniteElement& trial, const mfem::FiniteElement& test, mfem::ElementTransformation& trans) const
         {
            return m_intOrder(trial, test, trans);
         }

         void getElementMatrix(
               const mfem::FiniteElement& trial, const mfem::FiniteElement& test,
               mfem::ElementTransformation& trans, mfem::DenseMatrix& mat) const override;

         Integral* copy() const noexcept override
         {
            return new Integral(*this);
         }
      private:
         Integrand m_prod;
         std::function<int(
            const mfem::FiniteElement&, const mfem::FiniteElement&, mfem::ElementTransformation&)> m_intOrder;
   };
   Integral(const Dot<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>&)
      -> Integral<Dot<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>>;
   Integral(const ShapeFunctionBase<Trial>& lhs, const ShapeFunctionBase<Test>& rhs)
      -> Integral<Dot<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>>;


   /**
    * @brief Integral of a test function
    *
    * @f[
    *    \int_\Omega A(v) \ dx
    * @f]
    */
   template <>
   class Integral<ShapeFunctionBase<Test>>
      : public LinearFormDomainIntegrator
   {
      public:
         using Integrand = ShapeFunctionBase<Test>;

         /**
          * @f[
          *    \int \lambda A(v) \ dx
          * @f]
          */
         Integral(const ScalarCoefficientBase& lhs, const ShapeFunctionBase<Test>& rhs)
            : Integral(Dot(lhs, rhs))
         {}

         /**
          * @f[
          *    \int \vec{\lambda} \cdot \vec{A}(v) \ dx
          * @f]
          */
         Integral(const VectorCoefficientBase& lhs, const ShapeFunctionBase<Test>& rhs)
            : Integral(Dot(lhs, rhs))
         {}

         /**
          * @f[
          *    \int A(v) \ dx
          * @f]
          */
         Integral(const Integrand& integrand)
            : m_test(integrand.copy()),
              m_intOrder(
                    [](const mfem::FiniteElement& fe,
                       mfem::ElementTransformation& trans)
                    { return fe.GetOrder() + trans.OrderW(); })
         {}

         Integral(const Integral& other)
            : LinearFormDomainIntegrator(other),
              m_test(other.m_test->copy()),
              m_intOrder(other.m_intOrder)
         {}

         Integral(Integral&& other)
            : LinearFormDomainIntegrator(std::move(other)),
              m_test(std::move(other.m_test))
         {}

         Integral& setIntegrationOrder(
            std::function<
               int(const mfem::FiniteElement&, mfem::ElementTransformation&)> order)
         {
            m_intOrder = order;
            return *this;
         }

         int getIntegrationOrder(
               const mfem::FiniteElement& fe, mfem::ElementTransformation& trans) const
         {
            return m_intOrder(fe, trans);
         }

         void getElementVector(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans, mfem::Vector& vec) const override;

         Integral* copy() const noexcept override
         {
            return new Integral(*this);
         }
      private:
         std::unique_ptr<Integrand> m_test;
         std::function<int(
            const mfem::FiniteElement&, mfem::ElementTransformation&)> m_intOrder;
   };
   Integral(const ShapeFunctionBase<Test>&)
      -> Integral<ShapeFunctionBase<Test>>;
   Integral(const ScalarCoefficientBase&, const ShapeFunctionBase<Test>&)
      -> Integral<ShapeFunctionBase<Test>>;
   Integral(const VectorCoefficientBase&, const ShapeFunctionBase<Test>&)
      -> Integral<ShapeFunctionBase<Test>>;

   template <>
   class BoundaryIntegral<ShapeFunctionBase<Test>> : public LinearFormBoundaryIntegrator
   {
      public:
         using Integrand = ShapeFunctionBase<Test>;

         template <class Lhs>
         BoundaryIntegral(Lhs&& lhs, const ShapeFunctionBase<Test>& rhs)
            : m_integral(std::forward<Lhs>(lhs), rhs)
         {}

         BoundaryIntegral(const Integrand& integrand)
            : m_integral(integrand)
         {}

         BoundaryIntegral(const BoundaryIntegral& other)
            : LinearFormBoundaryIntegrator(other),
              m_integral(other.m_integral)
         {}

         BoundaryIntegral(BoundaryIntegral&& other)
            : LinearFormBoundaryIntegrator(std::move(other)),
              m_integral(std::move(other.m_integral))
         {}

         void getElementVector(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans, mfem::Vector& vec) const override
         {
            m_integral.getElementVector(fe, trans, vec);
         }

         BoundaryIntegral* copy() const noexcept override
         {
            return new BoundaryIntegral(*this);
         }
      private:
         Integral<Integrand> m_integral;
   };
   BoundaryIntegral(const ShapeFunctionBase<Test>&)
      -> BoundaryIntegral<ShapeFunctionBase<Test>>;

   template <class Lhs>
   BoundaryIntegral(Lhs&& lhs, const ShapeFunctionBase<Test>& rhs)
      -> BoundaryIntegral<ShapeFunctionBase<Test>>;
}

#endif
