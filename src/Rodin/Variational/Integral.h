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
   template <>
   class Integral<Dot<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>>
      : public BilinearFormDomainIntegrator
   {
      public:
         using Integrand = Dot<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>;

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
               const mfem::FiniteElement& trial, const mfem::FiniteElement& test, mfem::ElementTransformation& trans)
         {
            return m_intOrder(trial, test, trans);
         }

         void getElementMatrix(
               const mfem::FiniteElement& trial, const mfem::FiniteElement& test,
               mfem::ElementTransformation& trans, mfem::DenseMatrix& mat) override;

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

   template <>
   class Integral<ShapeFunctionBase<Test>>
      : public LinearFormDomainIntegrator
   {
      public:
         using Integrand = ShapeFunctionBase<Test>;

         Integral(const Integrand& test)
            : m_test(test.copy()),
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
               const mfem::FiniteElement& fe, mfem::ElementTransformation& trans)
         {
            return m_intOrder(fe, trans);
         }

         void getElementVector(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans, mfem::Vector& vec) override;

         Integral* copy() const noexcept override
         {
            return new Integral(*this);
         }
      private:
         std::unique_ptr<Integrand> m_test;
         std::function<int(
            const mfem::FiniteElement&, mfem::ElementTransformation&)> m_intOrder;
   };
   Integral(const ShapeFunctionBase<Test>&) -> Integral<ShapeFunctionBase<Test>>;
}

#endif
