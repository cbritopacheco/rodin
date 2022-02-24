#ifndef RODIN_VARIATIONAL_INTEGRAL_H
#define RODIN_VARIATIONAL_INTEGRAL_H

#include <set>
#include <utility>
#include <mfem.hpp>

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
   class Integral<FormLanguage::Product<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>>
      : public BilinearFormDomainIntegrator
   {
      public:
         using Integrand = FormLanguage::Product<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>;

         Integral(const Integrand& prod)
            : m_prod(prod)
         {}

         Integral(const Integral& other)
            : BilinearFormDomainIntegrator(other),
              m_prod(other.m_prod)
         {}

         Integral(Integral&& other)
            : BilinearFormDomainIntegrator(std::move(other)),
              m_prod(std::move(other.m_prod))
         {}

         void getElementMatrix(
               const mfem::FiniteElement& trial, const mfem::FiniteElement& test,
               mfem::ElementTransformation& trans, mfem::DenseMatrix& mat) override;

         Integral* copy() const noexcept override
         {
            return new Integral(*this);
         }
      private:
         Integrand m_prod;
   };
   Integral(const FormLanguage::Product<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>&)
      -> Integral<FormLanguage::Product<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>>;

   template <>
   class Integral<ShapeFunctionBase<Test>>
      : public LinearFormDomainIntegrator
   {
      public:
         using Integrand = ShapeFunctionBase<Test>;

         Integral(const Integrand& test)
            : m_test(test.copy())
         {}

         Integral(const Integral& other)
            : LinearFormDomainIntegrator(other),
              m_test(other.m_test->copy())
         {}

         Integral(Integral&& other)
            : LinearFormDomainIntegrator(std::move(other)),
              m_test(std::move(other.m_test))
         {}

         void getElementVector(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans, mfem::Vector& vec) override;

         Integral* copy() const noexcept override
         {
            return new Integral(*this);
         }
      private:
         std::unique_ptr<Integrand> m_test;
   };
   Integral(const ShapeFunctionBase<Test>&) -> Integral<ShapeFunctionBase<Test>>;
}

#endif
