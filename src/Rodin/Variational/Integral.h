#ifndef RODIN_VARIATIONAL_INTEGRAL_H
#define RODIN_VARIATIONAL_INTEGRAL_H

#include <set>
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
   template <class T>
   class Integral;

   template <>
   class Integral<FormLanguage::Product<TrialFunctionBase, TestFunctionBase>>
      : public BilinearFormDomainIntegrator
   {
      public:
         Integral(const FormLanguage::Product<TrialFunctionBase, TestFunctionBase>& prod)
            : m_prod(prod)
         {}

         Integral(const Integral& other) = default;

         void getElementMatrix(
               const mfem::FiniteElement& trial, const mfem::FiniteElement& test,
               mfem::ElementTransformation& trans, mfem::DenseMatrix& mat) override;

         IntegratorRegion getIntegratorRegion() const override
         {
            return Domain;
         }

         Integral* copy() const noexcept override
         {
            return new Integral(*this);
         }
      private:
         FormLanguage::Product<TrialFunctionBase, TestFunctionBase> m_prod;
   };
   Integral(const FormLanguage::Product<TrialFunctionBase, TestFunctionBase>&)
      -> Integral<FormLanguage::Product<TrialFunctionBase, TestFunctionBase>>;

   template <>
   class Integral<FormLanguage::Product<ScalarCoefficientBase, TestFunctionBase>>
      : public LinearFormDomainIntegrator
   {
      public:
         Integral(const FormLanguage::Product<ScalarCoefficientBase, TestFunctionBase>& prod)
            : m_prod(prod)
         {}

         Integral(const Integral& other) = default;

         void getElementVector(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans, mfem::Vector& vec) override;

         IntegratorRegion getIntegratorRegion() const override
         {
            return Domain;
         }

         Integral* copy() const noexcept override
         {
            return new Integral(*this);
         }
      private:
         FormLanguage::Product<ScalarCoefficientBase, TestFunctionBase> m_prod;
   };
   Integral(const FormLanguage::Product<ScalarCoefficientBase, TestFunctionBase>&)
      -> Integral<FormLanguage::Product<ScalarCoefficientBase, TestFunctionBase>>;
}

#endif