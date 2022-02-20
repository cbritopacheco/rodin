#ifndef RODIN_VARIATIONAL_INTEGRAL_H
#define RODIN_VARIATIONAL_INTEGRAL_H

#include <set>
#include <mfem.hpp>

#include "ForwardDecls.h"
#include "FormLanguage.h"
#include "TrialFunction.h"
#include "TestFunction.h"
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
         Integral(const FormLanguage::Product<TrialFunctionBase, TestFunctionBase>& prod);

         Integral(const Integral& other) = default;

         virtual void getElementMatrix(
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
}

#endif
