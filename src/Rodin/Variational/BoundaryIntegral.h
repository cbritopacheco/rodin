#ifndef RODIN_VARIATIONAL_BOUNDARYINTEGRAL_H
#define RODIN_VARIATIONAL_BOUNDARYINTEGRAL_H

#include "ForwardDecls.h"
#include "Integral.h"
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
   class BoundaryIntegral<ShapeFunctionBase<Test>> : public LinearFormBoundaryIntegrator
   {
      public:
         using Integrand = ShapeFunctionBase<Test>;

         BoundaryIntegral(const Integrand& op)
            : m_integral(op)
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
               mfem::ElementTransformation& trans, mfem::Vector& vec) override
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
   BoundaryIntegral(const ShapeFunctionBase<Test>&) -> BoundaryIntegral<ShapeFunctionBase<Test>>;
}

#endif

