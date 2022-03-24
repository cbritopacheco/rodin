#ifndef RODIN_VARIATIONAL_DERIVATIVE_H
#define RODIN_VARIATIONAL_DERIVATIVE_H

#include <cassert>
#include <cstdlib>

#include "ForwardDecls.h"
#include "H1.h"
#include "GridFunction.h"
#include "ScalarCoefficient.h"

namespace Rodin::Variational
{
   /**
    * @brief Derivative
    */
   class Derivative : public ScalarCoefficientBase
   {
      public:
         Derivative(int direction, int component, GridFunction<H1>& u)
            :  m_direction(direction),
               m_component(component),
               m_u(u)
         {}

         Derivative(const Derivative& other) = default;

         Derivative(Derivative&& other) = default;

         double getValue(
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip) const override
         {
            switch (trans.ElementType)
            {
               case mfem::ElementTransformation::ELEMENT:
               {
                  const mfem::FiniteElement* fe =
                     m_u.getFiniteElementSpace().getFES().GetFE(trans.ElementNo);
                  const mfem::IntegrationRule& ir = fe->GetNodes();
                  mfem::DenseMatrix dshape(fe->GetDof(), trans.GetSpaceDim());
                  fe->CalcPhysDShape(trans, dshape);
                  assert(false);
                  break;
               }
               case mfem::ElementTransformation::BDR_ELEMENT:
               {
                  break;
               }
               default:
                  Alert::Exception("Unhandled element type").raise();
            }
            return NAN;
         }

         Derivative* copy() const noexcept override
         {
            return new Derivative(*this);
         }
      private:
         int m_direction;
         int m_component;
         GridFunction<H1>& m_u;
   };
}

#endif
