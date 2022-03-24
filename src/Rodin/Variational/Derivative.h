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
               const mfem::IntegrationPoint&) const override
         {
            mfem::DenseMatrix grad;
            m_u.getHandle().GetVectorGradient(trans, grad);
            return grad(m_component, m_direction);
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

   Derivative Dx(GridFunction<H1>& u);
   Derivative Dy(GridFunction<H1>& u);
   Derivative Dz(GridFunction<H1>& u);
}

#endif
