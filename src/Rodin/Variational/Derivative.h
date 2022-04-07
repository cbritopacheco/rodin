#ifndef RODIN_VARIATIONAL_DERIVATIVE_H
#define RODIN_VARIATIONAL_DERIVATIVE_H

#include <cassert>
#include <cstdlib>

#include "ForwardDecls.h"
#include "H1.h"
#include "FiniteElementSpace.h"
#include "GridFunction.h"
#include "ScalarFunction.h"

namespace Rodin::Variational
{
   /**
    * @brief Function representing the Derivative of a GridFunction in H1 space
    *
    * Given a GridFunction @f$ u : \mathbb{R}^s \rightarrow \mathbb{R}^d @f$,
    * this class represents the derivative in the @f$ i @f$-th direction of the
    * @f$ j @f$-th component:
    * @f[
    *    \dfrac{\partial u_j}{ \partial x_i } ,
    * @f]
    * where @f$ \quad 0 \leq i < s @f$ and @f$ \ 0 \leq j < d @f$.
    */
   class Derivative : public ScalarFunctionBase
   {
      public:
         /**
          * @brief Constructs the derivative in the i-th direction of the j-th
          * component
          * @param[in] direction Spatial direction @f$ x_i @f$
          * @param[in] component Component @f$ u_j @f$ to differentiate
          * @param[in] u GridFunction in H1 space
          */
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

   /**
    * @brief %Utility function for computing @f$ \partial_x u @f$
    * @param[in] u GridFunction in H1 space
    *
    * Given a scalar function @f$ u : \mathbb{R}^s \rightarrow \mathbb{R} @f$,
    * this function constructs the derivative in the @f$ x @f$ direction
    * @f$
    *    \dfrac{\partial u}{\partial x}
    * @f$
    */
   Derivative Dx(GridFunction<H1>& u);

   /**
    * @brief %Utility function for computing @f$ \partial_y u @f$
    * @param[in] u GridFunction in H1 space
    *
    * Given a scalar function @f$ u : \mathbb{R}^s \rightarrow \mathbb{R} @f$,
    * this function constructs the derivative in the @f$ y @f$ direction
    * @f$
    *    \dfrac{\partial u}{\partial y}
    * @f$
    */
   Derivative Dy(GridFunction<H1>& u);

   /**
    * @brief %Utility function for computing @f$ \partial_z u @f$
    * @param[in] u GridFunction in H1 space
    *
    * Given a scalar function @f$ u : \mathbb{R}^s \rightarrow \mathbb{R} @f$,
    * this function constructs the derivative in the @f$ y @f$ direction
    * @f$
    *    \dfrac{\partial u}{\partial z}
    * @f$
    */
   Derivative Dz(GridFunction<H1>& u);
}

#endif
