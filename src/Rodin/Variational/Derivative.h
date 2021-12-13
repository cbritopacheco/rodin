#ifndef RODIN_VARIATIONAL_DERIVATIVE_H
#define RODIN_VARIATIONAL_DERIVATIVE_H

#include <cstddef>
#include <mfem.hpp>

#include "ForwardDecls.h"

namespace Rodin::Variational
{
   /**
    * @brief Represents the partial derivative with respect to the i-th
    * component and the j-th variable of the function
    * @f$ u \colon \mathbb{R}^n \rightarrow \mathbb{R}^m @f$.
    *
    * Denote
    * @f[
    *    u(x) = (u_1(x), \ldots, u_m(x))
    * @f]
    * the value of @f$ u @f$ at @f$ x = (x_1, \ldots, x_n) @f$. This class
    * represents the partial derivative with respect to the i-th component and
    * the j-th variable
    * @f[
    *    \dfrac{\partial u_i}{\partial x_j}
    * @f]
    *
    * @tparam i Index of the function component
    * @tparam j Index of the variable
    */
   template <size_t i, size_t j>
   class Derivative
   {
      private:
         GridFunctionBase& m_u;
         mfem::GridFunctionCoefficient m_mfemCoefficient;
   };
}

#endif
