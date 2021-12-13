/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DERIVATIVE_H
#define RODIN_VARIATIONAL_DERIVATIVE_H

#include <cstddef>
#include <mfem.hpp>

#include "ForwardDecls.h"
#include "ScalarCoefficient.h"

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
    * @tparam ComponentIndex Index of the function component @f$ u @f$
    * @tparam VariableIndex Index of the variable argument @f$ x @f$
    */
   template <int ComponentIndex, int VariableIndex>
   class Derivative : public ScalarCoefficientBase
   {
      static_assert(ComponentIndex >= 0 && VariableIndex >= 0,
            "ComponentIndex and VariableIndex must both be non-negative.");

      public:

         Derivative(GridFunctionBase& u);

         void buildMFEMCoefficient() override;

         mfem::Coefficient& getMFEMCoefficient() override;

         Derivative* copy() const noexcept override
         {
            return new Derivative(*this);
         }

      private:
         GridFunctionBase& m_u;
   };
}

#include "Derivative.hpp"

#endif
