/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DERIVATIVE_HPP
#define RODIN_VARIATIONAL_DERIVATIVE_HPP

#include "Derivative.h"

namespace Rodin::Variational
{
   template <int ComponentIndex, int VariableIndex>
   Derivative<ComponentIndex, VariableIndex>::Derivative(GridFunctionBase& u)
      : m_u(u)
   {}

   template <int ComponentIndex, int VariableIndex>
   void Derivative<ComponentIndex, VariableIndex>::buildMFEMCoefficient()
   {
   }

   template <int ComponentIndex, int VariableIndex>
   mfem::Coefficient& Derivative<ComponentIndex, VariableIndex>::getMFEMCoefficient()
   {
   }
}

#endif

