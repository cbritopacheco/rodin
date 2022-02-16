/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DERIVATIVE_HPP
#define RODIN_VARIATIONAL_DERIVATIVE_HPP

#include "GridFunction.h"

#include "Derivative.h"

namespace Rodin::Variational
{
   template <int DirectionIndex, int ComponentIndex>
   Derivative<DirectionIndex, ComponentIndex>::Derivative(GridFunctionBase& u)
      : m_u(u)
   {}

   template <int DirectionIndex, int ComponentIndex>
   void Derivative<DirectionIndex, ComponentIndex>::build()
   {
      m_mfemGridFunction.emplace(m_u.getHandle().FESpace());
      m_u.getHandle().GetDerivative(
            ComponentIndex, DirectionIndex - 1, *m_mfemGridFunction);
      m_mfemCoefficient.emplace(&(*m_mfemGridFunction));
   }

   template <int DirectionIndex, int ComponentIndex>
   mfem::Coefficient& Derivative<DirectionIndex, ComponentIndex>::get()
   {
      assert(m_mfemCoefficient);
      return *m_mfemCoefficient;
   }
}

#endif

