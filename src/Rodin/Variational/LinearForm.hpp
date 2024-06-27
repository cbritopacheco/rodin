/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_LINEARFORM_HPP
#define RODIN_VARIATIONAL_LINEARFORM_HPP

#include "Rodin/Alert.h"
#include "Rodin/Assembly/AssemblyBase.h"

#include "LinearForm.h"
#include "LinearFormIntegrator.h"
#include "FiniteElementSpace.h"


namespace Rodin::Variational
{
   template <class FES>
   void
   LinearForm<FES, Math::Vector<Scalar>>::assemble()
   {
      const auto& fes = getTestFunction().getFiniteElementSpace();
      m_vector = getAssembly().execute({ fes, getIntegrators() });
   }
}

#endif
