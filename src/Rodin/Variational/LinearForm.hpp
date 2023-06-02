/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_LINEARFORM_HPP
#define RODIN_VARIATIONAL_LINEARFORM_HPP

#include "Rodin/Alert.h"

#include "Assembly/AssemblyBase.h"

#include "LinearForm.h"
#include "LinearFormIntegrator.h"
#include "FiniteElementSpace.h"


namespace Rodin::Variational
{
   template <class FES>
   void
   LinearForm<FES, Context::Serial, Math::Vector>::assemble()
   {
      const auto& fes = getTestFunction().getFiniteElementSpace();
      const auto& mesh = getTestFunction().getFiniteElementSpace().getMesh();
      m_vector.reset(
            new VectorType(
               getAssembly().execute({mesh, fes, getIntegrators()})));
   }
}

#endif
