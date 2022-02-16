/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Alert.h"

#include "Problem.h"
#include "ScalarCoefficient.h"
#include "VectorCoefficient.h"

#include "NeumannBC.h"

namespace Rodin::Variational
{
   void NeumannBC<ScalarCoefficientBase>::imposeOn(ProblemBase& pb)
   {
      int maxBdrAttr = pb.getSolution()
                         .getFiniteElementSpace()
                         .getMesh()
                         .getHandle().bdr_attributes.Max();

      if (getBoundaryAttribute() > maxBdrAttr)
         Alert::Exception(
               "NeumannBC boundary attribute is out of range.").raise();

      // Project the coefficient onto the boundary
      m_nbcBdr = mfem::Array<int>(maxBdrAttr);
      m_nbcBdr = 0;
      m_nbcBdr[getBoundaryAttribute() - 1] = 1;

      getValue().build();
      pb.getLinearForm()
        .getHandle()
        .AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(getValue().get()), m_nbcBdr);
   }

   void NeumannBC<VectorCoefficientBase>::imposeOn(ProblemBase& pb)
   {
      int maxBdrAttr = pb.getSolution()
                         .getFiniteElementSpace()
                         .getMesh()
                         .getHandle().bdr_attributes.Max();

      if (getBoundaryAttribute() > maxBdrAttr)
         Alert::Exception(
               "NeumannBC boundary attribute is out of range.").raise();

      // Project the coefficient onto the boundary
      m_nbcBdr = mfem::Array<int>(maxBdrAttr);
      m_nbcBdr = 0;
      m_nbcBdr[getBoundaryAttribute() - 1] = 1;

      getValue().build();
      pb.getLinearForm()
        .getHandle()
        .AddBoundaryIntegrator(
              new mfem::VectorBoundaryLFIntegrator(
                 getValue().get()), m_nbcBdr);
   }
}
