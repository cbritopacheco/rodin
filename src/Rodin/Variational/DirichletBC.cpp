/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Alert.h"

#include "Problem.h"

#include "DirichletBC.h"

namespace Rodin::Variational
{
   void DirichletBC<ScalarCoefficientBase>::imposeOn(ProblemBase& pb)
   {
      int maxBdrAttr = pb.getSolution()
                         .getFiniteElementSpace()
                         .getMesh()
                         .getHandle().bdr_attributes.Max();

      if (getBoundaryAttribute() > maxBdrAttr)
         Rodin::Alert::Exception(
               "DirichletBC boundary attribute is out of range.").raise();

      // Project the coefficient onto the boundary
      m_essBdr = mfem::Array<int>(maxBdrAttr);
      m_essBdr = 0;
      m_essBdr[getBoundaryAttribute() - 1] = 1;

      getValue().buildMFEMCoefficient();
      pb.getSolution()
        .getHandle()
        .ProjectBdrCoefficient(
              getValue().getMFEMCoefficient(), m_essBdr);

      // Keep track of the boundary attributes that have been projected
      pb.getEssentialBoundary()[getBoundaryAttribute() - 1] = 1;
   }

   void DirichletBC<VectorCoefficientBase>::imposeOn(ProblemBase& pb)
   {
      int maxBdrAttr = pb.getSolution()
                         .getFiniteElementSpace()
                         .getMesh()
                         .getHandle().bdr_attributes.Max();

      if (getBoundaryAttribute() > maxBdrAttr)
         Rodin::Alert::Exception(
               "DirichletBC boundary attribute is out of range.").raise();

      // Project the coefficient onto the boundary
      m_essBdr = mfem::Array<int>(maxBdrAttr);
      m_essBdr = 0;
      m_essBdr[getBoundaryAttribute() - 1] = 1;

      getValue().buildMFEMVectorCoefficient();
      pb.getSolution()
        .getHandle()
        .ProjectBdrCoefficient(
              getValue().getMFEMVectorCoefficient(), m_essBdr);

      // Keep track of the boundary attributes that have been projected
      pb.getEssentialBoundary()[getBoundaryAttribute() - 1] = 1;
   }
}
