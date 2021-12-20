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

#include "DirichletBC.h"

namespace Rodin::Variational
{
   DirichletBC::DirichletBC(int bdrAttr, const ScalarCoefficientBase& value)
      :  m_bdrAttr(bdrAttr),
         m_value(std::unique_ptr<ScalarCoefficientBase>(value.copy()))
   {}

   DirichletBC::DirichletBC(int bdrAttr, const VectorCoefficientBase& value)
      :  m_bdrAttr(bdrAttr),
         m_value(std::unique_ptr<VectorCoefficientBase>(value.copy()))
   {}

   DirichletBC::DirichletBC(const DirichletBC& other)
      :  m_bdrAttr(other.m_bdrAttr)
   {
      std::visit(
         [this](auto&& arg)
         {
            using T = std::decay_t<decltype(arg)>;
            if constexpr (std::is_same_v<T, std::unique_ptr<VectorCoefficientBase>>)
            {
               m_value = std::unique_ptr<VectorCoefficientBase>(arg->copy());
            }
            else
            {
               m_value = std::unique_ptr<ScalarCoefficientBase>(arg->copy());
            }
         }, other.m_value);
   }

   int DirichletBC::getBoundaryAttribute() const
   {
      return m_bdrAttr;
   }

   void DirichletBC::imposeOn(ProblemBase& pb) const
   {
      int maxBdrAttr = pb.getSolution()
                         .getHandle()
                         .FESpace()
                         ->GetMesh()
                         ->bdr_attributes.Max();

      if (m_bdrAttr > maxBdrAttr)
         Rodin::Alert::Exception(
               "DirichletBC boundary attribute is out of range.").raise();

      // Project the coefficient onto the boundary
      mfem::Array<int> essBdr(maxBdrAttr);
      essBdr = 0;
      essBdr[m_bdrAttr - 1] = 1;

      std::visit(
         [&pb, &essBdr](auto&& arg)
         {
            using T = std::decay_t<decltype(arg)>;
            if constexpr (std::is_same_v<T, std::unique_ptr<VectorCoefficientBase>>)
            {
               arg->buildMFEMVectorCoefficient();
               pb.getSolution()
                 .getHandle()
                 .ProjectBdrCoefficient(
                       arg->getMFEMVectorCoefficient(), essBdr);
            }
            else
            {
               arg->buildMFEMCoefficient();
               pb.getSolution()
                 .getHandle()
                 .ProjectBdrCoefficient(
                       arg->getMFEMCoefficient(), essBdr);
            }
         }, m_value);

      // Keep track of the boundary attributes that have been projected
      pb.getEssentialBoundary()[m_bdrAttr - 1] = 1;
   }
}
