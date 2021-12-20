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
   NeumannBC::NeumannBC(int bdrAttr, const ScalarCoefficientBase& value)
      :  m_bdrAttr(bdrAttr),
         m_value(std::unique_ptr<ScalarCoefficientBase>(value.copy()))
   {}

   NeumannBC::NeumannBC(int bdrAttr, const VectorCoefficientBase& value)
      :  m_bdrAttr(bdrAttr),
         m_value(std::unique_ptr<VectorCoefficientBase>(value.copy()))
   {}

   NeumannBC::NeumannBC(const NeumannBC& other)
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

   int NeumannBC::getBoundaryAttribute() const
   {
      return m_bdrAttr;
   }

   void NeumannBC::imposeOn(ProblemBase& pb) const
   {
      int maxBdrAttr = pb.getSolution()
                         .getHandle()
                         .FESpace()
                         ->GetMesh()
                         ->bdr_attributes.Max();

      if (m_bdrAttr > maxBdrAttr)
         Rodin::Alert::Exception(
               "NeumannBC boundary attribute is out of range.").raise();

      // Project the coefficient onto the boundary
      mfem::Array<int> nbcBdr(maxBdrAttr);
      nbcBdr = 0;
      nbcBdr[m_bdrAttr - 1] = 1;

      std::visit(
         [&pb, &nbcBdr](auto&& arg)
         {
            using T = std::decay_t<decltype(arg)>;
            if constexpr (std::is_same_v<T, std::unique_ptr<VectorCoefficientBase>>)
            {
               arg->buildMFEMVectorCoefficient();
               pb.getLinearForm()
                  .getHandle()
                  .AddBoundaryIntegrator(
                        new mfem::VectorBoundaryLFIntegrator(
                           arg->getMFEMVectorCoefficient()), nbcBdr);
            }
            else
            {
               arg->buildMFEMCoefficient();
               pb.getLinearForm()
                 .getHandle()
                 .AddBoundaryIntegrator(
                       new mfem::BoundaryLFIntegrator(
                          arg->getMFEMCoefficient()), nbcBdr);
            }
         }, m_value);
      pb.getLinearForm().getHandle().Assemble();
   }
}
