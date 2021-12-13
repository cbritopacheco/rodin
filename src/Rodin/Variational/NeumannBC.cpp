/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Alert.h"

#include "ScalarCoefficient.h"
#include "VectorCoefficient.h"

#include "NeumannBC.h"

namespace Rodin::Variational
{
   NeumannBC::NeumannBC(int bdrAttr, const VectorCoefficient& value)
      :  m_bdrAttr(bdrAttr),
         m_value(std::unique_ptr<VectorCoefficient>(value.copy()))
   {}

   NeumannBC::NeumannBC(const NeumannBC& other)
      :  m_bdrAttr(other.m_bdrAttr),
         m_problem(other.m_problem)
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

   NeumannBC& NeumannBC::setProblem(ProblemBase& problem)
   {
      m_problem.emplace(problem);
      return *this;
   }

   void NeumannBC::eval()
   {
      assert(m_problem);

      int maxBdrAttr = m_problem->get()
                                 .getSolution()
                                 .getHandle()
                                 .FESpace()
                                 ->GetMesh()
                                 ->bdr_attributes.Max();

      if (m_bdrAttr > maxBdrAttr)
         Rodin::Alert::Exception(
               "NeumannBC boundary attribute is out of range.").raise();

      // Project the coefficient onto the boundary
      m_nbcBdr = mfem::Array<int>(maxBdrAttr);
      m_nbcBdr = 0;
      m_nbcBdr[m_bdrAttr - 1] = 1;

      std::visit(
         [this](auto&& arg)
         {
            using T = std::decay_t<decltype(arg)>;
            if constexpr (std::is_same_v<T, std::unique_ptr<VectorCoefficientBase>>)
            {
               arg->buildMFEMVectorCoefficient();
               m_problem->get()
                        .getLinearForm()
                        .getHandle()
                        .AddBoundaryIntegrator(
                              new mfem::VectorBoundaryLFIntegrator(
                                 arg->getMFEMVectorCoefficient()), m_nbcBdr);
            }
            else
            {
               arg->buildMFEMCoefficient();
               m_problem->get()
                        .getLinearForm()
                        .getHandle()
                        .AddBoundaryIntegrator(
                              new mfem::BoundaryLFIntegrator(
                                 arg->getMFEMCoefficient()), m_nbcBdr);
            }
         }, m_value);
   }

   NeumannBC* NeumannBC::copy() const noexcept
   {
      return new NeumannBC(*this);
   }
}
