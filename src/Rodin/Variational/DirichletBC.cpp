/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Alert.h"

#include "ScalarCoefficient.h"
#include "VectorCoefficient.h"

#include "DirichletBC.h"

namespace Rodin::Variational
{
   DirichletBC::DirichletBC(int bdrAttr, const ScalarCoefficientBase& value)
      :  m_bdrAttr(bdrAttr),
         m_value(std::unique_ptr<ScalarCoefficientBase>(value.copy()))
   {}

   DirichletBC::DirichletBC(int bdrAttr, const VectorCoefficient& value)
      :  m_bdrAttr(bdrAttr),
         m_value(std::unique_ptr<VectorCoefficient>(value.copy()))
   {}

   DirichletBC::DirichletBC(const DirichletBC& other)
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

   DirichletBC& DirichletBC::setProblem(ProblemBase& problem)
   {
      m_problem.emplace(problem);
      return *this;
   }

   void DirichletBC::eval()
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
               "DirichletBC boundary attribute is out of range.").raise();

      // Project the coefficient onto the boundary
      m_essBdr = mfem::Array<int>(maxBdrAttr);
      m_essBdr = 0;
      m_essBdr[m_bdrAttr - 1] = 1;

      std::visit(
         [this](auto&& arg)
         {
            using T = std::decay_t<decltype(arg)>;
            if constexpr (std::is_same_v<T, std::unique_ptr<VectorCoefficientBase>>)
            {
               arg->buildMFEMVectorCoefficient();
               m_problem->get().getSolution()
                               .getHandle()
                               .ProjectBdrCoefficient(
                                     arg->getMFEMVectorCoefficient(), m_essBdr);
            }
            else
            {
               arg->buildMFEMCoefficient();
               m_problem->get().getSolution()
                               .getHandle()
                               .ProjectBdrCoefficient(
                                     arg->getMFEMCoefficient(), m_essBdr);
            }
         }, m_value);

      // Keep track of the boundary attributes that have been projected
      m_problem->get().getEssentialBoundary()[m_bdrAttr - 1] = 1;
   }

   DirichletBC* DirichletBC::copy() const noexcept
   {
      return new DirichletBC(*this);
   }
}
