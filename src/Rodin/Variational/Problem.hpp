/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_PROBLEM_HPP
#define RODIN_VARIATIONAL_PROBLEM_HPP

#include "Rodin/Utility.h"

#include "FormLanguage/ProblemBody.h"

#include "GridFunction.h"
#include "DirichletBC.h"

#include "Problem.h"

namespace Rodin::Variational
{
   template <class TrialFEC, class TestFEC>
   Problem<TrialFEC, TestFEC>::Problem(TrialFunction<TrialFEC>& u, TestFunction<TestFEC>&)
      :  m_solution(u),
         m_bilinearForm(u.getFiniteElementSpace()),
         m_linearForm(u.getFiniteElementSpace())
   {}

   template <class TrialFEC, class TestFEC>
   Problem<TrialFEC, TestFEC>& Problem<TrialFEC, TestFEC>::operator=(const FormLanguage::ProblemBody& rhs)
   {
      m_pb.reset(rhs.copy());

      for (auto& bfi : m_pb->getBilinearFormDomainIntegratorList())
         m_bilinearForm.add(*bfi);

      // The LinearFormIntegrator instances have already been moved to the LHS
      for (auto& lfi : m_pb->getLinearFormDomainIntegratorList())
         m_linearForm.add(*lfi);
      for (auto& lfi : m_pb->getLinearFormBoundaryIntegratorList())
         m_linearForm.add(*lfi);

      m_solution.emplaceGridFunction();

      return *this;
   }

   template <class TrialFEC, class TestFEC>
   void Problem<TrialFEC, TestFEC>::assemble()
   {
      m_linearForm.assemble();
      m_bilinearForm.assemble();
   }

   template <class TrialFEC, class TestFEC>
   void Problem<TrialFEC, TestFEC>::update()
   {
      getSolution().getFiniteElementSpace().update();
      getSolution().update();
      getLinearForm().update();
      getBilinearForm().update();
      for (const auto& dbc : m_pb->getDirichletBCList())
      {
         // Project the coefficient onto the boundary
         if (m_solution.getFiniteElementSpace().getVectorDimension() == 1)
         {
            getSolution().projectOnBoundary(
                  dbc.getValue<ScalarCoefficientBase>(), dbc.getBoundaryAttributes());
         }
         else
         {
            assert(m_solution.getFiniteElementSpace().getVectorDimension() > 1);
            getSolution().projectOnBoundary(
                  dbc.getValue<VectorCoefficientBase>(), dbc.getBoundaryAttributes());
         }

         // Keep track of essential boundary attributes
         getEssentialBoundary().insert(dbc.getBoundaryAttributes().begin(), dbc.getBoundaryAttributes().end());
      }
   }

   template <class TrialFEC, class TestFEC>
   GridFunction<TrialFEC>& Problem<TrialFEC, TestFEC>::getSolution()
   {
      return m_solution.getGridFunction();
   }

   template <class TrialFEC, class TestFEC>
   std::set<int>& Problem<TrialFEC, TestFEC>::getEssentialBoundary()
   {
      return m_essBdr;
   }

   template <class TrialFEC, class TestFEC>
   BilinearForm<TrialFEC>& Problem<TrialFEC, TestFEC>::getBilinearForm()
   {
      return m_bilinearForm;
   }

   template <class TrialFEC, class TestFEC>
   LinearForm<TrialFEC>& Problem<TrialFEC, TestFEC>::getLinearForm()
   {
      return m_linearForm;
   }
}

#endif
