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
#include "FormLanguage/LinearFormIntegratorUnaryMinus.h"

#include "BoundaryCondition.h"
#include "GridFunction.h"
#include "DirichletBC.h"
#include "NeumannBC.h"

#include "Problem.h"

namespace Rodin::Variational
{
   template <class TrialFEC, class TestFEC>
   Problem<TrialFEC, TestFEC>::Problem(TrialFunction<TrialFEC>& u, TestFunction<TestFEC>&)
      :  m_solution(u),
         m_bilinearForm(u.getFiniteElementSpace()),
         m_linearForm(u.getFiniteElementSpace()),
         m_essBdr(u.getFiniteElementSpace()
                   .getMesh().getHandle().bdr_attributes.Max())
   {
      m_essBdr = 0;
   }

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

      // Neumann boundary conditions are imposed instantly
      for (auto& nbc : m_pb->getNeumannBCList())
         nbc->imposeOn(*this);

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
      m_solution.getFiniteElementSpace().update();
      m_solution.update();
      m_linearForm.update();
      m_bilinearForm.update();
      for (auto& dbc : m_pb->getDirichletBCList())
         dbc->imposeOn(*this);
   }

   template <class TrialFEC, class TestFEC>
   TrialFunction<TrialFEC>& Problem<TrialFEC, TestFEC>::getSolution()
   {
      return m_solution;
   }

   template <class TrialFEC, class TestFEC>
   mfem::Array<int>& Problem<TrialFEC, TestFEC>::getEssentialBoundary()
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
