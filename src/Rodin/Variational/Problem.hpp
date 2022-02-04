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
   template <class FEC>
   Problem<FEC>::Problem(GridFunction<FEC>& u)
      :  m_solution(u),
         m_bilinearForm(u.getFiniteElementSpace()),
         m_linearForm(u.getFiniteElementSpace()),
         m_essBdr(u.getFiniteElementSpace()
                   .getMesh().getHandle().bdr_attributes.Max())
   {
      m_essBdr = 0;
   }

   template <class FEC>
   Problem<FEC>& Problem<FEC>::operator=(const FormLanguage::ProblemBody& rhs)
   {
      m_pb.reset(rhs.copy());

      for (auto& bfi : m_pb->getBilinearFormDomainIntegratorList())
         m_bilinearForm.add(static_cast<BilinearFormDomainIntegrator&>(*bfi));

      // The LinearFormIntegrator instances have already been moved to the LHS
      for (auto& lfi : m_pb->getLinearFormDomainIntegratorList())
         m_linearForm.add(static_cast<LinearFormDomainIntegrator&>(*lfi));
      for (auto& lfi : m_pb->getLinearFormBoundaryIntegratorList())
         m_linearForm.add(static_cast<LinearFormBoundaryIntegrator&>(*lfi));

      // Neumann boundary conditions are imposed instantly
      for (auto& nbc : m_pb->getNeumannBCList())
         nbc->imposeOn(*this);

      return *this;
   }

   template <class FEC>
   void Problem<FEC>::assemble()
   {
      m_linearForm.assemble();
      m_bilinearForm.assemble();
   }

   template <class FEC>
   void Problem<FEC>::update()
   {
      m_solution.getFiniteElementSpace().update();
      m_solution.update();
      m_linearForm.update();
      m_bilinearForm.update();
      for (auto& dbc : m_pb->getDirichletBCList())
         dbc->imposeOn(*this);
   }

   template <class FEC>
   GridFunction<FEC>& Problem<FEC>::getSolution()
   {
      return m_solution;
   }

   template <class FEC>
   mfem::Array<int>& Problem<FEC>::getEssentialBoundary()
   {
      return m_essBdr;
   }

   template <class FEC>
   BilinearForm<FEC>& Problem<FEC>::getBilinearForm()
   {
      return m_bilinearForm;
   }

   template <class FEC>
   LinearForm<FEC>& Problem<FEC>::getLinearForm()
   {
      return m_linearForm;
   }
}

#endif
