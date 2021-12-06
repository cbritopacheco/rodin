/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_PROBLEM_HPP
#define RODIN_VARIATIONAL_PROBLEM_HPP

#include "GridFunction.h"

#include "Problem.h"

namespace Rodin::Variational
{
   template <class FEC>
   Problem<FEC>::Problem(GridFunction<FEC>& u, GridFunction<FEC>& v)
      :  m_solution(u),
         m_adjoint(v),
         m_bilinearForm(u.getFiniteElementSpace()),
         m_linearForm(u.getFiniteElementSpace()),
         m_essBdr(u.getFiniteElementSpace()
                   .getMesh().getHandle().bdr_attributes.Max())
   {
      // assert(u.getFiniteElementSpace() == v.getFiniteElementSpace());
      m_essBdr = 0;
   }

   template <class FEC>
   template <class Derived>
   Problem<FEC>& Problem<FEC>::operator=(
         const FormLanguage::ProblemBody<Derived>& rhs)
   {
      m_ast = std::unique_ptr<FormLanguage::ProblemBodyBase>(rhs.copy());

      // We first need to evaluate the problem body to build the problem
      m_ast->setProblem(*this).eval();

      // Check everything is fine
      if (m_bilinearForm.getHandle().GetDBFI()->Size() == 0)
         Alert::Exception("The number of bilinear form integrators is zero.").raise();

      return *this;
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
