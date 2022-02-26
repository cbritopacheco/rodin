/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "CG.h"

namespace Rodin::Solver
{
   CG& CG::printIterations(bool printIterations)
   {
      m_printIterations = printIterations;
      return *this;
   }

   CG& CG::setMaxIterations(int maxIterations)
   {
      m_maxIterations = maxIterations;
      return *this;
   }

   CG& CG::setRelativeTolerance(double rtol)
   {
      m_rtol = rtol;
      return *this;
   }

   CG& CG::setAbsoluteTolerance(double atol)
   {
      m_atol = atol;
      return *this;
   }

   void CG::solve(Variational::ProblemBase& problem)
   {
      auto& a = problem.getBilinearForm();
      auto& b = problem.getLinearForm();
      auto& u = problem.getSolution();

      problem.update();
      problem.assemble();

      // Compute essential true degrees of freedom
      int maxBdrAttr = u.getFiniteElementSpace().getMesh().getHandle().bdr_attributes.Max();
      mfem::Array<int> essTrueDofList;
      u.getFiniteElementSpace()
       .getFES()
       .GetEssentialTrueDofs(
             Utility::set2marker(problem.getEssentialBoundary(), maxBdrAttr), essTrueDofList);

      // Form the linear system
      mfem::SparseMatrix A;
      mfem::Vector B, X;
      a.getHandle()
       .FormLinearSystem(essTrueDofList, u.getHandle(), b.getHandle(), A, X, B);

      // Solve
      mfem::GSSmoother smoother(A);
      mfem::PCG(A, smoother, B, X, m_printIterations, m_maxIterations, m_rtol, m_atol);

      a.getHandle()
       .RecoverFEMSolution(X, b.getHandle(), u.getHandle());
   }
}

