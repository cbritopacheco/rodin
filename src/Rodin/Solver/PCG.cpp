/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "PCG.h"

namespace Rodin::Solver
{
   PCG& PCG::printIterations(bool printIterations)
   {
      m_printIterations = printIterations;
      return *this;
   }

   PCG& PCG::setMaxIterations(int maxIterations)
   {
      m_maxIterations = maxIterations;
      return *this;
   }

   PCG& PCG::setRelativeTolerance(double rtol)
   {
      m_rtol = rtol;
      return *this;
   }

   PCG& PCG::setAbsoluteTolerance(double atol)
   {
      m_atol = atol;
      return *this;
   }

   void PCG::solve(Variational::ProblemBase& problem) const
   {
      problem.assemble();

      auto& a = problem.getBilinearForm();
      auto& b = problem.getLinearForm();
      auto& u = problem.getSolution();

      // Compute essential true degrees of freedom
      mfem::Array<int> essTrueDofList;
      u.getHandle()
       .FESpace()
       ->GetEssentialTrueDofs(
                   problem.getEssentialBoundary(), essTrueDofList);

      // Form the linear system
      mfem::SparseMatrix A;
      mfem::Vector B, X;
      a.getHandle()
       .FormLinearSystem(essTrueDofList, u.getHandle(), b.getHandle(), A, X, B);

      mfem::GSSmoother M(A);
      mfem::PCG(A, M, B, X, m_printIterations, m_maxIterations, m_rtol, m_atol);

      a.getHandle()
       .RecoverFEMSolution(X, b.getHandle(), u.getHandle());
   }
}
