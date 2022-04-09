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
      problem.update().assemble();

      mfem::GSSmoother smoother;
      smoother.SetOperator(problem.getStiffnessMatrix());

      mfem::PCG(
            problem.getStiffnessMatrix(),
            smoother,
            problem.getMassVector(),
            problem.getInitialGuess(),
            m_printIterations, m_maxIterations, m_rtol, m_atol);

      problem.recoverSolution();
   }
}

