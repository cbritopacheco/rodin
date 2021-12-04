/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_PCG_H
#define RODIN_SOLVER_PCG_H

#include <optional>
#include <functional>

#include <mfem.hpp>

#include "Rodin/Variational/Problem.h"

#include "Solver.h"

namespace Rodin::Solver
{
   class PCG : public Solver
   {
      public:
         PCG() = default;
         ~PCG() = default;

         PCG& printIterations(bool printIterations);
         PCG& setMaxIterations(int maxIterations);
         PCG& setRelativeTolerance(double rtol);
         PCG& setAbsoluteTolerance(double atol);
         PCG& setSmoother(mfem::Solver& smoother);

         void solve(Variational::ProblemBase& problem) const override;

      private:
         bool m_printIterations;
         int  m_maxIterations;
         double m_rtol, m_atol;
         std::optional<std::reference_wrapper<mfem::Solver>> m_smoother;
   };
}

#endif
