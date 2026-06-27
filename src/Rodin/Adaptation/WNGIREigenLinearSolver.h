/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_WNGIREIGENLINEARSOLVER_H
#define RODIN_ADAPTATION_WNGIREIGENLINEARSOLVER_H

#include "WNGIRParameters.h"

namespace Rodin::Adaptation
{
  template <class ProblemT>
  class WNGIREigenLinearSolver
  {
    public:
      bool solve(
          ProblemT& problem,
          Math::Vector<Real>& solution,
          std::size_t& iterations,
          Real& error,
          const WNGIRParameters& parameters)
      {
        const auto& system = problem.getLinearSystem();
        using CGSolver = Eigen::ConjugateGradient<
          Eigen::SparseMatrix<Real>, Eigen::Lower | Eigen::Upper,
          Eigen::DiagonalPreconditioner<Real>>;
        CGSolver cg;
        const Eigen::Index ndofs = system.getOperator().rows();
        const std::size_t maxIterations = parameters.cgMaxIterations > 0
          ? parameters.cgMaxIterations
          : std::min<std::size_t>(
              2000,
              std::max<std::size_t>(100,
                2 * static_cast<std::size_t>(
                  std::max<Eigen::Index>(ndofs, 1))));
        cg.setMaxIterations(static_cast<int>(maxIterations));
        cg.setTolerance(parameters.cgRelativeTolerance > Real(0)
            ? parameters.cgRelativeTolerance : Real(1e-6));
        cg.compute(system.getOperator());
        if (cg.info() != Eigen::Success)
          return false;
        solution = cg.solve(system.getVector());
        iterations = static_cast<std::size_t>(
            std::max<Eigen::Index>(cg.iterations(), 0));
        error = cg.error();
        return cg.info() == Eigen::Success && solution.allFinite();
      }
  };
}

#endif
