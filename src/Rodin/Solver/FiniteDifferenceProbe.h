/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_FINITEDIFFERENCEPROBE_H
#define RODIN_SOLVER_FINITEDIFFERENCEPROBE_H

#include <cmath>
#include <functional>

#include "Rodin/Types.h"
#include "Rodin/Math/LinearSystem.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"

namespace Rodin::Solver
{
  /**
   * @brief Finite-difference utilities for probing nonlinear residuals/tangents.
   */
  class FiniteDifferenceProbe
  {
    public:
      using Scalar = Real;
      using DenseMatrix = Math::Matrix<Scalar>;
      using DenseVector = Math::Vector<Scalar>;

      struct Options
      {
        Scalar epsilon = 1e-7;
      };

      /**
       * @brief Returns the default Options instance.
       */
      static Options defaultOptions()
      {
        return Options{};
      }

      /**
       * @brief Builds FD Jacobian of a residual function R(x).
       *
       * @tparam ResidualEval Callable of signature `void(const DenseVector&, DenseVector&)`.
       */
      template <class ResidualEval>
      static DenseMatrix jacobian(const DenseVector& x, ResidualEval&& residual, const Options& options = defaultOptions())
      {
        DenseVector R0(x.size());
        residual(x, R0);

        DenseMatrix J(R0.size(), x.size());
        for (Eigen::Index j = 0; j < x.size(); ++j)
        {
          const Scalar h = options.epsilon * std::max(Scalar(1), std::abs(x[j]));
          DenseVector xp = x;
          DenseVector xm = x;
          xp[j] += h;
          xm[j] -= h;

          DenseVector Rp(R0.size());
          DenseVector Rm(R0.size());
          residual(xp, Rp);
          residual(xm, Rm);
          J.col(j) = (Rp - Rm) / (2.0 * h);
        }
        return J;
      }

      /**
       * @brief Computes relative Frobenius error between two Jacobians.
       */
      static Scalar relativeError(const DenseMatrix& reference, const DenseMatrix& approx)
      {
        const Scalar denom = std::max(reference.norm(), Scalar(1e-16));
        return (reference - approx).norm() / denom;
      }

      /**
       * @brief Probe Jacobian for a tangential problem assembled in Newton form.
       *
       * Uses
       * - @p setState(x): updates nonlinear state observed by the problem,
       * - `pb.assemble()` and `pb.getLinearSystem()` to read tangent operator and RHS,
       * - residual defined as @f$R(x) = -b(x)@f$.
       */
      template <class Problem, class SetStateFn>
      static DenseMatrix tangentialJacobian(
          Problem& pb,
          const DenseVector& x,
          SetStateFn&& setState,
          const Options& options = defaultOptions())
      {
        auto residual = [&](const DenseVector& state, DenseVector& out)
        {
          setState(state);
          pb.assemble();
          out = -pb.getLinearSystem().getVector();
        };
        return jacobian(x, residual, options);
      }
    };
}

#endif
