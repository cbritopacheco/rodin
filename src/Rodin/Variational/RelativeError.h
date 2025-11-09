/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_RELATIVEERROR_H
#define RODIN_VARIATIONAL_RELATIVEERROR_H

/**
 * @file
 * @brief Relative error computation utilities for finite element solutions.
 *
 * This file provides utilities to compute the relative error between a computed
 * finite element solution and an exact (analytical) solution. The relative error
 * in a norm @f$ \|\cdot\| @f$ is defined as:
 * @f[
 *   \text{RelErr} = \frac{\|u_h - u_{exact}\|}{\|u_{exact}\|}
 * @f]
 *
 * ## Supported Norms
 * - **L1**: @f$ \|f\|_{L^1} = \int_\Omega |f| \, dx @f$
 * - **L2**: @f$ \|f\|_{L^2} = \sqrt{\int_\Omega |f|^2 \, dx} @f$
 * - **L∞**: @f$ \|f\|_{L^\infty} = \max_{x \in \Omega} |f(x)| @f$
 */

#include "GridFunction.h"

namespace Rodin::Variational
{
    class RelativeError
    {
      public:
        enum class Norm
        {
          L1,
          L2,
          LInf
        };

        template <class FES, class FunctionType>
        static Real l1(const GridFunction<FES>& model, const FunctionType& exact)
        {
          return compute(model, exact, Norm::L1);
        }

        template <class FES, class FunctionType>
        static Real l2(const GridFunction<FES>& model, const FunctionType& exact)
        {
          return compute(model, exact, Norm::L2);
        }

        template <class FES, class FunctionType>
        static Real lInf(const GridFunction<FES>& model, const FunctionType& exact)
        {
          return compute(model, exact, Norm::LInf);
        }

        template <class FES, class FunctionType>
        static Real compute(const GridFunction<FES>& model, const FunctionType& exact, const Norm& norm)
        {
            switch (norm)
            {
                case Norm::L1:
                {
                  const auto& fes = model.getFiniteElementSpace();
                  GridFunction<FES> exactNorm(fes);
                  exactNorm = Abs(exact);
                  exactNorm.setWeights();
                  GridFunction<FES> diff(fes);
                  diff = Abs(model - exact);
                  diff.setWeights();
                  return Integral(diff).compute() / Integral(exactNorm).compute();
                }
                case Norm::L2:
                {
                  const auto& fes = model.getFiniteElementSpace();
                  GridFunction<FES> exactNorm(fes);
                  exactNorm = [&](const Geometry::Point& p) { auto v = exact(p); return Math::dot(v, v); };
                  exactNorm.setWeights();
                  GridFunction<FES> diff(fes);
                  diff = [&](const Geometry::Point& p) { auto v = exact(p) - model(p); return Math::dot(v, v); };
                  diff.setWeights();
                  return sqrt(Integral(diff).compute()) / sqrt(Integral(exactNorm).compute());
                }
                case Norm::LInf:
                {
                  const auto& fes = model.getFiniteElementSpace();
                  GridFunction<FES> exactNorm(fes);
                  exactNorm = Abs(exact);
                  GridFunction<FES> diff(fes);
                  diff = Abs(model - exact);
                  return diff.max() / exactNorm.max();
                }
            }
            return Math::Constants::nan();
        }
    };
}

#endif

