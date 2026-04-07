/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file RelativeError.h
 * @brief Relative error computation utilities for grid functions.
 *
 * This file provides utilities for computing relative errors between numerical
 * solutions (grid functions) and exact solutions in various norms. This is
 * essential for convergence analysis and solution verification in finite
 * element simulations.
 */
#ifndef RODIN_VARIATIONAL_RELATIVEERROR_H
#define RODIN_VARIATIONAL_RELATIVEERROR_H

#include "Rodin/Math/Common.h"

#include "Abs.h"
#include "GridFunction.h"
#include "Integral.h"
#include "RealFunction.h"
#include "Rodin/Assembly.h"

namespace Rodin::Variational
{
    /**
     * @ingroup RodinVariational
     * @brief Utility class for computing relative errors.
     *
     * RelativeError provides static methods to compute the relative error
     * between a computed grid function and an exact solution in various norms.
     * The relative error is defined as:
     * @f[
     *   E_{\text{rel}} = \frac{\|u_h - u_{\text{exact}}\|}{\|u_{\text{exact}}\|}
     * @f]
     * where @f$ u_h @f$ is the numerical solution and @f$ u_{\text{exact}} @f$
     * is the exact solution.
     *
     * ## Supported Norms
     * - **L1 norm**: @f$ \|u\|_{L^1} = \int_\Omega |u| \, dx @f$
     * - **L2 norm**: @f$ \|u\|_{L^2} = \sqrt{\int_\Omega u^2 \, dx} @f$
     * - **L∞ norm**: @f$ \|u\|_{L^\infty} = \max_{x \in \Omega} |u(x)| @f$
     *
     * ## Usage Example
     * ```cpp
     * // Compute L2 relative error
     * auto exact = [](const Point& p) { return sin(p.x()) * cos(p.y()); };
     * Real error_l2 = RelativeError::l2(uh, exact);
     * std::cout << "L2 relative error: " << error_l2 << std::endl;
     * 
     * // Compute L∞ relative error
     * Real error_linf = RelativeError::lInf(uh, exact);
     * ```
     *
     * ## Convergence Analysis
     * These error measures are used to verify convergence rates:
     * - P1 elements: @f$ E_{L^2} \sim O(h^2) @f$
     * - P2 elements: @f$ E_{L^2} \sim O(h^3) @f$
     */
    class RelativeError
    {
      public:
        /**
         * @brief Enumeration of supported norms.
         */
        enum class Norm
        {
          L1,    ///< L1 norm (integral of absolute value)
          L2,    ///< L2 norm (square root of integral of square)
          LInf   ///< L∞ norm (maximum absolute value)
        };

        /**
         * @brief Computes relative error in L1 norm.
         * @param[in] model Computed grid function
         * @param[in] exact Exact solution (function or grid function)
         * @returns Relative error @f$ \|u_h - u\|_{L^1} / \|u\|_{L^1} @f$
         */
        template <class FES, class Data, class FunctionType>
        static Real l1(const GridFunction<FES, Data>& model, const FunctionType& exact)
        {
          return compute(model, exact, Norm::L1);
        }

        /**
         * @brief Computes relative error in L2 norm.
         * @param[in] model Computed grid function
         * @param[in] exact Exact solution (function or grid function)
         * @returns Relative error @f$ \|u_h - u\|_{L^2} / \|u\|_{L^2} @f$
         */
        template <class FES, class Data, class FunctionType>
        static Real l2(const GridFunction<FES, Data>& model, const FunctionType& exact)
        {
          return compute(model, exact, Norm::L2);
        }

        /**
         * @brief Computes relative error in L∞ norm.
         * @param[in] model Computed grid function
         * @param[in] exact Exact solution (function or grid function)
         * @returns Relative error @f$ \|u_h - u\|_{L^\infty} / \|u\|_{L^\infty} @f$
         */
        template <class FES, class Data, class FunctionType>
        static Real lInf(const GridFunction<FES, Data>& model, const FunctionType& exact)
        {
          return compute(model, exact, Norm::LInf);
        }

        /**
         * @brief Computes relative error in specified norm.
         * @param[in] model Computed grid function
         * @param[in] exact Exact solution (function or grid function)
         * @param[in] norm Norm type to use
         * @returns Relative error in the specified norm
         */
        template <class FES, class Data, class FunctionType>
        static Real compute(const GridFunction<FES, Data>& model, const FunctionType& exact, const Norm& norm)
        {
            const auto& fes = model.getFiniteElementSpace();
            switch (norm)
            {
                case Norm::L1:
                {
                  GridFunction exactNorm(fes);
                  exactNorm = Abs(exact);
                  GridFunction diff(fes);
                  diff = Abs(model - exact);
                  return Integral(diff).compute() / Integral(exactNorm).compute();
                }
                case Norm::L2:
                {
                  GridFunction exactNorm(fes);
                  exactNorm = [&](const Geometry::Point& p) { auto v = exact(p); return Math::dot(v, v); };
                  GridFunction diff(fes);
                  diff = [&](const Geometry::Point& p) { auto v = exact(p) - model(p); return Math::dot(v, v); };
                  return Math::sqrt(Integral(diff).compute()) / Math::sqrt(Integral(exactNorm).compute());
                }
                case Norm::LInf:
                {
                  GridFunction exactNorm(fes);
                  exactNorm = Abs(exact);
                  GridFunction diff(fes);
                  diff = Abs(model - exact);
                  return diff.max() / exactNorm.max();
                }
            }
            return Math::nan<Real>();
        }
    };
}

#endif

