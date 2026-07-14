/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file NewtonRaphson.h
 * @brief Newton-Raphson method for root finding.
 *
 * This file implements a robust Newton-Raphson algorithm with automatic
 * bracketing and fallback mechanisms for finding roots of nonlinear equations.
 */
#ifndef RODIN_MATH_ROOTFINDING_NEWTONRAPHSON_H
#define RODIN_MATH_ROOTFINDING_NEWTONRAPHSON_H

#include <cmath>
#include <optional>

namespace Rodin::Math::RootFinding
{
  /**
   * @brief Newton-Raphson method for finding roots of scalar equations.
   *
   * This class implements a robust variant of the Newton-Raphson method for
   * finding roots of equations @f$ f(x) = 0 @f$. The method combines Newton's
   * method with interval bracketing and automatic fallbacks to ensure convergence.
   *
   * ## Algorithm
   * The Newton-Raphson iteration is:
   * @f[
   *   x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}
   * @f]
   *
   * ## Robustness Features
   * - **Interval bracketing**: Maintains a bracket @f$ [a, b] @f$ where @f$ f(a)f(b) \leq 0 @f$
   * - **Safeguards**: Falls back to bisection if Newton step goes outside bracket
   * - **Secant fallback**: Uses secant method if derivative is too small
   * - **Convergence checks**: Multiple criteria for termination
   *
   * ## Convergence Criteria
   * The algorithm terminates when any of these conditions are met:
   * 1. **Function tolerance**: @f$ |f(x)| < \varepsilon_{\text{abs}} @f$
   * 2. **Interval tolerance**: @f$ |b - a| < \varepsilon_{\text{abs}} + \varepsilon_{\text{rel}} |x| @f$
   * 3. **Maximum iterations**: Iteration count exceeds maximum
   *
   * ## Example Usage
   * ```cpp
   * // Find root of f(x) = x^2 - 2 in [1, 2]
   * NewtonRaphson<Real> solver(1e-12, 1e-9, 1e-12, 25);
   * auto f = [](Real x) {
   *   return std::make_pair(x*x - 2.0, 2.0*x);  // (f, f')
   * };
   * auto root = solver.solve(f, 1.5, 1.0, 2.0);
   * if (root) {
   *   std::cout << "Root: " << *root << std::endl;  // sqrt(2)
   * }
   * ```
   *
   * @tparam Scalar The scalar type (typically Real or double)
   */
  template <class Scalar>
  class NewtonRaphson
  {
    public:
      /**
       * @brief Constructs a Newton-Raphson solver with specified tolerances.
       *
       * @param[in] absTTol Absolute tolerance for interval width (default: @f$ 10^{-12} @f$)
       * @param[in] relTTol Relative tolerance for interval width (default: @f$ 10^{-9} @f$)
       * @param[in] absGTol Absolute tolerance for function value (default: @f$ 10^{-12} @f$)
       * @param[in] maxIter Maximum number of iterations (default: 25)
       */
      NewtonRaphson(Scalar absTTol = 1e-12, Scalar relTTol = 1e-9, Scalar absGTol = 1e-12,
        std::size_t maxIter = 25)
        : m_absTTol(absTTol),
          m_relTTol(relTTol),
          m_absGTol(absGTol),
          m_maxIter(maxIter)
      {}

      /**
       * @brief Solves for a root in the given interval.
       *
       * Finds a value @f$ x^* \in [a, b] @f$ such that @f$ f(x^*) \approx 0 @f$.
       * The function @f$ f @f$ must satisfy @f$ f(a) f(b) \leq 0 @f$ (root bracketing).
       *
       * @tparam F Function type that takes Scalar and returns std::pair<Scalar, Scalar>
       * @param[in] f Function object returning @f$ (f(x), f'(x)) @f$ as a pair
       * @param[in] t0 Initial guess for the root (should be in @f$ (a, b) @f$)
       * @param[in] lo Lower bound @f$ a @f$ of search interval
       * @param[in] hi Upper bound @f$ b @f$ of search interval
       * @return Optional containing the root if found, empty optional if failed
       *
       * @pre @f$ a < b @f$ and @f$ f(a) f(b) \leq 0 @f$
       * @pre Both @f$ a @f$ and @f$ b @f$ must be finite
       *
       * @note If derivative @f$ f'(x) @f$ is unavailable or unreliable, return NaN
       * for the derivative component and the method will use secant approximation.
       */
      template <class F>
      std::optional<Scalar> solve(const F& f, Scalar t0, Scalar lo, Scalar hi) const
      {
        // basic interval sanity
        assert(std::isfinite(lo) && std::isfinite(hi));
        assert(lo < hi);

        Scalar a = lo, b = hi;

        // Initial values and bracket check
        auto fa = f(a);
        Scalar fav = fa.first;

        auto fb = f(b);
        Scalar fbv = fb.first;

        // endpoints must be finite
        assert(std::isfinite(fav) && std::isfinite(fbv));

        // either already bracketed or caller misused the API
        assert(fav * fbv <= 0); // Root not bracketed at [lo, hi]

        if (!(fav * fbv <= 0))
          return {};

        if (!(t0 > a && t0 < b))
          t0 = 0.5 * (a + b);

        Scalar t = t0, tPrev = t;
        Scalar gPrev = std::numeric_limits<Scalar>::quiet_NaN();
        bool havePrev = false;

        for (std::size_t it = 0; it < m_maxIter; ++it)
        {
          const auto ft = f(t);
          const Scalar fv = ft.first;
          const Scalar fg = ft.second;

          // function and derivative must be finite
          assert(std::isfinite(fv));
          assert(std::isfinite(fg) || std::isnan(fg)); // allow NaN to trigger safeguards

          if (std::fabs(fv) < m_absGTol)
            return t;

          const Scalar at = m_absTTol + m_relTTol * std::fabs(t);
          if ((b - a) < at)
            return 0.5 * (a + b);

          // Maintain bracket using *current* endpoint values
          if (fv * fav < 0)
          {
            b = t;
            fbv = fv;
          }
          else
          {
            a = t;
            fav = fv;
          }

          if (std::fabs(fav) < m_absGTol)
            return a;

          if (std::fabs(fbv) < m_absGTol)
            return b;

          if (fav * fbv > 0)
            return {};

          // Newton step with scale-aware guard on gp
          Scalar dt;
          const Scalar gpTol = std::numeric_limits<Scalar>::epsilon() *
            (static_cast<Scalar>(1) + std::fabs(fg));
          if (std::fabs(fg) > gpTol)
          {
            dt = -fv / fg;
          }
          else if (havePrev)
          {
            const Scalar denom = (fv - gPrev);
            if (std::fabs(denom) > gpTol)
              dt = -fv * (t - tPrev) / denom; // secant
            else
              dt = 0;
          }
          else
          {
            dt = 0;
          }

          Scalar tNew = t + dt;

          // Safeguard: keep inside (a,b). If outside or too big, bisect.
          const Scalar mid = 0.5 * (a + b);
          const Scalar half = 0.5 * (b - a);

          if (tNew <= a || tNew >= b || std::fabs(tNew - t) > half)
            tNew = mid;

          if (std::fabs(tNew - t) < at)
            return tNew;

          tPrev = t;
          gPrev = fv;
          havePrev = true;
          t = tNew;
        }
        return {};
      }

    private:
      Scalar m_absTTol; ///< Absolute tolerance for interval width
      Scalar m_relTTol; ///< Relative tolerance for interval width
      Scalar m_absGTol; ///< Absolute tolerance for function value
      std::size_t m_maxIter; ///< Maximum number of iterations
  };
}

#endif
