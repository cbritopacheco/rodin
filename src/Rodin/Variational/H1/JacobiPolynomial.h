/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_JACOBIPOLYNOMIAL_H
#define RODIN_VARIATIONAL_H1_JACOBIPOLYNOMIAL_H

#include <cstddef>

#include "Rodin/Types.h"

namespace Rodin::Variational
{
  /**
   * @brief Jacobi polynomial evaluator.
   *
   * Evaluates the Jacobi polynomial @f$ P_K^{(\alpha, \beta)}(x) @f$ and its
   * derivative @f$ \partial_x P_K^{(\alpha, \beta)}(x) @f$ using:
   *
   *  - the three-term recurrence for @f$ P_n^{(\alpha,\beta)} @f$
   *  - the identity
   *    @f[
   *      \frac{d}{dx} P_K^{(\alpha, \beta)}(x)
   *      =
   *      \frac{1}{2}\bigl(K + \alpha + \beta + 1\bigr)
   *      P_{K - 1}^{(\alpha + 1, \beta + 1)}(x).
   *    @f]
   *
   * @tparam K  Polynomial degree.
   */
  template <size_t K>
  class JacobiPolynomial
  {
    public:
      /**
       * @brief Computes @f$ P_K^{(\Alpha, \Beta)}(x) @f$ and
       * @f$ \partial_x P_K^{(\Alpha, \Beta)}(x) @f$.
       *
       * @param[out] P   Value of @f$ P_K^{(\Alpha, \Beta)}(x) @f$.
       * @param[out] dP  Value of @f$ \partial_x P_K^{(\Alpha, \Beta)}(x) @f$.
       * @param[in]  x   Evaluation point.
       */
      static constexpr void getValue(Real& P, Real& dP, Real alpha, Real beta, Real x)
      {
        P = getValue(K, alpha, beta, x);

        if constexpr (K == 0)
        {
          dP = static_cast<Real>(0.0);
        }
        else
        {
          const Real Pm1 = getValue(
              K - 1,
              alpha + static_cast<Real>(1.0),
              beta  + static_cast<Real>(1.0),
              x);

          dP = static_cast<Real>(0.5)
              * (static_cast<Real>(K) + alpha + beta + static_cast<Real>(1.0))
              * Pm1;
        }
      }

    private:
      static constexpr Real getValue(size_t k, Real alpha, Real beta, Real x)
      {
        if (k == 0)
          return static_cast<Real>(1.0);

        const Real a_ab  = alpha + beta;
        const Real two   = static_cast<Real>(2.0);

        Real P0 = static_cast<Real>(1.0);
        Real P1 = static_cast<Real>(0.5) *
                  (alpha - beta + (a_ab + two) * x);

        if (k == 1)
          return P1;

        Real P2 = static_cast<Real>(0.0);

        for (size_t n = 1; n < k; ++n)
        {
          const Real nn    = static_cast<Real>(n);
          const Real two_n = two * nn;

          const Real a1 =
            two * (nn + static_cast<Real>(1.0))
                * (nn + a_ab + static_cast<Real>(1.0))
                * (two_n + a_ab);

          const Real a2 =
            (two_n + a_ab + static_cast<Real>(1.0))
            * (alpha * alpha - beta * beta);

          const Real a3 =
            (two_n + a_ab)
            * (two_n + a_ab + static_cast<Real>(1.0))
            * (two_n + a_ab + static_cast<Real>(2.0));

          const Real a4 =
            two * (nn + alpha) * (nn + beta)
            * (two_n + a_ab + static_cast<Real>(2.0));

          P2 = ((a2 + a3 * x) * P1 - a4 * P0) / a1;
          P0 = P1;
          P1 = P2;
        }

        return P1;
      }
  };
}

#endif
