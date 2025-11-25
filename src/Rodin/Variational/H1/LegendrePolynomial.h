/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_LEGENDREPOLYNOMIAL_H
#define RODIN_VARIATIONAL_H1_LEGENDREPOLYNOMIAL_H

#include "Rodin/Types.h"

namespace Rodin::Variational
{
  /**
   * @brief Compile-time degree Legendre polynomial evaluator.
   *
   * This class provides an efficient and numerically stable evaluation of the
   * Legendre polynomial @f$ P_K(x) @f$ and its derivative @f$ P'_K(x) @f$ on
   * the interval @f$ [-1, 1] @f$, using the standard three-term recurrence
   * relation and its differentiated form.
   *
   * The polynomial degree @f$ K @f$ is encoded as a template parameter, which
   * allows the compiler to optimize and possibly unroll the recurrence for
   * fixed degrees typically used in high-order finite element methods.
   *
   * This evaluator is intended for internal use in the construction of
   * modal bases and Gauss–Lobatto–Legendre nodes for H1-conforming elements.
   *
   * @tparam K Polynomial degree of the Legendre polynomial @f$ P_K @f$.
   */
  template <size_t K>
  class LegendrePolynomial
  {
    public:
      /**
       * @brief Computes @f$ P_K(x) @f$ and @f$ P'_K(x) @f$ using the
       * three-term recurrence relation.
       *
       * For @f$ K = 0 @f$ and @f$ K = 1 @f$, explicit closed forms are used:
       * @f$ P_0(x) = 1 @f$, @f$ P_1(x) = x @f$.
       * For @f$ K \ge 2 @f$, the standard recurrence
       * @f[
       *   (k + 1) P_{k + 1}(x)
       *   =
       *   (2k + 1)\,x\,P_k(x) - k\,P_{k - 1}(x)
       * @f]
       * and its derivative are used to advance from @f$ P_0, P_1 @f$ up to
       * @f$ P_K @f$.
       *
       * The evaluation is intended for arguments @f$ x \in [-1, 1] @f$.
       *
       * @param[out] P  Value of @f$ P_K(x) @f$.
       * @param[out] dP Value of @f$ P'_K(x) @f$.
       * @param x       Evaluation point (assumed in @f$ [-1, 1] @f$).
       */
      static constexpr void getValue(Real& P, Real& dP, Real x)
      {
        if constexpr (K == 0)
        {
          P  = 1.0;
          dP = 0.0;
        }
        else if constexpr (K == 1)
        {
          P  = x;
          dP = 1.0;
        }
        else
        {
          // Recurrence initialization:
          // P_0(x) = 1,   P_1(x) = x
          // P'_0(x) = 0,  P'_1(x) = 1
          Real P0  = 1.0;
          Real P1 = x;
          Real dP0 = 0.0;
          Real dP1 = 1.0;

          // 3-term recurrence:
          // (k + 1) P_{k + 1} = (2k + 1) x P_k - k P_{k-1}
          // differentiated:
          // (k + 1) P'_{k + 1} = (2k + 1)(P_k + x P'_k) - k P'_{k-1}
          for (size_t k = 1; k < K; ++k)
          {
            const Real kk = static_cast<Real>(k);
            const Real two_k_plus_one = 2.0 * kk + 1.0;
            const Real inv_kp1 = 1.0 / (kk + 1.0);
            const Real alpha = two_k_plus_one * inv_kp1; // (2k + 1) / (k + 1)
            const Real beta = kk * inv_kp1; // k / (k + 1)
            const Real P2 = alpha * x * P1 - beta * P0;
            const Real dP2 = alpha * (P1 + x * dP1) - beta * dP0;
            P0 = P1;
            P1 = P2;
            dP0 = dP1;
            dP1 = dP2;
          }

          P = P1;
          dP = dP1;
        }
      }
  };
}

#endif
