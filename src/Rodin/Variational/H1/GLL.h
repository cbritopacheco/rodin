/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_GLL_H
#define RODIN_VARIATIONAL_H1_GLL_H

#include <array>

#include "Rodin/Math/Constants.h"
#include "Rodin/Math/Common.h"

#include "LegendrePolynomial.h"

#define RODIN_VARIATIONAL_H1_GLL_TOLERANCE 1e-14

namespace Rodin::Variational
{
  /**
   * @brief Compile-time Gauss–Lobatto–Legendre (GLL) nodes on [-1, 1].
   *
   * This class provides the @f$ K+1 @f$ GLL nodes associated with the
   * Legendre polynomial degree @f$ K @f$. The nodes consist of the endpoints
   * @f$ \{-1, 1\} @f$ and the interior zeros of @f$ (1 - x^2) P'_K(x) @f$.
   *
   * The nodes are computed once at compile time (subject to compiler support
   * for `constexpr` floating-point operations) using a Newton iteration with
   * Chebyshev initial guesses, and cached in a static `constexpr` array.
   *
   * @tparam K Polynomial degree.
   * @tparam Tolerance Newton stopping tolerance (compile-time constant).
   * @tparam MaxIt Maximum number of Newton iterations (compile-time constant).
   */
  template <size_t K, size_t MaxIt = 25>
  class GLL
  {
    public:
      /// Number of GLL nodes (K + 1).
      static constexpr size_t getCount()
      {
        return K + 1;
      }

      /// i-th GLL node in [-1, 1], 0 <= i <= K.
      static constexpr Real getNode(size_t i)
      {
        return s_nodes[i];
      }

      /// Full array of GLL nodes in ascending order.
      static constexpr const std::array<Real, K + 1>& getNodes()
      {
        return s_nodes;
      }

    private:
      static constexpr Real constexpr_cos(Real x)
      {
        constexpr Real pi = Math::Constants::pi();
        const Real twoPi = static_cast<Real>(2.0) * pi;
        while (x > pi)
          x -= twoPi;
        while (x < -pi)
          x += twoPi;

        const Real x2 = x * x;

        Real term = static_cast<Real>(1.0);
        Real sum  = term;

        // N-term Taylor series for cos(x)
        // cos x ≈ 1 - x^2/2! + x^4/4! - ... up to x^{2N}/(2N)!
        for (size_t n = 1; n <= 14; ++n)
        {
          term *= -x2 / (static_cast<Real>(2 * n) * static_cast<Real>(2 * n - 1));
          sum  += term;
        }

        return sum;
      }

      static constexpr std::array<Real, K + 1> compute()
      {
        std::array<Real, K + 1> nodes{};

        if constexpr (K == 0)
        {
          nodes[0] = static_cast<Real>(0.0);
          return nodes;
        }

        nodes[0] = static_cast<Real>(-1.0);
        nodes[K] = static_cast<Real>(1.0);

        if constexpr (K == 1)
        {
          return nodes;
        }

        constexpr size_t maxIt = MaxIt;

        for (size_t i = 1; i < K; ++i)
        {
          // Initial guess: Chebyshev point
          const Real ii   = static_cast<Real>(i);
          const Real kReal = static_cast<Real>(K);
          Real x = static_cast<Real>(-constexpr_cos(Math::Constants::pi() * ii / kReal));

          for (size_t it = 0; it < maxIt; ++it)
          {
            Real P, dP;
            LegendrePolynomial<K>::getValue(P, dP, x); // P_K(x), P'_K(x)

            const Real one_minus_x2 = static_cast<Real>(1.0) - x * x;

            // f(x) = (1 - x^2) P'_K(x)
            const Real f = one_minus_x2 * dP;

            // From Legendre ODE: f'(x) = -K(K+1) P_K(x)
            const Real df = -static_cast<Real>(K) * static_cast<Real>(K + 1) * P;

            const Real dx = -f / df;
            x += dx;

            if (Math::abs(dx) < RODIN_VARIATIONAL_H1_GLL_TOLERANCE)
              break;
          }

          nodes[i] = x;
        }

        return nodes;
      }

      inline static const std::array<Real, K + 1> s_nodes = compute();
  };

  /**
   * @brief Gauss–Lobatto–Legendre nodes mapped from [-1, 1] to [0, 1].
   *
   * Given the standard GLL nodes @f$ x_i \in [-1, 1] @f$, this class provides
   * the corresponding mapped nodes
   * @f[
   *   \hat{x}_i = \frac{1}{2}\,(x_i + 1) \in [0, 1],
   * @f]
   * preserving the ordering and indexing.
   *
   * @tparam K Polynomial degree.
   */
  template <size_t K>
  class GLL01
  {
    public:
      /// Number of GLL nodes (K + 1).
      static constexpr size_t getCount()
      {
        return K + 1;
      }

      /// i-th mapped GLL node in [0, 1], 0 <= i <= K.
      static constexpr Real getNode(size_t i)
      {
        return s_nodes[i];
      }

      /// Full array of mapped GLL nodes in ascending order.
      static constexpr const std::array<Real, K + 1>& getNodes()
      {
        return s_nodes;
      }

    private:
      static constexpr std::array<Real, K + 1> compute()
      {
        std::array<Real, K + 1> nodes{};

        const auto& gll = GLL<K>::getNodes();

        for (size_t i = 0; i < K + 1; ++i)
        {
          // Map from [-1, 1] to [0, 1]: x̂ = (x + 1)/2
          nodes[i] = static_cast<Real>(0.5) * (gll[i] + static_cast<Real>(1.0));
        }

        return nodes;
      }

      inline static const std::array<Real, K + 1> s_nodes = compute();
  };
}

#endif
