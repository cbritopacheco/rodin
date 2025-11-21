/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_DUBINER_H
#define RODIN_VARIATIONAL_H1_DUBINER_H

#include <cstddef>

#include "Rodin/Math/Common.h"
#include "Rodin/Math/Matrix.h"

#include "Rodin/Utility/ForConstexpr.h"

#include "Fekete.h"
#include "JacobiPolynomial.h"

#define RODIN_VARIATIONAL_H1_DUBINER_TOLERANCE 1e-14

namespace Rodin::Variational
{
  template <size_t K, size_t P, size_t Q, size_t MaxItGLL = 25>
  class DubinerTriangle
  {
    static_assert(P + Q <= K, "DubinerTriangle: P + Q must be <= K.");

    public:
      void getBasis(Real& basis, Real r, Real s)
      {
        Real b = s;
        Real a;

        // Handle singularity at s = 1
        if (Math::abs(s - 1.0) < RODIN_VARIATIONAL_H1_DUBINER_TOLERANCE)
          a = -1.0;
        else
          a = 2.0 * (1.0 + r) / (1.0 - s) - 1.0;

        Real _unused;

        Real psi_p;
        JacobiPolynomial<P>::getValue(psi_p, _unused, 0.0, 0.0, a);

        Real psi_q;
        JacobiPolynomial<Q>::getValue(psi_q, _unused, 2.0 * P + 1.0, 0.0, b);

        basis = psi_p * psi_q * Math::pow(0.5 * (1.0 - b), P);
      }

      void getGradient(Real& dpsi_dr, Real& dpsi_ds, Real r, Real s)
      {
        Real b = s;
        Real a;

        if (Math::abs(s - 1.0) < RODIN_VARIATIONAL_H1_DUBINER_TOLERANCE)
          a = -1.0;
        else
          a = 2.0 * (1.0 + r) / (1.0 - s) - 1.0;

        Real Pa, dPa;
        JacobiPolynomial<P>::getValue(Pa, dPa, 0.0, 0.0, a);

        Real Pb, dPb;
        JacobiPolynomial<Q>::getValue(Pb, dPb, 2.0 * P + 1.0, 0.0, b);

        Real scale_b = std::pow(0.5 * (1.0 - b), P);

        // ∂ψ / ∂r: only a depends on r
        dpsi_dr = 0.0;
        if (Math::abs(s - 1.0) > RODIN_VARIATIONAL_H1_DUBINER_TOLERANCE)
        {
          Real da_dr = 2.0 / (1.0 - s);
          dpsi_dr = dPa * da_dr * Pb * scale_b;
        }

        // ∂ψ / ∂s: both a and b depend on s, and scale_b depends on b = s
        dpsi_ds = 0.0;
        if (Math::abs(s - 1.0) > RODIN_VARIATIONAL_H1_DUBINER_TOLERANCE)
        {
          Real da_ds = -2.0 * (1.0 + r) / ((1.0 - s) * (1.0 - s));
          dpsi_ds += dPa * da_ds * Pb * scale_b;
        }

        dpsi_ds += Pa * dPb * scale_b;  // db_ds = 1

        if (P > 0)
        {
          dpsi_ds += Pa * Pb * P * Math::pow(0.5 * (1.0 - b), P - 1) * (-0.5);
        }
      }

      constexpr
      static void getCollapsed(Real& r, Real& s, Real x, Real y)
      {
        // s = 2y - 1 always
        s = 2.0 * y - 1.0;

        // r = 2x / (1 - y) - 1, with collapse at the top edge
        if (1.0 - y > RODIN_VARIATIONAL_H1_DUBINER_TOLERANCE)
          r = 2.0 * (x / (1.0 - y)) - 1.0;
        else
          r = -1.0; // collapse at top vertex/edge
      }
  };

  template <size_t K, size_t MaxItGLL>
  class VandermondeTriangle
  {
    public:
      static const Math::Matrix<Real>& getMatrix()
      {
        static thread_local Math::Matrix<Real> s_vandermonde;

        constexpr size_t N = FeketeTriangle<K, MaxItGLL>::Count;

        if (s_vandermonde.size() == 0)
        {
          const auto& nodes = FeketeTriangle<K, MaxItGLL>::getNodes();
          s_vandermonde.resize(N, N);

          // Fill Vandermonde matrix
          size_t mode_idx = 0;
          Rodin::Utility::ForIndex<K + 1>([&](auto p_idx)
              {
                constexpr size_t P = p_idx.value;

                Rodin::Utility::ForIndex<K + 1 - P>([&](auto q_idx)
                {
                  constexpr size_t Q = q_idx.value;

                  for (size_t node_idx = 0; node_idx < N; ++node_idx)
                  {
                    Real r, s;

                    DubinerTriangle<K, P, Q, MaxItGLL>::getCollapsed(
                        r, s, nodes[node_idx].x(), nodes[node_idx].y());

                    DubinerTriangle<K, P, Q, MaxItGLL>::getBasis(
                        s_vandermonde(node_idx, mode_idx), r, s);
                  }

                  ++mode_idx;
                });
              });
        }
        return s_vandermonde;
      }
  };
}

#endif
