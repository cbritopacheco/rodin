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
  template <size_t K>
  class DubinerTriangle
  {
    public:
      template <size_t P, size_t Q>
      static constexpr void getBasis(Real& basis, Real r, Real s)
      {
        static_assert(P + Q <= K, "DubinerTriangle: P + Q must be <= K.");

        Real b = s;
        Real a;

        // Handle singularity at s = 1
        if (std::abs(s - 1.0) < RODIN_VARIATIONAL_H1_DUBINER_TOLERANCE)
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

      template <size_t P, size_t Q>
      static constexpr void getGradient(Real& dpsi_dr, Real& dpsi_ds, Real r, Real s)
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

        Real scale_b = Math::pow(0.5 * (1.0 - b), P);

        // ∂ψ / ∂r: only a depends on r
        dpsi_dr = 0.0;
        if (Math::abs(s - 1.0) > RODIN_VARIATIONAL_H1_DUBINER_TOLERANCE)
        {
          const Real da_dr = 2.0 / (1.0 - s);
          dpsi_dr = dPa * da_dr * Pb * scale_b;
        }

        // ∂ψ / ∂s: both a and b depend on s, and scale_b depends on b = s
        dpsi_ds = 0.0;
        if (Math::abs(s - 1.0) > RODIN_VARIATIONAL_H1_DUBINER_TOLERANCE)
        {
          const Real da_ds = -2.0 * (1.0 + r) / ((1.0 - s) * (1.0 - s));
          dpsi_ds += dPa * da_ds * Pb * scale_b;
        }

        dpsi_ds += Pa * dPb * scale_b;  // db_ds = 1

        if constexpr (P > 0)
          dpsi_ds += Pa * Pb * P * Math::pow(0.5 * (1.0 - b), P - 1) * (-0.5);
      }

      static constexpr void getCollapsed(Real& r, Real& s, Real x, Real y)
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

  template <size_t K, size_t MaxItGLL = 25>
  class VandermondeTriangle
  {
    public:
      static const Math::Matrix<Real>& getMatrix()
      {
        static thread_local Math::Matrix<Real> s_vandermonde;

        constexpr size_t N = FeketeTriangle<K>::Count;

        if (s_vandermonde.size() == 0)
        {
          const auto& nodes = FeketeTriangle<K>::getNodes();
          s_vandermonde.resize(N, N);

          // Fill Vandermonde matrix
          size_t mode_idx = 0;
          Rodin::Utility::ForIndex<K + 1>(
              [&](auto p_idx)
              {
                constexpr size_t P = p_idx.value;
                Rodin::Utility::ForIndex<K + 1 - P>(
                    [&](auto q_idx)
                    {
                      constexpr size_t Q = q_idx.value;
                      for (size_t node_idx = 0; node_idx < N; ++node_idx)
                      {
                        Real r, s;
                        DubinerTriangle<K>::getCollapsed(
                            r, s, nodes[node_idx].x(), nodes[node_idx].y());

                        DubinerTriangle<K>::template getBasis<P, Q>(
                            s_vandermonde(node_idx, mode_idx), r, s);
                      }
                      ++mode_idx;
                    });
              });
        }
        return s_vandermonde;
      }

      static const Math::Matrix<Real>& getInverse()
      {
        static thread_local Math::Matrix<Real> s_inv;

        if (s_inv.size() == 0)
        {
          const auto& V = VandermondeTriangle<K>::getMatrix();
          assert(V.rows() == V.cols() && "Triangle Vandermonde must be square.");
          Eigen::FullPivLU<Math::Matrix<Real>> lu(V);
          assert(lu.isInvertible());
          s_inv = lu.inverse();
        }

        return s_inv;
      }
  };

  // ---------------------------------------------------------------------------
  // Tetrahedron version
  // ---------------------------------------------------------------------------

  template <size_t K>
  class DubinerTetrahedron
  {
    public:
      // ψ_{P,Q,R}(a,b,c) on [-1,1]^3 (Dubiner tetrahedral basis)
      template <size_t P, size_t Q, size_t R>
      static constexpr void getBasis(Real& basis, Real a, Real b, Real c)
      {
        static_assert(P + Q + R <= K, "DubinerTetrahedron: P + Q + R must be <= K.");

        Real _unused;

        Real P_a;
        JacobiPolynomial<P>::getValue(P_a, _unused, 0.0, 0.0, a);

        Real P_b;
        JacobiPolynomial<Q>::getValue(P_b, _unused, 2.0 * P + 1.0, 0.0, b);

        Real P_c;
        JacobiPolynomial<R>::getValue(P_c, _unused, 2.0 * P + 2.0 * Q + 2.0, 0.0, c);

        Real scale_b = Math::pow(0.5 * (1.0 - b), P);
        Real scale_c = Math::pow(0.5 * (1.0 - c), P + Q);

        basis = P_a * P_b * P_c * scale_b * scale_c;
      }

      template <size_t P, size_t Q, size_t R>
      static constexpr void getGradient(Real& dpsi_da,
                                        Real& dpsi_db,
                                        Real& dpsi_dc,
                                        Real a, Real b, Real c)
      {
        Real P_a, dP_a;
        JacobiPolynomial<P>::getValue(P_a, dP_a, 0.0, 0.0, a);

        Real P_b, dP_b;
        JacobiPolynomial<Q>::getValue(P_b, dP_b, 2.0 * P + 1.0, 0.0, b);

        Real P_c, dP_c;
        JacobiPolynomial<R>::getValue(P_c, dP_c, 2.0 * P + 2.0 * Q + 2.0, 0.0, c);

        const Real scale_b = Math::pow(0.5 * (1.0 - b), P);
        const Real scale_c = Math::pow(0.5 * (1.0 - c), P + Q);

        // ∂ψ / ∂a
        dpsi_da = dP_a * P_b * P_c * scale_b * scale_c;

        // ∂ψ / ∂b
        Real dscale_b_db = 0.0;
        if constexpr (P > 0)
          dscale_b_db = P * Math::pow(0.5 * (1.0 - b), P - 1) * (-0.5);

        dpsi_db = P_a * (dP_b * scale_b + P_b * dscale_b_db) * P_c * scale_c;

        // ∂ψ / ∂c
        Real dscale_c_dc = 0.0;
        if constexpr (P + Q > 0)
          dscale_c_dc = (P + Q) * Math::pow(0.5 * (1.0 - c), P + Q - 1) * (-0.5);

        dpsi_dc = P_a * P_b * (dP_c * scale_c + P_c * dscale_c_dc) * scale_b;
      }

      // Map reference tetra (0,0,0)-(1,0,0)-(0,1,0)-(0,0,1) → (a,b,c) ∈ [-1,1]^3
      // Using a Duffy-type collapse:
      //
      //   c = 2 z - 1
      //   b = 2 y / (1 - z) - 1          if 1 - z > tol
      //   a = 2 x / (1 - y - z) - 1      if 1 - y - z > tol
      //
      static constexpr void getCollapsed(Real& a,
                                         Real& b,
                                         Real& c,
                                         Real x,
                                         Real y,
                                         Real z)
      {
        c = 2.0 * z - 1.0;

        const Real one_minus_z = 1.0 - z;
        if (one_minus_z > RODIN_VARIATIONAL_H1_DUBINER_TOLERANCE)
        {
          b = 2.0 * (y / one_minus_z) - 1.0;

          const Real one_minus_yz = 1.0 - y - z;
          if (one_minus_yz > RODIN_VARIATIONAL_H1_DUBINER_TOLERANCE)
            a = 2.0 * (x / one_minus_yz) - 1.0;
          else
            a = -1.0; // collapse along edge
        }
        else
        {
          // collapse along vertex line
          b = -1.0;
          a = -1.0;
        }
      }
  };

  template <size_t K>
  class VandermondeTetrahedron
  {
    public:
      static const Math::Matrix<Real>& getMatrix()
      {
        static thread_local Math::Matrix<Real> s_vandermonde;

        constexpr size_t N = FeketeTetrahedron<K>::Count;

        if (s_vandermonde.size() == 0)
        {
          const auto& nodes = FeketeTetrahedron<K>::getNodes();
          s_vandermonde.resize(N, N);

          size_t mode_idx = 0;
          Rodin::Utility::ForIndex<K + 1>(
              [&](auto p_idx)
              {
                constexpr size_t P = p_idx.value;
                Rodin::Utility::ForIndex<K + 1 - P>(
                    [&](auto q_idx)
                    {
                      constexpr size_t Q = q_idx.value;
                      Rodin::Utility::ForIndex<K + 1 - P - Q>(
                          [&](auto r_idx)
                          {
                            constexpr size_t R = r_idx.value;
                            for (size_t node_idx = 0; node_idx < N; ++node_idx)
                            {
                              Real a, b, c;
                              DubinerTetrahedron<K>::getCollapsed(
                                  a, b, c,
                                  nodes[node_idx].x(),
                                  nodes[node_idx].y(),
                                  nodes[node_idx].z());

                              DubinerTetrahedron<K>::template getBasis<P, Q, R>(
                                  s_vandermonde(node_idx, mode_idx), a, b, c);
                            }
                            ++mode_idx;
                          });
                    });
              });
        }

        return s_vandermonde;
      }

      static const Math::Matrix<Real>& getInverse()
      {
        static thread_local Math::Matrix<Real> s_inv;

        if (s_inv.size() == 0)
        {
          const auto& V = VandermondeTetrahedron<K>::getMatrix();
          assert(V.rows() == V.cols() && "Tetrahedron Vandermonde must be square.");
          Eigen::FullPivLU<Math::Matrix<Real>> lu(V);
          assert(lu.isInvertible());
          s_inv = lu.inverse();
        }

        return s_inv;
      }
  };
}

#endif
