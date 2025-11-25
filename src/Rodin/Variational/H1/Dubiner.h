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
  /**
   * @brief Dubiner orthogonal modal basis on the reference triangle.
   *
   * This class provides evaluation of the Dubiner basis functions
   * @f$ \psi_{p,q}(r,s) @f$ on the collapsed coordinate system
   * @f$ (r,s) \in [-1,1]^2 @f$, which maps from the reference triangle
   * with vertices @f$ (0,0), (1,0), (0,1) @f$.
   *
   * The Dubiner basis is constructed as a product of Jacobi polynomials:
   * @f[
   *   \psi_{p,q}(r,s) = P_p^{(0,0)}(a) \cdot P_q^{(2p+1,0)}(b)
   *     \cdot \left(\frac{1-b}{2}\right)^p
   * @f]
   * where @f$ (a,b) @f$ are collapsed coordinates derived from @f$ (r,s) @f$.
   *
   * This basis is L2-orthogonal on the reference triangle, making it
   * numerically stable for high-order polynomial approximations.
   *
   * @tparam K Maximum polynomial degree.
   *
   * @see JacobiPolynomial for the underlying polynomial evaluation.
   * @see VandermondeTriangle for the nodal-to-modal transformation.
   */
  template <size_t K>
  class DubinerTriangle
  {
    public:
      /**
       * @brief Evaluates the Dubiner basis function @f$ \psi_{P,Q}(r,s) @f$.
       *
       * @tparam P First modal index (0 ≤ P).
       * @tparam Q Second modal index (0 ≤ Q, P + Q ≤ K).
       * @param[out] basis The computed basis function value.
       * @param r First collapsed coordinate in [-1,1].
       * @param s Second collapsed coordinate in [-1,1].
       */
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

      /**
       * @brief Computes the gradient of @f$ \psi_{P,Q} @f$ w.r.t. collapsed coordinates.
       *
       * Uses the chain rule to compute:
       * @f[
       *   \frac{\partial \psi}{\partial r}, \quad \frac{\partial \psi}{\partial s}
       * @f]
       *
       * @tparam P First modal index.
       * @tparam Q Second modal index.
       * @param[out] dpsi_dr Derivative w.r.t. r.
       * @param[out] dpsi_ds Derivative w.r.t. s.
       * @param r First collapsed coordinate.
       * @param s Second collapsed coordinate.
       */
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

      /**
       * @brief Converts reference triangle coordinates to collapsed coordinates.
       *
       * Maps @f$ (x,y) @f$ from the reference triangle to @f$ (r,s) \in [-1,1]^2 @f$:
       * @f[
       *   s = 2y - 1, \quad r = \frac{2x}{1-y} - 1
       * @f]
       * with singularity handling at @f$ y = 1 @f$ (top vertex).
       *
       * @param[out] r First collapsed coordinate.
       * @param[out] s Second collapsed coordinate.
       * @param x First reference coordinate.
       * @param y Second reference coordinate.
       */
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

  /**
   * @brief Vandermonde matrix for nodal-to-modal transformation on triangles.
   *
   * Constructs and caches the Vandermonde matrix @f$ V @f$ where
   * @f$ V_{ij} = \psi_j(x_i) @f$ with @f$ \psi_j @f$ being the j-th Dubiner
   * mode and @f$ x_i @f$ the i-th Fekete node.
   *
   * The nodal basis functions @f$ \phi_i @f$ are recovered from the modal
   * basis via @f$ \phi_i(x) = \sum_j V^{-1}_{ij} \psi_j(x) @f$.
   *
   * Both @f$ V @f$ and @f$ V^{-1} @f$ are cached in thread-local storage.
   *
   * @tparam K Polynomial degree.
   *
   * @see DubinerTriangle for the modal basis functions.
   * @see FeketeTriangle for the nodal points.
   */
  template <size_t K>
  class VandermondeTriangle
  {
    public:
      /**
       * @brief Returns the Vandermonde matrix @f$ V @f$.
       *
       * The matrix is computed once and cached in thread-local storage.
       *
       * @return Reference to the N×N Vandermonde matrix where N = (K+1)(K+2)/2.
       */
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

  /**
   * @brief Dubiner orthogonal modal basis on the reference tetrahedron.
   *
   * Provides evaluation of the Dubiner basis functions
   * @f$ \psi_{p,q,r}(a,b,c) @f$ on collapsed coordinates
   * @f$ (a,b,c) \in [-1,1]^3 @f$, mapped from the reference tetrahedron
   * with vertices @f$ (0,0,0), (1,0,0), (0,1,0), (0,0,1) @f$.
   *
   * The basis uses a Duffy-type collapse and is constructed as:
   * @f[
   *   \psi_{p,q,r}(a,b,c) = P_p^{(0,0)}(a) \cdot P_q^{(2p+1,0)}(b)
   *     \cdot P_r^{(2p+2q+2,0)}(c) \cdot \left(\frac{1-b}{2}\right)^p
   *     \cdot \left(\frac{1-c}{2}\right)^{p+q}
   * @f]
   *
   * @tparam K Maximum polynomial degree.
   *
   * @see JacobiPolynomial, VandermondeTetrahedron
   */
  template <size_t K>
  class DubinerTetrahedron
  {
    public:
      /**
       * @brief Evaluates the Dubiner basis function @f$ \psi_{P,Q,R}(a,b,c) @f$.
       *
       * @tparam P First modal index.
       * @tparam Q Second modal index.
       * @tparam R Third modal index (P + Q + R ≤ K).
       * @param[out] basis The computed basis function value.
       * @param a First collapsed coordinate.
       * @param b Second collapsed coordinate.
       * @param c Third collapsed coordinate.
       */
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
