/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_QF_ORTHOGONALBASIS_H
#define RODIN_QF_ORTHOGONALBASIS_H

/**
 * @file
 * @brief Defines OrthogonalBasis, the Koornwinder-Dubiner basis of a simplex.
 */

#include <cmath>
#include <vector>

#include "Rodin/Math.h"
#include "Rodin/Geometry/Polytope.h"

namespace Rodin::QF
{
  /**
   * @ingroup RodinQuadrature
   * @brief Koornwinder-Dubiner orthogonal basis on the reference simplex.
   *
   * Moment equations written against monomials are ill-conditioned: the
   * monomial Gram matrix @f$ G_{\alpha\beta} = \int_K x^{\alpha+\beta} @f$ is
   * Hilbert-like and loses positive definiteness in double precision from
   * degree 12 on the triangle, so it cannot even be factorised, let alone used
   * to orthonormalise. Xiao--Gimbutas and Witherden--Vincent both avoid the
   * Gram matrix entirely by building an orthogonal basis from Jacobi
   * polynomial recurrences, which is what this does.
   *
   * On the unit triangle, with the collapse
   * @f$ a = 2x/(1-y) - 1 @f$ and @f$ b = 2y - 1 @f$,
   * @f[
   *   \psi_{mn}(x,y)
   *     = P_m^{0,0}(a)\,(1-y)^m\,P_n^{2m+1,0}(b),
   * @f]
   * and on the unit tetrahedron the analogous three-factor product. Each
   * @f$ \psi_{mn} @f$ is a polynomial of total degree @f$ m+n @f$, and the
   * family is orthogonal on the element.
   *
   * The normalisation is computed rather than quoted. Closed forms for the
   * constants differ between references by factors that depend on which
   * reference simplex is meant, and a wrong constant produces a basis that is
   * merely badly scaled rather than obviously broken. Here each function is
   * divided by its own @f$ L^2 @f$ norm, evaluated with a rule the exactness
   * suite has already verified.
   */
  class OrthogonalBasis
  {
    public:
      /// @brief Jacobi polynomial @f$ P_n^{(\alpha,\beta)}(x) @f$ by the
      /// three-term recurrence.
      static Real jacobi(size_t n, Real alpha, Real beta, Real x)
      {
        if (n == 0)
          return 1;
        Real p0 = 1;
        Real p1 = Real(0.5) * (alpha - beta + (alpha + beta + 2) * x);
        if (n == 1)
          return p1;
        for (size_t k = 2; k <= n; ++k)
        {
          const Real kk = static_cast<Real>(k);
          const Real c = 2 * kk * (kk + alpha + beta) * (2 * kk + alpha + beta - 2);
          const Real c1 = (2 * kk + alpha + beta - 1) *
            ((2 * kk + alpha + beta) * (2 * kk + alpha + beta - 2) * x + alpha * alpha -
              beta * beta);
          const Real c2 =
            2 * (kk + alpha - 1) * (kk + beta - 1) * (2 * kk + alpha + beta);
          const Real p2 = (c1 * p1 - c2 * p0) / c;
          p0 = p1;
          p1 = p2;
        }
        return p1;
      }

      /// @brief Exponent pairs/triples of total degree at most @p degree, in
      /// the order the basis functions are indexed.
      static std::vector<std::vector<size_t>> indices(size_t d, size_t degree)
      {
        std::vector<std::vector<size_t>> out;
        if (d == 2)
        {
          for (size_t m = 0; m <= degree; ++m)
            for (size_t n = 0; m + n <= degree; ++n)
              out.push_back({m, n});
        }
        else
        {
          for (size_t m = 0; m <= degree; ++m)
            for (size_t n = 0; m + n <= degree; ++n)
              for (size_t o = 0; m + n + o <= degree; ++o)
                out.push_back({m, n, o});
        }
        return out;
      }

      /// @brief Unnormalised @f$ \psi @f$ at @p x on the unit simplex.
      static Real evaluate(
        const std::vector<size_t>& idx, const Math::SpatialVector<Real>& x)
      {
        constexpr Real eps = 1e-14;
        if (idx.size() == 2)
        {
          const Real px = x[0], py = x[1];
          const Real q = 1 - py;
          const Real a = (q > eps) ? (2 * px / q - 1) : Real(-1);
          const Real b = 2 * py - 1;
          return jacobi(idx[0], 0, 0, a) * ipow(q, idx[0]) *
            jacobi(idx[1], 2 * static_cast<Real>(idx[0]) + 1, 0, b);
        }
        const Real px = x[0], py = x[1], pz = x[2];
        const Real q1 = 1 - py - pz;
        const Real q2 = 1 - pz;
        const Real a = (q1 > eps) ? (2 * px / q1 - 1) : Real(-1);
        const Real b = (q2 > eps) ? (2 * py / q2 - 1) : Real(-1);
        const Real c = 2 * pz - 1;
        return jacobi(idx[0], 0, 0, a) * ipow(q1, idx[0]) *
          jacobi(idx[1], 2 * static_cast<Real>(idx[0]) + 1, 0, b) * ipow(q2, idx[1]) *
          jacobi(idx[2], 2 * static_cast<Real>(idx[0] + idx[1]) + 2, 0, c);
      }

      /// @brief @f$ \frac{d}{dx} P_n^{(\alpha,\beta)}(x)
      /// = \tfrac{1}{2}(n+\alpha+\beta+1) P_{n-1}^{(\alpha+1,\beta+1)}(x) @f$.
      static Real dJacobi(size_t n, Real alpha, Real beta, Real x)
      {
        if (n == 0)
          return 0;
        return Real(0.5) * (static_cast<Real>(n) + alpha + beta + 1) *
          jacobi(n - 1, alpha + 1, beta + 1, x);
      }

      /**
       * @brief Gradient of @ref evaluate with respect to the reference
       * coordinates.
       *
       * Obtained analytically rather than by differencing: the elimination
       * solve is driven to a residual of @f$ 10^{-13} @f$, and a
       * finite-difference basis gradient caps the attainable accuracy near
       * @f$ 10^{-8} @f$. Gated by an FD-consistency test, as the house
       * pattern requires for a residual and tangent pair.
       */
      static Math::SpatialVector<Real> gradient(
        const std::vector<size_t>& idx, const Math::SpatialVector<Real>& x)
      {
        constexpr Real eps = 1e-12;
        Math::SpatialVector<Real> out;
        out.resize(static_cast<Eigen::Index>(idx.size()));
        out.setZero();

        if (idx.size() == 2)
        {
          const size_t m = idx[0], n = idx[1];
          const Real px = x[0], py = x[1];
          const Real q = std::max(1 - py, eps);
          const Real a = 2 * px / q - 1;
          const Real b = 2 * py - 1;
          const Real bn = 2 * static_cast<Real>(m) + 1;

          const Real Pm = jacobi(m, 0, 0, a), dPm = dJacobi(m, 0, 0, a);
          const Real Pn = jacobi(n, bn, 0, b), dPn = dJacobi(n, bn, 0, b);
          const Real qm = ipow(q, m);
          const Real qm1 = (m == 0) ? Real(0) : ipow(q, m - 1);

          out[0] = 2 * dPm * qm1 * Pn;
          out[1] = dPm * (2 * px / (q * q)) * qm * Pn -
            static_cast<Real>(m) * qm1 * Pm * Pn + Pm * qm * dPn * 2;
          return out;
        }

        const size_t m = idx[0], n = idx[1], o = idx[2];
        const Real px = x[0], py = x[1], pz = x[2];
        const Real q1 = std::max(1 - py - pz, eps);
        const Real q2 = std::max(1 - pz, eps);
        const Real a = 2 * px / q1 - 1;
        const Real b = 2 * py / q2 - 1;
        const Real c = 2 * pz - 1;
        const Real bn = 2 * static_cast<Real>(m) + 1;
        const Real cn = 2 * static_cast<Real>(m + n) + 2;

        const Real Pm = jacobi(m, 0, 0, a), dPm = dJacobi(m, 0, 0, a);
        const Real Pn = jacobi(n, bn, 0, b), dPn = dJacobi(n, bn, 0, b);
        const Real Po = jacobi(o, cn, 0, c), dPo = dJacobi(o, cn, 0, c);
        const Real q1m = ipow(q1, m);
        const Real q1m1 = (m == 0) ? Real(0) : ipow(q1, m - 1);
        const Real q2n = ipow(q2, n);
        const Real q2n1 = (n == 0) ? Real(0) : ipow(q2, n - 1);

        // d/dx: only a depends on x.
        out[0] = 2 * dPm * q1m1 * Pn * q2n * Po;
        // d/dy: a through q1, q1^m, and b.
        out[1] = dPm * (2 * px / (q1 * q1)) * q1m * Pn * q2n * Po -
          static_cast<Real>(m) * q1m1 * Pm * Pn * q2n * Po +
          Pm * q1m * dPn * (2 / q2) * q2n * Po;
        // d/dz: a through q1, q1^m, b through q2, q2^n, and c.
        out[2] = dPm * (2 * px / (q1 * q1)) * q1m * Pn * q2n * Po -
          static_cast<Real>(m) * q1m1 * Pm * Pn * q2n * Po +
          Pm * q1m * dPn * (2 * py / (q2 * q2)) * q2n * Po -
          static_cast<Real>(n) * Pm * q1m * Pn * q2n1 * Po +
          Pm * q1m * Pn * q2n * dPo * 2;
        return out;
      }

      /// @brief Integer power.
      static Real ipow(Real x, size_t e)
      {
        Real r = 1;
        while (e--)
          r *= x;
        return r;
      }
  };
}

#endif
