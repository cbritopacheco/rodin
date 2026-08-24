/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Defines CollapsedBasis, one orthogonal basis for every element.
 */
#ifndef RODIN_QF_COLLAPSEDBASIS_H
#define RODIN_QF_COLLAPSEDBASIS_H

#include <cassert>
#include <cmath>
#include <vector>

#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Types.h"

#include "GaussLegendre.h"

namespace Rodin::QF
{
  /**
   * @brief Orthogonal polynomial basis of any reference element, as a warped
   * tensor product of Jacobi polynomials.
   *
   * Every reference element here --- segment, triangle, quadrilateral,
   * tetrahedron, wedge, pyramid, hexahedron --- is the image of a cube under a
   * collapse, and on each of them the orthogonal basis is the same object:
   * @f[
   *   \psi_{\mathbf{i}}(x) = \prod_m
   *     \hat{P}^{(\alpha_m, 0)}_{i_m}(c_m)
   *     \left( \frac{1 - c_m}{2} \right)^{s_m},
   *   \qquad
   *   s_m = \sum_{l \in S_m} i_l,
   *   \qquad
   *   \alpha_m = 2 s_m + o_m .
   * @f]
   * The elements differ only in the collapse @f$ c_m(x) @f$ and in which
   * preceding indices raise each exponent --- two short tables, not seven code
   * paths. This is the basis of Karniadakis and Sherwin, and the one
   * Witherden and Vincent tabulate domain by domain
   * @cite witherden2015identification; writing it once is what lets a single
   * rule generator serve every element.
   *
   * The basis is orthogonal but not normalised: the normalising constants have
   * closed forms that differ per element and per index, and getting one of them
   * subtly wrong yields a basis that still looks plausible while quietly
   * spanning the wrong thing. They are instead measured with exactRule(), which
   * is exact and positive, and the measurement is checked by a test asserting
   * the Gram matrix is the identity.
   *
   * @see exactRule
   */
  class CollapsedBasis
  {
    public:
      /**
       * @brief How one collapsed direction is built.
       *
       * @p numerator is the coordinate the direction maps, @p denominator the
       * coordinates subtracted from one to form the shrinking factor, and
       * @p contributors the preceding directions whose indices raise this
       * one's Jacobi exponent by @p offset.
       */
      struct Axis
      {
          size_t numerator = 0;
          std::vector<size_t> denominator;
          std::vector<size_t> contributors;
          size_t offset = 0;
      };

      /// @brief The collapse of a reference element, one entry per direction.
      static const std::vector<Axis>& axes(Geometry::Polytope::Type g)
      {
        static const std::vector<Axis> segment = {{0, {}, {}, 0}};
        static const std::vector<Axis> triangle = {{0, {1}, {}, 0}, {1, {}, {0}, 1}};
        static const std::vector<Axis> quadrilateral = {{0, {}, {}, 0}, {1, {}, {}, 0}};
        static const std::vector<Axis> tetrahedron = {
          {0, {1, 2}, {}, 0}, {1, {2}, {0}, 1}, {2, {}, {0, 1}, 2}};
        static const std::vector<Axis> wedge = {
          {0, {1}, {}, 0}, {1, {}, {0}, 1}, {2, {}, {}, 0}};
        static const std::vector<Axis> pyramid = {
          {0, {2}, {}, 0}, {1, {2}, {}, 0}, {2, {}, {0, 1}, 2}};
        static const std::vector<Axis> hexahedron = {
          {0, {}, {}, 0}, {1, {}, {}, 0}, {2, {}, {}, 0}};
        static const std::vector<Axis> none;

        switch (g)
        {
          case Geometry::Polytope::Type::Segment:
            return segment;
          case Geometry::Polytope::Type::Triangle:
            return triangle;
          case Geometry::Polytope::Type::Quadrilateral:
            return quadrilateral;
          case Geometry::Polytope::Type::Tetrahedron:
            return tetrahedron;
          case Geometry::Polytope::Type::Wedge:
            return wedge;
          case Geometry::Polytope::Type::Pyramid:
            return pyramid;
          case Geometry::Polytope::Type::Hexahedron:
            return hexahedron;
          default:
            return none;
        }
      }

      /**
       * @brief Index tuples spanning the polynomials of total degree at most
       * @p degree.
       *
       * Each @f$ \psi_{\mathbf{i}} @f$ has total degree @f$ \sum_m i_m @f$ on
       * every one of these elements, the shrinking factors cancelling the
       * poles the collapse introduces, so the same bound serves all of them.
       */
      static std::vector<std::vector<size_t>> indices(
        Geometry::Polytope::Type g, size_t degree)
      {
        const size_t d = axes(g).size();
        std::vector<std::vector<size_t>> out;
        std::vector<size_t> current(d, 0);
        const std::function<void(size_t, size_t)> rec = [&](size_t axis, size_t budget) {
          if (axis == d)
          {
            out.push_back(current);
            return;
          }
          for (size_t i = 0; i <= budget; ++i)
          {
            current[axis] = i;
            rec(axis + 1, budget - i);
          }
          current[axis] = 0;
        };
        rec(0, degree);
        return out;
      }

      /**
       * @brief The collapsed coordinates at @p x, and their derivatives.
       *
       * @param[out] c Jacobi arguments, each in @f$ [-1, 1] @f$.
       * @param[out] shrink The factor @f$ (1 - c_m)/2 @f$ of each direction.
       * @param[out] dc Derivative of @p c with respect to the reference
       *   coordinates, or ignored when null.
       */
      static void collapse(Geometry::Polytope::Type g, const Math::SpatialVector<Real>& x,
        std::vector<Real>& c, std::vector<Real>& shrink, Math::Matrix<Real>* dc)
      {
        constexpr Real eps = 1e-14;
        const auto& all = axes(g);
        const size_t d = all.size();
        c.assign(d, 0);
        shrink.assign(d, 0);
        if (dc)
          *dc = Math::Matrix<Real>::Zero(
            static_cast<Eigen::Index>(d), static_cast<Eigen::Index>(d));

        for (size_t m = 0; m < d; ++m)
        {
          const auto& axis = all[m];
          Real q = 1;
          for (const size_t l : axis.denominator)
            q -= x[static_cast<Eigen::Index>(l)];
          const Real t = x[static_cast<Eigen::Index>(axis.numerator)];
          if (q <= eps)
          {
            // The apex of the collapse, where the direction degenerates. The
            // limit is the standard one; no quadrature point of a positive
            // interior rule lands here.
            c[m] = -1;
            shrink[m] = 1;
            continue;
          }
          c[m] = 2 * t / q - 1;
          shrink[m] = (1 - c[m]) / 2;
          if (!dc)
            continue;
          const Eigen::Index row = static_cast<Eigen::Index>(m);
          (*dc)(row, static_cast<Eigen::Index>(axis.numerator)) += 2 / q;
          for (const size_t l : axis.denominator)
            (*dc)(row, static_cast<Eigen::Index>(l)) += 2 * t / (q * q);
        }
      }

      /**
       * @brief Every Jacobi factor needed at one point, computed once.
       *
       * Evaluating a basis function on its own repeats the collapse and runs
       * the Jacobi recurrence from scratch, so evaluating a whole basis that
       * way costs @f$ O(n_e p) @f$ per point. The factors depend only on the
       * direction, the exponent, and the degree, so tabulating them once
       * leaves @f$ O(p^2) @f$ of setup and @f$ O(d) @f$ per mode --- which is
       * the difference between a search that finishes at degree twenty and one
       * that does not.
       *
       * @see tabulate
       */
      struct Tabulation
      {
          /// @brief value[m][s][i] is @f$ P^{(2s + o_m, 0)}_i(c_m) @f$.
          std::vector<std::vector<std::vector<Real>>> value;
          /// @brief Derivative of @ref value in @f$ c_m @f$.
          std::vector<std::vector<std::vector<Real>>> derivative;
          /// @brief shrink[m][s] is @f$ ((1 - c_m)/2)^s @f$.
          std::vector<std::vector<Real>> shrink;
          /// @brief Derivative of @ref shrink in @f$ c_m @f$.
          std::vector<std::vector<Real>> dshrink;
          /// @brief Derivative of the collapsed coordinates in @p x.
          Math::Matrix<Real> dc;
      };

      /// @brief Tabulates every factor needed to evaluate modes of total
      /// degree at most @p degree at @p x.
      static Tabulation tabulate(Geometry::Polytope::Type g, size_t degree,
        const Math::SpatialVector<Real>& x, bool wantDerivatives)
      {
        const auto& all = axes(g);
        const size_t d = all.size();
        std::vector<Real> c, shrink;
        Tabulation out;
        collapse(g, x, c, shrink, wantDerivatives ? &out.dc : nullptr);

        out.value.resize(d);
        out.derivative.resize(d);
        out.shrink.resize(d);
        out.dshrink.resize(d);
        for (size_t m = 0; m < d; ++m)
        {
          out.value[m].resize(degree + 1);
          out.derivative[m].resize(degree + 1);
          for (size_t s = 0; s <= degree; ++s)
          {
            const Real alpha =
              2 * static_cast<Real>(s) + static_cast<Real>(all[m].offset);
            out.value[m][s].resize(degree + 1 - s);
            out.derivative[m][s].assign(degree + 1 - s, 0);
            for (size_t i = 0; i + s <= degree; ++i)
            {
              out.value[m][s][i] = jacobi(i, alpha, 0, c[m]);
              if (wantDerivatives)
                out.derivative[m][s][i] = dJacobi(i, alpha, 0, c[m]);
            }
          }

          out.shrink[m].resize(degree + 1);
          out.dshrink[m].assign(degree + 1, 0);
          Real power = 1;
          for (size_t s = 0; s <= degree; ++s)
          {
            out.shrink[m][s] = power;
            if (wantDerivatives && s > 0)
              out.dshrink[m][s] =
                -static_cast<Real>(s) * ((s > 0) ? out.shrink[m][s - 1] : Real(0)) / 2;
            power *= shrink[m];
          }
        }
        return out;
      }

      /// @brief Value of the mode @p index from a tabulation.
      static Real evaluate(Geometry::Polytope::Type g, const Tabulation& table,
        const std::vector<size_t>& index)
      {
        const auto& all = axes(g);
        Real out = 1;
        for (size_t m = 0; m < all.size(); ++m)
        {
          const size_t s = power(all[m], index);
          out *= table.value[m][s][index[m]] * table.shrink[m][s];
        }
        return out;
      }

      /// @brief Gradient of the mode @p index from a tabulation.
      static Math::SpatialVector<Real> gradient(Geometry::Polytope::Type g,
        const Tabulation& table, const std::vector<size_t>& index)
      {
        const auto& all = axes(g);
        const size_t d = all.size();
        std::vector<Real> factor(d), dfactor(d);
        for (size_t m = 0; m < d; ++m)
        {
          const size_t s = power(all[m], index);
          const Real P = table.value[m][s][index[m]];
          const Real dP = table.derivative[m][s][index[m]];
          factor[m] = P * table.shrink[m][s];
          dfactor[m] = dP * table.shrink[m][s] + P * table.dshrink[m][s];
        }

        Math::SpatialVector<Real> out;
        out.resize(static_cast<Eigen::Index>(d));
        out.setZero();
        for (size_t m = 0; m < d; ++m)
        {
          Real rest = 1;
          for (size_t n = 0; n < d; ++n)
            if (n != m)
              rest *= factor[n];
          const Real scale = dfactor[m] * rest;
          for (size_t k = 0; k < d; ++k)
            out[static_cast<Eigen::Index>(k)] += scale *
              table.dc(static_cast<Eigen::Index>(m), static_cast<Eigen::Index>(k));
        }
        return out;
      }

      /// @brief Value of the basis function @p index at @p x.
      static Real evaluate(Geometry::Polytope::Type g, const std::vector<size_t>& index,
        const Math::SpatialVector<Real>& x)
      {
        const auto& all = axes(g);
        std::vector<Real> c, shrink;
        collapse(g, x, c, shrink, nullptr);
        Real out = 1;
        for (size_t m = 0; m < all.size(); ++m)
        {
          const Real alpha = exponent(all[m], index);
          out *= jacobi(index[m], alpha, 0, c[m]);
          out *= ipow(shrink[m], power(all[m], index));
        }
        return out;
      }

      /**
       * @brief Gradient of @ref evaluate with respect to the reference
       * coordinates.
       *
       * Analytic rather than differenced: the rule solvers this feeds are
       * driven to residuals near @f$ 10^{-15} @f$, and a differenced basis
       * gradient caps them several digits short of that.
       */
      static Math::SpatialVector<Real> gradient(Geometry::Polytope::Type g,
        const std::vector<size_t>& index, const Math::SpatialVector<Real>& x)
      {
        const auto& all = axes(g);
        const size_t d = all.size();
        std::vector<Real> c, shrink;
        Math::Matrix<Real> dc;
        collapse(g, x, c, shrink, &dc);

        // Factor of each direction, and its derivative in that direction's own
        // collapsed coordinate. The chain rule then maps them back to x.
        std::vector<Real> factor(d), dfactor(d);
        for (size_t m = 0; m < d; ++m)
        {
          const Real alpha = exponent(all[m], index);
          const size_t s = power(all[m], index);
          const Real P = jacobi(index[m], alpha, 0, c[m]);
          const Real dP = dJacobi(index[m], alpha, 0, c[m]);
          const Real w = ipow(shrink[m], s);
          // d/dc of ((1 - c)/2)^s is -s ((1 - c)/2)^(s - 1) / 2.
          const Real dw =
            (s == 0) ? Real(0) : -static_cast<Real>(s) * ipow(shrink[m], s - 1) / 2;
          factor[m] = P * w;
          dfactor[m] = dP * w + P * dw;
        }

        Math::SpatialVector<Real> out;
        out.resize(static_cast<Eigen::Index>(d));
        out.setZero();
        for (size_t m = 0; m < d; ++m)
        {
          Real rest = 1;
          for (size_t n = 0; n < d; ++n)
            if (n != m)
              rest *= factor[n];
          const Real scale = dfactor[m] * rest;
          for (size_t k = 0; k < d; ++k)
            out[static_cast<Eigen::Index>(k)] +=
              scale * dc(static_cast<Eigen::Index>(m), static_cast<Eigen::Index>(k));
        }
        return out;
      }

      /**
       * @brief A collapsed Gauss--Jacobi rule on @p g, exact for polynomials
       * of degree @f$ 2n - 1 @f$.
       *
       * The element is swept as a degenerated cube, the Jacobian of the
       * collapse carried exactly by the Jacobi weights rather than evaluated
       * at the nodes. Every weight is positive, which is what makes it a safe
       * instrument for measuring a basis --- unlike a rule of alternating
       * sign, whose cancellation at these degrees swamps what is being
       * measured.
       */
      static void exactRule(Geometry::Polytope::Type g, size_t n,
        std::vector<Math::SpatialVector<Real>>& points, std::vector<Real>& weights)
      {
        const auto& all = axes(g);
        const size_t d = all.size();

        // Each direction is integrated against the power of the shrinking
        // factor that the collapse leaves behind, which is exactly the number
        // of later directions that divide by it.
        std::vector<std::vector<Real>> node(d), weight(d);
        for (size_t m = 0; m < d; ++m)
        {
          size_t collapsed = 0;
          for (size_t l = 0; l < d; ++l)
            for (const size_t q : all[l].denominator)
              if (q == all[m].numerator)
                ++collapsed;
          if (collapsed == 0)
            GaussLegendre::gl1dUnit(n, node[m], weight[m]);
          else
            GaussLegendre::gj1dUnit(n, collapsed, node[m], weight[m]);
        }

        points.clear();
        weights.clear();
        std::vector<size_t> at(d, 0);
        const size_t total = static_cast<size_t>(std::pow(n, d) + 0.5);
        for (size_t flat = 0; flat < total; ++flat)
        {
          size_t rest = flat;
          for (size_t m = 0; m < d; ++m)
          {
            at[m] = rest % n;
            rest /= n;
          }
          // Map the cube coordinates back through the collapse. A direction is
          // scaled by every shrinking factor that divides it.
          Math::SpatialVector<Real> x;
          x.resize(static_cast<Eigen::Index>(d));
          Real w = 1;
          for (size_t m = 0; m < d; ++m)
            w *= weight[m][at[m]];
          for (size_t m = 0; m < d; ++m)
          {
            Real value = node[m][at[m]];
            for (const size_t l : all[m].denominator)
            {
              // Find the direction that owns coordinate l and shrink by it.
              for (size_t q = 0; q < d; ++q)
                if (all[q].numerator == l)
                  value *= 1 - node[q][at[q]];
            }
            x[static_cast<Eigen::Index>(all[m].numerator)] = value;
          }
          points.push_back(std::move(x));
          weights.push_back(w);
        }
      }

      /**
       * @brief Jacobi polynomial @f$ P_n^{(\alpha,\beta)}(x) @f$ by
       * recurrence.
       *
       * Both parameters are carried even though the basis only ever evaluates
       * @f$ \beta = 0 @f$: the derivative identity raises *both*, so a
       * routine fixed at @f$ \beta = 0 @f$ cannot express its own derivative
       * and silently returns the wrong polynomial instead.
       */
      static Real jacobi(size_t n, Real alpha, Real beta, Real x)
      {
        if (n == 0)
          return 1;
        Real previous = 1;
        Real current = Real(0.5) * (alpha - beta + (alpha + beta + 2) * x);
        for (size_t k = 2; k <= n; ++k)
        {
          const Real kk = static_cast<Real>(k);
          const Real a1 = 2 * kk * (kk + alpha + beta) * (2 * kk + alpha + beta - 2);
          const Real a2 = (2 * kk + alpha + beta - 1) *
            ((2 * kk + alpha + beta) * (2 * kk + alpha + beta - 2) * x + alpha * alpha -
              beta * beta);
          const Real a3 =
            2 * (kk + alpha - 1) * (kk + beta - 1) * (2 * kk + alpha + beta);
          const Real next = (a2 * current - a3 * previous) / a1;
          previous = current;
          current = next;
        }
        return current;
      }

      /**
       * @brief Derivative of @ref jacobi,
       * @f$ \tfrac{1}{2}(n + \alpha + \beta + 1) P^{(\alpha+1,\beta+1)}_{n-1} @f$.
       */
      static Real dJacobi(size_t n, Real alpha, Real beta, Real x)
      {
        if (n == 0)
          return 0;
        return Real(0.5) * (static_cast<Real>(n) + alpha + beta + 1) *
          jacobi(n - 1, alpha + 1, beta + 1, x);
      }

      /// @brief Integer power.
      static Real ipow(Real x, size_t e)
      {
        Real out = 1;
        for (size_t i = 0; i < e; ++i)
          out *= x;
        return out;
      }

    private:
      /// @brief Sum of the indices raising this direction's exponent.
      static size_t power(const Axis& axis, const std::vector<size_t>& index)
      {
        size_t out = 0;
        for (const size_t l : axis.contributors)
          out += index[l];
        return out;
      }

      /// @brief Jacobi exponent of this direction.
      static Real exponent(const Axis& axis, const std::vector<size_t>& index)
      {
        return 2 * static_cast<Real>(power(axis, index)) + static_cast<Real>(axis.offset);
      }
  };
}

#endif
