/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_QF_NODEELIMINATION_H
#define RODIN_QF_NODEELIMINATION_H

/**
 * @file
 * @brief Defines NodeElimination, which reduces a quadrature rule to a
 * near-minimal one by removing nodes.
 */

#include <Eigen/QR>

#include "SymmetricRuleSolver.h"

namespace Rodin::QF
{
  /**
   * @ingroup RodinQuadrature
   * @brief Reduces an exact rule to a near-minimal one by eliminating nodes.
   *
   * This is the construction of Xiao and Gimbutas @cite xiao2010numerical.
   * Begin from a rule that is already exact but has more points than
   * necessary; remove one node; the reduced rule no longer satisfies the
   * moment equations, but it is *close* to one that does, so it serves as the
   * initial guess for a Gauss--Newton solve on the remaining nodes. Repeat
   * until no node can be removed.
   *
   * The distinction from SymmetricRuleSolver matters. That class searches for
   * a rule of prescribed symmetry by guessing an orbit structure and
   * restarting Levenberg--Marquardt from random points; every solve is a
   * global search. Here every solve is warm-started beside a known solution,
   * so it converges in a few iterations, progress is monotone --- each success
   * removes exactly one point --- and the structure of the rule is discovered
   * rather than prescribed.
   *
   * The price is symmetry: eliminating a node breaks it. Xiao--Gimbutas rules
   * are consequently not symmetric, whereas the symmetric rules of Witherden
   * and Vincent @cite witherden2015identification are what
   * SymmetricRuleSolver produces. The two families need different generators,
   * and this is the one for the former.
   */
  class NodeElimination
  {
    public:
      /**
       * @brief Why an elimination round failed.
       *
       * reduce() reports only whether a candidate was accepted, and the two
       * ways of being rejected call for opposite responses: a solve that did
       * not converge wants more iterations or better conditioning, whereas one
       * that converged onto an inadmissible rule wants a stronger penalty.
       * Told apart, they are actionable; conflated, the second reads as the
       * first and sends you optimising the wrong thing.
       */
      struct Diagnostics
      {
          size_t candidatesTried = 0;    ///< Nodes attempted in the final round.
          size_t notConverged = 0;       ///< Rejected: residual above tolerance.
          size_t notAdmissible = 0;      ///< Rejected: converged but outside.
          Real bestResidual = std::numeric_limits<Real>::infinity();
          size_t worstNegativeWeights = 0; ///< In the best converged candidate.
          size_t worstExteriorPoints = 0;  ///< In the best converged candidate.
      };

      /// @brief A rule as free points and weights, with no symmetry imposed.
      struct Rule
      {
          std::vector<Math::SpatialVector<Real>> points;
          std::vector<Real> weights;

          size_t getSize() const
          {
            return weights.size();
          }
          Real getWeight(size_t i) const
          {
            return weights[i];
          }
          const Math::SpatialVector<Real>& getPoint(size_t i) const
          {
            return points[i];
          }
      };

      /**
       * @brief A collapsed Gauss product rule on the reference simplex.
       *
       * Used as the seed. The Duffy collapse contributes a Jacobian factor
       * that raises the degree in the collapsed directions, so each direction
       * carries the point count its own degree requires --- the same
       * bookkeeping that the wedge and pyramid dispatch got wrong.
       *
       * The rule is positive and exact by construction, which is what a seed
       * must be: elimination preserves exactness and positivity, it does not
       * establish them.
       */
      /// @brief Integer power.
      static Real ipow(Real x, size_t e)
      {
        Real r = 1;
        while (e--)
          r *= x;
        return r;
      }

      static Rule productSeed(Geometry::Polytope::Type g, size_t degree)
      {
        const size_t d = Geometry::Polytope::Traits(g).getDimension();
        // Gauss-Legendre on [0,1] with n points is exact to degree 2n - 1.
        const auto gauss = [](size_t n) {
          std::vector<Real> x(n), w(n);
          for (size_t i = 0; i < n; ++i)
          {
            // Newton on the Legendre polynomial, mapped to [0, 1].
            Real t = std::cos(
              M_PI * (static_cast<Real>(i) + 0.75) / (static_cast<Real>(n) + 0.5));
            for (int it = 0; it < 100; ++it)
            {
              Real p0 = 1, p1 = 0;
              for (size_t k = 0; k < n; ++k)
              {
                const Real p2 = p1;
                p1 = p0;
                p0 =
                  ((2 * static_cast<Real>(k) + 1) * t * p1 - static_cast<Real>(k) * p2) /
                  static_cast<Real>(k + 1);
              }
              const Real dp = static_cast<Real>(n) * (t * p0 - p1) / (t * t - 1);
              const Real dt = -p0 / dp;
              t += dt;
              if (std::abs(dt) < 1e-16)
                break;
            }
            Real p0 = 1, p1 = 0;
            for (size_t k = 0; k < n; ++k)
            {
              const Real p2 = p1;
              p1 = p0;
              p0 = ((2 * static_cast<Real>(k) + 1) * t * p1 - static_cast<Real>(k) * p2) /
                static_cast<Real>(k + 1);
            }
            const Real dp = static_cast<Real>(n) * (t * p0 - p1) / (t * t - 1);
            x[i] = 0.5 * (1 - t);
            w[i] = 1.0 / ((1 - t * t) * dp * dp);
          }
          return std::make_pair(x, w);
        };

        Rule rule;
        // Direction k of the collapse needs degree p + (d - 1 - k).
        std::vector<std::vector<Real>> nodes(d), wts(d);
        for (size_t k = 0; k < d; ++k)
        {
          const size_t need = degree + (d - 1 - k);
          const size_t n = (need + 2) / 2;
          auto gw = gauss(n);
          nodes[k] = std::move(gw.first);
          wts[k] = std::move(gw.second);
        }

        std::vector<size_t> idx(d, 0);
        for (;;)
        {
          // Duffy collapse: x_0 = u_0, x_k = u_k * prod_{j<k} (1 - u_j).
          // Its Jacobian is prod_{k < d-1} (1 - u_k)^{d-1-k}.
          Math::SpatialVector<Real> p;
          p.resize(static_cast<Eigen::Index>(d));
          Real w = 1, remaining = 1;
          for (size_t k = 0; k < d; ++k)
          {
            const Real u = nodes[k][idx[k]];
            p[static_cast<Eigen::Index>(k)] = remaining * u;
            w *= wts[k][idx[k]];
            remaining *= (1 - u);
          }
          for (size_t k = 0; k + 1 < d; ++k)
            w *= ipow(1 - nodes[k][idx[k]], d - 1 - k);
          rule.points.push_back(std::move(p));
          rule.weights.push_back(w);

          size_t k = 0;
          for (; k < d; ++k)
          {
            if (++idx[k] < nodes[k].size())
              break;
            idx[k] = 0;
          }
          if (k == d)
            break;
        }
        return rule;
      }

      /**
       * @brief Refines @p rule onto the moment equations by Gauss--Newton.
       *
       * Unknowns are the node coordinates together with @f$ \theta_q @f$,
       * where @f$ w_q = \theta_q^2 @f$, so weight positivity is a property of
       * the parameterisation. Points are free variables here rather than
       * functions of orbit parameters, which makes the Jacobian explicit:
       * @f[
       *   \frac{\partial r_\alpha}{\partial x_{q,k}}
       *     = s_\alpha\, w_q\, \alpha_k\, x_{q,k}^{\alpha_k - 1}
       *       \prod_{j \ne k} x_{q,j}^{\alpha_j},
       *   \qquad
       *   \frac{\partial r_\alpha}{\partial \theta_q}
       *     = 2 s_\alpha \theta_q \prod_k x_{q,k}^{\alpha_k}.
       * @f]
       * The finite-difference Jacobian SymmetricRuleSolver uses costs
       * @f$ O(\text{unknowns}) @f$ residual evaluations per iteration; this
       * costs one, which is what makes elimination affordable at high degree.
       */
      static bool refine(Geometry::Polytope::Type g, size_t degree, Rule& rule,
        size_t maxIterations = 200, Real tolerance = 1e-13)
      {
        const size_t d = Geometry::Polytope::Traits(g).getDimension();
        const SymmetricRuleSolver::MomentData moments(g, d, degree);
        const size_t n = rule.getSize();
        const Eigen::Index nu = static_cast<Eigen::Index>(n * (d + 1));
        const Eigen::Index ne = static_cast<Eigen::Index>(moments.alphas.size());

        Math::Vector<Real> x(nu);
        for (size_t q = 0; q < n; ++q)
        {
          for (size_t k = 0; k < d; ++k)
            x(static_cast<Eigen::Index>(q * (d + 1) + k)) =
              rule.points[q][static_cast<Eigen::Index>(k)];
          x(static_cast<Eigen::Index>(q * (d + 1) + d)) =
            std::sqrt(std::max(rule.weights[q], Real(0)));
        }

        Math::Vector<Real> r(ne);
        Math::Matrix<Real> J(ne, nu);
        const auto evaluate = [&](const Math::Vector<Real>& v, bool jacobian) {
          r.setZero();
          if (jacobian)
            J.setZero();
          for (size_t e = 0; e < moments.alphas.size(); ++e)
          {
            const auto& alpha = moments.alphas[e];
            const Real sc = moments.scale[e];
            Real sum = 0;
            for (size_t q = 0; q < n; ++q)
            {
              const Eigen::Index base = static_cast<Eigen::Index>(q * (d + 1));
              const Real theta = v(base + static_cast<Eigen::Index>(d));
              const Real w = theta * theta;
              Real mono = 1;
              for (size_t k = 0; k < d; ++k)
                mono *= ipow(v(base + static_cast<Eigen::Index>(k)), alpha[k]);
              sum += w * mono;
              if (!jacobian)
                continue;
              J(static_cast<Eigen::Index>(e), base + static_cast<Eigen::Index>(d)) =
                sc * 2 * theta * mono;
              for (size_t k = 0; k < d; ++k)
              {
                if (alpha[k] == 0)
                  continue;
                Real partial = 1;
                for (size_t j = 0; j < d; ++j)
                  partial *= ipow(v(base + static_cast<Eigen::Index>(j)),
                    (j == k) ? alpha[j] - 1 : alpha[j]);
                J(static_cast<Eigen::Index>(e), base + static_cast<Eigen::Index>(k)) =
                  sc * w * static_cast<Real>(alpha[k]) * partial;
              }
            }
            r(static_cast<Eigen::Index>(e)) = (sum - moments.exact[e]) * sc;
          }
        };

        // Penalised least-squares Newton, after Slobodkins and Tausch,
        // "A node elimination algorithm for cubature of high-dimensional
        // polytopes", arXiv:2207.10737.
        //
        // The moment system is underdetermined --- 231 equations against 363
        // unknowns at triangle degree 20 --- and that surplus is not a
        // nuisance but the resource that keeps nodes admissible. Writing
        // J = LQ for the economy LQ factorisation, the step splits as
        //
        //   dz_f = -Q^T L^{-1} f      (minimum-norm Newton step)
        //   dz_g = -(I - Q^T Q) g     (projection of -g onto null(J))
        //
        // with g the gradient of a logarithmic barrier on the constraints. The
        // first term restores exactness; the second slides along the solution
        // manifold away from the boundary without disturbing it. Screening
        // admissibility after the solve, as this did before, throws away
        // rules that this step would have kept.
        //
        // The barrier uses the half-space description Ax <= b published by
        // Geometry::Polytope::Traits, so the domain the solver respects is by
        // construction the reference element the rest of Rodin uses.
        const Geometry::Polytope::Traits traits(g);
        const auto& hs = traits.getHalfSpace();

        const auto barrierGradient = [&](const Math::Vector<Real>& v) {
          Math::Vector<Real> grad(nu);
          grad.setZero();
          for (size_t q = 0; q < n; ++q)
          {
            const Eigen::Index base = static_cast<Eigen::Index>(q * (d + 1));
            Math::Vector<Real> pt(static_cast<Eigen::Index>(d));
            for (size_t k = 0; k < d; ++k)
              pt(static_cast<Eigen::Index>(k)) = v(base + static_cast<Eigen::Index>(k));
            const Math::Vector<Real> slack = hs.vector - hs.matrix * pt;
            for (Eigen::Index l = 0; l < slack.size(); ++l)
            {
              const Real sl = std::max(slack(l), Real(1e-14));
              for (size_t k = 0; k < d; ++k)
                grad(base + static_cast<Eigen::Index>(k)) +=
                  hs.matrix(l, static_cast<Eigen::Index>(k)) / sl;
            }
            const Real theta = v(base + static_cast<Eigen::Index>(d));
            grad(base + static_cast<Eigen::Index>(d)) =
              -2 / (std::abs(theta) > 1e-14 ? theta : Real(1e-14));
          }
          return grad;
        };

        evaluate(x, true);
        Real cost = r.squaredNorm();
        for (size_t it = 0; it < maxIterations && std::sqrt(cost) > tolerance; ++it)
        {
          // J = LQ via the QR of J^T: J^T = Qt Rt gives L = Rt^T, Q = Qt^T.
          const Math::Matrix<Real> Jt = J.transpose();
          Eigen::HouseholderQR<Math::Matrix<Real>> qr(Jt);
          const Math::Matrix<Real> Qt =
            qr.householderQ() * Math::Matrix<Real>::Identity(nu, ne);
          const Math::Matrix<Real> Rt =
            qr.matrixQR().topLeftCorner(ne, ne).template triangularView<Eigen::Upper>();

          const Math::Vector<Real> y =
            Rt.transpose().template triangularView<Eigen::Lower>().solve(r);
          const Math::Vector<Real> dzf = -(Qt * y);

          const Math::Vector<Real> gb = barrierGradient(x);
          const Math::Vector<Real> dzg = -(gb - Qt * (Qt.transpose() * gb));

          // t is chosen to push the next iterate as far inside as possible.
          // max_k (beta_k + t alpha_k) is convex in t, so a ternary search
          // over a bounded interval finds its minimiser without enumerating
          // the intersections.
          const auto worstConstraint = [&](Real t) {
            Real worst = -std::numeric_limits<Real>::infinity();
            const Math::Vector<Real> xn = x + dzf + t * dzg;
            for (size_t q = 0; q < n; ++q)
            {
              const Eigen::Index base = static_cast<Eigen::Index>(q * (d + 1));
              Math::Vector<Real> pt(static_cast<Eigen::Index>(d));
              for (size_t k = 0; k < d; ++k)
                pt(static_cast<Eigen::Index>(k)) =
                  xn(base + static_cast<Eigen::Index>(k));
              const Math::Vector<Real> res = hs.matrix * pt - hs.vector;
              for (Eigen::Index l = 0; l < res.size(); ++l)
                worst = std::max(worst, res(l));
              const Real theta = xn(base + static_cast<Eigen::Index>(d));
              worst = std::max(worst, -theta * theta);
            }
            return worst;
          };

          Real lo = 0, hi = 1;
          for (int k = 0; k < 60; ++k)
          {
            const Real m1 = lo + (hi - lo) / 3;
            const Real m2 = hi - (hi - lo) / 3;
            if (worstConstraint(m1) < worstConstraint(m2))
              hi = m2;
            else
              lo = m1;
          }
          const Real tstar = 0.5 * (lo + hi);

          bool progressed = false;
          for (const Real t : {tstar, Real(0)})
          {
            Real alpha = 1;
            for (int back = 0; back < 24; ++back)
            {
              const Math::Vector<Real> xn = x + alpha * (dzf + t * dzg);
              evaluate(xn, false);
              const Real cn = r.squaredNorm();
              if (cn < cost)
              {
                x = xn;
                cost = cn;
                progressed = true;
                break;
              }
              alpha *= 0.5;
            }
            if (progressed)
              break;
          }
          if (!progressed)
            break;
          evaluate(x, true);
        }

        for (size_t q = 0; q < n; ++q)
        {
          const Eigen::Index base = static_cast<Eigen::Index>(q * (d + 1));
          for (size_t k = 0; k < d; ++k)
            rule.points[q][static_cast<Eigen::Index>(k)] =
              x(base + static_cast<Eigen::Index>(k));
          const Real theta = x(base + static_cast<Eigen::Index>(d));
          rule.weights[q] = theta * theta;
        }
        return std::sqrt(cost) <= tolerance;
      }

      /// @brief Whether every weight is positive and every point interior.
      static bool isAdmissible(
        Geometry::Polytope::Type g, const Rule& rule, Real tol = 1e-12)
      {
        const Geometry::Polytope::Traits traits(g);
        const auto& hs = traits.getHalfSpace();
        const size_t d = traits.getDimension();
        for (size_t q = 0; q < rule.getSize(); ++q)
        {
          if (!(rule.weights[q] > tol))
            return false;
          Math::Vector<Real> p(static_cast<Eigen::Index>(d));
          for (size_t k = 0; k < d; ++k)
            p(static_cast<Eigen::Index>(k)) =
              rule.points[q][static_cast<Eigen::Index>(k)];
          const Math::Vector<Real> res = hs.matrix * p - hs.vector;
          for (Eigen::Index i = 0; i < res.size(); ++i)
            if (res(i) > -tol)
              return false;
        }
        return true;
      }

      /**
       * @brief Eliminates nodes until none can be removed.
       *
       * Candidates are tried in increasing order of significance
       * @f$ w_q \lVert v_q \rVert @f$, where @f$ v_q @f$ is the node's column
       * of scaled monomial values: the node contributing least to the moments
       * is the one most likely to be removable, which is the selection Xiao
       * and Gimbutas use.
       */
      static Rule reduce(Geometry::Polytope::Type g, size_t degree, Rule rule,
        size_t maxCandidates = 0, Diagnostics* diagnostics = nullptr)
      {
        const size_t d = Geometry::Polytope::Traits(g).getDimension();
        const SymmetricRuleSolver::MomentData moments(g, d, degree);

        for (;;)
        {
          const size_t n = rule.getSize();
          if (n <= 1)
            return rule;

          std::vector<std::pair<Real, size_t>> significance;
          significance.reserve(n);
          for (size_t q = 0; q < n; ++q)
          {
            Real norm = 0;
            for (size_t e = 0; e < moments.alphas.size(); ++e)
            {
              Real mono = 1;
              for (size_t k = 0; k < d; ++k)
                mono *= ipow(
                  rule.points[q][static_cast<Eigen::Index>(k)], moments.alphas[e][k]);
              const Real v = mono * moments.scale[e];
              norm += v * v;
            }
            significance.emplace_back(rule.weights[q] * std::sqrt(norm), q);
          }
          std::sort(significance.begin(), significance.end());

          const size_t budget = maxCandidates ? std::min(maxCandidates, n) : n;
          bool removed = false;
          Diagnostics round;
          for (size_t c = 0; c < budget && !removed; ++c)
          {
            ++round.candidatesTried;
            const size_t victim = significance[c].second;
            Rule trial;
            trial.points.reserve(n - 1);
            trial.weights.reserve(n - 1);
            for (size_t q = 0; q < n; ++q)
            {
              if (q == victim)
                continue;
              trial.points.push_back(rule.points[q]);
              trial.weights.push_back(rule.weights[q]);
            }
            const bool converged = refine(g, degree, trial);
            const bool admissible = isAdmissible(g, trial);
            if (converged && admissible)
            {
              rule = std::move(trial);
              removed = true;
              break;
            }

            if (!converged)
            {
              ++round.notConverged;
            }
            else
            {
              ++round.notAdmissible;
              size_t neg = 0, out = 0;
              for (size_t q = 0; q < trial.getSize(); ++q)
              {
                if (!(trial.weights[q] > 0))
                  ++neg;
                Rule single;
                single.points.push_back(trial.points[q]);
                single.weights.push_back(1);
                if (!isAdmissible(g, single, 1e-12))
                  ++out;
              }
              if (neg + out < round.worstNegativeWeights + round.worstExteriorPoints ||
                round.notAdmissible == 1)
              {
                round.worstNegativeWeights = neg;
                round.worstExteriorPoints = out;
              }
            }
          }
          if (!removed)
          {
            if (diagnostics)
              *diagnostics = round;
            return rule;
          }
        }
      }
  };
}

#endif
