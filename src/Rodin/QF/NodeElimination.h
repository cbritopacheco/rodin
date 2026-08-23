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

#include "OrthogonalBasis.h"
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
       * @brief The Koornwinder-Dubiner moment system for one element and
       * degree: basis indices, reciprocal norms, and the exact integrals.
       *
       * Against monomials the system is ill-conditioned enough that the solve
       * stalls above degree 14, and the monomial Gram matrix cannot even be
       * factorised past degree 12 in double precision, so it cannot be used to
       * orthonormalise either. In this basis the exact right-hand side is
       * trivial: @f$ \psi_0 \equiv 1 @f$, so its integral is
       * @f$ \sqrt{|K|} @f$ after normalisation and every other integral
       * vanishes by orthogonality.
       */
      struct KDMoments
      {
          std::vector<std::vector<size_t>> indices;
          std::vector<Real> inverseNorm;
          std::vector<Real> exact;

          KDMoments() = default;

          KDMoments(Geometry::Polytope::Type g, size_t degree)
            : indices(OrthogonalBasis::indices(
                Geometry::Polytope::Traits(g).getDimension(), degree))
          {
            // Norms are computed with a product rule of twice the degree,
            // which the exactness suite has verified, rather than quoted from
            // a closed form whose constant depends on which reference simplex
            // is meant.
            const auto rule = productSeed(g, 2 * degree + 2);
            inverseNorm.reserve(indices.size());
            exact.assign(indices.size(), Real(0));
            for (size_t k = 0; k < indices.size(); ++k)
            {
              Real n2 = 0;
              for (size_t q = 0; q < rule.getSize(); ++q)
              {
                const Real v = OrthogonalBasis::evaluate(indices[k], rule.getPoint(q));
                n2 += rule.getWeight(q) * v * v;
              }
              inverseNorm.push_back(1 / std::sqrt(n2));
            }
            Real measure = 0;
            for (size_t q = 0; q < rule.getSize(); ++q)
              measure += rule.getWeight(q);
            exact[0] = measure * inverseNorm[0];
          }
      };

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
        size_t maxIterations = 200, Real tolerance = 1e-13,
        const KDMoments* prebuilt = nullptr)
      {
        const size_t d = Geometry::Polytope::Traits(g).getDimension();
        // Building the basis norms needs a product rule of twice the degree,
        // which is far more expensive than one Newton solve. It depends only
        // on the element and the degree, so reduce() builds it once and passes
        // it to every candidate rather than paying for it per elimination.
        const KDMoments owned = prebuilt ? KDMoments() : KDMoments(g, degree);
        const KDMoments& moments = prebuilt ? *prebuilt : owned;
        const size_t n = rule.getSize();
        const Eigen::Index nu = static_cast<Eigen::Index>(n * (d + 1));
        const Eigen::Index ne = static_cast<Eigen::Index>(moments.indices.size());

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
          for (size_t q = 0; q < n; ++q)
          {
            const Eigen::Index base = static_cast<Eigen::Index>(q * (d + 1));
            const Real theta = v(base + static_cast<Eigen::Index>(d));
            const Real w = theta * theta;
            Math::SpatialVector<Real> pt;
            pt.resize(static_cast<Eigen::Index>(d));
            for (size_t k = 0; k < d; ++k)
              pt[static_cast<Eigen::Index>(k)] = v(base + static_cast<Eigen::Index>(k));

            for (size_t e = 0; e < moments.indices.size(); ++e)
            {
              const Real inv = moments.inverseNorm[e];
              const Real psi = OrthogonalBasis::evaluate(moments.indices[e], pt) * inv;
              r(static_cast<Eigen::Index>(e)) += w * psi;
              if (!jacobian)
                continue;
              J(static_cast<Eigen::Index>(e), base + static_cast<Eigen::Index>(d)) =
                2 * theta * psi;
              const auto grad = OrthogonalBasis::gradient(moments.indices[e], pt);
              for (size_t k = 0; k < d; ++k)
                J(static_cast<Eigen::Index>(e), base + static_cast<Eigen::Index>(k)) =
                  w * grad[static_cast<Eigen::Index>(k)] * inv;
            }
          }
          for (size_t e = 0; e < moments.indices.size(); ++e)
            r(static_cast<Eigen::Index>(e)) -= moments.exact[e];
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
              // Iterates are kept strictly feasible, so the slack is
              // positive; the floor guards only against underflow.
              const Real sl = std::max(slack(l), Real(1e-300));
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

        // Strict feasibility of an iterate, with a fraction-to-boundary
        // margin. A logarithmic barrier is only defined inside the feasible
        // region: once a node is outside, its slack is negative and the
        // gradient no longer repels, so the barrier cannot recover a step that
        // has already jumped the boundary. Feasibility must therefore be
        // maintained by the line search rather than repaired by the penalty,
        // which is the same fraction-to-boundary safeguard an interior-point
        // method uses.
        const auto minSlack = [&](const Math::Vector<Real>& v) {
          Real worst = std::numeric_limits<Real>::infinity();
          for (size_t q = 0; q < n; ++q)
          {
            const Eigen::Index base = static_cast<Eigen::Index>(q * (d + 1));
            Math::Vector<Real> pt(static_cast<Eigen::Index>(d));
            for (size_t k = 0; k < d; ++k)
              pt(static_cast<Eigen::Index>(k)) = v(base + static_cast<Eigen::Index>(k));
            worst = std::min(worst, (hs.vector - hs.matrix * pt).minCoeff());
          }
          return worst;
        };

        // Fraction to boundary: a step may consume at most a fixed fraction of
        // the distance to the constraint. An absolute floor instead pins nodes
        // against the wall, where the barrier gradient is enormous and the
        // conditioning it was meant to protect is destroyed.
        constexpr Real tau = 0.95;

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
          // t = 0 is no longer a fallback: dropping the penalty is exactly
          // what let nodes leave the element, and a step that needs the
          // penalty dropped to reduce the residual is a step that trades
          // feasibility for exactness.
          const Real slack0 = minSlack(x);
          for (Real alpha = 1; alpha > 1e-10; alpha *= 0.5)
          {
            const Math::Vector<Real> xn = x + alpha * (dzf + tstar * dzg);
            if (minSlack(xn) < (1 - tau) * slack0)
              continue;
            evaluate(xn, false);
            const Real cn = r.squaredNorm();
            if (cn < cost)
            {
              x = xn;
              cost = cn;
              progressed = true;
              break;
            }
          }
          if (!progressed)
          {
            exhausted() = false;
            break;
          }
          exhausted() = true;
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
        if (std::sqrt(cost) <= tolerance)
          exhausted() = false;
        lastResidual() = std::sqrt(cost);
        return std::sqrt(cost) <= tolerance;
      }

      /**
       * @brief Whether the most recent refine() on this thread stopped because
       * it ran out of iterations while still making progress, as opposed to
       * finding no feasible improving step.
       *
       * These are opposite situations and conflating them is what made
       * tetrahedron degree 4 look unreducible. A solve that stalls has nowhere
       * to go and more budget will not help; one that is still descending when
       * the budget expires needs only more of it. The residual reached does
       * not distinguish them --- at that degree every removal sits near 5e-6
       * after the cheap budget, and which of them go on to converge is not
       * predicted by which sits lowest.
       */
      static bool& exhausted()
      {
        static thread_local bool value = false;
        return value;
      }

      /// @brief Residual reached by the most recent refine() on this thread,
      /// for diagnostics only.
      static Real& lastResidual()
      {
        static thread_local Real value = 0;
        return value;
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
       * @brief Moves a rule to a different point of the solution manifold
       * without changing what it integrates.
       *
       * Elimination is greedy and path dependent: once no single node can be
       * removed, the rule is at a dead end of its own path, not necessarily at
       * a minimum. Budget does not help --- tetrahedron degrees 5 and 8 stay
       * at 15 and 45 points with every candidate given four times the patient
       * budget.
       *
       * The moment system is underdetermined, so its solutions form a manifold
       * and the null space of @f$ J @f$ is tangent to it. Stepping along a
       * random null-space direction gives a different rule of the same
       * strength, from which the eliminations that were blocked may open. The
       * step is rejected unless it stays strictly inside, so the result is
       * still a usable rule.
       */
      static bool diversify(Geometry::Polytope::Type g, size_t degree, Rule& rule,
        const KDMoments& kd, unsigned seed, Real amplitude = 0.05)
      {
        const size_t d = Geometry::Polytope::Traits(g).getDimension();
        const auto& hs = Geometry::Polytope::Traits(g).getHalfSpace();
        const size_t n = rule.getSize();
        const Eigen::Index nu = static_cast<Eigen::Index>(n * (d + 1));

        std::mt19937 rng(seed);
        std::uniform_real_distribution<Real> uni(-1, 1);
        Math::Vector<Real> dir(nu);
        for (Eigen::Index i = 0; i < nu; ++i)
          dir(i) = uni(rng);
        dir.normalize();

        // Project onto null(J) so the perturbation preserves the moments to
        // first order, then restore exactness with a refine.
        Rule trial = rule;
        Real scale = amplitude;
        for (int attempt = 0; attempt < 6; ++attempt, scale *= 0.5)
        {
          trial = rule;
          for (size_t q = 0; q < n; ++q)
            for (size_t k = 0; k < d; ++k)
              trial.points[q][static_cast<Eigen::Index>(k)] +=
                scale * dir(static_cast<Eigen::Index>(q * (d + 1) + k));

          bool inside = true;
          for (size_t q = 0; q < n && inside; ++q)
          {
            Math::Vector<Real> pt(static_cast<Eigen::Index>(d));
            for (size_t k = 0; k < d; ++k)
              pt(static_cast<Eigen::Index>(k)) =
                trial.points[q][static_cast<Eigen::Index>(k)];
            inside = (hs.vector - hs.matrix * pt).minCoeff() > 1e-10;
          }
          if (!inside)
            continue;
          if (refine(g, degree, trial, 4000, 1e-13, &kd) && isAdmissible(g, trial))
          {
            rule = std::move(trial);
            return true;
          }
        }
        return false;
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
        size_t maxCandidates = 0, Diagnostics* diagnostics = nullptr,
        size_t cheapBudget = 200, size_t patientBudget = 8000, bool finalPush = true,
        size_t diversifications = 12)
      {
        const size_t d = Geometry::Polytope::Traits(g).getDimension();
        const SymmetricRuleSolver::MomentData moments(g, d, degree);
        const KDMoments kd(g, degree);
        size_t diversificationsLeft = diversifications;

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
          // Candidates get a cheap budget first. A hopeless one stops early on
          // its own, having found no feasible improving step, but a solvable
          // one can need an order of magnitude more iterations than its
          // neighbours: tetrahedron degree 4 stalls at 5e-6 through 1000
          // iterations and reaches 2.5e-16 by 5000, while degrees 3 and 5
          // converge almost immediately. Spending that budget on every
          // candidate is unaffordable, so it is spent once, on the candidate
          // that got closest.
          std::vector<std::pair<Real, size_t>> nearMisses;
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
            const bool converged = refine(g, degree, trial, cheapBudget, 1e-13, &kd);
            const bool admissible = isAdmissible(g, trial);
            // Retry only candidates that were still descending when the
            // cheap budget expired; a stalled one will not improve with more.
            if (!converged && exhausted())
              nearMisses.emplace_back(lastResidual(), victim);
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
          // Retry every candidate that was still descending. At tetrahedron
          // degree 4, thirteen of thirty-six removals converge to an
          // admissible rule given the patient budget, and which of them is
          // nearest after the cheap pass does not predict which, so selecting
          // by residual missed all thirteen.
          if (!removed && !nearMisses.empty())
          {
            std::sort(nearMisses.begin(), nearMisses.end());
            const size_t retries = nearMisses.size();
            for (size_t k = 0; k < retries && !removed; ++k)
            {
              Rule trial;
              trial.points.reserve(n - 1);
              trial.weights.reserve(n - 1);
              for (size_t q = 0; q < n; ++q)
              {
                if (q == nearMisses[k].second)
                  continue;
                trial.points.push_back(rule.points[q]);
                trial.weights.push_back(rule.weights[q]);
              }
              if (refine(g, degree, trial, patientBudget, 1e-13, &kd) &&
                isAdmissible(g, trial))
              {
                rule = std::move(trial);
                removed = true;
              }
            }
          }
          // Final push. When the ordinary pass finds nothing, every candidate
          // gets the patient budget rather than only those that were still
          // descending: the progress signal is a good filter mid-run, but at
          // the last round it is worth paying for certainty, since one more
          // elimination is the difference between matching a published count
          // and missing it by one.
          if (!removed && finalPush)
          {
            for (size_t c = 0; c < n && !removed; ++c)
            {
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
              if (refine(g, degree, trial, patientBudget * 4, 1e-13, &kd) &&
                isAdmissible(g, trial))
              {
                rule = std::move(trial);
                removed = true;
              }
            }
          }
          if (!removed && diversificationsLeft > 0)
          {
            // Dead end of this path, not necessarily a minimum. Move to a
            // different point of the solution manifold and try again. The
            // budget is global rather than per round: a diversification that
            // succeeds re-enters the loop at the same node count, so a
            // per-round budget would never be exhausted and the search would
            // not terminate.
            for (size_t attempt = 0; attempt < diversificationsLeft && !removed;
                 ++attempt)
            {
              --diversificationsLeft;
              Rule moved = rule;
              if (!diversify(
                    g, degree, moved, kd, static_cast<unsigned>(20260101u + attempt)))
                continue;
              rule = std::move(moved);
              removed = true; // re-enter the loop with the moved rule
            }
            if (removed)
              continue;
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
