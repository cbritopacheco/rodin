/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_QF_SYMMETRICRULESOLVER_H
#define RODIN_QF_SYMMETRICRULESOLVER_H

/**
 * @file
 * @brief Defines SymmetricRuleSolver, which solves the moment equations of a
 * fully symmetric simplex quadrature rule.
 */

#include <random>
#include <functional>
#include <numeric>

#include "SymmetricOrbit.h"

namespace Rodin::QF
{
  /**
   * @ingroup RodinQuadrature
   * @brief Solves the moment equations of a fully symmetric simplex rule.
   *
   * A rule of strength @f$ p @f$ is a set of orbits whose points and weights
   * integrate every polynomial of degree @f$ \le p @f$ exactly. Writing the
   * rule as symmetry orbits (SymmetricOrbit) turns that requirement into a
   * small nonlinear system: the unknowns are the free barycentric parameters
   * of each orbit together with one weight per orbit, and the equations are
   * the monomial moments of the reference simplex.
   *
   * This is the construction of Xiao and Gimbutas @cite xiao2010numerical and
   * of Witherden and Vincent @cite witherden2015identification. Neither
   * publishes a closed form; both search for solutions numerically and retain
   * those whose weights are positive and whose points lie inside the element.
   * The same is done here.
   *
   * ## Orbit configuration
   *
   * An orbit class is described by its *multiplicity pattern*: the partition
   * of @f$ d+1 @f$ giving how many barycentric coordinates share each value.
   * On the triangle @f$ \{3\} @f$, @f$ \{2,1\} @f$ and @f$ \{1,1,1\} @f$ are
   * the classes written @f$ S_3 @f$, @f$ S_{21}(a) @f$ and
   * @f$ S_{111}(a,b) @f$; on the tetrahedron @f$ \{4\} @f$, @f$ \{3,1\} @f$,
   * @f$ \{2,2\} @f$, @f$ \{2,1,1\} @f$ and @f$ \{1,1,1,1\} @f$ are
   * @f$ S_4 @f$ through @f$ S_{1111} @f$. A pattern with @f$ k @f$ parts
   * carries @f$ k-1 @f$ free parameters, the last being fixed by
   * @f$ \sum_i \lambda_i = 1 @f$, plus one weight.
   *
   * ## Method
   *
   * The residual is the normalised moment defect
   * @f[
   *   r_\alpha = \frac{1}{\max(|I_\alpha|, \varepsilon)}
   *     \left( \sum_q w_q x_q^\alpha - I_\alpha \right),
   *   \qquad |\alpha| \le p,
   * @f]
   * over the exact simplex moments
   * @f$ I_\alpha = \alpha! / (|\alpha| + d)! @f$. Normalising each equation by
   * its own moment keeps the rows of comparable size; the raw moments span
   * many orders of magnitude and the unnormalised system is badly scaled well
   * before the degrees of interest.
   *
   * The system is solved by Levenberg--Marquardt with a central-difference
   * Jacobian, from pseudo-random starting points. The generator is seeded
   * deterministically, so a given (geometry, degree, configuration) always
   * yields the same rule: coefficients that are regenerated must not drift.
   *
   * @see SymmetricOrbit
   */
  class SymmetricRuleSolver
  {
    public:
      /// @brief Partition of @f$ d+1 @f$ describing one orbit class.
      using Pattern = std::vector<size_t>;

      /// @brief The orbit classes making up a candidate rule.
      using Configuration = std::vector<Pattern>;

      /// @brief Outcome of a solve.
      struct Result
      {
        bool converged = false;             ///< Residual met the tolerance.
        bool admissible = false;            ///< Positive weights, interior points.
        std::vector<SymmetricOrbit> orbits; ///< The rule, as orbits.
        Real residual = std::numeric_limits<Real>::infinity(); ///< Final norm.
        size_t restarts = 0;                ///< Starting points consumed.
      };

      /// @brief Exact moment @f$ \int_K x^\alpha @f$ on the reference
      /// @f$ d @f$-simplex, @f$ \alpha! / (|\alpha| + d)! @f$.
      static Real simplexMoment(const std::vector<size_t>& alpha)
      {
        const auto fact = [](size_t n)
        {
          Real r = 1;
          for (size_t i = 2; i <= n; ++i) r *= static_cast<Real>(i);
          return r;
        };
        Real num = 1;
        size_t total = 0;
        for (const auto a : alpha)
        {
          num *= fact(a);
          total += a;
        }
        return num / fact(total + alpha.size());
      }

      /// @brief Every exponent tuple with @f$ |\alpha| \le p @f$ in
      /// @f$ d @f$ variables.
      static std::vector<std::vector<size_t>> monomials(size_t d, size_t p)
      {
        std::vector<std::vector<size_t>> out;
        std::vector<size_t> alpha(d, 0);
        const std::function<void(size_t, size_t)> rec =
          [&](size_t k, size_t budget)
          {
            if (k == d)
            {
              out.push_back(alpha);
              return;
            }
            for (size_t a = 0; a <= budget; ++a)
            {
              alpha[k] = a;
              rec(k + 1, budget - a);
            }
            alpha[k] = 0;
          };
        rec(0, p);
        return out;
      }

      /// @brief Number of unknowns implied by @p config: one weight per orbit
      /// plus one free barycentric parameter per part beyond the first.
      static size_t unknownCount(const Configuration& config)
      {
        size_t n = 0;
        for (const auto& pattern : config)
          n += pattern.size();  // (parts - 1) parameters + 1 weight
        return n;
      }

      /**
       * @brief Builds the orbits implied by a parameter vector.
       *
       * Layout per orbit: @f$ k-1 @f$ barycentric parameters followed by the
       * weight, where @f$ k @f$ is the number of parts. The value of the last
       * part is fixed by @f$ \sum_i \lambda_i = 1 @f$.
       */
      static std::vector<SymmetricOrbit> toOrbits(
        const Configuration& config, const Math::Vector<Real>& x)
      {
        std::vector<SymmetricOrbit> orbits;
        Eigen::Index cursor = 0;
        for (const auto& pattern : config)
        {
          const size_t k = pattern.size();
          std::vector<Real> values(k);
          Real used = 0;
          size_t assigned = 0;
          for (size_t i = 0; i + 1 < k; ++i)
          {
            values[i] = x(cursor++);
            used += static_cast<Real>(pattern[i]) * values[i];
            assigned += pattern[i];
          }
          const size_t remaining = pattern.back();
          values[k - 1] = (Real(1) - used) / static_cast<Real>(remaining);
          const Real weight = x(cursor++);

          SymmetricOrbit::Barycentric barycentric;
          for (size_t i = 0; i < k; ++i)
            for (size_t m = 0; m < pattern[i]; ++m)
              barycentric.push_back(values[i]);
          orbits.emplace_back(std::move(barycentric), weight);
        }
        return orbits;
      }

      /// @brief Normalised moment residual of the rule described by @p x.
      static Math::Vector<Real> residual(Geometry::Polytope::Type g,
        size_t degree, const Configuration& config, const Math::Vector<Real>& x)
      {
        const Geometry::Polytope::Traits traits(g);
        const size_t d = traits.getDimension();
        const auto alphas = monomials(d, degree);
        const auto orbits = toOrbits(config, x);

        Math::Vector<Real> r(static_cast<Eigen::Index>(alphas.size()));
        for (size_t e = 0; e < alphas.size(); ++e)
        {
          const auto& alpha = alphas[e];
          Real sum = 0;
          for (const auto& orbit : orbits)
          {
            for (const auto& b : orbit.expand())
            {
              const auto p = SymmetricOrbit::toReference(g, b);
              Real m = 1;
              for (size_t k = 0; k < d; ++k)
                m *= std::pow(p[static_cast<Eigen::Index>(k)],
                              static_cast<Real>(alpha[k]));
              sum += orbit.getWeight() * m;
            }
          }
          const Real exact = simplexMoment(alpha);
          r(static_cast<Eigen::Index>(e)) =
            (sum - exact) / std::max(std::abs(exact), Real(1e-14));
        }
        return r;
      }

      /**
       * @brief Searches for a rule of strength @p degree with the orbit
       * classes @p config.
       *
       * Restarts are drawn from a deterministically seeded generator, so the
       * same request always produces the same rule.
       */
      static Result solve(Geometry::Polytope::Type g, size_t degree,
        const Configuration& config, size_t maxRestarts = 256,
        Real tolerance = 1e-13, size_t maxIterations = 200,
        unsigned seed = 20260101u)
      {
        const Geometry::Polytope::Traits traits(g);
        const size_t d = traits.getDimension();
        const Eigen::Index n = static_cast<Eigen::Index>(unknownCount(config));
        const Real measure = simplexMoment(std::vector<size_t>(d, 0));

        std::mt19937 rng(seed);
        std::uniform_real_distribution<Real> bary(Real(0.02), Real(0.6));

        Result best;
        for (size_t restart = 0; restart < maxRestarts; ++restart)
        {
          // Start each orbit from random interior parameters and share the
          // element measure equally between the points, so the initial weights
          // are of the right order rather than of the right value.
          Math::Vector<Real> x(n);
          {
            size_t totalPoints = 0;
            for (const auto& pattern : config)
            {
              SymmetricOrbit::Barycentric probe;
              for (size_t i = 0; i < pattern.size(); ++i)
                for (size_t m = 0; m < pattern[i]; ++m)
                  probe.push_back(static_cast<Real>(i + 1));
              totalPoints += SymmetricOrbit(probe, 0).getSize();
            }
            Eigen::Index cursor = 0;
            for (const auto& pattern : config)
            {
              for (size_t i = 0; i + 1 < pattern.size(); ++i)
                x(cursor++) = bary(rng) / static_cast<Real>(pattern.size());
              x(cursor++) = measure / static_cast<Real>(totalPoints);
            }
          }

          Real lambda = 1e-3;
          Math::Vector<Real> r = residual(g, degree, config, x);
          Real cost = r.squaredNorm();
          for (size_t it = 0; it < maxIterations; ++it)
          {
            Math::Matrix<Real> J(r.size(), n);
            for (Eigen::Index j = 0; j < n; ++j)
            {
              const Real h = 1e-7 * std::max(std::abs(x(j)), Real(1));
              Math::Vector<Real> xp = x, xm = x;
              xp(j) += h;
              xm(j) -= h;
              J.col(j) = (residual(g, degree, config, xp)
                        - residual(g, degree, config, xm)) / (2 * h);
            }
            const Math::Matrix<Real> H =
              J.transpose() * J
              + lambda * Math::Matrix<Real>::Identity(n, n);
            const Math::Vector<Real> step = H.ldlt().solve(-J.transpose() * r);
            const Math::Vector<Real> xn = x + step;
            const Math::Vector<Real> rn = residual(g, degree, config, xn);
            const Real cn = rn.squaredNorm();
            if (cn < cost)
            {
              x = xn;
              r = rn;
              cost = cn;
              lambda = std::max(lambda * Real(0.3), Real(1e-12));
              if (std::sqrt(cost) < tolerance)
                break;
            }
            else
            {
              lambda *= 10;
              if (lambda > 1e12)
                break;
            }
          }

          Result candidate;
          candidate.residual = std::sqrt(cost);
          candidate.converged = candidate.residual < tolerance;
          candidate.orbits = toOrbits(config, x);
          candidate.admissible = isAdmissible(g, candidate.orbits);
          candidate.restarts = restart + 1;
          if (candidate.residual < best.residual)
            best = candidate;
          if (candidate.converged && candidate.admissible)
            return candidate;
        }
        return best;
      }

      /// @brief Whether every weight is positive and every point interior.
      static bool isAdmissible(Geometry::Polytope::Type g,
        const std::vector<SymmetricOrbit>& orbits, Real tol = 1e-12)
      {
        for (const auto& orbit : orbits)
        {
          if (!(orbit.getWeight() > Real(0)))
            return false;
          for (const auto& b : orbit.expand())
            for (const auto v : b)
              if (v < tol || v > Real(1) - tol)
                return false;
        }
        (void)g;
        return true;
      }
  };
}

#endif
