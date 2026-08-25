/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Defines SymmetricRuleGenerator, one symmetric rule search for every
 * element.
 */
#ifndef RODIN_QF_SYMMETRICRULEGENERATOR_H
#define RODIN_QF_SYMMETRICRULEGENERATOR_H

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <functional>
#include <atomic>
#include <thread>
#include <map>
#include <mutex>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/SVD>

#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Types.h"

#include "CollapsedBasis.h"
#include "SymmetryGroup.h"

namespace Rodin::QF
{
  /**
   * @brief Finds fully symmetric quadrature rules on any reference element.
   *
   * A symmetric rule is a union of orbits of the element's symmetry group. The
   * search therefore chooses how many orbits of each kind to use --- the
   * decomposition --- and then looks for seed points whose orbits integrate
   * exactly. This is the method of Witherden and Vincent
   * @cite witherden2015identification, and every part of it here is element
   * agnostic:
   *
   * - the orbit kinds come from SymmetryGroup, derived from the element rather
   *   than tabulated per domain;
   * - the moment equations are stated in CollapsedBasis, one orthogonal basis
   *   written once for all seven elements;
   * - the weights are eliminated, being linear once the points are fixed, so
   *   the search runs over seed positions alone;
   * - the Jacobian is analytic throughout, the elimination differentiated by
   *   variable projection.
   *
   * A seed is written as a strictly positive convex combination of the
   * vertices of the region its orbit kind allows, so it is interior for every
   * value of the search variables and the search never has to be told about
   * constraints. The combination is a softmax, which is smooth everywhere ---
   * unlike a clamp, whose derivative vanishes once it binds, freezing any
   * iterate that reaches the boundary.
   *
   * @see SymmetryGroup
   * @see CollapsedBasis
   */
  class SymmetricRuleGenerator
  {
    public:
      /// @brief Which orbit kinds a rule is built from, by index into
      /// SymmetryGroup::strata.
      using Decomposition = std::vector<size_t>;

      /// @brief A rule, and what the search made of it.
      struct Rule
      {
          std::vector<Math::SpatialVector<Real>> points;
          std::vector<Real> weights;
          std::vector<Real> orbitWeights;   ///< One per orbit.
          Decomposition decomposition;
          Real residual = std::numeric_limits<Real>::max();
          bool converged = false;
          bool admissible = false;   ///< Positive weights, interior points.
          size_t restarts = 0;

          size_t getSize() const
          {
            return points.size();
          }

          const Math::SpatialVector<Real>& getPoint(size_t i) const
          {
            return points[i];
          }

          Real getWeight(size_t i) const
          {
            return weights[i];
          }
      };

      /**
       * @brief The moment system of one decomposition, with the weights
       * eliminated.
       */
      class System
      {
        public:
          System(Geometry::Polytope::Type g, size_t degree, Decomposition decomposition)
            : m_g(g),
              m_degree(degree),
              m_d(Geometry::Polytope::Traits(g).getDimension()),
              m_decomposition(std::move(decomposition))
          {
            const auto& strata = SymmetryGroup::strata(m_g);
            m_no = static_cast<Eigen::Index>(m_decomposition.size());

            m_offset.clear();
            size_t variables = 0;
            size_t points = 0;
            for (const size_t which : m_decomposition)
            {
              m_offset.push_back(variables);
              variables += strata[which].domain.size();
              points += strata[which].orbitSize;
            }
            m_nz = static_cast<Eigen::Index>(std::max<size_t>(variables, 1));
            m_points = points;

            // The basis and its norms depend only on the element and the
            // degree, never on the decomposition, and measuring them costs a
            // full quadrature per mode. Recomputing that for each of the many
            // decompositions of a high degree dwarfs the search itself.
            const Normalisation& shared = normalisation(m_g, m_degree);
            m_basis = shared.basis;
            m_invNorm = shared.inverseNorm;
            m_ne = static_cast<Eigen::Index>(m_basis.size());
            const Real measure = shared.measure;
            m_average = measure / static_cast<Real>(points);

            m_A.resize(m_ne, m_no);
            m_dA.assign(static_cast<size_t>(m_nz), Math::Matrix<Real>(m_ne, m_no));
            m_omega.resize(m_no);
            m_b = Math::Vector<Real>::Zero(m_ne);
            m_b(0) = measure * m_invNorm[0];
            m_fullBasis = m_basis;
            m_fullInvNorm = m_invNorm;
            m_measure = measure;

            prune();
          }

          Eigen::Index getVariableCount() const
          {
            return m_nz;
          }

          size_t getPointCount() const
          {
            return m_points;
          }

          const Math::Vector<Real>& getWeights() const
          {
            return m_omega;
          }

          Real getRelativeFloor() const
          {
            Real least = std::numeric_limits<Real>::max();
            for (Eigen::Index j = 0; j < m_no; ++j)
              least = std::min(least, m_omega(j) / m_average);
            return least;
          }

          /// @brief A point of an orbit with its derivative in every variable.
          struct Node
          {
              Math::SpatialVector<Real> x;
              Math::Matrix<Real> dx;
          };

          /**
           * @brief Expands the search variables into the points of each orbit.
           *
           * The seed of an orbit is a softmax-weighted combination of the
           * vertices of the region its kind allows, hence strictly interior
           * whatever the variables are, and the orbit is that seed carried
           * through the stratum's coset representatives.
           */
          std::vector<std::vector<Node>> expand(const Math::Vector<Real>& z) const
          {
            const auto& strata = SymmetryGroup::strata(m_g);
            std::vector<std::vector<Node>> out(m_decomposition.size());
            for (size_t j = 0; j < m_decomposition.size(); ++j)
            {
              const auto& stratum = strata[m_decomposition[j]];
              const size_t n = stratum.domain.size();
              const size_t at = m_offset[j];

              // Softmax, shifted by its greatest entry so the exponentials
              // cannot overflow.
              Real largest = -std::numeric_limits<Real>::max();
              for (size_t i = 0; i < n; ++i)
                largest = std::max(largest, z(static_cast<Eigen::Index>(at + i)));
              std::vector<Real> share(n);
              Real total = 0;
              for (size_t i = 0; i < n; ++i)
              {
                share[i] = std::exp(z(static_cast<Eigen::Index>(at + i)) - largest);
                total += share[i];
              }
              for (size_t i = 0; i < n; ++i)
                share[i] /= total;

              Math::Vector<Real> seed =
                Math::Vector<Real>::Zero(static_cast<Eigen::Index>(m_d));
              for (size_t i = 0; i < n; ++i)
                seed += share[i] * stratum.domain[i];

              // d(seed)/d(z_i) = share_i (vertex_i - seed).
              Math::Matrix<Real> dseed =
                Math::Matrix<Real>::Zero(static_cast<Eigen::Index>(m_d), m_nz);
              for (size_t i = 0; i < n; ++i)
                dseed.col(static_cast<Eigen::Index>(at + i)) =
                  share[i] * (stratum.domain[i] - seed);

              for (const auto& map : stratum.representatives)
              {
                Node node;
                node.x.resize(static_cast<Eigen::Index>(m_d));
                const Math::Vector<Real> image = map(seed);
                for (size_t k = 0; k < m_d; ++k)
                  node.x[static_cast<Eigen::Index>(k)] =
                    image(static_cast<Eigen::Index>(k));
                node.dx = map.linear * dseed;
                out[j].push_back(std::move(node));
              }
            }
            return out;
          }

          /// @brief Builds the moment matrix, and its derivatives when asked.
          void assemble(const Math::Vector<Real>& z, bool wantDerivatives) const
          {
            const auto orbits = expand(z);
            m_A.setZero();
            for (auto& m : m_dA)
              m.setZero();

            for (Eigen::Index j = 0; j < m_no; ++j)
            {
              for (const auto& node : orbits[static_cast<size_t>(j)])
              {
                // Tabulated once for the point, then read off per mode: the
                // basis is evaluated in full at every point of every orbit on
                // every iteration, so repeating the recurrence per mode is the
                // dominant cost of the whole search.
                const auto table =
                  CollapsedBasis::tabulate(m_g, m_degree, node.x, wantDerivatives);
                for (Eigen::Index e = 0; e < m_ne; ++e)
                {
                  const size_t i = static_cast<size_t>(e);
                  m_A(e, j) +=
                    CollapsedBasis::evaluate(m_g, table, m_basis[i]) * m_invNorm[i];
                  if (!wantDerivatives)
                    continue;
                  const auto grad = CollapsedBasis::gradient(m_g, table, m_basis[i]);
                  for (Eigen::Index l = 0; l < m_nz; ++l)
                  {
                    Real chain = 0;
                    for (size_t r = 0; r < m_d; ++r)
                    {
                      const Eigen::Index rr = static_cast<Eigen::Index>(r);
                      chain += grad[rr] * node.dx(rr, l);
                    }
                    m_dA[static_cast<size_t>(l)](e, j) += chain * m_invNorm[i];
                  }
                }
              }
            }
          }

          /**
           * @brief Evaluates the residual, and the Jacobian when asked.
           *
           * Nothing here is a difference quotient.
           */
          void evaluate(const Math::Vector<Real>& z, Math::Vector<Real>& residual,
            Math::Matrix<Real>* jacobian) const
          {
            assemble(z, jacobian != nullptr);
            m_svd.compute(m_A, Eigen::ComputeThinU | Eigen::ComputeThinV);
            m_omega = m_svd.solve(m_b);

            residual.resize(m_ne + m_no);
            residual.head(m_ne) = m_A * m_omega - m_b;
            // Converging is not enough: a decomposition usually admits several
            // roots, and the ones a bare moment residual finds first routinely
            // carry a negative weight. Making negativity part of what is
            // minimised means convergence implies a usable rule.
            for (Eigen::Index j = 0; j < m_no; ++j)
              residual(m_ne + j) = std::max(Real(0), s_floor - m_omega(j) / m_average);
            if (!jacobian)
              return;

            // Variable projection. With omega itself a function of z,
            //   dr = (I - U U^T) dA w  -  U S^-1 V^T dA^T r,
            //   dw = -V S^-2 V^T dA^T r  -  V S^-1 U^T dA w.
            const auto& U = m_svd.matrixU();
            const auto& V = m_svd.matrixV();
            const auto& sv = m_svd.singularValues();
            const Math::Vector<Real> r = residual.head(m_ne);
            jacobian->resize(m_ne + m_no, m_nz);
            for (Eigen::Index l = 0; l < m_nz; ++l)
            {
              const auto& D = m_dA[static_cast<size_t>(l)];
              const Math::Vector<Real> Dw = D * m_omega;
              const Math::Vector<Real> Dtr = D.transpose() * r;
              const Math::Vector<Real> a = U.transpose() * Dw;
              const Math::Vector<Real> c = V.transpose() * Dtr;
              Math::Vector<Real> sa(m_no), sc(m_no), scc(m_no);
              for (Eigen::Index i = 0; i < m_no; ++i)
              {
                const Real si = (sv(i) > 0) ? sv(i) : Real(1);
                sa(i) = a(i) / si;
                sc(i) = c(i) / si;
                scc(i) = c(i) / (si * si);
              }
              jacobian->col(l).head(m_ne) = Dw - U * a - U * sc;
              const Math::Vector<Real> domega = -(V * scc) - (V * sa);
              for (Eigen::Index j = 0; j < m_no; ++j)
              {
                const bool active = (s_floor - m_omega(j) / m_average) > 0;
                (*jacobian)(m_ne + j, l) = active ? -domega(j) / m_average : Real(0);
              }
            }
          }

          /**
           * @brief Residual of the whole moment system, pruned modes included.
           *
           * The search works with the pruned system; this is what says the
           * pruning was sound.
           */
          Real fullResidual(const Math::Vector<Real>& z) const
          {
            Math::Vector<Real> residual;
            evaluate(z, residual, nullptr);
            const auto orbits = expand(z);
            Real worst = 0;
            for (size_t e = 0; e < m_fullBasis.size(); ++e)
            {
              Real got = 0;
              for (Eigen::Index j = 0; j < m_no; ++j)
                for (const auto& node : orbits[static_cast<size_t>(j)])
                  got += m_omega(j) *
                    CollapsedBasis::evaluate(m_g, m_fullBasis[e], node.x) *
                    m_fullInvNorm[e];
              const Real want = (e == 0) ? m_measure * m_fullInvNorm[0] : Real(0);
              worst = std::max(worst, std::abs(got - want));
            }
            return worst;
          }

        private:
          /**
           * @brief Drops the modes that symmetry satisfies identically.
           *
           * Summing a basis function over a symmetry orbit is its
           * symmetrisation, which vanishes for every mode carrying no
           * invariant component. Those rows are zero whatever the seeds are,
           * so they constrain nothing and merely enlarge a system that is
           * rebuilt at every point of every orbit on every iteration. On a
           * hexahedron of strength twenty they are the overwhelming majority,
           * the group having forty-eight elements.
           *
           * Which rows those are is discovered by probing rather than derived,
           * so no table of invariant modes has to be written per element.
           */
          void prune()
          {
            std::mt19937 rng(97531u);
            std::normal_distribution<Real> gauss(0, 1.5);
            Math::Vector<Real> scale = Math::Vector<Real>::Zero(m_ne);
            for (size_t probe = 0; probe < 4; ++probe)
            {
              Math::Vector<Real> z(m_nz);
              for (Eigen::Index i = 0; i < m_nz; ++i)
                z(i) = gauss(rng);
              assemble(z, false);
              for (Eigen::Index e = 0; e < m_ne; ++e)
                scale(e) = std::max(scale(e), m_A.row(e).cwiseAbs().maxCoeff());
            }

            std::vector<std::vector<size_t>> basis;
            std::vector<Real> inverseNorm;
            std::vector<Real> b;
            for (Eigen::Index e = 0; e < m_ne; ++e)
            {
              if (scale(e) <= 1e-11 && std::abs(m_b(e)) <= 1e-11)
                continue;
              basis.push_back(m_basis[static_cast<size_t>(e)]);
              inverseNorm.push_back(m_invNorm[static_cast<size_t>(e)]);
              b.push_back(m_b(e));
            }
            m_basis = std::move(basis);
            m_invNorm = std::move(inverseNorm);
            m_ne = static_cast<Eigen::Index>(m_basis.size());
            m_b.resize(m_ne);
            for (Eigen::Index e = 0; e < m_ne; ++e)
              m_b(e) = b[static_cast<size_t>(e)];
            m_A.resize(m_ne, m_no);
            m_dA.assign(static_cast<size_t>(m_nz), Math::Matrix<Real>(m_ne, m_no));
          }

          /// @brief Basis norms, shared by every decomposition of a degree.
          struct Normalisation
          {
              std::vector<std::vector<size_t>> basis;
              std::vector<Real> inverseNorm;
              Real measure = 1;
          };

          static const Normalisation& normalisation(
            Geometry::Polytope::Type g, size_t degree)
          {
            static std::mutex guard;
            static std::map<std::pair<Geometry::Polytope::Type, size_t>, Normalisation>
              cache;
            const std::lock_guard<std::mutex> lock(guard);
            const auto key = std::make_pair(g, degree);
            const auto found = cache.find(key);
            if (found != cache.end())
              return found->second;

            Normalisation built;
            built.basis = CollapsedBasis::indices(g, degree);
            std::vector<Math::SpatialVector<Real>> qp;
            std::vector<Real> qw;
            CollapsedBasis::exactRule(g, degree + 2, qp, qw);
            built.measure = std::accumulate(qw.begin(), qw.end(), Real(0));
            // Tabulated per quadrature point, so the Jacobi recurrences run
            // once for the whole basis rather than once per mode.
            std::vector<Real> square(built.basis.size(), 0);
            for (size_t q = 0; q < qw.size(); ++q)
            {
              const auto table = CollapsedBasis::tabulate(g, degree, qp[q], false);
              for (size_t e = 0; e < built.basis.size(); ++e)
              {
                const Real v = CollapsedBasis::evaluate(g, table, built.basis[e]);
                square[e] += qw[q] * v * v;
              }
            }
            built.inverseNorm.resize(built.basis.size());
            for (size_t e = 0; e < built.basis.size(); ++e)
              built.inverseNorm[e] = 1 / std::sqrt(std::max(square[e], Real(1e-300)));
            return cache.emplace(key, std::move(built)).first->second;
          }

          Geometry::Polytope::Type m_g;
          size_t m_degree;
          size_t m_d;
          Decomposition m_decomposition;
          std::vector<size_t> m_offset;
          Eigen::Index m_nz = 0;
          Eigen::Index m_ne = 0;
          Eigen::Index m_no = 0;
          size_t m_points = 0;
          Real m_average = 1;
          Real m_measure = 1;
          std::vector<std::vector<size_t>> m_basis;
          std::vector<Real> m_invNorm;
          std::vector<std::vector<size_t>> m_fullBasis;
          std::vector<Real> m_fullInvNorm;

          mutable Math::Matrix<Real> m_A;
          mutable std::vector<Math::Matrix<Real>> m_dA;
          mutable Math::Vector<Real> m_b;
          mutable Math::Vector<Real> m_omega;
          mutable Eigen::JacobiSVD<Math::Matrix<Real>> m_svd;
      };

      /// @brief Smallest weight, relative to the equal-weight value, a rule may
      /// carry.
      static constexpr Real s_floor = 1e-4;

      /// @brief How far short of the counting condition a decomposition may
      /// fall and still be tried.
      static constexpr size_t s_slack = 2;

      /**
       * @brief Solves one decomposition.
       *
       * Levenberg--Marquardt from pseudo-random seeds, restarted until it finds
       * a rule that is both exact and positive. Deterministically seeded, so a
       * given (element, degree, decomposition) always yields the same rule.
       */
      static Rule solve(Geometry::Polytope::Type g, size_t degree,
        const Decomposition& decomposition, size_t maxRestarts = 64,
        Real tolerance = 1e-12, size_t maxIterations = 400, unsigned seed = 20260101u)
      {
        const System problem(g, degree, decomposition);
        const Eigen::Index nz = problem.getVariableCount();

        std::mt19937 rng(seed);
        std::normal_distribution<Real> gauss(0, 1.0);

        Rule best;
        best.decomposition = decomposition;
        for (size_t restart = 0; restart < maxRestarts; ++restart)
        {
          Math::Vector<Real> z(nz);
          for (Eigen::Index i = 0; i < nz; ++i)
            z(i) = gauss(rng);

          Real lambda = 1e-3;
          Math::Vector<Real> r;
          Math::Matrix<Real> J;
          problem.evaluate(z, r, &J);
          Real cost = r.squaredNorm();
          for (size_t it = 0; it < maxIterations && std::sqrt(cost) > tolerance; ++it)
          {
            const Math::Matrix<Real> H =
              J.transpose() * J + lambda * Math::Matrix<Real>::Identity(nz, nz);
            const Math::Vector<Real> step = H.ldlt().solve(-J.transpose() * r);
            const Math::Vector<Real> zn = z + step;
            Math::Vector<Real> rn;
            problem.evaluate(zn, rn, nullptr);
            const Real cn = rn.squaredNorm();
            if (cn < cost)
            {
              z = zn;
              cost = cn;
              problem.evaluate(z, r, &J);
              lambda = std::max(lambda * Real(0.3), Real(1e-14));
            }
            else
            {
              lambda *= 10;
              if (lambda > 1e12)
                break;
            }
          }

          // Accepting at the tolerance leaves a rule that merely clears the
          // bar; a few more steps of the same iteration usually take it to
          // rounding, since the Jacobian is exact and the iteration converges
          // quadratically once close. The cost is negligible -- it happens
          // once, on a rule already found -- and the difference shows up
          // directly in how exactly the table integrates.
          if (std::sqrt(cost) < tolerance)
          {
            Math::Vector<Real> polished = z;
            Real best = cost;
            for (size_t it = 0; it < 24; ++it)
            {
              Math::Vector<Real> r2;
              Math::Matrix<Real> J2;
              problem.evaluate(polished, r2, &J2);
              const Math::Matrix<Real> H =
                J2.transpose() * J2 + Real(1e-14) * Math::Matrix<Real>::Identity(nz, nz);
              const Math::Vector<Real> next =
                polished + H.ldlt().solve(-J2.transpose() * r2);
              Math::Vector<Real> rn;
              problem.evaluate(next, rn, nullptr);
              const Real cn = rn.squaredNorm();
              if (!(cn < best))
                break;
              polished = next;
              best = cn;
            }
            if (best < cost)
            {
              z = polished;
              cost = best;
            }
          }

          Rule candidate = build(g, degree, decomposition, problem, z);
          candidate.restarts = restart + 1;
          candidate.residual = std::max(std::sqrt(cost), problem.fullResidual(z));
          candidate.converged = candidate.residual < tolerance;
          if (candidate.residual < best.residual)
            best = candidate;
          if (candidate.converged && candidate.admissible)
            return candidate;
        }
        return best;
      }

      /**
       * @brief Every decomposition using exactly @p points points.
       *
       * Witherden and Vincent record that the decomposition matters far more
       * than the initial guess, so the search enumerates them exhaustively for
       * one point count before considering a larger one --- rather than
       * capping the number of orbits, which silently puts the higher degrees
       * out of reach, since a rule of strength twenty needs a good many.
       */
      static std::vector<Decomposition> decompositions(
        Geometry::Polytope::Type g, size_t points)
      {
        const auto& strata = SymmetryGroup::strata(g);
        std::vector<Decomposition> out;
        Decomposition current;
        const std::function<void(size_t, size_t)> rec = [&](size_t from, size_t left) {
          if (left == 0)
          {
            out.push_back(current);
            return;
          }
          for (size_t i = from; i < strata.size(); ++i)
          {
            const size_t size = strata[i].orbitSize;
            if (size > left)
              continue;
            // A kind with no freedom contributes the same point however often
            // it is repeated, so it may appear at most once.
            if (strata[i].getDimension() == 0 &&
              std::count(current.begin(), current.end(), i) > 0)
              continue;
            current.push_back(i);
            rec(i, left - size);
            current.pop_back();
          }
        };
        rec(0, points);
        return out;
      }

      /**
       * @brief Number of independent conditions a symmetric rule must satisfy.
       *
       * A symmetric rule integrates every non-invariant mode correctly whatever
       * its parameters, so the conditions it actually has to meet number the
       * dimension of the invariant subspace of the polynomials of that degree,
       * not the dimension of the whole space.
       *
       * That dimension is given exactly by Burnside's lemma applied to the
       * action on polynomials,
       * @f[
       *   \dim \mathcal{P}_p^G = \frac{1}{|G|} \sum_{g \in G}
       *     \sum_{k=0}^{p} h_k\bigl(\lambda(g)\bigr),
       * @f]
       * where @f$ h_k @f$ is the complete homogeneous symmetric polynomial in
       * the eigenvalues of the linear part of @f$ g @f$ --- read off as the
       * power series of @f$ 1/\det(I - tA_g) @f$, which is Molien's series.
       * The translation part shifts only to lower degree and so leaves the
       * trace alone.
       *
       * Counting instead the modes whose symmetrisation is merely nonzero
       * would badly overcount: almost every mode has *some* invariant
       * component, and the resulting figure is nearly the dimension of the
       * whole space, which rejects exactly the decompositions that work.
       */
      static size_t invariantDimension(Geometry::Polytope::Type g, size_t degree)
      {
        const auto& group = SymmetryGroup::maps(g);
        if (group.empty())
          return 0;
        const Eigen::Index d = group.front().linear.rows();

        Real total = 0;
        for (const auto& map : group)
        {
          // Coefficients of det(I - tA), a polynomial of degree d, by
          // interpolation through d + 1 samples.
          const Eigen::Index n = d + 1;
          Math::Matrix<Real> vandermonde(n, n);
          Math::Vector<Real> sampled(n);
          for (Eigen::Index i = 0; i < n; ++i)
          {
            const Real t = Real(0.25) * static_cast<Real>(i + 1);
            Real power = 1;
            for (Eigen::Index j = 0; j < n; ++j)
            {
              vandermonde(i, j) = power;
              power *= t;
            }
            sampled(i) =
              (Math::Matrix<Real>::Identity(d, d) - t * map.linear).determinant();
          }
          const Math::Vector<Real> c = vandermonde.fullPivLu().solve(sampled);

          // Invert the series: h_k = -sum_{j>=1} c_j h_{k-j}, with c_0 = 1.
          std::vector<Real> h(degree + 1, 0);
          h[0] = 1;
          for (size_t k = 1; k <= degree; ++k)
          {
            Real sum = 0;
            for (Eigen::Index j = 1; j < n && static_cast<size_t>(j) <= k; ++j)
              sum += c(j) * h[k - static_cast<size_t>(j)];
            h[k] = -sum / c(0);
          }
          for (size_t k = 0; k <= degree; ++k)
            total += h[k];
        }
        return static_cast<size_t>(std::llround(total / static_cast<Real>(group.size())));
      }

      /**
       * @brief How hard a decomposition is worth trying.
       *
       * Each orbit contributes its seed's free parameters and one weight, and
       * a rule must satisfy every invariant mode. A decomposition with at
       * least as many unknowns as conditions can be expected to work if any
       * does, and is given the full budget of restarts. One with fewer is
       * overdetermined and can succeed only by a coincidence of the geometry,
       * so it is tried briefly rather than not at all.
       *
       * Not at all would be wrong. The six-point rule on a cube is its face
       * centres: a single orbit with no freedom whatever, one unknown against
       * two conditions, which works precisely because those points are
       * special. Treating the count as a necessary condition and discarding
       * such decompositions loses exactly the most economical rules, and the
       * search returned eight points there instead of six. Being
       * overdetermined also makes them rigid, so a few restarts settle them
       * and the full budget would only re-confirm the same outcome.
       */
      /// @brief Whether any orbit of @p decomposition lies on a face.
      static bool usesBoundary(
        Geometry::Polytope::Type g, const Decomposition& decomposition)
      {
        const auto& strata = SymmetryGroup::strata(g);
        for (const size_t which : decomposition)
          if (strata[which].boundary)
            return true;
        return false;
      }

      static size_t restartBudget(Geometry::Polytope::Type g,
        const Decomposition& decomposition, size_t conditions, size_t maxRestarts)
      {
        const auto& strata = SymmetryGroup::strata(g);
        size_t unknowns = decomposition.size(); // one weight per orbit
        for (const size_t which : decomposition)
          unknowns += strata[which].getFreedom();
        if (unknowns >= conditions)
          return maxRestarts;
        // Short of the count, but only just: a coincidence of the geometry can
        // still close a small gap, and this is where the economical rules live
        // --- the cube's six-point rule is one unknown short of two conditions.
        // Being overdetermined makes these rigid, so a few restarts settle
        // them. Further short than this nothing is ever found, and at the
        // higher strengths there are thousands of such decompositions at every
        // point count; trying them all is the difference between minutes and
        // hours.
        if (unknowns + s_slack >= conditions)
          return std::min<size_t>(maxRestarts, 4);
        // Still not nothing. The count is not a necessary condition, and the
        // rules that violate it are the economical ones: the cube's six-point
        // rule is one unknown short, and the pyramid's twenty-three point rule
        // at strength six is three short of twenty conditions and exists all
        // the same, because the conditions are not independent in the
        // parameters that decomposition offers. Skipping these outright lost
        // both. What keeps the search affordable is the ordering, which puts
        // them last, and the deadline, which stops it.
        return std::min<size_t>(maxRestarts, 2);
      }

      /**
       * @brief Searches for the smallest rule of strength @p degree.
       *
       * Point counts are tried upward from the counting bound, and every
       * decomposition of each is attempted before the count is increased, so
       * the first rule found is the smallest this method can express.
       */
      static Rule search(Geometry::Polytope::Type g, size_t degree,
        size_t maxPoints = 128, size_t maxRestarts = 64, Real tolerance = 1e-12,
        size_t threads = 0, Real seconds = 0)
      {
        const size_t d = Geometry::Polytope::Traits(g).getDimension();
        size_t moments = 1;
        for (size_t i = 1; i <= d; ++i)
          moments = moments * (degree + i) / i;
        const size_t lower = (moments + d) / (d + 1);

        const size_t conditions = invariantDimension(g, degree);
        const size_t workers = (threads > 0)
          ? threads
          : std::max<size_t>(1, std::min<size_t>(std::thread::hardware_concurrency(), 8));

        // Unattended sweeps must not disappear into a single strength, so the
        // search may be given a deadline. Reaching it returns the best rule
        // found rather than nothing, and the caller can tell the difference
        // because such a rule has not converged.
        const auto started = std::chrono::steady_clock::now();
        const auto expired = [&] {
          if (seconds <= 0)
            return false;
          const std::chrono::duration<Real> spent =
            std::chrono::steady_clock::now() - started;
          return spent.count() > seconds;
        };

        Rule best;
        for (size_t points = lower; points <= maxPoints && !expired(); ++points)
        {
          std::vector<Decomposition> candidates = decompositions(g, points);
          if (candidates.empty())
            continue;

          // Ordered by how likely each is to work: those with enough unknowns
          // to meet the conditions first, and among them the ones using fewest
          // orbits, which is what published rules look like and what solves
          // fastest. Since a success stops the others, the order decides how
          // much of the list is ever examined -- the difference between
          // settling a high strength in seconds and grinding through thousands
          // of hopeless decompositions first. The ordering is fixed, so the
          // rule found does not depend on it having been reached quickly.
          const auto unknownsOf = [&](const Decomposition& decomposition) {
            const auto& strata = SymmetryGroup::strata(g);
            size_t unknowns = decomposition.size();
            for (const size_t which : decomposition)
              unknowns += strata[which].getFreedom();
            return unknowns;
          };
          std::stable_sort(candidates.begin(), candidates.end(),
            [&](const Decomposition& a, const Decomposition& b) {
              const size_t ua = unknownsOf(a), ub = unknownsOf(b);
              const bool pa = ua >= conditions, pb = ub >= conditions;
              if (pa != pb)
                return pa;
              // Interior decompositions before those placing an orbit on a
              // face. Boundary orbits are admitted because a few published
              // rules need them, but they are the exception, and letting them
              // multiply the list ahead of the ordinary ones buries what
              // usually works: the pyramid at strength six found nothing until
              // the interior shapes were tried first.
              const bool ba = usesBoundary(g, a), bb = usesBoundary(g, b);
              if (ba != bb)
                return bb;
              // Among the workable ones, the tightest fit first. A rule that
              // exists generally has barely enough freedom to exist: published
              // rules sit just above the condition count, and decompositions
              // carrying far more freedom than that are the ones with slack to
              // wander in and rarely close.
              if (pa && ua != ub)
                return ua < ub;
              return a.size() < b.size();
            });

          // Decompositions are independent, so they are tried in parallel. The
          // rule returned is still the one from the earliest decomposition that
          // succeeds, not whichever thread happens to finish first, so the
          // result does not depend on the scheduling.
          std::vector<Rule> results(candidates.size());
          std::atomic<size_t> next{0};
          std::atomic<size_t> earliest{candidates.size()};
          std::vector<std::thread> pool;
          for (size_t w = 0; w < std::min(workers, candidates.size()); ++w)
          {
            pool.emplace_back([&] {
              for (;;)
              {
                const size_t i = next++;
                if (i >= candidates.size() || i > earliest.load() || expired())
                  return;
                size_t budget = restartBudget(g, candidates[i], conditions, maxRestarts);
                // Effort is concentrated on the decompositions the ordering
                // put first. A rule that exists is usually found in one of a
                // handful of shapes, and finding it is a matter of seeding the
                // iteration often enough to land in its basin; spreading the
                // same restarts thinly over hundreds of unlikely shapes
                // instead is how a strength ends up one point count too high.
                if (budget == maxRestarts && i < 8)
                  budget *= 8;
                results[i] = solve(g, degree, candidates[i], budget, tolerance);
                if (results[i].converged && results[i].admissible)
                {
                  size_t seen = earliest.load();
                  while (i < seen && !earliest.compare_exchange_weak(seen, i))
                    ;
                }
              }
            });
          }
          for (auto& worker : pool)
            worker.join();

          for (const auto& candidate : results)
          {
            if (candidate.converged && candidate.admissible)
              return candidate;
            if (candidate.residual < best.residual)
              best = candidate;
          }
        }
        return best;
      }

    private:
      /// @brief Assembles the rule a converged search describes.
      static Rule build(Geometry::Polytope::Type g, size_t degree,
        const Decomposition& decomposition, const System& problem,
        const Math::Vector<Real>& z)
      {
        Rule out;
        out.decomposition = decomposition;
        const auto orbits = problem.expand(z);
        const auto& omega = problem.getWeights();
        out.admissible = problem.getRelativeFloor() > s_floor;
        for (size_t j = 0; j < orbits.size(); ++j)
        {
          out.orbitWeights.push_back(omega(static_cast<Eigen::Index>(j)));
          for (const auto& node : orbits[j])
          {
            out.points.push_back(node.x);
            out.weights.push_back(omega(static_cast<Eigen::Index>(j)));
          }
        }
        return out;
      }
  };
}

#endif
