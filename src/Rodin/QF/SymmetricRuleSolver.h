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

#include <limits>
#include <algorithm>
#include <random>
#include <Eigen/SVD>
#include <Eigen/Cholesky>
#include <functional>
#include <numeric>

#include "SymmetricOrbit.h"
#include "GaussLegendre.h"
#include "CollapsedBasis.h"
#include "OrthogonalBasis.h"

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
   * @par Reach, and what this solver is not for
   * Every rule this class can produce is fully symmetric by construction, so
   * its reach is bounded by which point counts a symmetry decomposition can
   * express at all. That bound is real and it is not a defect of the search.
   * The published Xiao--Gimbutas tetrahedron rule of strength four, for
   * instance, has eleven points carrying eleven distinct weights and no
   * symmetry whatsoever --- it comes from node elimination, and no symmetric
   * decomposition of eleven points reproduces it. Measuring this solver
   * against such a table therefore compares it with a rule of a kind it
   * cannot represent.
   *
   * Use NodeElimination for asymmetric tables and this solver for symmetric
   * ones, and take whichever is smaller per degree; the two generators have
   * complementary reach, which is why the published tables of the two
   * families disagree on point counts.
   *
   * @see SymmetricOrbit
   * @see NodeElimination
   */
  class SymmetricRuleSolver
  {
    public:
      /// @brief Partition of the base simplex's vertex count, describing the
      /// multiplicity pattern of one orbit class.
      using Pattern = std::vector<size_t>;

      /**
       * @brief One orbit class of a candidate rule.
       *
       * On a simplex only @c barycentric is used. On a product element the
       * class additionally says how the orbit sits along the tensor
       * direction: @c midPlane places it on the single reflection-invariant
       * value @f$ 1/2 @f$, otherwise it is the reflected pair
       * @f$ \{c, 1-c\} @f$ with @f$ c @f$ a further unknown.
       */
      struct OrbitClass
      {
          Pattern barycentric;
          bool tensor = false;   ///< Element has a tensor direction.
          bool midPlane = true;  ///< Tensor part is {1/2} rather than {c, 1-c}.

          OrbitClass() = default;
          OrbitClass(Pattern p)
            : barycentric(std::move(p))
          {}
          OrbitClass(std::initializer_list<size_t> p)
            : barycentric(p)
          {}
          OrbitClass(Pattern p, bool mid)
            : barycentric(std::move(p)),
              tensor(true),
              midPlane(mid)
          {}

          friend bool operator==(const OrbitClass& a, const OrbitClass& b)
          {
            return a.barycentric == b.barycentric && a.tensor == b.tensor &&
              a.midPlane == b.midPlane;
          }
      };

      /// @brief The orbit classes making up a candidate rule.
      using Configuration = std::vector<OrbitClass>;

      /// @brief Outcome of a solve.
      struct Result
      {
          bool converged = false;             ///< Residual met the tolerance.
          bool admissible = false;            ///< Positive weights, interior points.
          std::vector<SymmetricOrbit> orbits; ///< The rule, as orbits.
          Real residual = std::numeric_limits<Real>::infinity(); ///< Final norm.
          size_t restarts = 0;                ///< Starting points consumed.
      };

      /// @brief Exact moment @f$ \int_K x^\alpha @f$ on the reference element
      /// @p g. The wedge is the unit triangle crossed with @f$ [0,1] @f$, so
      /// its moment is the triangle moment times @f$ 1/(c+1) @f$.
      static Real referenceMoment(
        Geometry::Polytope::Type g, const std::vector<size_t>& alpha)
      {
        if (g == Geometry::Polytope::Type::Pyramid)
        {
          // Reference pyramid {0 <= z <= 1, 0 <= x, y <= 1 - z}:
          //   int x^a y^b z^c = c! (a+b+2)! / [ (a+b+c+3)! (a+1)(b+1) ]
          const auto fact = [](size_t n) {
            Real r = 1;
            for (size_t i = 2; i <= n; ++i)
              r *= static_cast<Real>(i);
            return r;
          };
          return fact(alpha[2]) * fact(alpha[0] + alpha[1] + 2) /
            (fact(alpha[0] + alpha[1] + alpha[2] + 3) *
              static_cast<Real>((alpha[0] + 1) * (alpha[1] + 1)));
        }
        if (g != Geometry::Polytope::Type::Wedge)
          return simplexMoment(alpha);
        const std::vector<size_t> planar{alpha[0], alpha[1]};
        return simplexMoment(planar) / static_cast<Real>(alpha[2] + 1);
      }

      /// @brief Exact moment @f$ \int_K x^\alpha @f$ on the reference
      /// @f$ d @f$-simplex, @f$ \alpha! / (|\alpha| + d)! @f$.
      static Real simplexMoment(const std::vector<size_t>& alpha)
      {
        const auto fact = [](size_t n) {
          Real r = 1;
          for (size_t i = 2; i <= n; ++i)
            r *= static_cast<Real>(i);
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
        const std::function<void(size_t, size_t)> rec = [&](size_t k, size_t budget) {
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
        for (const auto& c : config)
        {
          n += c.barycentric.size();  // (parts - 1) parameters + 1 weight
          if (c.tensor && !c.midPlane)
            ++n;                      // the reflected pair's free coordinate
        }
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
        for (const auto& orbitClass : config)
        {
          const auto& pattern = orbitClass.barycentric;
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
          SymmetricOrbit::Tensor tensor;
          if (orbitClass.tensor)
          {
            if (orbitClass.midPlane)
            {
              tensor = {Real(0.5)};
            }
            else
            {
              const Real c = x(cursor++);
              tensor = {c, Real(1) - c};
            }
          }
          // Weights enter as w = theta^2. Positivity is then a property of
          // the parameterisation rather than something checked after the
          // solve: at higher degrees the moment equations have many exact
          // solutions, most of them with a negative weight, and a solver that
          // only screens afterwards converges to them and discards the result.
          const Real theta = x(cursor++);
          const Real weight = theta * theta;

          SymmetricOrbit::Barycentric barycentric;
          for (size_t i = 0; i < k; ++i)
            for (size_t m = 0; m < pattern[i]; ++m)
              barycentric.push_back(values[i]);
          if (tensor.empty())
            orbits.emplace_back(std::move(barycentric), weight);
          else
            orbits.emplace_back(std::move(barycentric), std::move(tensor), weight);
        }
        return orbits;
      }

      /**
       * @brief Monomial exponents, exact moments and row scalings for one
       * (element, degree) pair.
       *
       * The residual is evaluated O(unknowns) times per Levenberg-Marquardt
       * iteration and was rebuilding all three every time: allocating the
       * exponent list and recomputing factorials per call. They depend only on
       * the element and the degree, so they are built once per solve.
       */
      struct MomentData
      {
          std::vector<std::vector<size_t>> alphas;
          std::vector<Real> exact;
          std::vector<Real> scale;

          MomentData(Geometry::Polytope::Type g, size_t d, size_t degree)
            : alphas(monomials(d, degree))
          {
            exact.reserve(alphas.size());
            scale.reserve(alphas.size());
            for (const auto& alpha : alphas)
            {
              const Real m = referenceMoment(g, alpha);
              exact.push_back(m);
              scale.push_back(Real(1) / std::max(std::abs(m), Real(1e-14)));
            }

            const Eigen::Index n = static_cast<Eigen::Index>(alphas.size());
            Math::Matrix<Real> gram(n, n);
            std::vector<size_t> sum(d);
            for (Eigen::Index i = 0; i < n; ++i)
            {
              for (Eigen::Index j = 0; j < n; ++j)
              {
                for (size_t k = 0; k < d; ++k)
                  sum[k] =
                    alphas[static_cast<size_t>(i)][k] + alphas[static_cast<size_t>(j)][k];
                gram(i, j) = referenceMoment(g, sum);
              }
            }
            Eigen::LLT<Math::Matrix<Real>> llt(gram);
            if (llt.info() == Eigen::Success)
              chol = llt.matrixL();
          }

          /**
           * @brief Cholesky factor of the monomial Gram matrix
           * @f$ G_{\alpha\beta} = \int_K x^{\alpha+\beta} @f$.
           *
           * Normalising each moment equation by its own value fixes the row
           * scaling but not the conditioning of the monomial basis itself,
           * and that is what stalls the solve at higher degree. Multiplying
           * the residual by @f$ L^{-1} @f$ expresses it in an orthonormal
           * basis of the same space, which is what Xiao--Gimbutas and
           * Witherden--Vincent both do.
           *
           * Left empty when the factorisation fails, so the caller falls back
           * to the scaled monomial residual rather than using a corrupted
           * transform.
           */
          Math::Matrix<Real> chol;
      };

      /// @brief Integer power, avoiding std::pow in the residual hot loop.
      static Real ipow(Real x, size_t e)
      {
        Real r = 1;
        while (e--)
          r *= x;
        return r;
      }

      /// @brief Normalised moment residual of the rule described by @p x.
      ///
      /// The orbits are expanded once per evaluation, not once per monomial:
      /// the residual is called O(unknowns) times per Levenberg-Marquardt
      /// iteration for the finite-difference Jacobian, so an expansion inside
      /// the monomial loop dominates everything else.
      static Math::Vector<Real> residual(Geometry::Polytope::Type g,
        const Configuration& config, const Math::Vector<Real>& x,
        const MomentData& moments, size_t d)
      {
        std::vector<Math::SpatialVector<Real>> points;
        std::vector<Real> weights;
        for (const auto& orbit : toOrbits(config, x))
        {
          for (auto& p : orbit.expandPoints(g))
          {
            points.push_back(std::move(p));
            weights.push_back(orbit.getWeight());
          }
        }

        Math::Vector<Real> r(static_cast<Eigen::Index>(moments.alphas.size()));
        for (size_t e = 0; e < moments.alphas.size(); ++e)
        {
          const auto& alpha = moments.alphas[e];
          Real sum = 0;
          for (size_t q = 0; q < points.size(); ++q)
          {
            Real m = weights[q];
            for (size_t k = 0; k < d; ++k)
              m *= ipow(points[q][static_cast<Eigen::Index>(k)], alpha[k]);
            sum += m;
          }
          r(static_cast<Eigen::Index>(e)) = (sum - moments.exact[e]) * moments.scale[e];
        }
        return r;
      }

      /// @brief Convenience overload building the moment data on the spot.
      static Math::Vector<Real> residual(Geometry::Polytope::Type g, size_t degree,
        const Configuration& config, const Math::Vector<Real>& x)
      {
        const size_t d = Geometry::Polytope::Traits(g).getDimension();
        return residual(g, config, x, MomentData(g, d, degree), d);
      }

      /// Smallest weight, as a fraction of the equal-weight value, that a rule
      /// may carry. Published rules sit orders of magnitude above this.
      static constexpr Real s_floor = 1e-4;

      /// Total number of points a configuration expands to.
      static size_t expandOf(const Configuration& config, Geometry::Polytope::Type g)
      {
        size_t n = 0;
        for (const auto& c : config)
          n += classSize(c);
        return n;
      }

      /// Maps a search variable onto the open unit interval.
      static Real logistic(Real z)
      {
        return (z >= 0) ? 1 / (1 + std::exp(-z)) : std::exp(z) / (1 + std::exp(z));
      }

      /**
       * @brief The moment system of a configuration with the weights
       * eliminated.
       *
       * This is the formulation of Witherden and Vincent
       * @cite witherden2015identification, and it differs from solve() in the
       * decisive place. Once the points are prescribed the moment equations
       * are *linear* in the weights, so the weights need not be search
       * variables at all: they are recovered as the least squares solution of
       * a small linear system, and the nonlinear search runs over the orbit
       * shape parameters alone. Weights are roughly a third of the degrees of
       * freedom of a symmetric decomposition, and removing them greatly
       * reduces the iterations needed to converge --- more than paying for the
       * decomposition solved per evaluation.
       *
       * Writing @f$ A_{ij} = \sum_{x \in S_j} \psi_i(x) @f$ for the sum of
       * basis function @f$ i @f$ over orbit @f$ j @f$, and @f$ b_i @f$ for its
       * integral, the weights are @f$ \omega = A^+ b @f$ and the residual is
       * @f$ A\omega - b @f$.
       *
       * Three details decide whether this reaches machine precision.
       *
       * First, the residual is stated in an orthonormal basis, by applying the
       * inverse Cholesky factor of the monomial Gram matrix. Row scaling alone
       * leaves the monomial basis ill-conditioned, and that is what otherwise
       * limits the reachable degree.
       *
       * Second, orbit parameters are mapped from unconstrained search
       * variables by a stick-breaking transform rather than clamped. A clamp
       * has zero derivative once it binds, so a Jacobian column collapses and
       * the iterate is stuck where it stands; stick-breaking keeps every
       * iterate strictly inside the element and everywhere differentiable.
       *
       * Third, the Jacobian is analytic. The whole chain --- stick-breaking,
       * the barycentric expansion, the affine map to the reference element,
       * the monomials, and the elimination of @f$ \omega @f$ itself --- is
       * differentiable in closed form, the elimination by the variable
       * projection formula of Golub and Pereyra. A central-difference Jacobian
       * caps the attainable residual near @f$ 10^{-10} @f$ however the step is
       * tuned, which is a rule good to ten digits rather than sixteen.
       *
       * Convergence is judged on the residual *and* on positivity: a
       * configuration generally admits several roots, and the ones a bare
       * moment residual finds first routinely carry a negative weight. The
       * barrier makes negativity part of what the search minimises rather than
       * a verdict passed after the fact.
       *
       * @see verifyJacobian, which checks the analytic derivative against a
       * difference quotient.
       */
      class Eliminated
      {
        public:
          /// A point of an orbit together with its derivative in every
          /// search variable.
          struct Node
          {
              Math::SpatialVector<Real> x;
              Math::Matrix<Real> dx;   ///< Spatial dimension by variable count.
          };

          Eliminated(Geometry::Polytope::Type g, size_t degree, Configuration config)
            : m_g(g),
              m_d(Geometry::Polytope::Traits(g).getDimension()),
              m_config(std::move(config))
          {
            size_t shape = 0;
            size_t points = 0;
            for (const auto& c : m_config)
            {
              shape += c.barycentric.size() - 1;
              if (c.tensor && !c.midPlane)
                ++shape;
              points += classSize(c);
            }
            // A configuration with no free parameter still needs one variable
            // for the solver to hold; it simply has nothing to vary.
            m_degenerate = (shape == 0);
            m_nz = static_cast<Eigen::Index>(m_degenerate ? 1 : shape);
            m_no = static_cast<Eigen::Index>(m_config.size());

            // One orthogonal basis serves every element. Orthonormalising
            // monomials instead would mean a Cholesky of their Gram matrix,
            // which is Hilbert-like: already marginal near degree eight and
            // meaningless beyond, and that -- not the search -- is what caps
            // the attainable residual there. It is also the expensive part of
            // the setup, being quadratic in the number of moments.
            m_degree = degree;
            m_basis = CollapsedBasis::indices(m_g, degree);
            m_ne = static_cast<Eigen::Index>(m_basis.size());

            // Norms come from a collapsed Gauss--Jacobi rule, whose weights
            // are positive and which is exact well past the degree needed, so
            // there is no cancellation to amplify.
            std::vector<Math::SpatialVector<Real>> qp;
            std::vector<Real> qw;
            CollapsedBasis::exactRule(m_g, degree + 2, qp, qw);
            const Real measure = std::accumulate(qw.begin(), qw.end(), Real(0));
            m_average = measure / static_cast<Real>(points);

            m_invNorm.assign(m_basis.size(), 1);
            for (size_t e = 0; e < m_basis.size(); ++e)
            {
              Real square = 0;
              for (size_t q = 0; q < qw.size(); ++q)
              {
                const Real v = CollapsedBasis::evaluate(m_g, m_basis[e], qp[q]);
                square += qw[q] * v * v;
              }
              m_invNorm[e] = 1 / std::sqrt(std::max(square, Real(1e-300)));
            }

            // The integral of every non-constant mode is zero exactly, by
            // orthogonality to the constant one. Taking these numerically
            // instead would put quadrature noise where an exact zero belongs,
            // and the search would then solve equations that are quietly
            // wrong -- converging to a small residual against the wrong right
            // hand side, which no amount of searching detects.
            m_A.resize(m_ne, m_no);
            m_dA.assign(static_cast<size_t>(m_nz), Math::Matrix<Real>(m_ne, m_no));
            m_omega.resize(m_no);
            m_b = Math::Vector<Real>::Zero(m_ne);
            m_b(0) = measure * m_invNorm[0];

            prune();
          }

          /**
           * @brief Drops the modes that symmetry satisfies identically.
           *
           * Summing a basis function over a symmetry orbit is its
           * symmetrisation, which vanishes for every mode carrying no
           * invariant component. Those rows are zero whatever the orbit
           * parameters are, so they constrain nothing and merely enlarge the
           * system --- and the system is rebuilt at every point of every orbit
           * on every iteration. On a tetrahedron of strength twenty they are
           * the overwhelming majority: the symmetry group has twenty-four
           * elements, and the invariant subspace shrinks roughly in
           * proportion.
           *
           * Which rows those are is discovered rather than derived, by
           * probing at random parameters, so no per-element table of invariant
           * modes has to be written or maintained. A row that survives every
           * probe is kept. The pruning is verified before any rule is
           * accepted, by measuring the residual of the *full* system, so a row
           * dropped in error cannot pass unnoticed.
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
                z(i) = m_degenerate ? Real(0) : gauss(rng);
              assemble(z, false);
              for (Eigen::Index e = 0; e < m_ne; ++e)
                scale(e) = std::max(scale(e), m_A.row(e).cwiseAbs().maxCoeff());
            }

            m_keep.clear();
            for (Eigen::Index e = 0; e < m_ne; ++e)
              if (scale(e) > 1e-11 || std::abs(m_b(e)) > 1e-11)
                m_keep.push_back(e);

            m_full = m_ne;
            m_fullBasis = m_basis;
            m_fullInvNorm = m_invNorm;
            std::vector<std::vector<size_t>> basis;
            std::vector<Real> inverseNorm;
            Math::Vector<Real> b(static_cast<Eigen::Index>(m_keep.size()));
            for (size_t i = 0; i < m_keep.size(); ++i)
            {
              const size_t e = static_cast<size_t>(m_keep[i]);
              basis.push_back(m_basis[e]);
              inverseNorm.push_back(m_invNorm[e]);
              b(static_cast<Eigen::Index>(i)) = m_b(m_keep[i]);
            }
            m_basis = std::move(basis);
            m_invNorm = std::move(inverseNorm);
            m_b = std::move(b);
            m_ne = static_cast<Eigen::Index>(m_keep.size());
            m_A.resize(m_ne, m_no);
            m_dA.assign(static_cast<size_t>(m_nz), Math::Matrix<Real>(m_ne, m_no));
          }

          /**
           * @brief Residual of the *whole* moment system, pruned rows included.
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
              const Real want = (e == 0) ? m_b(0) : Real(0);
              worst = std::max(worst, std::abs(got - want));
            }
            return worst;
          }

          /// @brief Number of modes before pruning.
          Eigen::Index getFullCount() const
          {
            return m_full;
          }

          /// Number of search variables.
          Eigen::Index getVariableCount() const
          {
            return m_nz;
          }

          /// Number of residual entries.
          Eigen::Index getResidualCount() const
          {
            return m_ne + m_no;
          }

          /// Whether the configuration has any free parameter at all.
          bool isDegenerate() const
          {
            return m_degenerate;
          }

          /// The weights recovered by the most recent evaluation.
          const Math::Vector<Real>& getWeights() const
          {
            return m_omega;
          }

          /// The smallest weight the most recent evaluation produced,
          /// relative to the equal-weight value.
          Real getRelativeFloor() const
          {
            Real least = std::numeric_limits<Real>::max();
            for (Eigen::Index j = 0; j < m_no; ++j)
              least = std::min(least, m_omega(j) / m_average);
            return least;
          }

          /**
           * @brief Expands a search vector into the points of each orbit,
           * carrying the derivative of every point along the way.
           */
          std::vector<std::vector<Node>> expand(const Math::Vector<Real>& z) const
          {
            std::vector<std::vector<Node>> orbits(m_config.size());
            Eigen::Index cursor = 0;
            for (size_t c = 0; c < m_config.size(); ++c)
            {
              const auto& pattern = m_config[c].barycentric;
              const size_t k = pattern.size();

              // Stick-breaking: each part takes a fraction of the barycentric
              // budget still unspent, so the coordinates sum to one by
              // construction and stay positive for every value of z.
              std::vector<Real> value(k);
              std::vector<Math::Vector<Real>> dvalue(k, Math::Vector<Real>::Zero(m_nz));
              Real rest = 1;
              Math::Vector<Real> drest = Math::Vector<Real>::Zero(m_nz);
              for (size_t i = 0; i + 1 < k; ++i)
              {
                const Eigen::Index slot = cursor++;
                const Real fraction = logistic(z(slot));
                const Real dfraction = fraction * (1 - fraction);
                const Real m = static_cast<Real>(pattern[i]);
                value[i] = rest * fraction / m;
                dvalue[i] = drest * (fraction / m);
                dvalue[i](slot) += rest * dfraction / m;
                const Math::Vector<Real> keep = drest * (1 - fraction);
                drest = keep;
                drest(slot) -= rest * dfraction;
                rest *= (1 - fraction);
              }
              value[k - 1] = rest / static_cast<Real>(pattern.back());
              dvalue[k - 1] = drest / static_cast<Real>(pattern.back());

              // The tensor direction, when the element has one.
              std::vector<Real> tvalue;
              std::vector<Math::Vector<Real>> dtvalue;
              if (m_config[c].tensor)
              {
                if (m_config[c].midPlane)
                {
                  tvalue = {Real(0.5)};
                  dtvalue = {Math::Vector<Real>::Zero(m_nz)};
                }
                else
                {
                  const Eigen::Index slot = cursor++;
                  const Real t = logistic(z(slot));
                  Math::Vector<Real> row = Math::Vector<Real>::Zero(m_nz);
                  row(slot) = t * (1 - t);
                  tvalue = {t, 1 - t};
                  dtvalue = {row, Math::Vector<Real>(-row)};
                }
              }

              // Distinct permutations of the multiset of parts. Which part a
              // slot came from is tracked by index rather than by value, so
              // the derivative is attributed correctly even where two parts
              // happen to coincide.
              std::vector<size_t> index;
              for (size_t i = 0; i < k; ++i)
                for (size_t m = 0; m < pattern[i]; ++m)
                  index.push_back(i);
              std::sort(index.begin(), index.end());

              const bool product = m_config[c].tensor;
              const Geometry::Polytope::Type base =
                product ? SymmetricOrbit::baseSimplex(index.size()) : m_g;
              const Geometry::Polytope::Traits baseTraits(base);
              const size_t bd = baseTraits.getDimension();

              do
              {
                // The map to the reference element is affine in the
                // barycentric tuple, so the derivative maps the same way.
                Math::SpatialVector<Real> x;
                x.resize(static_cast<Eigen::Index>(m_d));
                x.setZero();
                Math::Matrix<Real> dx =
                  Math::Matrix<Real>::Zero(static_cast<Eigen::Index>(m_d), m_nz);
                for (size_t i = 0; i < index.size(); ++i)
                {
                  const auto& vertex = baseTraits.getVertex(i);
                  for (size_t r = 0; r < bd; ++r)
                  {
                    const Eigen::Index rr = static_cast<Eigen::Index>(r);
                    x[rr] += value[index[i]] * vertex[rr];
                    dx.row(rr) += dvalue[index[i]].transpose() * vertex[rr];
                  }
                }
                if (!product)
                {
                  orbits[c].push_back({x, dx});
                }
                else
                {
                  for (size_t t = 0; t < tvalue.size(); ++t)
                  {
                    Node node{x, dx};
                    const Eigen::Index last = static_cast<Eigen::Index>(m_d - 1);
                    node.x[last] = tvalue[t];
                    node.dx.row(last) = dtvalue[t].transpose();
                    orbits[c].push_back(std::move(node));
                  }
                }
              } while (std::next_permutation(index.begin(), index.end()));
            }
            return orbits;
          }

          /// @brief Builds the moment matrix, and its derivatives when asked.
          void assemble(const Math::Vector<Real>& z, bool wantDerivatives) const
          {
            const auto orbits = expand(z);
            m_A.setZero();
            for (auto& m : m_dA)
              m.setZero();
            const bool jacobian = wantDerivatives;
            for (Eigen::Index j = 0; j < m_no; ++j)
            {
              for (const auto& node : orbits[static_cast<size_t>(j)])
              {
                // Tabulated once for the point, then read off per mode: the
                // basis is evaluated in full at every point of every orbit on
                // every iteration, so repeating the recurrence per mode is the
                // dominant cost of the whole search.
                const auto table =
                  CollapsedBasis::tabulate(m_g, m_degree, node.x, jacobian);
                for (Eigen::Index e = 0; e < m_ne; ++e)
                {
                  const size_t i = static_cast<size_t>(e);
                  m_A(e, j) +=
                    CollapsedBasis::evaluate(m_g, table, m_basis[i]) * m_invNorm[i];
                  if (!jacobian)
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
           * @brief Evaluates the residual, and the Jacobian when @p jacobian
           * is given.
           *
           * Nothing in the Jacobian is a difference quotient.
           */
          void evaluate(const Math::Vector<Real>& z, Math::Vector<Real>& residual,
            Math::Matrix<Real>* jacobian) const
          {
            assemble(z, jacobian != nullptr);

            m_svd.compute(m_A, Eigen::ComputeThinU | Eigen::ComputeThinV);
            m_omega = m_svd.solve(m_b);

            residual.resize(m_ne + m_no);
            residual.head(m_ne) = m_A * m_omega - m_b;
            for (Eigen::Index j = 0; j < m_no; ++j)
              residual(m_ne + j) = std::max(Real(0), s_floor - m_omega(j) / m_average);
            if (!jacobian)
              return;

            // Variable projection. With omega itself a function of z,
            //   dr   = (I - U U^T) dA w  -  U S^-1 V^T dA^T r,
            //   dw   = -V S^-2 V^T dA^T r  -  V S^-1 U^T dA w.
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

          /// Rebuilds the orbits of the configuration, carrying the weights
          /// the most recent evaluation produced.
          std::vector<SymmetricOrbit> toOrbits(const Math::Vector<Real>& z) const
          {
            std::vector<SymmetricOrbit> out;
            Eigen::Index cursor = 0;
            for (size_t c = 0; c < m_config.size(); ++c)
            {
              const auto& pattern = m_config[c].barycentric;
              const size_t k = pattern.size();
              std::vector<Real> value(k);
              Real rest = 1;
              for (size_t i = 0; i + 1 < k; ++i)
              {
                const Real fraction = logistic(z(cursor++));
                value[i] = rest * fraction / static_cast<Real>(pattern[i]);
                rest *= (1 - fraction);
              }
              value[k - 1] = rest / static_cast<Real>(pattern.back());

              SymmetricOrbit::Tensor tensor;
              if (m_config[c].tensor)
              {
                if (m_config[c].midPlane)
                  tensor = {Real(0.5)};
                else
                {
                  const Real t = logistic(z(cursor++));
                  tensor = {t, 1 - t};
                }
              }

              SymmetricOrbit::Barycentric bary;
              for (size_t i = 0; i < k; ++i)
                for (size_t m = 0; m < pattern[i]; ++m)
                  bary.push_back(value[i]);
              const Real w = m_omega(static_cast<Eigen::Index>(c));
              out.push_back(tensor.empty()
                  ? SymmetricOrbit(std::move(bary), w)
                  : SymmetricOrbit(std::move(bary), std::move(tensor), w));
            }
            return out;
          }

        private:
          Geometry::Polytope::Type m_g;
          size_t m_d;
          size_t m_degree = 0;
          Configuration m_config;
          Eigen::Index m_nz = 0;
          Eigen::Index m_ne = 0;
          Eigen::Index m_no = 0;
          Real m_average = 1;
          bool m_degenerate = false;
          std::vector<std::vector<size_t>> m_basis;
          std::vector<Real> m_invNorm;
          std::vector<Eigen::Index> m_keep;
          std::vector<std::vector<size_t>> m_fullBasis;
          std::vector<Real> m_fullInvNorm;
          Eigen::Index m_full = 0;

          mutable Math::Matrix<Real> m_A;
          mutable std::vector<Math::Matrix<Real>> m_dA;
          mutable Math::Vector<Real> m_b;
          mutable Math::Vector<Real> m_omega;
          mutable Eigen::JacobiSVD<Math::Matrix<Real>> m_svd;
      };

      /**
       * @brief Largest relative disagreement between the analytic Jacobian of
       * Eliminated and a central difference of its residual.
       *
       * An analytic Jacobian that is quietly wrong does not announce itself:
       * the search still runs, converges more slowly, and returns rules a few
       * digits short. This is the check that would catch it.
       */
      static Real verifyJacobian(Geometry::Polytope::Type g, size_t degree,
        const Configuration& config, unsigned seed = 20260101u)
      {
        const Eliminated problem(g, degree, config);
        const Eigen::Index nz = problem.getVariableCount();
        std::mt19937 rng(seed);
        std::normal_distribution<Real> gauss(0, 1);

        Real worst = 0;
        for (size_t trial = 0; trial < 8; ++trial)
        {
          Math::Vector<Real> z(nz);
          for (Eigen::Index i = 0; i < nz; ++i)
            z(i) = gauss(rng);

          Math::Vector<Real> r;
          Math::Matrix<Real> J;
          problem.evaluate(z, r, &J);

          for (Eigen::Index l = 0; l < nz; ++l)
          {
            // A single step size proves nothing: too large and truncation
            // dominates, too small and roundoff does, and at high degree the
            // residual curves sharply enough that the usable window moves.
            // Sweeping and taking the best agreement compares the analytic
            // derivative against the difference quotient at its own optimum.
            Real closest = std::numeric_limits<Real>::max();
            for (const Real relative : {1e-3, 1e-4, 1e-5, 1e-6, 1e-7})
            {
              const Real h = relative * std::max(std::abs(z(l)), Real(1));
              Math::Vector<Real> zp = z, zm = z, rp, rm;
              zp(l) += h;
              zm(l) -= h;
              problem.evaluate(zp, rp, nullptr);
              problem.evaluate(zm, rm, nullptr);
              const Math::Vector<Real> fd = (rp - rm) / (2 * h);
              const Real scale = std::max(
                {fd.cwiseAbs().maxCoeff(), J.col(l).cwiseAbs().maxCoeff(), Real(1e-8)});
              closest = std::min(closest, (fd - J.col(l)).cwiseAbs().maxCoeff() / scale);
            }
            worst = std::max(worst, closest);
          }
        }
        return worst;
      }

      /**
       * @brief Solves a configuration with the weights eliminated.
       *
       * A Levenberg--Marquardt loop over Eliminated, restarted from
       * pseudo-random points until it finds a rule that is both exact and
       * positive. The generator is seeded deterministically, so a given
       * (geometry, degree, configuration) always yields the same rule:
       * coefficients that are regenerated must not drift.
       */
      static Result solveEliminated(Geometry::Polytope::Type g, size_t degree,
        const Configuration& config, size_t maxRestarts = 256, Real tolerance = 1e-12,
        size_t maxIterations = 600, unsigned seed = 20260101u)
      {
        const Eliminated problem(g, degree, config);
        const Eigen::Index nz = problem.getVariableCount();

        std::mt19937 rng(seed);
        // Search variables live on the whole line; this spread covers the unit
        // interval of stick fractions without crowding either end.
        std::normal_distribution<Real> gauss(0, 1.5);

        Result best;
        for (size_t restart = 0; restart < maxRestarts; ++restart)
        {
          Math::Vector<Real> z(nz);
          for (Eigen::Index i = 0; i < nz; ++i)
            z(i) = problem.isDegenerate() ? Real(0) : gauss(rng);

          Real lambda = 1e-3;
          Math::Vector<Real> r;
          Math::Matrix<Real> J;
          problem.evaluate(z, r, &J);
          Real cost = r.squaredNorm();
          for (size_t it = 0; it < maxIterations && std::sqrt(cost) > tolerance; ++it)
          {
            const Math::Matrix<Real> H =
              J.transpose() * J + lambda * Math::Matrix<Real>::Identity(nz, nz);
            const Math::Vector<Real> zn = z + H.ldlt().solve(-J.transpose() * r);
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

          Math::Vector<Real> settled;
          problem.evaluate(z, settled, nullptr);

          Result candidate;
          // Judged on the whole moment system rather than the pruned one, so
          // a mode dropped in error cannot be mistaken for a rule.
          candidate.residual = std::max(std::sqrt(cost), problem.fullResidual(z));
          candidate.converged = candidate.residual < tolerance;
          candidate.restarts = restart + 1;
          candidate.orbits = problem.toOrbits(z);
          candidate.admissible = problem.getRelativeFloor() > s_floor;

          if (candidate.residual < best.residual)
            best = candidate;
          if (candidate.converged && candidate.admissible)
            return candidate;
          if (problem.isDegenerate())
            break; // nothing left to re-seed
        }
        return best;
      }

      /**
       * @brief Searches for a rule of strength @p degree with the orbit
       * classes @p config.
       *
       * Restarts are drawn from a deterministically seeded generator, so the
       * same request always produces the same rule.
       */
      static Result solve(Geometry::Polytope::Type g, size_t degree,
        const Configuration& config, size_t maxRestarts = 256, Real tolerance = 1e-13,
        size_t maxIterations = 200, unsigned seed = 20260101u)
      {
        const Geometry::Polytope::Traits traits(g);
        const size_t d = traits.getDimension();
        const Eigen::Index n = static_cast<Eigen::Index>(unknownCount(config));
        const Real measure = referenceMoment(g, std::vector<size_t>(d, 0));

        const MomentData moments(g, d, degree);

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
            for (const auto& c : config)
              totalPoints += classSize(c);
            Eigen::Index cursor = 0;
            for (const auto& c : config)
            {
              const auto& pattern = c.barycentric;
              for (size_t i = 0; i + 1 < pattern.size(); ++i)
                x(cursor++) = bary(rng) / static_cast<Real>(pattern.size());
              if (c.tensor && !c.midPlane)
                x(cursor++) = bary(rng);
              x(cursor++) = std::sqrt(measure / static_cast<Real>(totalPoints));
            }
          }

          Real lambda = 1e-3;
          Math::Vector<Real> r = residual(g, config, x, moments, d);
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
              J.col(j) = (residual(g, config, xp, moments, d) -
                           residual(g, config, xm, moments, d)) /
                (2 * h);
            }
            const Math::Matrix<Real> H =
              J.transpose() * J + lambda * Math::Matrix<Real>::Identity(n, n);
            const Math::Vector<Real> step = H.ldlt().solve(-J.transpose() * r);
            const Math::Vector<Real> xn = x + step;
            const Math::Vector<Real> rn = residual(g, config, xn, moments, d);
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

      /// @brief The orbit classes available on a @f$ d @f$-simplex: every
      /// partition of @f$ d+1 @f$. On the triangle these are @f$ \{3\} @f$,
      /// @f$ \{2,1\} @f$, @f$ \{1,1,1\} @f$; on the tetrahedron
      /// @f$ \{4\} @f$, @f$ \{3,1\} @f$, @f$ \{2,2\} @f$, @f$ \{2,1,1\} @f$
      /// and @f$ \{1,1,1,1\} @f$.
      static std::vector<Pattern> patternsFor(Geometry::Polytope::Type g)
      {
        const size_t n = (g == Geometry::Polytope::Type::Wedge)
          ? 3 // the base simplex of the reference wedge is the triangle
          : Geometry::Polytope::Traits(g).getVertexCount();
        std::vector<Pattern> out;
        Pattern current;
        const std::function<void(size_t, size_t)> rec = [&](size_t remaining,
                                                          size_t largest) {
          if (remaining == 0)
          {
            out.push_back(current);
            return;
          }
          for (size_t part = std::min(remaining, largest); part >= 1; --part)
          {
            current.push_back(part);
            rec(remaining - part, part);
            current.pop_back();
          }
        };
        rec(n, n);
        return out;
      }

      /// @brief Number of points an orbit class contributes: the multinomial
      /// coefficient of its multiplicity pattern.
      static size_t patternSize(const Pattern& pattern)
      {
        SymmetricOrbit::Barycentric probe;
        for (size_t i = 0; i < pattern.size(); ++i)
          for (size_t m = 0; m < pattern[i]; ++m)
            probe.push_back(static_cast<Real>(i + 1));
        return SymmetricOrbit(probe, 0).getSize();
      }

      /// @brief Total number of points in a configuration.
      /// @brief Points contributed by one orbit class.
      static size_t classSize(const OrbitClass& c)
      {
        return patternSize(c.barycentric) * (c.tensor && !c.midPlane ? 2 : 1);
      }

      /// @brief Total number of points in a configuration.
      static size_t configurationSize(const Configuration& config)
      {
        size_t n = 0;
        for (const auto& c : config)
          n += classSize(c);
        return n;
      }

      /**
       * @note The Koornwinder-Dubiner basis, which lifted node elimination
       * from degree 14 to 20, was tried here and made this solver worse: the
       * triangle ceiling fell from 9 to 8, and degree 3 stopped converging
       * within 256 restarts where it had succeeded before.
       *
       * The two problems have opposite shape. Node elimination carries the
       * node coordinates as unknowns, so its moment system is heavily
       * underdetermined and its difficulty is conditioning; an orthogonal
       * basis addresses exactly that. Here the unknowns are a handful of orbit
       * parameters, so the system is strongly overdetermined -- ten equations
       * against four unknowns at triangle degree 3 -- and the difficulty is
       * finding the basin of an isolated root. Changing the basis rescales the
       * residual landscape without making that search easier, and empirically
       * makes it harder.
       *
       * Recorded so the experiment is not repeated.
       */

      /**
       * @brief Searches orbit configurations for the cheapest admissible rule
       * of strength @p degree.
       *
       * Configurations are enumerated as multisets of orbit classes, ordered
       * by total point count, and the first that converges with positive
       * weights and interior points is returned. The centroid class appears at
       * most once, since a second copy is the same point.
       *
       * This is the search Xiao and Gimbutas @cite xiao2010numerical and
       * Witherden and Vincent @cite witherden2015identification each perform;
       * the point counts it finds are therefore the figure of merit to compare
       * against theirs.
       */
      static Result search(Geometry::Polytope::Type g, size_t degree,
        size_t maxPoints = 32, size_t maxRestarts = 96, Configuration* found = nullptr,
        size_t maxOrbits = 6, bool eliminate = true, Real tolerance = 1e-12)
      {
        const bool product = (g == Geometry::Polytope::Type::Wedge);
        std::vector<OrbitClass> patterns;
        for (const auto& p : patternsFor(g))
        {
          if (!product)
          {
            patterns.emplace_back(p);
          }
          else
          {
            patterns.emplace_back(p, true); // mid-plane
            patterns.emplace_back(p, false); // reflected pair
          }
        }
        std::vector<Configuration> candidates;
        Configuration current;
        // Published symmetric rules use a handful of orbits; without a cap the
        // multiset enumeration explodes long before any of it is solved, which
        // on the wedge means no candidate is ever tried.
        const std::function<void(size_t, size_t)> rec = [&](size_t start, size_t points) {
          if (!current.empty())
            candidates.push_back(current);
          if (current.size() >= maxOrbits)
            return;
          for (size_t i = start; i < patterns.size(); ++i)
          {
            const size_t sz = classSize(patterns[i]);
            if (points + sz > maxPoints)
              continue;
            // A class with no free parameters contributes the same points
            // however often it is repeated, so it may appear at most once.
            const bool centroid = (patterns[i].barycentric.size() == 1) &&
              (!patterns[i].tensor || patterns[i].midPlane);
            if (centroid && std::count(current.begin(), current.end(), patterns[i]))
              continue;
            current.push_back(patterns[i]);
            rec(centroid ? i + 1 : i, points + sz);
            current.pop_back();
          }
        };
        rec(0, 0);

        std::stable_sort(candidates.begin(), candidates.end(),
          [](const Configuration& a, const Configuration& b) {
            return configurationSize(a) < configurationSize(b);
          });

        Result best;
        for (const auto& config : candidates)
        {
          const auto r = eliminate
            ? solveEliminated(g, degree, config, maxRestarts, tolerance)
            : solve(g, degree, config, maxRestarts, tolerance);
          if (r.converged && r.admissible)
          {
            if (found)
              *found = config;
            return r;
          }
          if (r.residual < best.residual)
            best = r;
        }
        return best;
      }

      /// @brief A candidate pyramid rule: one orbit class per entry.
      using PyramidConfiguration = std::vector<SymmetricOrbit::PyramidClass>;

      /// @brief Points contributed by a pyramid configuration.
      static size_t pyramidConfigurationSize(const PyramidConfiguration& config)
      {
        size_t n = 0;
        for (const auto c : config)
          n += SymmetricOrbit::pyramidClassSize(c);
        return n;
      }

      /// @brief Unknowns of a pyramid configuration: the shape parameters of
      /// each orbit, its height, and its weight.
      static size_t pyramidUnknownCount(const PyramidConfiguration& config)
      {
        size_t n = 0;
        for (const auto c : config)
          n += SymmetricOrbit::pyramidClassParameters(c) + 2;
        return n;
      }

      /// @brief Expands a pyramid parameter vector into points and weights.
      static void pyramidPoints(const PyramidConfiguration& config,
        const Math::Vector<Real>& x, std::vector<Math::SpatialVector<Real>>& points,
        std::vector<Real>& weights)
      {
        points.clear();
        weights.clear();
        Eigen::Index cursor = 0;
        for (const auto c : config)
        {
          const size_t np = SymmetricOrbit::pyramidClassParameters(c);
          const Real alpha = (np >= 1) ? x(cursor++) : Real(0);
          const Real beta = (np >= 2) ? x(cursor++) : Real(0);
          const Real z = x(cursor++);
          const Real theta = x(cursor++);
          const Real w = theta * theta;
          for (auto& p : SymmetricOrbit::expandPyramid(c, alpha, beta, z))
          {
            points.push_back(std::move(p));
            weights.push_back(w);
          }
        }
      }

      /**
       * @brief Searches pyramid orbit configurations for an admissible rule.
       *
       * Orbits live in the centred coordinates of SymmetricOrbit::expandPyramid
       * and are sheared onto the reference element, whose apex sits over a base
       * corner. Weights enter as theta^2, and the shape parameters are confined
       * to (-1, 1) by construction, so every iterate is a rule of the right
       * shape even before it is exact.
       */
      static Result searchPyramid(size_t degree, size_t maxPoints = 40,
        size_t maxRestarts = 400, PyramidConfiguration* found = nullptr,
        size_t maxOrbits = 5, unsigned seed = 20260101u)
      {
        using PC = SymmetricOrbit::PyramidClass;
        const auto g = Geometry::Polytope::Type::Pyramid;
        const MomentData moments(g, 3, degree);

        std::vector<PyramidConfiguration> candidates;
        PyramidConfiguration current;
        const std::vector<PC> classes = {PC::Centre, PC::Axis, PC::Diagonal, PC::General};
        const std::function<void(size_t, size_t)> rec = [&](size_t start, size_t points) {
          if (!current.empty())
            candidates.push_back(current);
          if (current.size() >= maxOrbits)
            return;
          for (size_t i = start; i < classes.size(); ++i)
          {
            const size_t sz = SymmetricOrbit::pyramidClassSize(classes[i]);
            if (points + sz > maxPoints)
              continue;
            current.push_back(classes[i]);
            rec(i, points + sz);
            current.pop_back();
          }
        };
        rec(0, 0);
        std::stable_sort(candidates.begin(), candidates.end(),
          [](const PyramidConfiguration& a, const PyramidConfiguration& b) {
            return pyramidConfigurationSize(a) < pyramidConfigurationSize(b);
          });

        Result best;
        for (const auto& config : candidates)
        {
          const Eigen::Index n = static_cast<Eigen::Index>(pyramidUnknownCount(config));
          std::mt19937 rng(seed);
          std::uniform_real_distribution<Real> shape(0.05, 0.95);
          std::uniform_real_distribution<Real> height(0.05, 0.9);

          std::vector<Math::SpatialVector<Real>> pts;
          std::vector<Real> wts;
          const auto residualOf = [&](const Math::Vector<Real>& v) {
            pyramidPoints(config, v, pts, wts);
            Math::Vector<Real> r(static_cast<Eigen::Index>(moments.alphas.size()));
            for (size_t e = 0; e < moments.alphas.size(); ++e)
            {
              const auto& a = moments.alphas[e];
              Real sum = 0;
              for (size_t q = 0; q < pts.size(); ++q)
                sum += wts[q] * ipow(pts[q][0], a[0]) * ipow(pts[q][1], a[1]) *
                  ipow(pts[q][2], a[2]);
              r(static_cast<Eigen::Index>(e)) =
                (sum - moments.exact[e]) * moments.scale[e];
            }
            return r;
          };

          for (size_t restart = 0; restart < maxRestarts; ++restart)
          {
            Math::Vector<Real> x(n);
            {
              Eigen::Index cursor = 0;
              const size_t total = pyramidConfigurationSize(config);
              for (const auto c : config)
              {
                for (size_t k = 0; k < SymmetricOrbit::pyramidClassParameters(c); ++k)
                  x(cursor++) = shape(rng);
                x(cursor++) = height(rng);
                x(cursor++) = std::sqrt(Real(1) / (3 * static_cast<Real>(total)));
              }
            }

            Real lambda = 1e-3;
            Math::Vector<Real> r = residualOf(x);
            Real cost = r.squaredNorm();
            for (size_t it = 0; it < 200 && std::sqrt(cost) > 1e-13; ++it)
            {
              Math::Matrix<Real> J(r.size(), n);
              for (Eigen::Index j = 0; j < n; ++j)
              {
                const Real h = 1e-7 * std::max(std::abs(x(j)), Real(1));
                Math::Vector<Real> xp = x, xm = x;
                xp(j) += h;
                xm(j) -= h;
                J.col(j) = (residualOf(xp) - residualOf(xm)) / (2 * h);
              }
              const Math::Matrix<Real> H =
                J.transpose() * J + lambda * Math::Matrix<Real>::Identity(n, n);
              const Math::Vector<Real> xn = x + H.ldlt().solve(-J.transpose() * r);
              const Math::Vector<Real> rn = residualOf(xn);
              const Real cn = rn.squaredNorm();
              if (cn < cost)
              {
                x = xn;
                r = rn;
                cost = cn;
                lambda = std::max(lambda * Real(0.3), Real(1e-12));
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
            candidate.converged = candidate.residual < 1e-13;
            candidate.restarts = restart + 1;
            pyramidPoints(config, x, pts, wts);
            candidate.admissible = true;
            const auto& hs = Geometry::Polytope::Traits(g).getHalfSpace();
            for (size_t q = 0; q < pts.size() && candidate.admissible; ++q)
            {
              if (!(wts[q] > 0))
                candidate.admissible = false;
              Math::Vector<Real> pt(3);
              for (int k = 0; k < 3; ++k)
                pt(k) = pts[q][k];
              if ((hs.matrix * pt - hs.vector).maxCoeff() > -1e-12)
                candidate.admissible = false;
            }
            if (candidate.residual < best.residual)
            {
              best = candidate;
              best.orbits.clear();
              for (size_t q = 0; q < pts.size(); ++q)
                best.orbits.emplace_back(
                  SymmetricOrbit::Barycentric{pts[q][0], pts[q][1], pts[q][2]}, wts[q]);
            }
            if (candidate.converged && candidate.admissible)
            {
              if (found)
                *found = config;
              Result out = candidate;
              out.orbits.clear();
              for (size_t q = 0; q < pts.size(); ++q)
                out.orbits.emplace_back(
                  SymmetricOrbit::Barycentric{pts[q][0], pts[q][1], pts[q][2]}, wts[q]);
              return out;
            }
          }
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
          for (const auto t : orbit.getTensor())
            if (t < tol || t > Real(1) - tol)
              return false;
        }
        (void)g;
        return true;
      }
  };
}

#endif
