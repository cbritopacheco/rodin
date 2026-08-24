/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Tests CollapsedBasis independently of any rule generator.
 *
 * The rule generators state their moment equations in this basis, so an error
 * here does not surface as a failure: the search converges happily against the
 * wrong equations and returns a rule that is confidently incorrect. These
 * tests therefore check the basis on its own terms --- that it is orthonormal,
 * that it spans what it claims to span, that its gradient is its own
 * derivative, and that the rule used to measure it is itself exact.
 */
#include <gtest/gtest.h>

#include <cmath>
#include <numeric>
#include <random>

#include "Rodin/QF/CollapsedBasis.h"

#include "QuadratureInvariants.h"

using namespace Rodin;
using namespace Rodin::QF;
using namespace Rodin::Geometry;
using namespace Rodin::Tests::QF;

namespace
{
  /// @brief Every element the basis claims to cover, with its known measure.
  struct Element
  {
      const char* name;
      Polytope::Type type;
      Real measure;
  };

  const std::vector<Element>& elements()
  {
    static const std::vector<Element> all = {
      {"segment", Polytope::Type::Segment, 1.0},
      {"triangle", Polytope::Type::Triangle, 0.5},
      {"quadrilateral", Polytope::Type::Quadrilateral, 1.0},
      {"tetrahedron", Polytope::Type::Tetrahedron, 1.0 / 6.0},
      {"wedge", Polytope::Type::Wedge, 0.5},
      {"pyramid", Polytope::Type::Pyramid, 1.0 / 3.0},
      {"hexahedron", Polytope::Type::Hexahedron, 1.0},
    };
    return all;
  }

  /// @brief Monomial exponents of total degree at most @p degree.
  std::vector<std::vector<size_t>> monomials(size_t d, size_t degree)
  {
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

  Real monomialAt(const std::vector<size_t>& alpha, const Math::SpatialVector<Real>& x)
  {
    Real out = 1;
    for (size_t k = 0; k < alpha.size(); ++k)
      out *= CollapsedBasis::ipow(x[static_cast<Eigen::Index>(k)], alpha[k]);
    return out;
  }
}

/**
 * @brief The measuring rule reproduces the measure and integrates monomials
 * exactly.
 *
 * Everything else here is measured with this rule, so it is checked first and
 * against a quantity known in closed form.
 */
TEST(CollapsedBasisTest, ExactRuleIsExactAndPositive)
{
  for (const auto& element : elements())
  {
    const size_t d = CollapsedBasis::axes(element.type).size();
    for (const size_t degree : {2u, 5u, 10u})
    {
      std::vector<Math::SpatialVector<Real>> points;
      std::vector<Real> weights;
      CollapsedBasis::exactRule(element.type, degree + 1, points, weights);

      const Real volume = std::accumulate(weights.begin(), weights.end(), Real(0));
      EXPECT_NEAR(volume, element.measure, 1e-13)
        << element.name << ", degree " << degree;

      for (const Real w : weights)
        EXPECT_GT(w, 0) << element.name << " weight should be positive";

      // Against the moments of the reference element, computed independently.
      for (const auto& alpha : monomials(d, degree))
      {
        Real got = 0;
        for (size_t q = 0; q < weights.size(); ++q)
          got += weights[q] * monomialAt(alpha, points[q]);
        const Real expected = exactMoment(
          element.type, alpha[0], (d > 1) ? alpha[1] : 0, (d > 2) ? alpha[2] : 0);
        EXPECT_NEAR(got, expected, 1e-13 * std::max(std::abs(expected), Real(1)))
          << element.name << ", monomial of degree " << degree;
      }
    }
  }
}

/**
 * @brief The basis is orthonormal once scaled, and every non-constant mode
 * integrates to zero.
 *
 * The eliminated rule solver takes the right hand side of its moment system to
 * be zero away from the constant mode. Were that false the search would
 * converge to a small residual against the wrong equations, which looks
 * exactly like success.
 */
TEST(CollapsedBasisTest, BasisIsOrthonormal)
{
  for (const auto& element : elements())
  {
    for (const size_t degree : {2u, 6u, 12u})
    {
      std::vector<Math::SpatialVector<Real>> points;
      std::vector<Real> weights;
      CollapsedBasis::exactRule(element.type, degree + 2, points, weights);
      const Real volume = std::accumulate(weights.begin(), weights.end(), Real(0));

      const auto index = CollapsedBasis::indices(element.type, degree);
      std::vector<std::vector<Real>> value(index.size());
      std::vector<Real> inverseNorm(index.size());
      for (size_t a = 0; a < index.size(); ++a)
      {
        value[a].resize(weights.size());
        Real square = 0;
        for (size_t q = 0; q < weights.size(); ++q)
        {
          value[a][q] = CollapsedBasis::evaluate(element.type, index[a], points[q]);
          square += weights[q] * value[a][q] * value[a][q];
        }
        ASSERT_GT(square, 0) << element.name << " mode " << a << " vanishes";
        inverseNorm[a] = 1 / std::sqrt(square);
      }

      for (size_t a = 0; a < index.size(); ++a)
      {
        Real mean = 0;
        for (size_t q = 0; q < weights.size(); ++q)
          mean += weights[q] * value[a][q];
        if (a == 0)
          EXPECT_NEAR(mean, volume, 1e-12) << element.name;
        else
          EXPECT_NEAR(mean, 0, 1e-12)
            << element.name << " mode " << a << " should integrate to zero";

        for (size_t b = 0; b <= a; ++b)
        {
          Real gram = 0;
          for (size_t q = 0; q < weights.size(); ++q)
            gram += weights[q] * value[a][q] * value[b][q];
          gram *= inverseNorm[a] * inverseNorm[b];
          EXPECT_NEAR(gram, (a == b) ? 1 : 0, 1e-11)
            << element.name << " Gram entry (" << a << ", " << b << ")";
        }
      }
    }
  }
}

/**
 * @brief The basis spans the polynomials of its stated degree.
 *
 * Orthonormality alone does not establish this: an orthonormal set can be
 * perfectly well conditioned and still span the wrong space, in which case the
 * moment conditions it expresses are not the ones intended. The check fits
 * every monomial of the degree in the basis and demands the fit be exact at
 * points that took no part in it.
 */
TEST(CollapsedBasisTest, BasisSpansPolynomialsOfItsDegree)
{
  std::mt19937 rng(20260824u);
  std::uniform_real_distribution<Real> uniform(0, 1);

  for (const auto& element : elements())
  {
    const size_t d = CollapsedBasis::axes(element.type).size();
    for (const size_t degree : {3u, 7u})
    {
      std::vector<Math::SpatialVector<Real>> points;
      std::vector<Real> weights;
      CollapsedBasis::exactRule(element.type, degree + 2, points, weights);
      const auto index = CollapsedBasis::indices(element.type, degree);
      const auto alphas = monomials(d, degree);
      ASSERT_EQ(index.size(), alphas.size())
        << element.name << ": basis and monomial counts must agree";

      // Coefficients by orthogonal projection, which is exact when the basis
      // spans and visibly wrong when it does not.
      for (const auto& alpha : alphas)
      {
        std::vector<Real> coefficient(index.size());
        for (size_t a = 0; a < index.size(); ++a)
        {
          Real inner = 0, square = 0;
          for (size_t q = 0; q < weights.size(); ++q)
          {
            const Real psi = CollapsedBasis::evaluate(element.type, index[a], points[q]);
            inner += weights[q] * psi * monomialAt(alpha, points[q]);
            square += weights[q] * psi * psi;
          }
          coefficient[a] = inner / square;
        }

        // Evaluated away from the nodes that produced the coefficients.
        for (size_t trial = 0; trial < 16; ++trial)
        {
          Math::SpatialVector<Real> x;
          x.resize(static_cast<Eigen::Index>(d));
          // Rejection onto the element keeps the sample strictly interior.
          for (size_t attempt = 0;; ++attempt)
          {
            for (size_t k = 0; k < d; ++k)
              x[static_cast<Eigen::Index>(k)] = uniform(rng);
            if (isInside(element.type, x, 0.0))
              break;
            ASSERT_LT(attempt, 10000u) << element.name;
          }

          Real fitted = 0;
          for (size_t a = 0; a < index.size(); ++a)
            fitted +=
              coefficient[a] * CollapsedBasis::evaluate(element.type, index[a], x);
          EXPECT_NEAR(fitted, monomialAt(alpha, x), 1e-10)
            << element.name << ": basis fails to represent a monomial of degree "
            << degree;
        }
      }
    }
  }
}

/**
 * @brief The analytic gradient is the derivative of the value.
 *
 * The step is swept rather than fixed: a single step proves nothing, since too
 * large a step is dominated by truncation and too small a one by roundoff, and
 * the usable window moves with the degree.
 */
TEST(CollapsedBasisTest, GradientMatchesDifferenceQuotient)
{
  std::mt19937 rng(20260824u);
  std::uniform_real_distribution<Real> uniform(0.05, 0.95);

  for (const auto& element : elements())
  {
    const size_t d = CollapsedBasis::axes(element.type).size();
    for (const size_t degree : {3u, 8u})
    {
      const auto index = CollapsedBasis::indices(element.type, degree);
      for (size_t trial = 0; trial < 24; ++trial)
      {
        Math::SpatialVector<Real> x;
        x.resize(static_cast<Eigen::Index>(d));
        for (size_t attempt = 0;; ++attempt)
        {
          for (size_t k = 0; k < d; ++k)
            x[static_cast<Eigen::Index>(k)] = uniform(rng);
          // Kept clear of the boundary, where the collapse degenerates and a
          // difference quotient straddles the singularity.
          if (isInside(element.type, x, 0.0))
            break;
          ASSERT_LT(attempt, 10000u) << element.name;
        }

        const auto& mode = index[trial % index.size()];
        const auto analytic = CollapsedBasis::gradient(element.type, mode, x);
        for (size_t k = 0; k < d; ++k)
        {
          Real closest = std::numeric_limits<Real>::max();
          for (const Real step : {1e-3, 1e-4, 1e-5, 1e-6, 1e-7})
          {
            Math::SpatialVector<Real> plus = x, minus = x;
            plus[static_cast<Eigen::Index>(k)] += step;
            minus[static_cast<Eigen::Index>(k)] -= step;
            if (!isInside(element.type, plus, 0.0) || !isInside(element.type, minus, 0.0))
              continue;
            const Real difference =
              (CollapsedBasis::evaluate(element.type, mode, plus) -
                CollapsedBasis::evaluate(element.type, mode, minus)) /
              (2 * step);
            const Real scale = std::max({std::abs(difference),
              std::abs(analytic[static_cast<Eigen::Index>(k)]), Real(1e-6)});
            closest = std::min(closest,
              std::abs(difference - analytic[static_cast<Eigen::Index>(k)]) / scale);
          }
          if (closest == std::numeric_limits<Real>::max())
            continue;
          EXPECT_LT(closest, 1e-5)
            << element.name << ", direction " << k << ", degree " << degree;
        }
      }
    }
  }
}

/**
 * @brief The tabulated evaluation agrees with the direct one.
 *
 * The tabulated path exists only for speed, and the two must therefore be the
 * same function. If they drift apart the generators quietly change what they
 * are solving, since they use the fast path and every check here uses the slow
 * one.
 */
TEST(CollapsedBasisTest, TabulatedEvaluationMatchesDirect)
{
  std::mt19937 rng(20260824u);
  std::uniform_real_distribution<Real> uniform(0.02, 0.98);

  for (const auto& element : elements())
  {
    const size_t d = CollapsedBasis::axes(element.type).size();
    for (const size_t degree : {3u, 9u})
    {
      const auto index = CollapsedBasis::indices(element.type, degree);
      for (size_t trial = 0; trial < 12; ++trial)
      {
        Math::SpatialVector<Real> x;
        x.resize(static_cast<Eigen::Index>(d));
        for (size_t attempt = 0;; ++attempt)
        {
          for (size_t k = 0; k < d; ++k)
            x[static_cast<Eigen::Index>(k)] = uniform(rng);
          if (isInside(element.type, x, 0.0))
            break;
          ASSERT_LT(attempt, 10000u) << element.name;
        }

        const auto table = CollapsedBasis::tabulate(element.type, degree, x, true);
        for (const auto& mode : index)
        {
          const Real direct = CollapsedBasis::evaluate(element.type, mode, x);
          const Real tabulated = CollapsedBasis::evaluate(element.type, table, mode);
          EXPECT_NEAR(tabulated, direct, 1e-11 * std::max(std::abs(direct), Real(1)))
            << element.name << ": tabulated value differs";

          const auto directGradient = CollapsedBasis::gradient(element.type, mode, x);
          const auto tabulatedGradient =
            CollapsedBasis::gradient(element.type, table, mode);
          for (size_t k = 0; k < d; ++k)
          {
            const Eigen::Index kk = static_cast<Eigen::Index>(k);
            EXPECT_NEAR(tabulatedGradient[kk], directGradient[kk],
              1e-10 * std::max(std::abs(directGradient[kk]), Real(1)))
              << element.name << ": tabulated gradient differs in direction " << k;
          }
        }
      }
    }
  }
}
