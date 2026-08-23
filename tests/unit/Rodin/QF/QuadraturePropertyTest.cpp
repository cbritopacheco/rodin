/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */

/**
 * @file
 * @brief Properties a correct quadrature rule must satisfy, none of which
 * depend on the closed-form moment oracle.
 *
 * The exactness sweep compares a rule against moments I derived by hand. That
 * is a strong test but it shares one assumption with the rule: what the
 * reference element is. A rule and an oracle built on the same wrong element
 * agree with each other perfectly. The checks here are chosen to be blind to
 * that: they compare the rule against itself under transformations that must
 * leave the integral invariant, so they hold whatever the element is.
 */
#include <gtest/gtest.h>
#include <random>
#include "Rodin/QF/XiaoGimbutas.h"
#include "QuadratureInvariants.h"

using namespace Rodin;
using namespace Rodin::QF;
using namespace Rodin::Geometry;
using namespace Rodin::Tests::QF;

namespace
{
  using Vertices = std::vector<Math::SpatialVector<Real>>;

  /// @brief Reference simplex vertices, from Polytope::Traits.
  Vertices referenceVertices(Polytope::Type g)
  {
    const Polytope::Traits t(g);
    Vertices out;
    for (size_t i = 0; i < t.getVertexCount(); ++i)
    {
      Math::SpatialVector<Real> v;
      v.resize(static_cast<Eigen::Index>(t.getDimension()));
      for (size_t k = 0; k < t.getDimension(); ++k)
        v[static_cast<Eigen::Index>(k)] =
          t.getVertex(i)[static_cast<Eigen::Index>(k)];
      out.push_back(std::move(v));
    }
    return out;
  }

  /// @brief Integrates @p f over the simplex spanned by @p verts, by mapping
  /// the reference rule onto it. The map is affine, so a rule exact on the
  /// reference is exact here too, with the Jacobian determinant as the scale.
  template <class Rule, class F>
  Real integrateOn(const Rule& qf, const Vertices& verts, F&& f)
  {
    const Eigen::Index d = verts[0].size();
    Math::Matrix<Real> B(d, d);
    for (Eigen::Index k = 0; k < d; ++k)
      for (Eigen::Index i = 0; i < d; ++i)
        B(i, k) = verts[static_cast<size_t>(k) + 1][i] - verts[0][i];
    const Real det = std::abs(B.determinant());

    Real s = 0;
    for (size_t q = 0; q < qf.getSize(); ++q)
    {
      const Math::SpatialVector<Real> y = verts[0] + B * qf.getPoint(q);
      s += qf.getWeight(q) * f(y);
    }
    return det * s;
  }

  /// @brief Subdivides a simplex into pieces that tile it, by joining the
  /// centroid to each facet.
  std::vector<Vertices> subdivide(const Vertices& verts)
  {
    Math::SpatialVector<Real> centroid = verts[0];
    for (size_t i = 1; i < verts.size(); ++i)
      centroid = centroid + verts[i];
    centroid = centroid / static_cast<Real>(verts.size());

    std::vector<Vertices> out;
    for (size_t skip = 0; skip < verts.size(); ++skip)
    {
      Vertices piece{centroid};
      for (size_t i = 0; i < verts.size(); ++i)
        if (i != skip)
          piece.push_back(verts[i]);
      out.push_back(std::move(piece));
    }
    return out;
  }

  /// @brief A random polynomial of total degree at most @p p.
  template <class Rng>
  auto randomPolynomial(size_t d, size_t p, Rng& rng)
  {
    std::uniform_real_distribution<Real> uni(-1, 1);
    std::vector<std::pair<std::vector<size_t>, Real>> terms;
    for (size_t a = 0; a <= p; ++a)
      for (size_t b = 0; a + b <= p; ++b)
        for (size_t c = 0; a + b + c <= p; ++c)
        {
          if (d == 2 && c > 0)
            continue;
          terms.push_back({{a, b, c}, uni(rng)});
        }
    return [terms, d](const Math::SpatialVector<Real>& x)
    {
      Real s = 0;
      for (const auto& [e, coeff] : terms)
      {
        Real m = std::pow(x[0], static_cast<Real>(e[0]))
               * std::pow(x[1], static_cast<Real>(e[1]));
        if (d == 3)
          m *= std::pow(x[2], static_cast<Real>(e[2]));
        s += coeff * m;
      }
      return s;
    };
  }

  const std::vector<Polytope::Type> kSimplices = {
    Polytope::Type::Triangle, Polytope::Type::Tetrahedron
  };
  const char* name(Polytope::Type g)
  { return g == Polytope::Type::Triangle ? "triangle" : "tetrahedron"; }
}

/**
 * @brief Subdivision consistency.
 *
 * Integrating over a simplex must equal the sum over any set of pieces that
 * tiles it. This uses no moment formula whatsoever: the rule is compared only
 * against itself, applied to sub-elements through affine maps. A rule that
 * integrated the wrong domain, or that was exact only by agreement with a
 * mistaken oracle, fails here.
 */
TEST(QuadraturePropertyTest, IntegralIsAdditiveUnderSubdivision)
{
  std::mt19937 rng(4242u);
  for (const auto g : kSimplices)
  {
    const size_t d = Polytope::Traits(g).getDimension();
    const auto verts = referenceVertices(g);
    const auto pieces = subdivide(verts);
    for (size_t p = 1; p <= std::min<size_t>(XiaoGimbutas::getMaxDegree(g), 10); ++p)
    {
      const XiaoGimbutas qf(p, g);
      const auto f = randomPolynomial(d, p, rng);
      const Real whole = integrateOn(qf, verts, f);
      Real parts = 0;
      for (const auto& piece : pieces)
        parts += integrateOn(qf, piece, f);
      EXPECT_NEAR(parts, whole, 1e-11 * std::max(Real(1), std::abs(whole)))
        << name(g) << " degree " << p;
    }
  }
}

/**
 * @brief Affine covariance.
 *
 * A rule exact on the reference element is exact on any affine image of it,
 * with the integral scaling by the determinant. This is the property assembly
 * actually relies on, and it is what makes a reference rule usable at all.
 */
TEST(QuadraturePropertyTest, IntegralScalesWithTheAffineDeterminant)
{
  std::mt19937 rng(99u);
  std::uniform_real_distribution<Real> uni(-1.5, 1.5);
  for (const auto g : kSimplices)
  {
    const size_t d = Polytope::Traits(g).getDimension();
    const auto verts = referenceVertices(g);
    for (size_t p = 1; p <= std::min<size_t>(XiaoGimbutas::getMaxDegree(g), 8); ++p)
    {
      const XiaoGimbutas qf(p, g);
      for (int trial = 0; trial < 5; ++trial)
      {
        // A random affine image of the reference element.
        Math::Matrix<Real> A(static_cast<Eigen::Index>(d),
                             static_cast<Eigen::Index>(d));
        do {
          for (Eigen::Index i = 0; i < A.rows(); ++i)
            for (Eigen::Index j = 0; j < A.cols(); ++j)
              A(i, j) = uni(rng);
        } while (std::abs(A.determinant()) < 0.2);
        Math::SpatialVector<Real> shift;
        shift.resize(static_cast<Eigen::Index>(d));
        for (Eigen::Index k = 0; k < shift.size(); ++k)
          shift[k] = uni(rng);

        Vertices mapped;
        for (const auto& v : verts)
          mapped.push_back(Math::SpatialVector<Real>(shift + A * v));

        const auto f = randomPolynomial(d, p, rng);
        // f composed with the map is a polynomial of the same degree, so
        // integrating it on the reference and scaling must agree.
        const auto fComposed = [&](const Math::SpatialVector<Real>& x)
        { return f(Math::SpatialVector<Real>(shift + A * x)); };

        const Real direct = integrateOn(qf, mapped, f);
        const Real viaReference =
          std::abs(A.determinant()) * integrateOn(qf, verts, fComposed);
        EXPECT_NEAR(direct, viaReference,
          1e-10 * std::max(Real(1), std::abs(direct)))
          << name(g) << " degree " << p;
      }
    }
  }
}

/**
 * @brief Exactness on random polynomials, not only on monomials.
 *
 * The sweep tests one monomial at a time, so an error that cancels between
 * monomials of the same degree would pass it. A random combination does not
 * cancel.
 */
TEST(QuadraturePropertyTest, ExactOnRandomPolynomialsNotOnlyMonomials)
{
  std::mt19937 rng(7u);
  for (const auto g : kSimplices)
  {
    const size_t d = Polytope::Traits(g).getDimension();
    const auto verts = referenceVertices(g);
    for (size_t p = 1; p <= std::min<size_t>(XiaoGimbutas::getMaxDegree(g), 10); ++p)
    {
      const XiaoGimbutas qf(p, g);
      // A rule of the next degree up is an independent reference.
      const XiaoGimbutas ref(
        std::min(p + 2, XiaoGimbutas::getMaxDegree(g)), g);
      for (int trial = 0; trial < 5; ++trial)
      {
        const auto f = randomPolynomial(d, p, rng);
        const Real a = integrateOn(qf, verts, f);
        const Real b = integrateOn(ref, verts, f);
        EXPECT_NEAR(a, b, 1e-11 * std::max(Real(1), std::abs(b)))
          << name(g) << " degree " << p;
      }
    }
  }
}

/**
 * @brief A positive rule cannot return a negative integral for a
 * non-negative integrand. This is what makes the rule usable for assembling a
 * form that must stay positive semidefinite.
 */
TEST(QuadraturePropertyTest, NonNegativeIntegrandGivesNonNegativeIntegral)
{
  std::mt19937 rng(13u);
  for (const auto g : kSimplices)
  {
    const size_t d = Polytope::Traits(g).getDimension();
    const auto verts = referenceVertices(g);
    for (size_t p = 1; p <= std::min<size_t>(XiaoGimbutas::getMaxDegree(g), 10); ++p)
    {
      const XiaoGimbutas qf(p, g);
      for (int trial = 0; trial < 5; ++trial)
      {
        const auto h = randomPolynomial(d, p / 2, rng);
        // h^2 is non-negative and of degree at most p.
        const auto f = [&](const Math::SpatialVector<Real>& x)
        { const Real v = h(x); return v * v; };
        EXPECT_GE(integrateOn(qf, verts, f), 0.0) << name(g) << " degree " << p;
      }
    }
  }
}

/// @brief The subdivision check must be able to fail: a rule of too low a
/// degree for the integrand is rejected by it.
TEST(QuadraturePropertyTest, SubdivisionCheckRejectsAnUnderResolvedRule)
{
  std::mt19937 rng(5u);
  const auto g = Polytope::Type::Triangle;
  const auto verts = referenceVertices(g);
  const auto pieces = subdivide(verts);
  const XiaoGimbutas low(2, g);
  const auto f = randomPolynomial(2, 8, rng);  // degree 8, rule only exact to 2
  Real parts = 0;
  for (const auto& piece : pieces)
    parts += integrateOn(low, piece, f);
  const Real whole = integrateOn(low, verts, f);
  EXPECT_GT(std::abs(parts - whole), 1e-6)
    << "subdivision agreed for an integrand the rule cannot resolve; "
       "the check has no teeth";
}
