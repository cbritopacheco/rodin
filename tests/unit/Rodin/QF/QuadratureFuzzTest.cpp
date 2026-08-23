/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */

/**
 * @file
 * @brief Randomised checks on the shipped quadrature rules.
 *
 * The property tests use a handful of fixed configurations each. A rule can
 * satisfy those and still fail on inputs nobody thought to write down, so the
 * same invariants are exercised here over many random ones: random
 * subdivisions rather than the centroid split, random affine images rather
 * than five, random polynomials rather than a fixed list.
 *
 * Every case is drawn from a fixed seed, so a failure is reproducible and a
 * regression cannot appear or vanish between runs.
 */
#include <gtest/gtest.h>
#include <random>
#include "Rodin/QF/XiaoGimbutas.h"
#include "Rodin/QF/WitherdenVincent.h"
#include "QuadratureInvariants.h"

using namespace Rodin;
using namespace Rodin::QF;
using namespace Rodin::Geometry;
using namespace Rodin::Tests::QF;

namespace
{
  using Vertices = std::vector<Math::SpatialVector<Real>>;

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
      s += qf.getWeight(q) * f(Math::SpatialVector<Real>(verts[0] + B * qf.getPoint(q)));
    return det * s;
  }

  template <class Rng>
  auto randomPolynomial(size_t d, size_t p, Rng& rng)
  {
    std::uniform_real_distribution<Real> uni(-1, 1);
    std::vector<std::pair<std::array<size_t, 3>, Real>> terms;
    for (size_t a = 0; a <= p; ++a)
      for (size_t b = 0; a + b <= p; ++b)
        for (size_t c = 0; a + b + c <= p; ++c)
        {
          if (d == 2 && c > 0) continue;
          terms.push_back({{a, b, c}, uni(rng)});
        }
    return [terms, d](const Math::SpatialVector<Real>& x)
    {
      Real s = 0;
      for (const auto& [e, k] : terms)
      {
        Real m = std::pow(x[0], (Real)e[0]) * std::pow(x[1], (Real)e[1]);
        if (d == 3) m *= std::pow(x[2], (Real)e[2]);
        s += k * m;
      }
      return s;
    };
  }

  /// @brief Splits a simplex by joining an arbitrary interior point to each
  /// facet. The centroid is one choice among infinitely many; any interior
  /// point tiles the element just as well, and a rule that only survives the
  /// symmetric split is not correct.
  std::vector<Vertices> subdivideAt(const Vertices& verts,
    const Math::SpatialVector<Real>& inner)
  {
    std::vector<Vertices> out;
    for (size_t skip = 0; skip < verts.size(); ++skip)
    {
      Vertices piece{inner};
      for (size_t i = 0; i < verts.size(); ++i)
        if (i != skip) piece.push_back(verts[i]);
      out.push_back(std::move(piece));
    }
    return out;
  }

  template <class Rng>
  Math::SpatialVector<Real> randomInterior(size_t d, Rng& rng)
  {
    std::uniform_real_distribution<Real> uni(0.12, 0.88);
    Math::SpatialVector<Real> x;
    x.resize(static_cast<Eigen::Index>(d));
    for (;;)
    {
      Real tot = 0;
      for (size_t k = 0; k < d; ++k)
      {
        x[static_cast<Eigen::Index>(k)] = uni(rng);
        tot += x[static_cast<Eigen::Index>(k)];
      }
      if (tot < 0.92) return x;
    }
  }

  const std::vector<Polytope::Type> kSimplices = {
    Polytope::Type::Triangle, Polytope::Type::Tetrahedron
  };
  const char* name(Polytope::Type g)
  { return g == Polytope::Type::Triangle ? "triangle" : "tetrahedron"; }
}

/// @brief Subdivision additivity over randomly placed split points.
TEST(QuadratureFuzzTest, AdditiveUnderRandomlyPlacedSubdivisions)
{
  std::mt19937 rng(20260823u);
  for (const auto g : kSimplices)
  {
    const size_t d = Polytope::Traits(g).getDimension();
    const auto verts = referenceVertices(g);
    for (size_t p = 1; p <= std::min<size_t>(XiaoGimbutas::getMaxDegree(g), 10); ++p)
    {
      const XiaoGimbutas qf(p, g);
      for (int trial = 0; trial < 25; ++trial)
      {
        const auto f = randomPolynomial(d, p, rng);
        const auto pieces = subdivideAt(verts, randomInterior(d, rng));
        Real parts = 0;
        for (const auto& piece : pieces)
          parts += integrateOn(qf, piece, f);
        const Real whole = integrateOn(qf, verts, f);
        ASSERT_NEAR(parts, whole, 1e-10 * std::max(Real(1), std::abs(whole)))
          << name(g) << " degree " << p << " trial " << trial;
      }
    }
  }
}

/// @brief Affine covariance over many random maps, including badly scaled and
/// orientation-reversing ones.
TEST(QuadratureFuzzTest, AffineCovariantUnderManyRandomMaps)
{
  std::mt19937 rng(555u);
  std::uniform_real_distribution<Real> uni(-2.5, 2.5);
  for (const auto g : kSimplices)
  {
    const size_t d = Polytope::Traits(g).getDimension();
    const auto verts = referenceVertices(g);
    for (size_t p = 1; p <= std::min<size_t>(XiaoGimbutas::getMaxDegree(g), 8); ++p)
    {
      const XiaoGimbutas qf(p, g);
      for (int trial = 0; trial < 30; ++trial)
      {
        Math::Matrix<Real> A((Eigen::Index)d, (Eigen::Index)d);
        do {
          for (Eigen::Index i = 0; i < A.rows(); ++i)
            for (Eigen::Index j = 0; j < A.cols(); ++j)
              A(i, j) = uni(rng);
        } while (std::abs(A.determinant()) < 0.05);
        Math::SpatialVector<Real> shift;
        shift.resize((Eigen::Index)d);
        for (Eigen::Index k = 0; k < shift.size(); ++k) shift[k] = uni(rng);

        Vertices mapped;
        for (const auto& v : verts)
          mapped.push_back(Math::SpatialVector<Real>(shift + A * v));
        const auto f = randomPolynomial(d, p, rng);
        const auto composed = [&](const Math::SpatialVector<Real>& x)
        { return f(Math::SpatialVector<Real>(shift + A * x)); };

        const Real direct = integrateOn(qf, mapped, f);
        const Real viaRef = std::abs(A.determinant()) * integrateOn(qf, verts, composed);
        ASSERT_NEAR(direct, viaRef, 1e-9 * std::max(Real(1), std::abs(direct)))
          << name(g) << " degree " << p << " trial " << trial
          << " det " << A.determinant();
      }
    }
  }
}

/// @brief The two families agree over many random polynomials.
TEST(QuadratureFuzzTest, FamiliesAgreeOverManyRandomPolynomials)
{
  std::mt19937 rng(31u);
  for (const auto g : kSimplices)
  {
    const size_t d = Polytope::Traits(g).getDimension();
    const auto verts = referenceVertices(g);
    const size_t top = std::min(XiaoGimbutas::getMaxDegree(g),
                                WitherdenVincent::getMaxDegree(g));
    for (size_t p = 1; p <= top; ++p)
    {
      const XiaoGimbutas xg(p, g);
      const WitherdenVincent wv(p, g);
      for (int trial = 0; trial < 40; ++trial)
      {
        const auto f = randomPolynomial(d, p, rng);
        const Real a = integrateOn(xg, verts, f), b = integrateOn(wv, verts, f);
        ASSERT_NEAR(a, b, 1e-10 * std::max(Real(1), std::abs(b)))
          << name(g) << " degree " << p << " trial " << trial;
      }
    }
  }
}

/**
 * @brief Randomised mutation: perturbing any single coefficient by any
 * appreciable amount must be caught.
 *
 * The fixed mutation test perturbs one weight by one part in 1e9. This varies
 * which coefficient, whether it is a weight or a coordinate, and by how much,
 * so the detection threshold is established across the table rather than at
 * one point of it.
 */
TEST(QuadratureFuzzTest, AnyPerturbedCoefficientIsDetected)
{
  std::mt19937 rng(6161u);
  std::uniform_real_distribution<Real> mag(1e-8, 1e-3);
  for (const auto g : kSimplices)
  {
    const size_t d = Polytope::Traits(g).getDimension();
    for (size_t p = 2; p <= std::min<size_t>(XiaoGimbutas::getMaxDegree(g), 10); ++p)
    {
      const XiaoGimbutas qf(p, g);
      ASSERT_LT(exactnessSweep(qf, g, p).worstRelativeError, 1e-11);
      std::uniform_int_distribution<size_t> node(0, qf.getSize() - 1);
      std::uniform_int_distribution<size_t> which(0, d);  // d = weight slot
      for (int trial = 0; trial < 20; ++trial)
      {
        const size_t q = node(rng), slot = which(rng);
        const Real delta = mag(rng);

        struct Perturbed
        {
          const XiaoGimbutas& base;
          size_t q, slot, d;
          Real delta;
          mutable Math::SpatialVector<Real> scratch;
          size_t getSize() const { return base.getSize(); }
          Real getWeight(size_t i) const
          { return base.getWeight(i) * ((i == q && slot == d) ? 1 + delta : 1); }
          const Math::SpatialVector<Real>& getPoint(size_t i) const
          {
            scratch = base.getPoint(i);
            if (i == q && slot < d)
              scratch[static_cast<Eigen::Index>(slot)] += delta;
            return scratch;
          }
        } bad{qf, q, slot, d, delta, {}};

        ASSERT_GT(exactnessSweep(bad, g, p).worstRelativeError, 1e-11)
          << name(g) << " degree " << p << " node " << q << " slot " << slot
          << " delta " << delta << " went undetected";
      }
    }
  }
}

/// @brief Weights and nodes survive many random reconstructions identically:
/// the tables are static data and constructing a rule must not mutate them.
TEST(QuadratureFuzzTest, RepeatedConstructionIsIdentical)
{
  std::mt19937 rng(77u);
  for (const auto g : kSimplices)
  {
    std::uniform_int_distribution<size_t> deg(1, XiaoGimbutas::getMaxDegree(g));
    for (int trial = 0; trial < 50; ++trial)
    {
      const size_t p = deg(rng);
      const XiaoGimbutas a(p, g), b(p, g);
      ASSERT_EQ(a.getSize(), b.getSize());
      for (size_t i = 0; i < a.getSize(); ++i)
      {
        ASSERT_EQ(a.getWeight(i), b.getWeight(i)) << name(g) << " degree " << p;
        for (Eigen::Index k = 0; k < a.getPoint(i).size(); ++k)
          ASSERT_EQ(a.getPoint(i)[k], b.getPoint(i)[k]);
      }
    }
  }
}
