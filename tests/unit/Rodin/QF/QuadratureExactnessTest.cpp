/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <string>
#include <limits>

#include "Rodin/Alert/Exception.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"
#include "Rodin/QF/Centroid.h"
#include "Rodin/QF/GrundmannMoller.h"
#include "Rodin/QF/GaussLegendre.h"
#include "Rodin/QF/TensorProduct.h"
#include "Rodin/QF/WitherdenVincent.h"
#include "Rodin/QF/XiaoGimbutas.h"

#include "QuadratureInvariants.h"

using namespace Rodin;
using namespace Rodin::QF;
using namespace Rodin::Geometry;
using namespace Rodin::Tests::QF;

namespace
{
  constexpr Real kExactTol = 1e-12;

  /**
   * @brief Lists every published reference geometry.
   *
   * Point carries the unit-mass rule on @f$\mathbb{R}^0@f$ and is included.
   * Taking types from the library also includes future published geometries.
   *
   * @returns All published reference geometry types.
   */
  std::vector<Polytope::Type> allGeometries()
  {
    std::vector<Polytope::Type> gs;
    for (const auto g : Polytope::Types)
      gs.push_back(g);
    return gs;
  }

  const char* name(Polytope::Type g)
  {
    switch (g)
    {
      case Polytope::Type::Point:
        return "Point";
      case Polytope::Type::Segment:
        return "Segment";
      case Polytope::Type::Triangle:
        return "Triangle";
      case Polytope::Type::Quadrilateral:
        return "Quadrilateral";
      case Polytope::Type::Tetrahedron:
        return "Tetrahedron";
      case Polytope::Type::Wedge:
        return "Wedge";
      case Polytope::Type::Pyramid:
        return "Pyramid";
      case Polytope::Type::Hexahedron:
        return "Hexahedron";
    }
    return "?";
  }
}

// ---------------------------------------------------------------------------
// Numerical: the rule integrates what it claims to integrate.
// ---------------------------------------------------------------------------

/**
 * @brief Checks every dispatched rule through its requested total degree.
 *
 * Closed-form moments provide an independent table-integrity gate on every
 * published reference geometry, including Point.
 */
TEST(QuadratureExactnessTest, DispatchedRuleIsExactToRequestedDegree)
{
  for (const auto g : allGeometries())
  {
    for (size_t order = 1; order <= 20; ++order)
    {
      const auto& qf = PolytopeQuadratureFormula::get(order, g);
      const auto report = exactnessSweep(qf, g, order);
      EXPECT_LT(report.worstRelativeError, kExactTol)
        << name(g) << " order " << order << ": worst monomial x^"
        << report.worstExponents[0] << " y^" << report.worstExponents[1] << " z^"
        << report.worstExponents[2] << " (" << report.monomialsTested
        << " monomials tested)";
    }
  }
}

/**
 * @brief Point rules have unit mass and the unique empty reference coordinate.
 *
 * A function on @f$\mathbb{R}^0@f$ is constant, so the one-point rule satisfies
 * @f$Q(c)=c@f$ at every requested order. Owned construction and canonical
 * reuse must both retain zero coordinate dimension, not a one-entry zero.
 */
TEST(QuadratureExactnessTest, PointRulesHaveUnitMassAndEmptyCoordinates)
{
  for (const size_t order : {0u, 1u, 15u, 16u, 17u, 51u})
  {
    SCOPED_TRACE(order);
    const auto owned = PolytopeQuadratureFormula::build(order, Polytope::Type::Point);
    const auto& canonical = PolytopeQuadratureFormula::get(order, Polytope::Type::Point);
    for (const auto* rule : {owned.get(), &canonical})
    {
      ASSERT_EQ(rule->getSize(), 1u);
      EXPECT_EQ(rule->getWeight(0), Real(1));
      EXPECT_EQ(rule->getPoint(0).size(), 0);
      EXPECT_TRUE(allPointsInside(*rule, Polytope::Type::Point));
    }
  }
}

/// @brief The arbitrary-order pyramid fallback remains exact beyond the table.
TEST(QuadratureExactnessTest, PyramidIsExactAtEvenOrders)
{
  for (size_t order = 12; order <= 20; order += 2)
  {
    const auto& qf = PolytopeQuadratureFormula::get(order, Polytope::Type::Pyramid);
    const auto report = exactnessSweep(qf, Polytope::Type::Pyramid, order);
    EXPECT_LT(report.worstRelativeError, kExactTol)
      << "Pyramid order " << order << ": worst at z^" << report.worstExponents[2];
  }
}

/// @brief The wedge tensor product is exact at odd orders.
TEST(QuadratureExactnessTest, WedgeIsExactAtOddOrders)
{
  for (size_t order = 1; order <= 19; order += 2)
  {
    const auto& qf = PolytopeQuadratureFormula::get(order, Polytope::Type::Wedge);
    const auto report = exactnessSweep(qf, Polytope::Type::Wedge, order);
    EXPECT_LT(report.worstRelativeError, kExactTol)
      << "Wedge order " << order << ": worst at x^" << report.worstExponents[0] << " y^"
      << report.worstExponents[1];
  }
}

/// @brief The dispatcher selects the documented family at every table boundary.
TEST(QuadratureExactnessTest, DispatchesDocumentedFamiliesAndFallbacks)
{
  using Type = Polytope::Type;

  EXPECT_NE(
    dynamic_cast<const Centroid*>(PolytopeQuadratureFormula::build(1, Type::Point).get()),
    nullptr);
  EXPECT_NE(dynamic_cast<const GaussLegendre*>(
              PolytopeQuadratureFormula::build(5, Type::Segment).get()),
    nullptr);

  EXPECT_NE(dynamic_cast<const XiaoGimbutas*>(
              PolytopeQuadratureFormula::build(50, Type::Triangle).get()),
    nullptr);
  EXPECT_NE(dynamic_cast<const GaussLegendre*>(
              PolytopeQuadratureFormula::build(51, Type::Triangle).get()),
    nullptr);
  EXPECT_NE(dynamic_cast<const TensorProduct*>(
              PolytopeQuadratureFormula::build(20, Type::Quadrilateral).get()),
    nullptr);

  EXPECT_NE(dynamic_cast<const XiaoGimbutas*>(
              PolytopeQuadratureFormula::build(15, Type::Tetrahedron).get()),
    nullptr);
  EXPECT_NE(dynamic_cast<const GaussLegendre*>(
              PolytopeQuadratureFormula::build(16, Type::Tetrahedron).get()),
    nullptr);
  EXPECT_NE(dynamic_cast<const TensorProduct*>(
              PolytopeQuadratureFormula::build(20, Type::Hexahedron).get()),
    nullptr);

  const auto wedge50 = PolytopeQuadratureFormula::build(50, Type::Wedge);
  const auto wedge51 = PolytopeQuadratureFormula::build(51, Type::Wedge);
  EXPECT_NE(dynamic_cast<const TensorProduct*>(wedge50.get()), nullptr);
  EXPECT_NE(dynamic_cast<const TensorProduct*>(wedge51.get()), nullptr);
  EXPECT_EQ(wedge50->getSize(), XiaoGimbutas(50, Type::Triangle).getSize() * 26);
  EXPECT_EQ(wedge51->getSize(), 26u * 26u * 26u);

  EXPECT_NE(dynamic_cast<const WitherdenVincent*>(
              PolytopeQuadratureFormula::build(10, Type::Pyramid).get()),
    nullptr);
  EXPECT_NE(dynamic_cast<const GaussLegendre*>(
              PolytopeQuadratureFormula::build(11, Type::Pyramid).get()),
    nullptr);
}

/**
 * @brief Checks positive conical-product exactness beyond the simplex tables.
 *
 * Full-degree moments, positive weights and reference-domain membership
 * are checked together at the first unsupported compact-table orders.
 * The degree-fifty-one line factor is checked independently as well: its
 * product with the certified triangle factor covers the wedge fallback.
 */
TEST(QuadratureExactnessTest, ConicalProductFallbacksAreExact)
{
  for (const auto& [g, order] : {std::pair{Polytope::Type::Segment, size_t{51}},
         std::pair{Polytope::Type::Triangle, size_t{51}},
         std::pair{Polytope::Type::Tetrahedron, size_t{16}}})
  {
    const auto qf = PolytopeQuadratureFormula::build(order, g);
    const auto report = exactnessSweep(*qf, g, order);
    EXPECT_LT(report.worstRelativeError, 1e-11)
      << name(g) << " order " << order << ": worst monomial x^"
      << report.worstExponents[0] << " y^" << report.worstExponents[1] << " z^"
      << report.worstExponents[2];
    EXPECT_TRUE(allWeightsPositive(*qf));
    EXPECT_TRUE(allPointsInside(*qf, g));
  }
}

// ---------------------------------------------------------------------------
// Invariant: properties that hold for every rule, at every order.
// ---------------------------------------------------------------------------

/**
 * @brief Checks that weights reproduce the reference measure.
 *
 * This isolates degree-zero exactness: a rescaled or renormalized rule may
 * satisfy it while failing higher moments.
 */
TEST(QuadratureExactnessTest, WeightsSumToReferenceMeasure)
{
  for (const auto g : allGeometries())
  {
    for (size_t order = 1; order <= 8; ++order)
    {
      const auto& qf = PolytopeQuadratureFormula::get(order, g);
      EXPECT_NEAR(weightSum(qf), referenceMeasure(g), kExactTol)
        << name(g) << " order " << order;
    }
  }
}

/**
 * @brief Checks that quadrature samples lie in their reference domains.
 *
 * Exterior nodes can request a geometric map outside its reference domain,
 * even when polynomial moments are exact.
 */
TEST(QuadratureExactnessTest, PointsLieInsideReferenceElement)
{
  for (const auto g : allGeometries())
  {
    for (size_t order = 1; order <= 8; ++order)
    {
      const auto& qf = PolytopeQuadratureFormula::get(order, g);
      EXPECT_TRUE(allPointsInside(qf, g)) << name(g) << " order " << order;
    }
  }
}

namespace
{
  /// @brief Asserts that two rules agree bitwise, node for node.
  template <class Rule>
  void expectIdentical(const Rule& a, const Rule& b, const std::string& where)
  {
    ASSERT_EQ(a.getSize(), b.getSize()) << where;
    for (size_t i = 0; i < a.getSize(); ++i)
    {
      EXPECT_EQ(a.getWeight(i), b.getWeight(i)) << where << " weight " << i;
      ASSERT_EQ(a.getPoint(i).size(), b.getPoint(i).size()) << where << " point " << i;
      for (Eigen::Index k = 0; k < a.getPoint(i).size(); ++k)
        EXPECT_EQ(a.getPoint(i)[k], b.getPoint(i)[k]) << where << " point " << i;
    }
  }
}

/**
 * @brief Determinism: two constructions of the same rule agree bitwise.
 *
 * The constructors are exercised directly rather than through
 * @ref Rodin::QF::PolytopeQuadratureFormula::get, which pools its rules and
 * would hand back the very same object twice -- comparing that with itself
 * would pass however the tables were read.
 */
TEST(QuadratureExactnessTest, ConstructionIsDeterministic)
{
  // Neither family covers every geometry -- Xiao-Gimbutas is simplices only,
  // and the segment is served by Gauss-Legendre -- so each is exercised only
  // where it ships rules, and the count guards against testing nothing.
  size_t compared = 0;
  for (const auto g : allGeometries())
  {
    for (size_t p = 1; p <= WitherdenVincent::getMaxDegree(g); ++p)
    {
      expectIdentical(WitherdenVincent(p, g), WitherdenVincent(p, g),
        std::string("Witherden-Vincent ") + name(g) + " degree " + std::to_string(p));
      ++compared;
    }

    for (size_t p = 1; p <= XiaoGimbutas::getMaxDegree(g); ++p)
    {
      expectIdentical(XiaoGimbutas(p, g), XiaoGimbutas(p, g),
        std::string("Xiao-Gimbutas ") + name(g) + " degree " + std::to_string(p));
      ++compared;
    }
  }
  EXPECT_GT(compared, 0u) << "no shipped rule was compared";
}

/**
 * @brief Checks that the moment oracle uses Rodin's reference domains.
 *
 * The hand-written moments are independent of the quadrature implementation.
 * Every published reference vertex and centroid must satisfy the domain's
 * half-space system. The constant moment must agree with the rule's mass.
 */
TEST(QuadratureExactnessTest, MomentOracleMatchesRodinReferenceElements)
{
  for (const auto g : allGeometries())
  {
    const Polytope::Traits traits(g);
    for (size_t i = 0; i < traits.getVertexCount(); ++i)
    {
      const auto& v = traits.getVertex(i);
      Math::SpatialVector<Real> x;
      x.resize(static_cast<Eigen::Index>(traits.getDimension()));
      for (size_t k = 0; k < traits.getDimension(); ++k)
        x[static_cast<Eigen::Index>(k)] = v[static_cast<Eigen::Index>(k)];
      EXPECT_TRUE(isInside(g, x, 1e-13)) << name(g) << ": reference vertex " << i
                                         << " does not satisfy its own half-space system";
    }

    EXPECT_TRUE(isInside(g, traits.getCentroid(), 1e-13))
      << name(g) << ": centroid outside its own element";

    const auto& qf = PolytopeQuadratureFormula::get(8, g);
    EXPECT_NEAR(weightSum(qf), referenceMeasure(g), kExactTol)
      << name(g) << ": measure from the moment oracle disagrees with the rule";
  }
}

// ---------------------------------------------------------------------------
// Characterisation: which rules are positive, and which are not.
// ---------------------------------------------------------------------------

/**
 * @brief Characterizes signed Grundmann--Möller weights above degree one.
 *
 * Negative weights do not preserve positive-semidefinite assembly in general.
 * The amplification @f$\sum|w|/\sum w@f$ distinguishes these rules from a
 * positive rule, independently of their polynomial exactness.
 */
TEST(QuadratureExactnessTest, GrundmannMollerHasSignedWeightsAboveDegreeOne)
{
  for (const auto g : {Polytope::Type::Triangle, Polytope::Type::Tetrahedron})
  {
    EXPECT_TRUE(allWeightsPositive(GrundmannMoller(0, g)))
      << name(g) << ": s = 0 (degree 1) is the one positive case";
    for (size_t s = 1; s <= 3; ++s)
    {
      const GrundmannMoller qf(s, g);
      EXPECT_FALSE(allWeightsPositive(qf)) << name(g) << " s = " << s;
      EXPECT_GT(weightAmplification(qf), Real(1)) << name(g) << " s = " << s;
      const size_t degree = 2 * s + 1;
      EXPECT_LT(exactnessSweep(qf, g, degree).worstRelativeError, kExactTol)
        << name(g) << " signed rule degree " << degree;
    }
  }
}

/**
 * @brief Checks positivity of tensor-product and collapsed rules.
 *
 * Their weights are products of positive Gaussian weights; collapse Jacobians
 * are carried by positive Gauss--Jacobi weights. This preserves the positivity
 * required by positive-semidefinite assembly. The amplification
 * @f$\sum|w|/\sum w@f$ must equal one up to the stated roundoff budget.
 */
TEST(QuadratureExactnessTest, TensorAndCollapsedRulesHaveStrictlyPositiveWeights)
{
  for (const auto g : {Polytope::Type::Segment, Polytope::Type::Quadrilateral,
         Polytope::Type::Hexahedron, Polytope::Type::Wedge, Polytope::Type::Pyramid})
  {
    for (size_t order = 1; order <= 8; ++order)
    {
      const auto& qf = PolytopeQuadratureFormula::get(order, g);
      EXPECT_TRUE(allWeightsPositive(qf)) << name(g) << " order " << order;
      EXPECT_NEAR(weightAmplification(qf), Real(1), 1e-14)
        << name(g) << " order " << order;
    }
  }
}

/// @brief The default dispatcher uses positive rules throughout the tabulated range.
TEST(QuadratureExactnessTest, DispatchedRulesHavePositiveWeights)
{
  for (const auto g : allGeometries())
  {
    const auto& qf = PolytopeQuadratureFormula::get(4, g);
    EXPECT_TRUE(allWeightsPositive(qf)) << name(g);
    EXPECT_NEAR(weightAmplification(qf), Real(1), 1e-14) << name(g);
  }
}

// ---------------------------------------------------------------------------
// Meta: the instrument has teeth.
// ---------------------------------------------------------------------------

namespace
{
  /// @brief Perturbs a selected sample to check the quadrature oracle.
  class PerturbedRule
  {
    public:
      /**
       * @brief Borrows a rule and owns the selected sample perturbation.
       * @param base Immutable source rule.
       * @param which Logical index of the altered sample.
       * @param delta Additive weight perturbation.
       * @param point Optional replacement reference coordinates.
       */
      PerturbedRule(const QuadratureFormulaBase& base, size_t which, Real delta,
        Optional<Math::SpatialVector<Real>> point = {})
        : m_base(base),
          m_which(which),
          m_delta(delta),
          m_point(std::move(point))
      {}

      /**
       * @brief Gets the unchanged source sample count.
       * @returns Number of source samples.
       */
      size_t getSize() const
      {
        return m_base.getSize();
      }

      /**
       * @brief Gets a source weight with the selected perturbation.
       * @param i Logical sample index.
       * @returns Original or perturbed weight.
       */
      Real getWeight(size_t i) const
      {
        return m_base.getWeight(i) + (i == m_which ? m_delta : Real(0));
      }

      /**
       * @brief Gets source or replacement reference coordinates.
       * @param i Logical sample index.
       * @returns Borrowed source coordinates or owned replacement coordinates.
       */
      const Math::SpatialVector<Real>& getPoint(size_t i) const
      {
        return i == m_which && m_point ? *m_point : m_base.getPoint(i);
      }

    private:
      const QuadratureFormulaBase& m_base;
      size_t m_which;
      Real m_delta;
      Optional<Math::SpatialVector<Real>> m_point;
  };
}

/**
 * @brief Checks sensitivity to a finite weight perturbation.
 *
 * A relative perturbation of @f$10^{-9}@f$, with the stated absolute floor,
 * must fail the moment budget that the unperturbed rule passes.
 */
TEST(QuadratureExactnessTest, ExactnessSweepRejectsAPerturbedWeight)
{
  for (const auto g : allGeometries())
  {
    constexpr size_t order = 4;
    const auto& qf = PolytopeQuadratureFormula::get(order, g);
    ASSERT_LT(exactnessSweep(qf, g, order).worstRelativeError, kExactTol)
      << name(g) << ": baseline must pass before the mutation means anything";

    const PerturbedRule bad(qf, 0, qf.getWeight(0) * Real(1e-9) + Real(1e-15));
    EXPECT_GT(exactnessSweep(bad, g, order).worstRelativeError, kExactTol)
      << name(g) << ": sweep failed to detect a perturbed weight";
  }
}

/// @brief The interiority check must reject a point outside the element.
TEST(QuadratureExactnessTest, InteriorityCheckRejectsAnExteriorPoint)
{
  Math::SpatialVector<Real> outside;
  outside.resize(2);
  outside[0] = 0.6;
  outside[1] = 0.6;   // x + y = 1.2 > 1
  EXPECT_FALSE(isInside(Polytope::Type::Triangle, outside, 1e-13));
  outside[0] = 0.25;
  outside[1] = 0.25;
  EXPECT_TRUE(isInside(Polytope::Type::Triangle, outside, 1e-13));
}

/**
 * @brief Domain membership rejects coordinate tuples of the wrong dimension.
 *
 * This is a logical shape condition, independent of half-space tolerances.
 * In particular, every nonempty coordinate tuple is invalid in
 * @f$\mathbb{R}^0@f$. Dimensions below the cell dimension must be rejected
 * before a component is read.
 */
TEST(QuadratureExactnessTest, InteriorityCheckRejectsWrongCoordinateDimension)
{
  for (const auto g : allGeometries())
  {
    SCOPED_TRACE(name(g));
    const Polytope::Traits traits(g);
    EXPECT_TRUE(isInside(g, traits.getCentroid(), Real(0)));
    for (size_t dimension = 0; dimension <= Math::SpatialVector<Real>::MaxSize;
         ++dimension)
    {
      if (dimension == traits.getDimension())
        continue;
      SCOPED_TRACE(dimension);
      Math::SpatialVector<Real> point;
      point.resize(static_cast<Eigen::Index>(dimension));
      point.setZero();
      EXPECT_FALSE(isInside(g, point, Real(0)));
    }
  }
}

/**
 * @brief Nonfinite coordinates do not belong to a real reference domain.
 *
 * Every active coordinate is checked. Point has no active coordinate;
 * nonempty tuples there are rejected by the dimension regression.
 */
TEST(QuadratureExactnessTest, InteriorityCheckRejectsNonfiniteCoordinates)
{
  for (const auto g : allGeometries())
  {
    SCOPED_TRACE(name(g));
    const Polytope::Traits traits(g);
    for (size_t axis = 0; axis < traits.getDimension(); ++axis)
    {
      for (const Real invalid :
        {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity(),
          -std::numeric_limits<Real>::infinity()})
      {
        auto point = traits.getCentroid();
        point[axis] = invalid;
        EXPECT_FALSE(isInside(g, point, Real(0))) << "axis " << axis;
      }
    }
  }
}

/**
 * @brief Exactness validation rejects malformed or nonfinite rule samples.
 *
 * NaN comparisons must not leave the worst-error accumulator at zero.
 * Finite data whose monomial evaluation overflows must also be rejected.
 * Valid finite signed weights remain admissible to the exactness instrument.
 */
TEST(QuadratureExactnessTest, ExactnessSweepRejectsInvalidSamples)
{
  constexpr size_t order = 4;
  for (const auto g : allGeometries())
  {
    SCOPED_TRACE(name(g));
    const Polytope::Traits traits(g);
    const auto& qf = PolytopeQuadratureFormula::get(order, g);
    ASSERT_LT(exactnessSweep(qf, g, order).worstRelativeError, kExactTol);
    for (const Real invalid : {std::numeric_limits<Real>::quiet_NaN(),
           std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity()})
    {
      EXPECT_THROW(
        exactnessSweep(PerturbedRule(qf, 0, invalid), g, order), Alert::Exception);
      for (size_t axis = 0; axis < traits.getDimension(); ++axis)
      {
        auto point = qf.getPoint(0);
        point[axis] = invalid;
        EXPECT_THROW(exactnessSweep(PerturbedRule(qf, 0, Real(0), point), g, order),
          Alert::Exception);
      }
    }
    Math::SpatialVector<Real> wrongDimension;
    wrongDimension.resize(static_cast<Eigen::Index>(
      (traits.getDimension() + 1) % (Math::SpatialVector<Real>::MaxSize + 1)));
    wrongDimension.setZero();
    EXPECT_THROW(exactnessSweep(PerturbedRule(qf, 0, Real(0), wrongDimension), g, order),
      Alert::Exception);
    if (traits.getDimension() > 0)
    {
      auto point = qf.getPoint(0);
      point[0] = std::numeric_limits<Real>::max();
      EXPECT_THROW(
        exactnessSweep(PerturbedRule(qf, 0, Real(0), point), g, order), Alert::Exception);
    }
  }
}
