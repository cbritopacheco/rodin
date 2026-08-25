/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

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

  /// @brief Every type Rodin publishes, minus Point, which carries no
  /// quadrature of its own. Taken from Polytope::Types so that a new element
  /// type is covered by these tests the moment it is declared, rather than
  /// when someone remembers to extend a list here.
  std::vector<Polytope::Type> allGeometries()
  {
    std::vector<Polytope::Type> gs;
    for (const auto g : Polytope::Types)
      if (g != Polytope::Type::Point)
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

/// @brief Every dispatched rule is exact to its requested degree, on every
/// element type. This is the table-integrity gate: a wrong coefficient shows
/// up here and nowhere else.
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

/// @brief The positive conical-product fallbacks are exact beyond the
/// published simplex tables.
TEST(QuadratureExactnessTest, ConicalProductFallbacksAreExact)
{
  for (const auto& [g, order] : {std::pair{Polytope::Type::Triangle, size_t{51}},
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

/// @brief The weights reproduce the measure of the reference element. This is
/// the degree-0 case of exactness, isolated because it is the one condition a
/// rescaled or renormalised rule can satisfy while failing every other.
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

/// @brief Every quadrature point lies in its reference element. A rule with
/// exterior points is unusable on curved cells even when it integrates
/// polynomials correctly, because the geometric map need not be defined there.
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

/// @brief Determinism: two constructions of the same rule agree bitwise.
TEST(QuadratureExactnessTest, ConstructionIsDeterministic)
{
  for (const auto g : allGeometries())
  {
    GrundmannMoller a(2, g == Polytope::Type::Triangle ? g : Polytope::Type::Triangle);
    GrundmannMoller b(2, Polytope::Type::Triangle);
    ASSERT_EQ(a.getSize(), b.getSize());
    for (size_t i = 0; i < a.getSize(); ++i)
    {
      EXPECT_EQ(a.getWeight(i), b.getWeight(i)) << "weight " << i;
      for (Eigen::Index k = 0; k < a.getPoint(i).size(); ++k)
        EXPECT_EQ(a.getPoint(i)[k], b.getPoint(i)[k]) << "point " << i;
    }
    break;
  }
}

/// @brief The closed-form moments are moments of *Rodin's* reference elements.
///
/// exactMoment() is written by hand so that the exactness sweep has an oracle
/// independent of any quadrature code. That independence is only worth having
/// if the domain it integrates over is the one Rodin actually uses, so this
/// test pins the two together: every reference vertex published by
/// Polytope::Traits must satisfy the half-space system, and the measure
/// implied by the moments must agree with the measure obtained by integrating
/// the constant 1 with a rule of the highest order under test.
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

/// @brief Grundmann-Möller has negative weights at every degree above one.
///
/// This is a characterisation test, not an aspiration: it documents why GM
/// cannot be used to assemble a form that must stay positive semidefinite, and
/// it will fail loudly if GM is ever silently swapped for a positive rule.
/// The amplification @f$ \sum|w|/\sum w @f$ is the quantity that matters.
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
    }
  }
}

/// @brief Every tensor-product and collapsed rule has positive weights.
///
/// This is load-bearing rather than incidental. Negative weights are the
/// reason Grundmann-Möller cannot assemble a form that must stay positive
/// semidefinite, and the reason a positive simplex family is being introduced
/// at all. That argument does not extend to the tensor-product and collapsed
/// elements, whose weights are products of positive Gaussian weights. The
/// collapse Jacobians are carried by positive Gauss--Jacobi weights.
///
/// Because that structural fact is what justifies leaving these elements on
/// their existing rules, it is asserted here rather than assumed. The
/// amplification @f$ \sum|w| / \sum w @f$ must be exactly one: any departure
/// means a sign appeared somewhere in the product.
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
  /// @brief A rule wrapper that perturbs one weight, used to prove the
  /// exactness sweep rejects a corrupted table.
  class PerturbedRule
  {
    public:
      PerturbedRule(const QuadratureFormulaBase& base, size_t which, Real delta)
        : m_base(base),
          m_which(which),
          m_delta(delta)
      {}

      size_t getSize() const
      {
        return m_base.getSize();
      }

      Real getWeight(size_t i) const
      {
        return m_base.getWeight(i) + (i == m_which ? m_delta : Real(0));
      }

      const Math::SpatialVector<Real>& getPoint(size_t i) const
      {
        return m_base.getPoint(i);
      }

    private:
      const QuadratureFormulaBase& m_base;
      size_t m_which;
      Real m_delta;
  };
}

/// @brief A tolerance test that cannot fail proves nothing. Perturbing a single
/// weight by one part in @f$ 10^{9} @f$ must be rejected by the sweep that the
/// unperturbed rule passes.
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
