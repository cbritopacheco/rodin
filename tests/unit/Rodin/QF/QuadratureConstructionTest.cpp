/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */

/**
 * @file
 * @brief Regression checks for quadrature construction and signed simplex weights.
 */
#include <gtest/gtest.h>
#include <limits>

#include "Rodin/Alert/Exception.h"
#include "Rodin/QF/GaussLegendre.h"
#include "Rodin/QF/GaussLobatto.h"
#include "Rodin/QF/GrundmannMoller.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"

using namespace Rodin;
using Type = Geometry::Polytope::Type;

TEST(QuadratureConstruction, PointSimplex)
{
  for (size_t s = 0; s <= RODIN_QF_GRUNDMANNMOLLER_MAX_S; ++s)
  {
    const QF::GrundmannMoller rule(s, Type::Point);
    ASSERT_EQ(rule.getSize(), 1u);
    EXPECT_EQ(rule.getPoint(0).size(), 0);
    EXPECT_EQ(rule.getWeight(0), 1);
  }
}

TEST(QuadratureConstruction, HighDegreeSimplexWeights)
{
  for (const auto geometry : {Type::Segment, Type::Triangle, Type::Tetrahedron})
  {
    const size_t dimension = Geometry::Polytope::Traits(geometry).getDimension();
    for (size_t s = 9; s <= RODIN_QF_GRUNDMANNMOLLER_MAX_S; ++s)
    {
      const QF::GrundmannMoller rule(s, geometry);
      for (size_t degree = 0; degree <= 2 * s + 1; ++degree)
      {
        long double sum = 0, magnitude = 0;
        for (size_t q = 0; q < rule.getSize(); ++q)
        {
          const long double term = rule.getWeight(q) *
            std::pow(static_cast<long double>(rule.getPoint(q)[0]), degree);
          sum += term;
          magnitude += std::abs(term);
        }
        long double exact = 1;
        for (size_t k = 1; k <= dimension; ++k)
          exact /= static_cast<long double>(degree + k);
        // Signed rules have growing cancellation; use an absolute forward-error bound.
        const long double tolerance =
          64 * std::numeric_limits<Real>::epsilon() * std::max(1.L, magnitude);
        EXPECT_NEAR(sum, exact, tolerance) << "s=" << s << " degree=" << degree;
      }
    }
  }
}

TEST(QuadratureConstruction, AnisotropicWedge)
{
  const QF::GaussLegendre gauss(Type::Wedge, 2, 3, 4);
  EXPECT_EQ(gauss.getSize(), 24u);
  const QF::GaussLobatto lobatto(Type::Wedge, 3, 4, 5);
  // Collapsed triangle endpoints merge; the extrusion retains all five nodes.
  EXPECT_EQ(lobatto.getSize(), 45u);
  const auto check = [](const auto& rule, size_t degree) {
    for (size_t a = 0; a <= degree; ++a)
    {
      for (size_t b = 0; b + a <= degree; ++b)
      {
        for (size_t c = 0; c <= degree; ++c)
        {
          long double sum = 0;
          for (size_t q = 0; q < rule.getSize(); ++q)
          {
            const auto& p = rule.getPoint(q);
            sum += rule.getWeight(q) * std::pow(p[0], a) * std::pow(p[1], b) *
              std::pow(p[2], c);
          }
          const long double exact = std::tgamma(a + 1.L) * std::tgamma(b + 1.L) /
            std::tgamma(a + b + 3.L) / (c + 1);
          EXPECT_NEAR(sum, exact, 2e-13);
        }
      }
    }
  };
  check(gauss, 3);
  check(lobatto, 3);
}

TEST(QuadratureConstruction, RejectsInvalidCountsAndOrders)
{
  const size_t huge = std::numeric_limits<size_t>::max();
  for (const auto geometry :
    {Type::Point, Type::Segment, Type::Triangle, Type::Quadrilateral, Type::Tetrahedron,
      Type::Wedge, Type::Hexahedron, Type::Pyramid})
  {
    EXPECT_THROW((QF::GaussLegendre(geometry, 0)), Alert::Exception);
    EXPECT_THROW((QF::GaussLobatto(geometry, 1)), Alert::Exception);
    EXPECT_THROW((QF::GaussLegendre(geometry, huge)), Alert::Exception);
    EXPECT_THROW((QF::GaussLobatto(geometry, huge)), Alert::Exception);
    if (geometry != Type::Point)
    {
      EXPECT_THROW(
        QF::PolytopeQuadratureFormula::build(huge, geometry), Alert::Exception);
    }
  }
  EXPECT_THROW((QF::GrundmannMoller(huge, Type::Triangle)), Alert::Exception);
  EXPECT_THROW((QF::GrundmannMoller(RODIN_QF_GRUNDMANNMOLLER_MAX_S + 1, Type::Triangle)),
    Alert::Exception);
  std::vector<Real> points, weights;
  EXPECT_THROW(QF::GaussLegendre::gl1dUnit(0, points, weights), Alert::Exception);
  EXPECT_THROW(QF::GaussLegendre::gj1dUnit(0, 1, points, weights), Alert::Exception);
  EXPECT_THROW(QF::GaussLegendre::gj1dUnit(huge, 1, points, weights), Alert::Exception);
}
