/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Tests for convergence-study infrastructure shared by all strategies.
 */

#include <cmath>
#include <limits>

#include <gtest/gtest.h>
#include <gtest/gtest-spi.h>

#include "Convergence.h"
#include "LiftedConvergence.h"
#include "FieldConvergence.h"
#include "Conductivity.h"
#include "LinearElasticity.h"

namespace Rodin::Tests::Convergence
{
  /** The analytic divergence-free shear field does not depend on lambda. */
  TEST(ConvergenceUtilities, DivergenceFreeElasticityHasKnownStressAndSource)
  {
    using Data = LinearElasticity::ManufacturedSolution;
    for (const auto geometry :
      {Geometry::Polytope::Type::Triangle, Geometry::Polytope::Type::Quadrilateral,
        Geometry::Polytope::Type::Tetrahedron, Geometry::Polytope::Type::Pyramid,
        Geometry::Polytope::Type::Hexahedron, Geometry::Polytope::Type::Wedge})
    {
      const auto mesh = UniformGrid(geometry).makeMesh(2);
      const size_t dim = mesh.getDimension();
      const Data soft(dim, 1.5, 1, Data::Field::DivergenceFree);
      const Data stiff(dim, 1e4, 1, Data::Field::DivergenceFree);
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        const auto rc = Geometry::Polytope::Traits(geometry).getCentroid();
        const Geometry::Point p(*cell, rc);
        const auto u = soft.exact(p), f = soft.forcing(p);
        const auto gradient = soft.jacobian(p), strain = soft.strain(p),
                   stress = soft.stress(p);
        const Real pi = Math::Constants::pi();
        EXPECT_EQ(u(0), std::sin(pi * p(1)));
        EXPECT_EQ(f(0), pi * pi * u(0));
        for (size_t i = 0; i < dim; ++i)
        {
          EXPECT_EQ(gradient(i, i), 0);
          EXPECT_EQ(u(i), stiff.exact(p)(i));
          EXPECT_EQ(f(i), stiff.forcing(p)(i));
          if (i > 0)
          {
            EXPECT_EQ(u(i), 0);
            EXPECT_EQ(f(i), 0);
          }
          for (size_t j = 0; j < dim; ++j)
          {
            const Real shear = (i == 0 && j == 1) || (i == 1 && j == 0)
              ? pi * std::cos(pi * p(1))
              : Real(0);
            EXPECT_EQ(gradient(i, j), i == 0 && j == 1 ? shear : Real(0));
            EXPECT_EQ(strain(i, j), shear / 2);
            EXPECT_EQ(stress(i, j), shear);
            EXPECT_EQ(stress(i, j), stiff.stress(p)(i, j));
          }
        }
      }
    }
  }

  TEST(ConductivityDataTest, ExponentialPhysicalDataAndSourcesOnAllGeometries)
  {
    // Floating algebra budget, unrelated to entity identification or topology.
    constexpr Real RoundoffBudget = 64 * std::numeric_limits<Real>::epsilon();
    using Type = Geometry::Polytope::Type;
    for (const auto geometry : {Type::Segment, Type::Triangle, Type::Quadrilateral,
           Type::Tetrahedron, Type::Pyramid, Type::Hexahedron, Type::Wedge})
    {
      SCOPED_TRACE(UniformGrid::getGeometryName(geometry));
      const auto mesh = UniformGrid(geometry).makeMesh(3);
      const size_t dim = mesh.getDimension();
      const auto cell = mesh.getPolytope(dim, 0);
      const Math::Vector<Real> rc = Math::Vector<Real>::Constant(dim, 0.25);
      const Geometry::Point point(*cell, rc);
      const auto& x = point.getPhysicalCoordinates();
      Real sum = 0;
      for (size_t j = 0; j < dim; ++j)
        sum += x(j);
      const Real exact = std::exp(sum);
      const ConductivityData data(dim, ConductivityData::Field::Exponential);
      EXPECT_DOUBLE_EQ(data.getSolution(x), exact);
      EXPECT_DOUBLE_EQ(data.getSolution()(point), exact);
      for (size_t j = 0; j < dim; ++j)
      {
        EXPECT_DOUBLE_EQ(data.getGradient(x)(j), exact);
        EXPECT_DOUBLE_EQ(data.getGradient()(point)(j), exact);
      }
      for (bool poisson : {false, true})
      {
        const Real gamma = poisson ? 1 : 1 + sum;
        const Real source = -Real(dim) * (gamma + (poisson ? 0 : 1)) * exact;
        EXPECT_DOUBLE_EQ(data.getCoefficient(poisson)(point), gamma);
        EXPECT_NEAR(data.getSource(poisson)(point), source,
          RoundoffBudget * std::max(Real(1), std::abs(source)));
      }
    }
  }

  TEST(FieldConvergenceTest, AcceptsEveryFieldOnNonuniformRefinementPaths)
  {
    FieldConvergence<2> algebraic, exponential;
    FieldConvergence<1> scalar;
    for (Real h : {0.5, 0.125, 0.0625})
    {
      algebraic.append(h, {ErrorNorms(h * h, h), ErrorNorms(h * h * h, h * h)});
      scalar.append(h, {ErrorNorms(h * h, h)});
    }
    for (Real p : {1, 3, 4})
      exponential.append(p,
        {ErrorNorms(std::exp(-2 * p), std::exp(-p)),
          ErrorNorms(std::exp(-3 * p), std::exp(-2 * p))});
    algebraic.expectAlgebraicFloor({1.9, 0.9});
    exponential.expectExponentialFloor({1.9, 0.9});
    scalar.expectAlgebraicFloor({1.9, 0.9});
    EXPECT_EQ(algebraic.getSize(), 3u);
    EXPECT_EQ(exponential.getSize(), 3u);
    EXPECT_EQ(scalar.getSize(), 3u);
  }

  TEST(FieldConvergenceTest, RejectsBadFinalIntervalInOtherField)
  {
    FieldConvergence<2> study;
    for (Real h : {0.5, 0.25, 0.125})
    {
      const Real other = h == 0.125 ? Real(0.25) : h;
      study.append(h, {ErrorNorms(h * h, h), ErrorNorms(2 * other * other, 2 * other)});
    }
    ::testing::TestPartResultArray failures;
    {
      ::testing::ScopedFakeTestPartResultReporter intercept(
        ::testing::ScopedFakeTestPartResultReporter::INTERCEPT_ONLY_CURRENT_THREAD,
        &failures);
      study.expectAlgebraicFloor({1.9, 0.9});
    }
    ASSERT_EQ(failures.size(), 4);
    for (int i = 0; i < failures.size(); ++i)
      EXPECT_NE(
        std::string(failures.GetTestPartResult(i).message()).find("field=1 interval=2"),
        std::string::npos);
  }

  TEST(FieldConvergenceTest, RejectsTwoLevelStudy)
  {
    FieldConvergence<1> study;
    study.append(0.5, {ErrorNorms(0.25, 0.5)}).append(0.25, {ErrorNorms(0.0625, 0.25)});
    ::testing::TestPartResultArray failures;
    {
      ::testing::ScopedFakeTestPartResultReporter intercept(
        ::testing::ScopedFakeTestPartResultReporter::INTERCEPT_ONLY_CURRENT_THREAD,
        &failures);
      study.expectAlgebraicFloor({1.9, 0.9});
    }
    ASSERT_EQ(failures.size(), 1);
    EXPECT_TRUE(failures.GetTestPartResult(0).fatally_failed());
  }

  TEST(FieldConvergenceTest, RejectsUnchangedFinalRefinementParameter)
  {
    FieldConvergence<1> algebraic, exponential;
    algebraic.append(0.5, {ErrorNorms(1, 1)})
      .append(0.25, {ErrorNorms(0.5, 0.5)})
      .append(0.25, {ErrorNorms(0.25, 0.25)});
    exponential.append(1, {ErrorNorms(1, 1)})
      .append(2, {ErrorNorms(0.5, 0.5)})
      .append(2, {ErrorNorms(0.25, 0.25)});
    ::testing::TestPartResultArray failures;
    {
      ::testing::ScopedFakeTestPartResultReporter intercept(
        ::testing::ScopedFakeTestPartResultReporter::INTERCEPT_ONLY_CURRENT_THREAD,
        &failures);
      algebraic.expectAlgebraicFloor({0.1, 0.1});
      exponential.expectExponentialFloor({0.1, 0.1});
    }
    ASSERT_EQ(failures.size(), 2);
    for (int i = 0; i < failures.size(); ++i)
    {
      EXPECT_TRUE(failures.GetTestPartResult(i).fatally_failed());
      EXPECT_NE(
        std::string(failures.GetTestPartResult(i).message()).find("field=0 interval=2"),
        std::string::npos);
    }
  }

  /** @brief Independent algebraic histories exercise every lifted rate interval. */
  TEST(LiftedConvergenceTest, AcceptsSeparateFieldAndGeometryOrders)
  {
    LiftedConvergence study;
    for (Real h : {0.25, 0.125, 0.0625})
    {
      const ErrorNorms field(h * h, h), geometry(h * h * h, h * h);
      const ErrorNorms total(field.getL2() + geometry.getL2(),
        field.getH1Seminorm() + geometry.getH1Seminorm());
      study.append(h, field, {field, geometry, total});
    }
    study.expectRates(1, 2);
  }

  TEST(LiftedConvergenceTest, RejectsBadFinalInterval)
  {
    LiftedConvergence study;
    for (Real h : {0.25, 0.125, 0.0625})
    {
      const Real scale = h == 0.0625 ? Real(0.25) : h;
      const ErrorNorms field(scale * scale, scale), geometry(h * h * h, h * h);
      const ErrorNorms total(field.getL2() + geometry.getL2(),
        field.getH1Seminorm() + geometry.getH1Seminorm());
      study.append(h, field, {field, geometry, total});
    }
    ::testing::TestPartResultArray failures;
    {
      ::testing::ScopedFakeTestPartResultReporter intercept(
        ::testing::ScopedFakeTestPartResultReporter::INTERCEPT_ONLY_CURRENT_THREAD,
        &failures);
      study.expectRates(1, 2);
    }
    ASSERT_GT(failures.size(), 0);
    for (int i = 0; i < failures.size(); ++i)
      EXPECT_NE(std::string(failures.GetTestPartResult(i).message()).find("interval=2"),
        std::string::npos);
  }

  TEST(LiftedConvergenceTest, RejectsInvalidNormDecomposition)
  {
    EXPECT_NONFATAL_FAILURE(
      LiftedConvergence::expectDecomposition({{1, 1}, {2, 2}, {4, 3}}), "values[2]");
    EXPECT_NONFATAL_FAILURE(
      LiftedConvergence::expectDecomposition({{1, 1}, {2, 2}, {0.5, 1}}), "values[2]");
  }

  /**
   * @brief Verifies that algebraic rates use the actual scale ratio.
   *
   * Here @f$h@f$ decreases by four rather than two. Errors proportional to
   * @f$h^2@f$ and @f$h@f$ must still report rates two and one respectively.
   */
  TEST(ErrorHistoryTest, ComputesAlgebraicRatesForNonDyadicSpacing)
  {
    ErrorHistory history;
    history.append(0.5, ErrorNorms(0.25, 0.5)).append(0.125, ErrorNorms(0.015625, 0.125));

    const auto rates = history.getAlgebraicRates(1);
    EXPECT_DOUBLE_EQ(rates.getL2(), 2);
    EXPECT_DOUBLE_EQ(rates.getH1Seminorm(), 1);
  }

  /**
   * @brief Verifies exponential rates for non-unit degree spacing.
   *
   * Errors @f$e^{-2p}@f$ and @f$e^{-p}@f$ must report decay constants two
   * and one even when two polynomial degrees separate the samples.
   */
  TEST(ErrorHistoryTest, ComputesExponentialRatesForNonUnitSpacing)
  {
    ErrorHistory history;
    history.append(1, ErrorNorms(std::exp(-2), std::exp(-1)))
      .append(3, ErrorNorms(std::exp(-6), std::exp(-3)));

    const auto rates = history.getExponentialRates(1);
    EXPECT_DOUBLE_EQ(rates.getL2(), 2);
    EXPECT_DOUBLE_EQ(rates.getH1Seminorm(), 1);
  }

  /** @brief Verifies geometry-independent unit-box boundary partitioning. */
  TEST(UnitBoxBoundaryTest, LabelsCoordinateSidesAndRemainder)
  {
    UniformGrid grid(Geometry::Polytope::Type::Triangle);
    auto mesh = grid.makeMesh(4);
    constexpr Geometry::Attribute lower = 11;
    constexpr Geometry::Attribute upper = 12;
    constexpr Geometry::Attribute remainder = 13;
    UnitBoxBoundary::labelCoordinatePartition(mesh, 0, lower, upper, remainder);

    size_t lowerCount = 0;
    size_t upperCount = 0;
    size_t remainderCount = 0;
    for (auto boundary = mesh.getBoundary(); boundary; ++boundary)
    {
      const auto attribute = boundary->getAttribute();
      ASSERT_TRUE(attribute.has_value());
      if (*attribute == lower)
        ++lowerCount;
      else if (*attribute == upper)
        ++upperCount;
      else if (*attribute == remainder)
        ++remainderCount;
      else
        FAIL() << "Unexpected boundary attribute " << *attribute;
    }

    EXPECT_GT(lowerCount, 0);
    EXPECT_GT(upperCount, 0);
    EXPECT_GT(remainderCount, 0);
  }
}
