/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief H-convergence validation for the Taylor--Hood Stokes discretization.
 */

#include <cstdint>
#include <functional>
#include <initializer_list>
#include <string>

#include <gtest/gtest.h>

#include "Convergence.h"
#include "Rodin/Assembly.h"
#include "Rodin/Solver/SparseLU.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Solver;
using namespace Rodin::Variational;

namespace Rodin::Tests::Convergence::H::Stokes
{
  using VectorCallable =
    std::function<Math::SpatialVector<Real>(const Point&)>;
  using ScalarCallable = std::function<Real(const Point&)>;
  using MatrixCallable =
    std::function<Math::SpatialMatrix<Real>(const Point&)>;
  using VectorField = VectorFunction<VectorCallable>;
  using ScalarField = RealFunction<ScalarCallable>;

  struct ManufacturedSolution
  {
    VectorField velocity;
    ScalarField pressure;
    VectorField forcing;
    MatrixCallable velocityJacobian;
    VectorField pressureGradient;
  };

  ManufacturedSolution makeAffineSolution(size_t dim)
  {
    return {
      VectorField(dim, VectorCallable([dim](const Point& p)
        {
          Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
          value.setZero();
          value(0) = p(1);
          return value;
        })),
      ScalarField(ScalarCallable(
        [](const Point& p) { return p(0) - 0.5; })),
      VectorField(dim, VectorCallable([dim](const Point&)
        {
          Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
          value.setZero();
          value(0) = 1;
          return value;
        })),
      MatrixCallable([dim](const Point&)
        {
          Math::SpatialMatrix<Real> value(
            static_cast<std::uint8_t>(dim), static_cast<std::uint8_t>(dim));
          value.setZero();
          value(0, 1) = 1;
          return value;
        }),
      VectorField(dim, VectorCallable([dim](const Point&)
        {
          Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
          value.setZero();
          value(0) = 1;
          return value;
        }))
    };
  }

  ManufacturedSolution makePolynomialSolution(size_t dim)
  {
    return {
      VectorField(dim, VectorCallable([dim](const Point& p)
        {
          Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
          value.setZero();
          value(0) = p(1) * p(1) * p(1);
          return value;
        })),
      ScalarField(ScalarCallable([](const Point& p)
        { return p(0) * p(0) - Real(1) / 3; })),
      VectorField(dim, VectorCallable([dim](const Point& p)
        {
          Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
          value.setZero();
          value(0) = -6 * p(1) + 2 * p(0);
          return value;
        })),
      MatrixCallable([dim](const Point& p)
        {
          Math::SpatialMatrix<Real> value(
            static_cast<std::uint8_t>(dim), static_cast<std::uint8_t>(dim));
          value.setZero();
          value(0, 1) = 3 * p(1) * p(1);
          return value;
        }),
      VectorField(dim, VectorCallable([dim](const Point& p)
        {
          Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
          value.setZero();
          value(0) = 2 * p(0);
          return value;
        }))
    };
  }

  struct StokesErrors
  {
    ErrorNorms velocity;
    ErrorNorms pressure;
    Real divergence;
  };

  StokesErrors solve(
    const UniformGridHierarchy& hierarchy,
    size_t pointsPerAxis,
    const ManufacturedSolution& data)
  {
    auto mesh = hierarchy.makeMesh(pointsPerAxis);
    const size_t dim = mesh.getSpaceDimension();
    constexpr size_t quadratureOrder = 12;

    H1 velocitySpace(std::integral_constant<size_t, 2>{}, mesh, dim);
    H1 pressureSpace(std::integral_constant<size_t, 1>{}, mesh);
    P0g meanSpace(mesh);

    TrialFunction u(velocitySpace);
    TrialFunction p(pressureSpace);
    TrialFunction lambda(meanSpace);
    TestFunction v(velocitySpace);
    TestFunction q(pressureSpace);
    TestFunction mu(meanSpace);

    auto viscosity = Integral(Jacobian(u), Jacobian(v));
    auto pressureVelocity = Integral(p, Div(v));
    auto incompressibility = Integral(Div(u), q);
    auto gaugePressure = Integral(lambda, q);
    auto gaugeMean = Integral(p, mu);
    auto body = Integral(data.forcing, v);
    viscosity.setOrder(quadratureOrder);
    pressureVelocity.setOrder(quadratureOrder);
    incompressibility.setOrder(quadratureOrder);
    gaugePressure.setOrder(quadratureOrder);
    gaugeMean.setOrder(quadratureOrder);
    body.setOrder(quadratureOrder);

    Problem stokes(u, p, lambda, v, q, mu);
    stokes = viscosity
           - pressureVelocity
           + incompressibility
           + gaugePressure
           + gaugeMean
           - body
           + DirichletBC(u, data.velocity);

    SparseLU solver(stokes);
    solver.solve();

    const auto velocity = ErrorNorm::computeVector(
      mesh, u.getSolution(), data.velocity, data.velocityJacobian,
      quadratureOrder);
    const auto pressure = ErrorNorm::compute(
      mesh, p.getSolution(), data.pressure, data.pressureGradient,
      quadratureOrder);
    const Real divergenceError = ErrorNorm::computeDivergenceL2(
      mesh, u.getSolution(), quadratureOrder);
    return {velocity, pressure, divergenceError};
  }

  void expectRates(
    const ErrorHistory& velocity,
    const ErrorHistory& pressure)
  {
    ASSERT_EQ(velocity.getSize(), pressure.getSize());
    ASSERT_GE(velocity.getSize(), 3);
    for (size_t i = 1; i < velocity.getSize(); ++i)
    {
      const auto& coarseVelocity = velocity.getSample(i - 1).error;
      const auto& fineVelocity = velocity.getSample(i).error;
      const auto& coarsePressure = pressure.getSample(i - 1).error;
      const auto& finePressure = pressure.getSample(i).error;
      ASSERT_TRUE(coarseVelocity.isFinite());
      ASSERT_TRUE(fineVelocity.isFinite());
      ASSERT_TRUE(coarsePressure.isFinite());
      ASSERT_TRUE(finePressure.isFinite());
      ASSERT_GT(coarseVelocity.getL2(), fineVelocity.getL2());
      ASSERT_GT(coarseVelocity.getH1Seminorm(), fineVelocity.getH1Seminorm());
      ASSERT_GT(coarsePressure.getL2(), finePressure.getL2());
      ASSERT_GT(coarsePressure.getH1Seminorm(), finePressure.getH1Seminorm());

      const auto velocityRates = velocity.getAlgebraicRates(i);
      const auto pressureRates = pressure.getAlgebraicRates(i);
      SCOPED_TRACE(::testing::Message()
        << "velocity L2 " << coarseVelocity.getL2() << " -> "
        << fineVelocity.getL2() << " (rate " << velocityRates.getL2()
        << "), velocity H1 " << coarseVelocity.getH1Seminorm() << " -> "
        << fineVelocity.getH1Seminorm() << " (rate "
        << velocityRates.getH1Seminorm() << "); pressure L2 "
        << coarsePressure.getL2() << " -> " << finePressure.getL2()
        << " (rate " << pressureRates.getL2() << "), pressure H1 "
        << coarsePressure.getH1Seminorm() << " -> "
        << finePressure.getH1Seminorm() << " (rate "
        << pressureRates.getH1Seminorm() << ')');
      EXPECT_GT(velocityRates.getL2(), 2.45);
      EXPECT_LT(velocityRates.getL2(), 3.55);
      EXPECT_GT(velocityRates.getH1Seminorm(), 1.55);
      EXPECT_LT(velocityRates.getH1Seminorm(), 2.45);
      EXPECT_GT(pressureRates.getL2(), 1.45);
      EXPECT_LT(pressureRates.getL2(), 2.55);
      EXPECT_GT(pressureRates.getH1Seminorm(), 0.55);
      EXPECT_LT(pressureRates.getH1Seminorm(), 1.45);
    }
  }

  class StokesHConvergenceTest
    : public ::testing::TestWithParam<Polytope::Type>
  {};

  /** @brief Verifies exact reproduction by the Taylor--Hood pair. */
  TEST_P(StokesHConvergenceTest, AffineDivergenceFreeSolutionIsExact)
  {
    const UniformGridHierarchy hierarchy(GetParam(), {3});
    const auto error = solve(hierarchy, 3,
      makeAffineSolution(hierarchy.getDimension()));
    EXPECT_LT(error.velocity.getL2(), 1e-10);
    EXPECT_LT(error.velocity.getH1Seminorm(), 1e-10);
    EXPECT_LT(error.pressure.getL2(), 1e-10);
    EXPECT_LT(error.pressure.getH1Seminorm(), 1e-10);
    EXPECT_LT(error.divergence, 1e-10);
  }

  /** @brief Verifies Taylor--Hood velocity and pressure h-rates. */
  TEST_P(StokesHConvergenceTest, PolynomialDivergenceFreeSolutionHasOptimalRates)
  {
    const UniformGridHierarchy hierarchy(GetParam(), {3, 5, 9});
    const auto data = makePolynomialSolution(hierarchy.getDimension());
    ErrorHistory velocity;
    ErrorHistory pressure;
    for (const size_t level : hierarchy.getLevels())
    {
      const auto error = solve(hierarchy, level, data);
      const Real h = hierarchy.getMeshSize(level);
      velocity.append(h, error.velocity);
      pressure.append(h, error.pressure);
    }
    expectRates(velocity, pressure);
  }

  std::string geometryName(
    const ::testing::TestParamInfo<Polytope::Type>& info)
  {
    return std::string(UniformGrid::getGeometryName(info.param));
  }

  INSTANTIATE_TEST_SUITE_P(
    AllUniformGridGeometries,
    StokesHConvergenceTest,
    ::testing::Values(
      Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron,
      Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron,
      Polytope::Type::Wedge),
    geometryName);
}
