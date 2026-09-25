/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief H-convergence validation for isotropic linear elasticity.
 *
 * The tests solve @f$-\nabla\cdot\sigma(u)=f@f$ with
 * @f$\sigma(u)=\lambda\,\nabla\cdot u\,I+2\mu\varepsilon(u)@f$ and measure
 * the displacement L2 error and the Frobenius norm of its Jacobian error.
 */

#include <cmath>
#include <cstdint>
#include <functional>
#include <initializer_list>
#include <string>

#include <gtest/gtest.h>

#include "Convergence.h"
#include "Rodin/Assembly.h"
#include "Rodin/Solver/CG.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Solver;
using namespace Rodin::Variational;

namespace Rodin::Tests::Convergence::H::LinearElasticity
{
  constexpr Attribute DirichletAttribute = 201;
  constexpr Attribute TractionAttribute = 202;

  using VectorCallable = std::function<Math::SpatialVector<Real>(const Point&)>;
  using MatrixCallable = std::function<Math::SpatialMatrix<Real>(const Point&)>;
  using VectorField = VectorFunction<VectorCallable>;

  struct ManufacturedSolution
  {
      VectorField exact;
      VectorField forcing;
      MatrixCallable jacobian;
  };

  ManufacturedSolution makeAffineSolution(size_t dim)
  {
    return {VectorField(dim, VectorCallable([dim](const Point& p) {
              Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
              for (size_t i = 0; i < dim; ++i)
                value(i) = Real(i + 1) * p(i) + Real(1 - Integer(i));
              return value;
            })),
      VectorField(dim, VectorCallable([dim](const Point&) {
        Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
        value.setZero();
        return value;
      })),
      MatrixCallable([dim](const Point&) {
        Math::SpatialMatrix<Real> value(
          static_cast<std::uint8_t>(dim), static_cast<std::uint8_t>(dim));
        value.setZero();
        for (size_t i = 0; i < dim; ++i)
          value(i, i) = Real(i + 1);
        return value;
      })};
  }

  ManufacturedSolution makeExponentialSolution(size_t dim, Real lambda, Real mu)
  {
    return {VectorField(dim, VectorCallable([dim](const Point& p) {
              Real sum = 0;
              for (size_t i = 0; i < dim; ++i)
                sum += p(i);
              const Real exponential = std::exp(sum);
              Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
              for (size_t i = 0; i < dim; ++i)
                value(i) = Real(i + 1) * exponential;
              return value;
            })),
      VectorField(dim, VectorCallable([dim, lambda, mu](const Point& p) {
        Real sum = 0;
        Real coefficientSum = 0;
        for (size_t i = 0; i < dim; ++i)
        {
          sum += p(i);
          coefficientSum += Real(i + 1);
        }
        const Real exponential = std::exp(sum);
        Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
        for (size_t i = 0; i < dim; ++i)
        {
          value(i) = -exponential *
            (mu * Real(dim) * Real(i + 1) + (lambda + mu) * coefficientSum);
        }
        return value;
      })),
      MatrixCallable([dim](const Point& p) {
        Real sum = 0;
        for (size_t i = 0; i < dim; ++i)
          sum += p(i);
        const Real exponential = std::exp(sum);
        Math::SpatialMatrix<Real> value(
          static_cast<std::uint8_t>(dim), static_cast<std::uint8_t>(dim));
        for (size_t i = 0; i < dim; ++i)
          for (size_t j = 0; j < dim; ++j)
            value(i, j) = Real(i + 1) * exponential;
        return value;
      })};
  }

  ManufacturedSolution makeDivergenceFreeSolution(size_t dim, Real mu)
  {
    const Real pi = Math::Constants::pi();
    const auto displacement = [dim, pi](const Point& p) {
      Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
      value.setZero();
      value(0) = std::sin(pi * p(1));
      return value;
    };

    return {VectorField(dim, VectorCallable(displacement)),
      VectorField(dim, VectorCallable([mu, pi, displacement](const Point& p) {
        return mu * pi * pi * displacement(p);
      })),
      MatrixCallable([dim, pi](const Point& p) {
        Math::SpatialMatrix<Real> value(
          static_cast<std::uint8_t>(dim), static_cast<std::uint8_t>(dim));
        value.setZero();
        value(0, 1) = pi * std::cos(pi * p(1));
        return value;
      })};
  }

  VectorField makeTraction(
    const LocalMesh& mesh, const MatrixCallable& exactJacobian, Real lambda, Real mu)
  {
    const size_t dim = mesh.getSpaceDimension();
    const BoundaryNormal normal(mesh);
    return VectorField(
      dim, VectorCallable([dim, exactJacobian, lambda, mu, normal](const Point& p) {
        const auto jacobian = exactJacobian(p);
        Real divergence = 0;
        for (size_t i = 0; i < dim; ++i)
          divergence += jacobian(i, i);

        Math::SpatialVector<Real> outward(static_cast<std::uint8_t>(dim));
        if (dim == 1)
          outward(0) = p(0) < 0.5 ? -1 : 1;
        else
          outward = normal(p);

        Math::SpatialVector<Real> traction(static_cast<std::uint8_t>(dim));
        traction.setZero();
        for (size_t i = 0; i < dim; ++i)
        {
          for (size_t j = 0; j < dim; ++j)
          {
            const Real stress =
              lambda * divergence * (i == j) + mu * (jacobian(i, j) + jacobian(j, i));
            traction(i) += stress * outward(j);
          }
        }
        return traction;
      }));
  }

  template <size_t K, class FES>
  ErrorNorms solveProblem(LocalMesh& mesh, FES& vh, const ManufacturedSolution& data,
    Real lambda, Real mu, bool mixedTraction)
  {
    TrialFunction u(vh);
    TestFunction v(vh);

    auto volumetric = Integral(lambda * Div(u), Div(v));
    auto shear = Integral(
      mu * (Jacobian(u) + Jacobian(u).T()), 0.5 * (Jacobian(v) + Jacobian(v).T()));
    auto body = Integral(data.forcing, v);
    const size_t quadratureOrder = K == 1 ? 8 : 12;
    volumetric.setOrder(quadratureOrder);
    shear.setOrder(quadratureOrder);
    body.setOrder(quadratureOrder);

    Problem linearElasticity(u, v);
    if (mixedTraction)
    {
      const auto traction = makeTraction(mesh, data.jacobian, lambda, mu);
      auto boundaryLoad = BoundaryIntegral(traction, v);
      boundaryLoad.setOrder(quadratureOrder);
      linearElasticity = volumetric + shear - body -
        boundaryLoad.over(TractionAttribute) +
        DirichletBC(u, data.exact).on(DirichletAttribute);
    }
    else
    {
      linearElasticity = volumetric + shear - body + DirichletBC(u, data.exact);
    }

    CG solver(linearElasticity);
    solver.setTolerance(1e-13).setMaxIterations(50000).solve();
    EXPECT_TRUE(solver.success());

    return ErrorNorm::computeVector(
      mesh, u.getSolution(), data.exact, data.jacobian, quadratureOrder);
  }

  template <size_t K>
  ErrorNorms solve(const UniformGridHierarchy& hierarchy, size_t pointsPerAxis,
    const ManufacturedSolution& data, Real lambda, Real mu, bool mixedTraction = false)
  {
    auto mesh = hierarchy.makeMesh(pointsPerAxis);
    if (mixedTraction)
    {
      UnitBoxBoundary::labelCoordinatePartition(
        mesh, 0, DirichletAttribute, TractionAttribute, TractionAttribute);
    }

    const size_t dim = mesh.getSpaceDimension();
    if constexpr (K == 1)
    {
      P1 vh(mesh, dim);
      return solveProblem<K>(mesh, vh, data, lambda, mu, mixedTraction);
    }
    else
    {
      H1 vh(std::integral_constant<size_t, K>{}, mesh, dim);
      return solveProblem<K>(mesh, vh, data, lambda, mu, mixedTraction);
    }
  }

  void expectRates(const ErrorHistory& history, Real minimumL2Rate, Real maximumL2Rate,
    Real minimumH1Rate, Real maximumH1Rate)
  {
    ASSERT_GE(history.getSize(), 3);
    for (size_t i = 1; i < history.getSize(); ++i)
    {
      const auto& coarse = history.getSample(i - 1).error;
      const auto& fine = history.getSample(i).error;
      ASSERT_TRUE(coarse.isFinite());
      ASSERT_TRUE(fine.isFinite());
      ASSERT_GT(coarse.getL2(), fine.getL2());
      ASSERT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());

      const auto rates = history.getAlgebraicRates(i);
      SCOPED_TRACE(::testing::Message()
        << "L2 errors " << coarse.getL2() << " -> " << fine.getL2() << ", rate "
        << rates.getL2() << "; H1-seminorm errors " << coarse.getH1Seminorm() << " -> "
        << fine.getH1Seminorm() << ", rate " << rates.getH1Seminorm());
      EXPECT_GT(rates.getL2(), minimumL2Rate);
      EXPECT_LT(rates.getL2(), maximumL2Rate);
      EXPECT_GT(rates.getH1Seminorm(), minimumH1Rate);
      EXPECT_LT(rates.getH1Seminorm(), maximumH1Rate);
    }
  }

  template <size_t K, class Factory>
  void testConvergence(Polytope::Type geometry, std::initializer_list<size_t> levels,
    Factory&& makeSolution, Real lambda, Real mu, bool mixedTraction, Real minimumL2Rate,
    Real maximumL2Rate, Real minimumH1Rate, Real maximumH1Rate)
  {
    const UniformGridHierarchy hierarchy(geometry, levels);
    const auto data = makeSolution(hierarchy.getDimension(), lambda, mu);
    ErrorHistory history;
    for (const size_t level : hierarchy.getLevels())
    {
      history.append(hierarchy.getMeshSize(level),
        solve<K>(hierarchy, level, data, lambda, mu, mixedTraction));
    }
    expectRates(history, minimumL2Rate, maximumL2Rate, minimumH1Rate, maximumH1Rate);
  }

  class LinearElasticityHConvergenceTest : public ::testing::TestWithParam<Polytope::Type>
  {};

  /** @brief Verifies exact reproduction of affine P1 displacements. */
  TEST_P(LinearElasticityHConvergenceTest, AffineP1PatchTestIsExact)
  {
    const UniformGridHierarchy hierarchy(GetParam(), {3, 5, 9});
    const auto data = makeAffineSolution(hierarchy.getDimension());
    const auto error = solve<1>(hierarchy, 5, data, 1.5, 0.5);
    EXPECT_LT(error.getL2(), 1e-10);
    EXPECT_LT(error.getH1Seminorm(), 1e-10);
  }

  /** @brief Verifies optimal P1 displacement rates with full Dirichlet data. */
  TEST_P(LinearElasticityHConvergenceTest, SmoothDirichletHasOptimalP1Rates)
  {
    testConvergence<1>(GetParam(), {5, 9, 17}, makeExponentialSolution, 1.5, 0.5, false,
      1.65, 2.35, 0.75, 1.25);
  }

  /** @brief Verifies optimal P2 displacement rates with full Dirichlet data. */
  TEST_P(LinearElasticityHConvergenceTest, SmoothDirichletHasOptimalP2Rates)
  {
    testConvergence<2>(GetParam(), {3, 5, 9}, makeExponentialSolution, 1.5, 0.5, false,
      2.45, 3.55, 1.55, 2.45);
  }

  /** @brief Verifies optimal P1 rates with mixed displacement and traction. */
  TEST_P(LinearElasticityHConvergenceTest, MixedTractionHasOptimalP1Rates)
  {
    if (GetParam() == Polytope::Type::Tetrahedron)
    {
      // The coarsest tetrahedral grid is pre-asymptotic for the dual L2
      // estimate at the displacement/traction interface.
      testConvergence<1>(GetParam(), {9, 17, 33}, makeExponentialSolution, 1.5, 0.5, true,
        1.65, 2.35, 0.75, 1.25);
    }
    else
    {
      testConvergence<1>(GetParam(), {5, 9, 17}, makeExponentialSolution, 1.5, 0.5, true,
        1.65, 2.35, 0.75, 1.25);
    }
  }

  /** @brief Verifies optimal P2 rates with mixed displacement and traction. */
  TEST_P(LinearElasticityHConvergenceTest, MixedTractionHasOptimalP2Rates)
  {
    testConvergence<2>(GetParam(), {3, 5, 9}, makeExponentialSolution, 1.5, 0.5, true,
      2.45, 3.55, 1.55, 2.45);
  }

  /**
   * @brief Verifies P2 convergence for a divergence-free displacement at
   * nearly incompressible material parameters.
   */
  TEST_P(
    LinearElasticityHConvergenceTest, NearlyIncompressibleDivergenceFreeHasOptimalP2Rates)
  {
    const size_t dim = Polytope::Traits(GetParam()).getDimension();
    if (dim == 1)
      GTEST_SKIP() << "A nonzero divergence-free displacement is impossible in 1D.";

    const auto factory = [](size_t dimension, Real, Real mu) {
      return makeDivergenceFreeSolution(dimension, mu);
    };
    testConvergence<2>(
      GetParam(), {3, 5, 9}, factory, 1e4, 1, false, 2.25, 3.75, 1.35, 2.65);
  }

  std::string geometryName(const ::testing::TestParamInfo<Polytope::Type>& info)
  {
    return std::string(UniformGrid::getGeometryName(info.param));
  }

  INSTANTIATE_TEST_SUITE_P(AllUniformGridGeometries, LinearElasticityHConvergenceTest,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    geometryName);
}
