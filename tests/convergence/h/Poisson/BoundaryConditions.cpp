/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Poisson h-convergence with Neumann and Robin boundary conditions.
 *
 * For @f$u=\exp(\sum_i x_i)@f$, the manufactured volume load is
 * @f$f=-d u@f$. On a natural boundary the exact flux is
 * @f$g=\nabla u\cdot n@f$, while Robin data with coefficient @f$\alpha@f$
 * is @f$r=g+\alpha u@f$. The pure-Neumann solution is shifted to have zero
 * mean and a global Lagrange multiplier enforces that normalization.
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
#include "Rodin/Solver/SparseLU.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Solver;
using namespace Rodin::Variational;

namespace Rodin::Tests::Convergence::H::Poisson
{
  constexpr Attribute DirichletAttribute = 101;
  constexpr Attribute NaturalAttribute = 102;
  constexpr Real RobinCoefficient = 2;

  using ScalarCallable = std::function<Real(const Point&)>;
  using VectorCallable = std::function<Math::SpatialVector<Real>(const Point&)>;
  using ScalarField = RealFunction<ScalarCallable>;
  using VectorField = VectorFunction<VectorCallable>;

  enum class BoundaryCondition
  {
    DirichletNeumann,
    DirichletRobin,
    PureNeumann
  };

  struct ExponentialManufacturedSolution
  {
      ExponentialManufacturedSolution(size_t dim, bool zeroMean)
        : exact([dim, zeroMean](const Point& p) {
            Real sum = 0;
            for (size_t i = 0; i < dim; ++i)
              sum += p(i);
            const Real mean =
              zeroMean ? std::pow(std::expm1(Real(1)), Real(dim)) : Real(0);
            return std::exp(sum) - mean;
          }),
          forcing([dim](const Point& p) {
            Real sum = 0;
            for (size_t i = 0; i < dim; ++i)
              sum += p(i);
            return -Real(dim) * std::exp(sum);
          }),
          gradient(dim, [dim](const Point& p) {
            Real sum = 0;
            for (size_t i = 0; i < dim; ++i)
              sum += p(i);
            Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
            value.setConstant(std::exp(sum));
            return value;
          })
      {}

      ScalarField exact;
      ScalarField forcing;
      VectorField gradient;
  };

  ScalarField makeNormalFlux(const LocalMesh& mesh, const VectorField& exactGradient)
  {
    const size_t dim = mesh.getSpaceDimension();
    const BoundaryNormal normal(mesh);
    return ScalarField(ScalarCallable([dim, exactGradient, normal](const Point& p) {
      const auto gradient = exactGradient(p);
      if (dim == 1)
        return (p(0) < 0.5 ? -1 : 1) * gradient(0);

      const auto outward = normal(p);
      Real flux = 0;
      for (size_t i = 0; i < dim; ++i)
        flux += gradient(i) * outward(i);
      return flux;
    }));
  }

  template <size_t K, class FES>
  ErrorNorms solveBoundaryProblem(LocalMesh& mesh, FES& vh, BoundaryCondition condition)
  {
    const bool pureNeumann = condition == BoundaryCondition::PureNeumann;
    const ExponentialManufacturedSolution data(mesh.getSpaceDimension(), pureNeumann);
    const ScalarField normalFlux = makeNormalFlux(mesh, data.gradient);

    TrialFunction u(vh);
    TestFunction v(vh);

    auto stiffness = Integral(Grad(u), Grad(v));
    auto load = Integral(data.forcing, v);
    auto naturalLoad = BoundaryIntegral(normalFlux, v);
    const size_t quadratureOrder = K == 1 ? 8 : 12;
    if constexpr (K > 1)
      stiffness.setOrder(quadratureOrder);
    load.setOrder(quadratureOrder);
    naturalLoad.setOrder(quadratureOrder);

    if (condition == BoundaryCondition::DirichletNeumann)
    {
      Problem poisson(u, v);
      poisson = stiffness - load - naturalLoad.over(NaturalAttribute) +
        DirichletBC(u, data.exact).on(DirichletAttribute);

      CG solver(poisson);
      solver.setTolerance(1e-13).setMaxIterations(20000).solve();
      EXPECT_TRUE(solver.success());
    }
    else if (condition == BoundaryCondition::DirichletRobin)
    {
      const ScalarField robinData(
        ScalarCallable([normalFlux, exact = data.exact](const Point& p) {
          return normalFlux(p) + RobinCoefficient * exact(p);
        }));
      auto robinMass = BoundaryIntegral(u, v);
      auto robinLoad = BoundaryIntegral(robinData, v);
      robinMass.setOrder(quadratureOrder);
      robinLoad.setOrder(quadratureOrder);

      Problem poisson(u, v);
      poisson = stiffness + RobinCoefficient * robinMass.over(NaturalAttribute) - load -
        robinLoad.over(NaturalAttribute) +
        DirichletBC(u, data.exact).on(DirichletAttribute);

      CG solver(poisson);
      solver.setTolerance(1e-13).setMaxIterations(20000).solve();
      EXPECT_TRUE(solver.success());
    }
    else
    {
      P0g constants(mesh);
      TrialFunction lambda(constants);
      TestFunction mu(constants);

      Problem poisson(u, lambda, v, mu);
      poisson = stiffness + Integral(lambda, v) + Integral(u, mu) - load -
        naturalLoad.over(NaturalAttribute);

      SparseLU solver(poisson);
      solver.solve();
    }

    return ErrorNorm::compute(
      mesh, u.getSolution(), data.exact, data.gradient, quadratureOrder);
  }

  template <size_t K>
  ErrorNorms solve(const UniformGridHierarchy& hierarchy, size_t pointsPerAxis,
    BoundaryCondition condition)
  {
    auto mesh = hierarchy.makeMesh(pointsPerAxis);
    if (condition == BoundaryCondition::PureNeumann)
    {
      UnitBoxBoundary::labelCoordinatePartition(
        mesh, 0, NaturalAttribute, NaturalAttribute, NaturalAttribute);
    }
    else
    {
      UnitBoxBoundary::labelCoordinatePartition(
        mesh, 0, DirichletAttribute, NaturalAttribute, NaturalAttribute);
    }

    if constexpr (K == 1)
    {
      P1 vh(mesh);
      return solveBoundaryProblem<K>(mesh, vh, condition);
    }
    else
    {
      H1 vh(std::integral_constant<size_t, K>{}, mesh);
      return solveBoundaryProblem<K>(mesh, vh, condition);
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

  template <size_t K>
  void testBoundaryCondition(Polytope::Type geometry, BoundaryCondition condition,
    std::initializer_list<size_t> levels, Real minimumL2Rate, Real maximumL2Rate,
    Real minimumH1Rate, Real maximumH1Rate)
  {
    const UniformGridHierarchy hierarchy(geometry, levels);
    ErrorHistory history;
    for (const size_t level : hierarchy.getLevels())
      history.append(hierarchy.getMeshSize(level), solve<K>(hierarchy, level, condition));
    expectRates(history, minimumL2Rate, maximumL2Rate, minimumH1Rate, maximumH1Rate);
  }

  class PoissonBoundaryHConvergenceTest : public ::testing::TestWithParam<Polytope::Type>
  {};

  /** @brief Verifies optimal P1 rates with Dirichlet--Neumann data. */
  TEST_P(PoissonBoundaryHConvergenceTest, DirichletNeumannHasOptimalP1Rates)
  {
    testBoundaryCondition<1>(GetParam(), BoundaryCondition::DirichletNeumann, {5, 9, 17},
      1.65, 2.35, 0.75, 1.25);
  }

  /** @brief Verifies optimal P2 rates with Dirichlet--Neumann data. */
  TEST_P(PoissonBoundaryHConvergenceTest, DirichletNeumannHasOptimalP2Rates)
  {
    testBoundaryCondition<2>(
      GetParam(), BoundaryCondition::DirichletNeumann, {3, 5, 9}, 2.45, 3.55, 1.55, 2.45);
  }

  /** @brief Verifies optimal P3 rates with Dirichlet--Neumann data. */
  TEST_P(PoissonBoundaryHConvergenceTest, DirichletNeumannHasOptimalP3Rates)
  {
    testBoundaryCondition<3>(
      GetParam(), BoundaryCondition::DirichletNeumann, {3, 5, 9}, 3.25, 4.75, 2.35, 3.65);
  }

  /** @brief Verifies optimal P1 rates with Dirichlet--Robin data. */
  TEST_P(PoissonBoundaryHConvergenceTest, DirichletRobinHasOptimalP1Rates)
  {
    testBoundaryCondition<1>(
      GetParam(), BoundaryCondition::DirichletRobin, {5, 9, 17}, 1.65, 2.35, 0.75, 1.25);
  }

  /** @brief Verifies optimal P2 rates with Dirichlet--Robin data. */
  TEST_P(PoissonBoundaryHConvergenceTest, DirichletRobinHasOptimalP2Rates)
  {
    testBoundaryCondition<2>(
      GetParam(), BoundaryCondition::DirichletRobin, {3, 5, 9}, 2.45, 3.55, 1.55, 2.45);
  }

  /** @brief Verifies optimal P3 rates with Dirichlet--Robin data. */
  TEST_P(PoissonBoundaryHConvergenceTest, DirichletRobinHasOptimalP3Rates)
  {
    testBoundaryCondition<3>(
      GetParam(), BoundaryCondition::DirichletRobin, {3, 5, 9}, 3.25, 4.75, 2.35, 3.65);
  }

  /** @brief Verifies optimal P1 pure-Neumann rates and mean normalization. */
  TEST_P(PoissonBoundaryHConvergenceTest, PureNeumannHasOptimalP1Rates)
  {
    testBoundaryCondition<1>(
      GetParam(), BoundaryCondition::PureNeumann, {5, 9, 17}, 1.65, 2.35, 0.75, 1.25);
  }

  /** @brief Verifies optimal P2 pure-Neumann rates and mean normalization. */
  TEST_P(PoissonBoundaryHConvergenceTest, PureNeumannHasOptimalP2Rates)
  {
    testBoundaryCondition<2>(
      GetParam(), BoundaryCondition::PureNeumann, {3, 5, 9}, 2.45, 3.55, 1.55, 2.45);
  }

  /** @brief Verifies optimal P3 pure-Neumann rates and mean normalization. */
  TEST_P(PoissonBoundaryHConvergenceTest, PureNeumannHasOptimalP3Rates)
  {
    testBoundaryCondition<3>(
      GetParam(), BoundaryCondition::PureNeumann, {3, 5, 9}, 3.25, 4.75, 2.35, 3.65);
  }

  std::string geometryName(const ::testing::TestParamInfo<Polytope::Type>& info)
  {
    return std::string(UniformGrid::getGeometryName(info.param));
  }

  INSTANTIATE_TEST_SUITE_P(AllUniformGridGeometries, PoissonBoundaryHConvergenceTest,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    geometryName);
}
