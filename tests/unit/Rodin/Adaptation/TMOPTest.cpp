/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cmath>
#include <initializer_list>
#include <array>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include <Rodin/Adaptation.h>
#include <Rodin/Geometry/LevelSetDiscretizerTriangles.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Adaptation::TMOP;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  static LocalMesh makeUnitTriangle()
  {
    auto mesh = LocalMesh::Builder()
      .initialize(2)
      .nodes(3)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .finalize();
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    return mesh;
  }

  static LocalMesh makeTriangle(
      const Math::SpatialPoint& a,
      const Math::SpatialPoint& b,
      const Math::SpatialPoint& c)
  {
    auto mesh = LocalMesh::Builder()
      .initialize(2)
      .nodes(3)
      .vertex(a)
      .vertex(b)
      .vertex(c)
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .finalize();
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    return mesh;
  }

  static LocalMesh makeTwoTriangleSquare()
  {
    auto mesh = LocalMesh::Builder()
      .initialize(2)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .polytope(Polytope::Type::Triangle, {1, 3, 2})
      .finalize();
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    return mesh;
  }

  static HighOrderTriangleGeometry upgradeUnitTriangle()
  {
    auto mesh = makeUnitTriangle();
    return HighOrderGeometryUpgrade().upgrade(mesh, 2);
  }

  static HighOrderTriangleGeometry makeInvertedGeometry()
  {
    HighOrderTriangleGeometry geometry;
    geometry.order = 2;
    geometry.nodes = {
      { Math::SpatialPoint({0, 0}), Math::SpatialPoint({0, 0}), true },
      { Math::SpatialPoint({0, 1}), Math::SpatialPoint({0, 1}), true },
      { Math::SpatialPoint({1, 0}), Math::SpatialPoint({1, 0}), true },
      { Math::SpatialPoint({0, 0.5}), Math::SpatialPoint({0, 0.5}), false },
      { Math::SpatialPoint({0.5, 0.5}), Math::SpatialPoint({0.5, 0.5}), false },
      { Math::SpatialPoint({0.5, 0}), Math::SpatialPoint({0.5, 0}), false } };
    geometry.cells.push_back({{0, 1, 2, 3, 4, 5}});
    return geometry;
  }

  static void distortMidpoint(HighOrderTriangleGeometry& geometry)
  {
    const Index midpoint = geometry.cells.front()[3];
    geometry.nodes[midpoint].x[1] += 0.25;
  }

  TEST(Rodin_Adaptation_TMOP, SquaredDistanceMetricIsZeroAtIdentity)
  {
    SquaredDistanceMetric metric;
    EXPECT_NEAR(metric.value(Matrix2::Identity()), 0, 1e-14);
  }

  TEST(Rodin_Adaptation_TMOP, DistortedMatrixHasLargerMetric)
  {
    SquaredDistanceMetric metric;
    Matrix2 A = Matrix2::Identity();
    A(0, 1) = 0.5;

    EXPECT_GT(metric.value(A), metric.value(Matrix2::Identity()));
  }

  TEST(Rodin_Adaptation_TMOP, InvertedAndNearSingularMatricesArePenalized)
  {
    SquaredDistanceMetric metric;
    Matrix2 inverted = Matrix2::Identity();
    inverted(0, 0) = -1;

    Matrix2 nearSingular = Matrix2::Identity();
    nearSingular(1, 1) = 1e-14;

    EXPECT_GT(metric.value(inverted), 1e10);
    EXPECT_GT(metric.value(nearSingular), 1e10);
  }

  TEST(Rodin_Adaptation_TMOP, ShapeAndAreaMetricsAreFiniteForIdentity)
  {
    ShapeDistortionMetric shape;
    AreaDistortionMetric area;
    ShapeSizeMetric shapeSize;

    EXPECT_NEAR(shape.value(Matrix2::Identity()), 0, 1e-14);
    EXPECT_NEAR(area.value(Matrix2::Identity()), 0, 1e-14);
    EXPECT_NEAR(shapeSize.value(Matrix2::Identity()), 0, 1e-14);
  }

  TEST(Rodin_Adaptation_TMOP, IdentityTargetWorks)
  {
    const auto target = IdentityTargetEvaluator().evaluate();
    EXPECT_NEAR((target.W - Matrix2::Identity()).norm(), 0, 1e-14);
  }

  TEST(Rodin_Adaptation_TMOP, P2UpgradeCreatesExpectedNodes)
  {
    auto geometry = upgradeUnitTriangle();

    ASSERT_EQ(geometry.nodes.size(), 6);
    ASSERT_EQ(geometry.cells.size(), 1);
    EXPECT_TRUE(geometry.nodes[0].fixed);
    EXPECT_TRUE(geometry.nodes[1].fixed);
    EXPECT_TRUE(geometry.nodes[2].fixed);

    const Index midpoint = geometry.cells.front()[3];
    EXPECT_FALSE(geometry.nodes[midpoint].fixed);
    EXPECT_NEAR(geometry.nodes[midpoint].x[0], 0.5, 1e-14);
    EXPECT_NEAR(geometry.nodes[midpoint].x[1], 0.0, 1e-14);
  }

  TEST(Rodin_Adaptation_TMOP, LinearlyInitializedP2HasLinearJacobian)
  {
    const auto geometry = upgradeUnitTriangle();
    CurvedTriangleJacobianEvaluator evaluator;

    for (const auto& point : std::vector<ReferencePoint>{
           { Real(1) / Real(3), Real(1) / Real(3) },
           { Real(0.2), Real(0.2) },
           { Real(0.6), Real(0.2) },
           { Real(0.2), Real(0.6) } })
    {
      const auto J = evaluator.jacobian(geometry, 0, point);
      EXPECT_NEAR((J - Matrix2::Identity()).norm(), 0, 1e-13);
      EXPECT_NEAR(J.determinant(), 1, 1e-13);
    }
  }

  TEST(Rodin_Adaptation_TMOP, MovingHighOrderNodeChangesJacobian)
  {
    auto geometry = upgradeUnitTriangle();
    CurvedTriangleJacobianEvaluator evaluator;
    const auto before = evaluator.jacobian(geometry, 0);

    distortMidpoint(geometry);
    const auto after = evaluator.jacobian(geometry, 0);

    EXPECT_GT((after - before).norm(), 1e-8);
  }

  TEST(Rodin_Adaptation_TMOP, InvalidCurvedGeometryIsDetected)
  {
    const auto geometry = makeInvertedGeometry();
    SquaredDistanceMetric metric;
    Objective objective(metric);

    EXPECT_LT(objective.minJacobian(geometry), 0);
    EXPECT_EQ(objective.invalidElementCount(geometry), 1);
  }

  TEST(Rodin_Adaptation_TMOP, ObjectiveIsFiniteAndIncreasesForDistortion)
  {
    auto geometry = upgradeUnitTriangle();
    SquaredDistanceMetric metric;
    Objective objective(metric);

    const Real clean = objective.value(geometry);
    distortMidpoint(geometry);
    const Real distorted = objective.value(geometry);

    EXPECT_TRUE(std::isfinite(clean));
    EXPECT_TRUE(std::isfinite(distorted));
    EXPECT_GT(distorted, clean);
  }

  TEST(Rodin_Adaptation_TMOP, DeviationPenaltyTracksNodeDisplacement)
  {
    auto geometry = upgradeUnitTriangle();
    SquaredDistanceMetric metric;
    Objective objective(metric);
    objective.setDeviationWeight(2.0);

    EXPECT_NEAR(objective.deviationValue(geometry), 0, 1e-14);
    distortMidpoint(geometry);
    EXPECT_GT(objective.deviationValue(geometry), 0);
  }

  TEST(Rodin_Adaptation_TMOP, OptimizerReducesObjective)
  {
    auto geometry = upgradeUnitTriangle();
    distortMidpoint(geometry);

    SquaredDistanceMetric metric;
    Objective objective(metric);
    objective.setDeviationWeight(0.1);
    const Real before = objective.value(geometry);

    OptimizerOptions options;
    options.maxIterations = 20;
    options.initialStepSize = 0.05;
    const auto report = Optimizer(geometry, objective)
      .setOptions(options)
      .optimize();

    EXPECT_LT(report.finalObjective, before);
    EXPECT_EQ(report.invalidElements, 0);
    EXPECT_GT(report.iterations, 0);
  }

  TEST(Rodin_Adaptation_TMOP, OptimizerDoesNotMoveFixedVertices)
  {
    auto geometry = upgradeUnitTriangle();
    const auto v0 = geometry.nodes[0].x;
    const auto v1 = geometry.nodes[1].x;
    const auto v2 = geometry.nodes[2].x;
    distortMidpoint(geometry);

    SquaredDistanceMetric metric;
    Objective objective(metric);
    Optimizer(geometry, objective).optimize();

    EXPECT_NEAR((geometry.nodes[0].x - v0).norm(), 0, 1e-14);
    EXPECT_NEAR((geometry.nodes[1].x - v1).norm(), 0, 1e-14);
    EXPECT_NEAR((geometry.nodes[2].x - v2).norm(), 0, 1e-14);
  }

  TEST(Rodin_Adaptation_TMOP, OptimizerRejectsInvalidInitialGeometry)
  {
    auto geometry = makeInvertedGeometry();
    SquaredDistanceMetric metric;
    Objective objective(metric);

    const auto report = Optimizer(geometry, objective).optimize();

    EXPECT_EQ(report.iterations, 0);
    EXPECT_EQ(report.invalidElements, 1);
    EXPECT_FALSE(report.converged);
    EXPECT_EQ(report.reason, "invalid initial geometry");
  }

  TEST(Rodin_Adaptation_TMOP, CutMeshCanBeUpgradedAndEvaluated)
  {
    auto mesh = makeTwoTriangleSquare();

    P1 space(mesh);
    GridFunction phi(space);
    phi[0] = -1;
    phi[1] = 1;
    phi[2] = -1;
    phi[3] = 1;

    auto cut = LevelSetDiscretizerTriangles(phi)
      .setSignTolerance(1e-12)
      .setInterfaceAttribute(99)
      .discretize();

    auto geometry = HighOrderGeometryUpgrade().upgrade(cut.mesh, 2);
    SquaredDistanceMetric metric;
    Objective objective(metric);

    EXPECT_GT(cut.interfaceGraph.edges.size(), 0);
    EXPECT_EQ(objective.invalidElementCount(geometry), 0);
    EXPECT_TRUE(std::isfinite(objective.value(geometry)));
  }

  struct MetricMatrixCase
  {
    Matrix2 A;
    Real expectedSquaredDistance;
  };

  class Rodin_Adaptation_TMOP_MetricMatrices
    : public testing::TestWithParam<MetricMatrixCase>
  {};

  TEST_P(Rodin_Adaptation_TMOP_MetricMatrices, SquaredDistanceMatchesFormula)
  {
    const auto c = GetParam();
    SquaredDistanceMetric metric;
    EXPECT_NEAR(metric.value(c.A), c.expectedSquaredDistance, 1e-14);
  }

  INSTANTIATE_TEST_SUITE_P(
      Values,
      Rodin_Adaptation_TMOP_MetricMatrices,
      testing::Values(
        MetricMatrixCase{(Matrix2() << 1, 0, 0, 1).finished(), 0},
        MetricMatrixCase{(Matrix2() << 2, 0, 0, 1).finished(), 0.5},
        MetricMatrixCase{(Matrix2() << 1, 0.5, 0, 1).finished(), 0.125},
        MetricMatrixCase{(Matrix2() << 1, 0, 0.25, 1).finished(), 0.03125},
        MetricMatrixCase{(Matrix2() << 1.5, 0, 0, 1.5).finished(), 0.25},
        MetricMatrixCase{(Matrix2() << 0.5, 0, 0, 0.5).finished(), 0.25},
        MetricMatrixCase{(Matrix2() << 1, 0.25, 0.25, 1).finished(), 0.0625},
        MetricMatrixCase{(Matrix2() << 0.75, 0.1, -0.2, 1.25).finished(), 0.0875}));

  class Rodin_Adaptation_TMOP_ShapeAreaMetrics
    : public testing::TestWithParam<Matrix2>
  {};

  TEST_P(Rodin_Adaptation_TMOP_ShapeAreaMetrics, ShapeSizeIsNonnegativeForValidMaps)
  {
    const auto A = GetParam();
    ShapeDistortionMetric shape;
    AreaDistortionMetric area;
    ShapeSizeMetric shapeSize;

    EXPECT_GT(A.determinant(), 0);
    EXPECT_GE(shape.value(A), -1e-14);
    EXPECT_GE(area.value(A), -1e-14);
    EXPECT_NEAR(shapeSize.value(A), shape.value(A) + area.value(A), 1e-14);
  }

  INSTANTIATE_TEST_SUITE_P(
      ValidMaps,
      Rodin_Adaptation_TMOP_ShapeAreaMetrics,
      testing::Values(
        (Matrix2() << 1, 0, 0, 1).finished(),
        (Matrix2() << 2, 0, 0, 2).finished(),
        (Matrix2() << 0.5, 0, 0, 0.5).finished(),
        (Matrix2() << 1, 0.25, 0, 1).finished(),
        (Matrix2() << 1.25, 0, 0.1, 0.9).finished(),
        (Matrix2() << 0.8, -0.2, 0.15, 1.3).finished(),
        (Matrix2() << 1.1, 0.4, -0.1, 1.0).finished(),
        (Matrix2() << 1.0, -0.3, 0.2, 1.2).finished()));

  struct AffineTriangleCase
  {
    Math::SpatialPoint a;
    Math::SpatialPoint b;
    Math::SpatialPoint c;
    Matrix2 J;
  };

  class Rodin_Adaptation_TMOP_AffineTriangles
    : public testing::TestWithParam<AffineTriangleCase>
  {};

  TEST_P(Rodin_Adaptation_TMOP_AffineTriangles, P2UpgradePreservesAffineJacobian)
  {
    const auto c = GetParam();
    auto mesh = makeTriangle(c.a, c.b, c.c);
    const auto geometry = HighOrderGeometryUpgrade().upgrade(mesh, 2);
    CurvedTriangleJacobianEvaluator evaluator;

    for (const auto& point : std::vector<ReferencePoint>{
           {0.2, 0.2}, {0.6, 0.2}, {0.2, 0.6}, {Real(1) / Real(3), Real(1) / Real(3)}})
    {
      const auto J = evaluator.jacobian(geometry, 0, point);
      EXPECT_NEAR((J - c.J).norm(), 0, 1e-13);
      EXPECT_GT(J.determinant(), 0);
    }
  }

  INSTANTIATE_TEST_SUITE_P(
      AffineMaps,
      Rodin_Adaptation_TMOP_AffineTriangles,
      testing::Values(
        AffineTriangleCase{{0, 0}, {1, 0}, {0, 1}, (Matrix2() << 1, 0, 0, 1).finished()},
        AffineTriangleCase{{1, 2}, {3, 2}, {1, 5}, (Matrix2() << 2, 0, 0, 3).finished()},
        AffineTriangleCase{{-1, 0}, {0, 0.5}, {-0.5, 2}, (Matrix2() << 1, 0.5, 0.5, 2).finished()},
        AffineTriangleCase{{0.2, 0.1}, {1.2, 0.3}, {0.4, 1.4}, (Matrix2() << 1, 0.2, 0.2, 1.3).finished()},
        AffineTriangleCase{{2, -1}, {2.5, -0.5}, {1.5, 1}, (Matrix2() << 0.5, -0.5, 0.5, 2).finished()},
        AffineTriangleCase{{-2, -1}, {-0.5, -1}, {-1.5, 0.25}, (Matrix2() << 1.5, 0.5, 0, 1.25).finished()}));

  class Rodin_Adaptation_TMOP_MidpointPerturbations
    : public testing::TestWithParam<std::pair<Index, Math::SpatialPoint>>
  {};

  TEST_P(Rodin_Adaptation_TMOP_MidpointPerturbations, PerturbingFreeNodeChangesObjective)
  {
    auto geometry = upgradeUnitTriangle();
    const auto c = GetParam();
    SquaredDistanceMetric metric;
    Objective objective(metric);
    const Real before = objective.value(geometry);

    const Index node = geometry.cells.front()[c.first];
    geometry.nodes[node].x[0] += c.second[0];
    geometry.nodes[node].x[1] += c.second[1];

    EXPECT_GT(objective.value(geometry), before);
    EXPECT_EQ(objective.invalidElementCount(geometry), 0);
  }

  INSTANTIATE_TEST_SUITE_P(
      FreeNodeMoves,
      Rodin_Adaptation_TMOP_MidpointPerturbations,
      testing::Values(
        std::pair<Index, Math::SpatialPoint>{3, Math::SpatialPoint({0.00, 0.10})},
        std::pair<Index, Math::SpatialPoint>{3, Math::SpatialPoint({0.10, 0.00})},
        std::pair<Index, Math::SpatialPoint>{4, Math::SpatialPoint({0.00, 0.10})},
        std::pair<Index, Math::SpatialPoint>{4, Math::SpatialPoint({-0.10, 0.00})},
        std::pair<Index, Math::SpatialPoint>{5, Math::SpatialPoint({0.08, 0.00})},
        std::pair<Index, Math::SpatialPoint>{5, Math::SpatialPoint({0.00, -0.08})}));

  class Rodin_Adaptation_TMOP_DeviationWeights
    : public testing::TestWithParam<Real>
  {};

  TEST_P(Rodin_Adaptation_TMOP_DeviationWeights, DeviationScalesWithWeight)
  {
    auto geometry = upgradeUnitTriangle();
    distortMidpoint(geometry);

    SquaredDistanceMetric metric;
    Objective objective(metric);
    objective.setDeviationWeight(GetParam());

    EXPECT_NEAR(objective.deviationValue(geometry), GetParam() * 0.25 * 0.25, 1e-14);
  }

  INSTANTIATE_TEST_SUITE_P(
      Weights,
      Rodin_Adaptation_TMOP_DeviationWeights,
      testing::Values(0.0, 0.1, 1.0, 10.0, 100.0));

  class Rodin_Adaptation_TMOP_OptimizerPerturbations
    : public testing::TestWithParam<std::pair<Index, Math::SpatialPoint>>
  {};

  TEST_P(Rodin_Adaptation_TMOP_OptimizerPerturbations, OptimizerReducesDifferentPerturbations)
  {
    auto geometry = upgradeUnitTriangle();
    const auto c = GetParam();
    const Index node = geometry.cells.front()[c.first];
    geometry.nodes[node].x[0] += c.second[0];
    geometry.nodes[node].x[1] += c.second[1];

    SquaredDistanceMetric metric;
    Objective objective(metric);
    objective.setDeviationWeight(0.1);
    const Real before = objective.value(geometry);

    OptimizerOptions options;
    options.maxIterations = 15;
    options.initialStepSize = 0.03;
    const auto report = Optimizer(geometry, objective)
      .setOptions(options)
      .optimize();

    EXPECT_LT(report.finalObjective, before);
    EXPECT_EQ(report.invalidElements, 0);
  }

  INSTANTIATE_TEST_SUITE_P(
      Moves,
      Rodin_Adaptation_TMOP_OptimizerPerturbations,
      testing::Values(
        std::pair<Index, Math::SpatialPoint>{3, Math::SpatialPoint({0.00, 0.20})},
        std::pair<Index, Math::SpatialPoint>{3, Math::SpatialPoint({0.15, 0.00})},
        std::pair<Index, Math::SpatialPoint>{4, Math::SpatialPoint({-0.15, 0.00})},
        std::pair<Index, Math::SpatialPoint>{4, Math::SpatialPoint({0.00, -0.12})},
        std::pair<Index, Math::SpatialPoint>{5, Math::SpatialPoint({0.12, 0.12})}));

  class Rodin_Adaptation_TMOP_LinearMeshObjective
    : public testing::TestWithParam<AffineTriangleCase>
  {};

  TEST_P(Rodin_Adaptation_TMOP_LinearMeshObjective, FunctionBaseMetricMatchesManualQuadrature)
  {
    const auto c = GetParam();
    auto mesh = makeTriangle(c.a, c.b, c.c);
    SquaredDistanceMetric metric;

    const Real objective = LinearMeshMetricObjective(metric).compute(mesh);
    const Real area = std::abs(c.J.determinant()) / Real(2);
    const Real expected = area * metric.value(c.J);

    EXPECT_NEAR(objective, expected, 1e-13);
  }

  INSTANTIATE_TEST_SUITE_P(
      FunctionBaseBridge,
      Rodin_Adaptation_TMOP_LinearMeshObjective,
      testing::Values(
        AffineTriangleCase{{0, 0}, {1, 0}, {0, 1}, (Matrix2() << 1, 0, 0, 1).finished()},
        AffineTriangleCase{{0, 0}, {2, 0}, {0, 1}, (Matrix2() << 2, 0, 0, 1).finished()},
        AffineTriangleCase{{0, 0}, {1, 0.2}, {0.1, 1.1}, (Matrix2() << 1, 0.1, 0.2, 1.1).finished()},
        AffineTriangleCase{{1, 1}, {1.5, 1}, {1, 1.5}, (Matrix2() << 0.5, 0, 0, 0.5).finished()}));

  TEST(Rodin_Adaptation_TMOP, ProblemRequiresInitialization)
  {
    auto mesh = makeUnitTriangle();
    Rodin::Adaptation::TMOP::Problem problem(mesh);

    try
    {
      (void) problem.value();
      FAIL() << "Expected uninitialized TMOP problem to throw.";
    }
    catch (const Alert::Exception& e)
    {
      EXPECT_NE(
        std::string(e.what()).find("has not been initialized"),
        std::string::npos);
    }
  }

  TEST(Rodin_Adaptation_TMOP, ProblemInitializesAndEvaluates)
  {
    auto mesh = makeUnitTriangle();
    SquaredDistanceMetric metric;
    Rodin::Adaptation::TMOP::Problem problem(mesh);
    problem
      .setMetric(metric)
      .setDeviationWeight(0.1)
      .initialize();

    EXPECT_TRUE(problem.isInitialized());
    EXPECT_EQ(problem.invalidElementCount(), 0);
    EXPECT_TRUE(std::isfinite(problem.value()));
  }

  TEST(Rodin_Adaptation_TMOP, ProblemOptimizerUpdatesGeometry)
  {
    auto mesh = makeUnitTriangle();
    SquaredDistanceMetric metric;
    Rodin::Adaptation::TMOP::Problem problem(mesh);
    problem
      .setMetric(metric)
      .setDeviationWeight(0.1)
      .initialize();
    distortMidpoint(problem.getGeometry());

    const Real before = problem.value();
    OptimizerOptions options;
    options.maxIterations = 10;
    options.initialStepSize = 0.03;
    const auto report = problem
      .setOptimizerOptions(options)
      .optimize();

    EXPECT_LT(report.finalObjective, before);
    EXPECT_EQ(report.finalInvalidElements, 0);
  }

  TEST(Rodin_Adaptation_TMOP, ProblemCanLeaveOriginalVerticesFree)
  {
    auto mesh = makeUnitTriangle();
    Rodin::Adaptation::TMOP::Problem problem(mesh);
    problem
      .setFixOriginalVertices(false)
      .initialize();

    EXPECT_FALSE(problem.getGeometry().nodes[0].fixed);
    EXPECT_FALSE(problem.getGeometry().nodes[1].fixed);
    EXPECT_FALSE(problem.getGeometry().nodes[2].fixed);
  }
}
