/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include <Rodin/Assembly/Default.h>
#include <Rodin/Adaptation.h>
#include <Rodin/Geometry/LevelSetDiscretizerTriangles.h>
#include <Rodin/Geometry.h>
#include <Rodin/Solver/NewtonSolver.h>
#include <Rodin/Solver/SparseLU.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Adaptation::TargetMatrixOptimization;
using namespace Rodin::Geometry;
using namespace Rodin::Solver;
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

  static LevelSetDiscretizerTrianglesResult cutTwoTriangleSquare()
  {
    auto mesh = makeTwoTriangleSquare();

    P1 space(mesh);
    GridFunction phi(space);
    phi[0] = -1;
    phi[1] = 1;
    phi[2] = -1;
    phi[3] = 1;

    return LevelSetDiscretizerTriangles(phi)
      .setSignTolerance(1e-12)
      .setInterfaceAttribute(99)
      .setNegativeCellAttribute(1)
      .setPositiveCellAttribute(2)
      .discretize();
  }

  struct MeshQualitySummary
  {
    Real minArea = std::numeric_limits<Real>::infinity();
    Real minQuality = std::numeric_limits<Real>::infinity();
    Index invertedCells = 0;
    Index degenerateCells = 0;
  };

  static Real signedTriangleArea(
      const Math::SpatialPoint& a,
      const Math::SpatialPoint& b,
      const Math::SpatialPoint& c)
  {
    return Real(0.5)
      * ((b[0] - a[0]) * (c[1] - a[1])
       - (b[1] - a[1]) * (c[0] - a[0]));
  }

  static Real triangleQuality(
      const Math::SpatialPoint& a,
      const Math::SpatialPoint& b,
      const Math::SpatialPoint& c)
  {
    const Real area = std::abs(signedTriangleArea(a, b, c));
    const Real l0 = (b - a).squaredNorm();
    const Real l1 = (c - b).squaredNorm();
    const Real l2 = (a - c).squaredNorm();
    const Real denom = l0 + l1 + l2;
    if (denom <= Real(0))
      return Real(0);
    return Real(4) * std::sqrt(Real(3)) * area / denom;
  }

  static MeshQualitySummary summarizeMesh(const LocalMesh& mesh)
  {
    MeshQualitySummary summary;
    const auto& conn = mesh.getConnectivity();
    for (Index c = 0; c < mesh.getCellCount(); ++c)
    {
      if (conn.getGeometry(2, c) != Polytope::Type::Triangle)
        continue;
      const auto& cell = conn.getPolytope(2, c);
      const auto x0 = mesh.getVertexCoordinates(cell(0));
      const auto x1 = mesh.getVertexCoordinates(cell(1));
      const auto x2 = mesh.getVertexCoordinates(cell(2));
      const Real area = signedTriangleArea(x0, x1, x2);
      const Real absArea = std::abs(area);
      summary.minArea = std::min(summary.minArea, absArea);
      summary.minQuality =
        std::min(summary.minQuality, triangleQuality(x0, x1, x2));
      if (area <= Real(0))
        summary.invertedCells++;
      if (absArea <= Real(1e-14))
        summary.degenerateCells++;
    }
    if (!std::isfinite(summary.minArea))
    {
      summary.minArea = 0;
      summary.minQuality = 0;
    }
    return summary;
  }

  static Real strictTMOPQualityEnergy(
      LocalMesh& mesh,
      const InitialElementTargetJacobian& target)
  {
    P1 space(mesh, 2);
    GridFunction displacement(space);
    displacement.getData().setZero();
    QualityTerm quality(SquaredDistanceMetric{}, target);
    return quality.energy(displacement);
  }

  static LocalMesh makeUniformSquare(size_t resolution)
  {
    auto mesh =
      LocalMesh::UniformGrid(Polytope::Type::Triangle, {resolution, resolution});
    mesh.scale(Real(1) / static_cast<Real>(resolution - 1));
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    return mesh;
  }

  using LevelSetFunction = std::function<Real(const Math::SpatialPoint&)>;

  struct LevelSetCutCase
  {
    std::string name;
    size_t resolution;
    LevelSetFunction phi;
  };

  static LevelSetDiscretizerTrianglesResult cutLevelSetCase(
      const LevelSetCutCase& c)
  {
    auto mesh = makeUniformSquare(c.resolution);
    P1 space(mesh);
    GridFunction phi(space);
    for (Index v = 0; v < mesh.getVertexCount(); ++v)
      phi[v] = c.phi(mesh.getVertexCoordinates(v));

    return LevelSetDiscretizerTriangles(phi)
      .setSignTolerance(1e-12)
      .setInterfaceAttribute(99)
      .setNegativeCellAttribute(1)
      .setPositiveCellAttribute(2)
      .discretize();
  }

  static std::vector<LevelSetCutCase> simpleCutCases()
  {
    return {
      {"vertical-line", 9,
        [](const Math::SpatialPoint& x) { return x[0] - Real(0.47); }},
      {"diagonal-line", 9,
        [](const Math::SpatialPoint& x) { return x[0] + x[1] - Real(0.93); }},
      {"shallow-line", 10,
        [](const Math::SpatialPoint& x) {
          return x[1] - (Real(0.28) + Real(0.37) * x[0]);
        }},
      {"circle", 11,
        [](const Math::SpatialPoint& x) {
          const Real dx = x[0] - Real(0.52);
          const Real dy = x[1] - Real(0.49);
          return std::sqrt(dx * dx + dy * dy) - Real(0.27);
        }},
      {"ellipse", 12,
        [](const Math::SpatialPoint& x) {
          const Real dx = (x[0] - Real(0.48)) / Real(0.31);
          const Real dy = (x[1] - Real(0.53)) / Real(0.19);
          return dx * dx + dy * dy - Real(1);
        }},
      {"sine-curve", 12,
        [](const Math::SpatialPoint& x) {
          return x[1]
            - (Real(0.50)
             + Real(0.08) * std::sin(Real(2) * Math::Constants::pi() * x[0]));
        }}
    };
  }

  static std::vector<LevelSetCutCase> complicatedCutCases()
  {
    return {
      {"two-disjoint-circles", 13,
        [](const Math::SpatialPoint& x) {
          const Real dx0 = x[0] - Real(0.34);
          const Real dy0 = x[1] - Real(0.44);
          const Real dx1 = x[0] - Real(0.68);
          const Real dy1 = x[1] - Real(0.61);
          const Real d0 = std::sqrt(dx0 * dx0 + dy0 * dy0) - Real(0.16);
          const Real d1 = std::sqrt(dx1 * dx1 + dy1 * dy1) - Real(0.13);
          return std::min(d0, d1);
        }},
      {"annulus", 14,
        [](const Math::SpatialPoint& x) {
          const Real dx = x[0] - Real(0.50);
          const Real dy = x[1] - Real(0.50);
          return std::abs(std::sqrt(dx * dx + dy * dy) - Real(0.27))
               - Real(0.055);
        }},
      {"wavy-circle", 14,
        [](const Math::SpatialPoint& x) {
          const Real dx = x[0] - Real(0.51);
          const Real dy = x[1] - Real(0.48);
          const Real theta = std::atan2(dy, dx);
          const Real radius = Real(0.25) + Real(0.035) * std::sin(Real(5) * theta);
          return std::sqrt(dx * dx + dy * dy) - radius;
        }},
      {"offset-crescent", 14,
        [](const Math::SpatialPoint& x) {
          const Real dx0 = x[0] - Real(0.48);
          const Real dy0 = x[1] - Real(0.50);
          const Real dx1 = x[0] - Real(0.57);
          const Real dy1 = x[1] - Real(0.50);
          const Real outer = std::sqrt(dx0 * dx0 + dy0 * dy0) - Real(0.29);
          const Real inner = Real(0.18) - std::sqrt(dx1 * dx1 + dy1 * dy1);
          return std::max(outer, inner);
        }}
    };
  }

  static void perturbMesh(
      LocalMesh& mesh,
      Real amplitude,
      const std::vector<bool>& protectedVertices = {})
  {
    for (Index v = 0; v < mesh.getVertexCount(); ++v)
    {
      if (!protectedVertices.empty()
          && protectedVertices[static_cast<size_t>(v)])
        continue;

      auto x = mesh.getVertexCoordinates(v);
      const bool boundary =
        x[0] <= Real(1e-12) || x[0] >= Real(1) - Real(1e-12) ||
        x[1] <= Real(1e-12) || x[1] >= Real(1) - Real(1e-12);
      if (boundary)
        continue;

      const Real sx = std::sin(Real(7) * x[0] + Real(3) * x[1]);
      const Real sy = std::cos(Real(5) * x[0] - Real(11) * x[1]);
      x[0] += amplitude * sx;
      x[1] += amplitude * sy;
      mesh.setVertexCoordinates(v, x);
    }
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
  }

  static std::vector<bool> interfaceProtectedVertices(
      const LevelSetDiscretizerTrianglesResult& cut)
  {
    std::vector<bool> protectedVertices(cut.mesh.getVertexCount(), false);
    for (Index v = 0; v < cut.report.vertexOrigins.size(); ++v)
    {
      if (cut.report.vertexOrigins[v].kind
          == OutputVertexOriginKind::InterfaceGraphVertex)
      {
        protectedVertices[static_cast<size_t>(v)] = true;
      }
    }

    const auto& conn = cut.mesh.getConnectivity();
    for (const auto& item : cut.report.interfaceEdgeProvenance)
    {
      const auto& edge = conn.getPolytope(1, item.first);
      protectedVertices[static_cast<size_t>(edge(0))] = true;
      protectedVertices[static_cast<size_t>(edge(1))] = true;
    }
    return protectedVertices;
  }

  static std::vector<bool> protectedVerticesForCase(
      const LevelSetDiscretizerTrianglesResult& cut,
      const LevelSetCutCase& c,
      Real interfaceBand)
  {
    auto protectedVertices = interfaceProtectedVertices(cut);
    for (Index v = 0; v < cut.mesh.getVertexCount(); ++v)
    {
      const auto x = cut.mesh.getVertexCoordinates(v);
      const bool boundary =
        x[0] <= Real(1e-12) || x[0] >= Real(1) - Real(1e-12) ||
        x[1] <= Real(1e-12) || x[1] >= Real(1) - Real(1e-12);
      if (boundary || std::abs(c.phi(x)) <= interfaceBand)
        protectedVertices[static_cast<size_t>(v)] = true;
    }
    return protectedVertices;
  }

  template <class FES, class Data>
  static void applyDisplacement(
      LocalMesh& mesh,
      const GridFunction<FES, Data>& displacement,
      Real scale)
  {
    const auto& data = displacement.getData();
    const Index vertexCount = static_cast<Index>(mesh.getVertexCount());
    const size_t sdim = mesh.getSpaceDimension();
    for (Index v = 0; v < vertexCount; ++v)
    {
      auto x = mesh.getVertexCoordinates(v);
      for (size_t c = 0; c < sdim; ++c)
        x[c] += scale * data(v + static_cast<Index>(c) * vertexCount);
      mesh.setVertexCoordinates(v, x);
    }
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
  }

  static void fillDisplacement(Eigen::VectorXd& data, Real amplitude)
  {
    for (Eigen::Index i = 0; i < data.size(); ++i)
    {
      const Real k = static_cast<Real>(i + 1);
      data(i) = amplitude * (std::sin(Real(0.37) * k)
               + Real(0.5) * std::cos(Real(0.19) * k));
    }
  }

  static void fillDirection(Eigen::VectorXd& data)
  {
    for (Eigen::Index i = 0; i < data.size(); ++i)
    {
      const Real k = static_cast<Real>(i + 1);
      data(i) = std::cos(Real(0.23) * k)
              + Real(0.25) * std::sin(Real(0.71) * k);
    }
    data /= data.norm();
  }

  static Real finiteDifferenceQualityDeviationError(
      LocalMesh& mesh,
      Real displacementAmplitude,
      Real qualityWeight,
      Real deviationWeight)
  {
    constexpr size_t vdim = 2;
    P1 space(mesh, vdim);
    GridFunction displacement(space);
    fillDisplacement(displacement.getData(), displacementAmplitude);

    TrialFunction du(space);
    TestFunction v(space);
    QualityTerm quality(qualityWeight);
    DeviationTerm deviation(deviationWeight);

    auto assembleResidual = [&]()
    {
      LinearForm residual(v);
      residual =
          quality.residual(displacement, v)
        + deviation.residual(displacement, v);
      residual.assemble();
      return residual.getVector();
    };

    BilinearForm tangent(du, v);
    tangent =
        quality.tangent(displacement, du, v)
      + deviation.tangent(du, v);
    tangent.assemble();

    const auto residual0 = assembleResidual();
    auto direction = displacement.getData();
    fillDirection(direction);

    const auto original = displacement.getData();
    const Real eps = 1e-7;
    displacement.getData() = original + eps * direction;
    const auto residual1 = assembleResidual();
    displacement.getData() = original;

    const auto fd = (residual1 - residual0) / eps;
    const auto jd = tangent.getOperator() * direction;
    const Real denom = std::max<Real>({Real(1), fd.norm(), jd.norm()});
    return (fd - jd).norm() / denom;
  }

  static Real finiteDifferenceFullTmopError(
      const LevelSetDiscretizerTrianglesResult& cut,
      Real displacementAmplitude,
      Real qualityWeight,
      Real deviationWeight,
      Real fitWeight)
  {
    constexpr size_t vdim = 2;
    P1 space(cut.mesh, vdim);
    GridFunction displacement(space);
    fillDisplacement(displacement.getData(), displacementAmplitude);

    TrialFunction du(space);
    TestFunction v(space);
    QualityTerm quality(qualityWeight);
    DeviationTerm deviation(deviationWeight);
    LevelSetFitTerm fit(cut.interfaceGraph, cut.report, fitWeight);

    auto assembleResidual = [&]()
    {
      LinearForm residual(v);
      residual =
          quality.residual(displacement, v)
        + deviation.residual(displacement, v)
        + fit.residual(displacement, v);
      residual.assemble();
      return residual.getVector();
    };

    BilinearForm tangent(du, v);
    tangent =
        quality.tangent(displacement, du, v)
      + deviation.tangent(du, v)
      + fit.tangent(displacement, du, v);
    tangent.assemble();

    const auto residual0 = assembleResidual();
    auto direction = displacement.getData();
    fillDirection(direction);

    const auto original = displacement.getData();
    const Real eps = 1e-7;
    displacement.getData() = original + eps * direction;
    const auto residual1 = assembleResidual();
    displacement.getData() = original;

    const auto fd = (residual1 - residual0) / eps;
    const auto jd = tangent.getOperator() * direction;
    const Real denom = std::max<Real>({Real(1), fd.norm(), jd.norm()});
    return (fd - jd).norm() / denom;
  }

  struct NativeTmopRun
  {
    bool accepted = false;
    Real scale = 1;
    MeshQualitySummary beforeQuality;
    MeshQualitySummary afterQuality;
    Real beforeEnergy = 0;
    Real afterEnergy = 0;
    Real sourceFitBefore = 0;
    Real sourceFitAfter = 0;
    Real initialResidual = 0;
    Real finalResidual = 0;
  };

  static NativeTmopRun runNativeTmop(
      LocalMesh& mesh,
      const LocalMesh& targetMesh,
      const InterfaceGraph& graph,
      const LevelSetDiscretizerTrianglesReport& report)
  {
    NativeTmopRun run;
    run.beforeQuality = summarizeMesh(mesh);
    InitialElementTargetJacobian target(targetMesh);
    run.beforeEnergy = strictTMOPQualityEnergy(mesh, target);

    constexpr size_t vdim = 2;
    P1 space(mesh, vdim);
    GridFunction displacement(space);
    displacement.getData().setZero();

    TrialFunction du(space);
    TestFunction v(space);
    QualityTerm quality(SquaredDistanceMetric{}, target, 1.0);
    DeviationTerm deviation(1.0);
    LevelSetFitTerm fit(graph, report, 1.0);
    LevelSetFitTerm fitDiagnostic(graph, report, 1.0);
    run.sourceFitBefore = fitDiagnostic.sourceSegmentDistanceEnergy(mesh);

    Variational::Problem tangentProblem(du, v);
    tangentProblem =
        quality.tangent(displacement, du, v)
      + deviation.tangent(du, v)
      + fit.tangent(displacement, du, v)
      + quality.residual(displacement, v)
      + deviation.residual(displacement, v)
      + fit.residual(displacement, v);

    SparseLU linearSolver{tangentProblem};
    NewtonSolver newton(linearSolver);
    // Deliberately full-step: these tests should fail if the model only works
    // with line search/damping instead of the Rodin Newton step itself.
    newton
      .setMaxIterations(1)
      .setDampingFactor(1.0)
      .setAbsoluteTolerance(1e-12)
      .setRelativeTolerance(1e-10);
    newton.solve(displacement);
    const auto newtonReport = newton.getReport();
    run.initialResidual = newtonReport.initial_residual;
    run.finalResidual = newtonReport.final_residual;

    const Real sourceFitLimit =
      std::max(Real(1e-8), Real(1.01) * run.sourceFitBefore + Real(1e-10));

    applyDisplacement(mesh, displacement, Real(1));
    run.afterQuality = summarizeMesh(mesh);
    run.afterEnergy = strictTMOPQualityEnergy(mesh, target);
    run.sourceFitAfter = fitDiagnostic.sourceSegmentDistanceEnergy(mesh);
    run.accepted =
      run.afterQuality.invertedCells == 0
      && run.afterQuality.degenerateCells == 0
      && std::isfinite(run.afterEnergy)
      && std::isfinite(run.sourceFitAfter)
      && run.sourceFitAfter <= sourceFitLimit
      && run.afterEnergy <= run.beforeEnergy * (Real(1) + Real(1e-6));
    return run;
  }

  static void expectCutHasUsableProvenance(
      const LevelSetDiscretizerTrianglesResult& cut)
  {
    ASSERT_GT(cut.interfaceGraph.edges.size(), 0);
    EXPECT_EQ(
      cut.report.interfaceEdgeProvenance.size(),
      cut.interfaceGraph.edges.size());
    EXPECT_EQ(cut.report.vertexOrigins.size(), cut.mesh.getVertexCount());
    EXPECT_EQ(cut.report.cellProvenance.size(), cut.mesh.getCellCount());

    const auto quality = summarizeMesh(cut.mesh);
    EXPECT_EQ(quality.invertedCells, 0);
    EXPECT_EQ(quality.degenerateCells, 0);

    LevelSetFitTerm fit(cut.interfaceGraph, cut.report);
    EXPECT_NEAR(fit.sourceSegmentDistanceEnergy(cut.mesh), 0, 1e-24);
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

  TEST(Rodin_Adaptation_TargetMatrixOptimization, SquaredDistanceMetricIsZeroAtIdentity)
  {
    SquaredDistanceMetric metric;
    EXPECT_NEAR(metric.value(Matrix2::Identity()), 0, 1e-14);
    EXPECT_NEAR(metric.gradient(Matrix2::Identity()).norm(), 0, 1e-14);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, DistortedMatrixHasLargerMetric)
  {
    SquaredDistanceMetric metric;
    Matrix2 A = Matrix2::Identity();
    A(0, 1) = 0.5;

    EXPECT_GT(metric.value(A), metric.value(Matrix2::Identity()));
    EXPECT_NEAR((metric.gradient(A) - (A - Matrix2::Identity())).norm(), 0, 1e-14);
    EXPECT_NEAR(
        (metric.hessianAction(A, Matrix2::Identity()) - Matrix2::Identity()).norm(),
        0,
        1e-14);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, SquaredDistanceMetricFollowsFormulaForSingularAndInvertedMatrices)
  {
    SquaredDistanceMetric metric;
    Matrix2 inverted = Matrix2::Identity();
    inverted(0, 0) = -1;

    Matrix2 nearSingular = Matrix2::Identity();
    nearSingular(1, 1) = 1e-14;

    EXPECT_NEAR(metric.value(inverted), 2, 1e-14);
    EXPECT_NEAR(metric.value(nearSingular), 0.5, 1e-13);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, ShapeAndAreaMetricsAreFiniteForIdentity)
  {
    ShapeDistortionMetric shape;
    AreaDistortionMetric area;
    ShapeSizeMetric shapeSize;

    EXPECT_NEAR(shape.value(Matrix2::Identity()), 0, 1e-14);
    EXPECT_NEAR(area.value(Matrix2::Identity()), 0, 1e-14);
    EXPECT_NEAR(shapeSize.value(Matrix2::Identity()), 0, 1e-14);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, IdentityTargetWorks)
  {
    auto mesh = makeUnitTriangle();
    auto cellIterator = mesh.getPolytope(2, 0);
    const auto& cell = *cellIterator;
    const Math::SpatialPoint rc({Real(1) / Real(3), Real(1) / Real(3)});

    IdentityTargetJacobian target;
    EXPECT_NEAR((target.evaluate(cell, rc) - Matrix2::Identity()).norm(), 0, 1e-14);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, P2UpgradeCreatesExpectedNodes)
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

  TEST(Rodin_Adaptation_TargetMatrixOptimization, LinearlyInitializedP2HasLinearJacobian)
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

  TEST(Rodin_Adaptation_TargetMatrixOptimization, MovingHighOrderNodeChangesJacobian)
  {
    auto geometry = upgradeUnitTriangle();
    CurvedTriangleJacobianEvaluator evaluator;
    const auto before = evaluator.jacobian(geometry, 0);

    distortMidpoint(geometry);
    const auto after = evaluator.jacobian(geometry, 0);

    EXPECT_GT((after - before).norm(), 1e-8);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, InvalidCurvedGeometryIsDetected)
  {
    const auto geometry = makeInvertedGeometry();
    SquaredDistanceMetric metric;
    Objective objective(metric);

    EXPECT_LT(objective.minJacobian(geometry), 0);
    EXPECT_EQ(objective.invalidElementCount(geometry), 1);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, ObjectiveIsFiniteAndIncreasesForDistortion)
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

  TEST(Rodin_Adaptation_TargetMatrixOptimization, DeviationPenaltyTracksNodeDisplacement)
  {
    auto geometry = upgradeUnitTriangle();
    SquaredDistanceMetric metric;
    Objective objective(metric);
    objective.setDeviationWeight(2.0);

    EXPECT_NEAR(objective.deviationValue(geometry), 0, 1e-14);
    distortMidpoint(geometry);
    EXPECT_GT(objective.deviationValue(geometry), 0);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, OptimizerReducesObjective)
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

  TEST(Rodin_Adaptation_TargetMatrixOptimization, OptimizerDoesNotMoveFixedVertices)
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

  TEST(Rodin_Adaptation_TargetMatrixOptimization, OptimizerRejectsInvalidInitialGeometry)
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

  TEST(Rodin_Adaptation_TargetMatrixOptimization, CutMeshCanBeUpgradedAndEvaluated)
  {
    auto cut = cutTwoTriangleSquare();
    auto geometry = HighOrderGeometryUpgrade().upgrade(cut.mesh, 2);
    SquaredDistanceMetric metric;
    Objective objective(metric);

    EXPECT_GT(cut.interfaceGraph.edges.size(), 0);
    EXPECT_EQ(objective.invalidElementCount(geometry), 0);
    EXPECT_TRUE(std::isfinite(objective.value(geometry)));
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, DeviationTermAssemblesResidualAndTangent)
  {
    auto mesh = makeTwoTriangleSquare();
    constexpr size_t vdim = 2;
    P1 space(mesh, vdim);

    GridFunction displacement(space);
    auto& data = displacement.getData();
    for (Eigen::Index i = 0; i < data.size(); ++i)
      data(i) = static_cast<Real>(i + 1) / static_cast<Real>(data.size());

    TrialFunction du(space);
    TestFunction v(space);
    DeviationTerm term;

    BilinearForm tangent(du, v);
    tangent = term.tangent(du, v);
    tangent.assemble();

    LinearForm residual(v);
    residual = term.residual(displacement, v);
    residual.assemble();

    const auto mismatch =
      tangent.getOperator() * displacement.getData() - residual.getVector();
    EXPECT_NEAR(mismatch.norm(), 0, 1e-13);
    EXPECT_GT(residual.getVector().dot(displacement.getData()), 0);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, QualityTermResidualAndTangentAreConsistent)
  {
    auto mesh = makeTwoTriangleSquare();
    constexpr size_t vdim = 2;
    P1 space(mesh, vdim);

    GridFunction displacement(space);
    auto& data = displacement.getData();
    for (Eigen::Index i = 0; i < data.size(); ++i)
      data(i) = 0.01 * std::sin(static_cast<Real>(i + 1));

    TrialFunction du(space);
    TestFunction v(space);
    QualityTerm quality;

    auto assembleResidual = [&]()
    {
      LinearForm residual(v);
      residual = quality.residual(displacement, v);
      residual.assemble();
      return residual.getVector();
    };

    BilinearForm tangent(du, v);
    tangent = quality.tangent(displacement, du, v);
    tangent.assemble();

    const auto residual0 = assembleResidual();
    auto direction = data;
    for (Eigen::Index i = 0; i < direction.size(); ++i)
      direction(i) = std::cos(static_cast<Real>(i + 1));
    direction /= direction.norm();

    const auto original = data;
    const Real eps = 1e-7;
    data = original + eps * direction;
    const auto residual1 = assembleResidual();
    data = original;

    const auto fd = (residual1 - residual0) / eps;
    const auto jd = tangent.getOperator() * direction;

    EXPECT_NEAR((fd - jd).norm(), 0, 1e-7);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, QualityTermResidualAndTangentAreConsistentForH1P2)
  {
    auto mesh = makeTwoTriangleSquare();
    VectorH1<2, LocalMesh> space(std::integral_constant<size_t, 2>{}, mesh, 2);

    GridFunction displacement(space);
    auto& data = displacement.getData();
    for (Eigen::Index i = 0; i < data.size(); ++i)
      data(i) = 0.003 * std::sin(static_cast<Real>(i + 1));

    TrialFunction du(space);
    TestFunction v(space);
    QualityTerm quality(SquaredDistanceMetric{}, IdentityTargetJacobian{}, 1.3);

    auto assembleResidual = [&]()
    {
      LinearForm residual(v);
      residual = quality.residual(displacement, v);
      residual.assemble();
      return residual.getVector();
    };

    BilinearForm tangent(du, v);
    tangent = quality.tangent(displacement, du, v);
    tangent.assemble();

    const auto residual0 = assembleResidual();
    auto direction = data;
    for (Eigen::Index i = 0; i < direction.size(); ++i)
      direction(i) = std::cos(Real(0.61) * static_cast<Real>(i + 1));
    direction /= direction.norm();

    const auto original = data;
    const Real eps = 1e-7;
    data = original + eps * direction;
    const auto residual1 = assembleResidual();
    data = original;

    const auto fd = (residual1 - residual0) / eps;
    const auto jd = tangent.getOperator() * direction;
    const Real denom = std::max<Real>({Real(1), fd.norm(), jd.norm()});

    EXPECT_LT((fd - jd).norm() / denom, 1e-7);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, LevelSetFitTermIsZeroOnFreshCutInterface)
  {
    const auto cut = cutTwoTriangleSquare();
    LevelSetFitTerm fit(cut.interfaceGraph, cut.report);

    EXPECT_GT(cut.report.interfaceEdgeProvenance.size(), 0);
    EXPECT_NEAR(fit.sourceSegmentDistanceEnergy(cut.mesh), 0, 1e-28);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, LevelSetFitTermIncreasesWhenInterfaceIsPerturbed)
  {
    auto cut = cutTwoTriangleSquare();
    LevelSetFitTerm fit(cut.interfaceGraph, cut.report);

    ASSERT_GT(cut.report.interfaceEdgeProvenance.size(), 0);
    const Index outputEdge = cut.report.interfaceEdgeProvenance.begin()->first;
    const auto& edge = cut.mesh.getConnectivity().getPolytope(1, outputEdge);
    auto x = cut.mesh.getVertexCoordinates(edge(0));
    x[0] += 0.05;
    cut.mesh.setVertexCoordinates(edge(0), x);

    EXPECT_GT(fit.sourceSegmentDistanceEnergy(cut.mesh), 1e-8);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, LevelSetFitTermResidualAndTangentAreConsistent)
  {
    const auto cut = cutTwoTriangleSquare();
    constexpr size_t vdim = 2;
    P1 space(cut.mesh, vdim);

    GridFunction displacement(space);
    auto& data = displacement.getData();
    for (Eigen::Index i = 0; i < data.size(); ++i)
      data(i) = 0.01 * std::sin(static_cast<Real>(i + 1));

    TrialFunction du(space);
    TestFunction v(space);
    LevelSetFitTerm fit(cut.interfaceGraph, cut.report, 3.0);

    auto assembleResidual = [&]()
    {
      LinearForm residual(v);
      residual = fit.residual(displacement, v);
      residual.assemble();
      return residual.getVector();
    };

    BilinearForm tangent(du, v);
    tangent = fit.tangent(displacement, du, v);
    tangent.assemble();

    const auto residual0 = assembleResidual();
    auto direction = data;
    for (Eigen::Index i = 0; i < direction.size(); ++i)
      direction(i) = std::cos(static_cast<Real>(i + 1));
    direction /= direction.norm();

    const auto original = data;
    const Real eps = 1e-7;
    data = original + eps * direction;
    const auto residual1 = assembleResidual();
    data = original;

    const auto fd = (residual1 - residual0) / eps;
    const auto jd = tangent.getOperator() * direction;
    const Real denom = std::max<Real>({Real(1), fd.norm(), jd.norm()});

    EXPECT_NEAR((fd - jd).norm() / denom, 0, 1e-7);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, CutMeshZeroDisplacementTermSmoke)
  {
    const auto cut = cutTwoTriangleSquare();
    constexpr size_t vdim = 2;
    P1 space(cut.mesh, vdim);
    GridFunction displacement(space);
    displacement = VectorFunction{ Zero(), Zero() };

    DeviationTerm deviation;
    LevelSetFitTerm fit(cut.interfaceGraph, cut.report);
    TestFunction v(space);
    LinearForm residual(v);
    residual = deviation.residual(displacement, v);
    residual.assemble();

    EXPECT_NEAR(residual.getVector().norm(), 0, 1e-14);
    EXPECT_NEAR(fit.sourceSegmentDistanceEnergy(cut.mesh), 0, 1e-28);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, NativeTermsFiniteDifferenceOnAffineMeshes)
  {
    std::vector<std::pair<std::string, LocalMesh>> meshes;
    meshes.emplace_back("unit-triangle", makeUnitTriangle());
    meshes.emplace_back("two-triangle-square", makeTwoTriangleSquare());
    meshes.emplace_back("skew-triangle", makeTriangle({0, 0}, {1.2, 0.1}, {0.2, 0.9}));
    meshes.emplace_back("large-skew-triangle", makeTriangle({0.1, 0.2}, {1.7, 0.4}, {0.3, 1.6}));

    for (auto& item : meshes)
    {
      SCOPED_TRACE(item.first);
      const Real error =
        finiteDifferenceQualityDeviationError(item.second, 0.013, 1.7, 4.0);
      EXPECT_LT(error, 1e-7) << "fd_error=" << error;
    }
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, NativeTermsFiniteDifferenceOnPerturbedAffineMeshes)
  {
    auto mesh = makeUniformSquare(6);
    perturbMesh(mesh, 0.015);

    const auto quality = summarizeMesh(mesh);
    ASSERT_EQ(quality.invertedCells, 0);
    ASSERT_EQ(quality.degenerateCells, 0);

    const Real error =
      finiteDifferenceQualityDeviationError(mesh, 0.01, 2.0, 0.75);
    EXPECT_LT(error, 1e-7);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, SimpleCutsHaveUsableProvenanceAndZeroSourceFit)
  {
    for (const auto& c : simpleCutCases())
    {
      SCOPED_TRACE(c.name);
      const auto cut = cutLevelSetCase(c);
      expectCutHasUsableProvenance(cut);
    }
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, ComplicatedCutsHaveUsableProvenanceAndZeroSourceFit)
  {
    for (const auto& c : complicatedCutCases())
    {
      SCOPED_TRACE(c.name);
      const auto cut = cutLevelSetCase(c);
      expectCutHasUsableProvenance(cut);
    }
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, FullTMOPFiniteDifferenceOnSimpleCuts)
  {
    for (const auto& c : simpleCutCases())
    {
      SCOPED_TRACE(c.name);
      const auto cut = cutLevelSetCase(c);
      const Real error =
        finiteDifferenceFullTmopError(cut, 0.004, 1.0, 8.0, 200.0);
      EXPECT_LT(error, 2e-6);
    }
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, FullTMOPFiniteDifferenceOnComplicatedCuts)
  {
    for (const auto& c : complicatedCutCases())
    {
      SCOPED_TRACE(c.name);
      const auto cut = cutLevelSetCase(c);
      const Real error =
        finiteDifferenceFullTmopError(cut, 0.003, 1.0, 10.0, 300.0);
      EXPECT_LT(error, 5e-6);
    }
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, NativeTMOPImprovesPerturbedSimpleCuts)
  {
    for (const auto& c : simpleCutCases())
    {
      SCOPED_TRACE(c.name);
      auto cut = cutLevelSetCase(c);
      auto mesh = cut.mesh;
      perturbMesh(mesh, 0.00075, protectedVerticesForCase(cut, c, 0.12));

      const auto perturbed = summarizeMesh(mesh);
      ASSERT_EQ(perturbed.invertedCells, 0);
      ASSERT_EQ(perturbed.degenerateCells, 0);

      const auto run =
        runNativeTmop(mesh, cut.mesh, cut.interfaceGraph, cut.report);
      EXPECT_TRUE(run.accepted);
      EXPECT_GT(run.scale, 0);
      EXPECT_EQ(run.afterQuality.invertedCells, 0);
      EXPECT_EQ(run.afterQuality.degenerateCells, 0);
      EXPECT_LE(run.afterEnergy, run.beforeEnergy * (Real(1) + Real(1e-8)));
      EXPECT_LE(run.sourceFitAfter, std::max(Real(1e-8), run.sourceFitBefore + Real(1e-10)));
      EXPECT_GT(run.initialResidual, 0);
      EXPECT_LE(run.finalResidual, run.initialResidual);
    }
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, NativeTMOPImprovesPerturbedComplicatedCuts)
  {
    for (const auto& c : complicatedCutCases())
    {
      SCOPED_TRACE(c.name);
      auto cut = cutLevelSetCase(c);
      auto mesh = cut.mesh;
      perturbMesh(mesh, 0.0005, protectedVerticesForCase(cut, c, 0.14));

      const auto perturbed = summarizeMesh(mesh);
      ASSERT_EQ(perturbed.invertedCells, 0);
      ASSERT_EQ(perturbed.degenerateCells, 0);

      const auto run =
        runNativeTmop(mesh, cut.mesh, cut.interfaceGraph, cut.report);
      EXPECT_TRUE(run.accepted);
      EXPECT_GT(run.scale, 0);
      EXPECT_EQ(run.afterQuality.invertedCells, 0);
      EXPECT_EQ(run.afterQuality.degenerateCells, 0);
      EXPECT_LE(run.afterEnergy, run.beforeEnergy * (Real(1) + Real(1e-8)));
      EXPECT_LE(run.sourceFitAfter, std::max(Real(1e-8), run.sourceFitBefore + Real(1e-10)));
      EXPECT_GT(run.initialResidual, 0);
      EXPECT_LE(run.finalResidual, run.initialResidual);
    }
  }

  // Strict target-matrix TMOP proof on the production path: residual/tangent
  // terms -> Variational::Problem -> SparseLU -> NewtonSolver, with the
  // identity target (W = I, T = A, mu(T) = 1/2||A - I||^2). A single interior
  // node is displaced; Newton must recover most of the distortion-induced
  // strict TMOP energy and must not regress worst-element quality or invert
  // cells. This is the test that demonstrates the method actually optimizes
  // mesh quality, not just that residual/tangent are FD-consistent.
  TEST(Rodin_Adaptation_TargetMatrixOptimization, StrictIdentityTargetNewtonReducesEnergyAndImprovesQuality)
  {
    auto strictEnergy = [](LocalMesh& m)
    {
      P1 s(m, 2);
      GridFunction z(s);
      z.getData().setZero();
      QualityTerm q(SquaredDistanceMetric{}, IdentityTargetJacobian{});
      return q.energy(z);
    };

    auto makeGrid = []()
    {
      // Raw uniform grid (integer coordinates 0..3): unit-scale triangles so
      // the identity target W = I is an appropriately scaled ideal element.
      auto m = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
      m.getConnectivity().compute(2, 1);
      m.getConnectivity().compute(1, 0);
      return m;
    };

    auto displaceCenter = [](LocalMesh& m) -> Index
    {
      for (Index v = 0; v < m.getVertexCount(); ++v)
      {
        const auto x = m.getVertexCoordinates(v);
        if (std::abs(x[0] - Real(1)) < 1e-9 && std::abs(x[1] - Real(1)) < 1e-9)
        {
          m.setVertexCoordinates(v, { Real(1.42), Real(1.31) });
          m.getConnectivity().compute(2, 1);
          m.getConnectivity().compute(1, 0);
          return v;
        }
      }
      return m.getVertexCount();
    };

    auto baselineMesh = makeGrid();
    const Real baselineEnergy = strictEnergy(baselineMesh);

    auto mesh = makeGrid();
    const Index moved = displaceCenter(mesh);
    ASSERT_LT(moved, mesh.getVertexCount());

    const auto beforeQuality = summarizeMesh(mesh);
    const Real beforeEnergy = strictEnergy(mesh);
    ASSERT_EQ(beforeQuality.invertedCells, 0);
    // The displacement genuinely raised the strict TMOP energy above the
    // (topology-fixed) irreducible baseline.
    ASSERT_GT(beforeEnergy, baselineEnergy * (Real(1) + Real(1e-6)));

    P1 space(mesh, 2);
    GridFunction displacement(space);
    displacement.getData().setZero();

    TrialFunction du(space);
    TestFunction v(space);
    QualityTerm quality(SquaredDistanceMetric{}, IdentityTargetJacobian{}, 1.0);
    DeviationTerm deviation(1e-3);

    Variational::Problem problem(du, v);
    problem =
        quality.tangent(displacement, du, v)
      + deviation.tangent(du, v)
      + quality.residual(displacement, v)
      + deviation.residual(displacement, v);

    SparseLU linearSolver{problem};
    NewtonSolver newton(linearSolver);
    newton
      .setMaxIterations(40)
      .setDampingFactor(1.0)
      .setAbsoluteTolerance(1e-12)
      .setRelativeTolerance(1e-10)
      .setStepTolerance(1e-12);
    newton.solve(displacement);
    const auto report = newton.getReport();

    applyDisplacement(mesh, displacement, Real(1));
    const auto afterQuality = summarizeMesh(mesh);
    const Real afterEnergy = strictEnergy(mesh);

    EXPECT_GT(report.initial_residual, 0);
    EXPECT_LE(report.final_residual, report.initial_residual);
    EXPECT_EQ(afterQuality.invertedCells, 0);
    EXPECT_EQ(afterQuality.degenerateCells, 0);
    // Real TMOP: the target-matrix quality energy must drop, recovering more
    // than half of the distortion-induced excess over the irreducible
    // (topology-fixed) baseline.
    EXPECT_LT(afterEnergy, beforeEnergy);
    EXPECT_LT(
        afterEnergy - baselineEnergy,
        Real(0.5) * (beforeEnergy - baselineEnergy));
    // Worst-element quality must not regress.
    EXPECT_GE(afterQuality.minQuality, beforeQuality.minQuality);
  }

  // Closed-form check of the production target: det(W_K)/2 = |K_0| and W_K is
  // the equilateral reference-to-physical Jacobian.
  TEST(Rodin_Adaptation_TargetMatrixOptimization, EquilateralSameAreaTargetClosedForm)
  {
    for (Real area : { Real(0.5), Real(0.013), Real(7.25) })
    {
      Matrix2 A0;
      A0 << std::sqrt(Real(2) * area), 0, 0, std::sqrt(Real(2) * area);
      const Matrix2 W =
        EquilateralSameAreaTargetJacobian::equilateralSameArea(A0);
      EXPECT_NEAR(Real(0.5) * std::abs(W.determinant()), area, 1e-12);
      const Real l = std::sqrt(Real(4) * area / std::sqrt(Real(3)));
      EXPECT_NEAR(W(0, 0), l, 1e-12);
      EXPECT_NEAR(W(0, 1), Real(0.5) * l, 1e-12);
      EXPECT_NEAR(W(1, 0), Real(0), 1e-12);
      EXPECT_NEAR(W(1, 1), l * std::sqrt(Real(3)) / Real(2), 1e-12);
      // det(T) magnitude = det(A0)/det(W) = 1 at u = 0, the scale-consistency
      // invariant the example's detT diagnostic checks for.
      EXPECT_NEAR(std::abs(A0.determinant() / W.determinant()), Real(1), 1e-12);
    }
  }

  // The production target on the regime where W = I fails: a [0,1]-scaled
  // (fine-element) distorted mesh. EquilateralSameAreaTargetJacobian is in
  // physical element scale, so strict TMOP recovers element shape, improves
  // worst-element quality, and inverts no cells.
  TEST(Rodin_Adaptation_TargetMatrixOptimization, EquilateralSameAreaTargetImprovesQualityAtFineScaleWithoutInversion)
  {
    auto mesh = makeUniformSquare(7);
    perturbMesh(mesh, 0.012);
    const auto before = summarizeMesh(mesh);
    ASSERT_EQ(before.invertedCells, 0);
    ASSERT_EQ(before.degenerateCells, 0);

    EquilateralSameAreaTargetJacobian target(mesh);
    auto strictEnergy = [&](LocalMesh& m)
    {
      P1 s(m, 2);
      GridFunction z(s);
      z.getData().setZero();
      QualityTerm q(SquaredDistanceMetric{}, target);
      return q.energy(z);
    };
    const Real beforeEnergy = strictEnergy(mesh);
    ASSERT_GT(beforeEnergy, 0);

    P1 space(mesh, 2);
    GridFunction displacement(space);
    displacement.getData().setZero();
    TrialFunction du(space);
    TestFunction v(space);
    QualityTerm quality(SquaredDistanceMetric{}, target, 1.0);
    DeviationTerm deviation(1e-3);

    Variational::Problem problem(du, v);
    problem =
        quality.tangent(displacement, du, v)
      + deviation.tangent(du, v)
      + quality.residual(displacement, v)
      + deviation.residual(displacement, v);

    SparseLU linearSolver{problem};
    NewtonSolver newton(linearSolver);
    newton
      .setMaxIterations(40)
      .setDampingFactor(1.0)
      .setAbsoluteTolerance(1e-12)
      .setRelativeTolerance(1e-10)
      .setStepTolerance(1e-12);
    newton.solve(displacement);
    const auto report = newton.getReport();

    applyDisplacement(mesh, displacement, Real(1));
    const auto after = summarizeMesh(mesh);
    const Real afterEnergy = strictEnergy(mesh);

    EXPECT_GT(report.initial_residual, 0);
    EXPECT_LE(report.final_residual, report.initial_residual);
    EXPECT_EQ(after.invertedCells, 0);
    EXPECT_EQ(after.degenerateCells, 0);
    EXPECT_LT(afterEnergy, beforeEnergy);
    EXPECT_GT(after.minQuality, before.minQuality);
  }

  static LocalMesh makeUnitQuad()
  {
    auto mesh = LocalMesh::Builder()
      .initialize(2)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({1, 1})
      .vertex({0, 1})
      .polytope(Polytope::Type::Quadrilateral, {0, 1, 2, 3})
      .finalize();
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    return mesh;
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, ConstantTargetJacobianFactoriesProduceExpectedMatrices)
  {
    EXPECT_NEAR(
        (ConstantTargetJacobian::identity().getMatrix() - Matrix2::Identity())
          .norm(), 0, 1e-14);

    const auto us = ConstantTargetJacobian::uniformScale(2.5).getMatrix();
    EXPECT_NEAR(us(0, 0), 2.5, 1e-14);
    EXPECT_NEAR(us(1, 1), 2.5, 1e-14);
    EXPECT_NEAR(us(0, 1), 0, 1e-14);
    EXPECT_NEAR(us(1, 0), 0, 1e-14);

    const auto dg = ConstantTargetJacobian::diagonal(3, 0.5).getMatrix();
    EXPECT_NEAR(dg(0, 0), 3, 1e-14);
    EXPECT_NEAR(dg(1, 1), 0.5, 1e-14);
    EXPECT_NEAR(std::abs(dg.determinant()), 1.5, 1e-14);

    const Real theta = 0.7;
    const auto rot = ConstantTargetJacobian::rotation(theta).getMatrix();
    EXPECT_NEAR(rot.determinant(), 1, 1e-13);
    EXPECT_NEAR((rot.transpose() * rot - Matrix2::Identity()).norm(), 0, 1e-13);

    const auto st = ConstantTargetJacobian::stretch(2, 3, theta).getMatrix();
    // det(R * diag(sx,sy)) = sx * sy.
    EXPECT_NEAR(st.determinant(), 6, 1e-12);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, ConstantTargetJacobianIsPolytopeIndependent)
  {
    const auto target = ConstantTargetJacobian::diagonal(2, 3);

    auto tri = makeTriangle({0, 0}, {1, 0}, {0, 1});
    auto quad = makeUnitQuad();
    auto triIt = tri.getPolytope(2, 0);
    auto quadIt = quad.getPolytope(2, 0);

    const Math::SpatialPoint rcA({Real(0.25), Real(0.25)});
    const Math::SpatialPoint rcB({Real(0.6), Real(0.1)});
    const Matrix2 wTri = target.evaluate(*triIt, rcA);
    const Matrix2 wQuad = target.evaluate(*quadIt, rcB);

    // Same matrix regardless of polytope type or sample point.
    EXPECT_NEAR((wTri - wQuad).norm(), 0, 1e-15);
    EXPECT_NEAR((wTri - target.getMatrix()).norm(), 0, 1e-15);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, ConstantTargetJacobianEnergyMatchesFormulaOnAffineTriangle)
  {
    // Affine triangle with map J = [b-a, c-a].
    auto mesh = makeTriangle({0.2, 0.1}, {1.2, 0.3}, {0.4, 1.4});
    Matrix2 J;
    J << 1.0, 0.2, 0.2, 1.3;

    P1 space(mesh, 2);
    GridFunction displacement(space);
    displacement.getData().setZero();

    SquaredDistanceMetric metric;
    const auto target = ConstantTargetJacobian::diagonal(1.5, 0.8);
    QualityTerm quality(metric, target);
    const Real energy = quality.energy(displacement);

    const Matrix2 W = target.getMatrix();
    const Matrix2 T = J * W.inverse();
    // Reference triangle area 1/2; affine integrand is constant.
    const Real expected = Real(0.5) * W.determinant() * metric.value(T);
    EXPECT_NEAR(energy, expected, 1e-12);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, AnalyticTargetJacobianEvaluatesAtPhysicalPoint)
  {
    auto mesh = makeTriangle({0, 0}, {2, 0}, {0, 2});
    auto it = mesh.getPolytope(2, 0);
    const auto& cell = *it;

    auto f = [](const Math::SpatialPoint& x)
    {
      Matrix2 W = Matrix2::Identity();
      W(0, 0) = Real(1) + x[0];
      W(1, 1) = Real(2) + x[1];
      return W;
    };
    auto target = makeAnalyticTargetJacobian(f);

    const Math::SpatialPoint rc({Real(0.25), Real(0.5)});
    const Geometry::Point point(cell, rc);
    const auto phys = point.getPhysicalCoordinates();
    const Matrix2 expected = f(phys);
    const Matrix2 got = target.evaluate(cell, rc);
    EXPECT_NEAR((got - expected).norm(), 0, 1e-13);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, AnalyticTargetJacobianResidualTangentFiniteDifferenceConsistent)
  {
    auto mesh = makeUniformSquare(6);
    perturbMesh(mesh, 0.01);

    P1 space(mesh, 2);
    GridFunction displacement(space);
    fillDisplacement(displacement.getData(), 0.02);

    TrialFunction du(space);
    TestFunction v(space);
    auto f = [](const Math::SpatialPoint& x)
    {
      Matrix2 W = Matrix2::Identity();
      W(0, 0) = Real(0.15) + Real(0.1) * x[0];
      W(1, 1) = Real(0.15) + Real(0.1) * x[1];
      return W;
    };
    QualityTerm quality(SquaredDistanceMetric{}, makeAnalyticTargetJacobian(f), 1.0);
    DeviationTerm deviation(1e-2);

    auto assembleResidual = [&]()
    {
      LinearForm residual(v);
      residual =
          quality.residual(displacement, v)
        + deviation.residual(displacement, v);
      residual.assemble();
      return residual.getVector();
    };

    BilinearForm tangent(du, v);
    tangent =
        quality.tangent(displacement, du, v)
      + deviation.tangent(du, v);
    tangent.assemble();

    const auto residual0 = assembleResidual();
    auto direction = displacement.getData();
    fillDirection(direction);
    const auto original = displacement.getData();
    const Real eps = 1e-7;
    displacement.getData() = original + eps * direction;
    const auto residual1 = assembleResidual();
    displacement.getData() = original;

    const auto fd = (residual1 - residual0) / eps;
    const auto jd = tangent.getOperator() * direction;
    const Real denom = std::max<Real>({Real(1), fd.norm(), jd.norm()});
    EXPECT_LT((fd - jd).norm() / denom, 1e-7);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, IsoparametricTargetMatchesAffineTriangleInitialElement)
  {
    for (const auto& tri : {
           makeTriangle({0, 0}, {1, 0}, {0, 1}),
           makeTriangle({0.2, 0.1}, {1.2, 0.3}, {0.4, 1.4}),
           makeTriangle({2, -1}, {2.5, -0.5}, {1.5, 1}) })
    {
      auto mesh = tri;
      IsoparametricTargetJacobian iso(mesh);
      InitialElementTargetJacobian initial(mesh);
      auto it = mesh.getPolytope(2, 0);
      const Math::SpatialPoint rc({Real(1) / 3, Real(1) / 3});
      // Affine triangle: isoparametric Jacobian == P1 element Jacobian.
      EXPECT_NEAR(
          (iso.evaluate(*it, rc) - initial.evaluate(*it, rc)).norm(),
          0, 1e-12);

      // Hence zero strict energy at u = 0 on the captured mesh, exactly like
      // InitialElementTargetJacobian.
      P1 space(mesh, 2);
      GridFunction displacement(space);
      displacement.getData().setZero();
      QualityTerm quality(SquaredDistanceMetric{}, iso);
      EXPECT_NEAR(quality.energy(displacement), 0, 1e-12);
    }
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, IsoparametricTargetWorksOnQuadrilateralPolytope)
  {
    auto mesh = makeUnitQuad();
    IsoparametricTargetJacobian target(mesh);
    auto it = mesh.getPolytope(2, 0);
    const auto& cell = *it;

    const Geometry::Polytope::Traits traits(cell.getGeometry());
    const Math::SpatialPoint rc = traits.getCentroid();
    const Geometry::Point point(cell, rc);
    const auto& J = point.getJacobian();
    Matrix2 expected;
    expected(0, 0) = J(0, 0);
    expected(0, 1) = J(0, 1);
    expected(1, 0) = J(1, 0);
    expected(1, 1) = J(1, 1);

    const Matrix2 got = target.evaluate(cell, rc);
    // The generic captured target runs on a non-triangle polytope and returns
    // that cell's isoparametric Jacobian (FES/polytope independent).
    EXPECT_NEAR((got - expected).norm(), 0, 1e-12);
    EXPECT_GT(std::abs(got.determinant()), 1e-12);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, ShapeDistortionMetricZeroOnRotationsAndScalings)
  {
    ShapeDistortionMetric metric;
    EXPECT_NEAR(metric.value(Matrix2::Identity()), 0, 1e-13);
    EXPECT_NEAR(metric.gradient(Matrix2::Identity()).norm(), 0, 1e-13);

    for (Real theta : { Real(0.3), Real(1.1), Real(-0.7) })
    {
      const Real c = std::cos(theta);
      const Real s = std::sin(theta);
      Matrix2 R;
      R << c, -s, s, c;
      // Scale-invariant: rotations and uniform scalings are exactly minimal.
      EXPECT_NEAR(metric.value(R), 0, 1e-12);
      EXPECT_NEAR(metric.value(Real(3.4) * R), 0, 1e-12);
      EXPECT_NEAR(metric.gradient(R).norm(), 0, 1e-12);
    }

    Matrix2 shear;
    shear << 1, 0.6, 0, 1;
    EXPECT_GT(metric.value(shear), 0);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, ShapeDistortionMetricGradientFiniteDifference)
  {
    ShapeDistortionMetric metric;
    std::vector<Matrix2> cases(4);
    cases[0] << 1.2, 0.3, 0.1, 0.9;
    cases[1] << 0.7, -0.2, 0.25, 1.4;
    cases[2] << 2.0, 0.0, 0.0, 0.5;
    cases[3] << 1.0, 0.8, -0.3, 1.1;

    const Real eps = 1e-7;
    for (const auto& T : cases)
    {
      const Matrix2 g = metric.gradient(T);
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
        {
          Matrix2 Tp = T, Tm = T;
          Tp(i, j) += eps;
          Tm(i, j) -= eps;
          const Real fd =
            (metric.value(Tp) - metric.value(Tm)) / (Real(2) * eps);
          EXPECT_NEAR(g(i, j), fd, 1e-5);
        }
    }
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, ShapeDistortionMetricHessianActionFiniteDifference)
  {
    ShapeDistortionMetric metric;
    std::vector<Matrix2> cases(3);
    cases[0] << 1.2, 0.3, 0.1, 0.9;
    cases[1] << 0.7, -0.2, 0.25, 1.4;
    cases[2] << 1.0, 0.8, -0.3, 1.1;
    Matrix2 H;
    H << 0.13, -0.21, 0.07, 0.31;

    const Real eps = 1e-7;
    for (const auto& T : cases)
    {
      const Matrix2 ha = metric.hessianAction(T, H);
      const Matrix2 fd =
        (metric.gradient(T + eps * H) - metric.gradient(T - eps * H))
        / (Real(2) * eps);
      EXPECT_LT((ha - fd).norm() / std::max(Real(1), fd.norm()), 1e-5);
    }
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, ShapeDistortionMetricIsBarrier)
  {
    ShapeDistortionMetric metric;
    Real previous = -1;
    for (Real s : { Real(0.5), Real(0.1), Real(0.02), Real(0.004) })
    {
      Matrix2 T;
      T << 1, 0, 0, s;  // det = s -> 0+
      const Real v = metric.value(T);
      EXPECT_GT(v, previous);  // strictly increasing as det -> 0+
      previous = v;
    }
    EXPECT_GT(previous, 100);  // unbounded growth near singularity
  }

  // The barrier metric correctly wired into the production assembly: the
  // residual is the energy gradient and the tangent is its Jacobian. This is
  // the prerequisite for using ShapeDistortionMetric in the Newton path; it is
  // deterministic and independent of solver step control.
  TEST(Rodin_Adaptation_TargetMatrixOptimization, ShapeDistortionResidualTangentFiniteDifferenceConsistent)
  {
    auto mesh = makeUniformSquare(6);
    perturbMesh(mesh, 0.008);

    P1 space(mesh, 2);
    GridFunction displacement(space);
    fillDisplacement(displacement.getData(), 0.01);

    TrialFunction du(space);
    TestFunction v(space);
    QualityTerm quality(
        ShapeDistortionMetric{}, EquilateralSameAreaTargetJacobian(mesh), 1.0);
    DeviationTerm deviation(1e-2);

    auto assembleResidual = [&]()
    {
      LinearForm residual(v);
      residual =
          quality.residual(displacement, v)
        + deviation.residual(displacement, v);
      residual.assemble();
      return residual.getVector();
    };

    BilinearForm tangent(du, v);
    tangent =
        quality.tangent(displacement, du, v)
      + deviation.tangent(du, v);
    tangent.assemble();

    const auto residual0 = assembleResidual();
    auto direction = displacement.getData();
    fillDirection(direction);
    const auto original = displacement.getData();
    const Real eps = 1e-7;
    displacement.getData() = original + eps * direction;
    const auto residual1 = assembleResidual();
    displacement.getData() = original;

    const auto fd = (residual1 - residual0) / eps;
    const auto jd = tangent.getOperator() * direction;
    const Real denom = std::max<Real>({Real(1), fd.norm(), jd.norm()});
    EXPECT_LT((fd - jd).norm() / denom, 2e-6);
  }

  // NOTE: a full Newton minimization with ShapeDistortionMetric is
  // intentionally NOT asserted here. The pure shape barrier is rank-deficient
  // (scale/rotation null space) and steers toward det <= 0 away from the
  // solution; minimizing it robustly needs a line search that rejects steps
  // crossing the barrier (or an MFEM-style worst-case/regularized metric).
  // NewtonSolver currently offers only constant damping, no line search, so
  // it diverges on this barrier regardless of damping. The barrier metric
  // itself is finished and verified by the metric-level and production
  // residual/tangent finite-difference tests above; safe step control is the
  // next pipeline item.

  TEST(Rodin_Adaptation_TargetMatrixOptimization, PathologicalNearVertexCutStillAssemblesConsistentTerms)
  {
    auto mesh = makeTwoTriangleSquare();
    P1 space(mesh);
    GridFunction phi(space);
    phi[0] = -1e-12;
    phi[1] = 1;
    phi[2] = 1;
    phi[3] = 1;

    const auto cut = LevelSetDiscretizerTriangles(phi)
      .setSignTolerance(0)
      .setDiagnosticTolerance(1e-8)
      .setInterfaceAttribute(99)
      .setNegativeCellAttribute(1)
      .setPositiveCellAttribute(2)
      .discretize();

    EXPECT_GT(cut.report.pathologicalCutCount, 0);
    expectCutHasUsableProvenance(cut);

    const Real error =
      finiteDifferenceFullTmopError(cut, 1e-5, 1.0, 10.0, 100.0);
    EXPECT_LT(error, 1e-4);
  }

  struct MetricMatrixCase
  {
    Matrix2 A;
    Real expectedSquaredDistance;
  };

  class Rodin_Adaptation_TargetMatrixOptimization_MetricMatrices
    : public testing::TestWithParam<MetricMatrixCase>
  {};

  TEST_P(Rodin_Adaptation_TargetMatrixOptimization_MetricMatrices, SquaredDistanceMatchesFormula)
  {
    const auto c = GetParam();
    SquaredDistanceMetric metric;
    EXPECT_NEAR(metric.value(c.A), c.expectedSquaredDistance, 1e-14);
  }

  INSTANTIATE_TEST_SUITE_P(
      Values,
      Rodin_Adaptation_TargetMatrixOptimization_MetricMatrices,
      testing::Values(
        MetricMatrixCase{(Matrix2() << 1, 0, 0, 1).finished(), 0},
        MetricMatrixCase{(Matrix2() << 2, 0, 0, 1).finished(), 0.5},
        MetricMatrixCase{(Matrix2() << 1, 0.5, 0, 1).finished(), 0.125},
        MetricMatrixCase{(Matrix2() << 1, 0, 0.25, 1).finished(), 0.03125},
        MetricMatrixCase{(Matrix2() << 1.5, 0, 0, 1.5).finished(), 0.25},
        MetricMatrixCase{(Matrix2() << 0.5, 0, 0, 0.5).finished(), 0.25},
        MetricMatrixCase{(Matrix2() << 1, 0.25, 0.25, 1).finished(), 0.0625},
        MetricMatrixCase{(Matrix2() << 0.75, 0.1, -0.2, 1.25).finished(), 0.0875}));

  class Rodin_Adaptation_TargetMatrixOptimization_ShapeAreaMetrics
    : public testing::TestWithParam<Matrix2>
  {};

  TEST_P(Rodin_Adaptation_TargetMatrixOptimization_ShapeAreaMetrics, ShapeSizeIsNonnegativeForValidMaps)
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
      Rodin_Adaptation_TargetMatrixOptimization_ShapeAreaMetrics,
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

  class Rodin_Adaptation_TargetMatrixOptimization_AffineTriangles
    : public testing::TestWithParam<AffineTriangleCase>
  {};

  TEST_P(Rodin_Adaptation_TargetMatrixOptimization_AffineTriangles, P2UpgradePreservesAffineJacobian)
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

  TEST_P(Rodin_Adaptation_TargetMatrixOptimization_AffineTriangles, StrictQualityEnergyMatchesIdentityTargetFormula)
  {
    const auto c = GetParam();
    auto mesh = makeTriangle(c.a, c.b, c.c);
    P1 space(mesh, 2);
    GridFunction displacement(space);
    displacement.getData().setZero();

    SquaredDistanceMetric metric;
    QualityTerm quality(metric, IdentityTargetJacobian{});
    const Real energy = quality.energy(displacement);

    // Identity target means W = I and T = A, integrated on the reference
    // triangle whose area is 1/2.
    EXPECT_NEAR(energy, Real(0.5) * metric.value(c.J), 1e-13);
  }

  TEST_P(Rodin_Adaptation_TargetMatrixOptimization_AffineTriangles, InitialElementTargetGivesZeroEnergyAtZeroDisplacement)
  {
    const auto c = GetParam();
    auto mesh = makeTriangle(c.a, c.b, c.c);
    P1 space(mesh, 2);
    GridFunction displacement(space);
    displacement.getData().setZero();

    QualityTerm quality(
        SquaredDistanceMetric{},
        InitialElementTargetJacobian(mesh));

    EXPECT_NEAR(quality.energy(displacement), 0, 1e-13);
  }

  TEST_P(Rodin_Adaptation_TargetMatrixOptimization_AffineTriangles, InitialElementTargetEnergyIncreasesWithDisplacement)
  {
    const auto c = GetParam();
    auto mesh = makeTriangle(c.a, c.b, c.c);
    P1 space(mesh, 2);
    GridFunction displacement(space);
    fillDisplacement(displacement.getData(), 0.01);

    QualityTerm quality(
        SquaredDistanceMetric{},
        InitialElementTargetJacobian(mesh));

    EXPECT_GT(quality.energy(displacement), 0);
  }

  INSTANTIATE_TEST_SUITE_P(
      AffineMaps,
      Rodin_Adaptation_TargetMatrixOptimization_AffineTriangles,
      testing::Values(
        AffineTriangleCase{{0, 0}, {1, 0}, {0, 1}, (Matrix2() << 1, 0, 0, 1).finished()},
        AffineTriangleCase{{1, 2}, {3, 2}, {1, 5}, (Matrix2() << 2, 0, 0, 3).finished()},
        AffineTriangleCase{{-1, 0}, {0, 0.5}, {-0.5, 2}, (Matrix2() << 1, 0.5, 0.5, 2).finished()},
        AffineTriangleCase{{0.2, 0.1}, {1.2, 0.3}, {0.4, 1.4}, (Matrix2() << 1, 0.2, 0.2, 1.3).finished()},
        AffineTriangleCase{{2, -1}, {2.5, -0.5}, {1.5, 1}, (Matrix2() << 0.5, -0.5, 0.5, 2).finished()},
        AffineTriangleCase{{-2, -1}, {-0.5, -1}, {-1.5, 0.25}, (Matrix2() << 1.5, 0.5, 0, 1.25).finished()}));

  class Rodin_Adaptation_TargetMatrixOptimization_MidpointPerturbations
    : public testing::TestWithParam<std::pair<Index, Math::SpatialPoint>>
  {};

  TEST_P(Rodin_Adaptation_TargetMatrixOptimization_MidpointPerturbations, PerturbingFreeNodeChangesObjective)
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
      Rodin_Adaptation_TargetMatrixOptimization_MidpointPerturbations,
      testing::Values(
        std::pair<Index, Math::SpatialPoint>{3, Math::SpatialPoint({0.00, 0.10})},
        std::pair<Index, Math::SpatialPoint>{3, Math::SpatialPoint({0.10, 0.00})},
        std::pair<Index, Math::SpatialPoint>{4, Math::SpatialPoint({0.00, 0.10})},
        std::pair<Index, Math::SpatialPoint>{4, Math::SpatialPoint({-0.10, 0.00})},
        std::pair<Index, Math::SpatialPoint>{5, Math::SpatialPoint({0.08, 0.00})},
        std::pair<Index, Math::SpatialPoint>{5, Math::SpatialPoint({0.00, -0.08})}));

  class Rodin_Adaptation_TargetMatrixOptimization_DeviationWeights
    : public testing::TestWithParam<Real>
  {};

  TEST_P(Rodin_Adaptation_TargetMatrixOptimization_DeviationWeights, DeviationScalesWithWeight)
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
      Rodin_Adaptation_TargetMatrixOptimization_DeviationWeights,
      testing::Values(0.0, 0.1, 1.0, 10.0, 100.0));

  class Rodin_Adaptation_TargetMatrixOptimization_OptimizerPerturbations
    : public testing::TestWithParam<std::pair<Index, Math::SpatialPoint>>
  {};

  TEST_P(Rodin_Adaptation_TargetMatrixOptimization_OptimizerPerturbations, OptimizerReducesDifferentPerturbations)
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
      Rodin_Adaptation_TargetMatrixOptimization_OptimizerPerturbations,
      testing::Values(
        std::pair<Index, Math::SpatialPoint>{3, Math::SpatialPoint({0.00, 0.20})},
        std::pair<Index, Math::SpatialPoint>{3, Math::SpatialPoint({0.15, 0.00})},
        std::pair<Index, Math::SpatialPoint>{4, Math::SpatialPoint({-0.15, 0.00})},
        std::pair<Index, Math::SpatialPoint>{4, Math::SpatialPoint({0.00, -0.12})},
        std::pair<Index, Math::SpatialPoint>{5, Math::SpatialPoint({0.12, 0.12})}));

  class Rodin_Adaptation_TargetMatrixOptimization_LinearMeshObjective
    : public testing::TestWithParam<AffineTriangleCase>
  {};

  TEST_P(Rodin_Adaptation_TargetMatrixOptimization_LinearMeshObjective, FunctionBaseMetricMatchesManualQuadrature)
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
      Rodin_Adaptation_TargetMatrixOptimization_LinearMeshObjective,
      testing::Values(
        AffineTriangleCase{{0, 0}, {1, 0}, {0, 1}, (Matrix2() << 1, 0, 0, 1).finished()},
        AffineTriangleCase{{0, 0}, {2, 0}, {0, 1}, (Matrix2() << 2, 0, 0, 1).finished()},
        AffineTriangleCase{{0, 0}, {1, 0.2}, {0.1, 1.1}, (Matrix2() << 1, 0.1, 0.2, 1.1).finished()},
        AffineTriangleCase{{1, 1}, {1.5, 1}, {1, 1.5}, (Matrix2() << 0.5, 0, 0, 0.5).finished()}));

  TEST(Rodin_Adaptation_TargetMatrixOptimization, ProblemRequiresInitialization)
  {
    auto mesh = makeUnitTriangle();
    Rodin::Adaptation::TargetMatrixOptimization::Problem problem(mesh);

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

  TEST(Rodin_Adaptation_TargetMatrixOptimization, ProblemInitializesAndEvaluates)
  {
    auto mesh = makeUnitTriangle();
    SquaredDistanceMetric metric;
    Rodin::Adaptation::TargetMatrixOptimization::Problem problem(mesh);
    problem
      .setMetric(metric)
      .setDeviationWeight(0.1)
      .initialize();

    EXPECT_TRUE(problem.isInitialized());
    EXPECT_EQ(problem.invalidElementCount(), 0);
    EXPECT_TRUE(std::isfinite(problem.value()));
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, ProblemOptimizerUpdatesGeometry)
  {
    auto mesh = makeUnitTriangle();
    SquaredDistanceMetric metric;
    Rodin::Adaptation::TargetMatrixOptimization::Problem problem(mesh);
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

  TEST(Rodin_Adaptation_TargetMatrixOptimization, ProblemCanLeaveOriginalVerticesFree)
  {
    auto mesh = makeUnitTriangle();
    Rodin::Adaptation::TargetMatrixOptimization::Problem problem(mesh);
    problem
      .setFixOriginalVertices(false)
      .initialize();

    EXPECT_FALSE(problem.getGeometry().nodes[0].fixed);
    EXPECT_FALSE(problem.getGeometry().nodes[1].fixed);
    EXPECT_FALSE(problem.getGeometry().nodes[2].fixed);
  }
}
