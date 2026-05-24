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
  static constexpr Attribute TestNegative = 1;
  static constexpr Attribute TestPositive = 2;
  static constexpr Attribute TestBoundary = 40;
  static constexpr Attribute TestInterface = 99;

  // These test-only helpers keep Eigen-style 2x2 literals working with
  // SpatialMatrix while preserving the suite's existing row-major semantics.
  inline Math::SpatialMatrix<Real> mat2(Real a, Real b, Real c, Real d)
  {
    Math::SpatialMatrix<Real> m(2, 2);
    m(0, 0) = a; m(0, 1) = b; m(1, 0) = c; m(1, 1) = d;
    return m;
  }

  struct Mat2Init
  {
    Math::SpatialMatrix<Real>& t;
    int i;
    explicit Mat2Init(Math::SpatialMatrix<Real>& m) : t(m), i(0) { t.resize(2, 2); }
    Mat2Init& operator,(Real v) { t(i / 2, i % 2) = v; ++i; return *this; }
  };

  inline Mat2Init operator<<(Math::SpatialMatrix<Real>& m, Real v)
  {
    Mat2Init s(m);
    s.t(0, 0) = v; s.i = 1;
    return s;
  }

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

  static LocalMesh makeTwoTriangleSquareWithVerticalInterface()
  {
    auto mesh = LocalMesh::Builder()
      .initialize(2)
      .nodes(6)
      .vertex({0, 0})
      .vertex({0.5, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({0.5, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 3})
      .polytope(Polytope::Type::Triangle, {1, 4, 3})
      .polytope(Polytope::Type::Triangle, {1, 2, 4})
      .polytope(Polytope::Type::Triangle, {2, 5, 4})
      .finalize();
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    mesh.getConnectivity().compute(1, 2);
    mesh.setAttribute({2, 0}, TestNegative);
    mesh.setAttribute({2, 1}, TestNegative);
    mesh.setAttribute({2, 2}, TestPositive);
    mesh.setAttribute({2, 3}, TestPositive);
    const auto& conn = mesh.getConnectivity();
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto& edge = conn.getPolytope(1, e);
      const auto a = mesh.getVertexCoordinates(edge(0));
      const auto b = mesh.getVertexCoordinates(edge(1));
      if (std::abs(a[0] - Real(0.5)) <= Real(1e-12)
       && std::abs(b[0] - Real(0.5)) <= Real(1e-12))
      {
        mesh.setAttribute({1, e}, TestInterface);
      }
      else
      {
        const bool boundary =
          (std::abs(a[0]) <= Real(1e-12) && std::abs(b[0]) <= Real(1e-12))
       || (std::abs(a[0] - Real(1)) <= Real(1e-12)
        && std::abs(b[0] - Real(1)) <= Real(1e-12))
       || (std::abs(a[1]) <= Real(1e-12) && std::abs(b[1]) <= Real(1e-12))
       || (std::abs(a[1] - Real(1)) <= Real(1e-12)
        && std::abs(b[1] - Real(1)) <= Real(1e-12));
        if (boundary)
          mesh.setAttribute({1, e}, TestBoundary);
      }
    }
    return mesh;
  }

  static LocalMesh makeCurvedP2Triangle()
  {
    auto mesh = makeUnitTriangle();
    Variational::RealH1Element<2> fe(Polytope::Type::Triangle);
    PointCloud pm(2, fe.getCount());
    for (size_t i = 0; i < fe.getCount(); ++i)
    {
      const auto& rc = fe.getNode(i);
      Math::SpatialPoint x{
        rc[0],
        rc[1]
      };

      const bool vertex =
        (std::abs(rc[0]) <= Real(1e-12) && std::abs(rc[1]) <= Real(1e-12))
     || (std::abs(rc[0] - Real(1)) <= Real(1e-12) && std::abs(rc[1]) <= Real(1e-12))
     || (std::abs(rc[0]) <= Real(1e-12) && std::abs(rc[1] - Real(1)) <= Real(1e-12));
      if (!vertex && std::abs(rc[1]) <= Real(1e-12))
        x[1] += Real(0.08);
      if (!vertex && std::abs(rc[0]) <= Real(1e-12))
        x[0] += Real(0.04);

      pm(0, i) = x[0];
      pm(1, i) = x[1];
    }
    mesh.setPolytopeTransformation(
        {2, 0},
        new ParametricTransformation<Variational::RealH1Element<2>>(
            pm, fe));
    return mesh;
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

  static LocalMesh makeUniformSquare(size_t resolution)
  {
    auto mesh =
      LocalMesh::UniformGrid(Polytope::Type::Triangle, {resolution, resolution});
    mesh.scale(Real(1) / static_cast<Real>(resolution - 1));
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    return mesh;
  }

  static void perturbMesh(LocalMesh& mesh, Real amplitude)
  {
    for (Index v = 0; v < mesh.getVertexCount(); ++v)
    {
      auto x = mesh.getVertexCoordinates(v);
      const bool boundary =
        x[0] <= Real(1e-12) || x[0] >= Real(1) - Real(1e-12)
     || x[1] <= Real(1e-12) || x[1] >= Real(1) - Real(1e-12);
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

  template <class Metric>
  static Real finiteDifferenceStrictQualityMetricError(Metric metric)
  {
    auto mesh = makeUniformSquare(6);
    perturbMesh(mesh, 0.006);

    P1 space(mesh, 2);
    GridFunction displacement(space);
    fillDisplacement(displacement.getData(), 0.003);

    TrialFunction du(space);
    TestFunction v(space);
    QualityTerm quality(
        metric,
        OrientedEquilateralSameAreaTargetJacobian(mesh),
        0.7);

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

  enum ObjectiveTermMask : int
  {
    ObjectiveQuality  = 1 << 0,
    ObjectiveGamma    = 1 << 1,
    ObjectivePhase    = 1 << 2,
    ObjectiveDeviation = 1 << 3,
    ObjectiveBoundary = 1 << 4,
    ObjectiveAll =
        ObjectiveQuality
      | ObjectiveGamma
      | ObjectivePhase
      | ObjectiveDeviation
      | ObjectiveBoundary
  };

  static const char* objectiveMaskName(int mask)
  {
    switch (mask)
    {
      case ObjectiveQuality: return "quality";
      case ObjectiveGamma: return "gamma";
      case ObjectivePhase: return "phase";
      case ObjectiveDeviation: return "deviation";
      case ObjectiveBoundary: return "boundary";
      case ObjectiveQuality | ObjectiveGamma | ObjectivePhase:
        return "quality+gamma+phase";
      case ObjectiveGamma | ObjectivePhase | ObjectiveBoundary:
        return "gamma+phase+boundary";
      case ObjectiveAll: return "all";
      default: return "custom";
    }
  }

  template <class Metric, class SpaceFactory>
  static Real finiteDifferenceObjectiveError(
      Metric metric,
      SpaceFactory makeSpace,
      int mask,
      Real finiteDifferenceStep)
  {
    auto mesh = makeTwoTriangleSquareWithVerticalInterface();
    auto space = makeSpace(mesh);
    GridFunction displacement(space);
    fillDisplacement(displacement.getData(), Real(0.002));

    TrialFunction du(space);
    TestFunction v(space);

    auto phi = [](const Math::SpatialPoint& x)
    {
      return x[0] - Real(0.5);
    };
    auto gradPhi = [](const Math::SpatialPoint&)
    {
      return Math::SpatialPoint{1, 0};
    };
    auto boundaryValue = [](const Math::SpatialPoint& x)
    {
      return x[1];
    };
    auto boundaryGradient = [](const Math::SpatialPoint&)
    {
      return Math::SpatialPoint{0, 1};
    };

    QualityTerm quality(metric, IdentityTargetJacobian{}, Real(0.8));
    quality.setQuadratureOrder(4);

    AnalyticLevelSetFitTerm gamma(
        phi, gradPhi, Optional<Attribute>(TestInterface), Real(1.1));
    gamma.setNormalization(Real(1));

    VolumetricPhaseConsistencyTerm phase(
        phi, gradPhi, TestNegative, TestPositive, Real(0.9));
    phase
      .setQuadratureOrder(4)
      .setEpsilon(Real(1))
      .setMargin(Real(1))
      .setNormalization(Real(1));

    DeviationTerm deviation(Real(0.7));

    AnalyticLevelSetFitTerm boundary(
        boundaryValue, boundaryGradient,
        Optional<Attribute>(TestBoundary), Real(0.6));
    boundary.setNormalization(Real(1));

    auto assembleResidual = [&]()
    {
      Math::Vector<Real> r =
        Math::Vector<Real>::Zero(displacement.getData().size());
      if (mask & ObjectiveQuality)
      {
        LinearForm form(v);
        form = quality.residual(displacement, v);
        form.assemble();
        r += form.getVector();
      }
      if (mask & ObjectiveGamma)
      {
        LinearForm form(v);
        form = gamma.residual(displacement, v);
        form.assemble();
        r += form.getVector();
      }
      if (mask & ObjectivePhase)
      {
        LinearForm form(v);
        form = phase.residual(displacement, v);
        form.assemble();
        r += form.getVector();
      }
      if (mask & ObjectiveDeviation)
      {
        LinearForm form(v);
        form = deviation.residual(displacement, v);
        form.assemble();
        r += form.getVector();
      }
      if (mask & ObjectiveBoundary)
      {
        LinearForm form(v);
        form = boundary.residual(displacement, v);
        form.assemble();
        r += form.getVector();
      }
      return r;
    };

    auto tangentAction = [&](const Math::Vector<Real>& direction)
    {
      Math::Vector<Real> jd =
        Math::Vector<Real>::Zero(displacement.getData().size());
      if (mask & ObjectiveQuality)
      {
        BilinearForm form(du, v);
        form = quality.tangent(displacement, du, v);
        form.assemble();
        jd += form.getOperator() * direction;
      }
      if (mask & ObjectiveGamma)
      {
        BilinearForm form(du, v);
        form = gamma.tangent(displacement, du, v);
        form.assemble();
        jd += form.getOperator() * direction;
      }
      if (mask & ObjectivePhase)
      {
        BilinearForm form(du, v);
        form = phase.tangent(displacement, du, v);
        form.assemble();
        jd += form.getOperator() * direction;
      }
      if (mask & ObjectiveDeviation)
      {
        BilinearForm form(du, v);
        form = deviation.tangent(du, v);
        form.assemble();
        jd += form.getOperator() * direction;
      }
      if (mask & ObjectiveBoundary)
      {
        BilinearForm form(du, v);
        form = boundary.tangent(displacement, du, v);
        form.assemble();
        jd += form.getOperator() * direction;
      }
      return jd;
    };

    const auto r0 = assembleResidual();
    auto direction = displacement.getData();
    fillDirection(direction);
    const auto u0 = displacement.getData();
    displacement.getData() = u0 + finiteDifferenceStep * direction;
    const auto r1 = assembleResidual();
    displacement.getData() = u0;

    const auto fd = (r1 - r0) / finiteDifferenceStep;
    const auto jd = tangentAction(direction);
    const Real denom = std::max<Real>({Real(1), fd.norm(), jd.norm()});
    return (fd - jd).norm() / denom;
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, SquaredDistanceMetricIsZeroAtIdentity)
  {
    SquaredDistanceMetric metric;
    EXPECT_NEAR(metric.value(Math::SpatialMatrix<Real>::Identity(2, 2)), 0, 1e-14);
    EXPECT_NEAR(metric.gradient(Math::SpatialMatrix<Real>::Identity(2, 2)).norm(), 0, 1e-14);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, DistortedMatrixHasLargerMetric)
  {
    SquaredDistanceMetric metric;
    Math::SpatialMatrix<Real> A = Math::SpatialMatrix<Real>::Identity(2, 2);
    A(0, 1) = 0.5;

    EXPECT_GT(metric.value(A), metric.value(Math::SpatialMatrix<Real>::Identity(2, 2)));
    EXPECT_NEAR((metric.gradient(A) - (A - Math::SpatialMatrix<Real>::Identity(2, 2))).norm(), 0, 1e-14);
    EXPECT_NEAR(
        (metric.hessianAction(A, Math::SpatialMatrix<Real>::Identity(2, 2)) - Math::SpatialMatrix<Real>::Identity(2, 2)).norm(),
        0,
        1e-14);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, SquaredDistanceMetricFollowsFormulaForSingularAndInvertedMatrices)
  {
    SquaredDistanceMetric metric;
    Math::SpatialMatrix<Real> inverted = Math::SpatialMatrix<Real>::Identity(2, 2);
    inverted(0, 0) = -1;

    Math::SpatialMatrix<Real> nearSingular = Math::SpatialMatrix<Real>::Identity(2, 2);
    nearSingular(1, 1) = 1e-14;

    EXPECT_NEAR(metric.value(inverted), 2, 1e-14);
    EXPECT_NEAR(metric.value(nearSingular), 0.5, 1e-13);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, ShapeAndAreaMetricsAreFiniteForIdentity)
  {
    ShapeDistortionMetric shape;
    AreaDistortionMetric area;
    ShapeSizeMetric shapeSize;

    EXPECT_NEAR(shape.value(Math::SpatialMatrix<Real>::Identity(2, 2)), 0, 1e-14);
    EXPECT_NEAR(area.value(Math::SpatialMatrix<Real>::Identity(2, 2)), 0, 1e-14);
    EXPECT_NEAR(shapeSize.value(Math::SpatialMatrix<Real>::Identity(2, 2)), 0, 1e-14);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, IdentityTargetWorks)
  {
    auto mesh = makeUnitTriangle();
    auto cellIterator = mesh.getPolytope(2, 0);
    const auto& cell = *cellIterator;
    const Math::SpatialPoint rc({Real(1) / Real(3), Real(1) / Real(3)});

    IdentityTargetJacobian target;
    EXPECT_NEAR((target.evaluate(cell, rc) - Math::SpatialMatrix<Real>::Identity(2, 2)).norm(), 0, 1e-14);
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
    EXPECT_NEAR(
        residual.getVector().dot(displacement.getData()),
        2 * term.energy(displacement),
        1e-13);
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

  TEST(Rodin_Adaptation_TargetMatrixOptimization, AnalyticLevelSetFitTermIsZeroOnFreshCutInterface)
  {
    auto cut = makeTwoTriangleSquareWithVerticalInterface();
    auto value = [](const Math::SpatialPoint& x)
    {
      return x[0] - Real(0.5);
    };
    auto gradient = [](const Math::SpatialPoint&)
    {
      return Math::SpatialPoint{1, 0};
    };

    AnalyticLevelSetFitTerm fit(
        value, gradient, Optional<Attribute>(TestInterface), 1.0);

    EXPECT_NEAR(fit.energy(cut), 0, 1e-28);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, AnalyticLevelSetFitTermIncreasesWhenInterfaceIsPerturbed)
  {
    auto cut = makeTwoTriangleSquareWithVerticalInterface();
    auto value = [](const Math::SpatialPoint& x)
    {
      return x[0] - Real(0.5);
    };
    auto gradient = [](const Math::SpatialPoint&)
    {
      return Math::SpatialPoint{1, 0};
    };
    AnalyticLevelSetFitTerm fit(
        value, gradient, Optional<Attribute>(TestInterface), 1.0);

    ASSERT_GT(Index(1), 0);
    Index outputEdge = 0;
    const auto& conn = cut.getConnectivity();
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
      if (const auto attr = cut.getAttribute(1, e);
          attr && *attr == TestInterface)
      {
        outputEdge = e;
        break;
      }
    const auto& edge = cut.getConnectivity().getPolytope(1, outputEdge);
    auto x = cut.getVertexCoordinates(edge(0));
    x[0] += 0.05;
    cut.setVertexCoordinates(edge(0), x);

    EXPECT_GT(fit.energy(cut), 1e-8);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, AnalyticLevelSetFitTermResidualAndTangentAreConsistent)
  {
    auto cut = makeTwoTriangleSquareWithVerticalInterface();
    constexpr size_t vdim = 2;
    P1 space(cut, vdim);

    GridFunction displacement(space);
    auto& data = displacement.getData();
    for (Eigen::Index i = 0; i < data.size(); ++i)
      data(i) = 0.01 * std::sin(static_cast<Real>(i + 1));

    TrialFunction du(space);
    TestFunction v(space);
    auto value = [](const Math::SpatialPoint& x)
    {
      return x[0] - Real(0.5);
    };
    auto gradient = [](const Math::SpatialPoint&)
    {
      return Math::SpatialPoint{1, 0};
    };
    AnalyticLevelSetFitTerm fit(
        value, gradient, Optional<Attribute>(TestInterface), 2.25);

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
      direction(i) = std::cos(Real(0.73) * static_cast<Real>(i + 1));
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

  TEST(Rodin_Adaptation_TargetMatrixOptimization, ObjectiveTermMatrixFiniteDifferenceP1)
  {
    const std::array<int, 8> masks = {{
      ObjectiveQuality,
      ObjectiveGamma,
      ObjectivePhase,
      ObjectiveDeviation,
      ObjectiveBoundary,
      ObjectiveQuality | ObjectiveGamma | ObjectivePhase,
      ObjectiveGamma | ObjectivePhase | ObjectiveBoundary,
      ObjectiveAll
    }};
    const std::array<Real, 4> steps = {{1e-5, 1e-6, 1e-7, 1e-8}};
    auto makeP1 = [](LocalMesh& mesh)
    {
      return P1(mesh, 2);
    };

    auto checkMetric = [&](const char* metricName, auto metric)
    {
      for (const int mask : masks)
      {
        Real best = std::numeric_limits<Real>::infinity();
        for (const Real eps : steps)
        {
          const Real error =
            finiteDifferenceObjectiveError(metric, makeP1, mask, eps);
          best = std::min(best, error);
        }
        SCOPED_TRACE(metricName);
        SCOPED_TRACE(objectiveMaskName(mask));
        EXPECT_LT(best, Real(5e-6));
      }
    };

    checkMetric("squared-distance", SquaredDistanceMetric{});
    checkMetric("area-distortion", AreaDistortionMetric{});
    checkMetric("shape-size-blend", ShapeSizeBlendMetric(Real(0.5)));
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, ObjectiveTermMatrixFiniteDifferenceP2)
  {
    const std::array<int, 6> masks = {{
      ObjectiveQuality,
      ObjectiveGamma,
      ObjectivePhase,
      ObjectiveDeviation,
      ObjectiveBoundary,
      ObjectiveAll
    }};
    const std::array<Real, 4> steps = {{1e-5, 1e-6, 1e-7, 1e-8}};
    auto makeP2 = [](LocalMesh& mesh)
    {
      return VectorH1<2, LocalMesh>(
          std::integral_constant<size_t, 2>{}, mesh, 2);
    };

    auto checkMetric = [&](const char* metricName, auto metric)
    {
      for (const int mask : masks)
      {
        Real best = std::numeric_limits<Real>::infinity();
        for (const Real eps : steps)
        {
          const Real error =
            finiteDifferenceObjectiveError(metric, makeP2, mask, eps);
          best = std::min(best, error);
        }
        SCOPED_TRACE(metricName);
        SCOPED_TRACE(objectiveMaskName(mask));
        EXPECT_LT(best, Real(2e-5));
      }
    };

    checkMetric("squared-distance", SquaredDistanceMetric{});
    checkMetric("area-distortion", AreaDistortionMetric{});
    checkMetric("shape-size-blend", ShapeSizeBlendMetric(Real(0.5)));
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
      Math::SpatialMatrix<Real> A0;
      A0 << std::sqrt(Real(2) * area), 0, 0, std::sqrt(Real(2) * area);
      const Math::SpatialMatrix<Real> W =
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

  TEST(Rodin_Adaptation_TargetMatrixOptimization, OrientedEquilateralSameAreaTargetPreservesAreaAndBestRotation)
  {
    Math::SpatialMatrix<Real> A0;
    A0 << 0.18, -0.04,
          0.03,  0.11;
    const Math::SpatialMatrix<Real> fixed =
      EquilateralSameAreaTargetJacobian::equilateralSameArea(A0);
    const Math::SpatialMatrix<Real> oriented =
      OrientedEquilateralSameAreaTargetJacobian::orientedEquilateralSameArea(A0);

    EXPECT_NEAR(
        std::abs(oriented.determinant()),
        std::abs(fixed.determinant()),
        1e-14);
    EXPECT_LE((A0 - oriented).squaredNorm(), (A0 - fixed).squaredNorm());
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, IdealElementTargetImprovesDistortedTriangleShape)
  {
    auto mesh = makeTriangle({0, 0}, {1.4, 0.0}, {0.08, 0.18});
    auto cell = mesh.getPolytope(2, 0);
    const Math::SpatialPoint rc({Real(1) / 3, Real(1) / 3});

    IdealElementTargetJacobian target(mesh);
    const Math::SpatialMatrix<Real> A0 = linearCellJacobian2D(*cell);
    const Math::SpatialMatrix<Real> W = target.evaluate(*cell, rc);
    const Math::SpatialMatrix<Real> oriented =
      OrientedEquilateralSameAreaTargetJacobian::orientedEquilateralSameArea(A0);

    EXPECT_NEAR((W - oriented).norm(), 0, 1e-14);
    EXPECT_NEAR(std::abs(W.determinant()), std::abs(A0.determinant()), 1e-14);

    const auto x0 = mesh.getVertexCoordinates(0);
    const auto x1 = mesh.getVertexCoordinates(1);
    const auto x2 = mesh.getVertexCoordinates(2);
    const Real beforeQuality = triangleQuality(x0, x1, x2);
    const Real targetQuality = triangleQuality(
        Math::SpatialPoint{0, 0},
        Math::SpatialPoint{W(0, 0), W(1, 0)},
        Math::SpatialPoint{W(0, 1), W(1, 1)});

    EXPECT_GT(targetQuality, beforeQuality);
    EXPECT_NEAR(targetQuality, 1, 1e-12);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, IdealElementTargetSupportsSquareSameAreaTarget)
  {
    Math::SpatialMatrix<Real> A0;
    A0 << 1.4, 0.2,
          0.1, 0.6;
    const Math::SpatialMatrix<Real> W =
      IdealElementTargetJacobian::orientedSquareSameArea(A0);
    const Real area = std::abs(A0.determinant());

    EXPECT_NEAR(std::abs(W.determinant()), area, 1e-14);
    EXPECT_NEAR(W.col(0).norm(), W.col(1).norm(), 1e-14);
    EXPECT_NEAR(W.col(0).dot(W.col(1)), 0, 1e-14);
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
        (ConstantTargetJacobian::identity().getMatrix() - Math::SpatialMatrix<Real>::Identity(2, 2))
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
    EXPECT_NEAR((rot.transpose() * rot - Math::SpatialMatrix<Real>::Identity(2, 2)).norm(), 0, 1e-13);

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
    const Math::SpatialMatrix<Real> wTri = target.evaluate(*triIt, rcA);
    const Math::SpatialMatrix<Real> wQuad = target.evaluate(*quadIt, rcB);

    // Same matrix regardless of polytope type or sample point.
    EXPECT_NEAR((wTri - wQuad).norm(), 0, 1e-15);
    EXPECT_NEAR((wTri - target.getMatrix()).norm(), 0, 1e-15);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, ConstantTargetJacobianEnergyMatchesFormulaOnAffineTriangle)
  {
    // Affine triangle with map J = [b-a, c-a].
    auto mesh = makeTriangle({0.2, 0.1}, {1.2, 0.3}, {0.4, 1.4});
    Math::SpatialMatrix<Real> J;
    J << 1.0, 0.2, 0.2, 1.3;

    P1 space(mesh, 2);
    GridFunction displacement(space);
    displacement.getData().setZero();

    SquaredDistanceMetric metric;
    const auto target = ConstantTargetJacobian::diagonal(1.5, 0.8);
    QualityTerm quality(metric, target);
    const Real energy = quality.energy(displacement);

    const Math::SpatialMatrix<Real> W = target.getMatrix();
    const Math::SpatialMatrix<Real> T = J * W.inverse();
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
      Math::SpatialMatrix<Real> W = Math::SpatialMatrix<Real>::Identity(2, 2);
      W(0, 0) = Real(1) + x[0];
      W(1, 1) = Real(2) + x[1];
      return W;
    };
    auto target = makeAnalyticTargetJacobian(f);

    const Math::SpatialPoint rc({Real(0.25), Real(0.5)});
    const Geometry::Point point(cell, rc);
    const auto phys = point.getPhysicalCoordinates();
    const Math::SpatialMatrix<Real> expected = f(phys);
    const Math::SpatialMatrix<Real> got = target.evaluate(cell, rc);
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
      Math::SpatialMatrix<Real> W = Math::SpatialMatrix<Real>::Identity(2, 2);
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

  TEST(Rodin_Adaptation_TargetMatrixOptimization, ParametricTargetMatchesAffineTriangleInitialElement)
  {
    for (const auto& tri : {
           makeTriangle({0, 0}, {1, 0}, {0, 1}),
           makeTriangle({0.2, 0.1}, {1.2, 0.3}, {0.4, 1.4}),
           makeTriangle({2, -1}, {2.5, -0.5}, {1.5, 1}) })
    {
      auto mesh = tri;
      ParametricTargetJacobian parametric(mesh);
      InitialElementTargetJacobian initial(mesh);
      auto it = mesh.getPolytope(2, 0);
      const Math::SpatialPoint rc({Real(1) / 3, Real(1) / 3});
      // Affine triangle: parametric Jacobian == P1 element Jacobian.
      EXPECT_NEAR(
          (parametric.evaluate(*it, rc) - initial.evaluate(*it, rc)).norm(),
          0, 1e-12);

      // Hence zero strict energy at u = 0 on the captured mesh, exactly like
      // InitialElementTargetJacobian.
      P1 space(mesh, 2);
      GridFunction displacement(space);
      displacement.getData().setZero();
      QualityTerm quality(SquaredDistanceMetric{}, parametric);
      EXPECT_NEAR(quality.energy(displacement), 0, 1e-12);
    }
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, ParametricTargetWorksOnQuadrilateralPolytope)
  {
    auto mesh = makeUnitQuad();
    ParametricTargetJacobian target(mesh);
    auto it = mesh.getPolytope(2, 0);
    const auto& cell = *it;

    const Geometry::Polytope::Traits traits(cell.getGeometry());
    const Math::SpatialPoint rc = traits.getCentroid();
    const Geometry::Point point(cell, rc);
    const auto& J = point.getJacobian();
    Math::SpatialMatrix<Real> expected(2, 2);
    expected(0, 0) = J(0, 0);
    expected(0, 1) = J(0, 1);
    expected(1, 0) = J(1, 0);
    expected(1, 1) = J(1, 1);

    const Math::SpatialMatrix<Real> got = target.evaluate(cell, rc);
    // The generic captured target runs on a non-triangle polytope and returns
    // that cell's parametric Jacobian (FES/polytope independent).
    EXPECT_NEAR((got - expected).norm(), 0, 1e-12);
    EXPECT_GT(std::abs(got.determinant()), 1e-12);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, QualityPreservingTargetIsZeroEnergyOnCurvedP2Geometry)
  {
    auto mesh = makeCurvedP2Triangle();
    VectorH1<2, LocalMesh> space(std::integral_constant<size_t, 2>{}, mesh, 2);
    GridFunction displacement(space);
    displacement.getData().setZero();

    QualityTerm quality(
        SquaredDistanceMetric{},
        QualityPreservingTargetJacobian{});
    quality.setQuadratureOrder(4);

    EXPECT_NEAR(quality.energy(displacement), 0, 1e-12);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, CurvedQualityTargetBlendsNaturalAndIdealJacobians)
  {
    auto mesh = makeCurvedP2Triangle();
    auto cellIterator = mesh.getPolytope(2, 0);
    const auto& cell = *cellIterator;
    const Math::SpatialPoint rc{Real(1) / Real(3), Real(1) / Real(3)};

    QualityPreservingTargetJacobian natural;
    IdealElementTargetJacobian ideal(mesh);
    CurvedQualityTargetJacobian target(mesh, 0.25);

    const Math::SpatialMatrix<Real> Wnatural = natural.evaluate(cell, rc);
    const Math::SpatialMatrix<Real> Wideal = ideal.evaluate(cell, rc);
    const Math::SpatialMatrix<Real> W = target.evaluate(cell, rc);

    EXPECT_GT(W.determinant(), 0);
    EXPECT_GT((W - Wnatural).norm(), 1e-12);
    EXPECT_LT((W - Wnatural).norm(), (Wideal - Wnatural).norm());
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, ProjectedInterfaceTargetCarriesCurvedMidsideGoal)
  {
    auto mesh = makeTriangle({0, 0}, {1, 0}, {0, 0.5});
    const auto& conn = mesh.getConnectivity();
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto& edge = conn.getPolytope(1, e);
      if ((edge(0) == 0 && edge(1) == 1)
          || (edge(0) == 1 && edge(1) == 0))
        mesh.setAttribute({1, e}, TestInterface);
    }

    auto project = [](const Math::SpatialPoint& x)
    {
      return Math::SpatialPoint{
        x[0],
        Real(0.1) * std::sin(Real(3.14159265358979323846) * x[0]) };
    };
    ProjectedInterfaceTargetJacobian target(mesh, TestInterface, project);
    ParametricTargetJacobian affine(mesh);

    auto cellIterator = mesh.getPolytope(2, 0);
    const auto& cell = *cellIterator;
    const Math::SpatialPoint rc{Real(1) / Real(3), Real(1) / Real(3)};
    const auto Wcurved = target.evaluate(cell, rc);
    const auto Waffine = affine.evaluate(cell, rc);

    EXPECT_GT(Wcurved.determinant(), 0);
    EXPECT_GT((Wcurved - Waffine).norm(), Real(1e-6));
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, ProjectedQualityTargetBlendsFitTargetTowardIdeal)
  {
    auto mesh = makeTriangle({0, 0}, {1, 0}, {0, 0.5});
    const auto& conn = mesh.getConnectivity();
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto& edge = conn.getPolytope(1, e);
      if ((edge(0) == 0 && edge(1) == 1)
          || (edge(0) == 1 && edge(1) == 0))
        mesh.setAttribute({1, e}, TestInterface);
    }

    auto project = [](const Math::SpatialPoint& x)
    {
      return Math::SpatialPoint{
        x[0],
        Real(0.1) * std::sin(Real(3.14159265358979323846) * x[0]) };
    };
    ProjectedInterfaceTargetJacobian projected(mesh, TestInterface, project);
    IdealElementTargetJacobian ideal(mesh);
    ProjectedQualityTargetJacobian target(
        mesh, TestInterface, project, Real(0.10));

    auto cellIterator = mesh.getPolytope(2, 0);
    const auto& cell = *cellIterator;
    const Math::SpatialPoint rc{Real(1) / Real(3), Real(1) / Real(3)};
    const auto Wprojected = projected.evaluate(cell, rc);
    const auto Wideal = ideal.evaluate(cell, rc);
    const auto W = target.evaluate(cell, rc);

    EXPECT_GT(W.determinant(), 0);
    EXPECT_GT((W - Wprojected).norm(), Real(1e-12));
    EXPECT_LT((W - Wprojected).norm(), (Wideal - Wprojected).norm());
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, ShapeDistortionMetricZeroOnRotationsAndScalings)
  {
    ShapeDistortionMetric metric;
    EXPECT_NEAR(metric.value(Math::SpatialMatrix<Real>::Identity(2, 2)), 0, 1e-13);
    EXPECT_NEAR(metric.gradient(Math::SpatialMatrix<Real>::Identity(2, 2)).norm(), 0, 1e-13);

    for (Real theta : { Real(0.3), Real(1.1), Real(-0.7) })
    {
      const Real c = std::cos(theta);
      const Real s = std::sin(theta);
      Math::SpatialMatrix<Real> R;
      R << c, -s, s, c;
      // Scale-invariant: rotations and uniform scalings are exactly minimal.
      EXPECT_NEAR(metric.value(R), 0, 1e-12);
      EXPECT_NEAR(metric.value(Real(3.4) * R), 0, 1e-12);
      EXPECT_NEAR(metric.gradient(R).norm(), 0, 1e-12);
    }

    Math::SpatialMatrix<Real> shear;
    shear << 1, 0.6, 0, 1;
    EXPECT_GT(metric.value(shear), 0);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, ShapeDistortionMetricGradientFiniteDifference)
  {
    ShapeDistortionMetric metric;
    std::vector<Math::SpatialMatrix<Real>> cases(4);
    cases[0] << 1.2, 0.3, 0.1, 0.9;
    cases[1] << 0.7, -0.2, 0.25, 1.4;
    cases[2] << 2.0, 0.0, 0.0, 0.5;
    cases[3] << 1.0, 0.8, -0.3, 1.1;

    const Real eps = 1e-7;
    for (const auto& T : cases)
    {
      const Math::SpatialMatrix<Real> g = metric.gradient(T);
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
        {
          Math::SpatialMatrix<Real> Tp = T, Tm = T;
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
    std::vector<Math::SpatialMatrix<Real>> cases(3);
    cases[0] << 1.2, 0.3, 0.1, 0.9;
    cases[1] << 0.7, -0.2, 0.25, 1.4;
    cases[2] << 1.0, 0.8, -0.3, 1.1;
    Math::SpatialMatrix<Real> H;
    H << 0.13, -0.21, 0.07, 0.31;

    const Real eps = 1e-7;
    for (const auto& T : cases)
    {
      const Math::SpatialMatrix<Real> ha = metric.hessianAction(T, H);
      const Math::SpatialMatrix<Real> fd =
        (Real(1) / (Real(2) * eps))
        * (metric.gradient(T + eps * H) - metric.gradient(T - eps * H));
      EXPECT_LT((ha - fd).norm() / std::max(Real(1), fd.norm()), 1e-5);
    }
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, AreaDistortionMetricDerivativesFiniteDifference)
  {
    AreaDistortionMetric metric;
    std::vector<Math::SpatialMatrix<Real>> cases(4);
    cases[0] << 1.2, 0.3, 0.1, 0.9;
    cases[1] << 0.7, -0.2, 0.25, 1.4;
    cases[2] << 2.0, 0.0, 0.0, 0.5;
    cases[3] << 1.0, 0.8, -0.3, 1.1;
    Math::SpatialMatrix<Real> H;
    H << 0.13, -0.21, 0.07, 0.31;

    const Real eps = 1e-7;
    for (const auto& T : cases)
    {
      const Math::SpatialMatrix<Real> g = metric.gradient(T);
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
        {
          Math::SpatialMatrix<Real> Tp = T, Tm = T;
          Tp(i, j) += eps;
          Tm(i, j) -= eps;
          const Real fd =
            (metric.value(Tp) - metric.value(Tm)) / (Real(2) * eps);
          EXPECT_NEAR(g(i, j), fd, 1e-6);
        }

      const Math::SpatialMatrix<Real> ha = metric.hessianAction(T, H);
      const Math::SpatialMatrix<Real> fd =
        (Real(1) / (Real(2) * eps))
        * (metric.gradient(T + eps * H) - metric.gradient(T - eps * H));
      EXPECT_LT((ha - fd).norm() / std::max(Real(1), fd.norm()), 1e-6);
    }
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, ShapeSizeMetricDerivativesFiniteDifference)
  {
    ShapeSizeMetric metric;
    std::vector<Math::SpatialMatrix<Real>> cases(4);
    cases[0] << 1.2, 0.3, 0.1, 0.9;
    cases[1] << 0.7, -0.2, 0.25, 1.4;
    cases[2] << 2.0, 0.0, 0.0, 0.5;
    cases[3] << 1.0, 0.8, -0.3, 1.1;
    Math::SpatialMatrix<Real> H;
    H << 0.13, -0.21, 0.07, 0.31;

    const Real eps = 1e-7;
    for (const auto& T : cases)
    {
      const Math::SpatialMatrix<Real> g = metric.gradient(T);
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
        {
          Math::SpatialMatrix<Real> Tp = T, Tm = T;
          Tp(i, j) += eps;
          Tm(i, j) -= eps;
          const Real fd =
            (metric.value(Tp) - metric.value(Tm)) / (Real(2) * eps);
          EXPECT_NEAR(g(i, j), fd, 1e-5);
        }

      const Math::SpatialMatrix<Real> ha = metric.hessianAction(T, H);
      const Math::SpatialMatrix<Real> fd =
        (Real(1) / (Real(2) * eps))
        * (metric.gradient(T + eps * H) - metric.gradient(T - eps * H));
      EXPECT_LT((ha - fd).norm() / std::max(Real(1), fd.norm()), 2e-5);
    }
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, SizeMetric77DerivativesFiniteDifference)
  {
    SizeMetric77 metric;
    std::vector<Math::SpatialMatrix<Real>> cases(4);
    cases[0] << 1.2, 0.3, 0.1, 0.9;
    cases[1] << 0.7, -0.2, 0.25, 1.4;
    cases[2] << 2.0, 0.0, 0.0, 0.5;
    cases[3] << 1.0, 0.8, -0.3, 1.1;
    Math::SpatialMatrix<Real> H;
    H << 0.13, -0.21, 0.07, 0.31;

    const Real eps = 1e-7;
    for (const auto& T : cases)
    {
      const Math::SpatialMatrix<Real> g = metric.gradient(T);
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
        {
          Math::SpatialMatrix<Real> Tp = T, Tm = T;
          Tp(i, j) += eps;
          Tm(i, j) -= eps;
          const Real fd =
            (metric.value(Tp) - metric.value(Tm)) / (Real(2) * eps);
          EXPECT_NEAR(g(i, j), fd, 1e-5);
        }
      const Math::SpatialMatrix<Real> ha = metric.hessianAction(T, H);
      const Math::SpatialMatrix<Real> fd =
        (Real(1) / (Real(2) * eps))
        * (metric.gradient(T + eps * H) - metric.gradient(T - eps * H));
      EXPECT_LT((ha - fd).norm() / std::max(Real(1), fd.norm()), 2e-5);
    }
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, ShapeSizeBlendMetric80DerivativesFiniteDifference)
  {
    ShapeSizeBlendMetric metric(0.5);
    std::vector<Math::SpatialMatrix<Real>> cases(4);
    cases[0] << 1.2, 0.3, 0.1, 0.9;
    cases[1] << 0.7, -0.2, 0.25, 1.4;
    cases[2] << 2.0, 0.0, 0.0, 0.5;
    cases[3] << 1.0, 0.8, -0.3, 1.1;
    Math::SpatialMatrix<Real> H;
    H << 0.13, -0.21, 0.07, 0.31;

    const Real eps = 1e-7;
    for (const auto& T : cases)
    {
      const Math::SpatialMatrix<Real> g = metric.gradient(T);
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
        {
          Math::SpatialMatrix<Real> Tp = T, Tm = T;
          Tp(i, j) += eps;
          Tm(i, j) -= eps;
          const Real fd =
            (metric.value(Tp) - metric.value(Tm)) / (Real(2) * eps);
          EXPECT_NEAR(g(i, j), fd, 1e-5);
        }
      const Math::SpatialMatrix<Real> ha = metric.hessianAction(T, H);
      const Math::SpatialMatrix<Real> fd =
        (Real(1) / (Real(2) * eps))
        * (metric.gradient(T + eps * H) - metric.gradient(T - eps * H));
      EXPECT_LT((ha - fd).norm() / std::max(Real(1), fd.norm()), 2e-5);
    }
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, ShapeDistortionMetricIsBarrier)
  {
    ShapeDistortionMetric metric;
    Real previous = -1;
    for (Real s : { Real(0.5), Real(0.1), Real(0.02), Real(0.004) })
    {
      Math::SpatialMatrix<Real> T;
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

  TEST(Rodin_Adaptation_TargetMatrixOptimization, AreaDistortionResidualTangentFiniteDifferenceConsistent)
  {
    EXPECT_LT(
      finiteDifferenceStrictQualityMetricError(AreaDistortionMetric{}),
      2e-6);
  }

  TEST(Rodin_Adaptation_TargetMatrixOptimization, ShapeSizeResidualTangentFiniteDifferenceConsistent)
  {
    EXPECT_LT(
      finiteDifferenceStrictQualityMetricError(ShapeSizeMetric{}),
      3e-6);
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

  struct MetricMatrixCase
  {
    Math::SpatialMatrix<Real> A;
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
        MetricMatrixCase{mat2(1, 0, 0, 1), 0},
        MetricMatrixCase{mat2(2, 0, 0, 1), 0.5},
        MetricMatrixCase{mat2(1, 0.5, 0, 1), 0.125},
        MetricMatrixCase{mat2(1, 0, 0.25, 1), 0.03125},
        MetricMatrixCase{mat2(1.5, 0, 0, 1.5), 0.25},
        MetricMatrixCase{mat2(0.5, 0, 0, 0.5), 0.25},
        MetricMatrixCase{mat2(1, 0.25, 0.25, 1), 0.0625},
        MetricMatrixCase{mat2(0.75, 0.1, -0.2, 1.25), 0.0875}));

  class Rodin_Adaptation_TargetMatrixOptimization_ShapeAreaMetrics
    : public testing::TestWithParam<Math::SpatialMatrix<Real>>
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
        mat2(1, 0, 0, 1),
        mat2(2, 0, 0, 2),
        mat2(0.5, 0, 0, 0.5),
        mat2(1, 0.25, 0, 1),
        mat2(1.25, 0, 0.1, 0.9),
        mat2(0.8, -0.2, 0.15, 1.3),
        mat2(1.1, 0.4, -0.1, 1.0),
        mat2(1.0, -0.3, 0.2, 1.2)));

  struct AffineTriangleCase
  {
    Math::SpatialPoint a;
    Math::SpatialPoint b;
    Math::SpatialPoint c;
    Math::SpatialMatrix<Real> J;
  };

  class Rodin_Adaptation_TargetMatrixOptimization_AffineTriangles
    : public testing::TestWithParam<AffineTriangleCase>
  {};

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
        AffineTriangleCase{{0, 0}, {1, 0}, {0, 1}, mat2(1, 0, 0, 1)},
        AffineTriangleCase{{1, 2}, {3, 2}, {1, 5}, mat2(2, 0, 0, 3)},
        AffineTriangleCase{{-1, 0}, {0, 0.5}, {-0.5, 2}, mat2(1, 0.5, 0.5, 2)},
        AffineTriangleCase{{0.2, 0.1}, {1.2, 0.3}, {0.4, 1.4}, mat2(1, 0.2, 0.2, 1.3)},
        AffineTriangleCase{{2, -1}, {2.5, -0.5}, {1.5, 1}, mat2(0.5, -0.5, 0.5, 2)},
        AffineTriangleCase{{-2, -1}, {-0.5, -1}, {-1.5, 0.25}, mat2(1.5, 0.5, 0, 1.25)}));

  // Increment 1 for isoparametric P2 TMOP: prove a curved
  // ParametricTransformation<RealH1Element<2>> set via
  // setPolytopeTransformation makes Geometry::Point::getJacobian() carry the
  // curvature (Jacobian varies across the element). This is the foundational
  // unknown the whole P2 plan rests on. Conformity is separately guaranteed
  // by computing edge-node coords from shared edges (affine-mapped reference
  // nodes), so no Fekete-order bookkeeping is needed.
  TEST(Rodin_Adaptation_TargetMatrixOptimization, P2ParametricTransformationCarriesCurvature)
  {
    auto mesh = makeUnitTriangle();              // V0(0,0) V1(1,0) V2(0,1)
    RealH1Element<2> fe(Polytope::Type::Triangle);
    ASSERT_EQ(fe.getCount(), 6u);
    const auto& ref = RealH1Element<2>::getNodes(Polytope::Type::Triangle);
    ASSERT_EQ(ref.size(), 6u);

    // On the unit reference triangle the affine geometry IS the ref coords.
    PointCloud affine(2, 6);
    PointCloud curved(2, 6);
    for (size_t j = 0; j < 6; ++j)
    {
      affine(0, j) = ref[j][0];
      affine(1, j) = ref[j][1];
      curved(0, j) = ref[j][0];
      curved(1, j) = ref[j][1];
      // Bump the edge node on V0-V1 (ref.y==0, 0<ref.x<1) off the chord.
      if (std::abs(ref[j][1]) < 1e-9
          && ref[j][0] > 1e-9 && ref[j][0] < Real(1) - 1e-9)
        curved(1, j) += Real(0.2);
    }

    const Math::SpatialPoint rc1{ Real(0.25), Real(0.25) };
    const Math::SpatialPoint rc2{ Real(0.50), Real(0.20) };

    mesh.setPolytopeTransformation({ size_t(2), Index(0) },
        new ParametricTransformation<RealH1Element<2>>(affine, fe));
    {
      auto it = mesh.getPolytope(2, 0);
      const Math::SpatialMatrix<Real> Ja1 = Geometry::Point(*it, rc1).getJacobian();
      const Math::SpatialMatrix<Real> Ja2 = Geometry::Point(*it, rc2).getJacobian();
      EXPECT_LT((Ja1 - Ja2).norm(), 1e-10);      // affine => constant Jacobian
    }

    mesh.setPolytopeTransformation({ size_t(2), Index(0) },
        new ParametricTransformation<RealH1Element<2>>(curved, fe));
    {
      auto it = mesh.getPolytope(2, 0);
      const Math::SpatialMatrix<Real> Jc1 = Geometry::Point(*it, rc1).getJacobian();
      const Math::SpatialMatrix<Real> Jc2 = Geometry::Point(*it, rc2).getJacobian();
      EXPECT_GT((Jc1 - Jc2).norm(), 1e-3);       // curvature is captured
      EXPECT_GT(Jc1.determinant(), 0.0);         // still a valid map
    }
  }

  // Increment 2: the existing TMOP QualityTerm residual/tangent must assemble
  // and stay finite-difference-consistent on an H1<2> (quadratic) vector
  // space over a P2 parametric mesh. Proves order-2 TMOP assembly works and
  // R = dE/du, J = d2E/du2 hold at P2 (the integrators were only ever run at
  // P1 before). Conformal P2 geometry via affine-mapped reference nodes.
  TEST(Rodin_Adaptation_TargetMatrixOptimization, P2QualityTermAssemblesAndIsFDConsistent)
  {
    auto mesh = makeTwoTriangleSquare();         // 2 tris sharing an edge
    RealH1Element<2> fe(Polytope::Type::Triangle);
    const auto& ref = RealH1Element<2>::getNodes(Polytope::Type::Triangle);
    const auto& conn = mesh.getConnectivity();
    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
    {
      const auto& cell = conn.getPolytope(2, c);
      const auto Va = mesh.getVertexCoordinates(cell(0));
      const auto Vb = mesh.getVertexCoordinates(cell(1));
      const auto Vc = mesh.getVertexCoordinates(cell(2));
      PointCloud pc(2, 6);
      for (size_t j = 0; j < 6; ++j)
      {
        const Real r = ref[j][0], s = ref[j][1];
        // affine map of the reference node through the cell vertices
        // (shared edges -> identical points by linearity => conformal).
        pc(0, j) = (1 - r - s) * Va[0] + r * Vb[0] + s * Vc[0];
        pc(1, j) = (1 - r - s) * Va[1] + r * Vb[1] + s * Vc[1];
      }
      mesh.setPolytopeTransformation({ size_t(2), c },
          new ParametricTransformation<RealH1Element<2>>(pc, fe));
    }

    constexpr size_t vdim = 2;
    H1 space(std::integral_constant<size_t, 2>{}, mesh, vdim);
    GridFunction displacement(space);
    fillDisplacement(displacement.getData(), 0.01);

    TrialFunction du(space);
    TestFunction v(space);
    QualityTerm quality(
        SquaredDistanceMetric{}, InitialElementTargetJacobian(mesh), 1.0);
    DeviationTerm deviation(1.0);

    auto assembleResidual = [&]()
    {
      LinearForm residual(v);
      residual = quality.residual(displacement, v)
               + deviation.residual(displacement, v);
      residual.assemble();
      return residual.getVector();
    };

    BilinearForm tangent(du, v);
    tangent = quality.tangent(displacement, du, v)
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
    const Real denom = std::max<Real>({ Real(1), fd.norm(), jd.norm() });
    EXPECT_LT((fd - jd).norm() / denom, 1e-6);
    EXPECT_GT(residual0.size(), 0);              // P2 space actually built
  }

  // Increment 3-4 of isoparametric P2 TMOP: an actual nonlinear Newton solve
  // in an H1<2> vector space over a P2 parametric mesh must reduce the strict
  // TMOP energy and keep the deformed map valid (det grad x_h > 0 at every
  // sampled quadrature point). Geometry is a conformal affine-P2
  // representation (reference nodes affine-mapped through each cell's vertices
  // => shared edges identical by linearity). One interior vertex is displaced
  // so the strict energy rises above the topology-fixed P2 baseline; Newton
  // on the now-FES-independent quality+deviation terms must recover most of
  // the excess. This is the end-to-end proof that TMOP optimizes in H1<2>,
  // not merely that residual/tangent are FD-consistent there.
  TEST(Rodin_Adaptation_TargetMatrixOptimization, P2StrictNewtonReducesEnergyOnParametricMesh)
  {
    RealH1Element<2> fe(Polytope::Type::Triangle);
    const auto& ref = RealH1Element<2>::getNodes(Polytope::Type::Triangle);

    auto makeGrid = []()
    {
      auto m = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
      m.getConnectivity().compute(2, 1);
      m.getConnectivity().compute(1, 0);
      return m;
    };

    auto displaceCenter = [](LocalMesh& m)
    {
      for (Index v = 0; v < m.getVertexCount(); ++v)
      {
        const auto x = m.getVertexCoordinates(v);
        if (std::abs(x[0] - Real(1)) < 1e-9 && std::abs(x[1] - Real(1)) < 1e-9)
        {
          m.setVertexCoordinates(v, { Real(1.42), Real(1.31) });
          m.getConnectivity().compute(2, 1);
          m.getConnectivity().compute(1, 0);
          return;
        }
      }
    };

    auto attachAffineP2 = [&](LocalMesh& m)
    {
      const auto& conn = m.getConnectivity();
      for (Index c = 0; c < static_cast<Index>(m.getCellCount()); ++c)
      {
        const auto& cell = conn.getPolytope(2, c);
        const auto Va = m.getVertexCoordinates(cell(0));
        const auto Vb = m.getVertexCoordinates(cell(1));
        const auto Vc = m.getVertexCoordinates(cell(2));
        PointCloud pc(2, 6);
        for (size_t j = 0; j < 6; ++j)
        {
          const Real r = ref[j][0], s = ref[j][1];
          pc(0, j) = (1 - r - s) * Va[0] + r * Vb[0] + s * Vc[0];
          pc(1, j) = (1 - r - s) * Va[1] + r * Vb[1] + s * Vc[1];
        }
        m.setPolytopeTransformation({ size_t(2), c },
            new ParametricTransformation<RealH1Element<2>>(pc, fe));
      }
    };

    auto p2StrictEnergy = [](LocalMesh& m)
    {
      H1 s(std::integral_constant<size_t, 2>{}, m, 2);
      GridFunction z(s);
      z.getData().setZero();
      QualityTerm q(SquaredDistanceMetric{}, IdentityTargetJacobian{}, 1.0);
      q.setQuadratureOrder(4);
      return q.energy(z);
    };

    auto baseline = makeGrid();
    attachAffineP2(baseline);
    const Real baselineEnergy = p2StrictEnergy(baseline);

    auto mesh = makeGrid();
    displaceCenter(mesh);
    attachAffineP2(mesh);
    const Real beforeEnergy = p2StrictEnergy(mesh);
    ASSERT_GT(beforeEnergy, baselineEnergy * (Real(1) + Real(1e-6)));

    H1 space(std::integral_constant<size_t, 2>{}, mesh, 2);
    GridFunction displacement(space);
    displacement.getData().setZero();

    TrialFunction du(space);
    TestFunction v(space);
    QualityTerm quality(SquaredDistanceMetric{}, IdentityTargetJacobian{}, 1.0);
    quality.setQuadratureOrder(4);
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

    const Real afterEnergy = quality.energy(displacement);

    // Validity: min det(grad x_h) over a quadrature stays positive at the
    // solution (the deformed P2 map is not inverted anywhere sampled).
    Real minDet = std::numeric_limits<Real>::infinity();
    const auto& qf =
      QF::PolytopeQuadratureFormula::get(4, Polytope::Type::Triangle);
    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
    {
      auto it = mesh.getPolytope(2, c);
      for (size_t q = 0; q < qf.getSize(); ++q)
      {
        const auto A = deformedCoordinateJacobian(
            displacement, *it, qf.getPoint(q));
        minDet = std::min(minDet, A.determinant());
      }
    }

    EXPECT_GT(report.initial_residual, 0);
    EXPECT_LE(report.final_residual, report.initial_residual);
    EXPECT_GT(minDet, 0);
    EXPECT_LT(afterEnergy, beforeEnergy);
    EXPECT_LT(
        afterEnergy - baselineEnergy,
        Real(0.5) * (beforeEnergy - baselineEnergy));
  }
}
