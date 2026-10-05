/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <Rodin/Assembly.h>
#include <Rodin/Solver/CG.h>

#include "examples/ShapeOptimization/KelvinBall/RotatedNitscheIntegrator.h"
#include "examples/ShapeOptimization/KelvinBall/TransportProjection.h"
#include "examples/ShapeOptimization/KelvinBall/HilbertIdentification.h"
#include "examples/ShapeOptimization/KelvinBall/Thickness.h"
#include "examples/ShapeOptimization/KelvinBall/RotatedCharacteristicContinuation.h"

namespace KelvinBall
{
  class CutCoverageTest : public ::testing::Test
  {
    protected:
      struct SamplingFlow
      {
          Real cutoff;
          struct Trace
          {
              const Geometry::Point& point;
              bool failed;
              bool exited() const
              {
                return failed;
              }
              const Geometry::Point& getPoint() const
              {
                return point;
              }
              Real getCorrection() const
              {
                return 0;
              }
          };
          Trace trace(const Geometry::Point& point) const
          {
            return {point, point.getReferenceCoordinates()(0) < cutoff};
          }
      };
      Mesh makeMesh(Real shift)
      {
        Mesh::Builder builder;
        builder.initialize(3).nodes(8);
        for (const auto& x : std::array<std::array<Real, 3>, 8>{
               {{0, 0, 0}, {1, 0, 0}, {0, 1, -1}, {0.2, 0.5, 0}, {shift, 0, 0},
                 {shift + 1, 0, 0}, {shift, 1, 1}, {shift + 0.2, 0.5, 0}}})
        {
          Math::SpatialPoint point(3);
          for (size_t i = 0; i < 3; ++i)
            point(i) = x[i];
          builder.vertex(point);
        }
        Index cell;
        IndexArray vertices(4);
        vertices << 0, 1, 2, 3;
        builder.polytope(Polytope::Type::Tetrahedron, vertices, cell);
        builder.attribute({3, cell}, Fluid);
        vertices << 5, 4, 6, 7;
        builder.polytope(Polytope::Type::Tetrahedron, vertices, cell);
        builder.attribute({3, cell}, Fluid);
        Mesh mesh = builder.finalize();
        mesh.getConnectivity().compute(2, 3);
        for (auto face = mesh.getBoundary(); face; ++face)
        {
          bool slave = true, master = true;
          for (Index v : face->getVertices())
          {
            const auto& x = mesh.getVertexCoordinates(v);
            slave = slave && std::abs(x(1) + x(2)) < 1e-12;
            master = master && std::abs(x(1) - x(2)) < 1e-12;
          }
          mesh.setAttribute({2, face->getIndex()},
            slave      ? SigmaMinus
              : master ? SigmaPlus
                       : Outer);
        }
        prepare(mesh);
        return mesh;
      }

      using System = Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;

      System makeSystem(size_t size)
      {
        System system;
        system.getOperator().resize(size, size);
        system.getVector().setZero(size);
        system.getSolution().setZero(size);
        return system;
      }
  };

  TEST_F(CutCoverageTest, MissingRotatedCharacteristicReturnsFailureWithoutThrowing)
  {
    auto mesh = makeMesh(2);
    AttributeFaceLocator locator(mesh, FlatSet<Attribute>{SigmaPlus, SigmaMinus});
    RotatedCharacteristicContinuation continuation(-0.1, mesh, locator, RotationPairs);
    bool checked = false;
    for (auto cell = mesh.getPolytope(3); cell; ++cell)
    {
      const auto& faces = mesh.getConnectivity().getIncidence({3, 2}, cell->getIndex());
      for (size_t local = 0; local < faces.size(); ++local)
      {
        const auto face = mesh.getPolytope(2, faces[local]);
        if (face->getAttribute() != SigmaMinus)
          continue;
        const auto& formula = QF::PolytopeQuadratureFormula::get(2, face->getGeometry());
        const auto& point = face->getQuadrature(formula).getPoint(0);
        Math::SpatialPoint reference;
        cell->getTransformation().inverse(reference, point.getPhysicalCoordinates());
        Index index = cell->getIndex();
        Real time = -0.1, correction = 0;
        BoundaryHit hit{time, index, reference, local, correction};
        EXPECT_FALSE(continuation(hit));
        checked = true;
      }
    }
    EXPECT_TRUE(checked);
  }

  TEST_F(CutCoverageTest, TransportCompleteCoverageKeepsAllSamples)
  {
    auto mesh = makeMesh(0);
    P1 space(mesh);
    GridFunction distance(space);
    distance = Real(3);
    TransportProjection projection(mesh, distance, SamplingFlow{-1}, 2);
    EXPECT_EQ(projection.getOmittedCount(), 0);
    EXPECT_EQ(projection.getFallbackCellCount(), 0);
    EXPECT_EQ(projection.getOmittedWeightFraction(), 0);
    for (auto cell = mesh.getPolytope(3); cell; ++cell)
    {
      const auto& formula = QF::PolytopeQuadratureFormula::get(2, cell->getGeometry());
      const auto& quadrature = cell->getQuadrature(formula);
      for (size_t q = 0; q < quadrature.getSize(); ++q)
      {
        EXPECT_EQ(projection.mask(quadrature.getPoint(q)), 1);
        EXPECT_NEAR(projection.value(quadrature.getPoint(q)), 3, 1e-14);
      }
    }
  }

  TEST_F(CutCoverageTest, TransportMissingSamplesRetainPreviousDistanceWhenRankDeficient)
  {
    auto mesh = makeMesh(0);
    P1 space(mesh);
    GridFunction distance(space);
    distance = RealFunction(
      [](const Geometry::Point& p) { return Real(2) + p.getPhysicalCoordinates()(0); });
    TransportProjection projection(mesh, distance, SamplingFlow{2}, 2);
    EXPECT_EQ(projection.getOmittedCount(), projection.getAttemptedCount());
    EXPECT_EQ(projection.getFallbackCellCount(), mesh.getCellCount());
    EXPECT_NEAR(projection.getOmittedWeightFraction(), 1, 1e-14);
    for (auto cell = mesh.getPolytope(3); cell; ++cell)
    {
      const auto& formula = QF::PolytopeQuadratureFormula::get(2, cell->getGeometry());
      const auto& quadrature = cell->getQuadrature(formula);
      for (size_t q = 0; q < quadrature.getSize(); ++q)
      {
        const auto& p = quadrature.getPoint(q);
        EXPECT_EQ(projection.mask(p), 1);
        EXPECT_NEAR(projection.value(p), distance.getValue(p), 1e-14);
      }
    }
  }

  TEST_F(CutCoverageTest, TransportPartialOmissionUsesIdenticalMassAndLoadMasks)
  {
    auto mesh = makeMesh(0);
    P1 space(mesh);
    GridFunction distance(space);
    distance = Real(3);
    TransportProjection projection(mesh, distance, SamplingFlow{0.03}, 4);
    EXPECT_GT(projection.getOmittedCount(), 0);
    EXPECT_LT(projection.getOmittedCount(), projection.getAttemptedCount());
    EXPECT_EQ(projection.getFallbackCellCount(), 0);
    TrialFunction u(space);
    TestFunction v(space);
    RealFunction mask([&](const Geometry::Point& p) { return projection.mask(p); });
    RealFunction value([&](const Geometry::Point& p) { return projection.value(p); });
    auto mass = Integral(mask * u, v);
    auto load = Integral(value, v);
    mass.setOrder(4);
    load.setOrder(4);
    Problem problem(u, v);
    problem = mass - load;
    problem.assemble();
    const auto& system = problem.getLinearSystem();
    const auto constant = Math::Vector<Real>::Constant(space.getSize(), 3);
    EXPECT_LT((system.getOperator() * constant - system.getVector()).norm(), 1e-12);
    EXPECT_GT(system.getOperator().norm(), 0);
  }

  TEST_F(CutCoverageTest, CompletelyMissingTraceAddsNeitherMatrixNorLoad)
  {
    auto mesh = makeMesh(2);
    P1 scalar(mesh);
    P1 vector(mesh, 3);
    GridFunction reference(scalar);
    reference = RealFunction(
      [](const Geometry::Point& p) { return p.getPhysicalCoordinates()(0); });
    RotatedNitscheIntegrator coupling(mesh);
    auto scalarSystem = makeSystem(scalar.getSize());
    EXPECT_NO_THROW(coupling.assembleScalarTracePenalty(scalar, scalarSystem, 1));
    EXPECT_EQ(scalarSystem.getOperator().norm(), 0);
    EXPECT_EQ(scalarSystem.getVector().norm(), 0);
    auto vectorSystem = makeSystem(vector.getSize());
    EXPECT_NO_THROW(coupling.assembleVector(vector, vectorSystem, 1, 10));
    EXPECT_EQ(vectorSystem.getOperator().norm(), 0);
    const size_t block = vector.getSize() + scalar.getSize();
    const std::array<size_t, 6> offsets{0, vector.getSize(), block,
      block + vector.getSize(), 2 * block, 2 * block + vector.getSize()};
    auto stokesSystem = makeSystem(3 * block);
    EXPECT_NO_THROW(
      coupling.assembleStokes<1>(vector, scalar, offsets, stokesSystem, 1, 10, 0.05));
    EXPECT_EQ(stokesSystem.getOperator().norm(), 0);
    EXPECT_EQ(stokesSystem.getVector().norm(), 0);
    GridFunction v(vector);
    v = VectorFunction(0.0, 0.0, 0.0);
    EXPECT_NO_THROW(coupling.scalarJump(reference));
    EXPECT_NO_THROW(coupling.vectorJump(v));
    EXPECT_NO_THROW(coupling.familyJump(v, v, v));
  }

  TEST_F(CutCoverageTest, ZeroStepProjectionRepairsExistingPeriodicMismatch)
  {
    auto mesh = makeMesh(0);
    P1 space(mesh);
    GridFunction previous(space);
    previous = RealFunction([](const Geometry::Point& p) { return p.z(); });
    RotatedNitscheIntegrator coupling(mesh);
    const Real initialJump = coupling.scalarJump(previous);
    ASSERT_GT(initialJump, Real(1e-3));
    Real lastJump = initialJump;
    for (const Real penalty : {Real(1), Real(100)})
    {
      TrialFunction projected(space);
      TestFunction test(space);
      Problem projection(projected, test);
      auto mass = Integral(projected, test);
      auto load = Integral(previous, test);
      mass.setOrder(2);
      load.setOrder(2);
      projection = mass + coupling.scalarTrace(projected, test, penalty) - load;
      projection.assemble();
      auto& system = projection.getLinearSystem();
      const Math::Vector<Real> originalLoad = system.getVector();
      // The old increment penalty instead adds the old jump to this load.
      EXPECT_EQ((system.getVector() - originalLoad).norm(), Real(0));
      Solver::CG(projection).solve();
      const auto& result = projected.getSolution();
      const Real jump = coupling.scalarJump(result);
      EXPECT_LT(jump, lastJump);
      EXPECT_LT(
        (system.getOperator() * result.getData() - originalLoad).norm(), Real(1e-10));
      EXPECT_GT((result.getData() - previous.getData()).norm(), Real(1e-3));
      lastJump = jump;
    }
    EXPECT_LT(lastJump, Real(0.01) * initialJump);
  }

  TEST_F(CutCoverageTest, ZeroJumpProjectionPreservesPeriodicConstant)
  {
    auto mesh = makeMesh(0);
    P1 space(mesh);
    GridFunction previous(space);
    previous = Real(3);
    TrialFunction projected(space);
    TestFunction test(space);
    Problem projection(projected, test);
    RotatedNitscheIntegrator coupling(mesh);
    projection = Integral(projected, test) + coupling.scalarTrace(projected, test, 10) -
      Integral(previous, test);
    projection.assemble();
    Solver::CG(projection).solve();
    EXPECT_LT(
      (projected.getSolution().getData() - previous.getData()).norm(), Real(1e-10));
    EXPECT_LT(coupling.scalarJump(projected.getSolution()), Real(1e-10));
  }

  TEST_F(CutCoverageTest, PartialOverlapRetainsSymmetricNonzeroCoupling)
  {
    auto mesh = makeMesh(0.3);
    P1 scalar(mesh);
    RotatedNitscheIntegrator coupling(mesh);
    size_t matched = 0, missed = 0;
    for (const auto& pair : RotationPairs)
    {
      for (auto face = mesh.getBoundary(); face; ++face)
      {
        if (face->getAttribute() != pair.slave)
          continue;
        const auto& formula = QF::PolytopeQuadratureFormula::get(4, face->getGeometry());
        const auto& quadrature = face->getQuadrature(formula);
        for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
        {
          if (coupling.getLocator().locate(
                pair.master, pair.rotation * quadrature.getPoint(qp).vector()))
            ++matched;
          else
            ++missed;
        }
      }
    }
    ASSERT_GT(matched, 0);
    ASSERT_GT(missed, 0);
    auto system = makeSystem(scalar.getSize());
    EXPECT_NO_THROW(coupling.assembleScalarTracePenalty(scalar, system, 1));
    EXPECT_GT(system.getOperator().norm(), 0);
    EXPECT_LT(
      (system.getOperator() - Math::SparseMatrix<Real>(system.getOperator().transpose()))
        .norm(),
      1e-12);
    Math::Vector<Real> ones = Math::Vector<Real>::Ones(scalar.getSize());
    EXPECT_LT((system.getOperator() * ones).norm(), 1e-12);
  }
  TEST_F(CutCoverageTest, NativeTraceEliminatesNonzeroDirichletColumns)
  {
    auto mesh = makeMesh(0);
    bool fixedFace = false;
    for (auto face = mesh.getBoundary(); face; ++face)
      if (face->getAttribute() == Outer)
      {
        mesh.setAttribute({2, face->getIndex()}, fixedFace ? Gamma : Outer);
        fixedFace = true;
      }
    P1 space(mesh);
    TrialFunction u(space);
    TestFunction v(space);
    RotatedNitscheIntegrator coupling(mesh);
    auto trace = coupling.scalarTrace(u, v, 10);
    Problem unconstrained(u, v);
    unconstrained = Integral(u, v) + trace - Integral(RealFunction{1}, v);
    unconstrained.assemble();
    Math::Matrix<Real> expected = unconstrained.getLinearSystem().getOperator();
    Math::Vector<Real> rhs = unconstrained.getLinearSystem().getVector();
    IndexSet fixed;
    for (auto face = mesh.getBoundary(); face; ++face)
      if (face->getAttribute() == Outer)
        for (const auto dof : space.getDOFs(2, face->getIndex()))
          fixed.insert(dof);
    ASSERT_FALSE(fixed.empty());
    for (const auto dof : fixed)
      rhs -= Real(2) * expected.col(dof);
    for (const auto dof : fixed)
    {
      expected.row(dof).setZero();
      expected.col(dof).setZero();
      expected(dof, dof) = 1;
      rhs(dof) = 2;
    }
    Problem constrained(u, v);
    constrained = Integral(u, v) + trace - Integral(RealFunction{1}, v) +
      DirichletBC(u, RealFunction{2}).on(Outer);
    constrained.assemble();
    const auto& system = constrained.getLinearSystem();
    EXPECT_LT((Math::Matrix<Real>(system.getOperator()) - expected).norm(), 1e-12);
    EXPECT_LT((system.getVector() - rhs).norm(), 1e-12);
  }

  TEST_F(CutCoverageTest, NativeVectorTracePreservesMappedOperator)
  {
    auto mesh = makeMesh(0.3);
    P1 space(mesh, 3);
    TrialFunction u(space);
    TestFunction v(space);
    RotatedNitscheIntegrator coupling(mesh);
    auto raw = makeSystem(space.getSize());
    coupling.assembleVector(space, raw, 4, 320);
    const auto trace = coupling.vectorTrace(u, v, 4, 320);
    EXPECT_GT(raw.getOperator().norm(), 0);
    EXPECT_LT((trace.getOperator() - raw.getOperator()).norm(), 1e-12);
    Problem problem(u, v);
    problem = trace;
    problem.assemble();
    EXPECT_LT(
      (problem.getLinearSystem().getOperator() - raw.getOperator()).norm(), 1e-12);
  }

  TEST_F(CutCoverageTest, HilbertIdentificationReusesOperatorAcrossDifferentLoads)
  {
    auto mesh = makeMesh(0);
    bool fixedFace = false;
    for (auto face = mesh.getBoundary(); face; ++face)
      if (face->getAttribute() == Outer)
      {
        mesh.setAttribute({2, face->getIndex()}, fixedFace ? Gamma : Outer);
        fixedFace = true;
      }
    VelocitySpace space(mesh, 3);
    RotatedNitscheIntegrator coupling(mesh);
    HilbertIdentification reused(space, coupling, Real(0.3), 320);
    HilbertIdentification::Field actual(space), expected(space);
    for (const Real density : {Real(-1), Real(2), Real(0)})
    {
      HilbertIdentification fresh(space, coupling, Real(0.3), 320);
      const auto diagnostics = reused.identify(RealFunction{density}, actual);
      fresh.identify(RealFunction{density}, expected);
      EXPECT_LT(diagnostics.residual, 1e-10);
      EXPECT_LT((actual.getData() - expected.getData()).norm(), 1e-10);
      if (density != 0)
        EXPECT_GT(actual.getData().norm(), 1e-5);
      else
        EXPECT_LT(actual.getData().norm(), 1e-12);
    }
    Math::Vector<Real> additional = Math::Vector<Real>::Zero(space.getSize());
    IndexSet fixed;
    for (auto face = mesh.getBoundary(); face; ++face)
      if (face->getAttribute() == Outer)
        for (const auto dof : space.getDOFs(2, face->getIndex()))
          fixed.insert(dof);
    for (Index dof = 0; dof < space.getSize(); ++dof)
      if (!fixed.contains(dof))
        additional(dof) = Real(0.1) * (dof + 1);
    ASSERT_GT(additional.norm(), 0);
    HilbertIdentification fresh(space, coupling, Real(0.3), 320);
    reused.identify(RealFunction{0}, actual, &additional);
    fresh.identify(RealFunction{0}, expected, &additional);
    EXPECT_LT((actual.getData() - expected.getData()).norm(), 1e-10);
    reused.identify(RealFunction{0}, actual);
    EXPECT_LT(actual.getData().norm(), 1e-12);
  }

  TEST_F(CutCoverageTest, ThicknessMeasurementsShareRaysWithoutDependingOnDiagnostics)
  {
    Mesh mesh = Mesh::UniformGrid(Polytope::Type::Tetrahedron, {5, 5, 5});
    mesh.scale(Real(0.5));
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      const auto x = centroid(mesh, *cell);
      const bool solid =
        x(0) > 0.5 && x(0) < 1.5 && x(1) > 0.5 && x(1) < 1.5 && x(2) > 0.5 && x(2) < 1.5;
      mesh.setAttribute({3, cell->getIndex()}, solid ? Obstacle : Fluid);
    }
    mesh.getConnectivity().compute(2, 3);
    for (auto face = mesh.getPolytope(2); face; ++face)
    {
      const auto& cells = mesh.getConnectivity().getIncidence({2, 3}, face->getIndex());
      if (cells.size() == 2 &&
        mesh.getPolytope(3, cells[0])->getAttribute() !=
          mesh.getPolytope(3, cells[1])->getAttribute())
        mesh.setAttribute({2, face->getIndex()}, Gamma);
    }
    VelocitySpace space(mesh, 3);
    GridFunction diagnosticNormal(space);
    diagnosticNormal = VectorFunction{0, 0, 0};
    Math::Vector<Real> load = Math::Vector<Real>::Zero(space.getSize());
    Math::Vector<Real> diagnosticLoad = load;
    Math::Vector<Real> nominalLoad = load;
    ThicknessPenalty guarded(mesh, Real(1.5));
    const auto result = guarded.evaluate(mesh, space, load, Real(1.2));
    const auto diagnostic =
      guarded.evaluate(mesh, space, diagnosticLoad, Real(1.2), &diagnosticNormal);
    ThicknessPenalty nominal(mesh, Real(1.2));
    const auto separate = nominal.evaluate(mesh, space, nominalLoad);
    ASSERT_GT(result.samples, 0);
    ASSERT_GT(result.violating, 0);
    EXPECT_GT(result.penalty, 0);
    EXPECT_GT(result.nominalPenalty, 0);
    EXPECT_LT((load - diagnosticLoad).norm(), 1e-12);
    EXPECT_DOUBLE_EQ(result.penalty, diagnostic.penalty);
    EXPECT_NEAR(result.nominalPenalty, separate.penalty, 1e-12);
    EXPECT_GT(load.norm(), 0);
    EXPECT_THROW(guarded.evaluate(mesh, space, load, Real(2)), std::runtime_error);
  }

}
