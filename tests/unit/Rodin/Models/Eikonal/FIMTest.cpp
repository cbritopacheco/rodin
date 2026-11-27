/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Eikonal/FIM.h"
#include "Rodin/Geometry/ParametricTransformation.h"
#include "Rodin/Variational/H1/H1Element.h"
#include "Rodin/Variational/MatrixFunction.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  /// @brief Native cell types with nonzero reference dimension.
  constexpr auto Geometries = std::array{Polytope::Type::Segment,
    Polytope::Type::Triangle, Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
    Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge};

  /**
   * @brief Builds one native reference cell with optional sheared vertex geometry.
   * @param geometry Native cell type.
   * @param shear Whether to shear the top vertices into the base interior.
   * @returns A mesh whose cell vertex indices match the native local indices.
   */
  LocalMesh referenceCell(Polytope::Type geometry, bool shear = false)
  {
    const Polytope::Traits traits(geometry);
    LocalMesh::Builder builder;
    builder.initialize(traits.getDimension()).nodes(traits.getVertexCount());
    IndexArray vertices(traits.getVertexCount());
    for (size_t vertex = 0; vertex < traits.getVertexCount(); ++vertex)
    {
      auto coordinates = traits.getVertex(vertex);
      if (shear && traits.getDimension() == 3)
      {
        coordinates(0) += Real(0.5) * coordinates(2);
        coordinates(1) += Real(0.5) * coordinates(2);
      }
      builder.vertex(coordinates);
      vertices[vertex] = vertex;
    }
    builder.polytope(geometry, vertices);
    auto mesh = builder.finalize();
    mesh.getConnectivity().compute(traits.getDimension(), 0);
    mesh.getConnectivity().compute(0, traits.getDimension());
    return mesh;
  }

  /**
   * @brief Selects native vertices on a reference coordinate plane.
   * @param geometry Native cell type.
   * @param coordinate Reference coordinate fixed to zero.
   * @returns Indices of the source vertices.
   */
  IndexVector planeSeeds(Polytope::Type geometry, size_t coordinate = 0)
  {
    const Polytope::Traits traits(geometry);
    IndexVector seeds;
    for (Index vertex = 0; vertex < traits.getVertexCount(); ++vertex)
    {
      if (traits.getVertex(vertex)(coordinate) == 0)
        seeds.push_back(vertex);
    }
    return seeds;
  }
}

/// @brief Every native family reproduces a coordinate-plane arrival at constant speed.
TEST(FIM, NativeCellsAndScalarTensorReduction)
{
  for (auto geometry : Geometries)
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    auto mesh = referenceCell(geometry);
    const size_t dimension = mesh.getDimension();
    P1 space(mesh);
    GridFunction scalar(space), tensor(space);
    const auto speed = [](const Point&) { return Real(2); };
    MatrixFunction coefficient(
      4 * Math::SpatialMatrix<Real>::Identity(dimension, dimension));
    Eikonal::FIM scalarSolver(scalar, speed);
    Eikonal::FIM tensorSolver(tensor, coefficient);
    const auto seeds = planeSeeds(geometry);
    scalarSolver.seed(seeds).solve();
    tensorSolver.seed(seeds).solve();
    const Polytope::Traits traits(geometry);
    for (Index vertex = 0; vertex < traits.getVertexCount(); ++vertex)
    {
      const Index dof = space.getDOFs(0, vertex)[0];
      EXPECT_NEAR(scalar[dof], traits.getVertex(vertex)(0) / 2, 1e-9);
      EXPECT_NEAR(tensor[dof], scalar[dof], 1e-12);
    }
  }
}

/// @brief Native quadrilateral face interiors supply the minimizing path without diagonals.
TEST(FIM, NativeQuadrilateralFaceInterior)
{
  auto mesh = referenceCell(Polytope::Type::Hexahedron, true);
  P1 space(mesh);
  GridFunction arrival(space);
  Eikonal::FIM solver(arrival, [](const Point&) { return Real(1); });
  solver.seed(planeSeeds(Polytope::Type::Hexahedron, 2)).solve();
  size_t checked = 0;
  for (auto vertex = mesh.getVertex(); vertex; ++vertex)
  {
    const auto& coordinates = vertex->getCoordinates();
    if (coordinates(0) == Real(0.5) && coordinates(1) == Real(0.5) && coordinates(2) == 1)
    {
      EXPECT_NEAR(arrival[space.getDOFs(0, vertex->getIndex())[0]], 1, 1e-9);
      ++checked;
    }
  }
  EXPECT_EQ(checked, 1);
}

/// @brief Native triangle and quadrilateral bases both reach an apex through their interior.
TEST(FIM, NativePyramidAndTetrahedronFaceInteriors)
{
  for (auto geometry : {Polytope::Type::Pyramid, Polytope::Type::Tetrahedron})
  {
    auto mesh = referenceCell(geometry);
    const Polytope::Traits traits(geometry);
    Index apex = 0;
    for (Index vertex = 0; vertex < traits.getVertexCount(); ++vertex)
    {
      if (traits.getVertex(vertex)(2) == 1)
      {
        apex = vertex;
        const Real center = geometry == Polytope::Type::Pyramid ? Real(0.5) : Real(1) / 3;
        mesh.setVertexCoordinates(vertex, Math::SpatialPoint{center, center, 1});
      }
    }
    P1 space(mesh);
    GridFunction arrival(space);
    Eikonal::FIM solver(arrival, [](const Point&) { return Real(1); });
    solver.setLocalTolerance(1e-7).seed(planeSeeds(geometry, 2)).solve();
    EXPECT_NEAR(arrival[space.getDOFs(0, apex)[0]], 1, 1e-10);
  }
}

/// @brief Invalid public inputs fail without relying on assertions.
TEST(FIM, InputContracts)
{
  auto mesh = referenceCell(Polytope::Type::Segment);
  P1 space(mesh);
  GridFunction arrival(space);
  Eikonal::FIM solver(arrival, [](const Point&) { return Real(1); });
  EXPECT_THROW(solver.solve(), Alert::Exception);
  EXPECT_THROW(solver.seed(IndexVector{}), Alert::Exception);
  EXPECT_THROW(solver.seed(IndexVector{2}), Alert::Exception);
  EXPECT_THROW(solver.setTolerance(0), Alert::Exception);
  EXPECT_THROW(
    solver.setTolerance(std::numeric_limits<Real>::infinity()), Alert::Exception);
  EXPECT_THROW(solver.setLocalTolerance(1), Alert::Exception);
  EXPECT_THROW(solver.setMaxIterations(0), Alert::Exception);
  for (Real speed : {Real(0), Real(-1), std::numeric_limits<Real>::quiet_NaN()})
  {
    Eikonal::FIM invalid(arrival, [speed](const Point&) { return speed; });
    invalid.seed(IndexVector{0});
    EXPECT_THROW(invalid.solve(), Alert::Exception);
  }
  for (const int problem : {0, 1, 2, 3})
  {
    auto surface = referenceCell(Polytope::Type::Triangle);
    P1 surfaceSpace(surface);
    GridFunction surfaceArrival(surfaceSpace);
    Math::SpatialMatrix<Real> tensor = Math::SpatialMatrix<Real>::Identity(2, 2);
    if (problem == 0)
      tensor(0, 0) = -1;
    else if (problem == 1)
      tensor(0, 1) = 1;
    else if (problem == 2)
      tensor(0, 0) = std::numeric_limits<Real>::infinity();
    else
      tensor = Math::SpatialMatrix<Real>::Identity(3, 3);
    Eikonal::FIM invalid(surfaceArrival, MatrixFunction(tensor));
    invalid.seed(IndexVector{0});
    EXPECT_THROW(invalid.solve(), Alert::Exception);
  }
}

/// @brief Seeds are order-independent and can be reset; disconnected components stay unreachable.
TEST(FIM, SeedLifecycleAndUnreachableComponents)
{
  auto mesh = LocalMesh::Builder()
                .initialize(1)
                .nodes(4)
                .vertex({0})
                .vertex({1})
                .vertex({2})
                .vertex({3})
                .polytope(Polytope::Type::Segment, {0, 1})
                .polytope(Polytope::Type::Segment, {2, 3})
                .finalize();
  mesh.getConnectivity().compute(1, 0);
  mesh.getConnectivity().compute(0, 1);
  P1 space(mesh);
  GridFunction arrival(space);
  Eikonal::FIM solver(arrival, [](const Point&) { return Real(1); });
  solver.seed(IndexVector{0, 0}).solve();
  EXPECT_NEAR(arrival[space.getDOFs(0, 1)[0]], 1, 1e-12);
  EXPECT_TRUE(std::isinf(arrival[space.getDOFs(0, 2)[0]]));
  solver.seed(IndexVector{1, 3}).solve();
  EXPECT_NEAR(arrival[space.getDOFs(0, 0)[0]], 1, 1e-12);
  EXPECT_NEAR(arrival[space.getDOFs(0, 2)[0]], 1, 1e-12);
  const auto previous = arrival.getData();
  solver.seed(IndexVector{3, 1, 3}).solve();
  EXPECT_LT((arrival.getData() - previous).norm(), 1e-12);
  solver.solve();
  EXPECT_EQ(solver.getIterations(), 0);
  EXPECT_LT((arrival.getData() - previous).norm(), 1e-12);
}

/// @brief A zero-dimensional mesh has prescribed times and no propagation paths.
TEST(FIM, PointGeometry)
{
  auto mesh = LocalMesh::Builder()
                .initialize(1)
                .nodes(1)
                .vertex({0})
                .polytope(Polytope::Type::Point, {0})
                .finalize();
  mesh.getConnectivity().compute(0, 0);
  P1 space(mesh);
  GridFunction arrival(space);
  Eikonal::FIM solver(arrival, [](const Point&) { return Real(1); });
  solver.seed(IndexVector{0}).solve();
  EXPECT_EQ(arrival[space.getDOFs(0, 0)[0]], 0);
  EXPECT_EQ(solver.getIterations(), 0);
}

/// @brief Native mixed cells evaluate discontinuous coefficients on the current cell trace.
TEST(FIM, MixedCellsAndCoefficientTraces)
{
  auto mesh = LocalMesh::Builder()
                .initialize(2)
                .nodes(6)
                .vertex({0, 0})
                .vertex({1, 0})
                .vertex({1, 1})
                .vertex({0, 1})
                .vertex({2, 0})
                .vertex({2, 1})
                .polytope(Polytope::Type::Quadrilateral, {0, 1, 2, 3})
                .polytope(Polytope::Type::Triangle, {1, 4, 2})
                .polytope(Polytope::Type::Triangle, {4, 5, 2})
                .finalize();
  mesh.getConnectivity().compute(2, 0);
  mesh.getConnectivity().compute(0, 2);
  P1 space(mesh);
  GridFunction arrival(space);
  Eikonal::FIM solver(arrival, [](const Point& point) {
    EXPECT_EQ(point.getPolytope().getDimension(), 2);
    return point.getPolytope().getIndex() == 0 ? Real(1) : Real(2);
  });
  solver.seed(IndexVector{0, 3}).solve();
  for (Index vertex : {Index(1), Index(2)})
    EXPECT_NEAR(arrival[space.getDOFs(0, vertex)[0]], 1, 1e-9);
  for (Index vertex : {Index(4), Index(5)})
    EXPECT_NEAR(arrival[space.getDOFs(0, vertex)[0]], 1.5, 1e-9);
}

/// @brief Active-list exhaustion is explicit and retains the latest iterate.
TEST(FIM, IterationLimit)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {4});
  mesh.getConnectivity().compute(1, 0);
  mesh.getConnectivity().compute(0, 1);
  P1 space(mesh);
  GridFunction arrival(space);
  Eikonal::FIM solver(arrival, [](const Point&) { return Real(1); });
  solver.setMaxIterations(1).seed(IndexVector{0});
  EXPECT_THROW(solver.solve(), Alert::Exception);
  EXPECT_EQ(solver.getIterations(), 1);
}

#ifdef RODIN_USE_OPENMP
/// @brief Dependency-level OpenMP updates preserve every serial iterate on native cells.
TEST(FIM, OpenMPMatchesSerialNativeCells)
{
  for (auto geometry : Geometries)
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    auto mesh = dimension == 1 ? LocalMesh::UniformGrid(geometry, {4})
      : dimension == 2         ? LocalMesh::UniformGrid(geometry, {4, 4})
                               : LocalMesh::UniformGrid(geometry, {3, 3, 3});
    mesh.getConnectivity().compute(dimension, 0);
    mesh.getConnectivity().compute(0, dimension);
    P1 space(mesh);
    GridFunction arrival(space);
    Eikonal::FIM solver(arrival, [](const Point& point) {
      return Real(1) + Real(0.1) * point.getPhysicalCoordinates()(0);
    });
    EXPECT_EQ(solver.getThreadCount(), 1);
    solver.seed(IndexVector{0}).solve();
    const auto reference = arrival.getData();
    const size_t iterations = solver.getIterations();
    for (size_t threads : {2u, 4u})
    {
      solver.setThreadCount(threads).seed(IndexVector{0}).solve();
      EXPECT_EQ(solver.getIterations(), iterations);
      EXPECT_TRUE((arrival.getData().array() == reference.array()).all());
    }
    Math::SpatialMatrix<Real> tensor =
      Math::SpatialMatrix<Real>::Identity(dimension, dimension);
    tensor(0, 0) = 4;
    if (dimension > 1)
      tensor(0, 1) = tensor(1, 0) = Real(0.2);
    Eikonal::FIM tensorSolver(arrival, MatrixFunction(tensor));
    tensorSolver.seed(IndexVector{0}).solve();
    const auto tensorReference = arrival.getData();
    const size_t tensorIterations = tensorSolver.getIterations();
    for (size_t threads : {2u, 4u})
    {
      tensorSolver.setThreadCount(threads).seed(IndexVector{0}).solve();
      EXPECT_EQ(tensorSolver.getIterations(), tensorIterations);
      EXPECT_TRUE((arrival.getData().array() == tensorReference.array()).all());
    }
  }
}

/// @brief Worker exceptions are propagated after joining, without terminating OpenMP workers.
TEST(FIM, OpenMPExceptions)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, {4, 4});
  mesh.getConnectivity().compute(2, 0);
  mesh.getConnectivity().compute(0, 2);
  P1 space(mesh);
  GridFunction arrival(space);
  // Collapse the source edge to make a worker's local travel cost invalid.
  mesh.setVertexCoordinates(1, mesh.getVertexCoordinates(0));
  Eikonal::FIM solver(arrival, [](const Point&) { return Real(1); });
  solver.setThreadCount(4).seed(IndexVector{0});
  EXPECT_THROW(solver.solve(), Alert::Exception);
  EXPECT_THROW(solver.setThreadCount(0), Alert::Exception);
}

#else
/// @brief Builds without OpenMP reject an unavailable execution backend explicitly.
TEST(FIM, UnavailableOpenMP)
{
  auto mesh = referenceCell(Polytope::Type::Segment);
  P1 space(mesh);
  GridFunction arrival(space);
  Eikonal::FIM solver(arrival, [](const Point&) { return Real(1); });
  EXPECT_THROW(solver.setThreadCount(2), Alert::Exception);
  solver.setThreadCount(1).seed(IndexVector{0}).solve();
  EXPECT_NEAR(arrival[space.getDOFs(0, 1)[0]], 1, 1e-12);
}
#endif

/// @brief Higher-order charts are rejected for every nonpoint native reference family.
TEST(FIM, RejectsHigherOrderGeometry)
{
  for (auto geometry : Geometries)
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    auto mesh = referenceCell(geometry);
    const size_t dimension = mesh.getDimension();
    RealH1Element<2> element(geometry);
    PointCloud controls(dimension, element.getCount());
    for (size_t node = 0; node < element.getCount(); ++node)
    {
      const auto& reference = element.getNode(node);
      for (size_t coordinate = 0; coordinate < dimension; ++coordinate)
        controls(coordinate, node) = reference(coordinate);
    }
    mesh.setPolytopeTransformation(
      {dimension, 0}, new ParametricTransformation(controls, element));
    P1 space(mesh);
    GridFunction arrival(space);
    Eikonal::FIM solver(arrival, [](const Point&) { return Real(1); });
    solver.seed(IndexVector{0});
    EXPECT_THROW(solver.solve(), Alert::Exception);
  }
}

/// @brief Flat embedded native cells use physical straight distances and ambient metrics.
TEST(FIM, EmbeddedFirstOrderCells)
{
  for (auto geometry : {Polytope::Type::Triangle, Polytope::Type::Quadrilateral})
  {
    const Polytope::Traits traits(geometry);
    LocalMesh::Builder builder;
    builder.initialize(3).nodes(traits.getVertexCount());
    IndexArray vertices(traits.getVertexCount());
    for (Index vertex = 0; vertex < traits.getVertexCount(); ++vertex)
    {
      const auto& r = traits.getVertex(vertex);
      builder.vertex({r(0), r(1), 0});
      vertices[vertex] = vertex;
    }
    builder.polytope(geometry, vertices);
    auto mesh = builder.finalize();
    mesh.getConnectivity().compute(2, 0);
    mesh.getConnectivity().compute(0, 2);
    P1 space(mesh);
    GridFunction arrival(space);
    auto tensor = Math::SpatialMatrix<Real>::Identity(3, 3);
    tensor(0, 0) = 4;
    tensor(2, 2) = 9;
    Eikonal::FIM solver(arrival, MatrixFunction(tensor));
    solver.seed(planeSeeds(geometry)).solve();
    for (Index vertex = 0; vertex < traits.getVertexCount(); ++vertex)
    {
      EXPECT_NEAR(
        arrival[space.getDOFs(0, vertex)[0]], traits.getVertex(vertex)(0) / 2, 1e-10);
    }
  }
}

/// @brief Each solve samples scalar and tensor inputs once per cell, outside worker regions.
TEST(FIM, CoefficientSamplingIsCellwise)
{
  for (auto geometry : Geometries)
  {
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    auto mesh = dimension == 1 ? LocalMesh::UniformGrid(geometry, {4})
      : dimension == 2         ? LocalMesh::UniformGrid(geometry, {3, 3})
                               : LocalMesh::UniformGrid(geometry, {3, 3, 3});
    mesh.getConnectivity().compute(dimension, 0);
    mesh.getConnectivity().compute(0, dimension);
    P1 space(mesh);
    GridFunction arrival(space);
    std::vector<size_t> visits(mesh.getCellCount(), 0);
    bool wrongPoint = false, workerCallback = false;
    Real speed = 1;
    auto sample = [&](const Point& point) {
      const auto& cell = point.getPolytope();
      ++visits.at(cell.getIndex());
      wrongPoint = wrongPoint || cell.getDimension() != dimension ||
        (point.getReferenceCoordinates() -
          Polytope::Traits(cell.getGeometry()).getCentroid())
            .squaredNorm() != 0;
#ifdef RODIN_USE_OPENMP
      workerCallback = workerCallback || omp_in_parallel();
#endif
      return speed;
    };
    Eikonal::FIM solver(arrival, sample);
#ifdef RODIN_USE_OPENMP
    solver.setThreadCount(4);
#endif
    solver.seed(IndexVector{0}).solve();
    EXPECT_FALSE(wrongPoint);
    EXPECT_FALSE(workerCallback);
    for (size_t count : visits)
      EXPECT_EQ(count, 1);
    const auto first = arrival.getData();
    std::fill(visits.begin(), visits.end(), 0);
    speed = 2;
    solver.seed(IndexVector{0}).solve();
    for (size_t count : visits)
      EXPECT_EQ(count, 1);
    EXPECT_LT((arrival.getData() - first / 2).norm(), 1e-8);
    std::fill(visits.begin(), visits.end(), 0);
    MatrixFunction tensor(dimension, dimension, [&](const Point& point) {
      sample(point);
      return Math::SpatialMatrix<Real>::Identity(dimension, dimension);
    });
    Eikonal::FIM tensorSolver(arrival, tensor);
#ifdef RODIN_USE_OPENMP
    tensorSolver.setThreadCount(4);
#endif
    tensorSolver.seed(IndexVector{0}).solve();
    for (size_t count : visits)
      EXPECT_EQ(count, 1);
    EXPECT_FALSE(wrongPoint);
    EXPECT_FALSE(workerCallback);
    EXPECT_LT((arrival.getData() - first).norm(), 1e-8);
  }
}

/// @brief A varying speed uses the explicit midpoint cell value, not a path integral.
TEST(FIM, SpatialCoefficientDiscretization)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {5});
  mesh.scale(Real(0.25));
  mesh.getConnectivity().compute(1, 0);
  mesh.getConnectivity().compute(0, 1);
  P1 space(mesh);
  GridFunction arrival(space);
  Eikonal::FIM solver(arrival, [](const Point& point) { return 1 + point.x(); });
  solver.seed(IndexVector{0}).solve();
  Real expected = 0;
  for (Index vertex = 1; vertex < mesh.getVertexCount(); ++vertex)
  {
    const Real midpoint = (Real(vertex) - Real(0.5)) / 4;
    expected += Real(0.25) / (1 + midpoint);
    EXPECT_NEAR(arrival[space.getDOFs(0, vertex)[0]], expected, 1e-12);
  }
}

/// @brief Analytic edge updates can precede a supporting time on a noncausal stencil.
TEST(FIM, AnalyticRefractionWithoutMarchingCutoff)
{
  auto mesh = LocalMesh::Builder()
                .initialize(2)
                .nodes(4)
                .vertex({0, 0})
                .vertex({-1, Real(0.2)})
                .vertex({2, Real(0.2)})
                .vertex({Real(0.5), Real(1.2)})
                .polytope(Polytope::Type::Triangle, {0, 2, 1})
                .polytope(Polytope::Type::Triangle, {1, 2, 3})
                .finalize();
  mesh.getConnectivity().compute(2, 0);
  mesh.getConnectivity().compute(0, 2);
  P1 space(mesh);
  GridFunction arrival(space);
  const Real upperSpeed = Real(1.2);
  const auto speed = [upperSpeed](const Point& point) {
    return point.getPolytope().getIndex() == 0 ? Real(1) : upperSpeed;
  };
  const Real exact =
    Real(1) / upperSpeed + Real(0.2) * std::sqrt(1 - Real(1) / (upperSpeed * upperSpeed));
  const auto verify = [&](const auto& coefficient) {
    Eikonal::FIM solver(arrival, coefficient);
    solver.setLocalTolerance(Real(0.25)).seed(IndexVector{1}).solve();
    const Real value = arrival[space.getDOFs(0, 0)[0]];
    EXPECT_NEAR(value, exact, 1e-12);
    EXPECT_LT(value, arrival[space.getDOFs(0, 2)[0]]);
    const auto coarse = arrival.getData();
    const size_t rounds = solver.getIterations();
    solver.setLocalTolerance(1e-10).seed(IndexVector{1}).solve();
    EXPECT_EQ(solver.getIterations(), rounds);
    EXPECT_TRUE((arrival.getData().array() == coarse.array()).all());
  };
  verify(speed);
  MatrixFunction tensor(2, 2, [&](const Point& point) {
    const Real value = speed(point);
    return (value * value) * Math::SpatialMatrix<Real>::Identity(2, 2);
  });
  verify(tensor);
}

/// @brief Triangular traces use the exact interior minimum independently of search resolution.
TEST(FIM, AnalyticTriangularTrace)
{
  auto mesh = LocalMesh::Builder()
                .initialize(3)
                .nodes(4)
                .vertex({0, 0, 0})
                .vertex({1, 0, 0})
                .vertex({0, 1, 0})
                .vertex({Real(0.2), Real(0.3), 1})
                .polytope(Polytope::Type::Tetrahedron, {0, 1, 2, 3})
                .finalize();
  mesh.getConnectivity().compute(3, 0);
  mesh.getConnectivity().compute(0, 3);
  const auto verify = [](const LocalMesh& current, const auto& coefficient, Real exact) {
    P1 space(current);
    GridFunction arrival(space);
    Eikonal::FIM solver(arrival, coefficient);
    solver.setLocalTolerance(Real(0.25)).seed(IndexVector{0, 1, 2}).solve();
    EXPECT_NEAR(arrival[space.getDOFs(0, 3)[0]], exact, 1e-12);
    const auto coarse = arrival.getData();
    solver.setLocalTolerance(1e-10).seed(IndexVector{0, 1, 2}).solve();
    EXPECT_TRUE((arrival.getData().array() == coarse.array()).all());
  };
  verify(mesh, [](const Point&) { return Real(2); }, Real(0.5));
  auto transform = Math::SpatialMatrix<Real>::Identity(3, 3);
  transform(0, 0) = Real(1.6);
  transform(1, 0) = Real(1.2);
  transform(0, 1) = Real(-0.6);
  transform(1, 1) = Real(0.8);
  transform(2, 2) = 3;
  for (auto vertex = mesh.getVertex(); vertex; ++vertex)
  {
    const Math::SpatialPoint physical = transform * vertex->getCoordinates();
    mesh.setVertexCoordinates(vertex->getIndex(), physical);
  }
  verify(mesh, MatrixFunction(transform * transform.transpose()), Real(1));
}
