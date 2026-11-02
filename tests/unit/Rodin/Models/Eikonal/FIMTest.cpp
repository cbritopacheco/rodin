/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <stdexcept>

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

/// @brief Higher-order maps are used along paths for all native reference families.
TEST(FIM, CurvedNativeCellsWithPullbackMetric)
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
      {
        controls(coordinate, node) = reference(coordinate) +
          Real(0.2) * reference(coordinate) * reference(coordinate);
      }
    }
    mesh.setPolytopeTransformation(
      {dimension, 0}, new ParametricTransformation(controls, element));
    P1 space(mesh);
    GridFunction arrival(space);
    MatrixFunction metric(dimension, dimension, [](const Point& point) {
      const auto& jacobian = point.getJacobian();
      return jacobian * jacobian.transpose();
    });
    Eikonal::FIM solver(arrival, metric);
    solver.seed(planeSeeds(geometry)).solve();
    const Polytope::Traits traits(geometry);
    for (Index vertex = 0; vertex < traits.getVertexCount(); ++vertex)
    {
      EXPECT_NEAR(
        arrival[space.getDOFs(0, vertex)[0]], traits.getVertex(vertex)(0), 1e-8);
    }
  }
}

/// @brief A curved embedded segment uses arc length rather than its endpoint chord.
TEST(FIM, EmbeddedCurvedSegmentAndQuadratureRefinement)
{
  auto mesh = LocalMesh::Builder()
                .initialize(2)
                .nodes(2)
                .vertex({0, 0})
                .vertex({1, Real(0.5)})
                .polytope(Polytope::Type::Segment, {0, 1})
                .finalize();
  mesh.getConnectivity().compute(1, 0);
  mesh.getConnectivity().compute(0, 1);
  RealH1Element<2> element(Polytope::Type::Segment);
  PointCloud controls(2, element.getCount());
  for (size_t node = 0; node < element.getCount(); ++node)
  {
    const Real coordinate = element.getNode(node)(0);
    controls(0, node) = coordinate;
    controls(1, node) = Real(0.5) * coordinate * coordinate;
  }
  mesh.setPolytopeTransformation({1, 0}, new ParametricTransformation(controls, element));
  P1 space(mesh);
  GridFunction arrival(space);
  Eikonal::FIM solver(arrival, [](const Point&) { return Real(1); });
  const IndexVector seeds{0};
  const Real exact = (std::sqrt(Real(2)) + std::asinh(Real(1))) / 2;
  solver.setIntegrationOrder(1).seed(seeds).solve();
  const Real coarseError = std::abs(arrival[space.getDOFs(0, 1)[0]] - exact);
  solver.setIntegrationOrder(20).seed(seeds).solve();
  const Real fineError = std::abs(arrival[space.getDOFs(0, 1)[0]] - exact);
  EXPECT_GT(coarseError, 1e-3);
  EXPECT_LT(fineError, 1e-12);
  EXPECT_GT(arrival[space.getDOFs(0, 1)[0]], std::sqrt(Real(1.25)));
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
  EXPECT_THROW(solver.setIntegrationOrder(0), Alert::Exception);
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

/// @brief Triangle and quadrilateral paths follow a curved surface with an ambient tensor.
TEST(FIM, EmbeddedCurvedSurfaces)
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
      builder.vertex({r(0), r(1), Real(0.2) * r(0) * r(1)});
      vertices[vertex] = vertex;
    }
    builder.polytope(geometry, vertices);
    auto mesh = builder.finalize();
    mesh.getConnectivity().compute(2, 0);
    mesh.getConnectivity().compute(0, 2);
    RealH1Element<2> element(geometry);
    PointCloud controls(3, element.getCount());
    for (size_t node = 0; node < element.getCount(); ++node)
    {
      const auto& r = element.getNode(node);
      controls(0, node) = r(0);
      controls(1, node) = r(1);
      controls(2, node) = Real(0.2) * r(0) * r(1);
    }
    mesh.setPolytopeTransformation(
      {2, 0}, new ParametricTransformation(controls, element));
    P1 space(mesh);
    GridFunction arrival(space);
    MatrixFunction coefficient(3, 3, [](const Point& point) {
      const auto& j = point.getJacobian();
      Math::SpatialVector<Real> normal{j(1, 0) * j(2, 1) - j(2, 0) * j(1, 1),
        j(2, 0) * j(0, 1) - j(0, 0) * j(2, 1), j(0, 0) * j(1, 1) - j(1, 0) * j(0, 1)};
      Math::SpatialMatrix<Real> tensor = j * j.transpose();
      for (size_t row = 0; row < 3; ++row)
      {
        for (size_t column = 0; column < 3; ++column)
          tensor(row, column) += normal(row) * normal(column);
      }
      return tensor;
    });
    Eikonal::FIM solver(arrival, coefficient);
    solver.seed(planeSeeds(geometry)).solve();
    for (Index vertex = 0; vertex < traits.getVertexCount(); ++vertex)
    {
      EXPECT_NEAR(
        arrival[space.getDOFs(0, vertex)[0]], traits.getVertex(vertex)(0), 1e-8);
    }
  }
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
  }
}

/// @brief Curved quadrilateral maps and tensor expressions remain deterministic across workers.
TEST(FIM, OpenMPCurvedTensorPaths)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, {4, 4});
  mesh.scale(Real(1) / 3);
  mesh.getConnectivity().compute(2, 0);
  mesh.getConnectivity().compute(0, 2);
  RealH1Element<2> element(Polytope::Type::Quadrilateral);
  for (auto cell = mesh.getCell(); cell; ++cell)
  {
    PointCloud controls(2, element.getCount());
    for (size_t node = 0; node < element.getCount(); ++node)
    {
      const auto& r = element.getNode(node);
      const Point point(*cell, r);
      const auto physical = point.getPhysicalCoordinates();
      controls(0, node) = physical(0);
      controls(1, node) =
        physical(1) + Real(0.05) * r(0) * (1 - r(0)) * r(1) * (1 - r(1));
    }
    mesh.setPolytopeTransformation(
      {2, cell->getIndex()}, new ParametricTransformation(controls, element));
  }
  P1 space(mesh);
  GridFunction arrival(space);
  Math::SpatialMatrix<Real> tensor = Math::SpatialMatrix<Real>::Identity(2, 2);
  tensor(0, 0) = 4;
  tensor(0, 1) = tensor(1, 0) = Real(0.5);
  Eikonal::FIM solver(arrival, MatrixFunction(tensor));
  solver.seed(IndexVector{0}).solve();
  const auto reference = arrival.getData();
  const size_t iterations = solver.getIterations();
  solver.setThreadCount(4).seed(IndexVector{0}).solve();
  EXPECT_EQ(solver.getIterations(), iterations);
  EXPECT_TRUE((arrival.getData().array() == reference.array()).all());
}

/// @brief Worker exceptions are propagated after joining, without terminating OpenMP workers.
TEST(FIM, OpenMPExceptions)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, {4, 4});
  mesh.getConnectivity().compute(2, 0);
  mesh.getConnectivity().compute(0, 2);
  P1 space(mesh);
  GridFunction arrival(space);
  Eikonal::FIM solver(arrival, [](const Point&) { return Real(-1); });
  solver.setThreadCount(4).seed(IndexVector{0});
  EXPECT_THROW(solver.solve(), Alert::Exception);
  EXPECT_THROW(solver.setThreadCount(0), Alert::Exception);
}

/// @brief A coefficient copy failure is caught even before a worker has an expression.
TEST(FIM, OpenMPCoefficientCopyException)
{
  struct Coefficient
  {
    bool* fail;
    explicit Coefficient(bool& flag)
      : fail(&flag)
    {}
    Coefficient(const Coefficient& other)
      : fail(other.fail)
    {
      if (*fail)
        throw std::runtime_error("User coefficient copy failed");
    }
    Real operator()(const Point&) const
    {
      return 1;
    }
  };
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
  mesh.getConnectivity().compute(2, 0);
  mesh.getConnectivity().compute(0, 2);
  P1 space(mesh);
  GridFunction arrival(space);
  bool fail = false;
  Eikonal::FIM solver(arrival, Coefficient(fail));
  solver.setThreadCount(4).seed(IndexVector{0});
  fail = true;
  EXPECT_THROW(solver.solve(), std::runtime_error);
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
