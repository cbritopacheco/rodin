/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <type_traits>
#include <Rodin/Assembly/Sequential.h>
#include <Rodin/Kokkos.h>

#include "../TestContext.h"

namespace Rodin::Tests::Unit
{
  using IterationHostMesh = Geometry::Mesh<Context::Local>;
  using IterationMesh = Geometry::Mesh<TestContext>;
  using ResidentIteration = Assembly::KokkosIteration<IterationMesh>;
  using IterationGeometry = Geometry::Polytope::Type;
  using Region = Geometry::Region;
  static_assert(std::is_trivially_copyable_v<ResidentIteration>);

  /**
   * @brief Compares resident traversal with the existing sequential strategy.
   * @param host Independent CPU mesh with matching entity numbering.
   * @param mesh Resident mesh to traverse through the iteration interface.
   */
  void checkIteration(const IterationHostMesh& host, const IterationMesh& mesh)
  {
    for (const auto region :
      {Region::Cells, Region::Faces, Region::Boundary, Region::Interface})
    {
      SCOPED_TRACE(static_cast<int>(region));
      const Assembly::SequentialIteration expected(host, region);
      const Assembly::KokkosIteration iteration(mesh, region);
      static_assert(std::is_same_v<decltype(iteration), const ResidentIteration>);
      ASSERT_EQ(iteration.getDimension(), expected.getDimension());
      ASSERT_EQ(iteration.getCount(), expected.getCount());

      constexpr Index vertexOffset = 6;
      ::Kokkos::View<Index**, TestMemorySpace> entities(
        "iteration entities", iteration.getCount(), vertexOffset + 8);
      ::Kokkos::View<Geometry::Attribute*, TestMemorySpace> attributes(
        "iteration attributes", iteration.getCount());
      ::Kokkos::parallel_for(
        "read iteration candidates",
        ::Kokkos::RangePolicy<TestExecutionSpace>(
          mesh.getContext().getExecutionSpace(), 0, iteration.getCount()),
        KOKKOS_LAMBDA(Index i) {
          const auto polytope = iteration.getPolytope(i);
          const auto vertices = polytope.getVertices();
          entities(i, 0) = iteration.filter(i);
          entities(i, 1) = polytope.getDimension();
          entities(i, 2) = polytope.getIndex();
          entities(i, 3) = static_cast<Index>(polytope.getGeometry());
          entities(i, 4) = vertices.size();
          entities(i, 5) = polytope.hasAttribute();
          if (polytope.hasAttribute())
            attributes(i) = polytope.getAttribute();
          for (Index v = 0; v < vertices.size(); ++v)
            entities(i, vertexOffset + v) = vertices[v];
        });
      mesh.getContext().getExecutionSpace().fence();
      const auto values =
        ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, entities);
      const auto labels =
        ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, attributes);
      for (Index i = 0; i < iteration.getCount(); ++i)
      {
        SCOPED_TRACE(i);
        const auto polytope = expected.getPolytope(i);
        const auto& vertices = polytope.getVertices();
        const auto attribute = polytope.getAttribute();
        EXPECT_EQ(values(i, 0), expected.filter(i));
        EXPECT_EQ(values(i, 1), polytope.getDimension());
        EXPECT_EQ(values(i, 2), polytope.getIndex());
        EXPECT_EQ(values(i, 3), static_cast<Index>(polytope.getGeometry()));
        ASSERT_EQ(values(i, 4), vertices.size());
        EXPECT_EQ(bool(values(i, 5)), bool(attribute));
        if (attribute)
        {
          EXPECT_EQ(labels(i), *attribute);
        }
        for (Index v = 0; v < vertices.size(); ++v)
          EXPECT_EQ(values(i, vertexOffset + v), vertices[v]);
      }

      Index selected = 0;
      ::Kokkos::parallel_reduce(
        "filter iteration candidates",
        ::Kokkos::RangePolicy<TestExecutionSpace>(
          mesh.getContext().getExecutionSpace(), 0, iteration.getCount()),
        KOKKOS_LAMBDA(Index i, Index & sum) {
          if (iteration.filter(i))
            ++sum;
        },
        selected);
      Index expectedSelected = 0;
      for (Index i = 0; i < expected.getCount(); ++i)
      {
        if (expected.filter(i))
          ++expectedSelected;
      }
      EXPECT_EQ(selected, expectedSelected);
    }
  }

  /**
   * @brief Discovers face relations independently and attaches optional labels.
   * @param host CPU mesh supplying the existing semantics.
   */
  void checkDiscoveredIteration(IterationHostMesh host)
  {
    IterationMesh mesh(TestContext{}, host);
    const size_t dimension = host.getDimension();
    host.getConnectivity().compute(dimension - 1, dimension);
    mesh.getConnectivity().compute(dimension - 1, dimension);
    for (size_t d = dimension - 1; d <= dimension; ++d)
    {
      for (Index i = 0; i < host.getPolytopeCount(d); ++i)
      {
        const Optional<Geometry::Attribute> attribute =
          i % 2 == 0 ? Optional<Geometry::Attribute>(i) : Optional<Geometry::Attribute>();
        host.setAttribute({d, i}, attribute);
        mesh.setAttribute({d, i}, attribute);
      }
    }
    checkIteration(host, mesh);
  }

  /// @brief Every reference family reuses the indexed region semantics after discovery.
  /**
   * @brief Counts candidates accepted by the iteration filter.
   * NVCC rejects extended lambdas in GoogleTest's private TestBody, so device
   * reductions used by test bodies live in free functions.
   * @param mesh Resident mesh supplying the execution instance.
   * @param iteration Resident traversal.
   * @returns Number of selected candidates.
   */
  Index countSelected(const IterationMesh& mesh, const ResidentIteration& iteration)
  {
    Index selected = 0;
    ::Kokkos::parallel_reduce(
      "count selected",
      ::Kokkos::RangePolicy<TestExecutionSpace>(
        mesh.getContext().getExecutionSpace(), 0, iteration.getCount()),
      KOKKOS_LAMBDA(Index i, Index & sum) {
        if (iteration.filter(i))
          ++sum;
      },
      selected);
    return selected;
  }

  /**
   * @brief Counts selected cells with three vertices.
   * @param mesh Resident mesh supplying the execution instance.
   * @param cells Resident cell traversal.
   * @returns Number of selected triangles.
   */
  Index countTriangles(const IterationMesh& mesh, const ResidentIteration& cells)
  {
    Index visited = 0;
    ::Kokkos::parallel_reduce(
      "count triangles",
      ::Kokkos::RangePolicy<TestExecutionSpace>(
        mesh.getContext().getExecutionSpace(), 0, cells.getCount()),
      KOKKOS_LAMBDA(Index i, Index & sum) {
        if (cells.filter(i) && cells.getPolytope(i).getVertices().size() == 3)
          ++sum;
      },
      visited);
    return visited;
  }

  /**
   * @brief Counts selected point cells whose index matches the candidate position.
   * @param mesh Resident mesh supplying the execution instance.
   * @param points Resident traversal of a point mesh.
   * @returns Number of consistent point cells.
   */
  Index countPoints(const IterationMesh& mesh, const ResidentIteration& points)
  {
    Index selected = 0;
    ::Kokkos::parallel_reduce(
      "count points",
      ::Kokkos::RangePolicy<TestExecutionSpace>(
        mesh.getContext().getExecutionSpace(), 0, points.getCount()),
      KOKKOS_LAMBDA(Index i, Index & sum) {
        const auto entity = points.getPolytope(i);
        if (points.filter(i) && entity.getDimension() == 0 && entity.getIndex() == i)
          ++sum;
      },
      selected);
    return selected;
  }

  TEST(Kokkos_Iteration, ReferenceFamiliesAndRegionFilters)
  {
    checkDiscoveredIteration(
      IterationHostMesh::UniformGrid(IterationGeometry::Segment, {4}));
    for (const auto geometry :
      {IterationGeometry::Triangle, IterationGeometry::Quadrilateral})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      checkDiscoveredIteration(IterationHostMesh::UniformGrid(geometry, {3, 3}));
    }
    for (const auto geometry :
      {IterationGeometry::Tetrahedron, IterationGeometry::Hexahedron,
        IterationGeometry::Wedge, IterationGeometry::Pyramid})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      checkDiscoveredIteration(IterationHostMesh::UniformGrid(geometry, {2, 2, 2}));
    }
  }

  /// @brief Embedded mixed cells preserve ordered vertices, candidate IDs and optional labels.
  TEST(Kokkos_Iteration, EmbeddedMixedMesh)
  {
    checkDiscoveredIteration(IterationHostMesh::Builder()
        .initialize(3)
        .nodes(5)
        .vertex({0, 0, 1})
        .vertex({1, 0, 1})
        .vertex({0, 1, 1})
        .vertex({1, 1, 1})
        .vertex({2, 0, 1})
        .polytope(IterationGeometry::Quadrilateral, {0, 1, 3, 2})
        .polytope(IterationGeometry::Triangle, {1, 4, 3})
        .finalize());
  }

  /// @brief Interface includes faces shared by more than two cells, as in LocalMesh.
  TEST(Kokkos_Iteration, NonmanifoldInterface)
  {
    auto host = IterationHostMesh::Builder()
                  .initialize(2)
                  .nodes(5)
                  .vertex({0, 0})
                  .vertex({1, 0})
                  .vertex({0, 1})
                  .vertex({0, -1})
                  .vertex({0.5, 2})
                  .polytope(IterationGeometry::Triangle, {0, 1, 2})
                  .polytope(IterationGeometry::Triangle, {1, 0, 3})
                  .polytope(IterationGeometry::Triangle, {0, 1, 4})
                  .finalize();
    host.getConnectivity().compute(1, 2);
    IterationMesh mesh(TestContext{}, host);
    checkIteration(host, mesh);
    const Assembly::KokkosIteration iteration(mesh, Region::Interface);
    EXPECT_EQ(iteration.getCount(), 7);
    const Index selected = countSelected(mesh, iteration);
    EXPECT_EQ(selected, 1);
  }

  /// @brief Cell and face traversal do not compute incidence; filtered regions require it explicitly.
  TEST(Kokkos_Iteration, ExplicitIncidenceRequirements)
  {
    auto host = IterationHostMesh::Builder()
                  .initialize(2)
                  .nodes(3)
                  .vertex({0, 0})
                  .vertex({1, 0})
                  .vertex({0, 1})
                  .polytope(IterationGeometry::Triangle, {0, 1, 2})
                  .finalize();
    IterationMesh mesh(TestContext{}, host);
    const Assembly::KokkosIteration cells(mesh, Region::Cells);
    EXPECT_EQ(cells.getCount(), 1);
    EXPECT_EQ(mesh.getFaceCount(), 0);
    EXPECT_FALSE(mesh.getConnectivity().hasIncidence(1, 2));
    const Assembly::KokkosIteration noFaces(mesh, Region::Boundary);
    EXPECT_EQ(noFaces.getCount(), 0);
    const Index visited = countTriangles(mesh, cells);
    EXPECT_EQ(visited, 1);

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().clear(1, 2);
    const Assembly::KokkosIteration faces(mesh, Region::Faces);
    EXPECT_EQ(faces.getCount(), 3);
    EXPECT_THROW((Assembly::KokkosIteration(mesh, Region::Boundary)), Alert::Exception);
    EXPECT_THROW((Assembly::KokkosIteration(mesh, Region::Interface)), Alert::Exception);
    EXPECT_FALSE(mesh.getConnectivity().hasIncidence(1, 2));
    mesh.getConnectivity().compute(1, 2);
    const Assembly::KokkosIteration boundary(mesh, Region::Boundary);
    const Assembly::KokkosIteration interface(mesh, Region::Interface);
    EXPECT_EQ(boundary.getCount(), faces.getCount());
    EXPECT_EQ(interface.getCount(), faces.getCount());
  }

  /// @brief Empty selections remain valid scans without eager boundary/interface iterator construction.
  TEST(Kokkos_Iteration, EmptySelection)
  {
    auto host = IterationHostMesh::Builder()
                  .initialize(2)
                  .nodes(3)
                  .vertex({0, 0})
                  .vertex({1, 0})
                  .vertex({0, 1})
                  .polytope(IterationGeometry::Triangle, {0, 1, 2})
                  .finalize();
    host.getConnectivity().compute(1, 2);
    IterationMesh mesh(TestContext{}, host);
    checkIteration(host, mesh);
    const Assembly::KokkosIteration interface(mesh, Region::Interface);
    EXPECT_GT(interface.getCount(), 0);
    const Index selected = countSelected(mesh, interface);
    EXPECT_EQ(selected, 0);
  }

  /// @brief Cell traversal handles empty and point meshes without face-dimension underflow.
  TEST(Kokkos_Iteration, ZeroDimensionalMeshesAndInvalidRegion)
  {
    IterationHostMesh empty;
    IterationMesh mesh(TestContext{}, empty);
    const Assembly::KokkosIteration cells(mesh, Region::Cells);
    EXPECT_EQ(cells.getDimension(), 0);
    EXPECT_EQ(cells.getCount(), 0);
    const Index emptyCount = countSelected(mesh, cells);
    EXPECT_EQ(emptyCount, 0);
    EXPECT_THROW((Assembly::KokkosIteration(mesh, Region::Faces)), Alert::Exception);
    EXPECT_THROW((Assembly::KokkosIteration(mesh, Region::Boundary)), Alert::Exception);
    EXPECT_THROW((Assembly::KokkosIteration(mesh, Region::Interface)), Alert::Exception);
    EXPECT_THROW(
      (Assembly::KokkosIteration(mesh, static_cast<Region>(-1))), Alert::Exception);

    auto point =
      IterationHostMesh::Builder().initialize(3).nodes(1).vertex({1, 2, 3}).finalize();
    mesh = IterationMesh(TestContext{}, point);
    const Assembly::KokkosIteration points(mesh, Region::Cells);
    EXPECT_EQ(points.getCount(), 1);
    const Index selected = countPoints(mesh, points);
    EXPECT_EQ(selected, 1);
  }
}
