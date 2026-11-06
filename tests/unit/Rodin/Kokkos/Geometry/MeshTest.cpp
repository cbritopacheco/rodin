/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <limits>
#include <type_traits>

#include <Rodin/Geometry.h>
#include <Rodin/IO.h>
#include <Rodin/Kokkos/Geometry.h>
#include <Rodin/Variational/H1/H1Element.h>
#include <Rodin/Variational/P1/P1Element.h>

#include "../TestContext.h"

namespace Rodin::Tests::Unit
{
  namespace DeviceGeometry = Rodin::Kokkos::Geometry;
  using HostMesh = Geometry::Mesh<Context::Local>;
  using DeviceMesh = Geometry::Mesh<TestContext>;
  using Polytope = Geometry::Polytope;
  static_assert(std::is_same_v<DeviceMesh::Context, TestContext>);

  /// @brief Host mesh whose virtual default cell map is independent of its vertices.
  class IdentityMesh final : public HostMesh
  {
    public:
      /**
       * @brief Reuses topology and coordinates without warming geometry caches.
       * @param source Uncached source mesh to move.
       */
      explicit IdentityMesh(HostMesh&& source)
        : HostMesh(std::move(source))
      {}

      /**
       * @brief Supplies an identity cell chart and the ordinary subentity charts.
       * @param dimension Entity reference dimension.
       * @param index Rank-local entity index.
       * @returns Newly allocated transformation owned by the caller.
       */
      Geometry::PolytopeTransformation* getDefaultPolytopeTransformation(
        size_t dimension, Index index) const override
      {
        if (dimension == getDimension())
          return new Geometry::IdentityTransformation(dimension);
        return HostMesh::getDefaultPolytopeTransformation(dimension, index);
      }
  };

  /**
   * @brief Constructs a nontrivial first-order physical image of a reference cell.
   * @param geometry Reference geometry family.
   * @param embedding Physical coordinate dimension.
   * @returns Single-cell host mesh, with vertices for point geometry.
   */
  HostMesh makeMesh(Polytope::Type geometry, size_t embedding)
  {
    const Polytope::Traits traits(geometry);
    const size_t dimension = traits.getDimension();
    HostMesh::Builder builder;
    builder.initialize(embedding).nodes(traits.getVertexCount());
    IndexArray vertices(traits.getVertexCount());
    for (Index a = 0; a < traits.getVertexCount(); ++a)
    {
      const auto reference = traits.getVertex(a);
      Math::SpatialPoint physical(static_cast<std::uint8_t>(embedding));
      for (size_t i = 0; i < embedding; ++i)
      {
        physical(i) = Real(i + 1);
        for (size_t j = 0; j < dimension; ++j)
        {
          physical(i) +=
            (i == j ? Real(1.5) : Real(0.1) * Real(i + j + 1)) * reference(j);
        }
      }
      if (dimension > 1)
        physical(0) += Real(0.1) * reference(0) * reference(1);
      builder.vertex(physical);
      vertices(a) = a;
    }
    if (dimension > 0)
      builder.polytope(geometry, vertices);
    return builder.finalize();
  }

  /**
   * @brief Compares numerical geometry on every entity through a device kernel.
   * @param source Host mesh supplying the independent geometry oracle.
   * @param mesh Resident mesh to evaluate.
   *
   * Two interior reference samples per entity exercise physical coordinates,
   * Jacobians and measures. Embedded entities use the host metric measure.
   */
  void checkGeometry(const HostMesh& source, const DeviceMesh& mesh)
  {
    constexpr Index sampleCount = 2;
    constexpr Real tolerance = 2e-12;
    const auto view = mesh.getDeviceView();
    const size_t embedding = source.getSpaceDimension();
    EXPECT_EQ(view.getDimension(), source.getDimension());
    EXPECT_EQ(view.getSpaceDimension(), embedding);
    ::Kokkos::View<Index*, TestMemorySpace> counts(
      "resident entity counts", source.getDimension() + 1);
    ::Kokkos::parallel_for(
      "read resident entity counts",
      ::Kokkos::RangePolicy<TestExecutionSpace>(
        mesh.getContext().getExecutionSpace(), 0, source.getDimension() + 1),
      KOKKOS_LAMBDA(Index dimension) { counts(dimension) = view.getCount(dimension); });
    const auto hostCounts =
      ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, counts);
    for (size_t dimension = 0; dimension <= source.getDimension(); ++dimension)
    {
      SCOPED_TRACE(dimension);
      const Index count = source.getPolytopeCount(dimension);
      ASSERT_EQ(hostCounts(dimension), count);
      ::Kokkos::View<Real**, TestMemorySpace> references(
        "geometry reference samples", count * sampleCount, 3);
      auto hostReferences = ::Kokkos::create_mirror_view(references);
      for (Index entity = 0; entity < count; ++entity)
      {
        const auto polytope = *source.getPolytope(dimension, entity);
        const Polytope::Traits traits(polytope.getGeometry());
        for (Index sample = 0; sample < sampleCount; ++sample)
        {
          auto reference = traits.getCentroid();
          if (sample != 0)
            reference = (reference + traits.getVertex(0)) / 2;
          for (size_t i = 0; i < dimension; ++i)
            hostReferences(entity * sampleCount + sample, i) = reference(i);
        }
      }
      ::Kokkos::deep_copy(references, hostReferences);
      const Index measureOffset = embedding + embedding * dimension;
      ::Kokkos::View<Real**, TestMemorySpace> results(
        "geometry sample results", count * sampleCount, measureOffset + 1);
      ::Kokkos::View<Index**, TestMemorySpace> extents(
        "geometry sample dimensions", count * sampleCount, 3);
      ::Kokkos::parallel_for(
        "evaluate resident mesh geometry",
        ::Kokkos::RangePolicy<TestExecutionSpace>(
          mesh.getContext().getExecutionSpace(), 0, count * sampleCount),
        KOKKOS_LAMBDA(Index item) {
          DeviceGeometry::Vector reference(dimension);
          for (size_t i = 0; i < dimension; ++i)
            reference(i) = references(item, i);
          const auto point =
            view.getPolytope(dimension, item / sampleCount).getPoint(reference);
          const auto physical = point.getPhysicalCoordinates();
          const auto jacobian = point.getJacobian();
          extents(item, 0) = physical.size();
          extents(item, 1) = jacobian.rows();
          extents(item, 2) = jacobian.cols();
          for (size_t i = 0; i < embedding; ++i)
          {
            results(item, i) = physical(i);
            for (size_t j = 0; j < dimension; ++j)
              results(item, embedding + i * dimension + j) = jacobian(i, j);
          }
          results(item, measureOffset) = point.getDistortion();
        });
      const auto hostResults =
        ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, results);
      const auto hostExtents =
        ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, extents);
      for (Index item = 0; item < count * sampleCount; ++item)
      {
        SCOPED_TRACE(item);
        Math::SpatialPoint reference(static_cast<std::uint8_t>(dimension));
        for (size_t i = 0; i < dimension; ++i)
          reference(i) = hostReferences(item, i);
        const auto entity = *source.getPolytope(dimension, item / sampleCount);
        const Geometry::Point expected(entity, reference);
        EXPECT_EQ(hostExtents(item, 0), embedding);
        EXPECT_EQ(hostExtents(item, 1), embedding);
        EXPECT_EQ(hostExtents(item, 2), dimension);
        for (size_t i = 0; i < embedding; ++i)
        {
          EXPECT_NEAR(
            hostResults(item, i), expected.getPhysicalCoordinates()(i), tolerance);
          for (size_t j = 0; j < dimension; ++j)
          {
            EXPECT_NEAR(hostResults(item, embedding + i * dimension + j),
              expected.getJacobian()(i, j), tolerance);
          }
        }
        EXPECT_NEAR(
          hostResults(item, measureOffset), expected.getDistortion(), tolerance);
      }
    }
  }

  /// @brief All supported cell families and embedding dimensions reproduce host maps.
  TEST(Kokkos_Mesh, GeometryMatchesEveryReferenceFamily)
  {
    const TestContext context;
    for (const auto geometry :
      {Polytope::Type::Point, Polytope::Type::Segment, Polytope::Type::Triangle,
        Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
        Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge})
    {
      const size_t dimension = Polytope::Traits(geometry).getDimension();
      for (size_t embedding = std::max(size_t(1), dimension); embedding <= 3; ++embedding)
      {
        SCOPED_TRACE(::testing::Message()
          << "geometry=" << int(geometry) << " embedding=" << embedding);
        auto source = makeMesh(geometry, embedding);
        for (size_t d = dimension; d > 0; --d)
          source.getConnectivity().compute(d, d - 1);
        DeviceMesh mesh(context, source);
        checkGeometry(source, mesh);
      }
    }
  }

  /// @brief Explicit identity maps remain identity regardless of vertex coordinates.
  TEST(Kokkos_Mesh, IdentityTransformationPreserved)
  {
    const TestContext context;
    auto source = makeMesh(Polytope::Type::Triangle, 2);
    source.setPolytopeTransformation({2, 0}, new Geometry::IdentityTransformation(2));
    DeviceMesh mesh(context, source);
    checkGeometry(source, mesh);
  }

  /// @brief Construction preserves uncached virtual geometry without slicing its source.
  TEST(Kokkos_Mesh, DerivedSourceTransformationPreserved)
  {
    const TestContext context;
    IdentityMesh source(makeMesh(Polytope::Type::Triangle, 2));
    DeviceMesh mesh(context, source);
    checkGeometry(source, mesh);
  }

  /// @brief Explicit first-order control coefficients are preserved independently of vertices.
  TEST(Kokkos_Mesh, ParametricControlPointsPreserved)
  {
    const TestContext context;
    auto source = makeMesh(Polytope::Type::Triangle, 2);
    Variational::RealP1Element element(Polytope::Type::Triangle);
    Geometry::PointCloud controls(2, element.getCount());
    for (Index a = 0; a < element.getCount(); ++a)
    {
      const auto reference = element.getNode(a);
      controls(0, a) = 5 + 3 * reference(0);
      controls(1, a) = 7 + 2 * reference(1) + Real(0.2) * reference(0);
    }
    source.setPolytopeTransformation({2, 0},
      new Geometry::ParametricTransformation<Variational::RealP1Element>(
        std::move(controls), std::move(element)));
    DeviceMesh mesh(context, source);
    checkGeometry(source, mesh);
  }

  /// @brief Packed incidence, ordered vertices, optional attributes and mixed geometry survive construction.
  void checkPreservesTopologyAndAttributes()
  {
    auto source = HostMesh::Builder()
                    .initialize(2)
                    .nodes(5)
                    .vertex({0, 0})
                    .vertex({1, 0})
                    .vertex({0, 1})
                    .vertex({1, 1})
                    .vertex({2, 0})
                    .polytope(Polytope::Type::Quadrilateral, {0, 1, 3, 2})
                    .polytope(Polytope::Type::Triangle, {1, 4, 3})
                    .finalize();
    for (size_t d = 0; d <= source.getDimension(); ++d)
    {
      for (size_t dp = 0; dp <= source.getDimension(); ++dp)
        source.getConnectivity().compute(d, dp);
    }
    for (size_t d = 0; d <= source.getDimension(); ++d)
    {
      for (Index i = 0; i < source.getPolytopeCount(d); ++i)
      {
        source.setAttribute({d, i},
          i % 2 == 0 ? Optional<Geometry::Attribute>{}
                     : Optional<Geometry::Attribute>{Geometry::Attribute(10 * d + i)});
      }
    }
    source.setAttribute({2, 1}, Geometry::Attribute(0));
    const TestContext context;
    DeviceMesh mesh(context, source);
    const auto view = mesh.getDeviceView();
    for (size_t d = 0; d <= source.getDimension(); ++d)
    {
      const Index count = source.getPolytopeCount(d);
      constexpr Index metadataSize = 5;
      constexpr Index vertexCapacity = RODIN_MAXIMUM_POLYTOPE_VERTICES;
      ::Kokkos::View<Index**, TestMemorySpace> entities(
        "resident entity metadata", count, metadataSize + vertexCapacity);
      ::Kokkos::parallel_for(
        "read resident entities",
        ::Kokkos::RangePolicy<TestExecutionSpace>(context.getExecutionSpace(), 0, count),
        KOKKOS_LAMBDA(Index i) {
          const auto entity = view.getPolytope(d, i);
          const auto vertices = entity.getVertices();
          entities(i, 0) = entity.getDimension();
          entities(i, 1) = entity.getIndex();
          entities(i, 2) = static_cast<Index>(entity.getGeometry());
          entities(i, 3) = entity.hasAttribute();
          entities(i, 4) = entity.hasAttribute() ? entity.getAttribute() : 0;
          for (Index a = 0; a < vertices.size(); ++a)
            entities(i, metadataSize + a) = vertices[a];
        });
      const auto hostEntities =
        ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, entities);
      for (Index i = 0; i < count; ++i)
      {
        const auto entity = *source.getPolytope(d, i);
        EXPECT_EQ(hostEntities(i, 0), d);
        EXPECT_EQ(hostEntities(i, 1), i);
        EXPECT_EQ(hostEntities(i, 2), static_cast<Index>(entity.getGeometry()));
        EXPECT_EQ(hostEntities(i, 3), bool(source.getAttribute(d, i)));
        if (const auto attribute = source.getAttribute(d, i))
          EXPECT_EQ(hostEntities(i, 4), *attribute);
        for (Index a = 0; a < entity.getVertices().size(); ++a)
          EXPECT_EQ(hostEntities(i, metadataSize + a), entity.getVertices()[a]);
      }
      for (size_t dp = 0; dp <= source.getDimension(); ++dp)
      {
        SCOPED_TRACE(::testing::Message() << "incidence=" << d << "," << dp);
        ::Kokkos::View<Index*, TestMemorySpace> presence("incidence presence", 1);
        ::Kokkos::parallel_for(
          "read resident incidence presence",
          ::Kokkos::RangePolicy<TestExecutionSpace>(context.getExecutionSpace(), 0, 1),
          KOKKOS_LAMBDA(
            Index) { presence(0) = view.getConnectivity().hasIncidence(d, dp); });
        const auto hostPresence =
          ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, presence);
        ASSERT_TRUE(hostPresence(0));
        Index stride = 1;
        for (Index i = 0; i < count; ++i)
        {
          stride = std::max(
            stride, Index(source.getConnectivity().getIncidence({d, dp}, i).size() + 1));
        }
        ::Kokkos::View<Index**, TestMemorySpace> incidences(
          "resident incidence results", count, stride);
        ::Kokkos::parallel_for(
          "read resident incidences",
          ::Kokkos::RangePolicy<TestExecutionSpace>(
            context.getExecutionSpace(), 0, count),
          KOKKOS_LAMBDA(Index i) {
            const auto incident = view.getConnectivity().getIncidence(d, dp, i);
            incidences(i, 0) = incident.size();
            for (Index a = 0; a < incident.size(); ++a)
              incidences(i, a + 1) = incident[a];
          });
        const auto hostIncidences =
          ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, incidences);
        for (Index i = 0; i < count; ++i)
        {
          const auto& expected = source.getConnectivity().getIncidence({d, dp}, i);
          ASSERT_EQ(hostIncidences(i, 0), expected.size());
          for (Index a = 0; a < expected.size(); ++a)
            EXPECT_EQ(hostIncidences(i, a + 1), expected[a]);
        }
      }
    }
    checkGeometry(source, mesh);
  }

  /// @brief Exercises PreservesTopologyAndAttributes with CUDA-compatible kernel scope.
  TEST(Kokkos_Mesh, PreservesTopologyAndAttributes)
  {
    checkPreservesTopologyAndAttributes();
  }

  /// @brief The specialization implements the existing virtual mesh interface directly.
  TEST(Kokkos_Mesh, ResidentMeshBaseAndImmediateUpdates)
  {
    static_assert(std::is_base_of_v<Geometry::MeshBase, DeviceMesh>);
    static_assert(!std::is_base_of_v<HostMesh, DeviceMesh>);
    const TestContext context;
    auto source = HostMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    source.getConnectivity().compute(1, 2);
    DeviceMesh mesh(context, source);
    Geometry::MeshBase& base = mesh;
    EXPECT_EQ(&base.getContext(), &mesh.getContext());
    EXPECT_EQ(base.getCellCount(), source.getCellCount());
    EXPECT_NEAR(base.getArea(), source.getArea(), 1e-12);
    EXPECT_NEAR(base.getPerimeter(), source.getPerimeter(), 1e-12);
    EXPECT_FALSE(base.isSubMesh());
    EXPECT_THROW(base.asSubMesh(), Alert::Exception);
    // Warm the ordinary host transformation/quadrature caches before mutation.
    const auto oldArea = base.getCell(0)->getMeasure();
    base.scale(2);
    source.scale(2);
    EXPECT_NEAR(base.getCell(0)->getMeasure(), 4 * oldArea, 1e-12);
    EXPECT_NEAR(base.getArea(), source.getArea(), 1e-12);
    checkGeometry(source, mesh);
    base.setVertexCoordinates(0, 0.25, 0);
    source.setVertexCoordinates(0, 0.25, 0);
    checkGeometry(source, mesh);
    base.setVertexCoordinates(1, Math::SpatialPoint({2.5, 0.1}));
    source.setVertexCoordinates(1, Math::SpatialPoint({2.5, 0.1}));
    checkGeometry(source, mesh);
    base.setAttribute({2, 0}, 0);
    source.setAttribute({2, 0}, 0);
    EXPECT_EQ(base.getAttribute(2, 0), Optional<Geometry::Attribute>(0));
    EXPECT_NEAR(base.getArea(0), source.getArea(0), 1e-12);
    EXPECT_NEAR(base.getArea(FlatSet<Geometry::Attribute>{0}),
      source.getArea(FlatSet<Geometry::Attribute>{0}), 1e-12);
    base.setAttribute({2, 0}, {});
    EXPECT_FALSE(base.getAttribute(2, 0));
    EXPECT_DOUBLE_EQ(base.getArea(0), 0);
  }

  /// @brief Copies allocate independent resident arrays; moves leave a usable empty source.
  TEST(Kokkos_Mesh, CopiesDoNotAliasResidentStorage)
  {
    const TestContext context;
    auto source = makeMesh(Polytope::Type::Triangle, 2);
    DeviceMesh mesh(context, source);
    DeviceMesh copied(mesh);
    EXPECT_NE(mesh.getDeviceView().coordinates, copied.getDeviceView().coordinates);
    EXPECT_NE(mesh.getDeviceView().vertexIndices, copied.getDeviceView().vertexIndices);
    mesh.scale(2);
    checkGeometry(source, copied);
    DeviceMesh assigned(context);
    assigned = copied;
    copied.scale(3);
    checkGeometry(source, assigned);
    DeviceMesh moved(std::move(assigned));
    checkGeometry(source, moved);
    EXPECT_TRUE(assigned.isEmpty());
    EXPECT_EQ(assigned.getDimension(), 0u);
    assigned = std::move(moved);
    checkGeometry(source, assigned);
    EXPECT_TRUE(moved.isEmpty());
  }

  /// @brief Construction owns its data independently of the source and builder.
  TEST(Kokkos_Mesh, BuilderAndSourceLifetime)
  {
    DeviceMesh mesh;
    {
      DeviceMesh::Builder builder;
      mesh = builder.initialize(2)
               .nodes(3)
               .vertex({0, 0})
               .vertex({1, 0})
               .vertex({0, 1})
               .polytope(Polytope::Type::Triangle, {0, 1, 2})
               .attribute({2, 0}, 7)
               .finalize();
    }
    EXPECT_NEAR(mesh.getArea(), 0.5, 1e-12);
    EXPECT_EQ(mesh.getAttribute(2, 0), Optional<Geometry::Attribute>(7));
    auto source = makeMesh(Polytope::Type::Triangle, 2);
    {
      HostMesh temporary(source);
      mesh = DeviceMesh(TestContext{}, temporary);
    }
    checkGeometry(source, mesh);
  }

  /// @brief All resident measure reductions reproduce the host quadrature policy.
  TEST(Kokkos_Mesh, MeasuresMatchEveryReferenceFamily)
  {
    for (const auto geometry :
      {Polytope::Type::Point, Polytope::Type::Segment, Polytope::Type::Triangle,
        Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
        Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge})
    {
      const size_t dimension = Polytope::Traits(geometry).getDimension();
      for (size_t embedding = std::max(size_t(1), dimension); embedding <= 3; ++embedding)
      {
        auto source = makeMesh(geometry, embedding);
        DeviceMesh mesh(TestContext{}, source);
        EXPECT_NEAR(mesh.getMeasure(dimension), source.getMeasure(dimension), 1e-11);
        mesh.scale(1.5);
        source.scale(1.5);
        EXPECT_NEAR(mesh.getMeasure(dimension), source.getMeasure(dimension), 1e-11);
      }
    }
  }

  /// @brief Unsupported maps and nonfinite values are rejected before resident data changes.
  TEST(Kokkos_Mesh, InvalidUpdatesPreserveResidentGeometry)
  {
    auto source = makeMesh(Polytope::Type::Triangle, 2);
    DeviceMesh mesh(TestContext{}, source);
    EXPECT_THROW(mesh.setVertexCoordinates(0, std::numeric_limits<Real>::quiet_NaN(), 0),
      Alert::Exception);
    Geometry::PointCloud points(2, 6);
    points.setZero();
    EXPECT_THROW(mesh.setPolytopeTransformation({2, 0},
                   new Geometry::ParametricTransformation(
                     points, Variational::H1Element<2, Real>(Polytope::Type::Triangle))),
      Alert::Exception);
    checkGeometry(source, mesh);
  }

  /// @brief Replacing compact control tables retains independent maps of other entities.
  TEST(Kokkos_Mesh, ResidentTransformationUpdates)
  {
    auto source = HostMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    DeviceMesh mesh(TestContext{}, source);
    for (Index i = 0; i < source.getCellCount(); ++i)
    {
      const auto& map = dynamic_cast<
        const Geometry::ParametricTransformation<Variational::RealP1Element>&>(
        source.getPolytopeTransformation(2, i));
      Geometry::PointCloud points(map.getPointMatrix());
      for (Index a = 0; a < points.getCount(); ++a)
        points(0, a) += Real(i + 1);
      mesh.setPolytopeTransformation({2, i},
        new Geometry::ParametricTransformation(
          points, Variational::RealP1Element(Polytope::Type::Triangle)));
      source.setPolytopeTransformation({2, i},
        new Geometry::ParametricTransformation(
          points, Variational::RealP1Element(Polytope::Type::Triangle)));
      checkGeometry(source, mesh);
      EXPECT_NEAR(mesh.getCell(i)->getMeasure(), source.getCell(i)->getMeasure(), 1e-12);
    }
    mesh.setPolytopeTransformation({2, 0}, new Geometry::IdentityTransformation(2));
    source.setPolytopeTransformation({2, 0}, new Geometry::IdentityTransformation(2));
    checkGeometry(source, mesh);
    mesh.setPolytopeTransformation({2, 1}, nullptr);
    source.setPolytopeTransformation({2, 1}, nullptr);
    checkGeometry(source, mesh);
    mesh.scale(2);
    source.scale(2);
    checkGeometry(source, mesh);
  }

  /// @brief Public file operations reuse Rodin formats and replace resident storage.
  TEST(Kokkos_Mesh, FileOperationsRoundTrip)
  {
    const auto directory = boost::filesystem::temp_directory_path() /
      boost::filesystem::unique_path("rodin-kokkos-mesh-%%%%-%%%%");
    ASSERT_TRUE(boost::filesystem::create_directories(directory));
    auto source = HostMesh::UniformGrid(Polytope::Type::Segment, {3});
    source.scale(0.5);
    source.setAttribute({1, 0}, 17);
    source.setAttribute({1, 1}, 23);
    source.getConnectivity().compute(0, 1);
    DeviceMesh mesh(TestContext{}, source);
    for (const auto format :
      {IO::FileFormat::MFEM, IO::FileFormat::MEDIT, IO::FileFormat::HDF5})
    {
      const auto path = directory / (std::to_string(static_cast<int>(format)) + ".mesh");
      mesh.save(path, format);
      DeviceMesh loaded(TestContext{}, makeMesh(Polytope::Type::Triangle, 2));
      loaded.load(path, format);
      EXPECT_EQ(loaded.getDimension(), 1u);
      EXPECT_EQ(loaded.getSpaceDimension(), 1u);
      EXPECT_NEAR(loaded.getMeasure(1), 1, 1e-12);
      EXPECT_NEAR(loaded.getMeasure(1, 17), 0.5, 1e-12);
      EXPECT_NEAR(loaded.getMeasure(1, 23), 0.5, 1e-12);
      const auto count = loaded.getCellCount();
      EXPECT_THROW(loaded.load(directory / "missing", format), Alert::Exception);
      EXPECT_EQ(loaded.getCellCount(), count);
      EXPECT_NEAR(loaded.getMeasure(1), 1, 1e-12);
    }
    EXPECT_GT(boost::filesystem::remove_all(directory), 0u);
  }

  /// @brief Device discovery creates usable mesh entities on every supported reference family.
  TEST(Kokkos_Mesh, DiscoveryBuildsResidentGeometry)
  {
    for (const auto geometry : {Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
           Polytope::Type::Tetrahedron, Polytope::Type::Hexahedron,
           Polytope::Type::Pyramid, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      const size_t dimension = Polytope::Traits(geometry).getDimension();
      auto source = makeMesh(geometry, 3);
      source.setAttribute({dimension, 0}, 19);
      DeviceMesh mesh(TestContext{}, source);
      for (size_t d = dimension; d > 1; --d)
      {
        source.getConnectivity().compute(d, d - 1);
        mesh.getConnectivity().compute(d, d - 1);
        checkGeometry(source, mesh);
        EXPECT_EQ(mesh.getAttribute(dimension, 0), Optional<Geometry::Attribute>(19));
        for (Index i = 0; i < mesh.getPolytopeCount(d - 1); ++i)
          EXPECT_FALSE(mesh.getAttribute(d - 1, i));
        EXPECT_NEAR(mesh.getMeasure(d - 1), source.getMeasure(d - 1), 1e-11);
      }
      source.getConnectivity().compute(dimension - 1, dimension);
      mesh.getConnectivity().compute(dimension - 1, dimension);
      EXPECT_NEAR(mesh.getPerimeter(), source.getPerimeter(), 1e-11);
      for (Index i = 0; i < mesh.getFaceCount(); ++i)
      {
        EXPECT_EQ(mesh.isBoundary(i), source.isBoundary(i));
        EXPECT_EQ(mesh.isInterface(i), source.isInterface(i));
      }
    }
  }

  /// @brief Growing edges and faces preserves attributes and independent maps at higher dimensions.
  TEST(Kokkos_Mesh, DiscoveryPreservesAttachedEntityData)
  {
    auto source = makeMesh(Polytope::Type::Hexahedron, 3);
    Geometry::PointCloud cellPoints(3, source.getVertexCount());
    for (Index i = 0; i < source.getVertexCount(); ++i)
    {
      for (size_t s = 0; s < 3; ++s)
        cellPoints(s, i) = 2 * source.getVertexCoordinates(i)(s);
    }
    source.setPolytopeTransformation({3, 0},
      new Geometry::ParametricTransformation(
        cellPoints, Variational::RealP1Element(Polytope::Type::Hexahedron)));
    source.setAttribute({3, 0}, 71);
    source.setAttribute({0, 2}, 0);
    DeviceMesh mesh(TestContext{}, source);
    // Warm caches before inserting dimensions ahead of the cell in resident storage.
    EXPECT_NEAR(mesh.getCell(0)->getMeasure(), source.getCell(0)->getMeasure(), 1e-11);
    source.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(3, 2);
    // A host transformation query also has to provision newly discovered slots.
    EXPECT_NEAR(mesh.getFace(0)->getMeasure(), source.getFace(0)->getMeasure(), 1e-11);
    const auto& faceMap =
      dynamic_cast<const Geometry::ParametricTransformation<Variational::RealP1Element>&>(
        source.getPolytopeTransformation(2, 0));
    Geometry::PointCloud facePoints(faceMap.getPointMatrix());
    for (Index i = 0; i < facePoints.getCount(); ++i)
      facePoints(0, i) *= 3;
    source.setPolytopeTransformation({2, 0},
      new Geometry::ParametricTransformation(
        facePoints, Variational::RealP1Element(Polytope::Type::Quadrilateral)));
    mesh.setPolytopeTransformation({2, 0},
      new Geometry::ParametricTransformation(
        facePoints, Variational::RealP1Element(Polytope::Type::Quadrilateral)));
    source.setAttribute({2, 0}, 42);
    mesh.setAttribute({2, 0}, 42);
    source.getConnectivity().compute(1, 3);
    mesh.getConnectivity().compute(1, 3);
    checkGeometry(source, mesh);
    EXPECT_EQ(mesh.getAttribute(3, 0), Optional<Geometry::Attribute>(71));
    EXPECT_EQ(mesh.getAttribute(2, 0), Optional<Geometry::Attribute>(42));
    EXPECT_EQ(mesh.getAttribute(0, 2), Optional<Geometry::Attribute>(0));
    mesh.setAttribute({1, 0}, 8);
    EXPECT_NEAR(mesh.getMeasure(1, 8), mesh.getPolytope(1, 0)->getMeasure(), 1e-11);
    mesh.scale(0.5);
    source.scale(0.5);
    checkGeometry(source, mesh);
  }

  /// @brief Copy and move materialize pending resident slots without requiring a user synchronization call.
  TEST(Kokkos_Mesh, CopyAndMoveAfterDiscovery)
  {
    auto source = makeMesh(Polytope::Type::Tetrahedron, 3);
    source.setAttribute({3, 0}, 7);
    DeviceMesh mesh(TestContext{}, source);
    source.getConnectivity().compute(2, 3);
    mesh.getConnectivity().compute(2, 3);
    DeviceMesh copied(mesh);
    checkGeometry(source, copied);
    EXPECT_EQ(copied.getAttribute(3, 0), Optional<Geometry::Attribute>(7));
    source.getConnectivity().compute(1, 2);
    mesh.getConnectivity().compute(1, 2);
    DeviceMesh moved(std::move(mesh));
    checkGeometry(source, moved);
    EXPECT_EQ(moved.getAttribute(3, 0), Optional<Geometry::Attribute>(7));
    EXPECT_TRUE(mesh.isEmpty());
  }

  /// @brief Empty meshes already have a resident topology and need no publication step.
  void checkEmptyMesh()
  {
    DeviceMesh mesh;
    EXPECT_TRUE(mesh.isEmpty());
    EXPECT_DOUBLE_EQ(mesh.getArea(), 0);
    EXPECT_DOUBLE_EQ(mesh.getVolume(), 0);
    EXPECT_DOUBLE_EQ(mesh.getPerimeter(), 0);
    const auto view = mesh.getDeviceView();
    ::Kokkos::View<Index*, TestMemorySpace> result("empty count", 1);
    ::Kokkos::parallel_for(
      "empty mesh query", ::Kokkos::RangePolicy<TestExecutionSpace>(0, 1),
      KOKKOS_LAMBDA(Index) { result(0) = view.getCount(0); });
    const auto host =
      ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, result);
    EXPECT_EQ(host(0), 0u);
  }

  /// @brief Exercises EmptyMesh with CUDA-compatible kernel scope.
  TEST(Kokkos_Mesh, EmptyMesh)
  {
    checkEmptyMesh();
  }
}
