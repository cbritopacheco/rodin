/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#include <cmath>
#include <array>
#include <cmath>
#include <type_traits>
#include <gtest/gtest.h>
#include <boost/filesystem.hpp>
#include <boost/mpi/collectives/all_reduce.hpp>
#include <Rodin/Geometry.h>
#include <Rodin/Kokkos/MPI.h>
#include <Rodin/MPI/IO/HDF5.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include "../TestContext.h"

#ifdef PETSC_VERSION_MAJOR
#error "Distributed Kokkos mesh and space support must not include PETSc."
#endif

namespace Rodin::Tests::Unit
{
  boost::mpi::environment* meshEnvironment = nullptr;
  boost::mpi::communicator* meshWorld = nullptr;

  /// @brief Owned cells counted through a borrowed device view.
  struct OwnedCells
  {
    /// @brief Borrowed local mesh and ownership flags.
    Kokkos::View<Geometry::Mesh<Context::Kokkos<Context::MPI, TestBackend>>> mesh;

    /// @param cell Local cell. @param count Number of owned cells.
    KOKKOS_INLINE_FUNCTION void operator()(Index cell, size_t& count) const
    {
      if (mesh.isOwned(mesh.dimension, cell))
        ++count;
    }
  };

  /// @brief Shard loading, numbering and device traversal need no PETSc runtime.
  TEST(Kokkos_MPI_Mesh, ShardRoundTripAndSpaces)
  {
    Context::MPI context(*meshEnvironment, *meshWorld);
    Geometry::Sharder<Context::MPI> sharder(context);
    if (meshWorld->rank() == 0)
    {
      auto local = Geometry::Mesh<Context::Local>::UniformGrid(
        Geometry::Polytope::Type::Triangle, {5, 5});
      local.getConnectivity().compute(2, 2);
      local.getConnectivity().compute(1, 2);
      Geometry::BalancedCompactPartitioner partitioner(local);
      partitioner.partition(static_cast<size_t>(meshWorld->size()));
      sharder.shard(partitioner);
      sharder.scatter(0);
    }
    const auto partition = sharder.gather(0);
    const auto file = boost::filesystem::temp_directory_path() /
      boost::filesystem::unique_path("rodin-kokkos-shard-%%%%-%%%%.h5");
    partition.save(file, IO::FileFormat::HDF5);
    Geometry::Mesh<Context::MPI> host(context);
    host.load(file, IO::FileFormat::HDF5);
    boost::filesystem::remove(file);
    host.getConnectivity().compute(2, 2);
    host.getConnectivity().compute(2, 0);
    host.getConnectivity().compute(2, 1);
    host.getConnectivity().compute(1, 2);
    host.getConnectivity().compute(1, 0);
    host.reconcile(1);
    const auto faceCells = host.getConnectivity().getIncidence(1, 2);
    host.getConnectivity().clear(1, 2);

    using DistributedContext = Context::Kokkos<Context::MPI, TestBackend>;
    Geometry::Mesh<DistributedContext> mesh(DistributedContext(context), host);
    EXPECT_NE(&mesh.getHost(), &host);
    Variational::P1 scalar(mesh);
    Variational::H1 quadratic(std::integral_constant<size_t, 2>{}, mesh);
    Variational::P0g constant(mesh);
    EXPECT_EQ(scalar.getSize(), 25u);
    EXPECT_EQ(quadratic.getSize(), 81u);
    EXPECT_EQ(constant.getSize(), 1u);
    const auto cellDOFs = quadratic.getDOFs(2, 0);
    EXPECT_FALSE(mesh.getConnectivity().hasIncidence(1, 2));
    mesh.getConnectivity().compute(1, 2);
    EXPECT_EQ(mesh.getConnectivity().getIncidence(1, 2), faceCells);
    EXPECT_TRUE((quadratic.getDOFs(2, 0) == cellDOFs).all());
    size_t localInterface = 0;
    for (auto face = mesh.getInterface(); face; ++face)
      ++localInterface;
    EXPECT_EQ(
      boost::mpi::all_reduce(*meshWorld, localInterface, std::plus<size_t>()), 40u);
    const auto query = scalar.getGhostedPositions();
    EXPECT_EQ(query.positions.extent(0), scalar.getLocalSize());
    const auto quadraticQuery = quadratic.getGhostedPositions();
    EXPECT_EQ(quadraticQuery.positions.extent(0), quadratic.getLocalSize());
    const auto constantQuery = constant.getGhostedPositions();
    EXPECT_EQ(constantQuery.positions.extent(0), constant.getLocalSize());
    const auto copy = scalar;
    EXPECT_NE(copy.getRevision(), scalar.getRevision());
    EXPECT_EQ(scalar.getLocalSize(), mesh.getLocal().getVertexCount());
    size_t localCount = 0;
    ::Kokkos::parallel_reduce("Count owned distributed cells",
      ::Kokkos::RangePolicy<TestExecutionSpace>(
        mesh.getContext().getExecutionSpace(), 0, mesh.getLocal().getCellCount()),
      OwnedCells{mesh.getView()}, localCount);
    const size_t count =
      boost::mpi::all_reduce(*meshWorld, localCount, std::plus<size_t>());
    EXPECT_EQ(count, 32u);
    EXPECT_EQ(mesh.getCellCount(), host.getShard().getCellCount());

    const auto original = host.getShard().getVertexCoordinates(0);
    const auto attribute = host.getShard().getAttribute(2, 0);
    const auto revision = mesh.getRevision();
    mesh.setAttribute({2, 0}, Geometry::Attribute(77));
    EXPECT_NE(mesh.getRevision(), revision);
    EXPECT_EQ(mesh.getAttribute(2, 0), Geometry::Attribute(77));
    EXPECT_EQ(mesh.getHost().getAttribute(2, 0), Geometry::Attribute(77));
    EXPECT_EQ(host.getShard().getAttribute(2, 0), attribute);
    mesh.scale(2);
    EXPECT_LE((mesh.getVertexCoordinates(0) - Real(2) * original).norm(), 1e-14);
    EXPECT_LE(
      (mesh.getHost().getVertexCoordinates(0) - Real(2) * original).norm(), 1e-14);
    EXPECT_EQ((host.getShard().getVertexCoordinates(0) - original).norm(), 0);
    for (Index i = 0; i < mesh.getVertexCount(); ++i)
    {
      auto point = mesh.getVertexCoordinates(i);
      point(0) += Real(0.125);
      mesh.setVertexCoordinates(i, point);
      mesh.setVertexCoordinates(i, point(1) - Real(0.25), 1);
      EXPECT_EQ(
        (mesh.getHost().getVertexCoordinates(i) - mesh.getVertexCoordinates(i)).norm(),
        0);
    }
    // Snapshot access is rank-local: the other ranks do not participate.
    if (meshWorld->rank() == 0)
    {
      auto snapshot = mesh.download();
      EXPECT_EQ(snapshot.getShard().getCellCount(), mesh.getCellCount());
      EXPECT_EQ(snapshot.getShard().getAttribute(2, 0), Geometry::Attribute(77));
      EXPECT_EQ(
        (snapshot.getVertexCoordinates(0) - mesh.getVertexCoordinates(0)).norm(), 0);
      snapshot.scale(3);
      EXPECT_NE(
        (snapshot.getVertexCoordinates(0) - mesh.getVertexCoordinates(0)).norm(), 0);
      auto copy = mesh;
      copy.setAttribute({2, 0}, Geometry::Attribute(88));
      EXPECT_EQ(mesh.getHost().getAttribute(2, 0), Geometry::Attribute(77));
      copy.getConnectivity().clear(1, 2);
      EXPECT_TRUE(mesh.getConnectivity().hasIncidence(1, 2));
      auto moved = std::move(copy);
      EXPECT_FALSE(moved.getConnectivity().hasIncidence(1, 2));
      EXPECT_EQ(moved.getAttribute(2, 0), Geometry::Attribute(88));
      moved.upload(host);
      EXPECT_EQ(moved.getAttribute(2, 0), attribute);
      EXPECT_EQ((moved.getVertexCoordinates(0) - original).norm(), 0);
      moved.getConnectivity().compute(1, 2);
      EXPECT_EQ(moved.getConnectivity().getIncidence(1, 2), faceCells);
    }
    meshWorld->barrier();
    host.save(file, IO::FileFormat::HDF5);
    Geometry::Mesh<DistributedContext> loaded{DistributedContext(context)};
    loaded.load(file, IO::FileFormat::HDF5);
    boost::filesystem::remove(file);
    EXPECT_EQ(loaded.getCellCount(), host.getShard().getCellCount());
    EXPECT_EQ((loaded.getVertexCoordinates(0) - original).norm(), 0);
  }

  /// @brief Discovery that needs new entities fails while a rank's existing ownership stays valid.
  TEST(Kokkos_MPI_Mesh, ResidentConnectivityPreservesOwnership)
  {
    using namespace Geometry;
    Context::MPI context(*meshEnvironment, *meshWorld);
    const auto rank = static_cast<Index>(meshWorld->rank());
    const auto x = static_cast<Real>(rank) * 2;
    const auto owned = Shard::State::Owned;
    Shard::Builder builder;
    builder.initialize(2, 2);
    IndexArray vertices(3);
    vertices << builder.vertex(3 * rank, Math::SpatialPoint{x, 0}, owned),
      builder.vertex(3 * rank + 1, Math::SpatialPoint{x + 1, 0}, owned),
      builder.vertex(3 * rank + 2, Math::SpatialPoint{x, 1}, owned);
    builder.polytope(2, rank, Polytope::Type::Triangle, vertices, owned);
    auto host =
      Mesh<Context::MPI>::Builder(context).initialize(builder.finalize()).finalize();
    using DistributedContext = Context::Kokkos<Context::MPI, TestBackend>;
    Mesh<DistributedContext> mesh(DistributedContext(context), host);
    using ResidentConnectivity = typename decltype(mesh)::LocalType::ConnectivityType;
    static_assert(
      !std::is_assignable_v<decltype(mesh.getConnectivity()), ResidentConnectivity>);
    static_assert(!std::is_move_assignable_v<typename decltype(mesh)::ConnectivityType>);
    static_assert(
      !std::is_move_constructible_v<typename decltype(mesh)::ConnectivityType>);
    static_assert(
      !std::is_convertible_v<decltype(mesh.getConnectivity()), ResidentConnectivity&>);

    auto& connectivity = mesh.getConnectivity();
    EXPECT_THROW(connectivity.compute(2, 1), Alert::Exception);
    EXPECT_FALSE(connectivity.hasIncidence(2, 1));
    EXPECT_EQ(mesh.getFaceCount(), 0u);
    size_t ownedCells = 0;
    ::Kokkos::parallel_reduce("Count cells after rejected discovery",
      ::Kokkos::RangePolicy<TestExecutionSpace>(
        mesh.getContext().getExecutionSpace(), 0, mesh.getCellCount()),
      OwnedCells{mesh.getView()}, ownedCells);
    EXPECT_EQ(ownedCells, 1u);
    connectivity.restrict(2, 1);
    EXPECT_TRUE(connectivity.getIncidence({2, 1}, 0).empty());

    host.getConnectivity().compute(2, 1);
    EXPECT_THROW(mesh.upload(host), Alert::Exception);
    EXPECT_EQ(mesh.getFaceCount(), 0u);
    host.reconcile(1);
    mesh.upload(host);
    EXPECT_EQ(mesh.getFaceCount(), 3u);
    mesh.getConnectivity().compute(1, 2);
    size_t boundary = 0;
    for (auto face = mesh.getBoundary(); face; ++face)
      ++boundary;
    EXPECT_EQ(boundary, 3u);
    EXPECT_FALSE(mesh.getInterface());
    EXPECT_NEAR(mesh.getPerimeter(), meshWorld->size() * (2 + std::sqrt(Real(2))), 1e-13);
  }

  /// @brief Facet traversal and integration use the global dimension on every rank.
  TEST(Kokkos_MPI_Mesh, FacetsUseDistributedDimension)
  {
    using namespace Geometry;
    using namespace Variational;
    Context::MPI context(*meshEnvironment, *meshWorld);
    const auto rank = static_cast<Index>(meshWorld->rank());
    const auto owned = Shard::State::Owned;
    Shard::Builder builder;
    builder.initialize(2, 2);
    if (rank == 0)
    {
      const Index a = builder.vertex(0, Math::SpatialPoint{0, 0}, owned);
      const Index b = builder.vertex(1, Math::SpatialPoint{1, 0}, owned);
      const Index c = builder.vertex(2, Math::SpatialPoint{0, 1}, owned);
      const Index d = builder.vertex(3, Math::SpatialPoint{1, 1}, owned);
      IndexArray triangle(3), edge(2);
      triangle << a, b, c;
      builder.polytope(2, 0, Polytope::Type::Triangle, triangle, owned);
      triangle << b, d, c;
      builder.polytope(2, 1, Polytope::Type::Triangle, triangle, owned);
      const std::array edges{std::pair{a, b}, std::pair{b, d}, std::pair{d, c},
        std::pair{c, a}, std::pair{b, c}};
      for (Index i = 0; i < edges.size(); ++i)
      {
        edge << edges[i].first, edges[i].second;
        builder.polytope(1, i, Polytope::Type::Segment, edge, owned);
      }
    }
    else if (rank == 1)
    {
      IndexArray segment(2);
      segment << builder.vertex(4, Math::SpatialPoint{2, 0}, owned),
        builder.vertex(5, Math::SpatialPoint{3, 0}, owned);
      builder.polytope(1, 5, Polytope::Type::Segment, segment, owned);
    }
    const auto host =
      Mesh<Context::MPI>::Builder(context).initialize(builder.finalize()).finalize();
    using DistributedContext = Context::Kokkos<Context::MPI, TestBackend>;
    Mesh<DistributedContext> mesh(DistributedContext(context), host);
    EXPECT_EQ(mesh.getGlobalDimension(), 2u);
    EXPECT_EQ(mesh.getDimension(), rank == 0 ? 2u : rank == 1 ? 1u : 0u);
    if (rank == 0)
      mesh.getConnectivity().compute(1, 2);
    size_t localBoundary = 0;
    for (auto face = mesh.getBoundary(); face; ++face)
      ++localBoundary;
    EXPECT_EQ(localBoundary, rank == 0 ? 4u : 0u);
    EXPECT_NEAR(mesh.getPerimeter(), 4, 1e-13);

    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    auto boundary = BoundaryIntegral(u, v);
    if (rank == 0)
    {
      mesh.getConnectivity().setDirty(1, 2);
      EXPECT_THROW(boundary.prepare(false), Alert::Exception);
      mesh.getConnectivity().compute(1, 2);
    }
    // Each selected facet contributes a two-by-two P1 matrix.
    EXPECT_EQ(boundary.assemble().values.extent(0), rank == 0 ? 16u : 0u);
    EXPECT_EQ(InterfaceIntegral(u, v).assemble().values.extent(0), rank == 0 ? 4u : 0u);
    EXPECT_EQ(FaceIntegral(u, v).assemble().values.extent(0), rank == 0 ? 20u : 0u);
  }
}

/// @brief Runs MPI/Kokkos tests without initializing or linking PETSc.
int main(int argc, char** argv)
{
  int status = 0;
  {
    boost::mpi::environment environment(argc, argv);
    boost::mpi::communicator world;
    Rodin::Tests::Unit::meshEnvironment = &environment;
    Rodin::Tests::Unit::meshWorld = &world;
    ::testing::InitGoogleTest(&argc, argv);
    ::Kokkos::ScopeGuard runtime(argc, argv);
    status = RUN_ALL_TESTS();
  }
  Rodin::Tests::Unit::exitTestProcess(status);
}
