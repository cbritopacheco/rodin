/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

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

    using DistributedContext = Context::Kokkos<Context::MPI, TestBackend>;
    Geometry::Mesh<DistributedContext> mesh(DistributedContext(context), host);
    EXPECT_NE(&mesh.getHost(), &host);
    Variational::P1 scalar(mesh);
    Variational::H1 quadratic(std::integral_constant<size_t, 2>{}, mesh);
    Variational::P0g constant(mesh);
    EXPECT_EQ(scalar.getSize(), 25u);
    EXPECT_EQ(quadratic.getSize(), 81u);
    EXPECT_EQ(constant.getSize(), 1u);
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
      copy.upload(host);
      EXPECT_EQ(copy.getAttribute(2, 0), attribute);
      EXPECT_EQ((copy.getVertexCoordinates(0) - original).norm(), 0);
      typename decltype(copy)::ConnectivityType emptyTopology(TestContext{});
      copy.getConnectivity() = std::move(emptyTopology);
      EXPECT_THROW(copy.getView(), Alert::Exception);
      EXPECT_THROW(copy.getHost(), Alert::Exception);
      copy.upload(host);
      EXPECT_NO_THROW(copy.getView());
    }
    meshWorld->barrier();
    host.save(file, IO::FileFormat::HDF5);
    Geometry::Mesh<DistributedContext> loaded{DistributedContext(context)};
    loaded.load(file, IO::FileFormat::HDF5);
    boost::filesystem::remove(file);
    EXPECT_EQ(loaded.getCellCount(), host.getShard().getCellCount());
    EXPECT_EQ((loaded.getVertexCoordinates(0) - original).norm(), 0);
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
