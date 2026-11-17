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
    const Geometry::Mesh<DistributedContext> mesh(DistributedContext(context), host);
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
    EXPECT_NE(copy.getNumberingRevision(), scalar.getNumberingRevision());
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
  }
}

/// @brief Runs MPI/Kokkos tests without initializing or linking PETSc.
int main(int argc, char** argv)
{
  boost::mpi::environment environment(argc, argv);
  boost::mpi::communicator world;
  Rodin::Tests::Unit::meshEnvironment = &environment;
  Rodin::Tests::Unit::meshWorld = &world;
  ::testing::InitGoogleTest(&argc, argv);
  ::Kokkos::ScopeGuard runtime(argc, argv);
  return RUN_ALL_TESTS();
}
