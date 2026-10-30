/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/// @file @brief Single-field boundary refreshes with empty distributed shards.
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#include <Rodin/MPI/Context/MPI.h>
#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/MPI/Variational/P1.h>
#include "SingleFieldTargetAssembly.h"
#include <Rodin/MPI/Variational/H1.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  boost::mpi::environment* environment = nullptr;
  boost::mpi::communicator* world = nullptr;

  /// @brief Retains all communicator ranks while leaving rank zero empty in parallel.
  class EmptyRootPartitioner final : public Partitioner
  {
    public:
      /// @brief Assigns single-cell partitions among the available non-root ranks.
      EmptyRootPartitioner(const LocalMesh& mesh, size_t count)
        : m_base(mesh),
          m_count(count)
      {
        m_base.partition(1);
      }

      /// @brief Returns the source interval mesh.
      const LocalMesh& getMesh() const override
      {
        return m_base.getMesh();
      }

      /// @brief Updates the cell grouping before assigning communicator ranks.
      void partition(size_t maximumSize, size_t dimension) override
      {
        m_base.partition(maximumSize, dimension);
      }

      /// @brief Maps each cell group to an available non-root rank.
      size_t getPartition(Index index) const override
      {
        return m_count == 1 ? 0 : 1 + m_base.getPartition(index) % (m_count - 1);
      }

      /// @brief Includes ranks with no cells in the shard layout.
      size_t getCount() const override
      {
        return m_count;
      }

    private:
      BalancedCompactPartitioner m_base;
      size_t m_count;
  };
}

/// @brief A distributed single-side update preserves the full constrained linear system.
template <bool Quadratic>
void checkDistributedBoundaryRefresh()
{
  Context::MPI context(*environment, *world);
  Sharder<Context::MPI> sharder(context);
  if (world->rank() == 0)
  {
    auto local = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    local.scale(Real(0.5));
    local.getConnectivity().compute(0, 1);
    local.getConnectivity().compute(1, 1);
    for (auto vertex = local.getVertex(); vertex; ++vertex)
    {
      const Real x = local.getVertexCoordinates(vertex->getIndex()).x();
      if (x == 0)
        local.setAttribute(vertex.key(), 11);
      if (x == 1)
        local.setAttribute(vertex.key(), 12);
    }
    EmptyRootPartitioner partitioner(local, static_cast<size_t>(world->size()));
    sharder.shard(partitioner);
    sharder.scatter(0);
  }
  auto mesh = sharder.gather(0);
  if (world->size() > 1 && world->rank() == 0)
  {
    EXPECT_EQ(mesh.getShard().getCellCount(), 0u);
  }
  Rodin::Tests::Unit::checkSingleFieldTargets<Assembly::MPI, Quadratic, 0>(
    mesh, MPI_COMM_WORLD, true);
}

/// @brief P1 full-to-partial transitions retain boundary lifting on every rank.
TEST(PETSc_MPIBoundaryRefresh, P1RefreshAfterFull)
{
  checkDistributedBoundaryRefresh<false>();
}

/// @brief H1 full-to-partial transitions retain boundary lifting on every rank.
TEST(PETSc_MPIBoundaryRefresh, H1RefreshAfterFull)
{
  checkDistributedBoundaryRefresh<true>();
}

int main(int argc, char** argv)
{
  boost::mpi::environment mpiEnvironment(argc, argv);
  boost::mpi::communicator communicator;
  environment = &mpiEnvironment;
  world = &communicator;
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  ::testing::InitGoogleTest(&argc, argv);
  const int localResult = RUN_ALL_TESTS();
  int result = 0;
  MPI_Allreduce(&localResult, &result, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  PetscFinalize();
  return result;
}
