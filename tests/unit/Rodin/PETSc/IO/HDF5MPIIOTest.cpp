/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <petsc.h>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include "Rodin/Geometry.h"
#include "Rodin/Geometry/BalancedCompactPartitioner.h"
#include "Rodin/Variational.h"
#include "Rodin/MPI/Context/MPI.h"
#include "Rodin/MPI/Geometry/Sharder.h"
#include "Rodin/MPI/Variational.h"
#include "Rodin/PETSc.h"
#include "Rodin/PETSc/IO.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

/// @brief Restoring owned coefficients also refreshes every local ghost value.
TEST(PETScHDF5MPI, CheckpointRefreshesGhosts)
{
  boost::mpi::environment environment;
  boost::mpi::communicator communicator;
  Context::MPI context(environment, communicator);
  Sharder<Context::MPI> sharder(context);
  if (communicator.rank() == 0)
  {
    auto local = LocalMesh::UniformGrid(Polytope::Type::Segment, {7});
    local.getConnectivity().compute(1, 1);
    local.getConnectivity().compute(0, 1);
    BalancedCompactPartitioner partitioner(local);
    partitioner.partition(communicator.size());
    sharder.shard(partitioner);
    sharder.scatter(0);
  }
  auto mesh = sharder.gather(0);
  P1<PetscScalar, Mesh<Context::MPI>> space(mesh);
  GridFunction<decltype(space), ::Vec> field(space);
  field = RealFunction(7.);
  field.sync();
  const auto path = boost::filesystem::temp_directory_path() /
    boost::filesystem::unique_path("rodin-hdf5-ghost-%%%%-%%%%.h5");
  field.save(path, IO::FileFormat::HDF5);
  field = RealFunction(0.);
  field.sync();
  field.load(path, IO::FileFormat::HDF5);
  field.sync();
  Index begin = 0, end = 0;
  space.getOwnershipRange(begin, end);
  size_t ghosts = 0;
  const auto& read = field;
  for (size_t local = 0; local < space.getShard().getSize(); ++local)
  {
    const Index global = space.getGlobalIndex(local);
    ghosts += global < begin || global >= end;
    EXPECT_LT(PetscAbsScalar(read[global] - PetscScalar(7)), 1e-12);
  }
  read.flush();
  if (communicator.size() > 1)
  {
    EXPECT_GT(ghosts, 0u);
  }
  boost::filesystem::remove(path);
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
