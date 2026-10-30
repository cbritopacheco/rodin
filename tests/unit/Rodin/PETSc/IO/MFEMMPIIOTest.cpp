/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Analytic and coefficient oracles for PETSc MFEM persistence. */
#include <gtest/gtest.h>
#include <petsc.h>
#include <memory>
#include <sstream>
#include <type_traits>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/PETSc.h>
#include <Rodin/PETSc/IO.h>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#include <Rodin/MPI/Context/MPI.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/MPI/Variational.h>
#include "MFEMTest.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::IO;

static boost::mpi::environment* g_environment = nullptr;
static boost::mpi::communicator* g_communicator = nullptr;

namespace
{
  template <size_t Family, size_t Rank>
  void checkExport(Polytope::Type geometry)
  {
    Context::MPI context(*g_environment, *g_communicator);
    Sharder<Context::MPI> sharder(context);
    if (g_communicator->rank() == 0)
    {
      const size_t dimension = Polytope::Traits(geometry).getDimension();
      auto local =
        LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 2));
      for (size_t from = 0; from <= dimension; ++from)
      {
        for (size_t to = 0; to <= dimension; ++to)
        {
          if (from != to || from == dimension)
            local.getConnectivity().compute(from, to);
        }
      }
      BalancedCompactPartitioner partitioner(local);
      partitioner.partition(static_cast<size_t>(g_communicator->size()));
      sharder.shard(partitioner);
      sharder.scatter(0);
    }
    auto mesh = sharder.gather(0);
    MFEMTest::checkExport<Family, Rank>(mesh);
  }

  class PETScMFEMMPI : public ::testing::TestWithParam<Polytope::Type>
  {};

  TEST_P(PETScMFEMMPI, CubicScalar)
  {
    checkExport<3, 0>(GetParam());
  }

  TEST_P(PETScMFEMMPI, CubicVector)
  {
    checkExport<3, 1>(GetParam());
  }

  TEST_P(PETScMFEMMPI, CubicMatrix)
  {
    checkExport<3, 2>(GetParam());
  }
  TEST_P(PETScMFEMMPI, LinearScalar)
  {
    checkExport<1, 0>(GetParam());
  }
  TEST_P(PETScMFEMMPI, LinearVector)
  {
    checkExport<1, 1>(GetParam());
  }
  TEST_P(PETScMFEMMPI, LinearMatrix)
  {
    checkExport<1, 2>(GetParam());
  }
  TEST_P(PETScMFEMMPI, ConstantScalar)
  {
    checkExport<0, 0>(GetParam());
  }
  TEST_P(PETScMFEMMPI, ConstantVector)
  {
    checkExport<0, 1>(GetParam());
  }
  TEST_P(PETScMFEMMPI, ConstantMatrix)
  {
    checkExport<0, 2>(GetParam());
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, PETScMFEMMPI,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge));
}

int main(int argc, char** argv)
{
  boost::mpi::environment environment(argc, argv);
  boost::mpi::communicator communicator;
  g_environment = &environment;
  g_communicator = &communicator;
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
