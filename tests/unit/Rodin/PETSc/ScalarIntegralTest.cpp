/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "ScalarIntegral.h"
#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#ifdef RODIN_USE_MPI
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/MPI/Variational/H1.h>
#include <Rodin/MPI/Variational/P0.h>
#include <Rodin/MPI/Variational/P0g.h>
#include <Rodin/MPI/Variational/P1.h>
#endif

using namespace Rodin;
using namespace Rodin::Geometry;

namespace
{
  constexpr Polytope::Type geometries[] = {Polytope::Type::Segment,
    Polytope::Type::Triangle, Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
    Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge};

  LocalMesh unitMesh(Polytope::Type geometry)
  {
    const size_t dim = Polytope::Traits(geometry).getDimension();
    auto mesh = LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dim, 2));
    for (size_t d = 0; d <= dim; ++d)
    {
      for (size_t dp = 0; dp <= dim; ++dp)
        mesh.getConnectivity().compute(d, dp);
    }
    return mesh;
  }

  TEST(PETSc_ScalarIntegral, Local)
  {
    for (auto geometry : geometries)
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      Tests::Unit::ScalarIntegral::check(unitMesh(geometry));
    }
  }

#ifdef RODIN_USE_MPI
  boost::mpi::environment* environment = nullptr;
  boost::mpi::communicator* world = nullptr;

  TEST(PETSc_ScalarIntegral, Distributed)
  {
    Context::MPI context(*environment, *world);
    for (auto geometry : geometries)
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      Sharder<Context::MPI> sharder(context);
      if (world->rank() == 0)
      {
        auto local = unitMesh(geometry);
        BalancedCompactPartitioner partitioner(local);
        partitioner.partition(world->size());
        sharder.shard(partitioner);
      }
      sharder.scatter(0);
      auto mesh = sharder.gather(0);
      // Incidence discovery and distributed ownership are separate operations.
      for (size_t d = 0; d <= mesh.getDimension(); ++d)
      {
        for (size_t dp = 0; dp <= mesh.getDimension(); ++dp)
          mesh.getConnectivity().compute(d, dp);
      }
      for (size_t d = 1; d < mesh.getDimension(); ++d)
        mesh.reconcile(d);
      Tests::Unit::ScalarIntegral::check(mesh);
    }
  }
#endif
}

int main(int argc, char** argv)
{
#ifdef RODIN_USE_MPI
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator comm;
  environment = &env;
  world = &comm;
#endif
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
