/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 *
 * Unit tests for the MPI assembly module.
 *
 * These tests verify that the MPI assembly infrastructure works correctly:
 * - MPIIteration iterates over local mesh elements
 * - Assembly::MPI<IndexMap<...>, DirichletBC<...>> assembles boundary conditions
 *   on a distributed mesh
 * - For 1-rank runs the assembled boundary DOF set is consistent with the
 *   sequential assembly
 *
 * Run with mpirun -n 1/2/3 as registered in CMakeLists.txt.
 */
#include <set>
#include <numeric>

#include <gtest/gtest.h>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>

#include <Rodin/Geometry.h>
#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#include <Rodin/MPI/Context/MPI.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/Assembly.h>
#include <Rodin/MPI/Variational/P1.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

// ---------------------------------------------------------------------------
// Global MPI handles (initialized in main())
// ---------------------------------------------------------------------------
static boost::mpi::environment* g_env   = nullptr;
static boost::mpi::communicator* g_world = nullptr;

namespace
{
  /**
   * @brief Creates a local mesh with all incidences required for sharding
   * and boundary-DOF assembly.
   */
  static Mesh<Context::Local> makeShardableMesh(
      Polytope::Type type,
      std::initializer_list<size_t> shape)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(type, shape);
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);
    mesh.getConnectivity().compute(D, 0);
    mesh.getConnectivity().compute(D, D - 1);
    mesh.getConnectivity().compute(D - 1, D);
    mesh.getConnectivity().compute(D - 1, 0);
    return mesh;
  }

  /**
   * @brief Distributes a uniform-grid mesh from rank 0.
   */
  static Mesh<Context::MPI> distributeFromRoot(
      const Context::MPI& ctx,
      Polytope::Type type,
      std::initializer_list<size_t> shape)
  {
    const auto& comm = ctx.getCommunicator();

    Sharder<Context::MPI> sharder(ctx);
    if (comm.rank() == 0)
    {
      auto localMesh = makeShardableMesh(type, shape);
      BalancedCompactPartitioner partitioner(localMesh);
      partitioner.partition(static_cast<size_t>(comm.size()));
      sharder.shard(partitioner);
      sharder.scatter(0);
    }
    return sharder.gather(0);
  }
}

namespace Rodin::Tests::Unit
{
  // =========================================================================
  // MPIIteration
  // =========================================================================

  /**
   * @brief MPIIteration over a Triangle mesh yields at least one cell.
   */
  TEST(Assembly_MPI_Iteration, TriangleMesh_HasCells)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 4, 4 });

    size_t localCount = 0;
    Assembly::MPIIteration iter(mpiMesh, Geometry::Region::Cells);
    for (auto it = iter.getIterator(); it; ++it)
      ++localCount;

    // At least one cell on every rank for a 4x4 mesh (32 triangles total)
    EXPECT_GT(localCount, 0u);
  }

  /**
   * @brief MPIIteration over a Quadrilateral mesh yields at least one cell.
   */
  TEST(Assembly_MPI_Iteration, QuadrilateralMesh_HasCells)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Quadrilateral, { 4, 4 });

    size_t localCount = 0;
    Assembly::MPIIteration iter(mpiMesh, Geometry::Region::Cells);
    for (auto it = iter.getIterator(); it; ++it)
      ++localCount;

    EXPECT_GT(localCount, 0u);
  }

  /**
   * @brief Global cell count equals the sum of local counts across all ranks.
   */
  TEST(Assembly_MPI_Iteration, GlobalCellCount_MatchesLocal_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 4, 4 });

    // Count only owned cells on this rank (ghost cells must be excluded to
    // avoid double-counting when reducing across ranks).
    const size_t D = mpiMesh.getDimension();
    const auto& shard = mpiMesh.getShard();
    size_t localCount = 0;
    {
      Assembly::MPIIteration iter(mpiMesh, Geometry::Region::Cells);
      for (auto it = iter.getIterator(); it; ++it)
      {
        if (shard.isOwned(D, it->getIndex()))
          ++localCount;
      }
    }

    // Reduce to root
    size_t totalCount = 0;
    boost::mpi::reduce(world, localCount, totalCount, std::plus<size_t>(), 0);

    if (world.rank() == 0)
    {
      // A 4x4 uniform triangle grid has 2*(4-1)*(4-1) = 18 triangles
      EXPECT_EQ(totalCount, 18u);
    }
  }

  // =========================================================================
  // Assembly::MPI<IndexMap<Real>, DirichletBC<...>>
  // =========================================================================

  /**
   * @brief MPI DirichletBC assembly yields at least one boundary DOF entry
   * on a Triangle mesh with a constant boundary condition.
   */
  TEST(Assembly_MPI_DirichletBC, ConstantBC_NonEmpty_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 4, 4 });

    P1<Real, Mesh<Context::MPI>> fes(mpiMesh);
    TrialFunction u(fes);

    DirichletBC dbc(u, RealFunction(1.0));
    dbc.assemble();

    // At least globally some DOFs must be fixed on the boundary
    const size_t localFixed = std::get<IndexMap<Real>>(dbc.getDOFs()).size();
    size_t globalFixed = 0;
    boost::mpi::reduce(world, localFixed, globalFixed, std::plus<size_t>(), 0);

    if (world.rank() == 0)
    {
      EXPECT_GT(globalFixed, 0u);
    }
  }

  /**
   * @brief MPI DirichletBC assembly: boundary DOF count is consistent with
   * the sequential assembly on the same mesh for a 1-rank run.
   */
  TEST(Assembly_MPI_DirichletBC, SingleRank_MatchesSequentialCount_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() != 1)
      GTEST_SKIP() << "This test is designed for exactly 1 MPI rank.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 4, 4 });

    P1<Real, Mesh<Context::MPI>> mpiFes(mpiMesh);
    TrialFunction uMPI(mpiFes);
    DirichletBC dbcMPI(uMPI, RealFunction(1.0));
    dbcMPI.assemble();
    const size_t mpiFixed = std::get<IndexMap<Real>>(dbcMPI.getDOFs()).size();

    // Sequential reference
    auto localMesh = makeShardableMesh(Polytope::Type::Triangle, { 4, 4 });
    P1 seqFes(localMesh);
    TrialFunction uSeq(seqFes);
    DirichletBC dbcSeq(uSeq, RealFunction(1.0));
    dbcSeq.assemble();
    const size_t seqFixed = std::get<IndexMap<Real>>(dbcSeq.getDOFs()).size();

    EXPECT_EQ(mpiFixed, seqFixed);
  }

  /**
   * @brief MPI DirichletBC assembly: boundary DOF values are correct
   * for a constant boundary condition.
   */
  TEST(Assembly_MPI_DirichletBC, ConstantBC_Values_AllEqualPrescribed_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    const Real gValue = 3.14;

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 4, 4 });

    P1<Real, Mesh<Context::MPI>> fes(mpiMesh);
    TrialFunction u(fes);

    DirichletBC dbc(u, RealFunction(gValue));
    dbc.assemble();

    for (const auto& [local, value]
           : std::get<IndexMap<Real>>(dbc.getDOFs()))
      EXPECT_NEAR(value, gValue, 1e-12);
  }

  // =========================================================================
  // All-geometry MPIIteration tests — Segment (1D), Tetrahedron/Hexahedron (3D)
  // =========================================================================

  /**
   * @brief MPIIteration over a 1D Segment mesh yields at least one cell.
   */
  TEST(Assembly_MPI_Iteration, SegmentMesh_HasCells)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Segment, { 10 });

    size_t localCount = 0;
    Assembly::MPIIteration iter(mpiMesh, Geometry::Region::Cells);
    for (auto it = iter.getIterator(); it; ++it)
      ++localCount;

    // At least one cell on every rank for a 10-segment mesh
    EXPECT_GT(localCount, 0u);
  }

  /**
   * @brief Global cell count for a 1D Segment mesh equals the expected value.
   */
  TEST(Assembly_MPI_Iteration, GlobalCellCount_MatchesLocal_Segment)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Segment, { 10 });

    // Count only owned cells on this rank (ghost cells must be excluded to
    // avoid double-counting when reducing across ranks).
    const size_t D = mpiMesh.getDimension();
    const auto& shard = mpiMesh.getShard();
    size_t localCount = 0;
    {
      Assembly::MPIIteration iter(mpiMesh, Geometry::Region::Cells);
      for (auto it = iter.getIterator(); it; ++it)
      {
        if (shard.isOwned(D, it->getIndex()))
          ++localCount;
      }
    }

    size_t totalCount = 0;
    boost::mpi::reduce(world, localCount, totalCount, std::plus<size_t>(), 0);

    if (world.rank() == 0)
    {
      // A UniformGrid Segment {10} has 9 cells
      EXPECT_EQ(totalCount, 9u);
    }
  }

  /**
   * @brief MPIIteration over a Tetrahedron mesh yields at least one cell.
   */
  TEST(Assembly_MPI_Iteration, TetrahedronMesh_HasCells)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Tetrahedron, { 3, 3, 3 });

    size_t localCount = 0;
    Assembly::MPIIteration iter(mpiMesh, Geometry::Region::Cells);
    for (auto it = iter.getIterator(); it; ++it)
      ++localCount;

    EXPECT_GT(localCount, 0u);
  }

  /**
   * @brief MPIIteration over a Hexahedron mesh yields at least one cell.
   */
  TEST(Assembly_MPI_Iteration, HexahedronMesh_HasCells)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Hexahedron, { 3, 3, 3 });

    size_t localCount = 0;
    Assembly::MPIIteration iter(mpiMesh, Geometry::Region::Cells);
    for (auto it = iter.getIterator(); it; ++it)
      ++localCount;

    EXPECT_GT(localCount, 0u);
  }

  // =========================================================================
  // All-geometry DirichletBC MPI tests
  // =========================================================================

  /**
   * @brief MPI DirichletBC assembly on Segment (1D) mesh: at least one
   * boundary DOF is found and its value matches the prescribed constant.
   */
  TEST(Assembly_MPI_DirichletBC, ConstantBC_Segment)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    const Real gValue = 2.71;

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Segment, { 10 });

    P1<Real, Mesh<Context::MPI>> fes(mpiMesh);
    TrialFunction u(fes);

    DirichletBC dbc(u, RealFunction(gValue));
    dbc.assemble();

    const size_t localFixed = std::get<IndexMap<Real>>(dbc.getDOFs()).size();
    size_t globalFixed = 0;
    boost::mpi::reduce(world, localFixed, globalFixed, std::plus<size_t>(), 0);

    if (world.rank() == 0)
    {
      EXPECT_GT(globalFixed, 0u);
    }

    for (const auto& [local, value]
           : std::get<IndexMap<Real>>(dbc.getDOFs()))
      EXPECT_NEAR(value, gValue, 1e-12);
  }

  /**
   * @brief MPI DirichletBC assembly on Tetrahedron (3D) mesh: at least one
   * boundary DOF found globally.
   */
  TEST(Assembly_MPI_DirichletBC, ConstantBC_NonEmpty_Tetrahedron)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Tetrahedron, { 3, 3, 3 });

    P1<Real, Mesh<Context::MPI>> fes(mpiMesh);
    TrialFunction u(fes);

    DirichletBC dbc(u, RealFunction(1.0));
    dbc.assemble();

    const size_t localFixed = std::get<IndexMap<Real>>(dbc.getDOFs()).size();
    size_t globalFixed = 0;
    boost::mpi::reduce(world, localFixed, globalFixed, std::plus<size_t>(), 0);

    if (world.rank() == 0)
    {
      EXPECT_GT(globalFixed, 0u);
    }
  }

  /**
   * @brief MPI DirichletBC assembly on Hexahedron (3D) mesh: at least one
   * boundary DOF found globally.
   */
  TEST(Assembly_MPI_DirichletBC, ConstantBC_NonEmpty_Hexahedron)
  {
    const auto& world = *g_world;
    if (world.size() > 3)
      GTEST_SKIP() << "Test designed for at most 3 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mpiMesh = distributeFromRoot(ctx, Polytope::Type::Hexahedron, { 3, 3, 3 });

    P1<Real, Mesh<Context::MPI>> fes(mpiMesh);
    TrialFunction u(fes);

    DirichletBC dbc(u, RealFunction(1.0));
    dbc.assemble();

    const size_t localFixed = std::get<IndexMap<Real>>(dbc.getDOFs()).size();
    size_t globalFixed = 0;
    boost::mpi::reduce(world, localFixed, globalFixed, std::plus<size_t>(), 0);

    if (world.rank() == 0)
    {
      EXPECT_GT(globalFixed, 0u);
    }
  }
}

// ---------------------------------------------------------------------------
// main() — initializes MPI environment used by all tests.
// ---------------------------------------------------------------------------
int main(int argc, char** argv)
{
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  g_env   = &env;
  g_world = &world;

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
