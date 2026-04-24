/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file PoissonTest.cpp
 * @brief Distributed PETSc tests for the Poisson problem on `Context::MPI`
 *        meshes.
 *
 * These tests distribute a mesh across MPI ranks, assemble the Poisson
 * variational problem with the PETSc backend, solve with `PETSc::Solver::CG`
 * (and `PETSc::Solver::GMRES`), and verify the solution against a known
 * manufactured solution.
 *
 * The L² error is computed per-rank and then reduced with
 * `boost::mpi::all_reduce` so that every rank can independently assert
 * correctness.
 *
 * ### Workflow
 *
 * 1. Rank 0 builds a uniform-grid mesh, partitions it with
 *    `BalancedCompactPartitioner`, and calls `sharder.scatter()`.
 * 2. Each rank saves its shard to HDF5.
 * 3. Each rank rereads its shard from HDF5 and rebuilds the shard metadata.
 * 4. Incidence data required for assembly and Dirichlet-BC enforcement is
 *    computed and the mesh is reconciled.
 * 5. `PETSc::Variational::TrialFunction / TestFunction` wrap the distributed
 *    `P1` space.
 * 6. `Problem<PETSc::Math::LinearSystem, U, V>` is assembled with the MPI
 *    assembly backend and solved with a PETSc KSP solver.
 * 7. The local per-rank L² error is summed globally and checked against the
 *    manufactured tolerance.
 */

#include <cmath>
#include <numeric>
#include <string>

#include <gtest/gtest.h>

#include <hdf5.h>
#include <petsc.h>
#include <boost/filesystem.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>

#include <Rodin/Configure.h>
#include <Rodin/Types.h>
#include <Rodin/Geometry.h>
#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#include <Rodin/Variational.h>
#include <Rodin/IO.h>
#include <Rodin/MPI.h>
#include <Rodin/MPI/Context/MPI.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/IO.h>
#include <Rodin/MPI/Variational/P1.h>
#include <Rodin/PETSc.h>

#if defined(_WIN32)
#include <process.h>
#define RODIN_GETPID _getpid
#else
#include <unistd.h>
#define RODIN_GETPID getpid
#endif

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

// ---------------------------------------------------------------------------
// Global MPI handles (initialized in main())
// ---------------------------------------------------------------------------
static boost::mpi::environment*  g_env   = nullptr;
static boost::mpi::communicator* g_world = nullptr;

namespace
{
  // -------------------------------------------------------------------------
  // Helper: build a mesh on rank 0 with full incidence data, then distribute.
  // -------------------------------------------------------------------------

  /**
   * @brief Builds a local mesh with all incidences needed for sharding and
   *        Dirichlet-BC assembly.
   */
  static Mesh<Context::Local> makeShardableMesh(
      Polytope::Type type,
      std::initializer_list<size_t> shape)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(type, shape);
    mesh.scale(1.0 / static_cast<Real>(*shape.begin() - 1));
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);
    mesh.getConnectivity().compute(D, 0);
    mesh.getConnectivity().compute(D, D - 1);
    mesh.getConnectivity().compute(D - 1, D);
    mesh.getConnectivity().compute(D - 1, 0);
    return mesh;
  }

  static Shard rebuildShardFromHDF5(
      const Mesh<Context::Local>& baseMesh,
      const boost::filesystem::path& filepath,
      size_t D)
  {
    using namespace Rodin::IO;

    const hid_t h5 = H5Fopen(filepath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (h5 < 0)
    {
      ADD_FAILURE() << "Cannot open HDF5 shard file " << filepath;
      return Shard{};
    }

    auto vtxState = HDF5::readVectorDataset<HDF5::U8>(h5, HDF5::shardStatePath(0));
    auto vtxLeft = HDF5::readVectorDataset<HDF5::U64>(h5, HDF5::shardPolytopeMapLeftPath(0));
    auto cellState = HDF5::readVectorDataset<HDF5::U8>(h5, HDF5::shardStatePath(D));
    auto cellLeft = HDF5::readVectorDataset<HDF5::U64>(h5, HDF5::shardPolytopeMapLeftPath(D));

    const std::string ownerPath = HDF5::shardOwnerGroupPath(D);
    auto ownerKeys = HDF5::readVectorDataset<HDF5::U64>(h5, ownerPath + "/Keys");
    auto ownerVals = HDF5::readVectorDataset<HDF5::U64>(h5, ownerPath + "/Values");

    const std::string haloPath = HDF5::shardHaloGroupPath(D);
    auto haloKeys = HDF5::readVectorDataset<HDF5::U64>(h5, haloPath + "/Keys");
    auto haloOffsets = HDF5::readVectorDataset<HDF5::U64>(h5, haloPath + "/Offsets");
    auto haloIndices = HDF5::readVectorDataset<HDF5::U64>(h5, haloPath + "/Indices");

    H5Fclose(h5);

    const auto stateOf = [](HDF5::U8 flag) -> Shard::State
    {
      if (flag == 1)
        return Shard::State::Owned;
      if (flag == 2)
        return Shard::State::Ghost;
      return Shard::State::Shared;
    };

    Shard::Builder builder;
    builder.initialize(baseMesh);

    for (Index vi = 0; vi < static_cast<Index>(baseMesh.getVertexCount()); vi++)
      builder.include({ 0, vi }, stateOf(vtxState[vi]));

    for (Index ci = 0; ci < static_cast<Index>(baseMesh.getCellCount()); ci++)
      builder.include({ D, ci }, stateOf(cellState[ci]));

    for (size_t i = 0; i < ownerKeys.size(); i++)
    {
      builder.setOwner(
          D,
          static_cast<Index>(ownerKeys[i]),
          static_cast<Index>(ownerVals[i]));
    }

    for (size_t i = 0; i < haloKeys.size(); i++)
    {
      const Index key = static_cast<Index>(haloKeys[i]);
      const size_t from = static_cast<size_t>(haloOffsets[i]);
      const size_t to = static_cast<size_t>(haloOffsets[i + 1]);
      for (size_t j = from; j < to; j++)
        builder.halo(D, key, static_cast<Index>(haloIndices[j]));
    }

    Shard shard = builder.finalize();

    {
      auto& vmap = shard.getPolytopeMap(0);
      vmap.right.clear();
      for (Index vi = 0; vi < static_cast<Index>(vmap.left.size()); vi++)
      {
        vmap.left[vi] = static_cast<Index>(vtxLeft[vi]);
        vmap.right[vtxLeft[vi]] = vi;
      }
    }

    {
      auto& cmap = shard.getPolytopeMap(D);
      cmap.right.clear();
      for (Index ci = 0; ci < static_cast<Index>(cmap.left.size()); ci++)
      {
        cmap.left[ci] = static_cast<Index>(cellLeft[ci]);
        cmap.right[cellLeft[ci]] = ci;
      }
    }

    return shard;
  }

  static Mesh<Context::MPI> reloadShardViaHDF5(
      const Context::MPI& ctx,
      Mesh<Context::MPI>&& mesh)
  {
    const auto& comm = ctx.getCommunicator();
    const size_t D = mesh.getDimension();

    const boost::filesystem::path rankFile =
        boost::filesystem::temp_directory_path()
        / ("rodin_petsc_mpi_poisson_" + std::to_string(RODIN_GETPID())
           + "_r" + std::to_string(comm.rank()) + ".h5");

    {
      IO::MeshPrinter<IO::FileFormat::HDF5, Context::MPI> printer(mesh);
      printer.print(rankFile);
    }
    comm.barrier();

    Mesh<Context::Local> baseMesh;
    baseMesh.load(rankFile, IO::FileFormat::HDF5);

    EXPECT_EQ(baseMesh.getCellCount(), mesh.getShard().getCellCount())
      << "Rank " << comm.rank() << ": cell count mismatch after HDF5 reload.";
    EXPECT_EQ(baseMesh.getVertexCount(), mesh.getShard().getVertexCount())
      << "Rank " << comm.rank() << ": vertex count mismatch after HDF5 reload.";

    Shard shard = rebuildShardFromHDF5(baseMesh, rankFile, D);

    Mesh<Context::MPI> reloaded =
        Mesh<Context::MPI>::Builder(ctx)
          .initialize(std::move(shard))
          .finalize();

    reloaded.getConnectivity().compute(D, D);
    reloaded.getConnectivity().compute(D, 0);
    reloaded.getConnectivity().compute(D, D - 1);
    reloaded.getConnectivity().compute(D - 1, D);
    reloaded.getConnectivity().compute(D - 1, 0);
    reloaded.getConnectivity().compute(1, 2);
    reloaded.reconcile(1);

    comm.barrier();
    boost::filesystem::remove(rankFile);
    return reloaded;
  }

  /**
   * @brief Distributes a uniform 2D grid from rank 0 to all ranks via
   *        `BalancedCompactPartitioner` + `Sharder<Context::MPI>`.
   *
   * After this function returns, every rank holds a shard that has been
   * persisted through HDF5, reread, and reconciled for FE assembly.
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
    auto mesh = sharder.gather(0);
    return reloadShardViaHDF5(ctx, std::move(mesh));
  }
}

namespace Rodin::Tests::Manufactured::PETSc::MPI
{
  // Alias to resolve the enclosing 'PETSc' namespace scope to Rodin::PETSc.
  namespace PETSc = ::Rodin::PETSc;

  // =========================================================================
  // Distributed Poisson tests
  // =========================================================================

  /**
   * @brief Distributes a 2D Poisson problem across MPI ranks and solves it
   *        with PETSc CG.
   *
   * Manufactured solution: @f$ u(x,y) = \sin(\pi x)\sin(\pi y) @f$.
   *
   * Right-hand side: @f$ f = 2\pi^2 \sin(\pi x)\sin(\pi y) @f$.
   *
   * The global @f$ L^2 @f$ error is obtained by summing local integrals
   * across all ranks.  Expected to stay below @ref RODIN_FUZZY_CONSTANT on a
   * 16×16 mesh.
   */
  TEST(PETSc_MPI_Poisson, CG_SimpleSine_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const auto pi = Math::Constants::pi();

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 16, 16 });

    P1<Real, Mesh<Context::MPI>> vh(mesh);

    auto f = 2 * pi * pi * sin(pi * F::x) * sin(pi * F::y);

    PETSc::Variational::TrialFunction u(vh);
    PETSc::Variational::TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());

    PETSc::Solver::CG solver(poisson);
    solver.solve();

    auto solution = sin(pi * F::x) * sin(pi * F::y);

    P1<Real, Mesh<Context::MPI>> sh(mesh);
    GridFunction<P1<Real, Mesh<Context::MPI>>, ::Vec> diff(sh);
    diff = Pow(u.getSolution() - solution, 2);

    const Real localError = Integral(diff).compute();

    Real globalError = 0;
    boost::mpi::all_reduce(world, localError, globalError, std::plus<Real>());

    EXPECT_NEAR(globalError, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @brief Same test as CG_SimpleSine_Triangle but on a Quadrilateral mesh.
   */
  TEST(PETSc_MPI_Poisson, CG_SimpleSine_Quadrilateral)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const auto pi = Math::Constants::pi();

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Quadrilateral, { 16, 16 });

    P1<Real, Mesh<Context::MPI>> vh(mesh);

    auto f = 2 * pi * pi * sin(pi * F::x) * sin(pi * F::y);

    PETSc::Variational::TrialFunction u(vh);
    PETSc::Variational::TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());

    PETSc::Solver::CG solver(poisson);
    solver.solve();

    auto solution = sin(pi * F::x) * sin(pi * F::y);

    P1<Real, Mesh<Context::MPI>> sh(mesh);
    GridFunction<P1<Real, Mesh<Context::MPI>>, ::Vec> diff(sh);
    diff = Pow(u.getSolution() - solution, 2);

    const Real localError = Integral(diff).compute();

    Real globalError = 0;
    boost::mpi::all_reduce(world, localError, globalError, std::plus<Real>());

    EXPECT_NEAR(globalError, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @brief Single-rank consistency check: distributed solution must match the
   *        sequential PETSc solution on the same mesh.
   *
   * When run with exactly one MPI rank, the distributed assembly should
   * produce the same result as the sequential (`Context::Local`) solve.
   * This test verifies that the distributed path does not introduce extra
   * error compared to the known manufactured-solution tolerance.
   */
  TEST(PETSc_MPI_Poisson, SingleRank_MatchesSequential_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() != 1)
      GTEST_SKIP() << "This test is designed for exactly 1 MPI rank.";

    const auto pi = Math::Constants::pi();

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 16, 16 });

    P1<Real, Mesh<Context::MPI>> vh(mesh);

    auto f = 2 * pi * pi * sin(pi * F::x) * sin(pi * F::y);

    PETSc::Variational::TrialFunction u(vh);
    PETSc::Variational::TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());

    PETSc::Solver::CG solver(poisson);
    solver.solve();

    auto solution = sin(pi * F::x) * sin(pi * F::y);

    P1<Real, Mesh<Context::MPI>> sh(mesh);
    GridFunction<P1<Real, Mesh<Context::MPI>>, ::Vec> diff(sh);
    diff = Pow(u.getSolution() - solution, 2);
    const Real localError = Integral(diff).compute();

    // On 1 rank the local error equals the global error.
    EXPECT_NEAR(localError, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @brief Distributed Poisson solved with PETSc GMRES.
   *
   * Same problem as CG_SimpleSine_Triangle but using the GMRES solver to
   * verify that the non-symmetric KSP path also converges correctly.
   */
  TEST(PETSc_MPI_Poisson, GMRES_SimpleSine_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const auto pi = Math::Constants::pi();

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 16, 16 });

    P1<Real, Mesh<Context::MPI>> vh(mesh);

    auto f = 2 * pi * pi * sin(pi * F::x) * sin(pi * F::y);

    PETSc::Variational::TrialFunction u(vh);
    PETSc::Variational::TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());

    PETSc::Solver::GMRES solver(poisson);
    solver.solve();

    auto solution = sin(pi * F::x) * sin(pi * F::y);

    P1<Real, Mesh<Context::MPI>> sh(mesh);
    GridFunction<P1<Real, Mesh<Context::MPI>>, ::Vec> diff(sh);
    diff = Pow(u.getSolution() - solution, 2);

    const Real localError = Integral(diff).compute();

    Real globalError = 0;
    boost::mpi::all_reduce(world, localError, globalError, std::plus<Real>());

    EXPECT_NEAR(globalError, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @brief P1-exact distributed Poisson: zero FE truncation error.
   *
   * Manufactured solution: @f$ u(x,y) = x + 2y + 1 @f$.
   *
   * Because the manufactured solution lies exactly in the P1 space, the
   * global L² error must be roundoff-zero.
   */
  TEST(PETSc_MPI_Poisson, CG_P1Exact_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 16, 16 });

    P1<Real, Mesh<Context::MPI>> vh(mesh);

    auto solution = F::x + 2 * F::y + 1;
    auto f = Zero();

    PETSc::Variational::TrialFunction u(vh);
    PETSc::Variational::TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, solution);

    PETSc::Solver::CG solver(poisson);
    solver.solve();

    P1<Real, Mesh<Context::MPI>> sh(mesh);
    GridFunction<P1<Real, Mesh<Context::MPI>>, ::Vec> diff(sh);
    diff = Pow(u.getSolution() - solution, 2);

    const Real localError = Integral(diff).compute();

    Real globalError = 0;
    boost::mpi::all_reduce(world, localError, globalError, std::plus<Real>());

    EXPECT_NEAR(globalError, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @brief Polynomial-manufactured distributed Poisson.
   *
   * Manufactured solution: @f$ u(x,y) = x(1-x)y(1-y) @f$.
   *
   * Right-hand side: @f$ f = 2y(1-y) + 2x(1-x) @f$.
   */
  TEST(PETSc_MPI_Poisson, CG_Polynomial_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 16, 16 });

    P1<Real, Mesh<Context::MPI>> vh(mesh);

    auto f = 2 * F::y * (1 - F::y) + 2 * F::x * (1 - F::x);

    PETSc::Variational::TrialFunction u(vh);
    PETSc::Variational::TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());

    PETSc::Solver::CG solver(poisson);
    solver.solve();

    auto solution = F::x * (1 - F::x) * F::y * (1 - F::y);

    P1<Real, Mesh<Context::MPI>> sh(mesh);
    GridFunction<P1<Real, Mesh<Context::MPI>>, ::Vec> diff(sh);
    diff = Pow(u.getSolution() - solution, 2);

    const Real localError = Integral(diff).compute();

    Real globalError = 0;
    boost::mpi::all_reduce(world, localError, globalError, std::plus<Real>());

    EXPECT_NEAR(globalError, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @brief Nonhomogeneous Dirichlet boundary condition, distributed.
   *
   * Manufactured solution: @f$ u(x,y) = \cos(\pi x)\cos(\pi y) @f$.
   *
   * Right-hand side: @f$ f = 2\pi^2\cos(\pi x)\cos(\pi y) @f$.
   */
  TEST(PETSc_MPI_Poisson, CG_NonhomogeneousDirichlet_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const auto pi = Math::Constants::pi();

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 16, 16 });

    P1<Real, Mesh<Context::MPI>> vh(mesh);

    auto solution = cos(pi * F::x) * cos(pi * F::y);
    auto f = 2 * pi * pi * solution;

    PETSc::Variational::TrialFunction u(vh);
    PETSc::Variational::TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, solution);

    PETSc::Solver::CG solver(poisson);
    solver.solve();

    P1<Real, Mesh<Context::MPI>> sh(mesh);
    GridFunction<P1<Real, Mesh<Context::MPI>>, ::Vec> diff(sh);
    diff = Pow(u.getSolution() - solution, 2);

    const Real localError = Integral(diff).compute();

    Real globalError = 0;
    boost::mpi::all_reduce(world, localError, globalError, std::plus<Real>());

    EXPECT_NEAR(globalError, 0, RODIN_FUZZY_CONSTANT);
  }
}

// ---------------------------------------------------------------------------
// main() — initializes MPI and PETSc before running all tests.
// ---------------------------------------------------------------------------
int main(int argc, char** argv)
{
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  g_env   = &env;
  g_world = &world;

  [[maybe_unused]] PetscErrorCode ierr = PetscInitialize(&argc, &argv, nullptr, nullptr);
  assert(ierr == PETSC_SUCCESS);

  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();

  PetscFinalize();
  return result;
}
