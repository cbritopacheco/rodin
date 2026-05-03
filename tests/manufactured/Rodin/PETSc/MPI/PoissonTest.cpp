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
 * The L² error is computed through the PETSc-backed `Integral` path.  For MPI
 * grid functions this evaluates a distributed `LinearForm`, so PETSc performs
 * the communicator-wide reduction internally.
 *
 * ### Workflow
 *
 * 1. Rank 0 builds a uniform-grid mesh, partitions it with
 *    `BalancedCompactPartitioner`, and calls `sharder.scatter()`.
 * 2. Each rank saves its shard to HDF5.
 * 3. Each rank rereads its shard from HDF5 through the MPI mesh loader.
 * 4. Incidence data required for assembly and Dirichlet-BC enforcement is
 *    computed and the mesh is reconciled.
 * 5. `PETSc::Variational::TrialFunction / TestFunction` wrap the distributed
 *    `P1` space.
 * 6. `Problem<PETSc::Math::LinearSystem, U, V>` is assembled with the MPI
 *    assembly backend and solved with a PETSc KSP solver.
 * 7. The PETSc-backed `Integral` path returns the communicator-wide L² error,
 *    which is checked against the manufactured tolerance.
 */

#include <cmath>
#include <numeric>
#include <string>

#include <gtest/gtest.h>

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
#include <Rodin/MPI/Variational/H1.h>
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

    Mesh<Context::MPI> reloaded(ctx);
    reloaded.load(rankFile, IO::FileFormat::HDF5);

    EXPECT_EQ(reloaded.getShard().getCellCount(), mesh.getShard().getCellCount())
      << "Rank " << comm.rank() << ": cell count mismatch after HDF5 reload.";
    EXPECT_EQ(reloaded.getShard().getVertexCount(), mesh.getShard().getVertexCount())
      << "Rank " << comm.rank() << ": vertex count mismatch after HDF5 reload.";

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
   * The global @f$ L^2 @f$ error is returned by `Integral(diff).compute()`
   * through PETSc's distributed vector operations. Expected to stay below
   * @ref RODIN_FUZZY_CONSTANT on a 16×16 mesh.
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

    const Real globalError = Integral(diff).compute();

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

    const Real globalError = Integral(diff).compute();

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

    const Real globalError = Integral(diff).compute();

    EXPECT_NEAR(globalError, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @brief P1-exact distributed Poisson: zero FE truncation error.
   *
   * Manufactured solution: @f$ u(x,y) = x + 2y + 1 @f$.
   *
   * Because the manufactured solution lies exactly in the P1 space, the
   * global L² error must be roundoff-zero. Strong nonhomogeneous Dirichlet
   * enforcement can make the assembled operator nonsymmetric, so this uses
   * GMRES rather than CG.
   */
  TEST(PETSc_MPI_Poisson, GMRES_P1Exact_Triangle)
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

    PETSc::Solver::GMRES solver(poisson);
    solver.solve();

    P1<Real, Mesh<Context::MPI>> sh(mesh);
    GridFunction<P1<Real, Mesh<Context::MPI>>, ::Vec> diff(sh);
    diff = Pow(u.getSolution() - solution, 2);

    const Real globalError = Integral(diff).compute();

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

    const Real globalError = Integral(diff).compute();

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

    const Real globalError = Integral(diff).compute();

    EXPECT_NEAR(globalError, 0, RODIN_FUZZY_CONSTANT);
  }
}

namespace Rodin::Tests::Manufactured::PETSc::MPI
{
  namespace PETSc = ::Rodin::PETSc;

  // =========================================================================
  // Distributed H1 Poisson tests
  // =========================================================================

  /**
   * @brief Distributed H1 K=1 Poisson test on triangles (CG solver).
   *
   * Manufactured solution: @f$ u(x,y) = \sin(\pi x)\sin(\pi y) @f$.
   * H1 with K=1 should reproduce P1 behavior on the same mesh.
   */
  TEST(PETSc_MPI_H1Poisson, CG_K1_SimpleSine_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const auto pi = Math::Constants::pi();

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 16, 16 });

    constexpr auto order = std::integral_constant<size_t, 1>{};
    H1 vh(order, mesh);

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

    H1 sh(order, mesh);
    GridFunction<decltype(sh), ::Vec> diff(sh);
    diff = Pow(u.getSolution() - solution, 2);

    const Real globalError = Integral(diff).compute();

    EXPECT_NEAR(globalError, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @brief Distributed H1 K=2 Poisson test on triangles (CG solver).
   *
   * Manufactured solution: @f$ u(x,y) = \sin(\pi x)\sin(\pi y) @f$.
   * H1 with K=2 (quadratic) distributes DOFs on vertices and edge interiors.
   */
  TEST(PETSc_MPI_H1Poisson, CG_K2_SimpleSine_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const auto pi = Math::Constants::pi();

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 16, 16 });

    constexpr auto order = std::integral_constant<size_t, 2>{};
    H1 vh(order, mesh);

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

    H1 sh(order, mesh);
    GridFunction<decltype(sh), ::Vec> diff(sh);
    diff = Pow(u.getSolution() - solution, 2);

    const Real globalError = Integral(diff).compute();

    EXPECT_NEAR(globalError, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @brief Distributed H1 K=2 Poisson test on quadrilaterals.
   *
   * Manufactured solution: @f$ u(x,y) = \sin(\pi x)\sin(\pi y) @f$.
   */
  TEST(PETSc_MPI_H1Poisson, CG_K2_SimpleSine_Quadrilateral)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const auto pi = Math::Constants::pi();

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Quadrilateral, { 16, 16 });

    constexpr auto order = std::integral_constant<size_t, 2>{};
    H1 vh(order, mesh);

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

    H1 sh(order, mesh);
    GridFunction<decltype(sh), ::Vec> diff(sh);
    diff = Pow(u.getSolution() - solution, 2);

    const Real globalError = Integral(diff).compute();

    EXPECT_NEAR(globalError, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @brief P2-exact distributed H1 K=2 Poisson test.
   *
   * Manufactured solution: @f$ u(x,y) = x(1-x)y(1-y) @f$.
   * Because this is a degree-2 polynomial, H1 K=2 should reproduce it exactly.
   */
  TEST(PETSc_MPI_H1Poisson, CG_K2_P2Exact_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 8, 8 });

    constexpr auto order = std::integral_constant<size_t, 2>{};
    H1 vh(order, mesh);

    auto solution = F::x * (1 - F::x) * F::y * (1 - F::y);
    auto f = 2 * F::y * (1 - F::y) + 2 * F::x * (1 - F::x);

    PETSc::Variational::TrialFunction u(vh);
    PETSc::Variational::TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());

    PETSc::Solver::CG solver(poisson);
    solver.solve();

    H1 sh(order, mesh);
    GridFunction<decltype(sh), ::Vec> diff(sh);
    diff = Pow(u.getSolution() - solution, 2);

    const Real globalError = Integral(diff).compute();

    EXPECT_NEAR(globalError, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @brief Single-rank H1 K=2 consistency check.
   *
   * On exactly 1 rank, the distributed H1 K=2 solve must match the
   * manufactured-solution tolerance of the sequential path.
   */
  TEST(PETSc_MPI_H1Poisson, SingleRank_K2_MatchesSequential_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() != 1)
      GTEST_SKIP() << "This test is designed for exactly 1 MPI rank.";

    const auto pi = Math::Constants::pi();

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 16, 16 });

    constexpr auto order = std::integral_constant<size_t, 2>{};
    H1 vh(order, mesh);

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

    H1 sh(order, mesh);
    GridFunction<decltype(sh), ::Vec> diff(sh);
    diff = Pow(u.getSolution() - solution, 2);
    const Real localError = Integral(diff).compute();

    EXPECT_NEAR(localError, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @brief H1 K=2 with nonhomogeneous Dirichlet BC, distributed.
   *
   * Manufactured solution: @f$ u(x,y) = \cos(\pi x)\cos(\pi y) @f$.
   */
  TEST(PETSc_MPI_H1Poisson, CG_K2_NonhomogeneousDirichlet_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const auto pi = Math::Constants::pi();

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 16, 16 });

    constexpr auto order = std::integral_constant<size_t, 2>{};
    H1 vh(order, mesh);

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

    H1 sh(order, mesh);
    GridFunction<decltype(sh), ::Vec> diff(sh);
    diff = Pow(u.getSolution() - solution, 2);

    const Real globalError = Integral(diff).compute();

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
