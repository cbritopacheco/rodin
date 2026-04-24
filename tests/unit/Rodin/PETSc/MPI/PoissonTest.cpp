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
 * 2. All ranks call `sharder.gather()` to obtain their local shard as a
 *    `Mesh<Context::MPI>`.
 * 3. Incidence data required for assembly and Dirichlet-BC enforcement is
 *    computed and reconciled.
 * 4. `PETSc::Variational::TrialFunction / TestFunction` wrap the distributed
 *    `P1` space.
 * 5. `Problem<PETSc::Math::LinearSystem, U, V>` is assembled with the MPI
 *    assembly backend and solved with a PETSc KSP solver.
 * 6. The local per-rank L² error is summed globally and checked against the
 *    manufactured tolerance.
 */

#include <cmath>
#include <numeric>

#include <gtest/gtest.h>

#include <petsc.h>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>

#include <Rodin/Configure.h>
#include <Rodin/Types.h>
#include <Rodin/Geometry.h>
#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#include <Rodin/Variational.h>
#include <Rodin/MPI.h>
#include <Rodin/MPI/Context/MPI.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/Variational/P1.h>
#include <Rodin/PETSc.h>

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

  /**
   * @brief Distributes a uniform 2D grid from rank 0 to all ranks via
   *        `BalancedCompactPartitioner` + `Sharder<Context::MPI>`.
   *
   * After this function returns, every rank holds its shard of the mesh.
   * The function also calls `compute(1,2)` and `reconcile(1)` on the
   * distributed mesh so that it is ready for FE assembly.
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
    mesh.getConnectivity().compute(1, 2);
    mesh.reconcile(1);
    return mesh;
  }
}

namespace Rodin::Tests::Unit::PETSc::MPI
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

    EXPECT_NEAR(globalError, 0, 1e-12);
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
