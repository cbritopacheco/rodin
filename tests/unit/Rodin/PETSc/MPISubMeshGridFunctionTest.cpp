/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 *
 * @file MPISubMeshGridFunctionTest.cpp
 * @brief PETSc distributed GridFunction tests for P0 FES on
 *        SubMesh<Context::MPI>.
 *
 * These tests verify GridFunction projection and numerical integration when
 * the P0 FES is defined on a distributed SubMesh (both boundary skin and
 * full-cell sub-regions).  They mirror the local-context tests in
 * tests/unit/Rodin/Variational/SubMeshGridFunctionTest.cpp but use
 * distributed meshes and the PETSc vector backend.
 *
 * NOTE: P1 on distributed SubMesh is not yet fully supported at np>=2 due
 * to a crash in the P1 distributed constructor when ghost vertex ownership
 * cannot be resolved across SubMesh shard boundaries.  P1 and cross-mesh
 * tests are excluded until that limitation is resolved.
 *
 * ## Test groups
 *
 *  1. P0 FES on boundary SubMesh — constant projection, owned DOF check.
 *  2. P0 FES on boundary SubMesh — integral ∫_∂Ω 1 dΓ = |∂Ω| = 4.
 *  3. P0 FES on cell SubMesh — constant projection, owned DOF check.
 *  4. P0 FES on cell SubMesh — integral ∫_Ω c dΩ = c.
 *  5. Quadrilateral mesh variant: P0 boundary SubMesh, integral = 4.
 *
 * Run with: mpirun -n 1/2/4 as registered in CMakeLists.txt.
 */

#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include <petsc.h>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>

#include <Rodin/Types.h>
#include <Rodin/Geometry.h>
#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#include <Rodin/Variational.h>
#include <Rodin/MPI/Context/MPI.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/Geometry/SubMesh.h>
#include <Rodin/MPI/Variational/P0.h>
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
  static Mesh<Context::Local> makeShardableMesh(
      Polytope::Type type,
      std::initializer_list<size_t> shape,
      Real scale = Real(1))
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(type, shape);
    if (scale != Real(1))
      mesh.scale(scale);
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);
    mesh.getConnectivity().compute(D, 0);
    mesh.getConnectivity().compute(D, D - 1);
    mesh.getConnectivity().compute(D - 1, D);
    mesh.getConnectivity().compute(D - 1, 0);
    return mesh;
  }

  static Mesh<Context::MPI> distributeFromRoot(
      const Context::MPI& ctx,
      Polytope::Type type,
      std::initializer_list<size_t> shape,
      Real scale = Real(1))
  {
    const auto& comm = ctx.getCommunicator();
    Sharder<Context::MPI> sharder(ctx);
    if (comm.rank() == 0)
    {
      auto localMesh = makeShardableMesh(type, shape, scale);
      BalancedCompactPartitioner partitioner(localMesh);
      partitioner.partition(static_cast<size_t>(comm.size()));
      sharder.shard(partitioner);
      sharder.scatter(0);
    }
    return sharder.gather(0);
  }

  static SubMesh<Context::MPI> makeBoundarySubMesh(
      const Mesh<Context::MPI>& mesh)
  {
    const size_t faceDim = mesh.getDimension() - 1;
    SubMesh<Context::MPI>::Builder builder;
    builder.initialize(mesh);
    for (auto it = mesh.getBoundary(); it; ++it)
      builder.include(faceDim, it->getIndex());
    return builder.finalize();
  }

  static SubMesh<Context::MPI> makeCellSubMesh(
      const Mesh<Context::MPI>& mesh)
  {
    const size_t cellDim = mesh.getDimension();
    const auto& shard = mesh.getShard();
    SubMesh<Context::MPI>::Builder builder;
    builder.initialize(mesh);
    for (Index i = 0; i < static_cast<Index>(shard.getCellCount()); ++i)
    {
      if (shard.isOwned(cellDim, i))
        builder.include(cellDim, i);
    }
    return builder.finalize();
  }
} // anonymous namespace

namespace Rodin::Tests::Unit::PETSc::MPI
{
  namespace PETSc = ::Rodin::PETSc;

  // ==========================================================================
  // Group 1 — P0 on boundary SubMesh: constant projection, owned DOF check
  // ==========================================================================

  /**
   * @brief Project a constant function onto distributed P0 on the boundary
   *        SubMesh.  Each owned DOF must equal the constant.
   */
  TEST(PETSc_MPI_P0_SubMesh, BoundarySubMesh_ConstantProjection_OwnedDOFs_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const Real c = 2.71828;

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub  = makeBoundarySubMesh(mesh);

    P0<Real, Mesh<Context::MPI>> fes(sub);
    GridFunction<decltype(fes), ::Vec> u(fes);
    u = RealFunction(c);

    Index begin = 0, end = 0;
    fes.getOwnershipRange(begin, end);

    for (Index i = begin; i < end; ++i)
      EXPECT_NEAR(PetscRealPart(u[i]), c, 1e-10);
    u.flush();
  }

  // ==========================================================================
  // Group 2 — P0 on boundary SubMesh: integral = boundary length
  // ==========================================================================

  /**
   * @brief Project constant f=1 onto P0 on the boundary SubMesh of a unit
   *        square and verify ∫_∂Ω u_h dΓ = |∂Ω| = 4.
   */
  TEST(PETSc_MPI_P0_SubMesh, BoundarySubMesh_Integral_EqualsBoundaryLength_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(
        ctx, Polytope::Type::Triangle, {6, 6}, Real(1) / Real(5));
    auto sub  = makeBoundarySubMesh(mesh);

    P0<Real, Mesh<Context::MPI>> fes(sub);
    GridFunction<decltype(fes), ::Vec> u(fes);
    u = RealFunction(1.0);

    const Real total = Integral(u).compute();
    EXPECT_NEAR(total, 4.0, 1e-9);
  }

  // ==========================================================================
  // Group 3 — P0 on cell SubMesh: constant projection, owned DOF check
  // ==========================================================================

  /**
   * @brief Project a constant onto P0 on the full-cell SubMesh.  Each owned
   *        DOF must equal the constant.
   */
  TEST(PETSc_MPI_P0_SubMesh, CellSubMesh_ConstantProjection_OwnedDOFs_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const Real c = 9.81;

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub  = makeCellSubMesh(mesh);

    P0<Real, Mesh<Context::MPI>> fes(sub);
    GridFunction<decltype(fes), ::Vec> u(fes);
    u = RealFunction(c);

    Index begin = 0, end = 0;
    fes.getOwnershipRange(begin, end);

    for (Index i = begin; i < end; ++i)
      EXPECT_NEAR(PetscRealPart(u[i]), c, 1e-10);
    u.flush();
  }

  // ==========================================================================
  // Group 4 — P0 on cell SubMesh: integral = c × |Ω|
  // ==========================================================================

  /**
   * @brief Project constant f=c onto P0 on a cell SubMesh of a unit square
   *        and verify ∫_Ω u_h dΩ = c × |Ω| = c.
   */
  TEST(PETSc_MPI_P0_SubMesh, CellSubMesh_Integral_EqualsArea_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const Real c = 3.0;

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(
        ctx, Polytope::Type::Triangle, {6, 6}, Real(1) / Real(5));
    auto sub  = makeCellSubMesh(mesh);

    P0<Real, Mesh<Context::MPI>> fes(sub);
    GridFunction<decltype(fes), ::Vec> u(fes);
    u = RealFunction(c);

    const Real total = Integral(u).compute();
    EXPECT_NEAR(total, c, 1e-9);
  }

  // ==========================================================================
  // Group 5 — Quadrilateral mesh variant
  // ==========================================================================

  /**
   * @brief P0 boundary SubMesh on Quadrilateral mesh: integral = |∂Ω| = 4.
   */
  TEST(PETSc_MPI_P0_SubMesh, BoundarySubMesh_Integral_EqualsBoundaryLength_Quadrilateral)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(
        ctx, Polytope::Type::Quadrilateral, {6, 6}, Real(1) / Real(5));
    auto sub  = makeBoundarySubMesh(mesh);

    P0<Real, Mesh<Context::MPI>> fes(sub);
    GridFunction<decltype(fes), ::Vec> u(fes);
    u = RealFunction(1.0);

    const Real total = Integral(u).compute();
    EXPECT_NEAR(total, 4.0, 1e-9);
  }

  /**
   * @brief P0 cell SubMesh on Quadrilateral mesh: integral = c × |Ω| = c.
   */
  TEST(PETSc_MPI_P0_SubMesh, CellSubMesh_Integral_EqualsArea_Quadrilateral)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const Real c = 2.5;

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(
        ctx, Polytope::Type::Quadrilateral, {6, 6}, Real(1) / Real(5));
    auto sub  = makeCellSubMesh(mesh);

    P0<Real, Mesh<Context::MPI>> fes(sub);
    GridFunction<decltype(fes), ::Vec> u(fes);
    u = RealFunction(c);

    const Real total = Integral(u).compute();
    EXPECT_NEAR(total, c, 1e-9);
  }
}

// ---------------------------------------------------------------------------
// main() — initializes PETSc and MPI environment used by all tests.
// ---------------------------------------------------------------------------
int main(int argc, char** argv)
{
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  g_env   = &env;
  g_world = &world;

  PetscErrorCode ierr = PetscInitialize(&argc, &argv, nullptr, nullptr);
  if (ierr != PETSC_SUCCESS)
    return static_cast<int>(ierr);

  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();

  PetscFinalize();
  return result;
}
