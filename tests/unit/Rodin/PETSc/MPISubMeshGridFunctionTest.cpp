/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 *
 * @file MPISubMeshGridFunctionTest.cpp
 * @brief PETSc distributed GridFunction tests for P0 and P1 FES on
 *        SubMesh<Context::MPI> and cross-mesh evaluation.
 *
 * These tests verify GridFunction projection and evaluation when the FES is
 * defined on a distributed SubMesh, or when a parent-mesh FES is evaluated
 * at points from a distributed SubMesh cell SubMesh (both meshes share the
 * same topological dimension).  They mirror the local-context tests in
 * tests/unit/Rodin/Variational/SubMeshGridFunctionTest.cpp but use
 * distributed meshes and the PETSc vector backend.
 *
 * ## Test groups
 *
 * ### Groups 1–4: SubMesh-FES GridFunction projection
 *  1. P0 FES on boundary SubMesh — constant projection, owned DOF check.
 *  2. P0 FES on boundary SubMesh — integral ∫_∂Ω 1 dΓ = |∂Ω| = 4.
 *  3. P1 FES on boundary SubMesh — constant projection, owned DOF check.
 *  4. P1 FES on boundary SubMesh — integral ∫_∂Ω 1 dΓ = 4.
 *
 * ### Groups 5–7: Cell SubMesh-FES GridFunction projection
 *  5. P0 FES on cell SubMesh — constant projection, owned DOF check.
 *  6. P0 FES on cell SubMesh — integral ∫_Ω c dΩ = c.
 *  7. P1 FES on cell SubMesh — constant projection, owned DOF check.
 *  8. P1 FES on cell SubMesh — integral ∫_Ω c dΩ = c.
 *
 * ### Groups 9–10: Cross-mesh evaluation (cell SubMesh ↔ parent mesh)
 *  9. Inclusion path: Parent-mesh P0 FES, TestFunction on cell SubMesh.
 *     Assembly iterates SubMesh cells and evaluates the parent GridFunction
 *     at SubMesh points → calls getValue(subMeshPoint) → inclusion path.
 * 10. Restriction path: SubMesh P1 FES, TestFunction on parent mesh.
 *     Assembly iterates parent cells and evaluates the SubMesh GridFunction
 *     at parent points → calls getValue(parentPoint) → restriction path.
 *
 * Run with: mpirun -n 1/2/3/4 as registered in CMakeLists.txt.
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
  // Helpers
  // -------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------

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
    // Scale to [0,1]×[0,1] so |∂Ω| = 4.
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
  // Group 3 — P1 on boundary SubMesh: constant projection, owned DOF check
  // ==========================================================================

  /**
   * @brief Project a constant onto distributed P1 on the boundary SubMesh.
   *        Each owned DOF must equal the constant.
   */
  TEST(PETSc_MPI_P1_SubMesh, BoundarySubMesh_ConstantProjection_OwnedDOFs_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const Real c = 3.14159;

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub  = makeBoundarySubMesh(mesh);

    P1<Real, Mesh<Context::MPI>> fes(sub);
    GridFunction<decltype(fes), ::Vec> u(fes);
    u = RealFunction(c);

    Index begin = 0, end = 0;
    fes.getOwnershipRange(begin, end);

    for (Index i = begin; i < end; ++i)
      EXPECT_NEAR(PetscRealPart(u[i]), c, 1e-10);
    u.flush();
  }

  // ==========================================================================
  // Group 4 — P1 on boundary SubMesh: integral = boundary length
  // ==========================================================================

  /**
   * @brief Project f=1 onto P1 on the boundary SubMesh of a unit square and
   *        verify ∫_∂Ω u_h dΓ = |∂Ω| = 4.
   */
  TEST(PETSc_MPI_P1_SubMesh, BoundarySubMesh_Integral_EqualsBoundaryLength_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(
        ctx, Polytope::Type::Triangle, {6, 6}, Real(1) / Real(5));
    auto sub  = makeBoundarySubMesh(mesh);

    P1<Real, Mesh<Context::MPI>> fes(sub);
    GridFunction<decltype(fes), ::Vec> u(fes);
    u = RealFunction(1.0);

    const Real total = Integral(u).compute();
    EXPECT_NEAR(total, 4.0, 1e-9);
  }

  // ==========================================================================
  // Group 5 — P0 on cell SubMesh: constant projection, owned DOF check
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
  // Group 6 — P0 on cell SubMesh: integral = c × |Ω|
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
  // Group 7 — P1 on cell SubMesh: constant projection, owned DOF check
  // ==========================================================================

  /**
   * @brief Project a constant onto P1 on the full-cell SubMesh.  Each owned
   *        DOF must equal the constant.
   */
  TEST(PETSc_MPI_P1_SubMesh, CellSubMesh_ConstantProjection_OwnedDOFs_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const Real c = 7.77;

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub  = makeCellSubMesh(mesh);

    P1<Real, Mesh<Context::MPI>> fes(sub);
    GridFunction<decltype(fes), ::Vec> u(fes);
    u = RealFunction(c);

    Index begin = 0, end = 0;
    fes.getOwnershipRange(begin, end);

    for (Index i = begin; i < end; ++i)
      EXPECT_NEAR(PetscRealPart(u[i]), c, 1e-10);
    u.flush();
  }

  // ==========================================================================
  // Group 8 — P1 on cell SubMesh: integral = c × |Ω|
  // ==========================================================================

  /**
   * @brief Project constant f=c onto P1 on a cell SubMesh of a unit square
   *        and verify ∫_Ω u_h dΩ = c.
   */
  TEST(PETSc_MPI_P1_SubMesh, CellSubMesh_Integral_EqualsArea_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const Real c = 1.41421;

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(
        ctx, Polytope::Type::Triangle, {6, 6}, Real(1) / Real(5));
    auto sub  = makeCellSubMesh(mesh);

    P1<Real, Mesh<Context::MPI>> fes(sub);
    GridFunction<decltype(fes), ::Vec> u(fes);
    u = RealFunction(c);

    const Real total = Integral(u).compute();
    EXPECT_NEAR(total, c, 1e-9);
  }

  // ==========================================================================
  // Group 9 — Inclusion path: parent P0 FES, TestFunction on cell SubMesh
  // ==========================================================================

  /**
   * @brief Project constant c onto P0 on the parent mesh.  Assemble a
   *        LinearForm with TestFunction on the cell SubMesh FES.
   *
   * Assembly iterates cell SubMesh cells and evaluates the parent
   * GridFunction at SubMesh quadrature points.  This calls
   * `parentGF.getValue(subPoint)` which goes through the inclusion path
   * (the SubMesh point's polytope is a SubMesh polytope, not a parent
   * polytope, so `fesMesh.inclusion(subPoint)` is used).
   *
   * Since the SubMesh covers the full domain and `parentGF = c`,
   * ∑_i b_i = c × |Ω|.  For a scaled [0,1]² mesh |Ω| = 1, so sum = c.
   */
  TEST(PETSc_MPI_P0_ParentMesh, InclusionPath_CellSubMeshTestFunction_VecSum_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const Real c = 5.5;

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(
        ctx, Polytope::Type::Triangle, {6, 6}, Real(1) / Real(5));
    auto sub  = makeCellSubMesh(mesh);

    // Parent-mesh GridFunction, constant c.
    P0<Real, Mesh<Context::MPI>> parentFes(mesh);
    GridFunction<decltype(parentFes), ::Vec> parentGF(parentFes);
    parentGF = RealFunction(c);

    // TestFunction on the cell SubMesh.
    P0<Real, Mesh<Context::MPI>> subFes(sub);
    PETSc::Variational::TestFunction v(subFes);

    // Build and assemble LinearForm: b_i = ∫_{K_i} parentGF · φ_i dK
    // Assembly evaluates parentGF at SubMesh quadrature points → inclusion.
    LinearForm lf(v);
    lf = Integral(parentGF, v);
    lf.assemble();

    // Sum all entries: ∑ b_i = ∫_Ω parentGF dΩ = c × 1 = c.
    PetscReal sum = 0.0;
    VecSum(lf.getVector(), &sum);
    EXPECT_NEAR(static_cast<Real>(sum), c, 1e-9);
  }

  /**
   * @brief Same inclusion-path test using P1 on the parent mesh.
   */
  TEST(PETSc_MPI_P1_ParentMesh, InclusionPath_CellSubMeshTestFunction_VecSum_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const Real c = 4.0;

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(
        ctx, Polytope::Type::Triangle, {6, 6}, Real(1) / Real(5));
    auto sub  = makeCellSubMesh(mesh);

    P1<Real, Mesh<Context::MPI>> parentFes(mesh);
    GridFunction<decltype(parentFes), ::Vec> parentGF(parentFes);
    parentGF = RealFunction(c);

    P1<Real, Mesh<Context::MPI>> subFes(sub);
    PETSc::Variational::TestFunction v(subFes);

    // Assembly iterates SubMesh cells, evaluates parentGF at SubMesh points.
    LinearForm lf(v);
    lf = Integral(parentGF, v);
    lf.assemble();

    PetscReal sum = 0.0;
    VecSum(lf.getVector(), &sum);
    EXPECT_NEAR(static_cast<Real>(sum), c, 1e-9);
  }

  // ==========================================================================
  // Group 10 — Restriction path: SubMesh P1 FES, TestFunction on parent mesh
  // ==========================================================================

  /**
   * @brief Project constant c onto P1 on the cell SubMesh.  Assemble a
   *        LinearForm with TestFunction on the parent mesh FES.
   *
   * Assembly iterates parent mesh cells and evaluates the SubMesh
   * GridFunction at parent quadrature points.  This calls
   * `subGF.getValue(parentPoint)` which goes through the restriction path:
   * the FES mesh is the SubMesh, the point lives in the parent mesh, so
   * `fesMesh.restriction(parentPoint)` is used.
   *
   * ∑_i b_i = c × |Ω| = c.
   */
  TEST(PETSc_MPI_P0_SubMesh, RestrictionPath_ParentMeshTestFunction_VecSum_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const Real c = 6.28;

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(
        ctx, Polytope::Type::Triangle, {6, 6}, Real(1) / Real(5));
    auto sub  = makeCellSubMesh(mesh);

    // SubMesh-FES GridFunction, constant c.
    P0<Real, Mesh<Context::MPI>> subFes(sub);
    GridFunction<decltype(subFes), ::Vec> subGF(subFes);
    subGF = RealFunction(c);

    // TestFunction on the parent mesh.
    P0<Real, Mesh<Context::MPI>> parentFes(mesh);
    PETSc::Variational::TestFunction v(parentFes);

    // Assembly iterates parent cells, evaluates subGF at parent points → restriction.
    LinearForm lf(v);
    lf = Integral(subGF, v);
    lf.assemble();

    PetscReal sum = 0.0;
    VecSum(lf.getVector(), &sum);
    EXPECT_NEAR(static_cast<Real>(sum), c, 1e-9);
  }

  /**
   * @brief Same restriction-path test with P1 on the SubMesh.
   */
  TEST(PETSc_MPI_P1_SubMesh, RestrictionPath_ParentMeshTestFunction_VecSum_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const Real c = 1.73205;

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(
        ctx, Polytope::Type::Triangle, {6, 6}, Real(1) / Real(5));
    auto sub  = makeCellSubMesh(mesh);

    P1<Real, Mesh<Context::MPI>> subFes(sub);
    GridFunction<decltype(subFes), ::Vec> subGF(subFes);
    subGF = RealFunction(c);

    P1<Real, Mesh<Context::MPI>> parentFes(mesh);
    PETSc::Variational::TestFunction v(parentFes);

    LinearForm lf(v);
    lf = Integral(subGF, v);
    lf.assemble();

    PetscReal sum = 0.0;
    VecSum(lf.getVector(), &sum);
    EXPECT_NEAR(static_cast<Real>(sum), c, 1e-9);
  }

  // ==========================================================================
  // Quadrilateral mesh variants
  // ==========================================================================

  /**
   * @brief P0 boundary SubMesh projection, Quadrilateral mesh variant.
   */
  TEST(PETSc_MPI_P0_SubMesh, BoundarySubMesh_ConstantProjection_OwnedDOFs_Quadrilateral)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const Real c = 1.618;

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Quadrilateral, {4, 4});
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

  /**
   * @brief P1 cell SubMesh integral, Quadrilateral mesh variant.
   */
  TEST(PETSc_MPI_P1_SubMesh, CellSubMesh_Integral_EqualsArea_Quadrilateral)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const Real c = 2.71828;

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(
        ctx, Polytope::Type::Quadrilateral, {6, 6}, Real(1) / Real(5));
    auto sub  = makeCellSubMesh(mesh);

    P1<Real, Mesh<Context::MPI>> fes(sub);
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

  PetscInitialize(&argc, &argv, nullptr, nullptr);
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
