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
 * These tests verify GridFunction projection and numerical integration when
 * the FES is defined on a distributed SubMesh (both boundary skin and
 * full-cell sub-regions) and when a parent-mesh FES is evaluated at points
 * from a distributed cell SubMesh (inclusion and restriction paths).  They
 * mirror the local-context tests in
 * tests/unit/Rodin/Variational/SubMeshGridFunctionTest.cpp but use
 * distributed meshes and the PETSc vector backend.
 *
 * ## Test groups
 *
 *  1. P0 FES on boundary SubMesh — constant projection, owned DOF check.
 *  2. P0 FES on boundary SubMesh — integral \int_\partial\Omega 1 d\Gamma = |\partial\Omega| = 4.
 *  3. P1 FES on boundary SubMesh — constant projection, owned DOF check.
 *  4. P1 FES on boundary SubMesh — integral \int_\partial\Omega 1 d\Gamma = 4.
 *  5. P0 FES on cell SubMesh — constant projection, owned DOF check.
 *  6. P0 FES on cell SubMesh — integral \int_\Omega c d\Omega = c.
 *  7. P1 FES on cell SubMesh — constant projection, owned DOF check.
 *  8. P1 FES on cell SubMesh — integral \int_\Omega c d\Omega = c.
 *  9. Inclusion path: parent-mesh P0/P1 GridFunction evaluated at SubMesh
 *     quadrature points via TestFunction on the cell SubMesh.
 * 10. Restriction path: SubMesh P0/P1 GridFunction evaluated at parent-mesh
 *     quadrature points via TestFunction on the parent mesh.
 * 11. Quadrilateral mesh variants.
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
#include <Rodin/MPI/Variational/P0g.h>
#include <Rodin/MPI/Variational/P0.h>
#include <Rodin/MPI/Variational/P1.h>
#include <Rodin/MPI/Variational/H1.h>
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

  TEST(PETSc_MPI_P0g_SubMesh, BoundarySubMesh_ConstantProjection_Integral_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const Real c = 2.25;

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(
        ctx, Polytope::Type::Triangle, {6, 6}, Real(1) / Real(5));
    auto sub  = makeBoundarySubMesh(mesh);

    P0g<Real, Mesh<Context::MPI>> fes(sub);
    GridFunction<decltype(fes), ::Vec> u(fes);
    u = RealFunction(c);

    const Real total = Integral(u).compute();
    EXPECT_NEAR(total, 4.0 * c, 1e-9);
  }

  TEST(PETSc_MPI_P0g_SubMesh, CellSubMesh_InclusionAndRestrictionPaths_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const Real parentValue = 1.75;
    const Real subValue = 3.5;

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(
        ctx, Polytope::Type::Triangle, {6, 6}, Real(1) / Real(5));
    auto sub  = makeCellSubMesh(mesh);

    P0g<Real, Mesh<Context::MPI>> parentFes(mesh);
    GridFunction<decltype(parentFes), ::Vec> parentGF(parentFes);
    parentGF = RealFunction(parentValue);

    P0g<Real, Mesh<Context::MPI>> subFes(sub);
    PETSc::Variational::TestFunction subTest(subFes);
    LinearForm inclusion(subTest);
    inclusion = Integral(parentGF, subTest);
    inclusion.assemble();

    PetscReal inclusionSum = 0.0;
    VecSum(inclusion.getVector(), &inclusionSum);
    EXPECT_NEAR(static_cast<Real>(inclusionSum), parentValue, 1e-9);

    GridFunction<decltype(subFes), ::Vec> subGF(subFes);
    subGF = RealFunction(subValue);

    PETSc::Variational::TestFunction parentTest(parentFes);
    LinearForm restriction(parentTest);
    restriction = Integral(subGF, parentTest);
    restriction.assemble();

    PetscReal restrictionSum = 0.0;
    VecSum(restriction.getVector(), &restrictionSum);
    EXPECT_NEAR(static_cast<Real>(restrictionSum), subValue, 1e-9);
  }

  TEST(PETSc_MPI_H1_SubMesh, BoundarySubMesh_ConstantProjection_Integral_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(
        ctx, Polytope::Type::Triangle, {6, 6}, Real(1) / Real(5));
    auto sub  = makeBoundarySubMesh(mesh);

    H1<1, Real, Mesh<Context::MPI>> fes(std::integral_constant<size_t, 1>{}, sub);
    GridFunction<decltype(fes), ::Vec> u(fes);
    u = RealFunction(1.0);

    const Real total = Integral(u).compute();
    EXPECT_NEAR(total, 4.0, 1e-9);
  }

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
   *        square and verify \int_\partial\Omega u_h d\Gamma = |\partial\Omega| = 4.
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
  // Group 4 — P0 on cell SubMesh: integral = c × |\Omega|
  // ==========================================================================

  /**
   * @brief Project constant f=c onto P0 on a cell SubMesh of a unit square
   *        and verify \int_\Omega u_h d\Omega = c × |\Omega| = c.
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
   * @brief P0 boundary SubMesh on Quadrilateral mesh: integral = |\partial\Omega| = 4.
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
   * @brief P0 cell SubMesh on Quadrilateral mesh: integral = c × |\Omega| = c.
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
   *        verify \int_\partial\Omega u_h d\Gamma = |\partial\Omega| = 4.
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
  // Group 8 — P1 on cell SubMesh: integral = c × |\Omega|
  // ==========================================================================

  /**
   * @brief Project constant f=c onto P1 on a cell SubMesh of a unit square
   *        and verify \int_\Omega u_h d\Omega = c.
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
  // Group 9 — Inclusion path: parent-mesh FES evaluated at SubMesh points
  // ==========================================================================

  /**
   * @brief Project constant c onto P0 on the parent mesh.  Assemble a
   *        LinearForm with TestFunction on the cell SubMesh FES.
   *
   * Assembly iterates cell SubMesh cells and evaluates the parent
   * GridFunction at SubMesh quadrature points.  This exercises the inclusion
   * path: `parentGF.getValue(subPoint)` calls `fesMesh.inclusion(subPoint)`
   * which walks the SubMesh ancestor chain to find the parent polytope.
   *
   * Since the SubMesh covers the full domain and `parentGF = c`,
   * ∑_i b_i = c × |\Omega| = c for a scaled [0,1]² mesh.
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

    P0<Real, Mesh<Context::MPI>> parentFes(mesh);
    GridFunction<decltype(parentFes), ::Vec> parentGF(parentFes);
    parentGF = RealFunction(c);

    P0<Real, Mesh<Context::MPI>> subFes(sub);
    PETSc::Variational::TestFunction v(subFes);

    LinearForm lf(v);
    lf = Integral(parentGF, v);
    lf.assemble();

    PetscReal sum = 0.0;
    VecSum(lf.getVector(), &sum);
    EXPECT_NEAR(static_cast<Real>(sum), c, 1e-9);
  }

  /**
   * @brief Same inclusion-path test using P1 on both the parent mesh and the
   *        cell SubMesh.
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

    LinearForm lf(v);
    lf = Integral(parentGF, v);
    lf.assemble();

    PetscReal sum = 0.0;
    VecSum(lf.getVector(), &sum);
    EXPECT_NEAR(static_cast<Real>(sum), c, 1e-9);
  }

  // ==========================================================================
  // Group 10 — Restriction path: SubMesh FES evaluated at parent-mesh points
  // ==========================================================================

  /**
   * @brief Project constant c onto P0 on the cell SubMesh.  Assemble a
   *        LinearForm with TestFunction on the parent mesh FES.
   *
   * Assembly iterates parent mesh cells and evaluates the SubMesh
   * GridFunction at parent quadrature points.  This exercises the restriction
   * path: `subGF.getValue(parentPoint)` sees that the FES mesh is a SubMesh
   * and calls `fesMesh.restriction(parentPoint)` to map the parent polytope
   * index to the SubMesh polytope index before interpolating.
   *
   * ∑_i b_i = c × |\Omega| = c for a scaled [0,1]² mesh.
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

    P0<Real, Mesh<Context::MPI>> subFes(sub);
    GridFunction<decltype(subFes), ::Vec> subGF(subFes);
    subGF = RealFunction(c);

    P0<Real, Mesh<Context::MPI>> parentFes(mesh);
    PETSc::Variational::TestFunction v(parentFes);

    LinearForm lf(v);
    lf = Integral(subGF, v);
    lf.assemble();

    PetscReal sum = 0.0;
    VecSum(lf.getVector(), &sum);
    EXPECT_NEAR(static_cast<Real>(sum), c, 1e-9);
  }

  /**
   * @brief Same restriction-path test with P1 on the cell SubMesh.
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
  // Group 11 — Quadrilateral mesh variants
  // ==========================================================================

  /**
   * @brief P0 boundary SubMesh constant projection, Quadrilateral mesh.
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

  PetscErrorCode ierr = PetscInitialize(&argc, &argv, nullptr, nullptr);
  if (ierr != PETSC_SUCCESS)
    return static_cast<int>(ierr);

  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();

  PetscFinalize();
  return result;
}
