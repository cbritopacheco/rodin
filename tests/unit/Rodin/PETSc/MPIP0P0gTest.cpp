/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 *
 * @file MPIP0P0gTest.cpp
 * @brief Distributed PETSc tests for P0 and P0g spaces on MPI meshes.
 *
 * These tests distribute a mesh across MPI ranks and exercise:
 *   - GridFunction projection onto distributed P0 and P0g.
 *   - Linear form assembly (\int f\cdotv d\Omega) with distributed P0 and P0g.
 *   - Bilinear form assembly (\int u\cdotv d\Omega) with distributed P0 and P0g.
 *   - Solving the L2 projection problem with distributed P0 via PETSc GMRES.
 *   - Solving the L2 projection problem with distributed P0g via PETSc GMRES.
 *
 * ### P0 tests
 *
 * For the P0 space (piecewise constant), the natural PDE is the L2 projection:
 * find u_h ∈ P0 such that \int u_h v d\Omega = \int f v d\Omega for all v ∈ P0.
 * The resulting system is diagonal with M_ii = |K_i| (cell area/volume),
 * so u_h = f_K (average of f over each cell K).
 *
 * ### P0g tests
 *
 * For the P0g space (global constant), the L2 projection gives:
 * |\Omega| \cdot c = \int f d\Omega, so c = (1/|\Omega|) \int f d\Omega.
 * For a constant f this reduces to c = f.
 */

#include <cmath>
#include <string>

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
#include <Rodin/MPI/Context/MPI.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/Variational/P0.h>
#include <Rodin/MPI/Variational/P0g.h>
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
} // anonymous namespace

// ---------------------------------------------------------------------------
// P0 — GridFunction projection
// ---------------------------------------------------------------------------
namespace Rodin::Tests::Unit::PETSc::MPI
{
  // Alias to resolve the enclosing 'PETSc' namespace scope to Rodin::PETSc.
  namespace PETSc = ::Rodin::PETSc;

  // =========================================================================
  // P0 — GridFunction projection
  // =========================================================================

  /**
   * @brief Project a constant function onto distributed P0 and verify each
   *        owned DOF equals the constant.
   *
   * For a constant f = c and P0 (piecewise constant), the L2 projection
   * trivially gives u_K = c on every cell K.
   */
  TEST(PETSc_MPI_P0, GridFunctionProjection_ConstantFunction_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const Real c = 3.14;

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 6, 6 });

    P0<Real, Mesh<Context::MPI>> fes(mesh);
    GridFunction<decltype(fes), ::Vec> u(fes);
    u = RealFunction(c);

    Index begin = 0;
    Index end   = 0;
    fes.getOwnershipRange(begin, end);

    // Acquire read access to the PETSc vector.
    for (Index i = begin; i < end; ++i)
    {
      const PetscScalar v = u[i];
      EXPECT_NEAR(PetscRealPart(v), c, 1e-10);
    }
    u.flush();
  }

  /**
   * @brief The global integral of the projected constant function equals
   *        the constant times the domain area.
   *
   * \int u_h d\Omega = c \cdot |\Omega|.
   */
  TEST(PETSc_MPI_P0, Integral_ProjectedConstant_EqualsArea_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const Real c = 2.0;

    Context::MPI ctx(*g_env, world);
    // Scale to [0,1]x[0,1] → area = 1.
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 6, 6 },
                                   Real(1) / Real(5));

    P0<Real, Mesh<Context::MPI>> fes(mesh);
    GridFunction<decltype(fes), ::Vec> u(fes);
    u = RealFunction(c);

    // Integral of the grid function over the domain
    const Real total = Integral(u).compute();

    // Expected: c * |\Omega| = c * 1 = c
    EXPECT_NEAR(total, c, 1e-9);
  }

  // =========================================================================
  // P0 — Linear form assembly
  // =========================================================================

  /**
   * @brief Assemble \int 1\cdotv d\Omega on distributed P0 and verify the global
   *        vector size equals the number of cells.
   *
   * For P0 each entry corresponds to one cell, so the vector size must equal
   * the global cell count.
   */
  TEST(PETSc_MPI_P0, LinearForm_GlobalVectorSize_EqualsCellCount_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 6, 6 });

    P0<Real, Mesh<Context::MPI>> fes(mesh);
    PETSc::Variational::TestFunction v(fes);

    LinearForm lf(v);
    lf = Integral(RealFunction(1.0), v);
    lf.assemble();

    ::Vec b = lf.getVector();
    PetscInt globalSize = 0;
    VecGetSize(b, &globalSize);

    // getCellCount() is collective — must be called on all ranks.
    const size_t globalCells = mesh.getCellCount();
    EXPECT_EQ(static_cast<size_t>(globalSize), globalCells);
  }

  /**
   * @brief Sum of \int 1\cdotv d\Omega over all P0 test functions equals domain area.
   *
   * For f = 1 and P0, each entry b_K = |K| (cell measure).
   * ∑_K b_K = |\Omega|.
   */
  TEST(PETSc_MPI_P0, LinearForm_SumEqualsDomainArea_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    // Scale to [0,1]x[0,1] → area = 1.
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 6, 6 },
                                   Real(1) / Real(5));

    P0<Real, Mesh<Context::MPI>> fes(mesh);
    PETSc::Variational::TestFunction v(fes);

    LinearForm lf(v);
    lf = Integral(RealFunction(1.0), v);
    lf.assemble();

    ::Vec b = lf.getVector();
    PetscReal sum = 0.0;
    VecSum(b, &sum);

    EXPECT_NEAR(static_cast<Real>(sum), 1.0, 1e-10);
  }

  // =========================================================================
  // P0 — Bilinear form assembly
  // =========================================================================

  /**
   * @brief Assemble the P0 mass matrix \int u\cdotv d\Omega and check global dimensions.
   *
   * For P0 the global matrix size is (N_cells × N_cells).
   */
  TEST(PETSc_MPI_P0, BilinearForm_MassMatrixDimensions_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 6, 6 });

    P0<Real, Mesh<Context::MPI>> fes(mesh);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction  v(fes);

    BilinearForm bf(u, v);
    bf = Integral(u, v);
    bf.assemble();

    ::Mat A = bf.getOperator();
    PetscInt globalRows = 0;
    PetscInt globalCols = 0;
    MatGetSize(A, &globalRows, &globalCols);

    // getCellCount() is collective — must be called on all ranks.
    const size_t globalCells = mesh.getCellCount();
    EXPECT_EQ(static_cast<size_t>(globalRows), globalCells);
    EXPECT_EQ(static_cast<size_t>(globalCols), globalCells);
  }

  // =========================================================================
  // P0 — Solve L2 projection
  // =========================================================================

  /**
   * @brief Solve the L2 projection of a constant onto P0 via PETSc GMRES.
   *
   * Problem: find u_h ∈ P0 s.t. \int u_h v d\Omega = \int c v d\Omega for all v ∈ P0.
   * Solution: u_h = c on every cell.  L2 error = 0.
   */
  TEST(PETSc_MPI_P0, SolveL2Projection_Constant_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const Real c = 5.0;

    Context::MPI ctx(*g_env, world);
    // Scale to [0,1]x[0,1] → area = 1.
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 6, 6 },
                                   Real(1) / Real(5));

    using FES = P0<Real, Mesh<Context::MPI>>;
    FES fes(mesh);

    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction  v(fes);

    Problem projection(u, v);
    projection = Integral(u, v) - Integral(RealFunction(c), v);

    PETSc::Solver::GMRES solver(projection);
    solver.solve();

    // Compute \int (u_h - c)^2 d\Omega
    FES sh(mesh);
    GridFunction<FES, ::Vec> diff(sh);
    diff = Pow(u.getSolution() - RealFunction(c), 2);

    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0.0, 1e-10);
  }

  /**
   * @brief Solve the L2 projection of a constant onto P0 via PETSc GMRES,
   *        Quadrilateral mesh variant.
   */
  TEST(PETSc_MPI_P0, SolveL2Projection_Constant_Quadrilateral)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const Real c = 2.71828;

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Quadrilateral, { 6, 6 },
                                   Real(1) / Real(5));

    using FES = P0<Real, Mesh<Context::MPI>>;
    FES fes(mesh);

    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction  v(fes);

    Problem projection(u, v);
    projection = Integral(u, v) - Integral(RealFunction(c), v);

    PETSc::Solver::GMRES solver(projection);
    solver.solve();

    FES sh(mesh);
    GridFunction<FES, ::Vec> diff(sh);
    diff = Pow(u.getSolution() - RealFunction(c), 2);

    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0.0, 1e-10);
  }

  /**
   * @brief Single-rank P0 L2 projection matches the sequential result.
   */
  TEST(PETSc_MPI_P0, SolveL2Projection_SingleRank_MatchesSequential_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() != 1)
      GTEST_SKIP() << "This test is designed for exactly 1 MPI rank.";

    const Real c = 7.0;

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 6, 6 },
                                   Real(1) / Real(5));

    using FES = P0<Real, Mesh<Context::MPI>>;
    FES fes(mesh);

    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction  v(fes);

    Problem projection(u, v);
    projection = Integral(u, v) - Integral(RealFunction(c), v);

    PETSc::Solver::GMRES solver(projection);
    solver.solve();

    FES sh(mesh);
    GridFunction<FES, ::Vec> diff(sh);
    diff = Pow(u.getSolution() - RealFunction(c), 2);

    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0.0, 1e-10);
  }

  // =========================================================================
  // P0g — GridFunction projection
  // =========================================================================

  /**
   * @brief Project a constant onto distributed P0g and verify the projection
   *        by checking \int u_h d\Omega = c \cdot |\Omega|.
   *
   * The test uses a [0,1]×[0,1] mesh (area = 1), so \int u_h d\Omega = c.
   */
  TEST(PETSc_MPI_P0g, GridFunctionProjection_ConstantFunction_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const Real c = 4.2;

    Context::MPI ctx(*g_env, world);
    // Scale to [0,1]x[0,1] → area = 1.
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 6, 6 },
                                   Real(1) / Real(5));

    P0g<Real, Mesh<Context::MPI>> fes(mesh);
    GridFunction<decltype(fes), ::Vec> u(fes);
    u = RealFunction(c);

    // Verify projection via collective integral: \int u_h d\Omega = c * 1 = c
    const Real total = Integral(u).compute();
    EXPECT_NEAR(total, c, 1e-10);
  }

  // =========================================================================
  // P0g — Linear form assembly
  // =========================================================================

  /**
   * @brief \int 1\cdotv d\Omega with P0g has global size 1.
   */
  TEST(PETSc_MPI_P0g, LinearForm_GlobalVectorSizeIsOne_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 6, 6 });

    P0g<Real, Mesh<Context::MPI>> fes(mesh);
    PETSc::Variational::TestFunction v(fes);

    LinearForm lf(v);
    lf = Integral(RealFunction(1.0), v);
    lf.assemble();

    ::Vec b = lf.getVector();
    PetscInt globalSize = 0;
    VecGetSize(b, &globalSize);

    EXPECT_EQ(globalSize, 1);
  }

  /**
   * @brief \int 1\cdotv d\Omega with P0g: the single entry equals |\Omega|.
   *
   * For f = 1 the single P0g DOF entry accumulates contributions from all
   * cells: ∑_K |K| = |\Omega|.
   */
  TEST(PETSc_MPI_P0g, LinearForm_SingleEntry_EqualsDomainArea_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    // Scale to [0,1]x[0,1] → area = 1.
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 6, 6 },
                                   Real(1) / Real(5));

    P0g<Real, Mesh<Context::MPI>> fes(mesh);
    PETSc::Variational::TestFunction v(fes);

    LinearForm lf(v);
    lf = Integral(RealFunction(1.0), v);
    lf.assemble();

    ::Vec b = lf.getVector();
    PetscReal domainArea = 0.0;
    VecSum(b, &domainArea);

    EXPECT_NEAR(static_cast<Real>(domainArea), 1.0, 1e-10);
  }

  // =========================================================================
  // P0g — Bilinear form assembly
  // =========================================================================

  /**
   * @brief Assemble \int u\cdotv d\Omega with P0g: global matrix is 1×1.
   *
   * For scalar P0g the system is 1×1 with the single entry equal to |\Omega|.
   */
  TEST(PETSc_MPI_P0g, BilinearForm_MassMatrix_1x1_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 6, 6 });

    P0g<Real, Mesh<Context::MPI>> fes(mesh);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction  v(fes);

    BilinearForm bf(u, v);
    bf = Integral(u, v);
    bf.assemble();

    ::Mat A = bf.getOperator();
    PetscInt globalRows = 0;
    PetscInt globalCols = 0;
    MatGetSize(A, &globalRows, &globalCols);

    EXPECT_EQ(globalRows, 1);
    EXPECT_EQ(globalCols, 1);
  }

  /**
   * @brief P0g mass matrix entry on [0,1]x[0,1] equals the domain area (1.0).
   */
  TEST(PETSc_MPI_P0g, BilinearForm_MassMatrixEntry_EqualsDomainArea_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    // Scale to [0,1]x[0,1] → area = 1.
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 6, 6 },
                                   Real(1) / Real(5));

    P0g<Real, Mesh<Context::MPI>> fes(mesh);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction  v(fes);

    BilinearForm bf(u, v);
    bf = Integral(u, v);
    bf.assemble();

    ::Mat A = bf.getOperator();

    // Only rank 0 owns the entry (0,0).
    if (world.rank() == 0)
    {
      PetscInt    row = 0;
      PetscInt    col = 0;
      PetscScalar val = 0.0;
      MatGetValues(A, 1, &row, 1, &col, &val);
      EXPECT_NEAR(PetscRealPart(val), 1.0, 1e-10);
    }
    world.barrier();
  }

  // =========================================================================
  // P0g — Solve L2 projection
  // =========================================================================

  /**
   * @brief Solve the L2 projection of a constant onto P0g via PETSc GMRES.
   *
   * Problem: find c ∈ P0g s.t. \int c v d\Omega = \int const v d\Omega for all v ∈ P0g.
   * Solution: c = const.  L2 error = 0.
   */
  TEST(PETSc_MPI_P0g, SolveL2Projection_Constant_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const Real c = 3.0;

    Context::MPI ctx(*g_env, world);
    // Scale to [0,1]x[0,1] → area = 1.
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 6, 6 },
                                   Real(1) / Real(5));

    using FES = P0g<Real, Mesh<Context::MPI>>;
    FES fes(mesh);

    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction  v(fes);

    Problem projection(u, v);
    projection = Integral(u, v) - Integral(RealFunction(c), v);

    PETSc::Solver::GMRES solver(projection);
    solver.solve();

    FES sh(mesh);
    GridFunction<FES, ::Vec> diff(sh);
    diff = Pow(u.getSolution() - RealFunction(c), 2);

    // Integral(diff) = \int (c_h - c)^2 d\Omega.  For a constant rhs this must be 0.
    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0.0, 1e-10);
  }

  /**
   * @brief Solve the P0g L2 projection on a Quadrilateral mesh.
   */
  TEST(PETSc_MPI_P0g, SolveL2Projection_Constant_Quadrilateral)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    const Real c = 1.5;

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Quadrilateral, { 6, 6 },
                                   Real(1) / Real(5));

    using FES = P0g<Real, Mesh<Context::MPI>>;
    FES fes(mesh);

    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction  v(fes);

    Problem projection(u, v);
    projection = Integral(u, v) - Integral(RealFunction(c), v);

    PETSc::Solver::GMRES solver(projection);
    solver.solve();

    FES sh(mesh);
    GridFunction<FES, ::Vec> diff(sh);
    diff = Pow(u.getSolution() - RealFunction(c), 2);

    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0.0, 1e-10);
  }

  /**
   * @brief P0g L2 projection of f(x,y) = x + y on [0,1]x[0,1].
   *
   * The global average of f(x,y) = x + y over [0,1]^2 is 1/2 + 1/2 = 1.
   * The L2 error \int (c_h - f)^2 d\Omega = \int (1 - x - y)^2 d\Omega = 1/6.
   */
  TEST(PETSc_MPI_P0g, SolveL2Projection_LinearF_GlobalAverageIsOne_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    // Scale to [0,1]x[0,1] → area = 1.
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 12, 12 },
                                   Real(1) / Real(11));

    using FES = P0g<Real, Mesh<Context::MPI>>;
    FES fes(mesh);

    auto f = F::x + F::y;

    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction  v(fes);

    Problem projection(u, v);
    projection = Integral(u, v) - Integral(f, v);

    PETSc::Solver::GMRES solver(projection);
    solver.solve();

    // Verify that the solution value (global average of f) is close to 1.
    // The solution is stored as a P0g GridFunction; DOF 0 = global average.
    FES sh(mesh);
    GridFunction<FES, ::Vec> diff(sh);
    diff = Pow(u.getSolution() - RealFunction(1.0), 2);

    // \int (c_h - 1)^2 d\Omega should be 0 (c_h = 1 = average of f over [0,1]^2)
    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0.0, 1e-6);
  }

} // namespace Rodin::Tests::Unit::PETSc::MPI

// ---------------------------------------------------------------------------
// main() — initializes MPI + PETSc used by all tests.
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
