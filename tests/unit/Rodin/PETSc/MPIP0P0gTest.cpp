/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file MPIP0P0gTest.cpp
 * @brief Distributed PETSc tests for P0 and P0g spaces on MPI meshes.
 *
 * Mesh distribution, interpolation, assembly and mass-projection solves use
 * the configured PETSc scalar field @f$\mathbb F\in\{\mathbb R,\mathbb C\}@f$.
 * Legacy projection data remain real-valued, including when embedded in complex
 * storage. Error integrands use @f$|u_h-c|^2=(u_h-c)\overline{(u_h-c)}@f$;
 * scalar comparisons retain both real and imaginary components.
 *
 * ## Piecewise constants
 *
 * Mass projection onto @f$P_0@f$ is defined by
 * @f[
 * \int_\Omega u_h\overline v\,dx=\int_\Omega f\overline v\,dx
 * \quad\text{for every }v\in P_0.
 * @f]
 * The diagonal mass matrix has @f$M_{ii}=|K_i|@f$, so each coefficient is the
 * cell average @f$u_h|_{K_i}=|K_i|^{-1}\int_{K_i}f\,dx@f$.
 *
 * ## Global constants
 *
 * Projection onto @f$P_{0g}@f$ yields the domain average
 * @f$c=|\Omega|^{-1}\int_\Omega f\,dx@f$. Constants are reproduced exactly
 * up to the stated algebraic tolerance. Interpolation through a DOF functional
 * is tested separately and is not asserted to compute an average.
 */

#include <cmath>
#include <string>
#include <type_traits>

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
  /// @brief Scalar data with the configured PETSc coefficient field.
  using PETScScalarFunction = std::conditional_t<std::is_same_v<PetscScalar, Complex>,
    ComplexFunction<PetscScalar>, RealFunction<PetscScalar>>;

  /// @brief Leaves rank zero empty while retaining the existing P0g DOF owner.
  class EmptyRootPartitioner final : public Partitioner
  {
    public:
      EmptyRootPartitioner(const Mesh<Context::Local>& mesh, size_t count)
        : m_partitioner(mesh),
          m_count(count)
      {
        m_partitioner.partition(count > 1 ? count - 1 : count);
      }

      const Mesh<Context::Local>& getMesh() const override
      {
        return m_partitioner.getMesh();
      }

      void partition(size_t count, size_t dimension) override
      {
        m_count = count;
        m_partitioner.partition(count > 1 ? count - 1 : count, dimension);
      }

      size_t getCount() const override
      {
        return m_count;
      }

      size_t getPartition(Index index) const override
      {
        return m_partitioner.getPartition(index) + (m_count > 1);
      }

    private:
      BalancedCompactPartitioner m_partitioner;
      size_t m_count;
  };

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

  /**
   * @brief Interpolation reaches the unique global-constant DOF owner even
   *        when that owner has no entities on which to evaluate the function.
   *
   * The two-cell segment is assigned only to nonzero ranks. With four ranks,
   * both the owner and one ghost holder have empty shards; the untouched ghost
   * must not replace the interpolated coefficient. Scalar and vector constants
   * are checked exactly on every copy. A subsequent empty selection must leave
   * both fields unchanged.
   */
  TEST(PETSc_MPI_P0g, InterpolationWithEmptyOwner)
  {
    const auto& world = *g_world;
    Context::MPI ctx(*g_env, world);
    Sharder<Context::MPI> sharder(ctx);
    if (world.rank() == 0)
    {
      auto local = makeShardableMesh(Polytope::Type::Segment, {3});
      EmptyRootPartitioner partitioner(local, static_cast<size_t>(world.size()));
      sharder.shard(partitioner);
      sharder.scatter(0);
    }
    auto mesh = sharder.gather(0);
    if (world.size() > 1 && world.rank() == 0)
    {
      EXPECT_EQ(mesh.getShard().getVertexCount(), 0);
    }

    P0g<PetscScalar, decltype(mesh)> scalarSpace(mesh);
    P0g<Math::SpatialVector<PetscScalar>, decltype(mesh)> vectorSpace(mesh, size_t(2));
    PETSc::Variational::GridFunction scalar(scalarSpace);
    PETSc::Variational::GridFunction vector(vectorSpace);
    scalar = PetscScalar(-7);
    vector = PetscScalar(-7);
    scalar = RealFunction(3);
    vector = VectorFunction(RealFunction(3), RealFunction(5));

    const auto& scalarRead = scalar;
    const auto& vectorRead = vector;
    EXPECT_EQ(scalarRead[0], PetscScalar(3));
    EXPECT_EQ(vectorRead[0], PetscScalar(3));
    EXPECT_EQ(vectorRead[1], PetscScalar(5));
    scalarRead.flush();
    vectorRead.flush();

    scalar.project(Region::Cells, RealFunction(11), [](const auto&) { return false; });
    vector.project(Region::Cells, VectorFunction(RealFunction(11), RealFunction(13)),
      [](const auto&) { return false; });
    EXPECT_EQ(scalarRead[0], PetscScalar(3));
    EXPECT_EQ(vectorRead[0], PetscScalar(3));
    EXPECT_EQ(vectorRead[1], PetscScalar(5));
    scalarRead.flush();
    vectorRead.flush();
  }

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

    P0<PetscScalar, Mesh<Context::MPI>> fes(mesh);
    GridFunction<decltype(fes), ::Vec> u(fes);
    u = RealFunction(c);

    Index begin = 0;
    Index end   = 0;
    fes.getOwnershipRange(begin, end);

    // Acquire read access to the PETSc vector.
    for (Index i = begin; i < end; ++i)
    {
      const PetscScalar v = u[i];
      EXPECT_NEAR(std::abs(v - PetscScalar(c)), 0, 1e-10);
    }
    u.flush();
  }

  /**
   * @brief The global integral of the projected constant function equals
   *        the constant times the domain area.
   *
   * @f$ \int_\Omega u_h\,d\Omega = c \cdot |\Omega| @f$.
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

    P0<PetscScalar, Mesh<Context::MPI>> fes(mesh);
    GridFunction<decltype(fes), ::Vec> u(fes);
    u = RealFunction(c);

    // Integral of the grid function over the domain
    const auto total = Integral(u).compute();

    // Expected: c * |\Omega| = c * 1 = c
    EXPECT_NEAR(std::abs(total - PetscScalar(c)), 0., 1e-9);
  }

  // =========================================================================
  // P0 — Linear form assembly
  // =========================================================================

  /**
   * @brief Assemble @f$ \int_\Omega 1\cdot v\,d\Omega @f$ on distributed P0 and verify the global
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

    P0<PetscScalar, Mesh<Context::MPI>> fes(mesh);
    PETSc::Variational::TestFunction v(fes);

    LinearForm lf(v);
    lf = Integral(PETScScalarFunction(1.0), v);
    lf.assemble();

    ::Vec b = lf.getVector();
    PetscInt globalSize = 0;
    VecGetSize(b, &globalSize);

    // getCellCount() is collective — must be called on all ranks.
    const size_t globalCells = mesh.getCellCount();
    EXPECT_EQ(static_cast<size_t>(globalSize), globalCells);
  }

  /**
   * @brief Sum of @f$ \int_\Omega 1\cdot v\,d\Omega @f$ over all P0 test functions equals domain area.
   *
   * For f = 1 and P0, each entry b_K = |K| (cell measure).
   * @f$ \sum_K b_K = |\Omega| @f$.
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

    P0<PetscScalar, Mesh<Context::MPI>> fes(mesh);
    PETSc::Variational::TestFunction v(fes);

    LinearForm lf(v);
    lf = Integral(PETScScalarFunction(1.0), v);
    lf.assemble();

    ::Vec b = lf.getVector();
    PetscScalar sum = 0.0;
    VecSum(b, &sum);

    EXPECT_NEAR(std::abs(sum - PetscScalar(1)), 0.0, 1e-10);
  }

  // =========================================================================
  // P0 — Bilinear form assembly
  // =========================================================================

  /**
   * @brief Assemble the P0 mass matrix @f$ \int_\Omega u\cdot v\,d\Omega @f$ and check global dimensions.
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

    P0<PetscScalar, Mesh<Context::MPI>> fes(mesh);
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
   * Problem: find @f$ u_h \in P0 @f$ such that @f$ \int_\Omega u_h v\,d\Omega = \int_\Omega c v\,d\Omega @f$ for all @f$ v \in P0 @f$.
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

    using FES = P0<PetscScalar, Mesh<Context::MPI>>;
    FES fes(mesh);

    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction  v(fes);

    Problem projection(u, v);
    projection = Integral(u, v) - Integral(PETScScalarFunction(c), v);

    PETSc::Solver::GMRES solver(projection);
    solver.solve();

    // The Hermitian square remains nonnegative for complex coefficient storage.
    FES sh(mesh);
    GridFunction<FES, ::Vec> diff(sh);
    diff = Dot(
      u.getSolution() - PETScScalarFunction(c), u.getSolution() - PETScScalarFunction(c));

    const auto error = Integral(diff).compute();
    EXPECT_NEAR(std::abs(error), 0.0, 1e-10);
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

    using FES = P0<PetscScalar, Mesh<Context::MPI>>;
    FES fes(mesh);

    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction  v(fes);

    Problem projection(u, v);
    projection = Integral(u, v) - Integral(PETScScalarFunction(c), v);

    PETSc::Solver::GMRES solver(projection);
    solver.solve();

    FES sh(mesh);
    GridFunction<FES, ::Vec> diff(sh);
    diff = Dot(
      u.getSolution() - PETScScalarFunction(c), u.getSolution() - PETScScalarFunction(c));

    const auto error = Integral(diff).compute();
    EXPECT_NEAR(std::abs(error), 0.0, 1e-10);
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

    using FES = P0<PetscScalar, Mesh<Context::MPI>>;
    FES fes(mesh);

    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction  v(fes);

    Problem projection(u, v);
    projection = Integral(u, v) - Integral(PETScScalarFunction(c), v);

    PETSc::Solver::GMRES solver(projection);
    solver.solve();

    FES sh(mesh);
    GridFunction<FES, ::Vec> diff(sh);
    diff = Dot(
      u.getSolution() - PETScScalarFunction(c), u.getSolution() - PETScScalarFunction(c));

    const auto error = Integral(diff).compute();
    EXPECT_NEAR(std::abs(error), 0.0, 1e-10);
  }

  // =========================================================================
  // P0g — GridFunction projection
  // =========================================================================

  /**
   * @brief Project a constant onto distributed P0g and verify the projection
   *        by checking @f$ \int_\Omega u_h\,d\Omega = c \cdot |\Omega| @f$.
   *
   * The test uses a @f$ [0,1]\times[0,1] @f$ mesh with area 1, so @f$ \int_\Omega u_h\,d\Omega = c @f$.
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

    P0g<PetscScalar, Mesh<Context::MPI>> fes(mesh);
    GridFunction<decltype(fes), ::Vec> u(fes);
    u = RealFunction(c);

    // Verify projection via collective integral: \int u_h d\Omega = c * 1 = c
    const auto total = Integral(u).compute();
    EXPECT_NEAR(std::abs(total - PetscScalar(c)), 0., 1e-10);
  }

  // =========================================================================
  // P0g — Linear form assembly
  // =========================================================================

  /**
   * @brief @f$ \int_\Omega 1\cdot v\,d\Omega @f$ with P0g has global size 1.
   */
  TEST(PETSc_MPI_P0g, LinearForm_GlobalVectorSizeIsOne_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 6, 6 });

    P0g<PetscScalar, Mesh<Context::MPI>> fes(mesh);
    PETSc::Variational::TestFunction v(fes);

    LinearForm lf(v);
    lf = Integral(PETScScalarFunction(1.0), v);
    lf.assemble();

    ::Vec b = lf.getVector();
    PetscInt globalSize = 0;
    VecGetSize(b, &globalSize);

    EXPECT_EQ(globalSize, 1);
  }

  /**
   * @brief @f$ \int_\Omega 1\cdot v\,d\Omega @f$ with P0g has a single entry equal to @f$ |\Omega| @f$.
   *
   * For @f$ f = 1 @f$ the single P0g DOF entry accumulates contributions
   * from all cells: @f$ \sum_K |K| = |\Omega| @f$.
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

    P0g<PetscScalar, Mesh<Context::MPI>> fes(mesh);
    PETSc::Variational::TestFunction v(fes);

    LinearForm lf(v);
    lf = Integral(PETScScalarFunction(1.0), v);
    lf.assemble();

    ::Vec b = lf.getVector();
    PetscScalar domainArea = 0.0;
    VecSum(b, &domainArea);

    EXPECT_NEAR(std::abs(domainArea - PetscScalar(1)), 0.0, 1e-10);
  }

  // =========================================================================
  // P0g — Bilinear form assembly
  // =========================================================================

  /**
   * @brief Assemble @f$ \int_\Omega u\cdot v\,d\Omega @f$ with P0g: global matrix is 1 by 1.
   *
   * For scalar P0g the system is 1 by 1 with the single entry equal to @f$ |\Omega| @f$.
   */
  TEST(PETSc_MPI_P0g, BilinearForm_MassMatrix_1x1_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, { 6, 6 });

    P0g<PetscScalar, Mesh<Context::MPI>> fes(mesh);
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

    P0g<PetscScalar, Mesh<Context::MPI>> fes(mesh);
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
      EXPECT_NEAR(std::abs(val - PetscScalar(1)), 0, 1e-10);
    }
    world.barrier();
  }

  // =========================================================================
  // P0g — Solve L2 projection
  // =========================================================================

  /**
   * @brief Solve the L2 projection of a constant onto P0g via PETSc GMRES.
   *
   * Problem: find @f$ c \in P0g @f$ such that @f$ \int_\Omega c v\,d\Omega = \int_\Omega \mathrm{const}\, v\,d\Omega @f$ for all @f$ v \in P0g @f$.
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

    using FES = P0g<PetscScalar, Mesh<Context::MPI>>;
    FES fes(mesh);

    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction  v(fes);

    Problem projection(u, v);
    projection = Integral(u, v) - Integral(PETScScalarFunction(c), v);

    PETSc::Solver::GMRES solver(projection);
    solver.solve();

    FES sh(mesh);
    GridFunction<FES, ::Vec> diff(sh);
    diff = Dot(
      u.getSolution() - PETScScalarFunction(c), u.getSolution() - PETScScalarFunction(c));

    // A represented constant has zero squared L2 error.
    const auto error = Integral(diff).compute();
    EXPECT_NEAR(std::abs(error), 0.0, 1e-10);
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

    using FES = P0g<PetscScalar, Mesh<Context::MPI>>;
    FES fes(mesh);

    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction  v(fes);

    Problem projection(u, v);
    projection = Integral(u, v) - Integral(PETScScalarFunction(c), v);

    PETSc::Solver::GMRES solver(projection);
    solver.solve();

    FES sh(mesh);
    GridFunction<FES, ::Vec> diff(sh);
    diff = Dot(
      u.getSolution() - PETScScalarFunction(c), u.getSolution() - PETScScalarFunction(c));

    const auto error = Integral(diff).compute();
    EXPECT_NEAR(std::abs(error), 0.0, 1e-10);
  }

  /**
   * @brief P0g L2 projection of f(x,y) = x + y on [0,1]x[0,1].
   *
   * The global average of f(x,y) = x + y over [0,1]^2 is 1/2 + 1/2 = 1.
   * The L2 error is @f$ \int_\Omega (c_h - f)^2\,d\Omega = \int_\Omega (1 - x - y)^2\,d\Omega = 1/6 @f$.
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

    using FES = P0g<PetscScalar, Mesh<Context::MPI>>;
    FES fes(mesh);

    auto f = PETScScalarFunction(1) * (F::x + F::y);

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
    diff = Dot(u.getSolution() - PETScScalarFunction(1.0),
      u.getSolution() - PETScScalarFunction(1.0));

    // The exact domain average of the manufactured source is one.
    const auto error = Integral(diff).compute();
    EXPECT_NEAR(std::abs(error), 0.0, 1e-6);
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
