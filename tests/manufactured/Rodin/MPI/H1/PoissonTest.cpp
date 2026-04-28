/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file PoissonTest.cpp
 * @brief Distributed MPI H1 Poisson manufactured tests — all three mesh
 *        construction workflows.
 *
 * Parametrized over (geometry, K) with K ∈ {1, 2, 3, 4, 6}.
 *
 * Geometries:
 *   2D — Triangle, Quadrilateral
 *   3D — Tetrahedron, Hexahedron, Wedge
 *
 * Each GTest parametrization is invoked by CTest with
 * MPIEXEC_NUMPROC_FLAG set to 1, 2, 3, 4, and 6 ranks.
 *
 * ### Mesh construction workflows under test
 *
 * **Workflow 1 — Shard → Save → Load → Reconcile → Assemble → Solve**
 * Rank 0 builds a local uniform-grid mesh, partitions it with
 * `BalancedCompactPartitioner`, and calls `Sharder::scatter()`.  Each rank
 * saves its shard to an HDF5 file, reloads it, recomputes all connectivity,
 * and reconciles intermediate entities before assembly.
 *
 * **Workflow 2 — Construct → Shard → Distribute → Reconcile → Assemble → Solve**
 * Rank 0 builds a local uniform-grid mesh and calls `Sharder::distribute()`,
 * which performs shard + scatter + gather in one step without any I/O
 * round-trip.  Connectivity and reconciliation are applied to the resulting
 * distributed mesh.
 *
 * **Workflow 3 — UniformGrid → Reconcile → Assemble → Solve**
 * Every rank constructs its portion of the mesh independently via
 * `Mesh<Context::MPI>::UniformGrid()`, which builds ghost layers and
 * ownership information without any inter-rank communication beyond what the
 * constructor already performs.  Connectivity and reconciliation are then
 * applied before assembly.
 *
 * ### Manufactured solutions
 *
 * **2D:** @f$ u = \sin(\pi x)\sin(\pi y) @f$ on a 16×16 grid (h = 1/15).
 * RHS: @f$ f = 2\pi^2 \sin(\pi x)\sin(\pi y) @f$.
 * Tolerance: @ref RODIN_FUZZY_CONSTANT.
 *
 * **3D:** @f$ u = x(1-x)\,y(1-y)\,z(1-z) @f$ on a 5×5×5 grid (h = 1/4).
 * Zero Dirichlet data on all six faces.  A polynomial solution is used to
 * avoid the @f$ \pi^{K+1} @f$ norm blow-up on coarse meshes.
 * Tolerance: @f$ 10 \times @f$ @ref RODIN_FUZZY_CONSTANT.
 */

#include <cmath>
#include <limits>
#include <string>
#include <tuple>

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
#include <Rodin/MPI/Variational/H1.h>
#include <Rodin/PETSc.h>

#if defined(_WIN32)
#  include <process.h>
#  define RODIN_GETPID _getpid
#else
#  include <unistd.h>
#  define RODIN_GETPID getpid
#endif

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

// ---------------------------------------------------------------------------
// Global MPI handles (initialized in main())
// ---------------------------------------------------------------------------
static boost::mpi::environment*  g_env   = nullptr;
static boost::mpi::communicator* g_world = nullptr;

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------
namespace
{
  /**
   * @brief Builds a local mesh with the incidence data required by the
   *        Sharder for 2D and 3D meshes.
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
   * @brief Computes all connectivity needed for H1 assembly at any polynomial
   *        degree and reconciles all intermediate entity dimensions.
   *
   * This helper is shared by all three mesh construction workflows.  It
   * assumes that cell-vertex connectivity (@f$ D \to 0 @f$) is already
   * present (which is true after shard-gather, HDF5 reload, and
   * `Mesh<Context::MPI>::UniformGrid`).
   *
   * For 3D meshes the edge chain (@f$ D \to 1 @f$, @f$ D-1 \to 1 @f$,
   * @f$ 1 \to 0 @f$) is computed in addition to the facet chain, and both
   * `reconcile(1)` and `reconcile(2)` are called.
   */
  static void setupConnectivityAndReconcile(Mesh<Context::MPI>& mesh)
  {
    const size_t D = mesh.getDimension();

    // Core chain — valid for D = 2 and D = 3.
    mesh.getConnectivity().compute(D, D);
    mesh.getConnectivity().compute(D, 0);
    mesh.getConnectivity().compute(D, D - 1);
    mesh.getConnectivity().compute(D - 1, D);
    mesh.getConnectivity().compute(D - 1, 0);

    // Edge chain — needed for K ≥ 2 interior edge DOFs in 3D.
    if (D >= 3)
    {
      mesh.getConnectivity().compute(D, 1);
      mesh.getConnectivity().compute(D - 1, 1);
      mesh.getConnectivity().compute(1, 0);
    }

    // Reconcile all intermediate entity dimensions (1 in 2D; 1 and 2 in 3D).
    for (size_t d = 1; d < D; ++d)
      mesh.reconcile(d);
  }

  // ---------------------------------------------------------------------------
  // Workflow 1: Shard → Save MPI mesh → Load MPI mesh → Reconcile
  // ---------------------------------------------------------------------------

  /**
   * @brief Persists a distributed shard to HDF5, reloads it, and applies
   *        full connectivity and reconciliation.
   *
   * This is the core of **Workflow 1**: after sharding the mesh is persisted
   * to disk and reloaded to exercise the complete I/O round-trip before FE
   * assembly.
   */
  static Mesh<Context::MPI> reloadShard(
      const Context::MPI& ctx,
      Mesh<Context::MPI>&& mesh)
  {
    const auto& comm = ctx.getCommunicator();
    const size_t D = mesh.getDimension();

    const boost::filesystem::path rankFile =
        boost::filesystem::temp_directory_path()
        / ("rodin_mpi_h1_poisson_full_"
           + std::to_string(RODIN_GETPID())
           + "_r" + std::to_string(comm.rank()) + ".h5");

    {
      IO::MeshPrinter<IO::FileFormat::HDF5, Context::MPI> printer(mesh);
      printer.print(rankFile);
    }
    comm.barrier();

    Mesh<Context::MPI> r(ctx);
    r.load(rankFile, IO::FileFormat::HDF5);

    EXPECT_EQ(r.getShard().getCellCount(), mesh.getShard().getCellCount())
      << "Rank " << comm.rank() << ": cell count mismatch after HDF5 reload.";
    EXPECT_EQ(r.getShard().getVertexCount(), mesh.getShard().getVertexCount())
      << "Rank " << comm.rank() << ": vertex count mismatch after HDF5 reload.";

    comm.barrier();
    boost::filesystem::remove(rankFile);

    setupConnectivityAndReconcile(r);
    comm.barrier();
    return r;
  }

  /**
   * @brief **Workflow 1** — Shard → Save MPI mesh → Load MPI mesh →
   *        Reconcile.
   *
   * Rank 0 builds a uniform-grid mesh, partitions it with
   * `BalancedCompactPartitioner`, scatters shards, each rank saves to HDF5
   * and reloads before assembly.
   */
  static Mesh<Context::MPI> meshViaShardSaveLoad(
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
    return reloadShard(ctx, std::move(mesh));
  }

  // ---------------------------------------------------------------------------
  // Workflow 2: Construct → Shard → Distribute → Reconcile
  // ---------------------------------------------------------------------------

  /**
   * @brief **Workflow 2** — Construct → Shard → Distribute → Reconcile.
   *
   * Rank 0 builds a local uniform-grid mesh, partitions it, shards and
   * scatters the shards to all ranks.  All ranks gather their local shard
   * directly — no HDF5 I/O round-trip occurs.  Connectivity and
   * reconciliation are then applied to the distributed mesh.
   *
   * This workflow differs from Workflow 1 only in the absence of the
   * HDF5 save/load step.
   */
  static Mesh<Context::MPI> meshViaSharderDistribute(
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
    setupConnectivityAndReconcile(mesh);
    comm.barrier();
    return mesh;
  }

  // ---------------------------------------------------------------------------
  // Workflow 3: UniformGrid → Reconcile
  // ---------------------------------------------------------------------------

  /**
   * @brief **Workflow 3** — UniformGrid → Reconcile → Assemble → Solve.
   *
   * All ranks construct their portion of the mesh independently via
   * `Mesh<Context::MPI>::UniformGrid()`, which builds ghost layers and
   * ownership information without a rank-0-centric sharding step.
   * Connectivity and reconciliation are applied before assembly.
   */
  static Mesh<Context::MPI> meshViaUniformGrid(
      const Context::MPI& ctx,
      Polytope::Type type,
      std::initializer_list<size_t> shape)
  {
    auto mesh = Mesh<Context::MPI>::UniformGrid(ctx, type, shape);
    mesh.scale(1.0 / static_cast<Real>(*shape.begin() - 1));
    setupConnectivityAndReconcile(mesh);
    ctx.getCommunicator().barrier();
    return mesh;
  }

  // -------------------------------------------------------------------------
  // Per-K Poisson solvers
  // -------------------------------------------------------------------------

  /**
   * @brief Solves the 2D Poisson problem with H1 of degree K and returns the
   *        globally reduced @f$ \int_\Omega (u - u_h)^2 @f$.
   *
   * Manufactured solution: @f$ u = \sin(\pi x)\sin(\pi y) @f$,
   * right-hand side: @f$ f = 2\pi^2 \sin(\pi x)\sin(\pi y) @f$.
   */
  template <size_t K>
  Real computeError2D(
      boost::mpi::communicator& world,
      Geometry::Mesh<Context::MPI>& mesh)
  {
    namespace PETSc = ::Rodin::PETSc;
    using FES = Variational::H1<K, Real, Geometry::Mesh<Context::MPI>>;

    constexpr auto order = std::integral_constant<size_t, K>{};
    const auto pi = Math::Constants::pi();

    FES vh(order, mesh);
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

    FES sh(order, mesh);
    GridFunction<FES, ::Vec> diff(sh);
    diff = Pow(u.getSolution() - solution, 2);

    const Real localError = Integral(diff).compute();
    Real globalError = 0;
    boost::mpi::all_reduce(world, localError, globalError, std::plus<Real>());
    return globalError;
  }

  /**
   * @brief Solves the 3D Poisson problem with H1 of degree K and returns the
   *        globally reduced @f$ \int_\Omega (u - u_h)^2 @f$.
   *
   * Manufactured solution: @f$ u = x(1-x)\,y(1-y)\,z(1-z) @f$.
   * Right-hand side:
   * @f$ f = 2\bigl[y(1-y)z(1-z) + x(1-x)z(1-z) + x(1-x)y(1-y)\bigr] @f$.
   * The solution satisfies homogeneous Dirichlet conditions on all six faces.
   */
  template <size_t K>
  Real computeError3D(
      boost::mpi::communicator& world,
      Geometry::Mesh<Context::MPI>& mesh)
  {
    namespace PETSc = ::Rodin::PETSc;
    using FES = Variational::H1<K, Real, Geometry::Mesh<Context::MPI>>;

    constexpr auto order = std::integral_constant<size_t, K>{};

    FES vh(order, mesh);

    auto f = 2 * F::y * (1 - F::y) * F::z * (1 - F::z)
           + 2 * F::x * (1 - F::x) * F::z * (1 - F::z)
           + 2 * F::x * (1 - F::x) * F::y * (1 - F::y);

    PETSc::Variational::TrialFunction u(vh);
    PETSc::Variational::TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());

    PETSc::Solver::CG solver(poisson);
    solver.solve();

    auto solution = F::x * (1 - F::x) * F::y * (1 - F::y) * F::z * (1 - F::z);

    FES sh(order, mesh);
    GridFunction<FES, ::Vec> diff(sh);
    diff = Pow(u.getSolution() - solution, 2);

    const Real localError = Integral(diff).compute();
    Real globalError = 0;
    boost::mpi::all_reduce(world, localError, globalError, std::plus<Real>());
    return globalError;
  }

  // -------------------------------------------------------------------------
  // Runtime K dispatchers
  // -------------------------------------------------------------------------

  Real run2D(
      size_t K,
      boost::mpi::communicator& world,
      Geometry::Mesh<Context::MPI>& mesh)
  {
    switch (K)
    {
      case 1: return computeError2D<1>(world, mesh);
      case 2: return computeError2D<2>(world, mesh);
      case 3: return computeError2D<3>(world, mesh);
      case 4: return computeError2D<4>(world, mesh);
      case 6: return computeError2D<6>(world, mesh);
      default:
        ADD_FAILURE() << "run2D: unsupported K=" << K;
        return std::numeric_limits<Real>::max();
    }
  }

  Real run3D(
      size_t K,
      boost::mpi::communicator& world,
      Geometry::Mesh<Context::MPI>& mesh)
  {
    switch (K)
    {
      case 1: return computeError3D<1>(world, mesh);
      case 2: return computeError3D<2>(world, mesh);
      case 3: return computeError3D<3>(world, mesh);
      case 4: return computeError3D<4>(world, mesh);
      case 6: return computeError3D<6>(world, mesh);
      default:
        ADD_FAILURE() << "run3D: unsupported K=" << K;
        return std::numeric_limits<Real>::max();
    }
  }

  // -------------------------------------------------------------------------
  // Human-readable GTest parameter names
  // -------------------------------------------------------------------------
  struct ParamName
  {
    std::string operator()(
        const testing::TestParamInfo<std::tuple<Polytope::Type, size_t>>& info)
        const
    {
      const auto& [geom, K] = info.param;
      std::string g;
      switch (geom)
      {
        case Polytope::Type::Triangle:      g = "Triangle";      break;
        case Polytope::Type::Quadrilateral: g = "Quadrilateral"; break;
        case Polytope::Type::Tetrahedron:   g = "Tetrahedron";   break;
        case Polytope::Type::Hexahedron:    g = "Hexahedron";    break;
        case Polytope::Type::Wedge:         g = "Wedge";         break;
        default:                            g = "Unknown";        break;
      }
      return g + "_K" + std::to_string(K);
    }
  };
} // anonymous namespace

// ---------------------------------------------------------------------------
// Test fixtures and parametrizations
// ---------------------------------------------------------------------------
namespace Rodin::Tests::Manufactured::MPI::H1Poisson
{
  using Param = std::tuple<Polytope::Type, size_t>;

  // =========================================================================
  // Workflow 1 — Shard → Save → Load → Reconcile → Assemble → Solve
  // =========================================================================

  class H1Poisson2D_Workflow1 : public ::testing::TestWithParam<Param> {};

  /**
   * @brief Workflow 1 — 2D Poisson with HDF5 mesh I/O round-trip.
   *
   * Mesh: 16×16 uniform grid, shard → save HDF5 → load HDF5 → reconcile.
   * Solution: sin(πx)sin(πy).
   * Tolerance: @ref RODIN_FUZZY_CONSTANT.
   */
  TEST_P(H1Poisson2D_Workflow1, SimpleSine)
  {
    const auto& [geom, K] = GetParam();
    auto& world = *g_world;

    Context::MPI ctx(*g_env, world);
    auto mesh = meshViaShardSaveLoad(ctx, geom, { 16, 16 });

    const Real globalError = run2D(K, world, mesh);
    EXPECT_NEAR(globalError, 0, RODIN_FUZZY_CONSTANT)
        << "np=" << world.size()
        << " geometry=" << static_cast<int>(geom)
        << " K=" << K;
  }

  INSTANTIATE_TEST_SUITE_P(
    AllGeometriesAndDegrees,
    H1Poisson2D_Workflow1,
    testing::Combine(
      testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral),
      testing::Values<size_t>(1, 2, 3, 4, 6)
    ),
    ParamName{}
  );

  class H1Poisson3D_Workflow1 : public ::testing::TestWithParam<Param> {};

  /**
   * @brief Workflow 1 — 3D Poisson with HDF5 mesh I/O round-trip.
   *
   * Mesh: 5×5×5 uniform grid, shard → save HDF5 → load HDF5 → reconcile.
   * Solution: x(1−x)y(1−y)z(1−z).
   * Tolerance: 10 × @ref RODIN_FUZZY_CONSTANT.
   */
  TEST_P(H1Poisson3D_Workflow1, PolynomialBC)
  {
    const auto& [geom, K] = GetParam();
    auto& world = *g_world;

    Context::MPI ctx(*g_env, world);
    auto mesh = meshViaShardSaveLoad(ctx, geom, { 5, 5, 5 });

    const Real globalError = run3D(K, world, mesh);
    EXPECT_NEAR(globalError, 0, 10 * RODIN_FUZZY_CONSTANT)
        << "np=" << world.size()
        << " geometry=" << static_cast<int>(geom)
        << " K=" << K;
  }

  INSTANTIATE_TEST_SUITE_P(
    AllGeometriesAndDegrees,
    H1Poisson3D_Workflow1,
    testing::Combine(
      testing::Values(
        Polytope::Type::Tetrahedron,
        Polytope::Type::Hexahedron,
        Polytope::Type::Wedge),
      testing::Values<size_t>(1, 2, 3, 4, 6)
    ),
    ParamName{}
  );

  // =========================================================================
  // Workflow 2 — Construct → Shard → Distribute → Reconcile → Assemble → Solve
  // =========================================================================

  class H1Poisson2D_Workflow2 : public ::testing::TestWithParam<Param> {};

  /**
   * @brief Workflow 2 — 2D Poisson with direct shard-gather (no I/O).
   *
   * Mesh: 16×16 uniform grid, shard → scatter → gather → reconcile (no HDF5).
   * Solution: sin(πx)sin(πy).
   * Tolerance: @ref RODIN_FUZZY_CONSTANT.
   */
  TEST_P(H1Poisson2D_Workflow2, SimpleSine)
  {
    const auto& [geom, K] = GetParam();
    auto& world = *g_world;

    Context::MPI ctx(*g_env, world);
    auto mesh = meshViaSharderDistribute(ctx, geom, { 16, 16 });

    const Real globalError = run2D(K, world, mesh);
    EXPECT_NEAR(globalError, 0, RODIN_FUZZY_CONSTANT)
        << "np=" << world.size()
        << " geometry=" << static_cast<int>(geom)
        << " K=" << K;
  }

  INSTANTIATE_TEST_SUITE_P(
    AllGeometriesAndDegrees,
    H1Poisson2D_Workflow2,
    testing::Combine(
      testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral),
      testing::Values<size_t>(1, 2, 3, 4, 6)
    ),
    ParamName{}
  );

  class H1Poisson3D_Workflow2 : public ::testing::TestWithParam<Param> {};

  /**
   * @brief Workflow 2 — 3D Poisson with direct shard-gather (no I/O).
   *
   * Mesh: 5×5×5 uniform grid, shard → scatter → gather → reconcile (no HDF5).
   * Solution: x(1−x)y(1−y)z(1−z).
   * Tolerance: 10 × @ref RODIN_FUZZY_CONSTANT.
   */
  TEST_P(H1Poisson3D_Workflow2, PolynomialBC)
  {
    const auto& [geom, K] = GetParam();
    auto& world = *g_world;

    Context::MPI ctx(*g_env, world);
    auto mesh = meshViaSharderDistribute(ctx, geom, { 5, 5, 5 });

    const Real globalError = run3D(K, world, mesh);
    EXPECT_NEAR(globalError, 0, 10 * RODIN_FUZZY_CONSTANT)
        << "np=" << world.size()
        << " geometry=" << static_cast<int>(geom)
        << " K=" << K;
  }

  INSTANTIATE_TEST_SUITE_P(
    AllGeometriesAndDegrees,
    H1Poisson3D_Workflow2,
    testing::Combine(
      testing::Values(
        Polytope::Type::Tetrahedron,
        Polytope::Type::Hexahedron,
        Polytope::Type::Wedge),
      testing::Values<size_t>(1, 2, 3, 4, 6)
    ),
    ParamName{}
  );

  // =========================================================================
  // Workflow 3 — UniformGrid → Reconcile → Assemble → Solve
  // =========================================================================

  class H1Poisson2D_Workflow3 : public ::testing::TestWithParam<Param> {};

  /**
   * @brief Workflow 3 — 2D Poisson built via distributed `UniformGrid`.
   *
   * Mesh: 16×16 grid constructed directly on all ranks via
   * `Mesh<Context::MPI>::UniformGrid()`, then reconciled.
   * Solution: sin(πx)sin(πy).
   * Tolerance: @ref RODIN_FUZZY_CONSTANT.
   */
  TEST_P(H1Poisson2D_Workflow3, SimpleSine)
  {
    const auto& [geom, K] = GetParam();
    auto& world = *g_world;

    Context::MPI ctx(*g_env, world);
    auto mesh = meshViaUniformGrid(ctx, geom, { 16, 16 });

    const Real globalError = run2D(K, world, mesh);
    EXPECT_NEAR(globalError, 0, RODIN_FUZZY_CONSTANT)
        << "np=" << world.size()
        << " geometry=" << static_cast<int>(geom)
        << " K=" << K;
  }

  INSTANTIATE_TEST_SUITE_P(
    AllGeometriesAndDegrees,
    H1Poisson2D_Workflow3,
    testing::Combine(
      testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral),
      testing::Values<size_t>(1, 2, 3, 4, 6)
    ),
    ParamName{}
  );

  class H1Poisson3D_Workflow3 : public ::testing::TestWithParam<Param> {};

  /**
   * @brief Workflow 3 — 3D Poisson built via distributed `UniformGrid`.
   *
   * Mesh: 5×5×5 grid constructed directly on all ranks via
   * `Mesh<Context::MPI>::UniformGrid()`, then reconciled.
   * Solution: x(1−x)y(1−y)z(1−z).
   * Tolerance: 10 × @ref RODIN_FUZZY_CONSTANT.
   */
  TEST_P(H1Poisson3D_Workflow3, PolynomialBC)
  {
    const auto& [geom, K] = GetParam();
    auto& world = *g_world;

    Context::MPI ctx(*g_env, world);
    auto mesh = meshViaUniformGrid(ctx, geom, { 5, 5, 5 });

    const Real globalError = run3D(K, world, mesh);
    EXPECT_NEAR(globalError, 0, 10 * RODIN_FUZZY_CONSTANT)
        << "np=" << world.size()
        << " geometry=" << static_cast<int>(geom)
        << " K=" << K;
  }

  INSTANTIATE_TEST_SUITE_P(
    AllGeometriesAndDegrees,
    H1Poisson3D_Workflow3,
    testing::Combine(
      testing::Values(
        Polytope::Type::Tetrahedron,
        Polytope::Type::Hexahedron,
        Polytope::Type::Wedge),
      testing::Values<size_t>(1, 2, 3, 4, 6)
    ),
    ParamName{}
  );
} // namespace Rodin::Tests::Manufactured::MPI::H1Poisson

// ---------------------------------------------------------------------------
// main() — initializes MPI and PETSc before running all tests.
// ---------------------------------------------------------------------------
int main(int argc, char** argv)
{
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  g_env   = &env;
  g_world = &world;

  [[maybe_unused]] PetscErrorCode ierr =
      PetscInitialize(&argc, &argv, nullptr, nullptr);
  assert(ierr == PETSC_SUCCESS);

  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();

  PetscFinalize();
  return result;
}
