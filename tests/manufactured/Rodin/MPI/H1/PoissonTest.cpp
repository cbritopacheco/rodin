/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file PoissonTest.cpp
 * @brief Distributed MPI H1 Poisson manufactured tests.
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
 * ### 2D tests
 * Manufactured solution @f$ u = \sin(\pi x)\sin(\pi y) @f$ on a 16×16
 * uniform grid (h = 1/15).  The global @f$ L^2 @f$ error
 * @f$ \int_\Omega (u - u_h)^2 @f$ must stay below
 * @ref RODIN_FUZZY_CONSTANT.
 *
 * ### 3D tests
 * Manufactured solution @f$ u = x(1-x)\,y(1-y)\,z(1-z) @f$ on a 5×5×5
 * uniform grid (h = 1/4, 4³ = 64 cells).  This degree-6 polynomial has
 * @f$ u = 0 @f$ on all six faces of @f$ [0,1]^3 @f$, so the Dirichlet
 * data is identically zero.  Using a polynomial avoids the large
 * @f$ \pi^{K+1} @f$ norm factors that make trigonometric solutions
 * inaccurate on coarse 3D meshes for K=1,2.  For K ≥ 6 the solution is
 * exactly representable; for smaller K the interpolation error is still
 * well within @f$ 10 \times @f$ @ref RODIN_FUZZY_CONSTANT.
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
   * @brief Persists a distributed shard to HDF5, reloads it, and recomputes
   *        all connectivity needed for H1 assembly at any polynomial degree.
   *
   * For 3D meshes the edge chain (D→1, D−1→1, 1→0) is computed in addition
   * to the facet chain, and both reconcile(1) and reconcile(2) are called.
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

    // Core connectivity — valid for D = 2 and D = 3.
    r.getConnectivity().compute(D, D);
    r.getConnectivity().compute(D, 0);
    r.getConnectivity().compute(D, D - 1);
    r.getConnectivity().compute(D - 1, D);
    r.getConnectivity().compute(D - 1, 0);

    // Edge chain — needed for K ≥ 2 interior edge DOFs in 3D.
    if (D >= 3)
    {
      r.getConnectivity().compute(D, 1);
      r.getConnectivity().compute(D - 1, 1);
      r.getConnectivity().compute(1, 0);
    }

    // Reconcile all intermediate entity dimensions (1 in 2D; 1 and 2 in 3D).
    for (size_t d = 1; d < D; ++d)
      r.reconcile(d);

    comm.barrier();
    boost::filesystem::remove(rankFile);
    return r;
  }

  /**
   * @brief Partitions and distributes a uniform-grid mesh from rank 0.
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
    return reloadShard(ctx, std::move(mesh));
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
  // 2D: Triangle and Quadrilateral × K ∈ {1,2,3,4,6}
  // =========================================================================

  class H1Poisson2D : public ::testing::TestWithParam<Param> {};

  /**
   * @brief Distributed 2D Poisson with H1 of degree K.
   *
   * Mesh: 16×16 uniform grid (h = 1/15).
   * Solution: sin(πx)sin(πy).
   * Tolerance: @ref RODIN_FUZZY_CONSTANT.
   */
  TEST_P(H1Poisson2D, SimpleSine)
  {
    const auto& [geom, K] = GetParam();
    auto& world = *g_world;

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, geom, { 16, 16 });

    const Real globalError = run2D(K, world, mesh);
    EXPECT_NEAR(globalError, 0, RODIN_FUZZY_CONSTANT)
        << "np=" << world.size()
        << " geometry=" << static_cast<int>(geom)
        << " K=" << K;
  }

  INSTANTIATE_TEST_SUITE_P(
    AllGeometriesAndDegrees,
    H1Poisson2D,
    testing::Combine(
      testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral),
      testing::Values<size_t>(1, 2, 3, 4, 6)
    ),
    ParamName{}
  );

  // =========================================================================
  // 3D: Tetrahedron, Hexahedron, and Wedge × K ∈ {1,2,3,4,6}
  // =========================================================================

  class H1Poisson3D : public ::testing::TestWithParam<Param> {};

  /**
   * @brief Distributed 3D Poisson with H1 of degree K.
   *
   * Mesh: 5×5×5 uniform grid (h = 1/4, 4³ cells).
   * Solution: x(1−x)y(1−y)z(1−z) — degree-6 polynomial, zero on ∂Ω.
   * Tolerance: 10 × @ref RODIN_FUZZY_CONSTANT.
   *
   * A polynomial solution is used instead of sin to avoid the
   * @f$ \pi^{K+1} @f$ norm blow-up that makes coarse-mesh trigonometric
   * solutions inaccurate for K = 1, 2.  For K ≥ 6 the solution is
   * representable exactly in the FE space.
   */
  TEST_P(H1Poisson3D, PolynomialBC)
  {
    const auto& [geom, K] = GetParam();
    auto& world = *g_world;

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, geom, { 5, 5, 5 });

    const Real globalError = run3D(K, world, mesh);
    EXPECT_NEAR(globalError, 0, 10 * RODIN_FUZZY_CONSTANT)
        << "np=" << world.size()
        << " geometry=" << static_cast<int>(geom)
        << " K=" << K;
  }

  INSTANTIATE_TEST_SUITE_P(
    AllGeometriesAndDegrees,
    H1Poisson3D,
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
