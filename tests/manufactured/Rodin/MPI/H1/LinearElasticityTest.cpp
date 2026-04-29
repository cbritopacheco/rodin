/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file LinearElasticityTest.cpp
 * @brief Distributed MPI H1 linear-elasticity manufactured tests.
 *
 * Parametrized over (geometry, K). 2D tests use K ∈ {1, 2, 3, 4, 6};
 * 3D tests use K ∈ {2, 3, 4, 6}.
 *
 * Geometries:
 *   2D — Triangle, Quadrilateral
 *   3D — Tetrahedron, Hexahedron, Wedge
 *
 * Each GTest parametrization is invoked by CTest with
 * MPIEXEC_NUMPROC_FLAG set to 1, 2, 3, 4, and 6 ranks.
 *
 * ### 2D tests
 * Manufactured solution @f$ \mathbf{u} = (\sin(\pi x)\sin(\pi y),
 * \sin(\pi x)\sin(\pi y)) @f$ on a 16×16 uniform grid (h = 1/15).
 * The global @f$ L^2 @f$ error
 * @f$ \int_\Omega \|\mathbf{u} - \mathbf{u}_h\|^2 @f$ must stay below
 * @ref RODIN_FUZZY_CONSTANT.
 *
 * Lamé parameters: @f$ \lambda = \mu = 1 @f$.
 * Body force (derived from strong form):
 * @f[
 *   f_1 = f_2 = \pi^2\bigl[(\lambda+3\mu)\sin(\pi x)\sin(\pi y)
 *                         - (\lambda+\mu)\cos(\pi x)\cos(\pi y)\bigr]
 * @f]
 *
 * ### 3D tests
 * Manufactured solution @f$ \mathbf{u} = (x(1-x),\, y(1-y),\, z(1-z)) @f$
 * on a 5×5×5 uniform grid (h = 1/4, 4³ = 64 cells).
 * The polynomial solution is imposed as inhomogeneous Dirichlet data.
 * For this solution the body force is constant:
 * @f[
 *   f_1 = f_2 = f_3 = 2(\lambda + 2\mu)
 * @f]
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
        / ("rodin_mpi_h1_le_"
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
  // Per-K linear-elasticity solvers
  // -------------------------------------------------------------------------

  /**
   * @brief Solves the 2D linear-elasticity problem with H1 of degree K and
   *        returns the globally reduced @f$ \int_\Omega \|\mathbf{u}-\mathbf{u}_h\|^2 @f$.
   *
   * Manufactured solution: @f$ \mathbf{u} = \sin(\pi x)\sin(\pi y)\,(1,1)^T @f$.
   * Lamé parameters: @f$ \lambda = \mu = 1 @f$.
   * Body force:
   * @f$ f_i = \pi^2[(\lambda+3\mu)\sin(\pi x)\sin(\pi y)
   *                  - (\lambda+\mu)\cos(\pi x)\cos(\pi y)] @f$
   */
  template <size_t K>
  Real computeError2D(
      boost::mpi::communicator& world,
      Geometry::Mesh<Context::MPI>& mesh)
  {
    namespace PETSc = ::Rodin::PETSc;
    using VFES = Variational::H1<K, Math::SpatialVector<Real>, Geometry::Mesh<Context::MPI>>;
    using SFES = Variational::H1<K, Real, Geometry::Mesh<Context::MPI>>;

    constexpr auto order = std::integral_constant<size_t, K>{};
    const auto pi = Math::Constants::pi();
    constexpr Real lambda = 1.0;
    constexpr Real mu     = 1.0;

    VFES vh(order, mesh, mesh.getSpaceDimension());

    auto A = sin(pi * F::x) * sin(pi * F::y);
    auto B = cos(pi * F::x) * cos(pi * F::y);
    VectorFunction f{
      pi * pi * ((lambda + 3 * mu) * A - (lambda + mu) * B),
      pi * pi * ((lambda + 3 * mu) * A - (lambda + mu) * B)
    };

    PETSc::Variational::TrialFunction u(vh);
    PETSc::Variational::TestFunction  v(vh);

    Problem elasticity(u, v);
    elasticity = Integral(lambda * Div(u), Div(v))
               + Integral(
                   mu * (Jacobian(u) + Jacobian(u).T()),
                   0.5 * (Jacobian(v) + Jacobian(v).T()))
               - Integral(f, v)
               + DirichletBC(u, Zero());

    PETSc::Solver::CG solver(elasticity);
    solver.solve();

    VectorFunction solution{ A, A };

    SFES sh(order, mesh);
    GridFunction<SFES, ::Vec> diff(sh);
    diff = Pow(Frobenius(u.getSolution() - solution), 2);

    const Real localError = Integral(diff).compute();
    Real globalError = 0;
    boost::mpi::all_reduce(world, localError, globalError, std::plus<Real>());
    return globalError;
  }

  /**
   * @brief Solves the 3D linear-elasticity problem with H1 of degree K and
   *        returns the globally reduced @f$ \int_\Omega \|\mathbf{u}-\mathbf{u}_h\|^2 @f$.
   *
   * Manufactured solution: @f$ \mathbf{u} = (x(1-x),\,y(1-y),\,z(1-z)) @f$.
   * Lamé parameters: @f$ \lambda = \mu = 1 @f$.
   * Body force: @f$ \mathbf{f} = (2(\lambda+2\mu),\, 2(\lambda+2\mu),\, 2(\lambda+2\mu)) @f$.
   *
   * The polynomial solution is imposed as inhomogeneous Dirichlet data.
   */
  template <size_t K>
  Real computeError3D(
      boost::mpi::communicator& world,
      Geometry::Mesh<Context::MPI>& mesh)
  {
    namespace PETSc = ::Rodin::PETSc;
    using VFES = Variational::H1<K, Math::SpatialVector<Real>, Geometry::Mesh<Context::MPI>>;
    using SFES = Variational::H1<K, Real, Geometry::Mesh<Context::MPI>>;

    constexpr auto order = std::integral_constant<size_t, K>{};
    constexpr Real lambda = 1.0;
    constexpr Real mu     = 1.0;

    VFES vh(order, mesh, mesh.getSpaceDimension());

    const Real c = 2.0 * (lambda + 2.0 * mu);
    VectorFunction f{ c, c, c };

    PETSc::Variational::TrialFunction u(vh);
    PETSc::Variational::TestFunction  v(vh);

    VectorFunction solution{
      F::x * (1 - F::x),
      F::y * (1 - F::y),
      F::z * (1 - F::z)
    };

    Problem elasticity(u, v);
    elasticity = Integral(lambda * Div(u), Div(v))
               + Integral(
                   mu * (Jacobian(u) + Jacobian(u).T()),
                   0.5 * (Jacobian(v) + Jacobian(v).T()))
               - Integral(f, v)
               + DirichletBC(u, solution);

    PETSc::Solver::CG solver(elasticity);
    solver.solve();

    SFES sh(order, mesh);
    GridFunction<SFES, ::Vec> diff(sh);
    diff = Pow(Frobenius(u.getSolution() - solution), 2);

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
namespace Rodin::Tests::Manufactured::MPI::H1LinearElasticity
{
  using Param = std::tuple<Polytope::Type, size_t>;

  // =========================================================================
  // 2D: Triangle and Quadrilateral × K ∈ {1,2,3,4,6}
  // =========================================================================

  class H1LinearElasticity2D : public ::testing::TestWithParam<Param> {};

  /**
   * @brief Distributed 2D linear elasticity with H1 of degree K.
   *
   * Mesh: 16×16 uniform grid (h = 1/15).
   * Solution: sin(πx)sin(πy)·(1,1)^T.
   * Tolerance: @ref RODIN_FUZZY_CONSTANT.
   */
  TEST_P(H1LinearElasticity2D, SineSolution)
  {
    const auto& [geom, K] = GetParam();
    auto& world = *g_world;

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, geom, { 16, 16 });

    const Real globalError = run2D(K, world, mesh);
    const Real tolerance = (K == 1 ? 5 : 1) * RODIN_FUZZY_CONSTANT;
    EXPECT_NEAR(globalError, 0, tolerance)
        << "np=" << world.size()
        << " geometry=" << static_cast<int>(geom)
        << " K=" << K;
  }

  INSTANTIATE_TEST_SUITE_P(
    AllGeometriesAndDegrees,
    H1LinearElasticity2D,
    testing::Combine(
      testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral),
      testing::Values<size_t>(1, 2, 3, 4, 6)
    ),
    ParamName{}
  );

  // =========================================================================
  // 3D: Tetrahedron, Hexahedron, and Wedge × K ∈ {2,3,4,6}
  // =========================================================================

  class H1LinearElasticity3D : public ::testing::TestWithParam<Param> {};

  /**
   * @brief Distributed 3D linear elasticity with H1 of degree K.
   *
   * Mesh: 5×5×5 uniform grid (h = 1/4, 4³ cells).
   * Solution: (x(1−x), y(1−y), z(1−z)) — per-component degree-2 polynomial,
   * imposed as inhomogeneous Dirichlet data. Exactly representable for K ≥ 2.
   * Tolerance: 10 × @ref RODIN_FUZZY_CONSTANT.
   *
   * A polynomial solution is used to avoid the @f$ \pi^{K+1} @f$ norm
   * blow-up that makes coarse-mesh trigonometric solutions inaccurate
   * for low degrees.
   */
  TEST_P(H1LinearElasticity3D, PolynomialSolution)
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
    H1LinearElasticity3D,
    testing::Combine(
      testing::Values(
        Polytope::Type::Tetrahedron,
        Polytope::Type::Hexahedron,
        Polytope::Type::Wedge),
      testing::Values<size_t>(2, 3, 4, 6)
    ),
    ParamName{}
  );
} // namespace Rodin::Tests::Manufactured::MPI::H1LinearElasticity

// ---------------------------------------------------------------------------
// main() — initializes MPI and PETSc before running all tests.
// ---------------------------------------------------------------------------
int main(int argc, char** argv)
{
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  g_env   = &env;
  g_world = &world;

  ::testing::InitGoogleTest(&argc, argv);

  [[maybe_unused]] PetscErrorCode ierr =
      PetscInitialize(&argc, &argv, nullptr, nullptr);
  assert(ierr == PETSC_SUCCESS);

  const int result = RUN_ALL_TESTS();

  PetscFinalize();
  return result;
}
