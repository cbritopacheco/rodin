/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file OseenScaling.cpp
 * @brief Times a distributed Oseen step assembled and solved on the
 * devices, one process per device.
 *
 * Rank 0 loads a MEDIT mesh of tetrahedra, or builds a tetrahedral grid of
 * the unit cube, partitions it and sends one shard to each process. Every
 * process builds its device mesh and the P1/P1 spaces, then assembles and
 * solves the linearized Navier-Stokes step of KokkosOseenAssembly,
 *
 *   rho/dt (u, v) + mu (Jac u, Jac v) + rho (Jac(u) w, v)
 *   + rho/2 (div(w) u, v) - (p, div v) + (div u, q) + eps (grad p, grad q)
 *   + delta (p, q) = rho/dt (w, v),
 *
 * with u = 0 on the wall, attribute 4 of the atrium meshes, and p = 0 on
 * the outlet, attribute 6, or u = 0 on the whole boundary of the grid. A MEDIT mesh is taken in millimeters and scaled to
 * meters, as examples/Heart/Atrium does. The matrices and vectors are PETSc
 * Kokkos objects, so the solver runs on the devices too: BiCGStab with
 * Jacobi by default, which PETSc options override.
 *
 * The option -oseen_partitioner greedy partitions the mesh with
 * GreedyPartitioner instead of BalancedCompactPartitioner, and the
 * first repetition prints the extreme local row and ghost column counts of
 * the operator. The options -oseen_eps, -oseen_scale and -oseen_wscale override the
 * pressure stabilization eps, the scale of a MEDIT mesh, 1e-3 by default,
 * and the factor of the advection field, 0.1 on a MEDIT mesh and 1 on the
 * grid. The option -oseen_refine k splits every cell of the mesh into
 * eight, k times, on rank 0 before the partition, and divides the time step
 * by 2^k to keep the CFL number; -oseen_dt overrides the step.
 *
 * Usage: KokkosMPIOseenScaling [mesh.mesh | n] [repeat] [renumber]. A number
 * n builds a grid of n cells per side. The word renumber orders the vertices
 * and the cells along space before the partition, see Mesh::renumber. Rank 0
 * prints one line of key=value timings, in seconds, per repetition; each one
 * is the time between two barriers after a fence, so the slowest process.
 * The first repetition also builds the sparsity pattern and the
 * preconditioner, and every solve starts from zero. With -log_view, the first
 * repetition and the next ones are reported in two stages.
 */

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include <petsc.h>
#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/broadcast.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#include <Rodin/Geometry/GreedyPartitioner.h>
#include <Rodin/IO.h>
#include <Rodin/Solver.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/Kokkos/MPI.h>

#if defined(KOKKOS_ENABLE_CUDA)
using Backend = Rodin::Context::CUDA;
#elif defined(KOKKOS_ENABLE_HIP)
using Backend = Rodin::Context::HIP;
#elif defined(KOKKOS_ENABLE_SYCL)
using Backend = Rodin::Context::SYCL;
#else
using Backend = ::Kokkos::DefaultHostExecutionSpace;
#endif

namespace Rodin::Examples
{
  using LocalMesh = Geometry::Mesh<Context::Local>;
  using HostMesh = Geometry::Mesh<Context::MPI>;
  using DeviceContext = Context::Kokkos<Context::MPI, Backend>;
  using DeviceMesh = Geometry::Mesh<DeviceContext>;

  /// @brief Physical and stabilization parameters of the Oseen step.
  struct Parameters
  {
    Real rho = 1060;
    Real dt = 1e-3;
    Real mu = 0.0049;
    Real eps = 1e-3;

    /// @brief Pressure penalty of examples/Heart/Atrium.
    Real delta = 1e-12;
  };

  /// @brief Attribute of the wall faces of the atrium meshes.
  constexpr Geometry::Attribute Wall = 4;

  /// @brief Attribute of the outlet faces of the atrium meshes.
  constexpr Geometry::Attribute Outlet = 6;

  /// @brief Wall time of a step on the slowest process.
  class Stopwatch
  {
    public:
      explicit Stopwatch(const boost::mpi::communicator& comm)
        : m_comm(comm)
      {}

      template <class Function>
      double operator()(Function&& f)
      {
        ::Kokkos::fence();
        m_comm.barrier();
        const auto t0 = std::chrono::steady_clock::now();
        f();
        ::Kokkos::fence();
        m_comm.barrier();
        const auto t1 = std::chrono::steady_clock::now();
        return std::chrono::duration<double>(t1 - t0).count();
      }

    private:
      const boost::mpi::communicator& m_comm;
  };

  /**
   * @brief PETSc log stage of repetition @p r, so that -log_view reports
   * the first repetition apart from the next ones.
   *
   * The first repetition also loads the device kernels and builds the
   * pattern and the preconditioner, which can take far longer than the
   * repetitions that follow.
   */
  PetscLogStage getRepetitionStage(size_t r)
  {
    const char* name = r == 0 ? "First repetition" : "Next repetitions";
    PetscLogStage stage = -1;
    PetscLogStageGetId(name, &stage);
    if (stage < 0)
      PetscLogStageRegister(name, &stage);
    return stage;
  }

  /// @returns Whether @p source is a cell count rather than a path.
  bool isGrid(const std::string& source)
  {
    return !source.empty() && std::all_of(source.begin(), source.end(), [](char c) {
      return c >= '0' && c <= '9';
    });
  }

  /**
   * @brief Splits every tetrahedron of @p mesh into eight, and every
   * triangle and segment of it into four and two.
   *
   * The new vertices are the midpoints of the edges. The interior octahedron
   * of each tetrahedron is cut along its shortest diagonal, and every child
   * keeps the orientation and the attribute of its parent, so the grading of
   * the mesh and its boundary attributes are kept with half the size.
   */
  LocalMesh refine(const LocalMesh& mesh)
  {
    using Geometry::Polytope;
    const auto& conn = mesh.getConnectivity();
    const size_t nv = mesh.getVertexCount();
    std::unordered_map<std::uint64_t, Index> midpoints;
    midpoints.reserve(2 * mesh.getCellCount());
    std::vector<Math::SpatialPoint> coordinates;
    const auto mid = [&](Index a, Index b) {
      const auto key =
        (static_cast<std::uint64_t>(std::min(a, b)) << 32) | std::max(a, b);
      const auto [it, inserted] = midpoints.emplace(key, nv + coordinates.size());
      if (inserted)
      {
        coordinates.push_back(
          0.5 * (mesh.getVertexCoordinates(a) + mesh.getVertexCoordinates(b)));
      }
      return it->second;
    };

    std::vector<std::array<Index, 4>> tets;
    std::vector<Optional<Geometry::Attribute>> tattrs;
    tets.reserve(8 * conn.getCount(3));
    tattrs.reserve(8 * conn.getCount(3));
    std::vector<Math::SpatialPoint> x(nv);
    for (Index i = 0; i < nv; ++i)
      x[i] = mesh.getVertexCoordinates(i);
    const auto volume = [&](const std::array<Index, 4>& t) {
      const auto at = [&](Index k) -> Math::SpatialPoint {
        return k < nv ? x[k] : coordinates[k - nv];
      };
      const Math::SpatialPoint a = at(t[1]) - at(t[0]);
      const Math::SpatialPoint b = at(t[2]) - at(t[0]);
      const Math::SpatialPoint c = at(t[3]) - at(t[0]);
      return a.x() * (b.y() * c.z() - b.z() * c.y()) -
        a.y() * (b.x() * c.z() - b.z() * c.x()) + a.z() * (b.x() * c.y() - b.y() * c.x());
    };
    const auto length = [&](Index a, Index b) {
      const auto at = [&](Index k) -> Math::SpatialPoint {
        return k < nv ? x[k] : coordinates[k - nv];
      };
      return (at(a) - at(b)).squaredNorm();
    };
    for (Index i = 0; i < conn.getCount(3); ++i)
    {
      const auto& v = conn.getPolytope(3, i);
      const Index m01 = mid(v[0], v[1]), m02 = mid(v[0], v[2]), m03 = mid(v[0], v[3]);
      const Index m12 = mid(v[1], v[2]), m13 = mid(v[1], v[3]), m23 = mid(v[2], v[3]);
      const bool positive = volume({v[0], v[1], v[2], v[3]}) > 0;
      std::array<std::array<Index, 4>, 8> children{{{v[0], m01, m02, m03},
        {m01, v[1], m12, m13}, {m02, m12, v[2], m23}, {m03, m13, m23, v[3]}}};
      // The three diagonals of the octahedron, each with the two pairs of
      // opposite vertices left around it.
      const std::array<std::array<Index, 6>, 3> diagonals{{{m01, m23, m02, m13, m03, m12},
        {m02, m13, m01, m23, m03, m12}, {m03, m12, m01, m23, m02, m13}}};
      const auto& d = *std::min_element(
        diagonals.begin(), diagonals.end(), [&](const auto& a, const auto& b) {
          return length(a[0], a[1]) < length(b[0], b[1]);
        });
      const std::array<Index, 4> ring{d[2], d[4], d[3], d[5]};
      for (size_t k = 0; k < 4; ++k)
        children[4 + k] = {d[0], d[1], ring[k], ring[(k + 1) % 4]};
      const auto attr = mesh.getAttribute(3, i);
      for (auto& c : children)
      {
        if ((volume(c) > 0) != positive)
          std::swap(c[0], c[1]);
        tets.push_back(c);
        tattrs.push_back(attr);
      }
    }

    // Every midpoint is made before the vertices are added.
    std::vector<std::pair<std::array<Index, 2>, Optional<Geometry::Attribute>>> segments;
    for (Index i = 0; i < conn.getCount(1); ++i)
    {
      const auto& v = conn.getPolytope(1, i);
      const Index m = mid(v[0], v[1]);
      const auto attr = mesh.getAttribute(1, i);
      segments.push_back({{v[0], m}, attr});
      segments.push_back({{m, v[1]}, attr});
    }
    std::vector<std::pair<std::array<Index, 3>, Optional<Geometry::Attribute>>> triangles;
    for (Index i = 0; i < conn.getCount(2); ++i)
    {
      const auto& v = conn.getPolytope(2, i);
      const Index m01 = mid(v[0], v[1]), m12 = mid(v[1], v[2]), m02 = mid(v[0], v[2]);
      const auto attr = mesh.getAttribute(2, i);
      triangles.push_back({{v[0], m01, m02}, attr});
      triangles.push_back({{m01, v[1], m12}, attr});
      triangles.push_back({{m02, m12, v[2]}, attr});
      triangles.push_back({{m01, m12, m02}, attr});
    }

    LocalMesh::Builder build;
    build.initialize(3).nodes(nv + coordinates.size());
    for (Index i = 0; i < nv; ++i)
    {
      build.vertex(x[i]);
      build.attribute({0, i}, mesh.getAttribute(0, i));
    }
    for (const auto& c : coordinates)
      build.vertex(c);
    build.reserve(1, segments.size());
    for (const auto& [s, attr] : segments)
    {
      Index index;
      build.polytope(Polytope::Type::Segment, {s[0], s[1]}, index);
      build.attribute({1, index}, attr);
    }
    build.reserve(2, triangles.size());
    for (const auto& [t, attr] : triangles)
    {
      Index index;
      build.polytope(Polytope::Type::Triangle, {t[0], t[1], t[2]}, index);
      build.attribute({2, index}, attr);
    }
    build.reserve(3, tets.size());
    for (size_t i = 0; i < tets.size(); ++i)
    {
      const auto& t = tets[i];
      Index index;
      build.polytope(Polytope::Type::Tetrahedron, {t[0], t[1], t[2], t[3]}, index);
      build.attribute({3, index}, tattrs[i]);
    }
    return build.finalize();
  }

  /**
   * @brief Loads or builds the mesh on rank 0 and sends a shard to every
   * process.
   * @param ctx MPI context.
   * @param source MEDIT path, or the number of cells per side of a grid.
   * @param renumber Whether to order the mesh along space first.
   * @param cells Set to the global number of cells.
   */
  HostMesh distribute(
    const Context::MPI& ctx, const std::string& source, bool renumber, size_t& cells)
  {
    const auto& comm = ctx.getCommunicator();
    Geometry::Sharder<Context::MPI> sharder(ctx);
    if (comm.rank() == 0)
    {
      LocalMesh mesh;
      if (isGrid(source))
      {
        const size_t n = std::stoul(source) + 1;
        mesh = LocalMesh::UniformGrid(Geometry::Polytope::Type::Tetrahedron, {n, n, n});
        mesh.scale(Real(1) / static_cast<Real>(n - 1));
      }
      else
      {
        mesh.load(source, IO::FileFormat::MEDIT);
        PetscReal scale = 1e-3;
        PetscOptionsGetReal(nullptr, nullptr, "-oseen_scale", &scale, nullptr);
        mesh.scale(scale);
      }
      PetscInt levels = 0;
      PetscOptionsGetInt(nullptr, nullptr, "-oseen_refine", &levels, nullptr);
      for (PetscInt l = 0; l < levels; ++l)
        mesh = refine(mesh);
      if (mesh.getPolytopeCount(Geometry::Polytope::Type::Tetrahedron) !=
        mesh.getCellCount())
        throw std::runtime_error(source + " has cells other than tetrahedra.");
      if (renumber)
        mesh = mesh.renumber();
      cells = mesh.getCellCount();
      auto& conn = mesh.getConnectivity();
      conn.compute(3, 3);
      conn.compute(3, 0);
      conn.compute(3, 2);
      conn.compute(2, 3);
      conn.compute(2, 0);
      char name[16] = "compact";
      PetscOptionsGetString(
        nullptr, nullptr, "-oseen_partitioner", name, sizeof(name), nullptr);
      if (std::string(name) == "greedy")
      {
        Geometry::GreedyPartitioner partitioner(mesh);
        partitioner.partition(static_cast<size_t>(comm.size()));
        sharder.shard(partitioner);
      }
      else
      {
        Geometry::BalancedCompactPartitioner partitioner(mesh);
        partitioner.partition(static_cast<size_t>(comm.size()));
        sharder.shard(partitioner);
      }
    }
    // The global mesh is released before the shards are sent.
    sharder.scatter(0);
    HostMesh mesh = sharder.gather(0);
    boost::mpi::broadcast(comm, cells, 0);
    auto& conn = mesh.getConnectivity();
    conn.compute(3, 3);
    conn.compute(3, 0);
    conn.compute(3, 2);
    conn.compute(2, 3);
    conn.compute(2, 0);
    mesh.reconcile(2);
    return mesh;
  }

  using VelocitySpace = Variational::P1<Math::SpatialVector<Real>, DeviceMesh>;
  using PressureSpace = Variational::P1<Real, DeviceMesh>;

  /**
   * @brief Prints the extreme local row and ghost column counts of the
   * operator over the processes, which bound the work and the halo of a
   * product.
   * @param comm World communicator.
   * @param A Assembled operator.
   * @param prefix Leading keys of the printed line.
   */
  void printLayout(
    const boost::mpi::communicator& comm, ::Mat A, const std::string& prefix)
  {
    PetscInt rows = 0, ghosts = 0;
    ::Mat Ad, Ao;
    const PetscInt* colmap = nullptr;
    MatGetLocalSize(A, &rows, nullptr);
    // A single process holds a sequential matrix, without ghost columns.
    if (comm.size() > 1 && MatMPIAIJGetSeqAIJ(A, &Ad, &Ao, &colmap) == PETSC_SUCCESS &&
      Ao)
      MatGetSize(Ao, nullptr, &ghosts);
    const std::array<long, 4> local{rows, -rows, ghosts, -ghosts};
    std::array<long, 4> extreme;
    boost::mpi::all_reduce(
      comm, local.data(), 4, extreme.data(), boost::mpi::maximum<long>());
    if (comm.rank() == 0)
    {
      std::cout << prefix << " rows=" << -extreme[1] << ".." << extreme[0]
                << " ghosts=" << -extreme[3] << ".." << extreme[2] << std::endl;
    }
  }

  /**
   * @brief Assembles and solves the Oseen step @p repeat times.
   * @param comm World communicator.
   * @param host Distributed host mesh, on which the advection is set.
   * @param vh Device velocity space.
   * @param qh Device pressure space.
   * @param grid Whether the mesh is the grid, whose whole boundary is a wall.
   * @param prefix Leading keys of the printed lines.
   * @param repeat Number of repetitions.
   */
  void run(const boost::mpi::communicator& comm, const HostMesh& host, VelocitySpace& vh,
    PressureSpace& qh, bool grid, const std::string& prefix, size_t repeat)
  {
    using namespace Rodin::Variational;
    Parameters P;
    PetscOptionsGetReal(nullptr, nullptr, "-oseen_eps", &P.eps, nullptr);
    // Each refinement halves the cells, and the step with them, so the CFL
    // number of the advection stays the same.
    PetscInt levels = 0;
    PetscOptionsGetInt(nullptr, nullptr, "-oseen_refine", &levels, nullptr);
    P.dt = std::ldexp(P.dt, -static_cast<int>(levels));
    PetscOptionsGetReal(nullptr, nullptr, "-oseen_dt", &P.dt, nullptr);

    // The advection field of KokkosOseenAssembly, interpolated on the host
    // and copied to the device space, whose degrees of freedom are laid out
    // the same way.
    P1 hvh(host, 3);
    PETSc::Variational::GridFunction hw(hvh);
    hw = VectorFunction{[](const Geometry::Point& p) { return std::sin(p.y()) + 0.5; },
      [](const Geometry::Point& p) { return std::cos(p.x()) * p.z(); },
      [](const Geometry::Point& p) { return p.x() - p.y(); }};
    PETSc::Variational::GridFunction w(vh);
    VecCopy(hw.getData(), w.getData());
    // Without stabilization, the solvers diverge once the advection rules
    // over rho/dt on the smallest cells of the atrium, 0.06 mm wide: the
    // field, about 0.75 m/s there, is scaled down to 7.5 cm/s by default.
    PetscReal wscale = grid ? 1 : 0.1;
    PetscOptionsGetReal(nullptr, nullptr, "-oseen_wscale", &wscale, nullptr);
    VecScale(w.getData(), wscale);
    VecGhostUpdateBegin(w.getData(), INSERT_VALUES, SCATTER_FORWARD);
    VecGhostUpdateEnd(w.getData(), INSERT_VALUES, SCATTER_FORWARD);

    PETSc::Variational::TrialFunction u(vh);
    PETSc::Variational::TrialFunction p(qh);
    PETSc::Variational::TestFunction v(vh);
    PETSc::Variational::TestFunction q(qh);
    const VectorFunction zero{Real(0), Real(0), Real(0)};

    Problem oseen(u, p, v, q);
    oseen = Integral(P.rho / P.dt * u, v) + Integral(P.mu * Jacobian(u), Jacobian(v)) +
      Integral(P.rho * Jacobian(u) * w, v) + Integral(0.5 * P.rho * Div(w) * u, v) -
      Integral(p, Div(v)) + Integral(Div(u), q) + Integral(P.eps * Grad(p), Grad(q)) +
      Integral(P.delta * p, q) - Integral(P.rho / P.dt * w, v);
    if (grid)
      oseen += DirichletBC(u, zero);
    else
    {
      // The constant pressure is held by the openings only, too weakly for
      // the iterative solvers, so the outlet fixes it as the pressure driven
      // openings of examples/Heart/Atrium do.
      oseen += DirichletBC(u, zero).on(Wall);
      oseen += DirichletBC(p, Zero()).on(Outlet);
    }

    PETSc::Solver::KSP solver(oseen);
    solver.setType(KSPBCGS);
    solver.setTolerances(Real(1e-8), Real(1e-50), Real(1e5), 10000);
    ::PC pc;
    KSPGetPC(solver.getHandle(), &pc);
    PCSetType(pc, PCJACOBI);

    Stopwatch time(comm);
    for (size_t r = 0; r < repeat; ++r)
    {
      PetscLogStagePush(getRepetitionStage(r));
      const double assemble = time([&] { oseen.assemble(); });
      // Each solve starts from zero, not from the previous solution.
      VecSet(oseen.getLinearSystem().getSolution(), 0);
      const double solve = time([&] { oseen.solve(solver); });
      PetscInt its = 0;
      KSPGetIterationNumber(solver.getHandle(), &its);
      KSPConvergedReason reason;
      KSPGetConvergedReason(solver.getHandle(), &reason);
      // The norm of the solution checks that every process count solves the
      // same problem.
      PetscReal norm = 0;
      VecNorm(oseen.getLinearSystem().getSolution(), NORM_2, &norm);
      if (r == 0)
        printLayout(comm, oseen.getLinearSystem().getOperator(), prefix);
      if (comm.rank() == 0)
      {
        std::cout << prefix << " rep=" << r << " assemble=" << assemble
                  << " solve=" << solve << " its=" << its << " reason=" << reason
                  << " norm=" << norm << std::endl;
      }
      PetscLogStagePop();
    }
  }
}

int main(int argc, char** argv)
{
  using namespace Rodin::Examples;
  namespace Context = Rodin::Context;
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  const std::string source = argc > 1 ? argv[1] : "32";
  const size_t repeat = argc > 2 ? std::stoul(argv[2]) : 3;
  const bool renumber = argc > 3 && std::string(argv[3]) == "renumber";
  {
    ::Kokkos::ScopeGuard runtime(argc, argv);
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    try
    {
      Context::MPI ctx(env, world);
      Stopwatch time(world);
      HostMesh host(ctx);
      size_t cells = 0;
      const double tmesh = time([&] { host = distribute(ctx, source, renumber, cells); });
      // A first copy warms the device up; the second one is timed.
      {
        const DeviceMesh warmup(DeviceContext(ctx), host);
      }
      std::optional<DeviceMesh> device;
      const double tdevice = time([&] { device.emplace(DeviceContext(ctx), host); });
      const DeviceMesh& mesh = *device;

      std::optional<VelocitySpace> vh;
      std::optional<PressureSpace> qh;
      const double tfes = time([&] {
        vh.emplace(mesh, 3);
        qh.emplace(mesh);
      });
      const size_t dofs = vh->getSize() + qh->getSize();
      const std::string name =
        isGrid(source) ? "grid" + source : source.substr(source.find_last_of('/') + 1);
      PetscInt levels = 0;
      PetscOptionsGetInt(nullptr, nullptr, "-oseen_refine", &levels, nullptr);
      const std::string prefix = "np=" + std::to_string(world.size()) + " mesh=" + name +
        (levels > 0 ? "+refine" + std::to_string(levels) : "") +
        (renumber ? "+renumber" : "") + " cells=" + std::to_string(cells) +
        " dofs=" + std::to_string(dofs);
      if (world.rank() == 0)
      {
        std::cout << prefix << " load=" << tmesh << " device=" << tdevice
                  << " fes=" << tfes << std::endl;
      }
      run(world, host, *vh, *qh, isGrid(source), prefix, repeat);
    }
    catch (const std::exception& e)
    {
      std::cerr << "rank " << world.rank() << ": " << e.what() << std::endl;
      world.abort(EXIT_FAILURE);
    }
    PetscFinalize();
  }
  return EXIT_SUCCESS;
}
