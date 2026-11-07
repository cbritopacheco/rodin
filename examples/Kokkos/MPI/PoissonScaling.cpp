/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file PoissonScaling.cpp
 * @brief Times a distributed Poisson problem assembled and solved on the
 * devices, one process per device.
 *
 * Every process generates its block of a uniform grid of the unit square or
 * cube, builds the device mesh and the H1 space of order k, then assembles
 * and solves
 * @f[
 *   -\Delta u = 1 \text{ in } \Omega, \quad u = 0 \text{ on } \partial\Omega
 * @f]
 * with PETSc conjugate gradients, algebraic multigrid by default. The
 * matrices and vectors are PETSc Kokkos objects, so the solver runs on the
 * devices too. PETSc options override the solver.
 *
 * Usage: KokkosMPIPoissonScaling [type] [n] [k] [repeat], with type one of
 * triangle, quadrilateral, tetrahedron or hexahedron, n the number of
 * vertices along each axis, k the order, 1 or 2. Rank 0 prints one line of
 * key=value timings, in seconds, per repetition; each one is the time between
 * two barriers after a fence, so the slowest process. The first repetition
 * also builds the sparsity pattern and the multigrid hierarchy, and every
 * solve starts from zero. With -log_view, the first repetition and the next
 * ones are reported in two stages.
 */

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <string>
#include <type_traits>

#include <petsc.h>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Solver.h>
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
  using HostMesh = Geometry::Mesh<Context::MPI>;
  using DeviceContext = Context::Kokkos<Context::MPI, Backend>;
  using DeviceMesh = Geometry::Mesh<DeviceContext>;

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

  Geometry::Polytope::Type parseType(const std::string& name)
  {
    using Type = Geometry::Polytope::Type;
    if (name == "triangle")
      return Type::Triangle;
    if (name == "quadrilateral")
      return Type::Quadrilateral;
    if (name == "tetrahedron")
      return Type::Tetrahedron;
    if (name == "hexahedron")
      return Type::Hexahedron;
    std::cerr << "Unknown cell type " << name << std::endl;
    std::exit(EXIT_FAILURE);
  }

  using P1Space = Variational::P1<Real, DeviceMesh>;
  using P2Space = Variational::H1<2, Real, DeviceMesh>;

  /**
   * @brief Assembles and solves @p pb with @p solver, @p repeat times.
   *
   * The first repetition builds the sparsity pattern and the multigrid
   * hierarchy; the next ones reuse the pattern, as a time loop would.
   */
  template <class Problem>
  void time(const boost::mpi::communicator& comm, Problem& pb, PETSc::Solver::CG& solver,
    const std::string& prefix, size_t repeat)
  {
    Stopwatch time(comm);
    solver.setTolerances(Real(1e-8), Real(1e-50), Real(1e5), 10000);
    ::PC pc;
    KSPGetPC(solver.getHandle(), &pc);
    PCSetType(pc, PCGAMG);
    for (size_t r = 0; r < repeat; ++r)
    {
      PetscLogStagePush(getRepetitionStage(r));
      const double assemble = time([&] { pb.assemble(); });
      // Each solve starts from zero, not from the previous solution.
      VecSet(pb.getLinearSystem().getSolution(), 0);
      const double solve = time([&] { pb.solve(solver); });
      PetscInt its = 0;
      KSPGetIterationNumber(solver.getHandle(), &its);
      KSPConvergedReason reason;
      KSPGetConvergedReason(solver.getHandle(), &reason);
      if (comm.rank() == 0)
      {
        std::cout << prefix << " rep=" << r << " assemble=" << assemble
                  << " solve=" << solve << " its=" << its << " reason=" << reason
                  << std::endl;
      }
      PetscLogStagePop();
    }
  }

  // The problems are built outside of templates: the deduction of the forms
  // from their operands does not parse in a dependent context.

  void run(const boost::mpi::communicator& comm, P1Space& fes, const std::string& prefix,
    size_t repeat)
  {
    using namespace Rodin::Variational;
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction v(fes);
    RealFunction f = 1;
    Problem pb(u, v);
    pb = Integral(Grad(u), Grad(v)) - Integral(f, v) + DirichletBC(u, Zero());
    PETSc::Solver::CG solver(pb);
    time(comm, pb, solver, prefix, repeat);
  }

  void run(const boost::mpi::communicator& comm, P2Space& fes, const std::string& prefix,
    size_t repeat)
  {
    using namespace Rodin::Variational;
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction v(fes);
    RealFunction f = 1;
    Problem pb(u, v);
    pb = Integral(Grad(u), Grad(v)) - Integral(f, v) + DirichletBC(u, Zero());
    PETSc::Solver::CG solver(pb);
    time(comm, pb, solver, prefix, repeat);
  }
}

int main(int argc, char** argv)
{
  using namespace Rodin::Examples;
  namespace Context = Rodin::Context;
  namespace Geometry = Rodin::Geometry;
  namespace Variational = Rodin::Variational;
  using Rodin::Real;
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  const std::string typeName = argc > 1 ? argv[1] : "triangle";
  const size_t n = argc > 2 ? std::stoul(argv[2]) : 1024;
  const size_t k = argc > 3 ? std::stoul(argv[3]) : 1;
  const size_t repeat = argc > 4 ? std::stoul(argv[4]) : 3;
  const auto type = parseType(typeName);
  {
    ::Kokkos::ScopeGuard runtime(argc, argv);
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    try
    {
      Context::MPI ctx(env, world);
      Stopwatch time(world);
      using Type = Geometry::Polytope::Type;
      const size_t d = (type == Type::Triangle || type == Type::Quadrilateral) ? 2 : 3;
      HostMesh host(ctx);
      const double tmesh = time([&] {
        if (d == 2)
          host = HostMesh::UniformGrid(ctx, type, {n, n});
        else
          host = HostMesh::UniformGrid(ctx, type, {n, n, n});
        host.scale(Real(1) / static_cast<Real>(n - 1));
        auto& conn = host.getConnectivity();
        conn.compute(d, d);
        conn.compute(d, 0);
        conn.compute(d, d - 1);
        conn.compute(d - 1, d);
        conn.compute(d - 1, 0);
        host.reconcile(d - 1);
        if (k > 1 && d == 3)
        {
          conn.compute(d, 1);
          conn.compute(1, 0);
          conn.compute(2, 1);
          host.reconcile(1);
        }
      });
      // A first copy warms the device up; the second one is timed.
      {
        const DeviceMesh warmup(DeviceContext(ctx), host);
      }
      std::optional<DeviceMesh> device;
      const double tdevice = time([&] { device.emplace(DeviceContext(ctx), host); });
      const DeviceMesh& mesh = *device;

      std::string prefix;
      const auto report = [&](size_t dofs, double tfes) {
        prefix = "np=" + std::to_string(world.size()) + " type=" + typeName +
          " n=" + std::to_string(n) + " k=" + std::to_string(k) +
          " dofs=" + std::to_string(dofs);
        if (world.rank() == 0)
        {
          std::cout << prefix << " mesh=" << tmesh << " device=" << tdevice
                    << " fes=" << tfes << std::endl;
        }
      };
      if (k == 1)
      {
        std::optional<P1Space> fes;
        const double tfes = time([&] { fes.emplace(mesh); });
        report(fes->getSize(), tfes);
        run(world, *fes, prefix, repeat);
      }
      else if (k == 2)
      {
        std::optional<P2Space> fes;
        const double tfes = time([&] { fes.emplace(std::integral_constant<size_t, 2>(), mesh); });
        report(fes->getSize(), tfes);
        run(world, *fes, prefix, repeat);
      }
      else
      {
        std::cerr << "Unsupported order " << k << std::endl;
      }
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
