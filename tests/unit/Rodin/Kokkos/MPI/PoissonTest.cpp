/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file PoissonTest.cpp
 * @brief Distributed device assembly of PETSc problems against the host.
 *
 * Rank 0 builds and partitions a grid, every rank reloads its shard, and the
 * same problem is solved with the distributed host assembly and with the
 * distributed device assembly. The owned coefficients of both solutions must
 * agree to the solver tolerance.
 */

#include <algorithm>
#include <cmath>
#include <string>

#include <gtest/gtest.h>

#include <petsc.h>
#include <boost/filesystem.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>

#include <Rodin/Geometry.h>
#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#include <Rodin/IO.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/MPI/IO.h>
#include <Rodin/Kokkos/MPI.h>

#include "../TestContext.h"

#include <unistd.h>

namespace
{
  boost::mpi::environment* g_env = nullptr;
  boost::mpi::communicator* g_world = nullptr;
}

namespace Rodin::Tests::Unit
{
  namespace PETSc = ::Rodin::PETSc;
  namespace DeviceGeometry = Rodin::Kokkos::Geometry;
  using Polytope = Geometry::Polytope;
  using HostMesh = Geometry::Mesh<Context::MPI>;
  using DeviceContext = Context::Kokkos<Context::MPI, TestBackend>;
  using DeviceMesh = Geometry::Mesh<DeviceContext>;

  /**
   * @brief Distributes a unit grid from rank 0, through an HDF5 reload of
   * each shard as the host distributed tests do.
   */
  HostMesh distributeFromRoot(const Context::MPI& ctx, Polytope::Type type, size_t n)
  {
    const auto& comm = ctx.getCommunicator();
    Geometry::Sharder<Context::MPI> sharder(ctx);
    if (comm.rank() == 0)
    {
      auto mesh = Geometry::Mesh<Context::Local>::UniformGrid(type, {n, n});
      mesh.scale(Real(1) / static_cast<Real>(n - 1));
      const size_t d = mesh.getDimension();
      mesh.getConnectivity().compute(d, d);
      mesh.getConnectivity().compute(d, 0);
      mesh.getConnectivity().compute(d, d - 1);
      mesh.getConnectivity().compute(d - 1, d);
      mesh.getConnectivity().compute(d - 1, 0);
      Geometry::BalancedCompactPartitioner partitioner(mesh);
      partitioner.partition(static_cast<size_t>(comm.size()));
      sharder.shard(partitioner);
      sharder.scatter(0);
    }
    auto mesh = sharder.gather(0);

    const boost::filesystem::path file = boost::filesystem::temp_directory_path() /
      ("rodin_kokkos_mpi_poisson_" + std::to_string(getpid()) + "_r" +
        std::to_string(comm.rank()) + ".h5");
    {
      IO::MeshPrinter<IO::FileFormat::HDF5, Context::MPI> printer(mesh);
      printer.print(file);
    }
    comm.barrier();
    HostMesh res(ctx);
    res.load(file, IO::FileFormat::HDF5);
    const size_t d = res.getDimension();
    res.getConnectivity().compute(d, d);
    res.getConnectivity().compute(d, 0);
    res.getConnectivity().compute(d, d - 1);
    res.getConnectivity().compute(d - 1, d);
    res.getConnectivity().compute(d - 1, 0);
    res.getConnectivity().compute(1, 2);
    res.reconcile(1);
    comm.barrier();
    boost::filesystem::remove(file);
    return res;
  }

  /// @brief Solves a problem with conjugate gradients to a tight tolerance.
  template <class Problem>
  void solve(Problem& problem)
  {
    PETSc::Solver::CG solver(problem);
    solver.setTolerances(Real(1e-13), Real(1e-50), Real(1e5), 10000);
    solver.solve();
  }

  /**
   * @brief Expects two distributed vectors with the same layout to agree.
   * @param host Reference vector. @param device Compared vector.
   */
  void expectSameVector(::Vec host, ::Vec device)
  {
    PetscInt nh = 0, nd = 0;
    VecGetLocalSize(host, &nh);
    VecGetLocalSize(device, &nd);
    ASSERT_EQ(nh, nd);
    const PetscScalar* x = nullptr;
    const PetscScalar* y = nullptr;
    VecGetArrayRead(host, &x);
    VecGetArrayRead(device, &y);
    PetscReal local[2] = {0, 0};
    for (PetscInt i = 0; i < nh; ++i)
    {
      local[0] = std::max(local[0], PetscAbsScalar(x[i] - y[i]));
      local[1] = std::max(local[1], PetscAbsScalar(x[i]));
    }
    VecRestoreArrayRead(host, &x);
    VecRestoreArrayRead(device, &y);
    PetscReal global[2];
    MPI_Allreduce(local, global, 2, MPIU_REAL, MPI_MAX, *g_world);
    ASSERT_GT(global[1], PetscReal(0));
    EXPECT_LE(global[0], PetscReal(1e-9) * global[1]);
  }

  /**
   * @brief Expects a distributed device system to match the host one.
   * @param host Problem assembled on the host. @param device Problem
   * assembled on the device.
   */
  template <class Host, class Device>
  void expectSameSystem(Host& host, Device& device)
  {
    const auto& h = host.getLinearSystem();
    const auto& d = device.getLinearSystem();
    PetscBool kokkos = PETSC_FALSE;
    PetscObjectTypeCompareAny(reinterpret_cast<PetscObject>(d.getOperator()), &kokkos,
      MATSEQAIJKOKKOS, MATMPIAIJKOKKOS, "");
    EXPECT_TRUE(kokkos);
    PetscInt mh = 0, md = 0;
    MatGetLocalSize(h.getOperator(), &mh, nullptr);
    MatGetLocalSize(d.getOperator(), &md, nullptr);
    EXPECT_EQ(mh, md);
    PetscReal ah = 0, ad = 0, bh = 0, bd = 0;
    MatNorm(h.getOperator(), NORM_FROBENIUS, &ah);
    MatNorm(d.getOperator(), NORM_FROBENIUS, &ad);
    VecNorm(h.getVector(), NORM_2, &bh);
    VecNorm(d.getVector(), NORM_2, &bd);
    ASSERT_GT(ah, PetscReal(0));
    ASSERT_GT(bh, PetscReal(0));
    EXPECT_LE(std::abs(ah - ad), PetscReal(1e-12) * ah);
    EXPECT_LE(std::abs(bh - bd), PetscReal(1e-12) * bh);
  }

  /// @brief A Poisson problem with lifted Dirichlet values matches the host.
  void checkPoisson(Polytope::Type type)
  {
    using namespace Rodin::Variational;
    Context::MPI ctx(*g_env, *g_world);
    const HostMesh host = distributeFromRoot(ctx, type, 16);
    const DeviceMesh mesh(DeviceContext(ctx), host);
    P1<Real, HostMesh> hfes(host);
    P1<Real, DeviceMesh> dfes(mesh);
    EXPECT_EQ(hfes.getSize(), dfes.getSize());

    PETSc::Variational::TrialFunction hu(hfes);
    PETSc::Variational::TestFunction hv(hfes);
    PETSc::Variational::TrialFunction du(dfes);
    PETSc::Variational::TestFunction dv(dfes);

    RealFunction hf([](const Geometry::Point& p) { return 1 + p.x() * p.y(); });
    RealFunction hg([](const Geometry::Point& p) { return 2 * p.x() - p.y() + 1; });
    RealFunction df(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 1 + p.x() * p.y();
    });
    RealFunction dg(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 2 * p.x() - p.y() + 1;
    });

    Problem hp(hu, hv);
    hp = Integral(Grad(hu), Grad(hv)) + Integral(hu, hv) - Integral(hf, hv) +
      DirichletBC(hu, hg);
    Problem dp(du, dv);
    dp = Integral(Grad(du), Grad(dv)) + Integral(du, dv) - Integral(df, dv) +
      DirichletBC(du, dg);

    solve(hp);
    solve(dp);
    expectSameSystem(hp, dp);
    expectSameVector(hu.getSolution().getData(), du.getSolution().getData());

    // A second solve reuses the cached pattern and coordinates.
    solve(dp);
    expectSameVector(hu.getSolution().getData(), du.getSolution().getData());
  }

  TEST(Kokkos_MPI_Poisson, Triangles)
  {
    checkPoisson(Polytope::Type::Triangle);
  }

  TEST(Kokkos_MPI_Poisson, Quadrilaterals)
  {
    checkPoisson(Polytope::Type::Quadrilateral);
  }
}

/**
 * @brief Runs tests within the lifetimes of the MPI, Kokkos and PETSc runtimes.
 * @param argc Number of command-line arguments.
 * @param argv Command-line argument array.
 * @returns GoogleTest exit status.
 */
int main(int argc, char** argv)
{
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  g_env = &env;
  g_world = &world;
  ::testing::InitGoogleTest(&argc, argv);
  int res = 0;
  {
    ::Kokkos::ScopeGuard runtime(argc, argv);
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    res = RUN_ALL_TESTS();
    PetscFinalize();
  }
  return res;
}
