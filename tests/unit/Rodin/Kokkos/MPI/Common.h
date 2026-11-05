/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TESTS_UNIT_KOKKOS_MPI_COMMON_H
#define RODIN_TESTS_UNIT_KOKKOS_MPI_COMMON_H

/**
 * @file Common.h
 * @brief Meshes, solvers and comparisons shared by the distributed device
 * tests.
 */

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include <petsc.h>
#include <boost/filesystem.hpp>
#include <boost/mpi/collectives/all_gather.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/serialization/vector.hpp>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#include <Rodin/IO.h>
#include <Rodin/Solver.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/MPI/IO.h>
#include <Rodin/Kokkos/MPI.h>

#include "../TestContext.h"

#include <unistd.h>

namespace Rodin::Tests::Unit
{
  /// @brief MPI environment of the test executable.
  inline boost::mpi::environment* g_env = nullptr;

  /// @brief World communicator of the test executable.
  inline boost::mpi::communicator* g_world = nullptr;

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
  inline HostMesh distributeFromRoot(const Context::MPI& ctx, Polytope::Type type, size_t n)
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
   * @brief Solves a saddle point problem with a redundant LU factorization.
   *
   * Block Jacobi ILU, the default preconditioner, breaks down on the rows of
   * the multipliers once the system is split over several processes.
   */
  template <class Problem>
  void solveDirect(Problem& problem)
  {
    PETSc::Solver::GMRES solver(problem);
    ::PC pc;
    KSPGetPC(solver.getHandle(), &pc);
    PCSetType(pc, PCREDUNDANT);
    // A distributed device system lays out the rows of the multipliers among
    // those of each process, where the default nested dissection ordering of
    // the factorization meets a zero pivot.
    ::KSP inner;
    PCRedundantGetKSP(pc, &inner);
    ::PC lu;
    KSPGetPC(inner, &lu);
    PCSetType(lu, PCLU);
    PCFactorSetMatOrderingType(lu, MATORDERINGRCM);
    solver.solve();
  }

  /**
   * @brief Owned rows of a field in the host and device systems.
   *
   * The host numbers the rows field after field over all processes, and a
   * distributed device assembly of several fields lays out the owned rows of
   * a process field after field, so that a field owns the rows
   * @f$ [h, h + n) @f$ of the host system and @f$ [d, d + n) @f$ of the
   * device one on each process.
   */
  struct FieldRows
  {
    PetscInt host;   ///< First owned row in the host system.
    PetscInt device; ///< First owned row in the device system.
    PetscInt count;  ///< Number of owned rows.
  };

  /**
   * @brief Gets the owned rows of each field.
   * @param v Device vector. @param starts First row of each field on the
   * process in the device system, empty when both systems share the layout.
   */
  inline std::vector<FieldRows> getFieldRows(::Vec v, const std::vector<PetscInt>& starts)
  {
    PetscInt begin = 0, end = 0;
    VecGetOwnershipRange(v, &begin, &end);
    if (starts.empty())
      return { FieldRows{ begin, begin, end - begin } };
    std::vector<FieldRows> res;
    PetscInt offset = 0;
    for (size_t k = 0; k < starts.size(); ++k)
    {
      const PetscInt count = (k + 1 < starts.size() ? starts[k + 1] : end) - starts[k];
      PetscInt first = 0, total = 0;
      MPI_Exscan(&count, &first, 1, MPIU_INT, MPI_SUM, *g_world);
      if (g_world->rank() == 0)
        first = 0;
      MPI_Allreduce(&count, &total, 1, MPIU_INT, MPI_SUM, *g_world);
      res.push_back(FieldRows{ offset + first, starts[k], count });
      offset += total;
    }
    return res;
  }

  /**
   * @brief Calls a function on the subvectors of each field.
   * @param host Host vector. @param device Device vector. @param fields Owned
   * rows of each field. @param f Called with the host and device subvectors
   * and the first owned row of the field in the host system.
   */
  template <class F>
  void forEachField(::Vec host, ::Vec device, const std::vector<FieldRows>& fields, F&& f)
  {
    for (const auto& field : fields)
    {
      ::IS ih, id;
      ISCreateStride(*g_world, field.count, field.host, 1, &ih);
      ISCreateStride(*g_world, field.count, field.device, 1, &id);
      ::Vec sh, sd;
      VecGetSubVector(host, ih, &sh);
      VecGetSubVector(device, id, &sd);
      f(sh, sd, field.host);
      VecRestoreSubVector(host, ih, &sh);
      VecRestoreSubVector(device, id, &sd);
      ISDestroy(&ih);
      ISDestroy(&id);
    }
  }

  /**
   * @brief Expects two distributed vectors to agree.
   * @param host Reference vector. @param device Compared vector. @param
   * fields Owned rows of each field, see getFieldRows.
   */
  inline void expectSameVector(::Vec host, ::Vec device, const std::vector<FieldRows>& fields)
  {
    PetscInt nh = 0, nd = 0;
    VecGetLocalSize(host, &nh);
    VecGetLocalSize(device, &nd);
    ASSERT_EQ(nh, nd);
    PetscReal local[2] = {0, 0};
    forEachField(host, device, fields, [&](::Vec sh, ::Vec sd, PetscInt) {
      PetscInt n = 0;
      VecGetLocalSize(sh, &n);
      const PetscScalar* x = nullptr;
      const PetscScalar* y = nullptr;
      VecGetArrayRead(sh, &x);
      VecGetArrayRead(sd, &y);
      for (PetscInt i = 0; i < n; ++i)
      {
        local[0] = std::max(local[0], PetscAbsScalar(x[i] - y[i]));
        local[1] = std::max(local[1], PetscAbsScalar(x[i]));
      }
      VecRestoreArrayRead(sh, &x);
      VecRestoreArrayRead(sd, &y);
    });
    PetscReal global[2];
    MPI_Allreduce(local, global, 2, MPIU_REAL, MPI_MAX, *g_world);
    ASSERT_GT(global[1], PetscReal(0));
    EXPECT_LE(global[0], PetscReal(1e-9) * global[1]);
  }

  /**
   * @brief Expects two distributed vectors with the same layout to agree.
   * @param host Reference vector. @param device Compared vector.
   */
  inline void expectSameVector(::Vec host, ::Vec device)
  {
    expectSameVector(host, device, getFieldRows(device, {}));
  }

  /**
   * @brief Expects a distributed device operator to match the host one.
   * @param host Operator assembled on the host. @param device Operator
   * assembled on the device. @param fields Owned rows of each field of a
   * square system, see getFieldRows, or empty when the operators share their
   * layout.
   */
  inline void expectSameOperator(::Mat host, ::Mat device,
                                 const std::vector<FieldRows>& fields = {})
  {
    PetscBool kokkos = PETSC_FALSE;
    PetscObjectTypeCompareAny(reinterpret_cast<PetscObject>(device), &kokkos,
      MATSEQAIJKOKKOS, MATMPIAIJKOKKOS, "");
    EXPECT_TRUE(kokkos);
    PetscInt mh = 0, md = 0, nh = 0, nd = 0;
    MatGetLocalSize(host, &mh, &nh);
    MatGetLocalSize(device, &md, &nd);
    ASSERT_EQ(mh, md);
    ASSERT_EQ(nh, nd);
    PetscReal ah = 0, ad = 0;
    MatNorm(host, NORM_FROBENIUS, &ah);
    MatNorm(device, NORM_FROBENIUS, &ad);
    ASSERT_GT(ah, PetscReal(0));
    EXPECT_LE(std::abs(ah - ad), PetscReal(1e-12) * ah);

    // The norm does not see a permutation of the rows: the operators must
    // also agree on a vector that differs at each host row.
    ::Vec xh, yh, xd, yd;
    MatCreateVecs(host, &xh, &yh);
    MatCreateVecs(device, &xd, &yd);
    const auto columns = fields.empty() ? getFieldRows(xd, {}) : fields;
    const auto rows = fields.empty() ? getFieldRows(yd, {}) : fields;
    forEachField(xh, xd, columns, [](::Vec sh, ::Vec sd, PetscInt first) {
      PetscInt n = 0;
      VecGetLocalSize(sh, &n);
      PetscScalar* x = nullptr;
      PetscScalar* y = nullptr;
      VecGetArray(sh, &x);
      VecGetArray(sd, &y);
      for (PetscInt i = 0; i < n; ++i)
        x[i] = y[i] = std::sin(PetscReal(1 + first + i));
      VecRestoreArray(sh, &x);
      VecRestoreArray(sd, &y);
    });
    MatMult(host, xh, yh);
    MatMult(device, xd, yd);
    expectSameVector(yh, yd, rows);
    VecDestroy(&xh);
    VecDestroy(&yh);
    VecDestroy(&xd);
    VecDestroy(&yd);
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
    const auto fields = getFieldRows(d.getVector(), d.getFieldStarts());
    expectSameOperator(h.getOperator(), d.getOperator(), fields);
    expectSameVector(h.getVector(), d.getVector(), fields);
  }
}

#endif
