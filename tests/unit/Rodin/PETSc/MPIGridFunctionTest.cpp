/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <petsc.h>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/operations.hpp>

#include <cassert>
#include <cmath>
#include <functional>

#include <Rodin/Geometry.h>
#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#include <Rodin/Variational.h>
#include <Rodin/MPI/Context/MPI.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/Variational/P1.h>
#include <Rodin/PETSc/Variational/GridFunction.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

static boost::mpi::environment*  g_env   = nullptr;
static boost::mpi::communicator* g_world = nullptr;

namespace
{
  /// @brief Returns the PETSc object id of a vector.
  ///
  /// Unlike the raw handle, an id is never reused by a later object, which
  /// makes it a reliable witness of whether a vector was reallocated.
  PetscObjectId objectId(const ::Vec& vec)
  {
    PetscObjectId id = 0;
    const PetscErrorCode ierr = PetscObjectGetId((PetscObject)vec, &id);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
    return id;
  }

  Mesh<Context::Local> makeShardableMesh()
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, { 5, 5 });
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);
    mesh.getConnectivity().compute(D, 0);
    return mesh;
  }

  Mesh<Context::MPI> distributeFromRoot(const Context::MPI& ctx)
  {
    const auto& comm = ctx.getCommunicator();
    Sharder<Context::MPI> sharder(ctx);

    if (comm.rank() == 0)
    {
      auto localMesh = makeShardableMesh();
      BalancedCompactPartitioner partitioner(localMesh);
      partitioner.partition(static_cast<size_t>(comm.size()));
      sharder.shard(partitioner);
      sharder.scatter(0);
    }

    return sharder.gather(0);
  }

  /// @brief Writes @p value into every DOF of the local shard, owned and ghost
  ///        alike, through `operator[]`.
  ///
  /// Writing only the owned range is not equivalent: `flush()` scatters ghost
  /// entries back to their owners with `INSERT_VALUES`, so untouched ghost
  /// slots would overwrite the values their owners just wrote.
  template <class FES, class GF, class Value>
  void writeShardDOFs(
    const Mesh<Context::MPI>& mesh, const FES& fes, GF& gf, const Value& value)
  {
    const size_t D = mesh.getDimension();
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      const auto& fe = fes.getFiniteElement(D, cell->getIndex());
      const auto& dofs = fes.getDOFs(D, cell->getIndex());
      for (size_t local = 0; local < fe.getCount(); ++local)
        gf[dofs[local]] = value(static_cast<Index>(dofs[local]));
    }
  }

  /// @brief Verifies rank filtered const read does not deadlock for PET sc MPI grid function by checking tolerance-based numerical results, MPI behavior.
  TEST(PETSc_MPI_GridFunction, RankFilteredConstReadDoesNotDeadlock)
  {
    auto& world = *g_world;
    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx);
    P1 fes(mesh);
    Rodin::PETSc::Variational::GridFunction gf(fes);

    gf = static_cast<PetscScalar>(3.0);

    if (world.rank() == 0)
    {
      Index begin = 0;
      Index end = 0;
      fes.getOwnershipRange(begin, end);
      ASSERT_LT(begin, end);

      const auto& cgf = gf;
      EXPECT_DOUBLE_EQ(static_cast<double>(PetscRealPart(cgf[begin])), 3.0);
      cgf.flush();
    }

    world.barrier();
    SUCCEED();
  }

  /// @brief Verifies rank filtered mutable access does not deadlock for PET sc MPI grid function by checking MPI behavior.
  TEST(PETSc_MPI_GridFunction, RankFilteredMutableAccessDoesNotDeadlock)
  {
    auto& world = *g_world;
    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx);
    P1 fes(mesh);
    Rodin::PETSc::Variational::GridFunction gf(fes);

    if (world.rank() == 0)
    {
      Index begin = 0;
      Index end = 0;
      fes.getOwnershipRange(begin, end);
      ASSERT_LT(begin, end);

      gf[begin] = static_cast<PetscScalar>(5.0);
    }

    world.barrier();
    SUCCEED();
  }

  TEST(PETSc_MPI_GridFunction, PointEvaluationUsesOwnedAndGhostCoefficients)
  {
    auto& world = *g_world;
    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx);
    P1 fes(mesh);
    P1 vectorFES(mesh, size_t(2));
    Rodin::PETSc::Variational::GridFunction gf(fes);
    Rodin::PETSc::Variational::GridFunction vector(vectorFES);
    gf = static_cast<PetscScalar>(3.25);
    vector = static_cast<PetscScalar>(-1.75);

    size_t evaluated = 0;
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      const Point p(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
      EXPECT_NEAR(static_cast<Real>(PetscRealPart(gf(p))), 3.25, 1e-14);
      const auto value = vector(p);
      EXPECT_NEAR(static_cast<Real>(PetscRealPart(value(0))), -1.75, 1e-14);
      EXPECT_NEAR(static_cast<Real>(PetscRealPart(value(1))), -1.75, 1e-14);
      ++evaluated;
    }
    EXPECT_GT(evaluated, 0);
    static_cast<const decltype(gf)&>(gf).flush();
    static_cast<const decltype(vector)&>(vector).flush();
    world.barrier();
  }

  /// @brief Verifies that axpy accumulates a scaled grid function on the owned
  ///        DOFs and refreshes the ghost layer, by checking exact expected
  ///        values and MPI behavior.
  TEST(PETSc_MPI_GridFunction, AxpyAccumulatesScaledGridFunctionAndRefreshesGhosts)
  {
    auto& world = *g_world;
    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx);
    P1 fes(mesh);
    Rodin::PETSc::Variational::GridFunction y(fes);
    Rodin::PETSc::Variational::GridFunction x(fes);

    y = static_cast<PetscScalar>(2.0);
    x = static_cast<PetscScalar>(3.0);

    y.axpy(static_cast<PetscScalar>(0.5), x);

    Index begin = 0;
    Index end = 0;
    fes.getOwnershipRange(begin, end);
    ASSERT_LT(begin, end);

    const auto& cy = y;
    for (Index i = begin; i < end; ++i)
      EXPECT_DOUBLE_EQ(static_cast<double>(PetscRealPart(cy[i])), 3.5);
    cy.flush();

    // Point evaluation reads owned and ghost coefficients alike, so a stale
    // ghost layer would show up here.
    size_t evaluated = 0;
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      const Point p(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
      EXPECT_NEAR(static_cast<Real>(PetscRealPart(y(p))), 3.5, 1e-14);
      ++evaluated;
    }
    EXPECT_GT(evaluated, 0);
    cy.flush();
    world.barrier();
  }

  /// @brief Verifies that axpy flushes DOFs written through operator[] and
  ///        propagates them to the ghost layer, by checking exact expected
  ///        values and MPI behavior.
  TEST(PETSc_MPI_GridFunction, AxpyFlushesPendingWritesBeforeGhostUpdate)
  {
    auto& world = *g_world;
    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx);
    P1 fes(mesh);
    Rodin::PETSc::Variational::GridFunction y(fes);
    Rodin::PETSc::Variational::GridFunction zero(fes);

    zero = static_cast<PetscScalar>(0.0);

    Index begin = 0;
    Index end = 0;
    fes.getOwnershipRange(begin, end);
    ASSERT_LT(begin, end);

    // Every rank writes the same value into every DOF of its shard -- owned
    // and ghost alike, as element-wise assembly does -- without flushing.
    // Adding zero must still land every DOF on that value: the write has to
    // be restored to the vector first, and the ghost layer refreshed
    // afterwards.
    writeShardDOFs(mesh, fes, y, [](Index) { return static_cast<PetscScalar>(7.0); });

    y.axpy(static_cast<PetscScalar>(1.0), zero);

    size_t evaluated = 0;
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      const Point p(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
      EXPECT_NEAR(static_cast<Real>(PetscRealPart(y(p))), 7.0, 1e-14);
      ++evaluated;
    }
    EXPECT_GT(evaluated, 0);
    static_cast<const decltype(y)&>(y).flush();
    world.barrier();
  }

  /// @brief Verifies that norm is a global reduction over the owned DOFs by
  ///        comparing against an independent MPI reduction.
  TEST(PETSc_MPI_GridFunction, NormReducesOverAllRanks)
  {
    auto& world = *g_world;
    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx);
    P1 fes(mesh);
    Rodin::PETSc::Variational::GridFunction gf(fes);

    Index begin = 0;
    Index end = 0;
    fes.getOwnershipRange(begin, end);
    ASSERT_LT(begin, end);

    // A value that depends only on the global DOF index, so every rank agrees
    // on the shared entries.
    const auto dofValue = [](Index i) {
      return static_cast<PetscScalar>(1.0 + static_cast<Real>(i));
    };
    writeShardDOFs(mesh, fes, gf, dofValue);
    gf.flush();

    Real localSquares = 0.0;
    Real localAbs = 0.0;
    Real localMax = 0.0;
    for (Index i = begin; i < end; ++i)
    {
      const Real value = 1.0 + static_cast<Real>(i);
      localSquares += value * value;
      localAbs += std::abs(value);
      localMax = std::max(localMax, std::abs(value));
    }

    Real globalSquares = 0.0;
    Real globalAbs = 0.0;
    Real globalMax = 0.0;
    boost::mpi::all_reduce(world, localSquares, globalSquares, std::plus<Real>());
    boost::mpi::all_reduce(world, localAbs, globalAbs, std::plus<Real>());
    boost::mpi::all_reduce(world, localMax, globalMax, boost::mpi::maximum<Real>());

    EXPECT_NEAR(gf.norm(), std::sqrt(globalSquares), 1e-10 * std::sqrt(globalSquares));
    EXPECT_NEAR(gf.norm(NORM_1), globalAbs, 1e-10 * globalAbs);
    EXPECT_NEAR(gf.norm(NORM_INFINITY), globalMax, 1e-10 * globalMax);
    world.barrier();
  }

  /// @brief Verifies that copy assignment reuses the destination PETSc vector
  ///        and leaves owned and ghost DOFs consistent, by checking exact
  ///        expected values and MPI behavior.
  TEST(PETSc_MPI_GridFunction, CopyAssignmentReusesVectorHandleAndSyncsGhosts)
  {
    auto& world = *g_world;
    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx);
    P1 fes(mesh);
    Rodin::PETSc::Variational::GridFunction source(fes);
    Rodin::PETSc::Variational::GridFunction destination(fes);

    Index begin = 0;
    Index end = 0;
    fes.getOwnershipRange(begin, end);
    ASSERT_LT(begin, end);

    writeShardDOFs(
      mesh, fes, source, [](Index) { return static_cast<PetscScalar>(-2.75); });
    source.flush();

    // Deliberately no explicit flush() on the destination: the assignment
    // must restore its array before overwriting the vector.
    writeShardDOFs(
      mesh, fes, destination, [](Index) { return static_cast<PetscScalar>(11.0); });

    const ::Vec before = destination.getData();
    const PetscObjectId beforeId = objectId(before);
    destination = source;
    EXPECT_EQ(before, destination.getData());
    EXPECT_EQ(beforeId, objectId(destination.getData()));
    EXPECT_NE(destination.getData(), source.getData());

    const auto& cDestination = destination;
    for (Index i = begin; i < end; ++i)
      EXPECT_DOUBLE_EQ(static_cast<double>(PetscRealPart(cDestination[i])), -2.75);
    cDestination.flush();

    size_t evaluated = 0;
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      const Point p(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
      EXPECT_NEAR(static_cast<Real>(PetscRealPart(destination(p))), -2.75, 1e-14);
      ++evaluated;
    }
    EXPECT_GT(evaluated, 0);
    cDestination.flush();
    world.barrier();
  }
}

int main(int argc, char** argv)
{
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  g_env = &env;
  g_world = &world;

  PetscInitialize(&argc, &argv, nullptr, nullptr);
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
