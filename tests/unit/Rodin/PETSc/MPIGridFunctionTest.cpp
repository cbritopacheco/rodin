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
#include <limits>

#include <Rodin/Geometry.h>
#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#include <Rodin/Variational.h>
#include <Rodin/MPI/Context/MPI.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/Variational.h>
#include <Rodin/PETSc/Variational/GridFunction.h>
#include <Rodin/Serialization/Vector.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

static boost::mpi::environment*  g_env   = nullptr;
static boost::mpi::communicator* g_world = nullptr;

namespace
{
  /// @brief Partitions over nonzero ranks so every geometry has an empty root.
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

  Mesh<Context::MPI> distributeWithEmptyRoot(const Context::MPI& ctx, Polytope::Type type)
  {
    const auto& comm = ctx.getCommunicator();
    Sharder<Context::MPI> sharder(ctx);
    if (comm.rank() == 0)
    {
      const size_t dim = Polytope::Traits(type).getDimension();
      auto local = dim == 1 ? Mesh<Context::Local>::UniformGrid(type, {3})
        : dim == 2          ? Mesh<Context::Local>::UniformGrid(type, {3, 3})
                            : Mesh<Context::Local>::UniformGrid(type, {3, 3, 3});
      local.getConnectivity().compute(dim, dim);
      for (size_t d = 1; d <= dim; ++d)
      {
        local.getConnectivity().compute(dim, d - 1);
        local.getConnectivity().compute(d, d - 1);
      }
      EmptyRootPartitioner partitioner(local, static_cast<size_t>(comm.size()));
      sharder.shard(partitioner);
      sharder.scatter(0);
    }
    return sharder.gather(0);
  }

  /**
   * @brief Checks owner and ghost coefficients against a global source oracle.
   *
   * Each owned eligible cell contributes its actual finite-element functionals
   * to an independent all-gather. For each global DOF, the reference chooses
   * the smallest distributed cell index. The expression depends on that index,
   * so competing incident-cell sources produce distinguishable coefficients.
   * No coordinate matching or nodal interpretation of H1 coefficients is used.
   * Full, odd-cell and empty selections also check preservation outside the
   * selected region. The const read is deliberately kept open before project.
   */
  template <class FES, class Function>
  void checkInterpolation(const FES& fes, size_t dim, const Function& function)
  {
    const auto& comm = *g_world;
    const auto& shard = fes.getMesh().getShard();
    Rodin::PETSc::Variational::GridFunction field(fes);
    const auto& read = field;
    const Index none = std::numeric_limits<Index>::max();
    std::vector<PetscScalar> aliasExpected(fes.getSize(), PetscScalar(-7));
    for (size_t filter = 0; filter < 3; ++filter)
    {
      SCOPED_TRACE(::testing::Message() << "filter=" << filter);
      const auto eligible = [&](const Polytope& cell) {
        const Index gid = shard.getPolytopeMap(dim).left.at(cell.getIndex());
        return filter == 0 || (filter == 1 && gid % 2 == 1);
      };
      std::vector<Index> localIndices;
      std::vector<PetscScalar> localValues;
      for (Index cell = 0; cell < shard.getPolytopeCount(dim); ++cell)
      {
        if (!shard.isOwned(dim, cell) || !eligible(*fes.getMesh().getPolytope(dim, cell)))
          continue;
        const auto& fe = fes.getFiniteElement(dim, cell);
        const auto pullback = fes.getPullback({dim, cell}, function);
        for (Index local = 0; local < fe.getCount(); ++local)
        {
          localIndices.push_back(fes.getGlobalIndex({dim, cell}, local));
          localIndices.push_back(shard.getPolytopeMap(dim).left.at(cell));
          localValues.push_back(fe.getLinearForm(local)(pullback));
        }
      }
      std::vector<std::vector<Index>> allIndices;
      std::vector<std::vector<PetscScalar>> allValues;
      boost::mpi::all_gather(comm, localIndices, allIndices);
      boost::mpi::all_gather(comm, localValues, allValues);
      std::vector<Index> source(fes.getSize(), none);
      std::vector<PetscScalar> expected(fes.getSize(), PetscScalar(-7));
      for (size_t rank = 0; rank < allIndices.size(); ++rank)
      {
        for (size_t i = 0; i < allValues[rank].size(); ++i)
        {
          const Index dof = allIndices[rank][2 * i];
          const Index cell = allIndices[rank][2 * i + 1];
          if (cell < source[dof])
          {
            source[dof] = cell;
            expected[dof] = allValues[rank][i];
          }
        }
      }
      if (filter == 0)
        for (Index dof = 0; dof < fes.getSize(); ++dof)
        {
          if (source[dof] != none)
            aliasExpected[dof] = PetscScalar(-14);
        }

      field = PetscScalar(-7);
      if (fes.getShard().getSize() > 0)
        EXPECT_EQ(read[fes.getGlobalIndex(0)], PetscScalar(-7));
      field.project(Region::Cells, function, eligible);
      Index begin = 0;
      Index end = 0;
      fes.getOwnershipRange(begin, end);
      for (Index dof = begin; dof < end; ++dof)
        EXPECT_EQ(read[dof], expected[dof]) << "owned global DOF=" << dof;
      for (Index local = 0; local < fes.getShard().getSize(); ++local)
      {
        const Index dof = fes.getGlobalIndex(local);
        EXPECT_EQ(read[dof], expected[dof]) << "represented global DOF=" << dof;
      }
      read.flush();
    }

    // Every source coefficient is read before any destination value changes.
    // This checks aliasing through an expression, including higher-order
    // functionals whose evaluation reads several coefficients at once.
    field.project(Region::Cells, field + field, [](const auto&) { return true; });
    for (Index local = 0; local < fes.getShard().getSize(); ++local)
    {
      const Index dof = fes.getGlobalIndex(local);
      EXPECT_LE(PetscAbsScalar(read[dof] - aliasExpected[dof]), Real(1e-10));
    }
    read.flush();
  }

  void checkInterpolationSpaces(const Mesh<Context::MPI>& mesh, size_t dim)
  {
    const auto& shard = mesh.getShard();
#if defined(PETSC_USE_COMPLEX)
    ComplexFunction scalar([&](const Point& point) {
#else
    RealFunction scalar([&](const Point& point) {
#endif
      const Index gid = shard.getPolytopeMap(dim).left.at(point.getPolytope().getIndex());
      Real value = Real(2) + Real(gid);
      const auto& coordinates = point.getPhysicalCoordinates();
      for (size_t d = 0; d < coordinates.size(); ++d)
        value += coordinates(d);
#if defined(PETSC_USE_COMPLEX)
      return Complex(value, Real(2) * value);
#else
      return value;
#endif
    });
    const auto vector = VectorFunction(scalar, Real(2) * scalar);
    const auto checkScalar = [&](const auto& fes, const char* name) {
      SCOPED_TRACE(name);
      checkInterpolation(fes, dim, scalar);
    };
    const auto checkVector = [&](const auto& fes, const char* name) {
      SCOPED_TRACE(name);
      checkInterpolation(fes, dim, vector);
    };
    using Scalar = PetscScalar;
    using Vector = Math::SpatialVector<Scalar>;
    using MeshType = Mesh<Context::MPI>;
    checkScalar(P0<Scalar, MeshType>(mesh), "P0 scalar");
    checkVector(P0<Vector, MeshType>(mesh, size_t(2)), "P0 vector");
    checkScalar(P0g<Scalar, MeshType>(mesh), "P0g scalar");
    checkVector(P0g<Vector, MeshType>(mesh, size_t(2)), "P0g vector");
    checkScalar(P1<Scalar, MeshType>(mesh), "P1 scalar");
    checkVector(P1<Vector, MeshType>(mesh, size_t(2)), "P1 vector");
    checkScalar(
      H1<1, Scalar, MeshType>(std::integral_constant<size_t, 1>{}, mesh), "H1 K1 scalar");
    checkVector(
      H1<1, Vector, MeshType>(std::integral_constant<size_t, 1>{}, mesh, size_t(2)),
      "H1 K1 vector");
    checkScalar(
      H1<2, Scalar, MeshType>(std::integral_constant<size_t, 2>{}, mesh), "H1 K2 scalar");
    checkVector(
      H1<2, Vector, MeshType>(std::integral_constant<size_t, 2>{}, mesh, size_t(2)),
      "H1 K2 vector");
    checkScalar(
      H1<3, Scalar, MeshType>(std::integral_constant<size_t, 3>{}, mesh), "H1 K3 scalar");
    checkVector(
      H1<3, Vector, MeshType>(std::integral_constant<size_t, 3>{}, mesh, size_t(2)),
      "H1 K3 vector");
  }

  class PETSc_MPI_Interpolation : public ::testing::TestWithParam<Polytope::Type>
  {};

  TEST_P(PETSc_MPI_Interpolation, SourcesAndStorageAcrossSpaces)
  {
    Context::MPI ctx(*g_env, *g_world);
    auto mesh = distributeWithEmptyRoot(ctx, GetParam());
    if (g_world->size() > 1 && g_world->rank() == 0)
      EXPECT_EQ(mesh.getShard().getVertexCount(), 0);
    checkInterpolationSpaces(mesh, Polytope::Traits(GetParam()).getDimension());
  }

  /**
   * @brief Zero-dimensional interpolation distinguishes actual point entities
   *        from an empty domain, including PETSc vectors with no owned entries.
   *
   * The point belongs only to the last rank, leaving the P0g coefficient owner
   * without any source entity whenever more than one rank participates.
   */
  TEST(PETSc_MPI_GridFunction, EmptyAndPointMeshInterpolation)
  {
    Context::MPI ctx(*g_env, *g_world);
    for (const bool includePoint : {false, true})
    {
      SCOPED_TRACE(includePoint);
      Shard::Builder builder;
      builder.initialize(0, 3);
      if (includePoint && g_world->rank() == g_world->size() - 1)
        builder.vertex(0, Math::SpatialPoint{0, 0, 0}, Shard::State::Owned);
      auto mesh =
        Mesh<Context::MPI>::Builder(ctx).initialize(builder.finalize()).finalize();
      EXPECT_EQ(mesh.getDimension(), 0);
      EXPECT_EQ(mesh.getShard().getVertexCount(),
        size_t(includePoint && g_world->rank() == g_world->size() - 1));
      checkInterpolationSpaces(mesh, 0);
    }
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, PETSc_MPI_Interpolation,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge));

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

  /// @brief Verifies that sync releases pending read/write access and refreshes
  ///        ghosts before direct PETSc use and point evaluation.
  TEST(PETSc_MPI_GridFunction, SyncReleasesAccessAndRefreshesGhosts)
  {
    auto& world = *g_world;
    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx);
    P1 fes(mesh);
    Rodin::PETSc::Variational::GridFunction gf(fes);

    writeShardDOFs(mesh, fes, gf, [](Index) { return static_cast<PetscScalar>(5.0); });
    ASSERT_TRUE(gf.getArrayWrite().acquired);
    gf.sync();
    EXPECT_FALSE(gf.getArrayWrite().acquired);

    size_t evaluated = 0;
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      const Point p(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
      EXPECT_NEAR(static_cast<Real>(PetscRealPart(gf(p))), 5.0, 1e-14);
      ++evaluated;
    }
    EXPECT_GT(evaluated, 0);

    const auto& cgf = gf;
    Index begin = 0;
    Index end = 0;
    fes.getOwnershipRange(begin, end);
    ASSERT_LT(begin, end);
    EXPECT_DOUBLE_EQ(static_cast<double>(PetscRealPart(cgf[begin])), 5.0);
    ASSERT_TRUE(cgf.getArrayRead().acquired);
    gf.sync();
    EXPECT_FALSE(cgf.getArrayRead().acquired);

    PetscReal norm = 0.0;
    PetscErrorCode ierr = VecNorm(gf.getData(), NORM_INFINITY, &norm);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
    EXPECT_DOUBLE_EQ(norm, 5.0);
    world.barrier();
  }

  /// @brief Verifies that reductions synchronize pending shard writes without
  ///        requiring callers to flush first.
  TEST(PETSc_MPI_GridFunction, ReductionsSynchronizePendingWrites)
  {
    auto& world = *g_world;
    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx);
    P1 fes(mesh);
    Rodin::PETSc::Variational::GridFunction gf(fes);

    const auto dofValue = [](Index i) {
      return static_cast<PetscScalar>(static_cast<Real>(i) - 4.0);
    };
    writeShardDOFs(mesh, fes, gf, dofValue);

    Real localSquares = 0.0;
    Real localAbs = 0.0;
    Index begin = 0;
    Index end = 0;
    fes.getOwnershipRange(begin, end);
    ASSERT_LT(begin, end);
    for (Index i = begin; i < end; ++i)
    {
      const Real value = static_cast<Real>(i) - 4.0;
      localSquares += value * value;
      localAbs += std::abs(value);
    }

    Real globalSquares = 0.0;
    Real globalAbs = 0.0;
    boost::mpi::all_reduce(world, localSquares, globalSquares, std::plus<Real>());
    boost::mpi::all_reduce(world, localAbs, globalAbs, std::plus<Real>());

    EXPECT_NEAR(gf.norm(), std::sqrt(globalSquares), 1e-10 * std::sqrt(globalSquares));
    EXPECT_NEAR(gf.norm(NORM_1), globalAbs, 1e-10 * globalAbs);

    Index idx = fes.getSize();
    EXPECT_DOUBLE_EQ(static_cast<double>(PetscRealPart(gf.min(idx))), -4.0);
    EXPECT_EQ(idx, 0);

    const Index maxIdx = static_cast<Index>(fes.getSize() - 1);
    EXPECT_DOUBLE_EQ(
      static_cast<double>(PetscRealPart(gf.max(idx))), static_cast<double>(maxIdx) - 4.0);
    EXPECT_EQ(idx, maxIdx);
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
