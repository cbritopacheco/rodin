/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Exact mixed PETSc target assembly with empty distributed shards. */
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#include <Rodin/MPI/Context/MPI.h>
#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/MPI/Variational/P1.h>
#include <Rodin/MPI/Variational/P0.h>
#include "MixedTargetAssembly.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  boost::mpi::environment* environment = nullptr;
  boost::mpi::communicator* world = nullptr;

  /// @brief Retains all communicator ranks while leaving rank zero empty in parallel.
  class EmptyRootPartitioner final : public Partitioner
  {
    public:
      /// @brief Assigns single-cell partitions among the available non-root ranks.
      EmptyRootPartitioner(const LocalMesh& mesh, size_t count)
        : m_base(mesh),
          m_count(count)
      {
        m_base.partition(1);
      }

      /// @brief Returns the source interval mesh.
      const LocalMesh& getMesh() const override
      {
        return m_base.getMesh();
      }

      /// @brief Updates the cell grouping before assigning communicator ranks.
      void partition(size_t maximumSize, size_t dimension) override
      {
        m_base.partition(maximumSize, dimension);
      }

      /// @brief Maps each cell group to an available non-root rank.
      size_t getPartition(Index index) const override
      {
        return m_count == 1 ? 0 : 1 + m_base.getPartition(index) % (m_count - 1);
      }

      /// @brief Includes ranks with no cells in the shard layout.
      size_t getCount() const override
      {
        return m_count;
      }

    private:
      BalancedCompactPartitioner m_base;
      size_t m_count;
  };
}

/// @brief Distributed targeted assembly matches exact blocks with nonzero boundary data.
template <bool Affine>
void checkDistributedMixedTargets()
{
  Context::MPI context(*environment, *world);
  Sharder<Context::MPI> sharder(context);
  if (world->rank() == 0)
  {
    auto local = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    local.scale(Real(0.5));
    local.getConnectivity().compute(0, 1);
    local.getConnectivity().compute(1, 1);
    for (auto vertex = local.getVertex(); vertex; ++vertex)
    {
      const Real x = local.getVertexCoordinates(vertex->getIndex()).x();
      if (x == 0)
        local.setAttribute(vertex.key(), 11);
      if (x == 1)
        local.setAttribute(vertex.key(), 12);
    }
    EmptyRootPartitioner partitioner(local, static_cast<size_t>(world->size()));
    sharder.shard(partitioner);
    sharder.scatter(0);
  }
  auto mesh = sharder.gather(0);
  if (world->size() > 1 && world->rank() == 0)
  {
    EXPECT_EQ(mesh.getShard().getCellCount(), 0u);
  }
  P1<PetscScalar, Mesh<Context::MPI>> space(mesh);
  std::array<Real, 3> coordinates{}, counts{};
  for (auto vertex = mesh.getShard().getVertex(); vertex; ++vertex)
  {
    const auto dof = space.getDOFs(0, vertex->getIndex())[0];
    ASSERT_LT(dof, coordinates.size());
    const Point point(*vertex, Math::SpatialPoint{});
    coordinates[dof] = point.x();
    counts[dof] = 1;
  }
  ASSERT_EQ(MPI_Allreduce(
              MPI_IN_PLACE, coordinates.data(), 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD),
    MPI_SUCCESS);
  ASSERT_EQ(
    MPI_Allreduce(MPI_IN_PLACE, counts.data(), 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD),
    MPI_SUCCESS);
  std::vector<std::pair<Real, Index>> positions;
  for (Index dof = 0; dof < coordinates.size(); ++dof)
  {
    ASSERT_GT(counts[dof], 0);
    positions.emplace_back(coordinates[dof] / counts[dof], dof);
  }
  std::sort(positions.begin(), positions.end());
  Rodin::Tests::Unit::checkMixedTargetAssembly<Assembly::MPI, Affine>(
    space, positions, MPI_COMM_WORLD);
}

/// @brief Distributed value constraints preserve exact targeted blocks with an empty root shard.
TEST(PETSc_MPIMixedTargetAssembly, EmptyShardsPreserveExactTargets)
{
  checkDistributedMixedTargets<false>();
}

/// @brief Distributed affine constraints obey conjugate-transpose projection, including ranks with no cells.
TEST(PETSc_MPIMixedTargetAssembly, EmptyShardsAffineIdentification)
{
  checkDistributedMixedTargets<true>();
}

int main(int argc, char** argv)
{
  boost::mpi::environment mpiEnvironment(argc, argv);
  boost::mpi::communicator mpiWorld;
  environment = &mpiEnvironment;
  world = &mpiWorld;
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  ::testing::InitGoogleTest(&argc, argv);
  const int local = RUN_ALL_TESTS();
  int result = 0;
  MPI_Allreduce(&local, &result, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  PetscFinalize();
  return result;
}

namespace
{
  /// @brief Checks each heterogeneous block against exact segment mass integrals.
  template <bool FirstConstant, bool SecondConstant>
  void checkHeterogeneousMassBlocks()
  {
    Context::MPI context(*environment, *world);
    Sharder<Context::MPI> sharder(context);
    if (world->rank() == 0)
    {
      auto local = LocalMesh::UniformGrid(Polytope::Type::Segment, {5});
      local.scale(Real(.25));
      local.getConnectivity().compute(0, 1);
      local.getConnectivity().compute(1, 1);
      BalancedCompactPartitioner partitioner(local);
      partitioner.partition(static_cast<size_t>(world->size()));
      sharder.shard(partitioner);
      sharder.scatter(0);
    }
    auto mesh = sharder.gather(0);
    using MeshType = decltype(mesh);
    using First = std::conditional_t<FirstConstant, P0<PetscScalar, MeshType>,
      P1<PetscScalar, MeshType>>;
    using Second = std::conditional_t<SecondConstant, P0<PetscScalar, MeshType>,
      P1<PetscScalar, MeshType>>;
    First first(mesh);
    Second second(mesh);
    PETSc::Variational::TrialFunction u(first);
    PETSc::Variational::TrialFunction p(second);
    PETSc::Variational::TestFunction v(first);
    PETSc::Variational::TestFunction q(second);
    Problem problem(u, v, p, q);
    problem = Integral(2 * u, v) + Integral(3 * p, v) + Integral(5 * u, q) +
      Integral(7 * p, q) - Integral(11 * v) - Integral(13 * q);
    const size_t firstSize = first.getSize(), secondSize = second.getSize(),
                 totalSize = firstSize + secondSize;
    Math::Matrix<PetscScalar> expected(totalSize, totalSize);
    expected.setZero();
    Math::Vector<PetscScalar> rhs(totalSize);
    rhs.setZero();
    auto block = [&](const auto& rows, const auto& columns, size_t rowOffset,
                   size_t columnOffset, Real coefficient, Real length) {
      for (Index i = 0; i < Index(rows.size()); ++i)
      {
        for (Index j = 0; j < Index(columns.size()); ++j)
        {
          Real weight = (rows.size() == 2 && columns.size() == 2)
            ? (i == j ? length / 3 : length / 6)
            : length / Real(rows.size() * columns.size());
          expected(rowOffset + rows[i], columnOffset + columns[j]) +=
            coefficient * weight;
        }
      }
    };
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      if (!mesh.getShard().isOwned(1, cell->getIndex()))
        continue;
      auto firstDOFs = first.getDOFs(1, cell->getIndex());
      auto secondDOFs = second.getDOFs(1, cell->getIndex());
      const Real length = cell->getMeasure();
      block(firstDOFs, firstDOFs, 0, 0, 2, length);
      block(firstDOFs, secondDOFs, 0, firstSize, 3, length);
      block(secondDOFs, firstDOFs, firstSize, 0, 5, length);
      block(secondDOFs, secondDOFs, firstSize, firstSize, 7, length);
      for (auto i : firstDOFs)
        rhs[i] += 11 * length / Real(firstDOFs.size());
      for (auto i : secondDOFs)
        rhs[firstSize + i] += 13 * length / Real(secondDOFs.size());
    }
    MPI_Allreduce(MPI_IN_PLACE, expected.data(), int(totalSize * totalSize), MPIU_SCALAR,
      MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(
      MPI_IN_PLACE, rhs.data(), int(totalSize), MPIU_SCALAR, MPI_SUM, MPI_COMM_WORLD);
    problem.assemble();
    problem.setFieldSplits();
    auto& system = problem.getLinearSystem();
    const auto matrix = system.getOperator();
    PetscInt ownershipBegin, ownershipEnd;
    ASSERT_EQ(
      MatGetOwnershipRange(matrix, &ownershipBegin, &ownershipEnd), PETSC_SUCCESS);
    const auto& splits = system.getFieldSplits().getSplits();
    ASSERT_EQ(splits.size(), 2u);
    for (size_t field = 0; field < splits.size(); ++field)
    {
      PetscInt globalSize, localSize;
      ASSERT_EQ(ISGetSize(splits[field].is, &globalSize), PETSC_SUCCESS);
      ASSERT_EQ(ISGetLocalSize(splits[field].is, &localSize), PETSC_SUCCESS);
      EXPECT_EQ(globalSize, PetscInt(field == 0 ? firstSize : secondSize));
      const PetscInt offset = field == 0 ? 0 : PetscInt(firstSize);
      const PetscInt* indices;
      ASSERT_EQ(ISGetIndices(splits[field].is, &indices), PETSC_SUCCESS);
      for (PetscInt i = 0; i < localSize; ++i)
      {
        EXPECT_GE(indices[i], ownershipBegin);
        EXPECT_LT(indices[i], ownershipEnd);
        EXPECT_GE(indices[i], offset);
        EXPECT_LT(indices[i], offset + globalSize);
      }
      ASSERT_EQ(ISRestoreIndices(splits[field].is, &indices), PETSC_SUCCESS);
    }
    const auto vector = system.getVector();
    Math::Matrix<PetscScalar> actual(totalSize, totalSize);
    actual.setZero();
    Math::Vector<PetscScalar> load(totalSize);
    load.setZero();
    PetscInt begin, end;
    ASSERT_EQ(MatGetOwnershipRange(matrix, &begin, &end), PETSC_SUCCESS);
    for (PetscInt i = begin; i < end; ++i)
    {
      for (PetscInt j = 0; j < PetscInt(totalSize); ++j)
      {
        PetscScalar value;
        ASSERT_EQ(MatGetValues(matrix, 1, &i, 1, &j, &value), PETSC_SUCCESS);
        actual(i, j) = value;
      }
    }
    ASSERT_EQ(VecGetOwnershipRange(vector, &begin, &end), PETSC_SUCCESS);
    for (PetscInt i = begin; i < end; ++i)
    {
      PetscScalar value;
      ASSERT_EQ(VecGetValues(vector, 1, &i, &value), PETSC_SUCCESS);
      load[i] = value;
    }
    MPI_Allreduce(MPI_IN_PLACE, actual.data(), int(totalSize * totalSize), MPIU_SCALAR,
      MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(
      MPI_IN_PLACE, load.data(), int(totalSize), MPIU_SCALAR, MPI_SUM, MPI_COMM_WORLD);
    const Real matrixError = (actual - expected).norm(), rhsError = (load - rhs).norm();
    EXPECT_LT(matrixError, 1e-12);
    EXPECT_LT(rhsError, 1e-12);
  }
}

/// @brief P1/P0 block assembly uses each space's own distributed DOF map.
TEST(PETSc_MPIMixedTargetAssembly, P1P0MassBlocks)
{
  checkHeterogeneousMassBlocks<false, true>();
}

/// @brief Reversing heterogeneous block order retains exact matrices and loads.
TEST(PETSc_MPIMixedTargetAssembly, P0P1MassBlocks)
{
  checkHeterogeneousMassBlocks<true, false>();
}

/**
 * @brief Local assembly preserves unequal component maps of the same space type.
 *
 * Only component zero is coupled. The independent segment mass matrix therefore
 * has zero rows/columns in all other components. Covers standalone forms, single
 * fields and mixed problems, whose index lookups share thread-local P1 storage.
 */
TEST(PETSc_MPIMixedTargetAssembly, UnequalVectorComponents)
{
  Context::MPI context(*environment, *world);
  Sharder<Context::MPI> sharder(context);
  if (world->rank() == 0)
  {
    auto local = LocalMesh::UniformGrid(Polytope::Type::Segment, {9});
    local.getConnectivity().compute(0, 1);
    local.getConnectivity().compute(1, 1);
    BalancedCompactPartitioner partitioner(local);
    partitioner.partition(static_cast<size_t>(world->size()));
    sharder.shard(partitioner);
    sharder.scatter(0);
  }
  auto mesh = sharder.gather(0);
  using MeshType = decltype(mesh);
  P1<Math::SpatialVector<PetscScalar>, MeshType> trial(mesh, 2), test(mesh, 3);
  P0<PetscScalar, MeshType> auxiliary(mesh);
  PETSc::Variational::TrialFunction u(trial);
  PETSc::Variational::TrialFunction p(auxiliary);
  PETSc::Variational::TestFunction v(test);
  PETSc::Variational::TestFunction q(auxiliary);
  Math::Matrix<PetscScalar> expected(test.getSize(), trial.getSize());
  expected.setZero();
  for (auto cell = mesh.getCell(); cell; ++cell)
  {
    if (!mesh.getShard().isOwned(cell->getDimension(), cell->getIndex()))
      continue;
    const auto& vertices = cell->getVertices();
    for (size_t a = 0; a < vertices.size(); ++a)
    {
      const auto rows = test.getDOFs(0, vertices[a]);
      for (size_t b = 0; b < vertices.size(); ++b)
      {
        const auto columns = trial.getDOFs(0, vertices[b]);
        expected(rows[0], columns[0]) +=
          cell->getMeasure() * Real(a == b ? 2 : 1) / Real(6);
      }
    }
  }
  ASSERT_EQ(MPI_Allreduce(MPI_IN_PLACE, expected.data(), int(expected.size()),
              MPIU_SCALAR, MPI_SUM, MPI_COMM_WORLD),
    MPI_SUCCESS);
  const auto check = [&](::Mat matrix) {
    PetscInt begin, end, columns;
    ASSERT_EQ(MatGetOwnershipRange(matrix, &begin, &end), PETSC_SUCCESS);
    ASSERT_EQ(MatGetSize(matrix, nullptr, &columns), PETSC_SUCCESS);
    for (PetscInt row = begin; row < end; ++row)
    {
      for (PetscInt column = 0; column < columns; ++column)
      {
        PetscScalar actual;
        ASSERT_EQ(MatGetValues(matrix, 1, &row, 1, &column, &actual), PETSC_SUCCESS);
        const PetscScalar value = row < expected.rows() && column < expected.cols()
          ? expected(row, column)
          : PetscScalar(0);
        EXPECT_LT(PetscAbsScalar(actual - value), 1e-12)
          << "entry=" << row << ',' << column;
      }
    }
  };
  BilinearForm form(u, v);
  form = Integral(Component(u, 0), Component(v, 0));
  form.assemble();
  check(form.getOperator());
  Problem single(u, v);
  single = Integral(Component(u, 0), Component(v, 0));
  single.assemble();
  check(single.getLinearSystem().getOperator());
  Problem mixed(u, v, p, q);
  mixed = Integral(Component(u, 0), Component(v, 0));
  mixed.assemble();
  check(mixed.getLinearSystem().getOperator());
}
