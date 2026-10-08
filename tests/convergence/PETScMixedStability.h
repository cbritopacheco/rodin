/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TESTS_CONVERGENCE_PETSC_MIXED_STABILITY_H
#define RODIN_TESTS_CONVERGENCE_PETSC_MIXED_STABILITY_H

#include <limits>
#include <petscmat.h>
#include "MixedStability.h"

namespace Rodin::Tests::Convergence
{
  /** @brief Explicitly global PETSc transport for the finite pressure-spectrum oracle.
   * All matrix-communicator ranks participate. Complete sparse matrices and
   * logical essential-DOF indices and interpolated constant coefficients are
   * replicated solely for this small-mesh
   * verification operation; this is not a scalable solver or a mesh accessor.
   * Algebraic spectra are measured by the backend-independent MixedStability.
   * This adapter supports real PETSc and does not infer ownership from values.
   * @par Architecture
   * Integer boundary indices are united globally. Matrix rows and pressure
   * constant coefficients are collected in their global algebraic ordering.
   * The unchanged backend-independent oracle then eliminates the essential
   * indices, constructs the physical zero-mean pressure basis and compares
   * Schur eigenvalues with the independently whitened divergence spectrum.
   * This global spectral calculation is performed once, on communicator rank
   * zero. Its dimensions, eigenvalues and consistency defects are broadcast
   * to all participants; dense spectral workspaces are not replicated.
   * The caller synchronizes pending pressure-vector writes before entry.
   */
  class PETScMixedStability
  {
    public:
      static void compute(MixedStability::Result& result, Mat A, Mat B, Mat M,
        const IndexMap<Real>& constrained, Vec pressureConstant)
      {
        const MPI_Comm comm = PetscObjectComm(reinterpret_cast<PetscObject>(A));
        int ranks = 0, rank = 0;
        ASSERT_EQ(MPI_Comm_size(comm, &ranks), MPI_SUCCESS);
        ASSERT_EQ(MPI_Comm_rank(comm, &rank), MPI_SUCCESS);
        ASSERT_LE(constrained.size(), size_t(std::numeric_limits<int>::max()));
        std::vector<PetscInt> local;
        for (const auto& [index, value] : constrained)
        {
          ASSERT_LE(index, Index(std::numeric_limits<PetscInt>::max()));
          local.push_back(static_cast<PetscInt>(index));
        }
        int count = static_cast<int>(local.size());
        std::vector<int> counts(ranks), offsets(ranks);
        ASSERT_EQ(MPI_Allgather(&count, 1, MPI_INT, counts.data(), 1, MPI_INT, comm), MPI_SUCCESS);
        int total = 0;
        for (int rank = 0; rank < ranks; ++rank)
        {
          ASSERT_LE(counts[rank], std::numeric_limits<int>::max() - total);
          offsets[rank] = total;
          total += counts[rank];
        }
        std::vector<PetscInt> indices(total);
        ASSERT_EQ(MPI_Allgatherv(local.data(), count, MPIU_INT, indices.data(),
          counts.data(), offsets.data(), MPIU_INT, comm), MPI_SUCCESS);
        IndexMap<Real> global;
        for (PetscInt index : indices)
        {
          ASSERT_GE(index, 0);
          global.emplace(static_cast<Index>(index), 0);
        }
        Math::SparseMatrix<Real> energy, divergence, mass;
        copyGlobal(energy, A, ranks);
        copyGlobal(divergence, B, ranks);
        copyGlobal(mass, M, ranks);
        ASSERT_FALSE(::testing::Test::HasFatalFailure());
        // Collection is an explicit global operation, not a field accessor.
        // VecScatterCreateToAll copies owned coefficients in global ordering;
        // ghost entries do not become additional pressure coefficients.
        PetscInt size = 0;
        ASSERT_EQ(VecGetSize(pressureConstant, &size), PETSC_SUCCESS);
        ASSERT_EQ(size, mass.rows());
        Vec complete = nullptr;
        VecScatter scatter = nullptr;
        ASSERT_EQ(VecScatterCreateToAll(pressureConstant, &scatter, &complete), PETSC_SUCCESS);
        EXPECT_EQ(VecScatterBegin(scatter, pressureConstant, complete,
          INSERT_VALUES, SCATTER_FORWARD), PETSC_SUCCESS);
        EXPECT_EQ(VecScatterEnd(scatter, pressureConstant, complete,
          INSERT_VALUES, SCATTER_FORWARD), PETSC_SUCCESS);
        const PetscScalar* values = nullptr;
        EXPECT_EQ(VecGetArrayRead(complete, &values), PETSC_SUCCESS);
        Math::Vector<Real> constant(size);
        for (PetscInt i = 0; i < size; ++i) constant(i) = values[i];
        EXPECT_EQ(VecRestoreArrayRead(complete, &values), PETSC_SUCCESS);
        EXPECT_EQ(VecScatterDestroy(&scatter), PETSC_SUCCESS);
        EXPECT_EQ(VecDestroy(&complete), PETSC_SUCCESS);
        if (rank == 0)
          MixedStability::compute(result, energy, divergence, mass, global, constant);
        // A global oracle measures one global spectrum. Repeating the same
        // dense factorization/SVD on every participant adds no independent
        // evidence and multiplies its peak workspace by the rank count.
        int failed = rank == 0 && ::testing::Test::HasFatalFailure();
        ASSERT_EQ(MPI_Bcast(&failed, 1, MPI_INT, 0, comm), MPI_SUCCESS);
        ASSERT_EQ(failed, 0) << "Global pressure-spectrum computation failed";
        std::array<PetscInt, 2> dimensions{
          static_cast<PetscInt>(result.freeVelocity),
          static_cast<PetscInt>(result.zeroMeanPressure)};
        ASSERT_EQ(MPI_Bcast(dimensions.data(), dimensions.size(), MPIU_INT, 0, comm),
          MPI_SUCCESS);
        ASSERT_GE(dimensions[0], 0);
        ASSERT_GT(dimensions[1], 0);
        ASSERT_LE(dimensions[1], std::numeric_limits<int>::max());
        result.freeVelocity = static_cast<Index>(dimensions[0]);
        result.zeroMeanPressure = static_cast<Index>(dimensions[1]);
        result.eigenvalues.resize(dimensions[1]);
        static_assert(std::is_same_v<Real, PetscReal>);
        ASSERT_EQ(MPI_Bcast(result.eigenvalues.data(), static_cast<int>(dimensions[1]),
          MPIU_REAL, 0, comm), MPI_SUCCESS);
        std::array<Real, 4> defects{result.meanBasisDefect, result.eigenResidual,
          result.spectralDifference, result.velocityResidual};
        ASSERT_EQ(MPI_Bcast(defects.data(), defects.size(), MPIU_REAL, 0, comm), MPI_SUCCESS);
        result.meanBasisDefect = defects[0];
        result.eigenResidual = defects[1];
        result.spectralDifference = defects[2];
        result.velocityResidual = defects[3];
      }

    private:
      static void copyGlobal(Math::SparseMatrix<Real>& destination, Mat source, int ranks)
      {
        Mat complete = nullptr;
        if (ranks == 1)
          ASSERT_EQ(MatDuplicate(source, MAT_COPY_VALUES, &complete), PETSC_SUCCESS);
        else
          ASSERT_EQ(MatCreateRedundantMatrix(source, ranks, MPI_COMM_NULL,
            MAT_INITIAL_MATRIX, &complete), PETSC_SUCCESS);
        PetscInt rows = 0, columns = 0;
        EXPECT_EQ(MatGetSize(complete, &rows, &columns), PETSC_SUCCESS);
        PetscInt first = 0, last = 0;
        EXPECT_EQ(MatGetOwnershipRange(complete, &first, &last), PETSC_SUCCESS);
        EXPECT_EQ(first, 0);
        EXPECT_EQ(last, rows);
        std::vector<Eigen::Triplet<Real>> entries;
        for (PetscInt row = 0; row < rows; ++row)
        {
          PetscInt count = 0;
          const PetscInt* indices = nullptr;
          const PetscScalar* values = nullptr;
          EXPECT_EQ(MatGetRow(complete, row, &count, &indices, &values), PETSC_SUCCESS);
          for (PetscInt i = 0; i < count; ++i)
            entries.emplace_back(row, indices[i], values[i]);
          EXPECT_EQ(MatRestoreRow(complete, row, &count, &indices, &values), PETSC_SUCCESS);
        }
        EXPECT_EQ(MatDestroy(&complete), PETSC_SUCCESS);
        destination.resize(rows, columns);
        destination.setFromTriplets(entries.begin(), entries.end());
      }
  };
}
#endif
