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
   * logical essential-DOF indices are replicated solely for this small-mesh
   * verification operation; this is not a scalable solver or a mesh accessor.
   * Algebraic spectra are measured by the backend-independent MixedStability.
   * This adapter supports real PETSc and does not infer ownership from values.
   */
  class PETScMixedStability
  {
    public:
      static void compute(MixedStability::Result& result, Mat A, Mat B, Mat M,
        const IndexMap<Real>& constrained)
      {
        const MPI_Comm comm = PetscObjectComm(reinterpret_cast<PetscObject>(A));
        int ranks = 0;
        ASSERT_EQ(MPI_Comm_size(comm, &ranks), MPI_SUCCESS);
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
        // Pressure is P1: its nodal coefficient vector for the constant one
        // is exactly one. No higher-order coefficient-layout claim is made.
        const Math::Vector<Real> constant = Math::Vector<Real>::Ones(mass.rows());
        MixedStability::compute(result, energy, divergence, mass, global, constant);
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
