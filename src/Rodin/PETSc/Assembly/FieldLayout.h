/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_ASSEMBLY_FIELDLAYOUT_H
#define RODIN_PETSC_ASSEMBLY_FIELDLAYOUT_H

/**
 * @file
 * @brief Rows of the fields of a distributed system in the PETSc layout.
 */
#include <algorithm>
#include <cassert>
#include <vector>

#include <mpi.h>
#include <petscsys.h>

namespace Rodin::PETSc::Assembly
{
  /**
   * @brief Rows of the fields of a distributed system in the PETSc layout.
   *
   * A degree of freedom of a system of several fields is known by its index
   * in the system, the offset of its field plus its index in the space of the
   * field, and the constraints and the preassembled operators are resolved
   * in this numbering. PETSc gives each process a contiguous range of rows.
   * The range of a process holds its owned degrees of freedom of the first
   * field, then those of the second one, and so on, so that the rows of a
   * process are the degrees of freedom of its part of the mesh. With the
   * offsets alone, the rows of a process would not be its degrees of
   * freedom, and the halo of a product would be about as large as its local
   * part.
   */
  class FieldLayout
  {
    public:
      /**
       * @brief Field of the system.
       */
      struct Field
      {
        /// @brief Index in the system of the first degree of freedom of the field.
        PetscInt offset;

        /// @brief First owned index of this process in the space of the field.
        PetscInt first;

        /// @brief Number of degrees of freedom of the field owned by this process.
        PetscInt count;
      };

      /// @brief Layout without fields, in which an index is its row.
      FieldLayout() = default;

      /**
       * @brief Lays out the rows of the processes in rank order, and the
       * owned degrees of freedom of a process field after field; collective.
       * @param comm Communicator of the processes.
       * @param fields Fields of the system, by increasing offset.
       */
      FieldLayout(MPI_Comm comm, const std::vector<Field>& fields)
      {
        int rank = 0, processes = 1;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &processes);
        const size_t n = fields.size();
        std::vector<PetscInt> local(2 * n);
        for (size_t k = 0; k < n; ++k)
        {
          assert(k == 0 || fields[k - 1].offset <= fields[k].offset);
          m_offsets.push_back(fields[k].offset);
          local[2 * k] = fields[k].first;
          local[2 * k + 1] = fields[k].count;
        }
        std::vector<PetscInt> all(local.size() * static_cast<size_t>(processes));
        MPI_Allgather(local.data(), static_cast<int>(local.size()), MPIU_INT, all.data(),
          static_cast<int>(local.size()), MPIU_INT, comm);
        m_firsts.assign(n, std::vector<PetscInt>(static_cast<size_t>(processes)));
        m_rows.assign(n, std::vector<PetscInt>(static_cast<size_t>(processes)));
        PetscInt row = 0;
        for (size_t r = 0; r < static_cast<size_t>(processes); ++r)
        {
          for (size_t k = 0; k < n; ++k)
          {
            m_firsts[k][r] = all[r * local.size() + 2 * k];
            m_rows[k][r] = row;
            row += all[r * local.size() + 2 * k + 1];
          }
        }
        m_starts.resize(n);
        for (size_t k = 0; k < n; ++k)
          m_starts[k] = m_rows[k][static_cast<size_t>(rank)];
      }

      /**
       * @brief Gets the row of a degree of freedom.
       * @param index Index of the degree of freedom in the system.
       * @return Row of the degree of freedom in the PETSc layout.
       */
      PetscInt operator()(PetscInt index) const
      {
        if (m_offsets.empty())
          return index;
        const size_t k = static_cast<size_t>(
          std::upper_bound(m_offsets.begin(), m_offsets.end(), index) - m_offsets.begin() - 1);
        const PetscInt g = index - m_offsets[k];
        const auto& first = m_firsts[k];
        const size_t r = static_cast<size_t>(
          std::upper_bound(first.begin(), first.end(), g) - first.begin() - 1);
        return m_rows[k][r] + g - first[r];
      }

      /**
       * @brief Gets the first row of each field on this process.
       * @return First rows, empty for a layout without fields.
       */
      const std::vector<PetscInt>& getFieldStarts() const
      {
        return m_starts;
      }

    private:
      std::vector<PetscInt> m_offsets;
      std::vector<std::vector<PetscInt>> m_firsts;
      std::vector<std::vector<PetscInt>> m_rows;
      std::vector<PetscInt> m_starts;
  };
}

#endif
