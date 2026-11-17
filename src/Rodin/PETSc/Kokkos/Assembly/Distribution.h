/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_KOKKOS_ASSEMBLY_DISTRIBUTION_H
#define RODIN_PETSC_KOKKOS_ASSEMBLY_DISTRIBUTION_H

/**
 * @file
 * @brief Translation from finite-element numbering to PETSc row ownership.
 */

#include <algorithm>
#include <cassert>
#include <unordered_map>
#include <utility>
#include <vector>
#include <petscsys.h>
#include "Rodin/Types.h"

namespace Rodin::Assembly::Internal
{
  /**
   * @brief Distribution of the degrees of freedom of a system over MPI
   * processes.
   *
   * The local degrees of freedom, ghosts included, are those the device
   * integrators number. Each one has a global index and is owned by one
   * process, the one owning it in its space. The process holds the
   * contiguous range @f$ [b, e) @f$ of the rows of the PETSc layout, whose
   * length is its number of owned degrees of freedom. With a single field
   * this range is the owned range of the space.
   *
   * With several fields, a degree of freedom is known by its index in the
   * system as on the host, the offset of its field plus its index in the
   * space of the field, and the constraints and the preassembled operators
   * are resolved in this numbering. The range of the process holds its
   * owned degrees of freedom of the first field, then those of the second
   * one, and so on, so that the rows of a process are the degrees of
   * freedom of its part of the mesh. getPETScIndex() maps the index in the
   * system to the row of the PETSc layout, and getFieldStarts() gives the
   * first row of each field on the process.
   */
  struct COODistribution
  {
    /// @brief Communicator of the processes.
    MPI_Comm comm = MPI_COMM_NULL;

    /// @brief First row of this process in the PETSc layout.
    PetscInt begin = 0;

    /// @brief One past the last row of this process in the PETSc layout.
    PetscInt end = 0;

    /// @brief Global number of degrees of freedom.
    PetscInt size = 0;

    /// @brief Global index of each local degree of freedom.
    std::vector<PetscInt> globals;

    /// @brief Whether this process owns each local degree of freedom.
    std::vector<unsigned char> owned;

    /// @brief Local index of the first degree of freedom of each field.
    std::vector<Index> starts;

    /// @brief Global index of the first degree of freedom of each field.
    std::vector<PetscInt> offsets;

    /// @brief First owned index and owned count of each field in its space.
    std::vector<std::pair<PetscInt, PetscInt>> ranges;

    /**
     * @brief First owned index in its space of each field on each process,
     * indexed by field then rank; empty when the PETSc layout is the numbering
     * of the system.
     */
    std::vector<std::vector<PetscInt>> firsts;

    /**
     * @brief First row in the PETSc layout of each field on each process,
     * indexed by field then rank.
     */
    std::vector<std::vector<PetscInt>> rows;

    /**
     * @brief Whether the local degrees of freedom have the same rows in the
     * PETSc layout as those of another distribution.
     * @param other Distribution compared.
     * @return Whether the communicators, the global and local degrees of
     * freedom and the layouts of the fields are the same.
     */
    bool hasLayout(const COODistribution& other) const
    {
      return comm == other.comm && begin == other.begin && end == other.end &&
        size == other.size && globals == other.globals && offsets == other.offsets &&
        firsts == other.firsts && rows == other.rows;
    }

    /// @brief Gets the number of owned degrees of freedom. @return Owned count.
    PetscInt getOwnedSize() const
    {
      return end - begin;
    }

    /**
     * @brief Appends the local degrees of freedom of a field.
     * @param fes Space numbering the local degrees of freedom of the field.
     * @param offset Global index of the first degree of freedom of the field.
     */
    template <class FES>
    void append(const FES& fes, PetscInt offset)
    {
      Index first = 0, last = 0;
      fes.getOwnershipRange(first, last);
      const size_t count = fes.getLocalSize();
      starts.push_back(static_cast<Index>(globals.size()));
      offsets.push_back(offset);
      ranges.emplace_back(
        static_cast<PetscInt>(first), static_cast<PetscInt>(last - first));
      globals.reserve(globals.size() + count);
      owned.reserve(owned.size() + count);
      for (size_t i = 0; i < count; ++i)
      {
        const Index g = fes.getGlobalIndex(static_cast<Index>(i));
        globals.push_back(offset + static_cast<PetscInt>(g));
        owned.push_back(g >= first && g < last);
      }
    }

    /**
     * @brief Gets the global offset of a field.
     * @param start Local index of the first degree of freedom of the field.
     * @return Global index of the first degree of freedom of the field.
     */
    PetscInt getOffset(Index start) const
    {
      const auto it = std::find(starts.begin(), starts.end(), start);
      assert(it != starts.end());
      return offsets[it - starts.begin()];
    }

    /// @brief Gets the local index of each global index held. @return Map from global to local indices.
    std::unordered_map<PetscInt, Index> getLocals() const
    {
      std::unordered_map<PetscInt, Index> res;
      res.reserve(globals.size());
      for (size_t i = 0; i < globals.size(); ++i)
        res.emplace(globals[i], static_cast<Index>(i));
      return res;
    }

    /**
     * @brief Lays out the rows of the processes in rank order, and the owned
     * degrees of freedom of a process field after field; collective.
     * @param communicator Communicator of the processes.
     * @param total Global number of degrees of freedom.
     */
    void close(MPI_Comm communicator, PetscInt total)
    {
      comm = communicator;
      size = total;
      const PetscInt count =
        static_cast<PetscInt>(std::count(owned.begin(), owned.end(), 1));
      PetscInt first = 0;
      MPI_Exscan(&count, &first, 1, MPIU_INT, MPI_SUM, comm);
      int rank = 0, processes = 1;
      MPI_Comm_rank(comm, &rank);
      MPI_Comm_size(comm, &processes);
      begin = rank == 0 ? 0 : first;
      end = begin + count;

      // First owned index and owned count of each field on each process.
      const size_t fields = ranges.size();
      std::vector<PetscInt> local(2 * fields);
      for (size_t k = 0; k < fields; ++k)
      {
        local[2 * k] = ranges[k].first;
        local[2 * k + 1] = ranges[k].second;
      }
      std::vector<PetscInt> all(local.size() * static_cast<size_t>(processes));
      MPI_Allgather(local.data(), static_cast<int>(local.size()), MPIU_INT, all.data(),
        static_cast<int>(local.size()), MPIU_INT, comm);
      firsts.assign(fields, std::vector<PetscInt>(static_cast<size_t>(processes)));
      rows.assign(fields, std::vector<PetscInt>(static_cast<size_t>(processes)));
      PetscInt row = 0;
      for (size_t r = 0; r < static_cast<size_t>(processes); ++r)
      {
        for (size_t k = 0; k < fields; ++k)
        {
          firsts[k][r] = all[r * local.size() + 2 * k];
          rows[k][r] = row;
          row += all[r * local.size() + 2 * k + 1];
        }
      }
      assert(rows.empty() || rows[0][static_cast<size_t>(rank)] == begin);
    }

    /**
     * @brief Gets the row of a degree of freedom in the PETSc layout.
     * @param global Index of the degree of freedom in the system.
     * @return Row of the degree of freedom.
     */
    PetscInt getPETScIndex(PetscInt global) const
    {
      if (firsts.empty())
        return global;
      const size_t k = static_cast<size_t>(
        std::upper_bound(offsets.begin(), offsets.end(), global) - offsets.begin() - 1);
      const PetscInt g = global - offsets[k];
      const auto& first = firsts[k];
      const size_t r = static_cast<size_t>(
        std::upper_bound(first.begin(), first.end(), g) - first.begin() - 1);
      return rows[k][r] + g - first[r];
    }

    /// @brief Gets the row in the PETSc layout of each local degree of freedom. @return Rows.
    std::vector<PetscInt> getPETScIndices() const
    {
      std::vector<PetscInt> res(globals.size());
      for (size_t i = 0; i < globals.size(); ++i)
        res[i] = getPETScIndex(globals[i]);
      return res;
    }

    /// @brief Gets the first row of each field on this process. @return Rows.
    std::vector<PetscInt> getFieldStarts() const
    {
      int rank = 0;
      MPI_Comm_rank(comm, &rank);
      std::vector<PetscInt> res(rows.size());
      for (size_t k = 0; k < rows.size(); ++k)
        res[k] = rows[k][static_cast<size_t>(rank)];
      return res;
    }

    /**
     * @brief Distribution of a distributed device space.
     * @param fes Space numbering the local degrees of freedom.
     * @return Distribution of the space.
     */
    template <class FES>
    static COODistribution of(const FES& fes)
    {
      COODistribution res;
      res.comm = fes.getMesh().getContext().getCommunicator();
      res.append(fes, 0);
      Index begin = 0, end = 0;
      fes.getOwnershipRange(begin, end);
      res.begin = static_cast<PetscInt>(begin);
      res.end = static_cast<PetscInt>(end);
      res.size = static_cast<PetscInt>(fes.getSize());
      return res;
    }
  };

}

#endif
