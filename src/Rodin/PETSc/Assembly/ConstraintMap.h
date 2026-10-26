/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_ASSEMBLY_CONSTRAINTMAP_H
#define RODIN_PETSC_ASSEMBLY_CONSTRAINTMAP_H

/**
 * @file
 * @brief Collective closure of distributed affine identification constraints.
 */
#include <petscsys.h>
#include <boost/mpi/collectives/all_gather.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include "Rodin/Assembly/ConstraintMap.h"

namespace Rodin::PETSc::Assembly
{
  /**
   * @brief Affine constraint map with collective dependency closure.
   *
   * # Architecture
   *
   * Identification rows specify @f$x_s=\sum_m C_{sm}x_m+d_s@f$.
   * A master can itself be constrained on another rank. Every rank therefore
   * needs the complete constraint dependency graph before expanding a local
   * matrix entry. Only prescribed values and identification rows are exchanged;
   * matrix entries and unconstrained rows are not communicated here.
   *
   * # Implementation
   *
   * Two collective gathers publish fixed rows and affine identification rows.
   * Replicated boundary rows must agree, independently of master ordering.
   * The native @ref Rodin::Assembly::ConstraintMap then resolves transitive
   * dependencies, propagates affine offsets and rejects cycles. Every rank
   * receives the same graph, including ranks with no locally constrained DOFs.
   * Assemblers invoke this collective only for identification conditions;
   * ordinary value-prescribing boundaries retain their existing communication.
   */
  class ConstraintMap final : public Rodin::Assembly::ConstraintMap<PetscScalar>
  {
    public:
      /// @brief Backend-independent constraint algebra.
      using Parent = Rodin::Assembly::ConstraintMap<PetscScalar>;

      /**
       * @brief Creates a map indexed by global degrees of freedom.
       * @param count Global degree-of-freedom extent.
       */
      explicit ConstraintMap(size_t count)
        : Parent(count)
      {}

      /**
       * @brief Collectively resolves constraints across all communicator ranks.
       * @param communicator Communicator of the assembled linear system.
       * @throws Alert::Exception For inconsistent replicated rows or cycles.
       *
       * All ranks must call with the same global extent, before local
       * finalization or projection of any matrix/vector contributions.
       */
      void finalize(MPI_Comm communicator)
      {
        using Coefficient = std::pair<Index, PetscScalar>;
        using Identification =
          std::pair<Index, std::pair<PetscScalar, std::vector<Coefficient>>>;
        std::vector<Coefficient> fixed;
        std::vector<Identification> identified;
        for (Index row = 0; row < size(); ++row)
        {
          if (isFixed(row))
            fixed.emplace_back(row, getFixedValue(row));
        }
        for (const Index row : getIdentifiedRows())
        {
          std::vector<Coefficient> entries;
          entries.reserve(expand(row).size());
          for (const auto& entry : expand(row))
            entries.emplace_back(entry.index, entry.coefficient);
          identified.emplace_back(
            row, std::pair{getIdentificationValue(row), std::move(entries)});
        }
        const boost::mpi::communicator world(communicator, boost::mpi::comm_attach);
        std::vector<std::vector<Coefficient>> allFixed;
        std::vector<std::vector<Identification>> allIdentified;
        boost::mpi::all_gather(world, fixed, allFixed);
        boost::mpi::all_gather(world, identified, allIdentified);
        for (const auto& rows : allFixed)
        {
          for (const auto& [row, value] : rows)
          {
            if (isFixed(row) && getFixedValue(row) != value)
            {
              Alert::MemberFunctionException(*this, __func__)
                << "Inconsistent distributed prescribed value at DOF " << row << "."
                << Alert::Raise;
            }
            if (!isFixed(row))
              setFixed(row, value);
          }
        }
        for (const auto& rows : allIdentified)
        {
          for (const auto& [row, definition] : rows)
          {
            const auto& [value, coefficients] = definition;
            std::vector<Entry> entries;
            entries.reserve(coefficients.size());
            for (const auto& [master, coefficient] : coefficients)
              entries.push_back({master, coefficient});
            if (isIdentified(row))
            {
              const auto& previous = expand(row);
              bool matches =
                previous.size() == entries.size() && getIdentificationValue(row) == value;
              for (const auto& entry : entries)
              {
                const auto found = std::find_if(
                  previous.begin(), previous.end(), [&](const Entry& candidate) {
                    return candidate.index == entry.index &&
                      candidate.coefficient == entry.coefficient;
                  });
                matches = matches && found != previous.end();
              }
              if (!matches)
              {
                Alert::MemberFunctionException(*this, __func__)
                  << "Inconsistent distributed identification at DOF " << row << "."
                  << Alert::Raise;
              }
            }
            else
              setIdentification(row, entries, value);
          }
        }
        if (!getIdentifiedRows().empty())
          Parent::finalize();
      }
  };
}

#endif
