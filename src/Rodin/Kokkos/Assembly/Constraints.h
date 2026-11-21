/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_ASSEMBLY_CONSTRAINTS_H
#define RODIN_KOKKOS_ASSEMBLY_CONSTRAINTS_H

/**
 * @file
 * @brief Resolution and device expansion of Rodin assembly constraints.
 */

#include <algorithm>
#include <string>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

#include <Kokkos_Core.hpp>

#include "Rodin/Assembly/ConstraintMap.h"
#include "Rodin/Variational/DirichletBC.h"
#include "Rodin/Variational/PeriodicBC.h"
#include "Rodin/Kokkos/Context/Kokkos.h"

namespace Rodin::Assembly
{
  /**
   * @brief Constraints of a device system, resolved on the host.
   *
   * ## Architecture
   *
   * Value constraints are recorded as they are assembled. When an
   * identification Dirichlet condition or a periodic condition is present,
   * resolve() replays every condition, in order, in a host ConstraintMap,
   * as the host assemblies do, and uploads the expansion of every degree
   * of freedom on its masters.
   *
   * @tparam Backend Kokkos backend of the mesh context.
   */
  template <class Backend>
  class Constraints
  {
    public:
      /// @brief Mesh context.
      using ContextType = ::Rodin::Context::Kokkos<::Rodin::Context::Local, Backend>;

      /// @brief Execution space of the kernels.
      using ExecutionSpace = typename ContextType::ExecutionSpaceType;

      /// @brief Memory space of the expansions.
      using MemorySpace = typename ContextType::MemorySpaceType;

      /// @brief Device array.
      template <class T>
      using Array = ::Kokkos::View<T*, MemorySpace>;

      /// @brief Expansions of the degrees of freedom on their masters.
      struct Expansions
      {
        /// @brief Start of every expansion, followed by the entry count.
        Array<Index> starts;

        /// @brief Masters.
        Array<Index> indices;

        /// @brief Coefficients of the masters.
        Array<Real> coefficients;

        /**
         * @brief Affine values of the identified degrees of freedom, zero
         * for the others.
         */
        Array<Real> values;

        /**
         * @brief Number of masters of a degree of freedom.
         * @param i Index of the degree of freedom.
         * @return Number of entries of the expansion of @p i.
         */
        KOKKOS_INLINE_FUNCTION size_t count(Index i) const
        {
          return starts(i + 1) - starts(i);
        }
      };

      /// @brief Identification constraints resolved on the host.
      struct Resolved
      {
        /// @brief Expansions, empty when no degree of freedom stays identified.
        Expansions expansions;

        /**
         * @brief Rows of the constraints of the identified rows, by
         * increasing row.
         */
        std::vector<Index> rows;

        /// @brief Columns of the constraint rows.
        std::vector<Index> cols;

        /// @brief Values of the constraint rows.
        std::vector<Real> values;

        /**
         * @brief Whether the right-hand side of each row is imposed, the
         * rows being fixed or identified.
         */
        std::vector<unsigned char> imposed;

        /**
         * @brief Whether some degree of freedom stays identified.
         * @return True when the expansions are uploaded.
         */
        bool isExpanded() const
        {
          return expansions.starts.extent(0) > 0;
        }
      };

      /**
       * @brief Creates empty constraints.
       * @param rows Number of rows of the system.
       * @param size Number of degrees of freedom, the larger of the row
       * and column counts.
       */
      Constraints(size_t rows, size_t size)
        : m_rows(rows),
          m_fixed(size, 0),
          m_values(size, Real(0)),
          m_constrained(false),
          m_identified(false)
      {}

      /**
       * @brief Assembles a Dirichlet condition and records its constraints.
       * @param dbc Value or identification Dirichlet condition.
       * @param offset Offset of the constrained space in the system.
       * @param masterOffset Offset of the space of the masters of an
       * identification.
       */
      void fix(
        Variational::DirichletBCBase<Real>& dbc, size_t offset, size_t masterOffset)
      {
        using ValueDOFs = typename Variational::DirichletBCBase<Real>::ValueDOFs;
        dbc.assemble();
        m_conditions.push_back({&dbc, nullptr, offset, masterOffset});
        const auto* dofs = std::get_if<ValueDOFs>(&dbc.getDOFs());
        if (!dofs)
        {
          m_identified = true;
          return;
        }
        for (const auto& [local, value] : *dofs)
        {
          m_fixed[offset + local] = 1;
          m_values[offset + local] = value;
          m_constrained = true;
        }
      }

      /**
       * @brief Assembles a periodic condition and records its constraints.
       * @param pbc Periodic condition.
       * @param offset Offset of the constrained space in the system.
       */
      void identify(Variational::PeriodicBCBase<Real>& pbc, size_t offset)
      {
        pbc.assemble();
        m_conditions.push_back({nullptr, &pbc, offset, offset});
        m_identified = true;
      }

      /**
       * @brief Whether some degree of freedom is fixed.
       * @return True when at least one degree of freedom is fixed.
       */
      bool isConstrained() const
      {
        return m_constrained;
      }

      /**
       * @brief Whether some condition identifies degrees of freedom.
       * @return True when an identification or periodic condition was
       * recorded.
       */
      bool isIdentified() const
      {
        return m_identified;
      }

      /**
       * @brief Whether each degree of freedom is fixed.
       * @return Flag of each degree of freedom, nonzero when fixed.
       */
      const std::vector<unsigned char>& getFixed() const
      {
        return m_fixed;
      }

      /**
       * @brief Values of the fixed degrees of freedom, and of the
       * identified ones once resolved.
       * @return Value of each degree of freedom, zero when unconstrained.
       */
      const std::vector<Real>& getValues() const
      {
        return m_values;
      }

      /**
       * @brief Replays the conditions in a host constraint map and uploads
       * the expansions of the degrees of freedom identified after
       * flattening.
       *
       * The fixed degrees of freedom are those of the flattened map, which
       * may fix degrees of freedom identified on fixed masters.
       *
       * @param execution Execution space of the uploads.
       * @return Resolved constraints, with empty expansions when no degree
       * of freedom stays identified.
       */
      Resolved resolve(const ExecutionSpace& execution)
      {
        using DBC = Variational::DirichletBCBase<Real>;
        using Entry = typename ConstraintMap<Real>::Entry;
        ConstraintMap<Real> map(m_fixed.size());
        for (const auto& c : m_conditions)
        {
          if (c.pbc)
          {
            for (const auto& [slave, relation] : c.pbc->getDOFs())
            {
              std::vector<Entry> entries;
              for (Index k = 0; k < static_cast<Index>(relation.first.size()); ++k)
                entries.push_back({c.offset + relation.first[k], relation.second[k]});
              map.setIdentification(c.offset + slave, entries);
            }
            continue;
          }
          if (const auto* dofs = std::get_if<typename DBC::ValueDOFs>(&c.dbc->getDOFs()))
          {
            for (const auto& [local, value] : *dofs)
              map.setFixed(c.offset + local, value);
            continue;
          }
          const auto& affine = c.dbc->getIdentificationValues();
          for (const auto& [slave, relation] :
            std::get<typename DBC::IdentifiedDOFs>(c.dbc->getDOFs()))
          {
            std::vector<Entry> entries;
            for (Index k = 0; k < static_cast<Index>(relation.first.size()); ++k)
              entries.push_back({c.masterOffset + relation.first[k], relation.second[k]});
            const auto it = affine.find(slave);
            map.setIdentification(
              c.offset + slave, entries, it == affine.end() ? Real(0) : it->second);
          }
        }
        if (!map.getIdentifiedRows().empty())
          map.finalize();

        const size_t size = map.size();
        m_constrained = false;
        for (size_t i = 0; i < size; ++i)
        {
          const bool fixed = map.isFixed(i);
          m_fixed[i] = fixed;
          m_values[i] = fixed ? map.getFixedValue(i) : Real(0);
          m_constrained = m_constrained || fixed;
        }
        Resolved res;
        if (map.getIdentifiedRows().empty())
          return res;
        res.imposed.assign(m_fixed.begin(), m_fixed.begin() + m_rows);
        std::vector<Index> starts(size + 1, 0), indices;
        std::vector<Real> coefficients, values(size, Real(0));
        for (size_t i = 0; i < size; ++i)
        {
          for (const auto& e : map.expand(i))
          {
            indices.push_back(e.index);
            coefficients.push_back(e.coefficient);
          }
          starts[i + 1] = indices.size();
        }
        for (const Index gs : map.getIdentifiedRows())
        {
          values[gs] = map.getIdentificationValue(gs);
          // The value of an identified row is imposed with the fixed ones.
          m_values[gs] = values[gs];
          if (static_cast<size_t>(gs) >= m_rows)
            continue;
          res.imposed[gs] = 1;
          res.rows.push_back(gs);
          res.cols.push_back(gs);
          res.values.push_back(Real(1));
          for (const auto& e : map.expand(gs))
          {
            res.rows.push_back(gs);
            res.cols.push_back(e.index);
            res.values.push_back(-e.coefficient);
          }
        }
        res.expansions.starts = upload(execution, starts, "Rodin expansion starts");
        res.expansions.indices = upload(execution, indices, "Rodin expansion masters");
        res.expansions.coefficients =
          upload(execution, coefficients, "Rodin expansion coefficients");
        res.expansions.values = upload(execution, values, "Rodin expansion values");
        execution.fence("Rodin expansions");
        return res;
      }

      /**
       * @brief Resolves the conditions of a distributed system in the global
       * numbering, then expands the local degrees of freedom on local
       * indices.
       *
       * The periodic and identification conditions relate global indices
       * of their spaces, which may be held by other processes; each
       * process may give any of them, and the map gathers them all. The fixed degrees of freedom are
       * those already synchronized over the processes. The masters the
       * process does not hold are appended to the local degrees of
       * freedom, so that their entries are sent to their owners, and only
       * the owner of an identified degree of freedom writes its row.
       *
       * @param execution Execution space of the uploads.
       * @param map Empty constraint map of the global degrees of freedom.
       * @param globals Global index of each local degree of freedom,
       * appended with the masters the process does not hold.
       * @param owned Whether the process owns each local degree of freedom,
       * appended likewise.
       * @param offset Callable giving the global index of the first degree
       * of freedom of a field from its local one.
       * @param finalize Callable resolving the map over the processes;
       * collective.
       * @return Resolved constraints in local indices, with empty
       * expansions when no degree of freedom stays identified.
       */
      template <class Map, class Offset, class Finalize>
      Resolved resolve(const ExecutionSpace& execution, Map& map,
        std::vector<Index>& globals, std::vector<unsigned char>& owned, Offset&& offset,
        Finalize&& finalize)
      {
        using Entry = typename Map::Entry;
        const size_t local = globals.size();
        assert(owned.size() == local && m_fixed.size() >= local);
        for (size_t i = 0; i < local; ++i)
        {
          if (m_fixed[i] && owned[i])
            map.setFixed(globals[i], m_values[i]);
        }
        using DBC = Variational::DirichletBCBase<Real>;
        for (const auto& c : m_conditions)
        {
          const Index first = offset(c.offset);
          if (c.dbc)
          {
            // The value conditions are already among the fixed degrees of
            // freedom.
            const auto* dofs =
              std::get_if<typename DBC::IdentifiedDOFs>(&c.dbc->getDOFs());
            if (!dofs)
              continue;
            const Index masters = offset(c.masterOffset);
            const auto& affine = c.dbc->getIdentificationValues();
            for (const auto& [slave, relation] : *dofs)
            {
              std::vector<Entry> entries;
              for (Index k = 0; k < static_cast<Index>(relation.first.size()); ++k)
                entries.push_back({masters + relation.first[k], relation.second[k]});
              const auto it = affine.find(slave);
              map.setIdentification(
                first + slave, entries, it == affine.end() ? Real(0) : it->second);
            }
            continue;
          }
          for (const auto& [slave, relation] : c.pbc->getDOFs())
          {
            std::vector<Entry> entries;
            for (Index k = 0; k < static_cast<Index>(relation.first.size()); ++k)
              entries.push_back({first + relation.first[k], relation.second[k]});
            map.setIdentification(first + slave, entries);
          }
        }
        finalize(map);

        std::unordered_map<Index, Index> locals;
        locals.reserve(local);
        for (size_t i = 0; i < local; ++i)
          locals.emplace(globals[i], static_cast<Index>(i));
        const auto find = [&](Index global) {
          const auto [it, inserted] =
            locals.emplace(global, static_cast<Index>(globals.size()));
          if (inserted)
          {
            globals.push_back(global);
            owned.push_back(0);
          }
          return it->second;
        };
        Resolved res;
        std::vector<Index> starts(1, 0), indices;
        std::vector<Real> coefficients;
        std::vector<std::vector<std::pair<Index, Real>>> rows(local);
        for (size_t i = 0; i < local; ++i)
        {
          for (const auto& e : map.expand(globals[i]))
          {
            const Index master = find(e.index);
            indices.push_back(master);
            coefficients.push_back(e.coefficient);
            if (map.isIdentified(globals[i]))
              rows[i].emplace_back(master, e.coefficient);
          }
          starts.push_back(indices.size());
        }
        // The appended masters, never identified once flattened, expand on
        // themselves.
        const size_t size = globals.size();
        for (size_t i = local; i < size; ++i)
        {
          indices.push_back(static_cast<Index>(i));
          coefficients.push_back(Real(1));
          starts.push_back(indices.size());
        }
        m_rows = size;
        m_fixed.assign(size, 0);
        m_values.assign(size, Real(0));
        m_constrained = false;
        std::vector<Real> values(size, Real(0));
        res.imposed.assign(size, 0);
        for (size_t i = 0; i < local; ++i)
        {
          const Index g = globals[i];
          if (map.isFixed(g))
          {
            m_fixed[i] = 1;
            m_values[i] = map.getFixedValue(g);
            m_constrained = true;
            // The right-hand side of a fixed row holds its value.
            res.imposed[i] = 1;
          }
          if (!map.isIdentified(g))
            continue;
          values[i] = map.getIdentificationValue(g);
          m_values[i] = values[i];
          res.imposed[i] = 1;
          if (!owned[i])
            continue;
          res.rows.push_back(i);
          res.cols.push_back(i);
          res.values.push_back(Real(1));
          for (const auto& [master, coefficient] : rows[i])
          {
            res.rows.push_back(i);
            res.cols.push_back(master);
            res.values.push_back(-coefficient);
          }
        }
        // The map, and so the expansion, is the same on every process.
        if (map.getIdentifiedRows().empty())
          return res;
        res.expansions.starts = upload(execution, starts, "Rodin expansion starts");
        res.expansions.indices = upload(execution, indices, "Rodin expansion masters");
        res.expansions.coefficients =
          upload(execution, coefficients, "Rodin expansion coefficients");
        res.expansions.values = upload(execution, values, "Rodin expansion values");
        execution.fence("Rodin expansions");
        return res;
      }

      /**
       * @brief Replaces the fixed degrees of freedom, before a distributed
       * resolution.
       * @param fixed Whether each local degree of freedom is fixed.
       * @param values Fixed value of each local degree of freedom.
       */
      void setFixed(std::vector<unsigned char> fixed, std::vector<Real> values)
      {
        assert(fixed.size() == values.size());
        m_constrained = std::find(fixed.begin(), fixed.end(),
                          static_cast<unsigned char>(1)) != fixed.end();
        m_fixed = std::move(fixed);
        m_values = std::move(values);
      }

    private:
      struct Condition
      {
        Variational::DirichletBCBase<Real>* dbc;
        Variational::PeriodicBCBase<Real>* pbc;
        size_t offset;
        size_t masterOffset;
      };

      /**
       * @brief Copies a host vector into a new device array.
       * @param execution Execution space of the copy.
       * @param data Host values.
       * @param label Label of the device array.
       * @return Device array holding a copy of @p data.
       */
      template <class T>
      static Array<T> upload(
        const ExecutionSpace& execution, const std::vector<T>& data, const char* label)
      {
        const Array<T> res(::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing,
                             std::string(label)),
          data.size());
        ::Kokkos::deep_copy(execution, res,
          ::Kokkos::View<const T*, ::Kokkos::HostSpace,
            ::Kokkos::MemoryTraits<::Kokkos::Unmanaged>>(data.data(), data.size()));
        return res;
      }

      size_t m_rows;
      std::vector<unsigned char> m_fixed;
      std::vector<Real> m_values;
      bool m_constrained;
      bool m_identified;
      std::vector<Condition> m_conditions;
  };
}

#endif
