/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ASSEMBLY_CONSTRAINTMAP_H
#define RODIN_ASSEMBLY_CONSTRAINTMAP_H

#include <cassert>
#include <stdexcept>
#include <vector>

#include "Rodin/Types.h"

namespace Rodin::Assembly
{
  /**
   * @brief Read-only expansion map for value and identification constraints.
   *
   * Identification constraints are stored as rows of the expansion map
   * @f$ x_s = \sum_k c_k x_{m_k} @f$.  For unconstrained DOFs,
   * @ref expand returns the identity row `{ { i, 1 } }`.
   *
   * Once populated, this object is read-only and can be shared by OpenMP
   * assembly loops.
   */
  template <class Scalar>
  class ConstraintMap
  {
    public:
      struct Entry
      {
        Index index;
        Scalar coefficient;
      };

      using Expansion = std::vector<Entry>;

      ConstraintMap() = default;

      explicit ConstraintMap(size_t size)
      {
        reset(size);
      }

      void reset(size_t size)
      {
        m_fixed.assign(size, Optional<Scalar>{});
        m_expansions.clear();
        m_expansions.resize(size);
        m_isIdentified.assign(size, false);
        m_identifiedRows.clear();
        for (Index i = 0; i < size; i++)
          m_expansions[static_cast<size_t>(i)].push_back({ i, Scalar(1) });
      }

      size_t size() const
      {
        return m_expansions.size();
      }

      bool isFixed(Index i) const
      {
        return static_cast<size_t>(i) < m_fixed.size()
            && m_fixed[static_cast<size_t>(i)].has_value();
      }

      bool isIdentified(Index i) const
      {
        return static_cast<size_t>(i) < m_expansions.size()
            && m_isIdentified[static_cast<size_t>(i)];
      }

      Scalar getFixedValue(Index i) const
      {
        assert(isFixed(i));
        return *m_fixed[static_cast<size_t>(i)];
      }

      const Expansion& expand(Index i) const
      {
        assert(static_cast<size_t>(i) < m_expansions.size());
        return m_expansions[static_cast<size_t>(i)];
      }

      const std::vector<Index>& getIdentifiedRows() const
      {
        return m_identifiedRows;
      }

      void addFixed(Index i, Scalar value)
      {
        checkIndex(i);
        if (isIdentified(i))
          throw std::invalid_argument(
            "DirichletBC conflict: a DOF cannot be both value-prescribed and identified.");
        m_fixed[static_cast<size_t>(i)] = value;
      }

      template <class Entries>
      void addIdentification(Index slave, const Entries& entries)
      {
        checkIndex(slave);
        if (isFixed(slave))
          throw std::invalid_argument(
            "DirichletBC conflict: a DOF cannot be both value-prescribed and identified.");
        if (isIdentified(slave))
          throw std::invalid_argument(
            "DirichletBC conflict: a DOF has multiple identification constraints.");

        Expansion expansion;
        for (const auto& e : entries)
        {
          checkIndex(e.index);
          if (e.coefficient != Scalar(0))
            expansion.push_back(e);
        }

        if (expansion.empty())
          return;

        if (expansion.size() == 1
            && expansion.front().index == slave
            && expansion.front().coefficient == Scalar(1))
          return;

        m_expansions[static_cast<size_t>(slave)] = std::move(expansion);
        m_isIdentified[static_cast<size_t>(slave)] = true;
        m_identifiedRows.push_back(slave);
      }

    private:
      void checkIndex(Index i) const
      {
        if (static_cast<size_t>(i) >= m_expansions.size())
          throw std::out_of_range("ConstraintMap index out of range.");
      }

      std::vector<Optional<Scalar>> m_fixed;
      std::vector<Expansion> m_expansions;
      std::vector<bool> m_isIdentified;
      std::vector<Index> m_identifiedRows;
  };
}

#endif
