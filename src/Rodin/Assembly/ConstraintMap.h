/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ASSEMBLY_CONSTRAINTMAP_H
#define RODIN_ASSEMBLY_CONSTRAINTMAP_H

#include <algorithm>
#include <cassert>
#include <vector>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Alert/Raise.h"
#include "Rodin/Types.h"

namespace Rodin::Assembly
{
  /**
   * @brief Read-only expansion map for value and identification constraints.
   *
   * Identification constraints are stored as rows of the expansion map
   * @f$ x_s = \sum_k c_k x_{m_k} + d_s @f$. For unconstrained DOFs,
   * @ref expand returns the identity row `{ { i, 1 } }`.
   *
   * The map supports canonicalization of identification constraints:
   * - duplicate masters are merged;
   * - zero coefficients are pruned;
   * - self-references are removed algebraically;
   * - identity constraints are ignored;
   * - constraints such as @f$ x_s = c x_s + d_s @f$, @f$ c \ne 1 @f$,
   *   are converted to fixed constraints;
   * - transitive identifications are flattened by @ref finalize.
   *
   * After calling @ref finalize, this object can be shared by OpenMP assembly
   * loops as a read-only expansion map.
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
        m_identificationValues.assign(size, Scalar(0));
        m_identifiedRows.clear();

        for (Index i = 0; i < static_cast<Index>(size); i++)
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

      Scalar getIdentificationValue(Index i) const
      {
        check(i);
        return m_identificationValues[static_cast<size_t>(i)];
      }

      void setFixed(Index i, Scalar value)
      {
        check(i);

        if (isIdentified(i))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Conflict detected. A DOF cannot be both value-prescribed and identified."
            << Alert::Raise;
        }

        m_fixed[static_cast<size_t>(i)] = value;
      }

      template <class Entries>
      void setIdentification(Index slave, const Entries& entries)
      {
        setIdentification(slave, entries, Scalar(0));
      }

      template <class Entries>
      void setIdentification(Index slave, const Entries& entries, Scalar value)
      {
        check(slave);

        if (isFixed(slave))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Conflict detected. A DOF cannot be both value-prescribed and identified."
            << Alert::Raise;
        }

        if (isIdentified(slave))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Conflict detected. A DOF cannot have multiple identification constraints."
            << Alert::Raise;
        }

        Expansion canonical = canonicalize(slave, entries, value);

        if (canonical.empty())
          return;

        if (isFixed(slave))
          return;

        m_expansions[static_cast<size_t>(slave)] = std::move(canonical);
        m_identificationValues[static_cast<size_t>(slave)] = value;
        m_isIdentified[static_cast<size_t>(slave)] = true;
        m_identifiedRows.push_back(slave);
      }

      /**
       * @brief Flattens transitive identification constraints.
       *
       * Example:
       * @f[
       *   x_1 = x_2, \qquad x_2 = x_3
       * @f]
       * becomes:
       * @f[
       *   x_1 = x_3, \qquad x_2 = x_3.
       * @f]
       *
       * Cycles such as @f$ x_1 = x_2 @f$, @f$ x_2 = x_1 @f$ are rejected.
       */
      void finalize()
      {
        enum class State
        {
          Unvisited,
          Visiting,
          Done
        };

        std::vector<State> state(size(), State::Unvisited);

        auto flatten = [&] (auto&& self, Index i) -> Expansion
        {
          check(i);

          const size_t idx = static_cast<size_t>(i);

          if (!m_isIdentified[idx])
            return Expansion{ Entry{ i, Scalar(1) } };

          if (state[idx] == State::Visiting)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Invalid identification constraints: cyclic dependency involving DOF "
              << i << "." << Alert::Raise;
          }

          if (state[idx] == State::Done)
            return m_expansions[idx];

          state[idx] = State::Visiting;

          Expansion flattened;
          for (const auto& e : m_expansions[idx])
          {
            check(e.index);

            if (e.coefficient == Scalar(0))
              continue;

            const Expansion masterExpansion = self(self, e.index);
            for (const auto& me : masterExpansion)
            {
              if (me.coefficient == Scalar(0))
                continue;

              accumulate(
                  flattened,
                  me.index,
                  e.coefficient * me.coefficient);
            }
          }

          canonicalize(i, flattened);

          state[idx] = State::Done;
          return m_expansions[idx];
        };

        const std::vector<Index> rows = m_identifiedRows;
        for (const Index row : rows)
        {
          if (isIdentified(row))
            flatten(flatten, row);
        }

        rebuild();
      }

    private:
      void check(Index i) const
      {
        if (static_cast<size_t>(i) >= m_expansions.size())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "ConstraintMap index out of range: " << i
            << " is not in [0, " << m_expansions.size() << ")."
            << Alert::Raise;
        }
      }

      static void accumulate(Expansion& expansion, Index index, Scalar coefficient)
      {
        if (coefficient == Scalar(0))
          return;

        auto it = std::find_if(
            expansion.begin(),
            expansion.end(),
            [&] (const Entry& entry)
            {
              return entry.index == index;
            });

        if (it == expansion.end())
        {
          expansion.push_back({ index, coefficient });
        }
        else
        {
          it->coefficient += coefficient;
        }
      }

      static void prune(Expansion& expansion)
      {
        expansion.erase(
            std::remove_if(
              expansion.begin(),
              expansion.end(),
              [] (const Entry& e)
              {
                return e.coefficient == Scalar(0);
              }),
            expansion.end());
      }

      template <class Entries>
      Expansion merge(const Entries& entries)
      {
        Expansion expansion;

        for (const auto& e : entries)
        {
          check(e.index);

          if (e.coefficient == Scalar(0))
            continue;

          accumulate(expansion, e.index, e.coefficient);
        }

        prune(expansion);
        return expansion;
      }

      /**
       * @brief Canonicalizes the row
       *        @f$ x_s = \sum_k c_k x_{m_k} + d_s @f$.
       *
       * If the row contains a self term,
       * @f[
       *   x_s = c_s x_s + \sum_m c_m x_m + d_s,
       * @f]
       * it is rewritten as
       * @f[
       *   x_s = \sum_m \frac{c_m}{1 - c_s} x_m
       *       + \frac{d_s}{1 - c_s}.
       * @f]
       *
       * Special cases:
       * - @f$ x_s = x_s @f$ is ignored when @f$ d_s=0 @f$;
       * - @f$ x_s = x_s + d_s @f$ with @f$ d_s\ne 0 @f$ is rejected;
       * - @f$ x_s = c x_s + d_s @f$, @f$ c \ne 1 @f$, becomes
       *   @f$ x_s = d_s/(1-c) @f$;
       * - @f$ x_s = x_s + \sum_m c_m x_m + d_s @f$ is rejected.
       */
      template <class Entries>
      Expansion canonicalize(Index slave, const Entries& entries, Scalar& value)
      {
        Expansion expansion = merge(entries);

        Expansion reduced;
        reduced.reserve(expansion.size());

        Scalar selfCoefficient = Scalar(0);

        for (const auto& e : expansion)
        {
          if (e.coefficient == Scalar(0))
            continue;

          if (e.index == slave)
            selfCoefficient += e.coefficient;
          else
            reduced.push_back(e);
        }

        prune(reduced);

        const Scalar scale = Scalar(1) - selfCoefficient;

        if (scale == Scalar(0))
        {
          if (reduced.empty())
          {
            if (value == Scalar(0))
            {
              // x_s = x_s. Identity constraint: ignore.
              return {};
            }

            Alert::MemberFunctionException(*this, __func__)
              << "Invalid affine identification: slave DOF " << slave
              << " is constrained by x_s = x_s + d with nonzero d."
              << Alert::Raise;
          }

          Alert::MemberFunctionException(*this, __func__)
            << "Invalid identification: slave DOF " << slave
            << " appears as a master with coefficient exactly 1 alongside "
            << "other masters or a nonzero defect." << Alert::Raise;
        }

        if (selfCoefficient != Scalar(0))
        {
          for (auto& e : reduced)
            e.coefficient /= scale;
          value /= scale;
        }

        prune(reduced);

        if (reduced.empty())
        {
          // x_s = c x_s with c != 1, or all non-self terms cancelled:
          // (1 - c) x_s = d, hence x_s = d / (1 - c).
          setFixed(slave, value);
          return {};
        }

        return reduced;
      }

      /**
       * @brief Canonicalizes and stores an already merged expansion row.
       *
       * Used by @ref finalize after recursive flattening.
       */
      void canonicalize(Index slave, Expansion& expansion)
      {
        Scalar value = m_identificationValues[static_cast<size_t>(slave)];
        Expansion canonical = canonicalize(slave, expansion, value);

        if (canonical.empty())
        {
          m_expansions[static_cast<size_t>(slave)].clear();
          m_expansions[static_cast<size_t>(slave)].push_back({ slave, Scalar(1) });
          m_isIdentified[static_cast<size_t>(slave)] = false;
          m_identificationValues[static_cast<size_t>(slave)] = Scalar(0);
          return;
        }

        m_expansions[static_cast<size_t>(slave)] = std::move(canonical);
        m_identificationValues[static_cast<size_t>(slave)] = value;
        m_isIdentified[static_cast<size_t>(slave)] = true;
      }

      void rebuild()
      {
        m_identifiedRows.clear();

        for (Index i = 0; i < static_cast<Index>(m_isIdentified.size()); i++)
        {
          if (m_isIdentified[static_cast<size_t>(i)])
            m_identifiedRows.push_back(i);
        }
      }

      std::vector<Optional<Scalar>> m_fixed;
      std::vector<Expansion> m_expansions;
      std::vector<bool> m_isIdentified;
      std::vector<Scalar> m_identificationValues;
      std::vector<Index> m_identifiedRows;
  };
}

#endif
