/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_VARIATIONAL_DOFMAPPING_H
#define RODIN_VARIATIONAL_DOFMAPPING_H

/**
 * @file DOFMapping.h
 * @brief Union and local numbering of two finite-element supports.
 */
#include <algorithm>
#include <vector>
#include "Rodin/Types.h"
#include "Rodin/Alert/MemberFunctionException.h"

namespace Rodin::Variational
{
  /**
   * @brief Local numbering of the union of two finite-element DOF sets.
   * @tparam Enabled Whether an expression requires a union of distinct supports.
   *
   * The disabled specialization is empty. Ordinary expressions therefore
   * retain their existing storage and basis path. Enabled mappings retain
   * scratch capacity across bindings and sort only at entity binding time.
   */
  template <bool Enabled = true>
  class DOFMapping
  {};

  /// @brief Union numbering and side-local indices for distinct supports.
  template <>
  class DOFMapping<true>
  {
    public:
      /**
       * @brief Captures the first element-local DOF map.
       * @param lhs First support, in its own local basis order.
       *
       * Each input contains distinct indices, as in the native finite element
       * spaces. Equal indices across inputs denote the same global basis.
       * Absent side-local indices represent a zero extension on that side.
       */
      void setLHS(const IndexArray& lhs)
      {
        m_records.clear();
        for (size_t local = 0; local < static_cast<size_t>(lhs.size()); ++local)
          m_records.push_back({lhs[local], local, false});
      }

      /**
       * @brief Completes the binding with the second support and constructs the union.
       * @param rhs Second support in its local basis order.
       *
       * The first support is captured by setLHS before querying the second.
       * This permits finite element spaces to return reusable scratch arrays.
       */
      void setRHS(const IndexArray& rhs)
      {
        for (size_t local = 0; local < static_cast<size_t>(rhs.size()); ++local)
          m_records.push_back({rhs[local], local, true});
        std::sort(
          m_records.begin(), m_records.end(), [](const auto& first, const auto& second) {
            return first.index < second.index;
          });
        size_t count = 0;
        for (size_t position = 0; position < m_records.size(); ++position)
        {
          if (position == 0 || m_records[position - 1].index != m_records[position].index)
            ++count;
        }
        m_indices.resize(count);
        m_lhs.assign(count, {});
        m_rhs.assign(count, {});
        size_t local = 0;
        for (size_t position = 0; position < m_records.size(); ++position)
        {
          const auto& record = m_records[position];
          if (position > 0 && m_records[position - 1].index != record.index)
            ++local;
          m_indices[local] = record.index;
          auto& side = record.rhs ? m_rhs[local] : m_lhs[local];
          if (side)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "An element-local DOF map contains a repeated index." << Alert::Raise;
          }
          side = record.local;
        }
      }

      /// @returns The union in ascending global DOF order.
      const IndexArray& getIndices() const
      {
        return m_indices;
      }

      /**
       * @brief Looks up a union basis in the first support.
       * @param local Position in the union numbering.
       * @returns First-side local index, or no index for a zero extension.
       */
      const Optional<size_t>& getLHS(size_t local) const
      {
        return m_lhs[local];
      }

      /**
       * @brief Looks up a union basis in the second support.
       * @param local Position in the union numbering.
       * @returns Second-side local index, or no index for a zero extension.
       */
      const Optional<size_t>& getRHS(size_t local) const
      {
        return m_rhs[local];
      }

    private:
      /// @brief One input entry before the supports are merged.
      struct Record
      {
        Index index;
        size_t local;
        bool rhs;
      };

      IndexArray m_indices;
      std::vector<Optional<size_t>> m_lhs, m_rhs;
      std::vector<Record> m_records;
  };
}

#endif
