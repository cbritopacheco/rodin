/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_SIMPLEXCOUNT_H
#define RODIN_GEOMETRY_SIMPLEXCOUNT_H

#include <vector>

namespace Rodin::Geometry
{
  class SimplexCount
  {
    public:
      SimplexCount() = default;

      SimplexCount(size_t meshDimension)
        : m_counts(meshDimension + 1, 0)
      {}

      SimplexCount(std::initializer_list<size_t> l)
        : m_counts(l)
      {}

      SimplexCount(const SimplexCount&) = default;

      SimplexCount(SimplexCount&&) = default;

      SimplexCount& operator=(const SimplexCount&) = default;

      SimplexCount& operator=(SimplexCount&&) = default;

      inline
      SimplexCount& initialize(size_t meshDimension)
      {
        m_counts.resize(meshDimension + 1, 0);
        return *this;
      }

      inline
      size_t& at(size_t d)
      {
        return m_counts[d];
      }

      inline
      const size_t& at(size_t d) const
      {
        return m_counts[d];
      }

    private:
      std::vector<size_t> m_counts;
  };
}

#endif
