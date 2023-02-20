/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_RANGESHAPE_H
#define RODIN_VARIATIONAL_RANGESHAPE_H

#include <vector>
#include <ostream>
#include <cassert>

#include "ForwardDecls.h"

namespace Rodin::Variational
{
  class RangeShape
  {
    public:
      constexpr
      RangeShape(std::initializer_list<int> l)
        : m_height(l.begin()[0]), m_width(l.begin()[1])
      {
        assert(l.size() == 2);
      }

      constexpr
      RangeShape(size_t height, size_t width)
        : m_height(height), m_width(width)
      {}

      constexpr
      RangeShape(int height, int width)
        : m_height(height), m_width(width)
      {
        assert(height > 0);
        assert(width > 0);
      }

      constexpr
      RangeShape(const RangeShape&) = default;

      constexpr
      RangeShape(RangeShape&&) = default;

      inline
      constexpr
      int height() const
      {
        assert(m_height > 0);
        return m_height;
      }

      inline
      constexpr
      int width() const
      {
        assert(m_width > 0);
        return m_width;
      }

      inline
      constexpr
      RangeShape product(const RangeShape& rhs) const
      {
        assert(m_height == rhs.m_width);
        return { m_width, rhs.m_height };
      }

      inline
      constexpr
      RangeShape transpose() const
      {
        return { m_width, m_height };
      }

      inline
      constexpr
      bool operator==(const RangeShape& other)
      {
        return m_height == other.m_height && m_width == other.m_width;
      }

      inline
      constexpr
      bool operator!=(const RangeShape& other)
      {
        return !operator==(other);
      }

    private:
      const int m_height;
      const int m_width;
  };

  std::ostream& operator<<(std::ostream& os, const RangeShape& obj);
}

#endif
