#ifndef RODIN_MATH_UNIT_H
#define RODIN_MATH_UNIT_H

#include "Rodin/Types.h"

namespace Rodin::Math
{
  /**
   * @brief Base class for units.
   */
  template <class Derived, class T>
  class Unit
  {
    public:
      using Type = T;

      static Unit One()
      {
        return Unit(T(1));
      }

      static Unit Zero()
      {
        return Unit(T(0));
      }

      constexpr
      Unit(T v)
        : m_v(v)
      {}

      constexpr
      Unit(const Unit&) = default;

      constexpr
      Unit(Unit&&) = default;

      constexpr
      Unit& operator=(const Unit&) = default;

      constexpr
      Unit& operator=(Unit&&) = default;

      explicit
      operator T() const
      {
        return m_v;
      }

      inline
      constexpr
      bool operator==(const Unit& other) const
      {
        return m_v == other.m_v;
      }

      inline
      constexpr
      bool operator!=(const Unit& other) const
      {
        return !operator==(other);
      }

      inline
      constexpr
      bool operator<(const Unit& other) const
      {
        return m_v < other.m_v;
      }

      inline
      constexpr
      bool operator>(const Unit& other) const
      {
        return m_v > other.m_v;
      }

      inline
      constexpr
      bool operator<=(const Unit& other) const
      {
        return !operator>(other);
      }

      inline
      constexpr
      bool operator>=(const Unit& other) const
      {
        return !operator<(other);
      }

      inline
      constexpr
      auto operator+(const Unit& other) const
      {
        return Unit(m_v + other.m_v);
      }

      inline
      constexpr
      auto operator-(const Unit& other) const
      {
        return Unit(m_v - other.m_v);
      }

      inline
      constexpr
      auto operator*(const Unit& other) const
      {
        return Unit(m_v * other.m_v);
      }


      inline
      constexpr
      auto operator/(const Unit& other) const
      {
        return Unit(m_v / other.m_v);
      }

      inline
      constexpr
      Unit operator+() const
      {
        return *this;
      }

      inline
      constexpr
      Unit operator-() const
      {
        return Unit(-m_v);
      }

      inline
      constexpr
      Unit& operator+=(const Unit& other)
      {
        m_v += other.m_v;
        return *this;
      }

      inline
      constexpr
      Unit& operator-=(const Unit& other)
      {
        m_v -= other.m_v;
        return *this;
      }

      inline
      constexpr
      Unit& operator*=(const Unit& other)
      {
        m_v *= other.m_v;
        return *this;
      }

      inline
      constexpr
      Unit& operator/=(const Unit& other)
      {
        m_v /= other.m_v;
        return *this;
      }

    private:
      T m_v;
  };
}

#endif

