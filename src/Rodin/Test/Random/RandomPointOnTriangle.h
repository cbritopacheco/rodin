/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TEST_RANDOM_RANDOMPOINTONTRIANGLE_H
#define RODIN_TEST_RANDOM_RANDOMPOINTONTRIANGLE_H

#include "Rodin/Math/Vector.h"
#include "RandomFloat.h"

namespace Rodin::Test::Random
{
  class PointOnTriangle
  {
    public:
      PointOnTriangle()
        : m_a({{0, 0}}),
          m_b({{1, 0}}),
          m_c({{0, 1}}),
          m_gen(0, 1)
      {}

      Math::Vector2<Scalar> operator()()
      {
        auto x = m_gen();
        auto y = m_gen();
        auto q = std::abs(x - y);
        auto t = 0.5 * (x + y - q);
        auto u = 1 - 0.5 * (q + x + y);
        return
          Math::Vector2<Scalar>{
            { q * m_a.x() + t * m_b.x() + u * m_c.x(),
              q * m_b.y() + t * m_b.y() + u * m_c.y() } };
      }

    private:
      Math::Vector2<Scalar> m_a;
      Math::Vector2<Scalar> m_b;
      Math::Vector2<Scalar> m_c;
      RandomScalar m_gen;
  };
}

#endif
