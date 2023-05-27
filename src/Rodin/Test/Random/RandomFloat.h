/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TEST_RANDOM_RANDOMFLOAT_H
#define RODIN_TEST_RANDOM_RANDOMFLOAT_H

#include <random>
#include <limits>
#include <cassert>
#include <type_traits>

namespace Rodin::Test::Random
{
  template <class T = float>
  class RandomFloat
  {
   static_assert(std::is_floating_point_v<T>, "Template parameter T must be of floating point type");
   public:
    constexpr
    RandomFloat(
       T a = std::numeric_limits<T>::min(),
       T b = std::numeric_limits<T>::max(),
       unsigned int seed = std::random_device()())
      : m_distrib(a, b), m_seed(seed)
    {
      assert(a <= b);
    }

    RandomFloat& setSeed(unsigned int seed)
    {
      m_gen.seed(seed);
    }

    constexpr
    unsigned int getSeed() const
    {
      return m_seed;
    }

    T operator()()
    {
      return m_distrib(m_gen);
    }

   private:
    std::mt19937 m_gen;
    std::uniform_real_distribution<T> m_distrib;
    unsigned int m_seed;
  };
}

#endif

