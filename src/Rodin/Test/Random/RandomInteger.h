/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TEST_RANDOM_RANDOMINTEGER_H
#define RODIN_TEST_RANDOM_RANDOMINTEGER_H

#include <random>
#include <limits>
#include <cassert>
#include <type_traits>

namespace Rodin::Test::Random
{
  namespace detail
  {
    template <class T>
    struct UniformIntDistType
    {
      // libc++ rejects char / signed char / unsigned char as IntType.
      using type = std::conditional_t<
        std::is_same_v<T, char> || std::is_same_v<T, signed char> || std::is_same_v<T, unsigned char>,
        int,
        T>;
    };

    template <class T>
    using UniformIntDistTypeT = typename UniformIntDistType<T>::type;
  }

  template <class T = int>
  class RandomInteger
  {
    static_assert(std::is_integral_v<T>, "Template parameter T must be an integral type");

  public:
    using dist_type = detail::UniformIntDistTypeT<T>;

    RandomInteger(
      T a = std::numeric_limits<T>::min(),
      T b = std::numeric_limits<T>::max(),
      unsigned int seed = std::random_device()())
      : m_gen(seed),
        m_distrib(static_cast<dist_type>(a), static_cast<dist_type>(b)),
        m_seed(seed)
    {
      assert(a <= b);
    }

    RandomInteger& setSeed(unsigned int seed)
    {
      m_seed = seed;
      m_gen.seed(seed);
      return *this;
    }

    unsigned int getSeed() const { return m_seed; }

    T operator()()
    {
      dist_type x = m_distrib(m_gen);
      return static_cast<T>(x);
    }

  private:
    std::mt19937 m_gen;
    std::uniform_int_distribution<dist_type> m_distrib;
    unsigned int m_seed;
  };
}

#endif
