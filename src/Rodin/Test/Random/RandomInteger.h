/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file RandomInteger.h
 * @brief Uniform random integer generators for randomized tests.
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
    /**
     * @brief Maps integral types to a distribution-compatible type.
     * @tparam T Requested integral type.
     */
    template <class T>
    struct UniformIntDistType
    {
      // libc++ rejects char / signed char / unsigned char as IntType.
      /// @brief Integral type accepted by std::uniform_int_distribution.
        using type = std::conditional_t<std::is_same_v<T, char> ||
            std::is_same_v<T, signed char> || std::is_same_v<T, unsigned char>,
          int, T>;
    };

    /// @brief Helper alias for UniformIntDistType::type.
    template <class T>
    using UniformIntDistTypeT = typename UniformIntDistType<T>::type;
  }

  /**
   * @brief Uniform random integer generator for tests.
   * @tparam T Integral value type.
   */
  template <class T = int>
  class RandomInteger
  {
    static_assert(std::is_integral_v<T>, "Template parameter T must be an integral type");

  public:
    /// @brief Distribution-compatible integer type.
    using dist_type = detail::UniformIntDistTypeT<T>;

    /**
     * @brief Constructs a generator over an interval.
     * @param a Lower bound of the distribution.
     * @param b Upper bound of the distribution.
     * @param seed Seed for the underlying generator.
     */
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

    /**
     * @brief Sets the generator seed.
     * @param seed Seed for the underlying generator.
     * @return Reference to this generator.
     */
    RandomInteger& setSeed(unsigned int seed)
    {
      m_seed = seed;
      m_gen.seed(seed);
      return *this;
    }

    /**
     * @brief Gets the current seed.
     * @return Seed value.
     */
    unsigned int getSeed() const { return m_seed; }

    /**
     * @brief Generates the next random integer.
     * @return Random value sampled from the distribution.
     */
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
