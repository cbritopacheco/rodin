/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file RandomFloat.h
 * @brief Uniform random floating-point generators for randomized tests.
 */
#ifndef RODIN_TEST_RANDOM_RANDOMFLOAT_H
#define RODIN_TEST_RANDOM_RANDOMFLOAT_H

#include <random>
#include <limits>
#include <cassert>
#include <type_traits>

#include "Rodin/Types.h"

namespace Rodin::Test::Random
{
  /**
   * @brief Uniform random floating-point generator for tests.
   * @tparam T Floating-point value type.
   */
  template <class T = float>
  class RandomFloat
  {
   static_assert(std::is_floating_point_v<T>, "Template parameter T must be of floating point type");
   public:
    /**
     * @brief Constructs a generator over an interval.
     * @param a Lower bound of the distribution.
     * @param b Upper bound of the distribution.
     * @param seed Seed for the underlying generator.
     */
     constexpr RandomFloat(T a = std::numeric_limits<T>::min(),
       T b = std::numeric_limits<T>::max(), unsigned int seed = std::random_device()())
       : m_distrib(a, b),
         m_seed(seed)
     {
       assert(a <= b);
     }

    /**
     * @brief Sets the generator seed.
     * @param seed Seed for the underlying generator.
     * @return Reference to this generator.
     */
     RandomFloat& setSeed(unsigned int seed)
     {
       m_gen.seed(seed);
       return *this;
     }

    /**
     * @brief Gets the seed recorded by this generator.
     * @return Seed value.
     */
     constexpr unsigned int getSeed() const
     {
       return m_seed;
     }

    /**
     * @brief Generates the next random value.
     * @return Random value sampled from the distribution.
     */
     T operator()()
     {
       return m_distrib(m_gen);
     }

   private:
    std::mt19937 m_gen;
    std::uniform_real_distribution<T> m_distrib;
    unsigned int m_seed;
  };

  /// @brief Uniform random generator for Rodin::Real values.
  using RandomReal = RandomFloat<Real>;
}

#endif
