/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#include <array>
#include <limits>
#include <sstream>
#include <type_traits>
#include <utility>

#include <gtest/gtest.h>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/complex.hpp>

#include <Rodin/Serialization/Array.h>

namespace Rodin::Tests::Unit
{
  template <class Archives>
  class ArraySerialization : public ::testing::Test
  {};

  using ArrayArchives = ::testing::Types<
    std::pair<boost::archive::text_oarchive, boost::archive::text_iarchive>,
    std::pair<boost::archive::binary_oarchive, boost::archive::binary_iarchive>>;
  TYPED_TEST_SUITE(ArraySerialization, ArrayArchives);

  /// @brief Archives replace stale arrays and preserve ordered scalar values and tails.
  TYPED_TEST(ArraySerialization, ReplacesStaleValues)
  {
    using Output = typename TypeParam::first_type;
    using Input = typename TypeParam::second_type;
    const auto check = []<class Scalar>() {
      const std::array<Real, 7> real{3, -7, 0, 11, -2, 19, -5};
      const std::array<Real, 7> imaginary{-2, 3, 7, 0, -11, 5, 13};
      const std::array<Index, 7> indices{
        0, 1, std::numeric_limits<Index>::max(), 7, 11, 19, 5};
      const auto expected = [&](Index index) {
        if constexpr (std::is_same_v<Scalar, Complex>)
          return Scalar(real[index], imaginary[index]);
        else if constexpr (std::is_same_v<Scalar, Index>)
          return indices[index];
        else
          return Scalar(real[index]);
      };
      constexpr Index Tail = 193;
      for (const Index count : {0, 1, 4, 7})
      {
        SCOPED_TRACE(count);
        Array<Scalar> original(count);
        for (Index index = 0; index < count; ++index)
          original[index] = expected(index);
        std::ostringstream encoded;
        {
          Output archive(encoded);
          archive << original << Tail;
        }
        for (const Index initialSize : {0, 2, 9})
        {
          SCOPED_TRACE(initialSize);
          Array<Scalar> loaded(initialSize);
          loaded.setConstant(Scalar(101));
          std::istringstream stream(encoded.str());
          Index tail = 0;
          {
            Input archive(stream);
            archive >> loaded >> tail;
          }
          ASSERT_EQ(static_cast<size_t>(loaded.size()), count);
          EXPECT_EQ(tail, Tail);
          for (Index index = 0; index < count; ++index)
          {
            EXPECT_EQ(loaded[index], expected(index));
            EXPECT_EQ(original[index], expected(index));
          }
        }
      }
    };
    check.template operator()<Real>();
    check.template operator()<Complex>();
    check.template operator()<Index>();
  }
}
