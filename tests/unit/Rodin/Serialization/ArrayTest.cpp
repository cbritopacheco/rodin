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
#include <Rodin/Serialization/EigenMatrix.h>

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

namespace Rodin::Tests::Unit
{
  /// @brief Rodin vectors preserve Eigen's archive representation and ownership.
  TYPED_TEST(ArraySerialization, EigenBackedVectorRoundTrip)
  {
    using Output = typename TypeParam::first_type;
    using Input = typename TypeParam::second_type;
    Math::Vector<Complex> dynamic(2);
    dynamic << Complex(1, 2), Complex(-3, 4);
    const Math::FixedSizeVector<Complex, 2> fixed = dynamic;
    std::stringstream stream;
    {
      Output archive(stream);
      archive << dynamic << fixed;
    }
    Math::Vector<Complex> loaded(7);
    Math::FixedSizeVector<Complex, 2> loadedFixed;
    {
      Input archive(stream);
      archive >> loaded >> loadedFixed;
    }
    EXPECT_EQ(loaded.size(), 2);
    EXPECT_EQ(loaded, dynamic);
    EXPECT_EQ(loadedFixed, fixed);
    EXPECT_EQ(loaded.dot(loadedFixed), Complex(30));
    // Legacy Eigen vectors use the same coefficient and dimension encoding.
    std::stringstream legacy;
    {
      const Eigen::VectorX<Complex> native = dynamic;
      Output archive(legacy);
      archive << native;
    }
    {
      Input archive(legacy);
      archive >> loaded;
    }
    EXPECT_EQ(loaded, dynamic);
  }
}

namespace Rodin::Tests::Unit
{
  /// @brief Matrix owners retain Eigen archive compatibility for rectangular and empty shapes.
  TYPED_TEST(ArraySerialization, EigenBackedMatrixRoundTrip)
  {
    using Output = typename TypeParam::first_type;
    using Input = typename TypeParam::second_type;
    for (int rows = 0; rows <= 3; ++rows)
    {
      for (int cols = 0; cols <= 3; ++cols)
      {
        Math::Matrix<Complex> original(rows, cols);
        for (int i = 0; i < rows; ++i)
        {
          for (int j = 0; j < cols; ++j)
            original(i, j) = Complex(i + 3 * j + 1, i - j - 2);
        }
        std::stringstream stream;
        {
          Output archive(stream);
          archive << original;
        }
        Math::Matrix<Complex> loaded = Math::Matrix<Complex>::Constant(5, 7, Complex(99));
        {
          Input archive(stream);
          archive >> loaded;
        }
        EXPECT_EQ(loaded.rows(), rows);
        EXPECT_EQ(loaded.cols(), cols);
        EXPECT_EQ(loaded, original);
        stream.clear();
        stream.seekg(0);
        {
          Eigen::MatrixX<Complex> native;
          Input archive(stream);
          archive >> native;
          EXPECT_EQ(native, original);
        }
        std::stringstream legacy;
        {
          const Eigen::MatrixX<Complex> native = original;
          Output archive(legacy);
          archive << native;
        }
        {
          Input archive(legacy);
          archive >> loaded;
        }
        EXPECT_EQ(loaded, original);
      }
    }
    const Math::FixedSizeMatrix<Complex, 2, 3> rectangle =
      Math::FixedSizeMatrix<Complex, 2, 3>::Constant(Complex(2, 3));
    const Math::FixedSizeMatrix<Complex, 1, 3> row =
      Math::FixedSizeMatrix<Complex, 1, 3>::Constant(Complex(1, -2));
    std::stringstream stream;
    {
      Output archive(stream);
      archive << rectangle << row;
    }
    Math::FixedSizeMatrix<Complex, 2, 3> loaded;
    Math::FixedSizeMatrix<Complex, 1, 3> loadedRow;
    {
      Input archive(stream);
      archive >> loaded >> loadedRow;
    }
    EXPECT_EQ(loaded, rectangle);
    EXPECT_EQ(loadedRow, row);
    EXPECT_EQ(loadedRow.dot(row), Complex(15));
  }
}
