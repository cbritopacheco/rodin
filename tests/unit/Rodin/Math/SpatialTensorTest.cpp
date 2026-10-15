/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Exhaustive spatial tensor extents, arithmetic and contractions. */
#include <cmath>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <gtest/gtest.h>
#include "Rodin/Alert/Exception.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Math/SpatialTensor.h"

using namespace Rodin;
using namespace Rodin::Math;

namespace Rodin::Tests::Unit
{
  /// @brief Invalid axis diagnostics identify the member and preserve tensor state.
  TEST(SpatialTensor, InvalidAxisReportsMemberContext)
  {
    SpatialTensor<Real> tensor(1, 2, 3);
    const auto extents = tensor.getExtents();
    try
    {
      tensor.getDimension(3);
      FAIL() << "An invalid tensor axis must be rejected.";
    }
    catch (const Alert::Exception& error)
    {
      EXPECT_NE(std::string(error.what()).find("getDimension"), std::string::npos);
    }
    EXPECT_EQ(tensor.getExtents(), extents);
  }

  /// @brief Invalid resize diagnostics identify the member and preserve tensor state.
  TEST(SpatialTensor, InvalidResizeReportsMemberContext)
  {
    SpatialTensor<Real> tensor(1, 2, 3);
    const auto extents = tensor.getExtents();
    try
    {
      tensor.resize({1, 4, 3});
      FAIL() << "An extent above spatial capacity must be rejected.";
    }
    catch (const Alert::Exception& error)
    {
      EXPECT_NE(std::string(error.what()).find("resize"), std::string::npos);
    }
    EXPECT_EQ(tensor.getExtents(), extents);
    EXPECT_EQ(tensor.size(), 6u);
  }

  template <class Scalar, size_t Rank>
  void checkTensorExtents()
  {
    using Tensor = SpatialTensor<Scalar, Rank>;
    size_t combinations = 1;
    for (size_t axis = 0; axis < Rank; ++axis)
      combinations *= 4;
    for (size_t shape = 0; shape < combinations; ++shape)
    {
      SCOPED_TRACE(shape);
      typename Tensor::Extents extents;
      size_t code = shape, count = 1;
      for (size_t axis = 0; axis < Rank; ++axis)
      {
        extents[axis] = code % 4;
        code /= 4;
        count *= extents[axis];
      }
      Tensor a(extents), b(extents);
      EXPECT_EQ(a.size(), count);
      EXPECT_EQ(a.getExtents(), extents);
      for (size_t axis = 0; axis < Rank; ++axis)
        EXPECT_EQ(a.getDimension(axis), extents[axis]);
      EXPECT_THROW(a.getDimension(Rank), Alert::Exception);
      Scalar expectedDot = 0;
      Real expectedSquaredNorm = 0;
      for (size_t i = 0; i < count; ++i)
      {
        EXPECT_EQ(a[i], Scalar(0));
        if constexpr (std::is_same_v<Scalar, Complex>)
        {
          a[i] = Scalar(Real(i + 1), Real(i % 3 + 1));
          b[i] = Scalar(Real(i % 5 + 1), -Real(i % 7 + 1));
        }
        else
        {
          a[i] = Scalar(i + 1);
          b[i] = Scalar(i % 5 + 1);
        }
        expectedDot += a[i] * Math::conj(b[i]);
        expectedSquaredNorm += std::real(a[i] * Math::conj(a[i]));
        typename Tensor::Extents indices;
        size_t offset = i;
        for (size_t axis = Rank; axis-- > 0;)
        {
          indices[axis] = offset % extents[axis];
          offset /= extents[axis];
        }
        EXPECT_EQ(std::apply([&](auto... index) { return a(index...); }, indices), a[i]);
        const auto& readOnly = a;
        EXPECT_EQ(
          std::apply([&](auto... index) { return readOnly(index...); }, indices), a[i]);
      }
      EXPECT_EQ(a.dot(b), expectedDot);
      EXPECT_EQ(a.squaredNorm(), expectedSquaredNorm);
      EXPECT_EQ(a.norm(), std::sqrt(expectedSquaredNorm));
      const auto sum = a + b, difference = a - b, negative = -a;
      const auto scaled = 2 * a, divided = a / 2;
      const auto conjugated = a.conjugate();
      auto added = a, subtracted = a;
      added += b;
      subtracted -= b;
      for (size_t i = 0; i < count; ++i)
      {
        EXPECT_EQ(sum[i], a[i] + b[i]);
        EXPECT_EQ(added[i], sum[i]);
        EXPECT_EQ(difference[i], a[i] - b[i]);
        EXPECT_EQ(subtracted[i], difference[i]);
        EXPECT_EQ(negative[i], -a[i]);
        EXPECT_EQ(scaled[i], Scalar(2) * a[i]);
        EXPECT_EQ(divided[i], a[i] / Scalar(2));
        EXPECT_EQ(conjugated[i], Math::conj(a[i]));
      }
      EXPECT_EQ(sum.getExtents(), extents);
      const auto copy = a;
      Tensor assigned;
      assigned = copy;
      EXPECT_EQ(assigned.getExtents(), extents);
      EXPECT_EQ(assigned.size(), count);
      for (size_t i = 0; i < count; ++i)
        EXPECT_EQ(assigned[i], a[i]);
      assigned.setConstant(Scalar(3));
      for (size_t i = 0; i < count; ++i)
        EXPECT_EQ(assigned[i], Scalar(3));
      assigned.setZero();
      EXPECT_EQ(assigned.squaredNorm(), 0);
      assigned.resize(extents);
      EXPECT_EQ(assigned.size(), count);
      for (size_t axis = 0; axis < Rank; ++axis)
      {
        auto invalid = extents;
        invalid[axis] = Tensor::MaxSize + 1;
        EXPECT_THROW(assigned.resize(invalid), Alert::Exception);
        EXPECT_EQ(assigned.getExtents(), extents);
        EXPECT_EQ(assigned.size(), count);
        EXPECT_THROW(Tensor bad(invalid), Alert::Exception);
      }

      // Independent nested sums check every free and contracted axis, including zero.
      if constexpr (Rank == 3)
      {
        SpatialVector<Complex> v(extents[2]);
        for (size_t k = 0; k < extents[2]; ++k)
          v[k] = Complex(Real(k + 1), -Real(k + 2));
        SpatialVector<Real> realVector(extents[2]);
        for (size_t k = 0; k < extents[2]; ++k)
          realVector[k] = std::real(v[k]);
        const auto realProduct = a * realVector;
        const auto product = a * v;
        EXPECT_EQ(product.rows(), extents[0]);
        EXPECT_EQ(product.cols(), extents[1]);
        for (size_t i = 0; i < extents[0]; ++i)
        {
          for (size_t j = 0; j < extents[1]; ++j)
          {
            Complex expected = 0;
            for (size_t k = 0; k < extents[2]; ++k)
              expected += a[(i * extents[1] + j) * extents[2] + k] * v[k];
            EXPECT_EQ(product(i, j), expected);
            Scalar realExpected = 0;
            for (size_t k = 0; k < extents[2]; ++k)
              realExpected += a[(i * extents[1] + j) * extents[2] + k] * realVector[k];
            EXPECT_EQ(realProduct(i, j), realExpected);
          }
        }
      }
      else if constexpr (Rank == 4)
      {
        SpatialMatrix<Complex> m(extents[2], extents[3]);
        for (size_t k = 0; k < extents[2]; ++k)
        {
          for (size_t l = 0; l < extents[3]; ++l)
            m(k, l) = Complex(Real(k + l + 1), -Real(k + 2 * l + 1));
        }
        SpatialMatrix<Real> realMatrix(extents[2], extents[3]);
        for (size_t k = 0; k < extents[2]; ++k)
        {
          for (size_t l = 0; l < extents[3]; ++l)
            realMatrix(k, l) = std::real(m(k, l));
        }
        const auto realProduct = a * realMatrix;
        const auto product = a * m;
        EXPECT_EQ(product.rows(), extents[0]);
        EXPECT_EQ(product.cols(), extents[1]);
        for (size_t i = 0; i < extents[0]; ++i)
        {
          for (size_t j = 0; j < extents[1]; ++j)
          {
            Complex expected = 0;
            for (size_t k = 0; k < extents[2]; ++k)
            {
              for (size_t l = 0; l < extents[3]; ++l)
              {
                expected +=
                  a[((i * extents[1] + j) * extents[2] + k) * extents[3] + l] * m(k, l);
              }
            }
            EXPECT_EQ(product(i, j), expected);
            Scalar realExpected = 0;
            for (size_t k = 0; k < extents[2]; ++k)
            {
              for (size_t l = 0; l < extents[3]; ++l)
              {
                realExpected +=
                  a[((i * extents[1] + j) * extents[2] + k) * extents[3] + l] *
                  realMatrix(k, l);
              }
            }
            EXPECT_EQ(realProduct(i, j), realExpected);
          }
        }
      }
    }
  }

  TEST(SpatialTensor, RankThreeRealAllAxisExtents)
  {
    checkTensorExtents<Real, 3>();
  }
  TEST(SpatialTensor, RankThreeComplexAllAxisExtents)
  {
    checkTensorExtents<Complex, 3>();
  }
  TEST(SpatialTensor, RankFourRealAllAxisExtents)
  {
    checkTensorExtents<Real, 4>();
  }
  TEST(SpatialTensor, RankFourComplexAllAxisExtents)
  {
    checkTensorExtents<Complex, 4>();
  }
  TEST(SpatialTensor, RankFiveRealAllAxisExtents)
  {
    checkTensorExtents<Real, 5>();
  }
  TEST(SpatialTensor, RankFiveComplexAllAxisExtents)
  {
    checkTensorExtents<Complex, 5>();
  }

  TEST(SpatialTensor, DefaultAndResizeTransitions)
  {
    SpatialTensor<Real> a;
    EXPECT_EQ(a.size(), 0);
    EXPECT_EQ(a.getExtents(), (SpatialTensor<Real>::Extents{0, 0, 0}));
    a.resize({3, 3, 3}).setConstant(2);
    a.resize({1, 1, 1});
    EXPECT_EQ(a.size(), 1);
    EXPECT_EQ(a[0], 2);
    a.resize({3, 0, 2});
    EXPECT_EQ(a.size(), 0);
    EXPECT_EQ(a.norm(), 0);
    a.resize({3, 3, 3}).setZero();
    EXPECT_EQ(a.size(), 27);
    EXPECT_EQ(a.norm(), 0);
  }

#ifndef NDEBUG
  TEST(SpatialTensor, RejectsInvalidAccessAndMismatchedExtents)
  {
    SpatialTensor<Real> a(2, 2, 2), b(1, 2, 2);
    EXPECT_DEATH((void)a[8], "Assertion|assertion");
    EXPECT_DEATH((void)a(0, 2, 0), "Assertion|assertion");
    EXPECT_DEATH((void)a.dot(b), "Assertion|assertion");
    EXPECT_DEATH((void)(a + b), "Assertion|assertion");
    EXPECT_DEATH((void)(a * SpatialVector<Real>(3)), "Assertion|assertion");
    SpatialTensor<Real, 4> c(2, 2, 2, 2);
    EXPECT_DEATH((void)(c * SpatialMatrix<Real>(2, 3)), "Assertion|assertion");
  }
#endif
}
