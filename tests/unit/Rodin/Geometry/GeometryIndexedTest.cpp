/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */

/// @file @brief Exception safety of geometry-indexed value construction.
#include <gtest/gtest.h>
#include "Rodin/Geometry/GeometryIndexed.h"

using namespace Rodin;
using namespace Rodin::Geometry;

namespace
{
  /// @brief Counts completed constructions while a selected construction throws.
  struct CountedValue
  {
    static inline size_t live = 0;
    static inline size_t remaining = 100;

    CountedValue()
    {
      if (remaining-- == 0)
        Alert::Exception() << "Injected value construction failure." << Alert::Raise;
      ++live;
    }

    CountedValue(const CountedValue&)
      : CountedValue()
    {}
    CountedValue(CountedValue&&)
      : CountedValue()
    {}
    ~CountedValue()
    {
      --live;
    }
  };
}

TEST(GeometryIndexed, ConstructorUnwinding)
{
  for (size_t fail = 0; fail < Polytope::Types.size(); ++fail)
  {
    CountedValue::remaining = fail;
    EXPECT_THROW((GeometryIndexed<CountedValue>()), Alert::Exception);
    EXPECT_EQ(CountedValue::live, 0u);
  }
  CountedValue::remaining = 100;
  {
    GeometryIndexed<CountedValue> source;
    const size_t count = Polytope::Types.size();
    for (size_t fail = 0; fail < count; ++fail)
    {
      CountedValue::remaining = fail;
      EXPECT_THROW((GeometryIndexed<CountedValue>(source)), Alert::Exception);
      EXPECT_EQ(CountedValue::live, count);
      CountedValue::remaining = fail;
      EXPECT_THROW((GeometryIndexed<CountedValue>(std::move(source))), Alert::Exception);
      EXPECT_EQ(CountedValue::live, count);
    }
  }
  EXPECT_EQ(CountedValue::live, 0u);
}

TEST(GeometryIndexed, InitializerMapping)
{
  using Type = Polytope::Type;
  GeometryIndexed<int> values{{Type::Point, 0}, {Type::Segment, 1}, {Type::Triangle, 2},
    {Type::Quadrilateral, 3}, {Type::Tetrahedron, 4}, {Type::Hexahedron, 5},
    {Type::Wedge, 6}, {Type::Pyramid, 7}};
  EXPECT_EQ(values[Type::Pyramid], 7);
  EXPECT_EQ(values[Type::Triangle], 2);
  EXPECT_THROW((GeometryIndexed<int>{{Type::Point, 1}}), Alert::Exception);
  EXPECT_THROW((GeometryIndexed<int>{{Type::Point, 0}, {Type::Segment, 1},
                 {Type::Triangle, 2}, {Type::Quadrilateral, 3}, {Type::Tetrahedron, 4},
                 {Type::Hexahedron, 5}, {Type::Wedge, 6}, {Type::Wedge, 7}}),
    Alert::Exception);
}
