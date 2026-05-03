/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <sstream>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/array.hpp>

#include <Rodin/Geometry/Polytope.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  // ---- Construction ----

  TEST(Geometry_PolytopeKey, DefaultConstruction)
  {
    Polytope::Key key;
    EXPECT_EQ(key.size(), 0);
  }

  TEST(Geometry_PolytopeKey, SizedConstruction)
  {
    Polytope::Key key(3);
    EXPECT_EQ(key.size(), 3);
  }

  TEST(Geometry_PolytopeKey, InitializerListConstruction)
  {
    Polytope::Key key{0, 1, 2};
    EXPECT_EQ(key.size(), 3);
    EXPECT_EQ(key(0), 0);
    EXPECT_EQ(key(1), 1);
    EXPECT_EQ(key(2), 2);
  }

  TEST(Geometry_PolytopeKey, InitializerListConstruction_Segment)
  {
    Polytope::Key key{5, 10};
    EXPECT_EQ(key.size(), 2);
    EXPECT_EQ(key(0), 5);
    EXPECT_EQ(key(1), 10);
  }

  TEST(Geometry_PolytopeKey, InitializerListConstruction_Hex)
  {
    Polytope::Key key{0, 1, 2, 3, 4, 5, 6, 7};
    EXPECT_EQ(key.size(), 8);
    for (std::uint8_t i = 0; i < 8; ++i)
      EXPECT_EQ(key(i), static_cast<Index>(i));
  }

  // ---- Access ----

  TEST(Geometry_PolytopeKey, ParenthesisOperator)
  {
    Polytope::Key key{10, 20, 30};
    EXPECT_EQ(key(0), 10);
    EXPECT_EQ(key(1), 20);
    EXPECT_EQ(key(2), 30);
  }

  TEST(Geometry_PolytopeKey, BracketOperator)
  {
    Polytope::Key key{10, 20, 30};
    EXPECT_EQ(key[0], 10);
    EXPECT_EQ(key[1], 20);
    EXPECT_EQ(key[2], 30);
  }

  TEST(Geometry_PolytopeKey, MutableAccess)
  {
    Polytope::Key key(3);
    key(0) = 100;
    key(1) = 200;
    key(2) = 300;
    EXPECT_EQ(key(0), 100);
    EXPECT_EQ(key(1), 200);
    EXPECT_EQ(key(2), 300);
  }

  TEST(Geometry_PolytopeKey, MutableBracketAccess)
  {
    Polytope::Key key(2);
    key[0] = 42;
    key[1] = 99;
    EXPECT_EQ(key[0], 42);
    EXPECT_EQ(key[1], 99);
  }

  // ---- Resize ----

  TEST(Geometry_PolytopeKey, Resize)
  {
    Polytope::Key key{1, 2, 3};
    EXPECT_EQ(key.size(), 3);
    key.resize(5);
    EXPECT_EQ(key.size(), 5);
    // Original values preserved
    EXPECT_EQ(key(0), 1);
    EXPECT_EQ(key(1), 2);
    EXPECT_EQ(key(2), 3);
  }

  TEST(Geometry_PolytopeKey, ResizeShrink)
  {
    Polytope::Key key{1, 2, 3, 4};
    key.resize(2);
    EXPECT_EQ(key.size(), 2);
    EXPECT_EQ(key(0), 1);
    EXPECT_EQ(key(1), 2);
  }

  // ---- Iteration ----

  TEST(Geometry_PolytopeKey, IteratorRange)
  {
    Polytope::Key key{10, 20, 30};
    std::vector<Index> values;
    for (auto it = key.begin(); it != key.end(); ++it)
      values.push_back(*it);
    ASSERT_EQ(values.size(), 3);
    EXPECT_EQ(values[0], 10);
    EXPECT_EQ(values[1], 20);
    EXPECT_EQ(values[2], 30);
  }

  TEST(Geometry_PolytopeKey, ConstIteratorRange)
  {
    const Polytope::Key key{5, 15, 25};
    std::vector<Index> values;
    for (auto it = key.begin(); it != key.end(); ++it)
      values.push_back(*it);
    ASSERT_EQ(values.size(), 3);
    EXPECT_EQ(values[0], 5);
    EXPECT_EQ(values[1], 15);
    EXPECT_EQ(values[2], 25);
  }

  TEST(Geometry_PolytopeKey, RangeForLoop)
  {
    Polytope::Key key{1, 2, 3, 4};
    Index sum = 0;
    for (auto v : key)
      sum += v;
    EXPECT_EQ(sum, 10);
  }

  TEST(Geometry_PolytopeKey, EmptyIteration)
  {
    Polytope::Key key;
    EXPECT_EQ(key.begin(), key.end());
  }

  // ---- getVertices ----

  TEST(Geometry_PolytopeKey, GetVertices)
  {
    Polytope::Key key{3, 7, 11};
    const auto& verts = key.getVertices();
    EXPECT_EQ(verts[0], 3);
    EXPECT_EQ(verts[1], 7);
    EXPECT_EQ(verts[2], 11);
  }

  // ---- SymmetricHash ----

  TEST(Geometry_PolytopeKey, SymmetricHash_SameForPermutations)
  {
    Polytope::Key::SymmetricHash hasher;
    Polytope::Key k1{0, 1, 2};
    Polytope::Key k2{2, 1, 0};
    Polytope::Key k3{1, 0, 2};
    EXPECT_EQ(hasher(k1), hasher(k2));
    EXPECT_EQ(hasher(k1), hasher(k3));
  }

  TEST(Geometry_PolytopeKey, SymmetricHash_DifferentForDistinctKeys)
  {
    Polytope::Key::SymmetricHash hasher;
    Polytope::Key k1{0, 1, 2};
    Polytope::Key k2{0, 1, 3};
    // Not guaranteed but extremely likely for distinct vertex sets
    EXPECT_NE(hasher(k1), hasher(k2));
  }

  TEST(Geometry_PolytopeKey, SymmetricHash_SegmentPermutations)
  {
    Polytope::Key::SymmetricHash hasher;
    Polytope::Key k1{5, 10};
    Polytope::Key k2{10, 5};
    EXPECT_EQ(hasher(k1), hasher(k2));
  }

  TEST(Geometry_PolytopeKey, SymmetricHash_TetPermutations)
  {
    Polytope::Key::SymmetricHash hasher;
    Polytope::Key k1{0, 1, 2, 3};
    Polytope::Key k2{3, 2, 1, 0};
    Polytope::Key k3{1, 3, 0, 2};
    EXPECT_EQ(hasher(k1), hasher(k2));
    EXPECT_EQ(hasher(k1), hasher(k3));
  }

  // ---- SymmetricEquality ----

  TEST(Geometry_PolytopeKey, SymmetricEquality_SamePermutations)
  {
    Polytope::Key::SymmetricEquality eq;
    Polytope::Key k1{0, 1, 2};
    Polytope::Key k2{2, 1, 0};
    EXPECT_TRUE(eq(k1, k2));
  }

  TEST(Geometry_PolytopeKey, SymmetricEquality_Different)
  {
    Polytope::Key::SymmetricEquality eq;
    Polytope::Key k1{0, 1, 2};
    Polytope::Key k2{0, 1, 3};
    EXPECT_FALSE(eq(k1, k2));
  }

  TEST(Geometry_PolytopeKey, SymmetricEquality_DifferentSizes)
  {
    Polytope::Key::SymmetricEquality eq;
    Polytope::Key k1{0, 1};
    Polytope::Key k2{0, 1, 2};
    EXPECT_FALSE(eq(k1, k2));
  }

  TEST(Geometry_PolytopeKey, SymmetricEquality_Segments)
  {
    Polytope::Key::SymmetricEquality eq;
    Polytope::Key k1{5, 10};
    Polytope::Key k2{10, 5};
    EXPECT_TRUE(eq(k1, k2));
  }

  // ---- Serialization ----

  TEST(Geometry_PolytopeKey, BoostSerialization)
  {
    Polytope::Key original{7, 14, 21};
    std::stringstream ss;
    {
      boost::archive::text_oarchive oa(ss);
      oa << original;
    }
    Polytope::Key loaded;
    {
      boost::archive::text_iarchive ia(ss);
      ia >> loaded;
    }
    EXPECT_EQ(loaded.size(), original.size());
    for (std::uint8_t i = 0; i < original.size(); ++i)
      EXPECT_EQ(loaded(i), original(i));
  }

  TEST(Geometry_PolytopeKey, BoostSerialization_Empty)
  {
    Polytope::Key original;
    std::stringstream ss;
    {
      boost::archive::text_oarchive oa(ss);
      oa << original;
    }
    Polytope::Key loaded(5);  // non-empty initially
    {
      boost::archive::text_iarchive ia(ss);
      ia >> loaded;
    }
    EXPECT_EQ(loaded.size(), 0);
  }

  TEST(Geometry_PolytopeKey, BoostSerialization_MaxSize)
  {
    Polytope::Key original{0, 1, 2, 3, 4, 5, 6, 7};
    std::stringstream ss;
    {
      boost::archive::text_oarchive oa(ss);
      oa << original;
    }
    Polytope::Key loaded;
    {
      boost::archive::text_iarchive ia(ss);
      ia >> loaded;
    }
    EXPECT_EQ(loaded.size(), 8);
    for (std::uint8_t i = 0; i < 8; ++i)
      EXPECT_EQ(loaded(i), static_cast<Index>(i));
  }
}
