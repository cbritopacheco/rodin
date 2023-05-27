/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <fstream>
#include <gtest/gtest.h>

#include <Rodin/IO.h>
#include <Rodin/IO/MEDIT.h>

#include <Rodin/Geometry.h>
#include <Rodin/Configure.h>

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;

TEST(Rodin_Geometry_Mesh, 2D_PolytopeTransformation_1)
{
  constexpr const size_t rdim = 2;
  constexpr const size_t sdim = 2;
  constexpr const size_t meshDim = 2;

  Mesh mesh =
    Mesh<Rodin::Context::Serial>::Builder()
    .initialize(sdim)
    .nodes(4)
    .vertex({0, 0})
    .vertex({1, 0})
    .vertex({0, 1})
    .vertex({1, 1})
    .polytope(Polytope::Geometry::Triangle, {0, 1, 2})
    .polytope(Polytope::Geometry::Triangle, {1, 3, 2})
    .finalize();

  {
    const auto& trans = mesh.getPolytopeTransformation(meshDim, 0);
    Math::Vector r(rdim);
    r << 0.5, 0.5;
    Math::Vector p(sdim);
    p << 0.5, 0.5;
    trans.transform(r);
    const auto res = trans.transform(r);
    EXPECT_NEAR((res - p).norm(), 0.0, RODIN_FUZZY_CONSTANT);
  }

  {
    const auto& trans = mesh.getPolytopeTransformation(meshDim, 0);

    {
      Math::Vector r(rdim);
      r << 0.0, 0.0;
      Math::Vector p(sdim);
      p << 0.0, 0.0;
      const auto res = trans.transform(r);
      EXPECT_NEAR((res - p).norm(), 0.0, RODIN_FUZZY_CONSTANT);
    }

    {
      Math::Vector r(rdim);
      r << 1.0, 0.0;
      Math::Vector p(sdim);
      p << 1.0, 0.0;
      const auto res = trans.transform(r);
      EXPECT_NEAR((res - p).norm(), 0.0, RODIN_FUZZY_CONSTANT);
    }

    {
      Math::Vector r(rdim);
      r << 0.0, 1.0;
      Math::Vector p(sdim);
      p << 0.0, 1.0;
      trans.transform(r);
      const auto res = trans.transform(r);
      EXPECT_NEAR((res - p).norm(), 0.0, RODIN_FUZZY_CONSTANT);
    }

    {
      Math::Vector r(rdim);
      r << (1.0 / 3.0), (1.0 / 3.0);
      Math::Vector p(sdim);
      p << (1.0 / 3.0), (1.0 / 3.0);
      const auto res = trans.transform(r);
      EXPECT_NEAR((res - p).norm(), 0.0, RODIN_FUZZY_CONSTANT);
    }
  }

  {
    const auto& trans = mesh.getPolytopeTransformation(meshDim, 1);

    {
      Math::Vector r(rdim);
      r << 0, 0;
      Math::Vector p(sdim);
      p << 1, 0;
      const auto res = trans.transform(r);
      EXPECT_NEAR((res - p).norm(), 0.0, RODIN_FUZZY_CONSTANT);
    }

    {
      Math::Vector r(rdim);
      r << 1, 0;
      Math::Vector p(sdim);
      p << 1, 1;
      const auto res = trans.transform(r);
      EXPECT_NEAR((res - p).norm(), 0.0, RODIN_FUZZY_CONSTANT);
    }

    {
      Math::Vector r(rdim);
      r << 0, 1;
      Math::Vector p(sdim);
      p << 0, 1;
      const auto res = trans.transform(r);
      EXPECT_NEAR((res - p).norm(), 0.0, RODIN_FUZZY_CONSTANT);
    }

    {
      Math::Vector r(rdim);
      r << (1.0 / 3.0), (1.0 / 3.0);
      Math::Vector p(sdim);
      p << (2.0 / 3.0), (2.0 / 3.0);
      const auto res = trans.transform(r);
      EXPECT_NEAR((res - p).norm(), 0.0, RODIN_FUZZY_CONSTANT);
    }
  }
}

