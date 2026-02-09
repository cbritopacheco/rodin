/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Variational
{
  class GetOrderOpsTest : public ::testing::Test
  {
    protected:
      void SetUp() override
      {
        m_mesh = Mesh().UniformGrid(Polytope::Type::Triangle, {2, 2});
        m_mesh.scale(1.0);
        m_mesh.getConnectivity().compute(1, 2);
        m_mesh.getConnectivity().compute(2, 1);
      }

      const Mesh<Context::Local>& mesh() const { return m_mesh; }

    private:
      Mesh<Context::Local> m_mesh;
  };

  TEST_F(GetOrderOpsTest, BasicOrdersAndCompositions)
  {
    auto cell = *mesh().getCell();

    P1 vh(mesh());
    TrialFunction u(vh);
    TestFunction  v(vh);

    GridFunction gf(vh);
    gf = Zero();

    // Base functions
    ASSERT_TRUE(u.getOrder(cell).has_value());
    EXPECT_EQ(*u.getOrder(cell), 1);
    ASSERT_TRUE(v.getOrder(cell).has_value());
    EXPECT_EQ(*v.getOrder(cell), 1);
    ASSERT_TRUE(Zero().getOrder(cell).has_value());
    EXPECT_EQ(*Zero().getOrder(cell), 0);

    // Sum propagates max order
    auto sumOp = u + v;
    ASSERT_TRUE(sumOp.getOrder(cell).has_value());
    EXPECT_EQ(*sumOp.getOrder(cell), 1);

    // Multiplication adds polynomial orders
    auto multOp = u * u;
    ASSERT_TRUE(multOp.getOrder(cell).has_value());
    EXPECT_EQ(*multOp.getOrder(cell), 2);

    // Division only polynomial if denominator is constant
    auto div = u / RealFunction(2.0);
    ASSERT_TRUE(div.getOrder(cell).has_value());
    EXPECT_EQ(*div.getOrder(cell), 1);
    auto divNonConst = u / v;
    EXPECT_FALSE(divNonConst.getOrder(cell).has_value());

    // Pow with non-negative integer exponent
    auto sq = Pow(u, 2);
    ASSERT_TRUE(sq.getOrder(cell).has_value());
    EXPECT_EQ(*sq.getOrder(cell), 2);
    auto powNonInt = Pow(u, 2.5);
    EXPECT_FALSE(powNonInt.getOrder(cell).has_value());

    // abs, sqrt, and exp return order 0 for constant operands
    auto absVal = abs(RealFunction(3.0));
    ASSERT_TRUE(absVal.getOrder(cell).has_value());
    EXPECT_EQ(*absVal.getOrder(cell), 0);
    auto sqrtConst = sqrt(RealFunction(4.0));
    ASSERT_TRUE(sqrtConst.getOrder(cell).has_value());
    EXPECT_EQ(*sqrtConst.getOrder(cell), 0);
    auto expConst = exp(RealFunction(1.0));
    ASSERT_TRUE(expConst.getOrder(cell).has_value());
    EXPECT_EQ(*expConst.getOrder(cell), 0);
    auto absPoly = abs(u);
    EXPECT_FALSE(absPoly.getOrder(cell).has_value());
    auto sqrtPoly = sqrt(u);
    EXPECT_FALSE(sqrtPoly.getOrder(cell).has_value());
    auto expPoly = exp(u);
    EXPECT_FALSE(expPoly.getOrder(cell).has_value());
  }
}
