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
    ASSERT_TRUE(gf.getOrder(cell).has_value());
    EXPECT_EQ(*gf.getOrder(cell), 1);
    ASSERT_TRUE(Zero().getOrder(cell).has_value());
    EXPECT_EQ(*Zero().getOrder(cell), 0);

    // Sum propagates max order
    auto sum = u + v;
    ASSERT_TRUE(sum.getOrder(cell).has_value());
    EXPECT_EQ(*sum.getOrder(cell), 1);

    // Mult adds polynomial orders
    auto mult = u * u;
    ASSERT_TRUE(mult.getOrder(cell).has_value());
    EXPECT_EQ(*mult.getOrder(cell), 2);
  }
}
