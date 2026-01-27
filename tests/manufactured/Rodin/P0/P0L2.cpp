/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Variational.h"
#include "Rodin/Solver/CG.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

namespace Rodin::Tests::Manufactured::P0
{
  class Manufactured_P0_L2_Test : public ::testing::Test
  {
  protected:
    void SetUp() override
    {
      m_mesh = Mesh().UniformGrid(Polytope::Type::Quadrilateral, { 10, 10 });
      m_mesh.scale(1.0 / 9.0);
      m_mesh.getConnectivity().compute(1, 2);
      m_mesh.getConnectivity().compute(2, 1);
      m_mesh.getConnectivity().compute(1, 0);
    }

    const Mesh<Context::Local>& getMesh() const { return m_mesh; }

    static RealFunction rhs()
    {
      return [](const Geometry::Point& p) -> double
      {
        return 2.0 * (p.x() * (1.0 - p.x()) + p.y() * p.y());
      };
    }

  private:
    Mesh<Context::Local> m_mesh;
  };

  TEST_F(Manufactured_P0_L2_Test, P0_L2Projection)
  {
    const auto& mesh = getMesh();

    P0 p0h(mesh);
    TrialFunction p(p0h);
    TestFunction  q(p0h);

    const auto f = rhs();

    Problem p_l2(p, q);
    p_l2 = Integral(p, q) - Integral(f, q);
    CG(p_l2).solve();

    GridFunction diff(p0h);
    diff = Pow(p.getSolution() - f, 2);
    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }
}
