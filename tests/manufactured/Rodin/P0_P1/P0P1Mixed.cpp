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
#include "Rodin/Solver/BiCGSTAB.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

namespace Rodin::Tests::Manufactured::P0P1
{
  class Manufactured_P0P1_Mixed_Test : public ::testing::Test
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

  TEST_F(Manufactured_P0P1_Mixed_Test, P0P1_MixedProblem)
  {
    const auto& mesh = getMesh();

    P0 p0h(mesh);
    P1 p1h(mesh);

    TrialFunction u(p1h);
    TestFunction  v(p1h);

    TrialFunction p(p0h);
    TestFunction  q(p0h);

    const auto f = rhs();

    // First, L2 projection of f onto P0: solve (p, q) = (f, q)
    Problem p_l2(p, q);
    p_l2 = Integral(p, q) - Integral(f, q);
    CG(p_l2).solve();

    // L2 projection of p onto P1: solve (u, v) = (p, v)
    Problem u_l2(u, v);
    u_l2 = Integral(u, v) - Integral(p.getSolution(), v);
    CG(u_l2).solve();

    // Reset for mixed solve
    p.getSolution() = 0.0;
    u.getSolution() = 0.0;

    // Mixed system:
    // (u, v) - (p, v) + (p, q) - (f, q) = 0
    Problem mixed(u, v, p, q);
    mixed = Integral(u, v)
          - Integral(p, v)
          + Integral(p, q)
          - Integral(f, q);

    BiCGSTAB(mixed).solve();

    // Check that u matches the L2 projection of f onto P1 (from u_l2)
    GridFunction diff(p1h);
    diff = Pow(u.getSolution() - u_l2.getSolution(), 2);
    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }
}
