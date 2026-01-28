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

namespace Rodin::Tests::Manufactured::P0L2
{
  template <size_t NX, size_t NY, size_t NZ = 1>
  class Manufactured_P0_L2_Test : public ::testing::TestWithParam<Polytope::Type>
  {
  protected:
    void SetUp() override
    {
      const auto geom = GetParam();
      if (geom == Polytope::Type::Tetrahedron || geom == Polytope::Type::Hexahedron || geom == Polytope::Type::Wedge)
      {
        m_mesh = Mesh().UniformGrid(geom, { NX, NY, NZ });
        m_mesh.scale(1.0 / (NX - 1));
        m_mesh.getConnectivity().compute(2, 3);
        m_mesh.getConnectivity().compute(3, 2);
      }
      else
      {
        m_mesh = Mesh().UniformGrid(geom, { NX, NY });
        m_mesh.scale(1.0 / (NX - 1));
        m_mesh.getConnectivity().compute(1, 2);
        m_mesh.getConnectivity().compute(2, 1);
      }
      m_mesh.getConnectivity().compute(1, 0);
    }

    const Mesh<Context::Local>& getMesh() const { return m_mesh; }

    static auto rhs_polynomial()
    {
      return RealFunction([](const Geometry::Point& p) -> double
      {
        return 2.0 * (p.x() * (1.0 - p.x()) + p.y() * p.y());
      });
    }

    static auto rhs_sine()
    {
      auto pi = Rodin::Math::Constants::pi();
      return RealFunction([pi](const Geometry::Point& p) -> double
      {
        return std::sin(pi * p.x()) * std::sin(pi * p.y()) * std::sin(pi * p.z());
      });
    }

  private:
    Mesh<Context::Local> m_mesh;
  };

  using Manufactured_P0_L2_Test_10x10 =
    Manufactured_P0_L2_Test<10, 10>;
  using Manufactured_P0_L2_Test_20x20x20 =
    Manufactured_P0_L2_Test<20, 20, 20>;

  TEST_P(Manufactured_P0_L2_Test_10x10, P0_L2Projection_PolynomialRHS)
  {
    const auto& mesh = getMesh();

    P0 p0h(mesh);
    TrialFunction p(p0h);
    TestFunction  q(p0h);

    const auto f = rhs_polynomial();

    Problem p_l2(p, q);
    p_l2 = Integral(p, q) - Integral(f, q);
    CG(p_l2).solve();

    GridFunction diff(p0h);
    diff = Pow(p.getSolution() - f, 2);
    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  TEST_P(Manufactured_P0_L2_Test_20x20x20, P0_L2Projection_SineRHS)
  {
    const auto& mesh = getMesh();

    P0 p0h(mesh);
    TrialFunction p(p0h);
    TestFunction  q(p0h);

    const auto f = rhs_sine();

    Problem p_l2(p, q);
    p_l2 = Integral(p, q) - Integral(f, q);
    CG(p_l2).solve();

    GridFunction diff(p0h);
    diff = Pow(p.getSolution() - f, 2);
    const Real l2_err = Integral(diff).compute();

    GridFunction ref(p0h);
    ref = Pow(f, 2);
    const Real l2_ref = Integral(ref).compute();

    const Real rel_err = l2_ref > 0 ? l2_err / l2_ref : 0;
    EXPECT_LT(rel_err, 1e-3);
  }

  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage2D,
    Manufactured_P0_L2_Test_10x10,
    ::testing::Values(
      Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral)
  );

  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage3D,
    Manufactured_P0_L2_Test_20x20x20,
    ::testing::Values(
      Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron,
      Polytope::Type::Wedge
      )
  );
}
