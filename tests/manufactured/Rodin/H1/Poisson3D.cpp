/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <algorithm>
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Solver/ForwardDecls.h"
#include "Rodin/Variational.h"
#include "Rodin/Variational/H1.h"
#include "Rodin/Solver/CG.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

namespace Rodin::Tests::Manufactured::H1Poisson3D
{
  template <size_t M>
  class Manufactured_Poisson3D_H1_Test
    : public ::testing::TestWithParam<Polytope::Type>
  {
    protected:
      void SetUp() override
      {
        m_mesh = Mesh().UniformGrid(GetParam(), { M, M, M });
        m_mesh.scale(1.0 / (M - 1));

        // Full 3D connectivity chain needed by H1 on all supported cell types.
        m_mesh.getConnectivity().compute(2, 3);
        m_mesh.getConnectivity().compute(3, 2);
        m_mesh.getConnectivity().compute(2, 1);
        m_mesh.getConnectivity().compute(1, 0);
      }

      const auto& getMesh() const
      {
        return m_mesh;
      }

      Polytope::Type getGeometry() const
      {
        return GetParam();
      }

    private:
      Mesh<Context::Local> m_mesh;
  };

  using Manufactured_Poisson3D_H1_Test_8 =
    Manufactured_Poisson3D_H1_Test<8>;

  using Manufactured_Poisson3D_H1_Test_16 =
    Manufactured_Poisson3D_H1_Test<16>;

  TEST_P(Manufactured_Poisson3D_H1_Test_8, Poisson3D_P1ExactResidual_H1)
  {
    constexpr auto order = std::integral_constant<size_t, 1>{};
    const auto& mesh = this->getMesh();
    H1 vh(order, mesh);

    auto solution = F::x + F::y + F::z + 1;
    auto f = Zero();

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, solution);

    CG(poisson).setTolerance(1e-12).solve();

    GridFunction u_exact(vh);
    u_exact = solution;

    auto& A = poisson.getLinearSystem().getOperator();
    auto& b = poisson.getLinearSystem().getVector();
    auto& x = poisson.getLinearSystem().getSolution();

    auto r  = A * x - b;
    auto re = A * u_exact.getData() - b;

    const Real scale = std::max<Real>(b.norm(), 1);
    EXPECT_NEAR(r.norm() / scale,  0, 1e-10)
      << "Geometry = " << this->getGeometry();
    EXPECT_NEAR(re.norm() / scale, 0, 1e-12)
      << "Geometry = " << this->getGeometry();

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    EXPECT_NEAR(Integral(diff).compute(), 0, 1e-12)
      << "Geometry = " << this->getGeometry();
  }

  TEST_P(Manufactured_Poisson3D_H1_Test_8, Manufactured_Poisson3D_H1_2)
  {
    const auto pi = Rodin::Math::Constants::pi();
    constexpr auto order = std::integral_constant<size_t, 2>{};

    const auto& mesh = this->getMesh();
    H1 vh(order, mesh);

    TrialFunction u(vh);
    TestFunction  v(vh);

    auto solution = sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z);
    auto f = 3 * pi * pi * solution;

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, solution);

    CG(poisson).setTolerance(1e-10).solve();

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);

    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT)
      << "Geometry = " << this->getGeometry();
  }

  TEST_P(Manufactured_Poisson3D_H1_Test_8, Manufactured_Poisson3D_H1_4)
  {
    const auto pi = Rodin::Math::Constants::pi();
    constexpr auto order = std::integral_constant<size_t, 4>{};

    const auto& mesh = this->getMesh();
    H1 vh(order, mesh);

    TrialFunction u(vh);
    TestFunction  v(vh);

    auto solution = sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z);
    auto f = 3 * pi * pi * solution;

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, solution);

    CG(poisson).setTolerance(1e-10).solve();

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);

    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT)
      << "Geometry = " << this->getGeometry();
  }

  TEST_P(Manufactured_Poisson3D_H1_Test_8, Manufactured_Poisson3D_H1_5)
  {
    const auto pi = Rodin::Math::Constants::pi();
    constexpr auto order = std::integral_constant<size_t, 5>{};

    const auto& mesh = this->getMesh();
    H1 vh(order, mesh);

    TrialFunction u(vh);
    TestFunction  v(vh);

    auto solution = sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z);
    auto f = 3 * pi * pi * solution;

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, solution);

    CG(poisson).setTolerance(1e-10).solve();

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);

    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT)
      << "Geometry = " << this->getGeometry();
  }

  TEST_P(Manufactured_Poisson3D_H1_Test_16, Manufactured_Poisson3D_H1_2_RefinedMesh)
  {
    const auto pi = Rodin::Math::Constants::pi();
    constexpr auto order = std::integral_constant<size_t, 2>{};

    const auto& mesh = this->getMesh();
    H1 vh(order, mesh);

    TrialFunction u(vh);
    TestFunction  v(vh);

    auto solution = sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z);
    auto f = 3 * pi * pi * solution;

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, solution);

    CG(poisson).setTolerance(1e-10).solve();

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);

    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT)
      << "Geometry = " << this->getGeometry();
  }

  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage3D,
    Manufactured_Poisson3D_H1_Test_8,
    ::testing::Values(
      Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron,
      Polytope::Type::Wedge
    )
  );

  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage3D,
    Manufactured_Poisson3D_H1_Test_16,
    ::testing::Values(
      Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron,
      Polytope::Type::Wedge
    )
  );
}
