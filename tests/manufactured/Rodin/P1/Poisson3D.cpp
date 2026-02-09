/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Context/Local.h"
#include "Rodin/Variational.h"
#include "Rodin/Solver/CG.h"
#include "Rodin/Test/Random.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;
using namespace Rodin::Test::Random;


/**
 * @brief Manufactured solutions for the 3D Poisson problem on Tetrahedron meshes.
 *
 * The system is given by:
 * @f[
 * \left\{
 * \begin{aligned}
 *   -\Delta u &= f \quad \text{in } \Omega,\\
 *  u &= g \quad \text{on } \partial\Omega.
 * \end{aligned}
 * \right.
 * @f]
 *
 * The weak formulation is: Find @f$ u \in V @f$ such that
 * @f[
 *   \int_\Omega \nabla u \cdot \nabla v \,dx = \int_\Omega f \, v \,dx,
 * @f]
 * for all @f$ v \in V @f$, with the essential boundary condition
 * @f[
 *   u = g \quad \text{on } \partial\Omega.
 * @f]
 */
namespace Rodin::Tests::Manufactured::Poisson3D
{
  template <size_t M>
  class Manufactured_Poisson3D_Test : public ::testing::Test
  {
  protected:
    void SetUp() override
    {
      m_mesh = Mesh().UniformGrid(Polytope::Type::Tetrahedron, {M, M, M});
      m_mesh.scale(1.0 / (M - 1));
      m_mesh.getConnectivity().compute(2, 3);
    }

    const auto& getMesh() const { return m_mesh; }

    template <class FES, class Expr, class GF>
    static Real relL2Error(const FES& vh, const GF& uh, const Expr& uex)
    {
      GridFunction err(vh);
      err = Pow(uh - uex, 2);
      const Real l2e = Math::sqrt(Integral(err).compute());

      GridFunction sol(vh);
      sol = Pow(uex, 2);
      const Real l2u = Math::sqrt(Integral(sol).compute());

      return l2e / (l2u + 1e-30);
    }

  private:
    Mesh<Context::Local> m_mesh;
  };

  using Manufactured_Poisson3D_Test_8 =
    Rodin::Tests::Manufactured::Poisson3D::Manufactured_Poisson3D_Test<8>;
  using Manufactured_Poisson3D_Test_16 =
    Rodin::Tests::Manufactured::Poisson3D::Manufactured_Poisson3D_Test<16>;
  using Manufactured_Poisson3D_Test_32 =
    Rodin::Tests::Manufactured::Poisson3D::Manufactured_Poisson3D_Test<32>;

  TEST_F(Manufactured_Poisson3D_Test_16, MeshVolumeIsOne)
  {
    const auto& mesh = this->getMesh();
    P1 vh(mesh);

    GridFunction one(vh);
    one = 1.0;

    const Real vol = Integral(one).compute();
    EXPECT_NEAR(vol, 1.0, 1e-12);
  }

  TEST_F(Manufactured_Poisson3D_Test_8, Poisson3D_P1ExactResidual)
  {
    const auto& mesh = this->getMesh();
    P1 vh(mesh);

    auto solution = F::x + F::y + F::z + 1;
    auto f = Zero();

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, solution);

    CG(poisson).solve();

    GridFunction u_exact(vh);
    u_exact = solution;

    auto& A = poisson.getLinearSystem().getOperator();
    auto& b = poisson.getLinearSystem().getVector();
    auto& x = poisson.getLinearSystem().getSolution();

    auto r = A * x - b;
    auto re = A * u_exact.getData() - b;

    const Real scale = std::max<Real>(b.norm(), 1);
    EXPECT_NEAR(r.norm() / scale, 0, 1e-10);
    EXPECT_NEAR(re.norm() / scale, 0, 1e-12);

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    EXPECT_NEAR(Integral(diff).compute(), 0, 1e-12);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1] \times [0, 1]
   * @f]
   *
   * @f[
   *  \Gamma = \partial \Omega
   * @f]
   *
   * @f[
   *  f(x, y, z) = 3 \pi^2 \sin(\pi x) \sin(\pi y) \sin(\pi z)
   * @f]
   *
   * @f[
   *  g(x, y, z) = 0
   * @f]
   *
   * @f[
   *  u(x, y, z) = \sin(\pi x) \sin(\pi y) \sin(\pi z)
   * @f]
   *
   */
  TEST_F(Manufactured_Poisson3D_Test_16, Poisson3D_SimpleSine)
  {
    auto pi = Rodin::Math::Constants::pi();

    const auto& mesh = this->getMesh();

    P1 vh(mesh);
    auto f = 3 * pi * pi * sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z);

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());
    CG(poisson).solve();

    auto solution = sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z);

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);

    Real error = Integral(diff).compute();
    EXPECT_NEAR(0.5 * error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1] \times [0, 1]
   * @f]
   *
   * @f[
   *  u(x, y, z) = x (1-x) y (1-y) z (1-z)
   * @f]
   *
   * @f[
   *  f(x, y, z) = 2 y (1-y) z (1-z) + 2 x (1-x) z (1-z) + 2 x (1-x) y (1-y)
   * @f]
   *
   * @f[
   *  g(x, y, z) = 0
   * @f]
   *
   */
  TEST_F(Manufactured_Poisson3D_Test_8, Poisson3D_Polynomial)
  {
    const auto& mesh = this->getMesh();

    P1 vh(mesh);
    auto f = 2 * F::y * (1 - F::y) * F::z * (1 - F::z)
           + 2 * F::x * (1 - F::x) * F::z * (1 - F::z)
           + 2 * F::x * (1 - F::x) * F::y * (1 - F::y);

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());
    CG(poisson).solve();

    auto solution = F::x * (1 - F::x) * F::y * (1 - F::y) * F::z * (1 - F::z);

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1] \times [0, 1]
   * @f]
   *
   * @f[
   *  u(x, y, z) = x (1-x) \sin(\pi y) \sin(\pi z)
   * @f]
   *
   * @f[
   *  f(x, y, z) = 2 \sin(\pi y) \sin(\pi z) + 2 \pi^2 x (1-x) \sin(\pi y) \sin(\pi z)
   * @f]
   *
   * @f[
   *  g(x, y, z) = 0
   * @f]
   *
   */
  TEST_F(Manufactured_Poisson3D_Test_16, Poisson3D_TrigonometricPolynomial)
  {
    auto pi = Rodin::Math::Constants::pi();

    const auto& mesh = this->getMesh();

    P1 vh(mesh);
    auto f = 2 * sin(pi * F::y) * sin(pi * F::z)
           + 2 * pi * pi * F::x * (1 - F::x) * sin(pi * F::y) * sin(pi * F::z);

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());
    CG(poisson).solve();

    auto solution = F::x * (1 - F::x) * sin(pi * F::y) * sin(pi * F::z);

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1] \times [0, 1]
   * @f]
   *
   * @f[
   *  u(x, y, z) = \cos(\pi x) \cos(\pi y) \cos(\pi z)
   * @f]
   *
   * @f[
   *  f(x, y, z) = 3 \pi^2 \cos(\pi x) \cos(\pi y) \cos(\pi z)
   * @f]
   *
   * @f[
   *  g(x, y, z) = u(x,y,z)
   * @f]
   *
   */
  TEST_F(Manufactured_Poisson3D_Test_16, Poisson3D_NonhomogeneousDirichlet)
  {
    auto pi = Rodin::Math::Constants::pi();

    const auto& mesh = this->getMesh();

    P1 vh(mesh);
    auto f = 3 * pi * pi * cos(pi * F::x) * cos(pi * F::y) * cos(pi * F::z);

    TrialFunction u(vh);
    TestFunction  v(vh);

    // Use the manufactured solution as Dirichlet data.
    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, cos(pi * F::x) * cos(pi * F::y) * cos(pi * F::z));
    CG(poisson).solve();

    auto solution = cos(pi * F::x) * cos(pi * F::y) * cos(pi * F::z);

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1] \times [0, 1]
   * @f]
   *
   * @f[
   *  u(x, y, z) = \sin(\pi x) \sin(\pi y) e^z
   * @f]
   *
   * @f[
   *  f(x, y, z) = (2\pi^2-1) \sin(\pi x) \sin(\pi y) e^z
   * @f]
   *
   * @f[
   *  g(x, y, z) = u(x,y,z)
   * @f]
   *
   */
  TEST_F(Manufactured_Poisson3D_Test_32, Poisson3D_MixedBoundary)
  {
    auto pi = Rodin::Math::Constants::pi();

    const auto& mesh = this->getMesh();

    P1 vh(mesh);
    auto f = (2 * pi * pi - 1) * sin(pi * F::x) * sin(pi * F::y) * exp(F::z);

    TrialFunction u(vh);
    TestFunction  v(vh);

    // Apply Dirichlet conditions on the entire boundary.
    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, sin(pi * F::x) * sin(pi * F::y) * exp(F::z));
    CG(poisson).solve();

    auto solution = sin(pi * F::x) * sin(pi * F::y) * exp(F::z);

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1] \times [0, 1]
   * @f]
   *
   * @f[
   *  u(x,y,z) = x+y+z
   * @f]
   *
   * @f[
   *  f(x,y,z) = 0
   * @f]
   *
   * @f[
   *  g(x,y,z) = u(x,y,z)
   * @f]
   */
  TEST_F(Manufactured_Poisson3D_Test_8, Poisson3D_LinearNonhomogeneous)
  {
    const auto& mesh = this->getMesh();
    P1 vh(mesh);
    auto f = Zero();
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, F::x + F::y + F::z);
    CG(poisson).solve();
    auto solution = F::x + F::y + F::z;
    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1] \times [0, 1]
   * @f]
   *
   * Vector Poisson problem:
   *
   * @f[
   *  u(x, y, z) = \begin{pmatrix}
   *    \sin(\pi x)\sin(\pi y)\sin(\pi z)\\
   *    \sin(\pi x)\sin(\pi y)\sin(\pi z)\\
   *    \sin(\pi x)\sin(\pi y)\sin(\pi z)
   *  \end{pmatrix}
   * @f]
   *
   * @f[
   *  f(x, y, z) = 3\pi^2 \begin{pmatrix}
   *    \sin(\pi x)\sin(\pi y)\sin(\pi z)\\
   *    \sin(\pi x)\sin(\pi y)\sin(\pi z)\\
   *    \sin(\pi x)\sin(\pi y)\sin(\pi z)
   *  \end{pmatrix}
   * @f]
   *
   * @f[
   *  g(x, y, z) = 0
   * @f]
   *
   */
  TEST_F(Manufactured_Poisson3D_Test_32, VectorPoisson3D_SimpleSine)
  {
    auto pi = Rodin::Math::Constants::pi();

    const auto& mesh = this->getMesh();

    P1 vh(mesh, mesh.getSpaceDimension());
    auto sol_expr = sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z);
    VectorFunction f{
      3 * pi * pi * sol_expr,
      3 * pi * pi * sol_expr,
      3 * pi * pi * sol_expr
    };

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Jacobian(u), Jacobian(v))
            - Integral(f, v)
            + DirichletBC(u, Zero(mesh.getSpaceDimension()));
    CG(poisson).solve();

    VectorFunction solution{
      sol_expr,
      sol_expr,
      sol_expr
    };

    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(Frobenius(u.getSolution() - solution), 2);

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1] \times [0, 1]
   * @f]
   *
   * Vector Poisson problem:
   *
   * @f[
   *  u(x, y, z) = \begin{pmatrix}
   *    x(1-x)y(1-y)z(1-z)\\
   *    x(1-x)y(1-y)z(1-z)\\
   *    x(1-x)y(1-y)z(1-z)
   *  \end{pmatrix}
   * @f]
   *
   * @f[
   *  f(x, y, z) = \begin{pmatrix}
   *    2y(1-y)z(1-z) + 2x(1-x)z(1-z) + 2x(1-x)y(1-y)\\
   *    2y(1-y)z(1-z) + 2x(1-x)z(1-z) + 2x(1-x)y(1-y)\\
   *    2y(1-y)z(1-z) + 2x(1-x)z(1-z) + 2x(1-x)y(1-y)
   *  \end{pmatrix}
   * @f]
   *
   * @f[
   *  g(x, y, z) = 0
   * @f]
   *
   */
  TEST_F(Manufactured_Poisson3D_Test_8, VectorPoisson3D_Polynomial)
  {
    const auto& mesh = this->getMesh();

    P1 vh(mesh, mesh.getSpaceDimension());
    auto f_expr = 2 * F::y * (1 - F::y) * F::z * (1 - F::z)
                + 2 * F::x * (1 - F::x) * F::z * (1 - F::z)
                + 2 * F::x * (1 - F::x) * F::y * (1 - F::y);
    VectorFunction f{
      f_expr,
      f_expr,
      f_expr
    };

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Jacobian(u), Jacobian(v))
            - Integral(f, v)
            + DirichletBC(u, Zero(mesh.getSpaceDimension()));
    CG(poisson).solve();

    auto sol_expr = F::x * (1 - F::x) * F::y * (1 - F::y) * F::z * (1 - F::z);
    VectorFunction solution{
      sol_expr,
      sol_expr,
      sol_expr
    };

    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(Frobenius(u.getSolution() - solution), 2);

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1] \times [0, 1]
   * @f]
   *
   * Vector Poisson problem with different components:
   *
   * @f[
   *  u(x, y, z) = \begin{pmatrix}
   *    x(1-x)\sin(\pi y)\sin(\pi z)\\
   *    y(1-y)\sin(\pi x)\sin(\pi z)\\
   *    z(1-z)\sin(\pi x)\sin(\pi y)
   *  \end{pmatrix}
   * @f]
   *
   * @f[
   *  f(x, y, z) = \begin{pmatrix}
   *    2\sin(\pi y)\sin(\pi z) + 2\pi^2 x(1-x)\sin(\pi y)\sin(\pi z)\\
   *    2\sin(\pi x)\sin(\pi z) + 2\pi^2 y(1-y)\sin(\pi x)\sin(\pi z)\\
   *    2\sin(\pi x)\sin(\pi y) + 2\pi^2 z(1-z)\sin(\pi x)\sin(\pi y)
   *  \end{pmatrix}
   * @f]
   *
   * @f[
   *  g(x, y, z) = 0
   * @f]
   *
   */
  TEST_F(Manufactured_Poisson3D_Test_16, VectorPoisson3D_TrigonometricPolynomial)
  {
    auto pi = Rodin::Math::Constants::pi();

    const auto& mesh = this->getMesh();

    P1 vh(mesh, mesh.getSpaceDimension());
    VectorFunction f{
      2 * sin(pi * F::y) * sin(pi * F::z)
        + 2 * pi * pi * F::x * (1 - F::x) * sin(pi * F::y) * sin(pi * F::z),
      2 * sin(pi * F::x) * sin(pi * F::z)
        + 2 * pi * pi * F::y * (1 - F::y) * sin(pi * F::x) * sin(pi * F::z),
      2 * sin(pi * F::x) * sin(pi * F::y)
        + 2 * pi * pi * F::z * (1 - F::z) * sin(pi * F::x) * sin(pi * F::y)
    };

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Jacobian(u), Jacobian(v))
            - Integral(f, v)
            + DirichletBC(u, Zero(mesh.getSpaceDimension()));
    CG(poisson).solve();

    VectorFunction solution{
      F::x * (1 - F::x) * sin(pi * F::y) * sin(pi * F::z),
      F::y * (1 - F::y) * sin(pi * F::x) * sin(pi * F::z),
      F::z * (1 - F::z) * sin(pi * F::x) * sin(pi * F::y)
    };

    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(Frobenius(u.getSolution() - solution), 2);

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0,1] \times [0,1] \times [0,1]
   * @f]
   *
   * Vector Poisson problem:
   *
   * @f[
   *  u(x, y, z) = \begin{pmatrix}
   *    \cos(\pi x) \cos(\pi y) \cos(\pi z)\\
   *    \sin(\pi x) \sin(\pi y) \sin(\pi z)\\
   *    \sin(\pi x) \cos(\pi y) \sin(\pi z)
   *  \end{pmatrix}
   * @f]
   *
   * @f[
   *  f(x, y, z) = 3\pi^2 \begin{pmatrix}
   *    \cos(\pi x) \cos(\pi y) \cos(\pi z)\\
   *    \sin(\pi x) \sin(\pi y) \sin(\pi z)\\
   *    \sin(\pi x) \cos(\pi y) \sin(\pi z)
   *  \end{pmatrix}
   * @f]
   *
   * @f[
   *  g(x, y, z) = u(x,y,z)
   * @f]
   *
   */
  TEST_F(Manufactured_Poisson3D_Test_32, VectorPoisson3D_NonhomogeneousDirichlet)
  {
    auto pi = Rodin::Math::Constants::pi();

    const auto& mesh = this->getMesh();

    P1 vh(mesh, mesh.getSpaceDimension());
    VectorFunction f{
      3 * pi * pi * cos(pi * F::x) * cos(pi * F::y) * cos(pi * F::z),
      3 * pi * pi * sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z),
      3 * pi * pi * sin(pi * F::x) * cos(pi * F::y) * sin(pi * F::z)
    };

    TrialFunction u(vh);
    TestFunction  v(vh);

    // Use the manufactured solution as Dirichlet data.
    Problem poisson(u, v);
    poisson = Integral(Jacobian(u), Jacobian(v))
            - Integral(f, v)
            + DirichletBC(u, VectorFunction{
                cos(pi * F::x) * cos(pi * F::y) * cos(pi * F::z),
                sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z),
                sin(pi * F::x) * cos(pi * F::y) * sin(pi * F::z)
              });
    CG(poisson).solve();

    VectorFunction solution{
      cos(pi * F::x) * cos(pi * F::y) * cos(pi * F::z),
      sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z),
      sin(pi * F::x) * cos(pi * F::y) * sin(pi * F::z)
    };

    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(Frobenius(u.getSolution() - solution), 2);

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }
}
