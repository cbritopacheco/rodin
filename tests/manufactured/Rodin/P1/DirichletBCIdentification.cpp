/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Variational.h"
#include "Rodin/Solver/SparseLU.h"
#include "Rodin/Solver/CG.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

/**
 * @brief Manufactured solutions for the identification Dirichlet BC
 *        `DirichletBC(u, A(v))`, where `A(v)` is a linear shape-function
 *        expression.
 *
 * The tests in this file verify the *algebraic* identification by exploiting
 * specific reductions:
 *
 *  - `DirichletBC(u, -u).on(\Gamma)` reduces to `u = 0 on \Gamma`, since the constraint
 *    @f$ u_s = -u_s @f$ implies @f$ 2 u_s = 0 @f$, hence @f$ u_s = 0 @f$.
 *
 *  - `DirichletBC(u, c\cdot u).on(\Gamma)` with @f$ c \neq 1 @f$ also reduces to
 *    @f$ u = 0 @f$ on @f$ \Gamma @f$, since
 *    @f$ u_s = c\,u_s \;\Leftrightarrow\; (1-c)\,u_s = 0 @f$.
 *
 *  - `DirichletBC(u, f\cdot v).on(\Gamma)` with `v` itself constrained to a known
 *    function `g_v` via a value-prescribing BC produces an effective
 *    Dirichlet of value @f$ f(x_s) g_v(x_s) @f$ on each slave node @f$ x_s @f$
 *    (Lagrange same-FES). Comparing solutions confirms the constraint
 *    coefficients are computed exactly.
 *
 * Each test compares two Poisson solves: a *reference* solve using a plain
 * value-prescribing `DirichletBC(u, g)`, and an *identification* solve using
 * `DirichletBC(u, A(u))` (or `A(v)` for the cross-block case). The two
 * resulting solutions must agree at the FE-discretization level (ideally
 * bitwise on Lagrange-exact data, otherwise within numerical roundoff /
 * iterative-solver tolerance).
 */
namespace Rodin::Tests::Manufactured::DirichletBCIdentification
{
  template <size_t M>
  class Manufactured_DirichletBCIdentification_Test
    : public ::testing::TestWithParam<Polytope::Type>
  {
    protected:
      Mesh<Context::Local> getMesh()
      {
        Mesh<Context::Local> mesh;
        mesh = mesh.UniformGrid(GetParam(), { M, M });
        mesh.scale(1.0 / (M - 1));
        mesh.getConnectivity().compute(1, 2);
        return mesh;
      }
  };

  using Test_16x16 =
    Manufactured_DirichletBCIdentification_Test<16>;
  using Test_32x32 =
    Manufactured_DirichletBCIdentification_Test<32>;

  /**
   * @brief `DirichletBC(u, -u).on(\Gamma)` is algebraically equivalent to
   *        `DirichletBC(u, 0).on(\Gamma)` for the Poisson problem.
   *
   * Solves
   * @f[
   *   -\Delta u = 2\pi^2 \sin(\pi x)\sin(\pi y) \;\;\text{in}\;\; \Omega,
   *   \qquad u = 0 \;\;\text{on}\;\; \partial\Omega
   * @f]
   * twice — once with the plain value-prescribing `DirichletBC(u, Zero())`,
   * once with the identification `DirichletBC(u, -u)` — and compares the
   * two solutions vertex-by-vertex.
   */
  TEST_P(Test_32x32, NegativeSelfIdentificationEquivalentToZero)
  {
    const Real pi = Rodin::Math::Constants::pi();
    const auto f = 2 * pi * pi * sin(pi * F::x) * sin(pi * F::y);

    // ---- Reference: plain homogeneous Dirichlet ----
    Mesh meshRef = this->getMesh();
    P1 vhRef(meshRef);
    TrialFunction uRef(vhRef);
    TestFunction  vRef(vhRef);
    Problem refProblem(uRef, vRef);
    refProblem = Integral(Grad(uRef), Grad(vRef))
               - Integral(f, vRef)
               + DirichletBC(uRef, Zero());
    SparseLU(refProblem).solve();
    const auto& uRefData = uRef.getSolution().getData();

    // ---- Identification path: DirichletBC(u, -u) ----
    Mesh meshId = this->getMesh();
    P1 vhId(meshId);
    TrialFunction uId(vhId);
    TestFunction  vId(vhId);
    Problem idProblem(uId, vId);
    idProblem = Integral(Grad(uId), Grad(vId))
              - Integral(f, vId)
              + DirichletBC(uId, -uId);
    SparseLU(idProblem).solve();
    const auto& uIdData = uId.getSolution().getData();

    // The two solutions should agree to within roundoff / direct-solver
    // tolerance. Using SparseLU (direct) for both rules out iterative
    // residual differences.
    ASSERT_EQ(uRefData.size(), uIdData.size());
    for (Eigen::Index i = 0; i < uRefData.size(); i++)
      EXPECT_NEAR(uRefData(i), uIdData(i), 1e-10);
  }

  /**
   * @brief `DirichletBC(u, c\cdot u).on(\Gamma)` with @f$ c \neq 1 @f$ and @f$ c \neq 0 @f$
   *        reduces to @f$ u=0 @f$ on @f$ \Gamma @f$ for the Poisson problem.
   *
   * For each @f$ c \in \{ -1, -2, 0.5, 3 \} @f$ the identification
   * solution must coincide with the homogeneous-Dirichlet solution.
   *
   * @note `c = 1` is excluded because @f$ u_s = u_s @f$ is vacuous
   *       (singular row). `c = 0` is excluded because the assembler's
   *       strict-zero pruning drops the empty constraint row entirely;
   *       to pin DOFs to zero use `DirichletBC(u, Zero())`, not
   *       `DirichletBC(u, 0\cdot v)`.
   */
  TEST_P(Test_32x32, ScalarSelfIdentificationEquivalentToZero)
  {
    const Real pi = Rodin::Math::Constants::pi();
    const auto f = 2 * pi * pi * sin(pi * F::x) * sin(pi * F::y);

    // Reference: plain homogeneous Dirichlet ---
    Mesh meshRef = this->getMesh();
    P1 vhRef(meshRef);
    TrialFunction uRef(vhRef);
    TestFunction  vRef(vhRef);
    Problem refProblem(uRef, vRef);
    refProblem = Integral(Grad(uRef), Grad(vRef))
               - Integral(f, vRef)
               + DirichletBC(uRef, Zero());
    SparseLU(refProblem).solve();
    const auto& uRefData = uRef.getSolution().getData();

    for (const Real c : { Real(-1), Real(-2), Real(0.5), Real(3) })
    {
      Mesh meshId = this->getMesh();
      P1 vhId(meshId);
      TrialFunction uId(vhId);
      TestFunction  vId(vhId);
      Problem idProblem(uId, vId);
      idProblem = Integral(Grad(uId), Grad(vId))
                - Integral(f, vId)
                + DirichletBC(uId, RealFunction(c) * uId);
      SparseLU(idProblem).solve();
      const auto& uIdData = uId.getSolution().getData();

      ASSERT_EQ(uRefData.size(), uIdData.size());
      for (Eigen::Index i = 0; i < uRefData.size(); i++)
      {
        EXPECT_NEAR(uRefData(i), uIdData(i), 1e-10)
          << "Mismatch at i=" << i << " for c=" << c;
      }
    }
  }

  /**
   * @brief Verifies that the L2 error of the manufactured solution under
   *        `DirichletBC(u, -u)` is the same as under `DirichletBC(u, 0)`.
   *
   * Same as above but expressed as an integral comparison rather than
   * vertex-wise comparison — guards against any silent rescaling of the
   * solution vector.
   */
  TEST_P(Test_16x16, NegativeSelfIdentificationFEError)
  {
    const Real pi = Rodin::Math::Constants::pi();
    const auto f        = 2 * pi * pi * sin(pi * F::x) * sin(pi * F::y);
    const auto solution = sin(pi * F::x) * sin(pi * F::y);

    // Identification path
    Mesh mesh = this->getMesh();
    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction  v(vh);
    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, -u);
    SparseLU(poisson).solve();

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    const Real error = Integral(diff).compute();

    // Same tolerance as Manufactured_Poisson_Test_16x16.Poisson_SimpleSine
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @brief P1-exact solution test: `u = 1 + x + 2y`, harmonic, with
   *        `DirichletBC(u, A(u))` reducing the value of u on the boundary.
   *
   * For an affine solution the FE solution is exact at the nodes. Using
   * `DirichletBC(u, 2*u + Zero())` does not type-check (mixed types), so
   * we test `DirichletBC(u, 2*u)` and verify the boundary collapses to 0
   * (since (1-2) u_s = 0 ⇒ u_s = 0). This means the solution under this
   * BC is the same as Poisson with homogeneous boundary, which can be
   * computed independently and compared.
   */
  TEST_P(Test_16x16, ScalarSelfIdentificationP1Exact)
  {
    Mesh meshRef = this->getMesh();
    P1 vhRef(meshRef);
    TrialFunction uRef(vhRef);
    TestFunction  vRef(vhRef);
    Problem refProblem(uRef, vRef);
    refProblem = Integral(Grad(uRef), Grad(vRef))
               + DirichletBC(uRef, Zero());
    SparseLU(refProblem).solve();
    const auto& uRefData = uRef.getSolution().getData();

    Mesh meshId = this->getMesh();
    P1 vhId(meshId);
    TrialFunction uId(vhId);
    TestFunction  vId(vhId);
    Problem idProblem(uId, vId);
    idProblem = Integral(Grad(uId), Grad(vId))
              + DirichletBC(uId, RealFunction(2.0) * uId);
    SparseLU(idProblem).solve();
    const auto& uIdData = uId.getSolution().getData();

    ASSERT_EQ(uRefData.size(), uIdData.size());
    for (Eigen::Index i = 0; i < uRefData.size(); i++)
      EXPECT_NEAR(uRefData(i), uIdData(i), 1e-10);
  }

  /**
   * @brief The identification BC produces the documented
   *        @c DirichletBCBase::IdentifiedDOFs variant alternative, never
   *        the @c DirichletBCBase::ValueDOFs alternative — used as a
   *        runtime invariant during the manufactured tests.
   */
  TEST_P(Test_16x16, IdentificationAlternativeIsActive)
  {
    Mesh mesh = this->getMesh();
    P1 vh(mesh);
    TrialFunction u(vh);
    TrialFunction v(vh);

    DirichletBC dbc(u, RealFunction(-1.0) * v);
    dbc.assemble();

    using IdentifiedDOFs = DirichletBCBase<Real>::IdentifiedDOFs;
    EXPECT_TRUE(std::holds_alternative<IdentifiedDOFs>(dbc.getDOFs()));
    EXPECT_FALSE(std::holds_alternative<IndexMap<Real>>(dbc.getDOFs()));

    ASSERT_TRUE(dbc.getValueUUID().has_value());
    EXPECT_EQ(*dbc.getValueUUID(), v.getUUID());
  }

  TEST_P(Test_16x16, CrossFieldVectorMasterMatchesValueReference)
  {
    const auto etaX = 1.0 + F::x;
    const auto etaY = 2.0 - F::y;
    const auto uBoundary = 2.0 * etaX - 0.5 * etaY;

    Mesh meshRef = this->getMesh();
    P1 scalarRef(meshRef);
    P1 vectorRef(meshRef, meshRef.getSpaceDimension());
    TrialFunction uRef(scalarRef);
    TestFunction  vRef(scalarRef);
    TrialFunction etaRef(vectorRef);
    TestFunction  zetaRef(vectorRef);

    Problem refProblem(uRef, vRef, etaRef, zetaRef);
    refProblem =
      Integral(uRef, vRef)
      + Integral(etaRef, zetaRef)
      + DirichletBC(uRef, uBoundary)
      + DirichletBC(etaRef, VectorFunction{ etaX, etaY });
    SparseLU(refProblem).solve();

    Mesh meshId = this->getMesh();
    P1 scalarId(meshId);
    P1 vectorId(meshId, meshId.getSpaceDimension());
    TrialFunction uId(scalarId);
    TestFunction  vId(scalarId);
    TrialFunction etaId(vectorId);
    TestFunction  zetaId(vectorId);

    Problem idProblem(uId, vId, etaId, zetaId);
    idProblem =
      Integral(uId, vId)
      + Integral(etaId, zetaId)
      + DirichletBC(
          uId,
          RealFunction(2.0) * etaId.x() + RealFunction(-0.5) * etaId.y())
      + DirichletBC(etaId, VectorFunction{ etaX, etaY });
    SparseLU(idProblem).solve();

    const auto& uRefData = uRef.getSolution().getData();
    const auto& uIdData  = uId.getSolution().getData();
    ASSERT_EQ(uRefData.size(), uIdData.size());
    for (Eigen::Index i = 0; i < uRefData.size(); i++)
      EXPECT_NEAR(uRefData(i), uIdData(i), 1e-10) << "u dof " << i;

    const auto& etaRefData = etaRef.getSolution().getData();
    const auto& etaIdData  = etaId.getSolution().getData();
    ASSERT_EQ(etaRefData.size(), etaIdData.size());
    for (Eigen::Index i = 0; i < etaRefData.size(); i++)
      EXPECT_NEAR(etaRefData(i), etaIdData(i), 1e-10) << "eta dof " << i;
  }

  INSTANTIATE_TEST_SUITE_P(
    MeshParams16x16,
    Test_16x16,
    ::testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral));

  INSTANTIATE_TEST_SUITE_P(
    MeshParams32x32,
    Test_32x32,
    ::testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral));
}
