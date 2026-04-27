/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 *
 * Manufactured tests for the Assembly module.
 *
 * These tests verify that the Default assembler (Sequential or OpenMP
 * depending on the build configuration) produces numerically correct
 * solutions across all supported mesh geometry types:
 *   - 1D: Segment
 *   - 2D: Triangle, Quadrilateral
 *   - 3D: Tetrahedron, Hexahedron
 *
 * Each test uses a manufactured (P1-exact or polynomial) solution so that
 * the error in the assembled discrete system is roundoff, regardless of mesh
 * refinement level.
 */
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Variational.h"
#include "Rodin/Solver/CG.h"
#include "Rodin/Solver/SparseLU.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

namespace Rodin::Tests::Manufactured::Assembly
{
  // =========================================================================
  // 2-D manufactured tests — Triangle and Quadrilateral
  // =========================================================================

  /**
   * @brief 2-D fixture parameterised on polytope type.
   *
   * All tests use an 8×8 uniform grid on [0,1]² so that they are fast.
   */
  template <size_t N>
  class Assembly_2D_Test : public ::testing::TestWithParam<Polytope::Type>
  {
    protected:
      void SetUp() override
      {
        m_mesh = Mesh<Context::Local>::UniformGrid(GetParam(), { N, N });
        m_mesh.scale(1.0 / (N - 1));
        m_mesh.getConnectivity().compute(1, 2);
      }

      const Mesh<Context::Local>& getMesh() const { return m_mesh; }

    private:
      Mesh<Context::Local> m_mesh;
  };

  using Assembly_2D_Test_8x8 = Assembly_2D_Test<8>;

  // -----------------------------------------------------------------------
  // Test 1: Poisson with P1-exact affine manufactured solution.
  //
  // -Δu = 0,  u = x + 2y + 1  on ∂Ω.
  //
  // With f = 0 and an affine Dirichlet BC the discrete solution equals the
  // manufactured solution to machine precision, regardless of mesh size or
  // type.
  // -----------------------------------------------------------------------
  TEST_P(Assembly_2D_Test_8x8, Poisson_P1ExactSolution_ZeroError)
  {
    const auto& mesh = getMesh();
    P1 vh(mesh);

    const auto solution = F::x + 2 * F::y + RealFunction(1.0);
    const auto f = RealFunction(0.0);

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, solution);

    CG(poisson).solve();

    const auto& A = poisson.getLinearSystem().getOperator();
    const auto& b = poisson.getLinearSystem().getVector();
    const auto& x = poisson.getLinearSystem().getSolution();

    // Residual of the computed solution
    const Real res = (A * x - b).norm();

    // Residual of the exact solution coefficients
    GridFunction u_exact(vh);
    u_exact = solution;
    const Real res_exact = (A * u_exact.getData() - b).norm();

    const Real scale = std::max(b.norm(), Real(1));
    EXPECT_NEAR(res / scale,       0.0, 1e-8);
    EXPECT_NEAR(res_exact / scale, 0.0, 1e-12);

    // L2 error
    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(u.getSolution() - solution, 2);
    EXPECT_NEAR(Integral(diff).compute(), 0.0, 1e-12);
  }

  // -----------------------------------------------------------------------
  // Test 2: Reaction-diffusion with quadratic manufactured solution.
  //
  // -Δu + u = f,  u = x*(1-x)*y*(1-y)  on ∂Ω (= 0).
  // f = 2*y*(1-y) + 2*x*(1-x) + x*(1-x)*y*(1-y).
  // -----------------------------------------------------------------------
  TEST_P(Assembly_2D_Test_8x8, ReactionDiffusion_QuadraticSolution_LowError)
  {
    const auto& mesh = getMesh();
    P1 vh(mesh);

    const auto solution = F::x * (RealFunction(1.0) - F::x)
                        * F::y * (RealFunction(1.0) - F::y);

    const auto f = RealFunction(2.0) * F::y * (RealFunction(1.0) - F::y)
                 + RealFunction(2.0) * F::x * (RealFunction(1.0) - F::x)
                 + solution;

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem rxn(u, v);
    rxn = Integral(Grad(u), Grad(v))
        + Integral(u, v)
        - Integral(f, v)
        + DirichletBC(u, RealFunction(0.0));

    CG(rxn).solve();

    // L2 error should be O(h^2) which for 8×8 mesh is small but not roundoff
    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(u.getSolution() - solution, 2);
    const Real error = Integral(diff).compute();
    EXPECT_LT(error, 1e-4);
  }

  // -----------------------------------------------------------------------
  // Test 3: Preassembled BilinearForm produces same system as inline integral.
  //
  // Stiffness assembled separately then passed as a preassembled BF must
  // yield an identical matrix and the same solution as the inline integral.
  // -----------------------------------------------------------------------
  TEST_P(Assembly_2D_Test_8x8, PreassembledBF_SameSystemAsInline)
  {
    const auto& mesh = getMesh();
    P1 vh(mesh);

    const auto solution = F::x + 2 * F::y + RealFunction(1.0);

    // --- inline ---
    TrialFunction u1(vh);
    TestFunction  v1(vh);
    Problem inline_prob(u1, v1);
    inline_prob = Integral(Grad(u1), Grad(v1))
                + DirichletBC(u1, solution);
    inline_prob.assemble();

    // --- preassembled ---
    TrialFunction u2(vh);
    TestFunction  v2(vh);
    BilinearForm stiff(u2, v2);
    stiff = Integral(Grad(u2), Grad(v2));
    stiff.assemble();

    Problem pre_prob(u2, v2);
    pre_prob = stiff + DirichletBC(u2, solution);
    pre_prob.assemble();

    const auto& A1 = inline_prob.getLinearSystem().getOperator();
    const auto& A2 = pre_prob.getLinearSystem().getOperator();
    const auto& b1 = inline_prob.getLinearSystem().getVector();
    const auto& b2 = pre_prob.getLinearSystem().getVector();

    ASSERT_EQ(A1.rows(), A2.rows());
    const Math::SparseMatrix<Real> Adiff = A1 - A2;
    EXPECT_NEAR(Adiff.norm(), 0.0, 1e-12);
    EXPECT_NEAR((b1 - b2).norm(), 0.0, 1e-12);
  }

  // -----------------------------------------------------------------------
  // Test 4: Preassembled LinearForm produces same RHS as inline integral.
  // -----------------------------------------------------------------------
  TEST_P(Assembly_2D_Test_8x8, PreassembledLF_SameRHSAsInline)
  {
    const auto& mesh = getMesh();
    P1 vh(mesh);

    // --- inline ---
    TrialFunction u1(vh);
    TestFunction  v1(vh);
    Problem inline_prob(u1, v1);
    inline_prob = Integral(Grad(u1), Grad(v1))
                - Integral(RealFunction(1.0), v1);
    inline_prob.assemble();

    // --- preassembled ---
    TrialFunction u2(vh);
    TestFunction  v2(vh);
    LinearForm load(v2);
    load = Integral(RealFunction(1.0), v2);
    load.assemble();

    Problem pre_prob(u2, v2);
    pre_prob = Integral(Grad(u2), Grad(v2)) - load;
    pre_prob.assemble();

    const auto& b1 = inline_prob.getLinearSystem().getVector();
    const auto& b2 = pre_prob.getLinearSystem().getVector();
    ASSERT_EQ(b1.size(), b2.size());
    EXPECT_NEAR((b1 - b2).norm(), 0.0, 1e-12);
  }

  INSTANTIATE_TEST_SUITE_P(
    AllGeometries2D,
    Assembly_2D_Test_8x8,
    ::testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral)
  );

  // =========================================================================
  // 3-D manufactured tests — Tetrahedron
  // =========================================================================

  /**
   * @brief 3-D fixture on a small Tetrahedron mesh (4×4×4 grid on [0,1]³).
   */
  class Assembly_Tet_Test : public ::testing::Test
  {
    protected:
      void SetUp() override
      {
        m_mesh = Mesh<Context::Local>::UniformGrid(
          Polytope::Type::Tetrahedron, { 4, 4, 4 });
        m_mesh.scale(1.0 / 3.0);
        m_mesh.getConnectivity().compute(2, 3);
      }

      const Mesh<Context::Local>& getMesh() const { return m_mesh; }

    private:
      Mesh<Context::Local> m_mesh;
  };

  // -----------------------------------------------------------------------
  // Test 5: 3-D Poisson with affine P1-exact manufactured solution.
  //
  // u = x + 2y + 3z + 1,  f = 0 (affine → -Δu = 0).
  // -----------------------------------------------------------------------
  TEST_F(Assembly_Tet_Test, Poisson3D_P1ExactSolution_ZeroError)
  {
    const auto& mesh = getMesh();
    P1 vh(mesh);

    const auto solution = F::x + RealFunction(2.0) * F::y
                        + RealFunction(3.0) * F::z + RealFunction(1.0);
    const auto f = RealFunction(0.0);

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, solution);

    CG(poisson).solve();

    const auto& A = poisson.getLinearSystem().getOperator();
    const auto& b = poisson.getLinearSystem().getVector();
    const auto& x = poisson.getLinearSystem().getSolution();

    GridFunction u_exact(vh);
    u_exact = solution;

    const Real res_exact = (A * u_exact.getData() - b).norm();
    const Real scale     = std::max(b.norm(), Real(1));
    EXPECT_NEAR(res_exact / scale, 0.0, 1e-12);

    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(u.getSolution() - solution, 2);
    EXPECT_NEAR(Integral(diff).compute(), 0.0, 1e-12);
  }

  // -----------------------------------------------------------------------
  // Test 6: 3-D BilinearForm dimensions correct on Tetrahedron mesh.
  // -----------------------------------------------------------------------
  TEST_F(Assembly_Tet_Test, StiffnessMatrix_CorrectDimensions)
  {
    const auto& mesh = getMesh();
    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction  v(vh);

    BilinearForm stiff(u, v);
    stiff = Integral(Grad(u), Grad(v));
    stiff.assemble();

    const auto& A = stiff.getOperator();
    EXPECT_EQ(A.rows(), static_cast<Eigen::Index>(vh.getSize()));
    EXPECT_EQ(A.cols(), static_cast<Eigen::Index>(vh.getSize()));
  }

  // -----------------------------------------------------------------------
  // Test 7: 3-D LinearForm size correct on Tetrahedron mesh.
  // -----------------------------------------------------------------------
  TEST_F(Assembly_Tet_Test, LoadVector_CorrectSize)
  {
    const auto& mesh = getMesh();
    P1 vh(mesh);
    TestFunction v(vh);

    LinearForm load(v);
    load = Integral(RealFunction(1.0), v);
    load.assemble();

    EXPECT_EQ(load.getVector().size(), static_cast<Eigen::Index>(vh.getSize()));
  }

  // =========================================================================
  // P0 manufactured tests — Triangle and Quadrilateral
  // =========================================================================

  /**
   * @brief P0 fixture for 2-D meshes.
   */
  template <size_t N>
  class Assembly_P0_2D_Test : public ::testing::TestWithParam<Polytope::Type>
  {
    protected:
      void SetUp() override
      {
        m_mesh = Mesh<Context::Local>::UniformGrid(GetParam(), { N, N });
        m_mesh.scale(1.0 / (N - 1));
        m_mesh.getConnectivity().compute(1, 2);
      }

      const Mesh<Context::Local>& getMesh() const { return m_mesh; }

    private:
      Mesh<Context::Local> m_mesh;
  };

  using Assembly_P0_2D_Test_8x8 = Assembly_P0_2D_Test<8>;

  /**
   * @brief P0 mass matrix row sums equal element areas (∫_K 1 dK = |K|).
   */
  TEST_P(Assembly_P0_2D_Test_8x8, P0MassMatrix_RowSumsEqualAreas)
  {
    const auto& mesh = getMesh();
    P0 vh(mesh);
    TrialFunction u(vh);
    TestFunction  v(vh);

    BilinearForm mass(u, v);
    mass = Integral(u, v);
    mass.assemble();

    const auto& M = mass.getOperator();

    // For P0, the mass matrix is diagonal: M_ii = |K_i|
    // Sum of all diagonal entries = total mesh area = 1
    Math::Vector<Real> ones(M.cols());
    ones.setOnes();
    const Real total = (M * ones).sum();
    EXPECT_NEAR(total, 1.0, 1e-12);
  }

  /**
   * @brief P0 load vector sums to domain area (∫_Ω 1 dΩ = 1).
   */
  TEST_P(Assembly_P0_2D_Test_8x8, P0LoadVector_SumEqualsArea)
  {
    const auto& mesh = getMesh();
    P0 vh(mesh);
    TestFunction v(vh);

    LinearForm load(v);
    load = Integral(RealFunction(1.0), v);
    load.assemble();

    EXPECT_NEAR(load.getVector().sum(), 1.0, 1e-12);
  }

  INSTANTIATE_TEST_SUITE_P(
    AllGeometries2D,
    Assembly_P0_2D_Test_8x8,
    ::testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral)
  );

  // =========================================================================
  // 1-D manufactured tests — Segment
  // =========================================================================

  /**
   * @brief 1-D fixture on a Segment mesh (16 cells on [0,1]).
   */
  class Assembly_Segment_Test : public ::testing::Test
  {
    protected:
      void SetUp() override
      {
        m_mesh = Mesh<Context::Local>::UniformGrid(
          Polytope::Type::Segment, { 16 });
        m_mesh.scale(1.0 / 15.0);
        m_mesh.getConnectivity().compute(1, 0);
        m_mesh.getConnectivity().compute(0, 1);
      }

      const Mesh<Context::Local>& getMesh() const { return m_mesh; }

    private:
      Mesh<Context::Local> m_mesh;
  };

  /**
   * @brief 1-D P1 LinearForm vector size equals the FES DOF count.
   */
  TEST_F(Assembly_Segment_Test, LF_VectorSize_P1)
  {
    const auto& mesh = getMesh();
    P1 vh(mesh);
    TestFunction v(vh);

    LinearForm load(v);
    load = Integral(RealFunction(1.0), v);
    load.assemble();

    EXPECT_EQ(load.getVector().size(), static_cast<Eigen::Index>(vh.getSize()));
  }

  /**
   * @brief 1-D P1 load vector sums to domain length = 1.
   */
  TEST_F(Assembly_Segment_Test, LF_LoadVector_SumsToLength)
  {
    const auto& mesh = getMesh();
    P1 vh(mesh);
    TestFunction v(vh);

    LinearForm load(v);
    load = Integral(RealFunction(1.0), v);
    load.assemble();

    EXPECT_NEAR(load.getVector().sum(), 1.0, 1e-12);
  }

  /**
   * @brief 1-D P1 stiffness matrix has correct square dimensions.
   */
  TEST_F(Assembly_Segment_Test, BF_StiffnessMatrix_CorrectDimensions)
  {
    const auto& mesh = getMesh();
    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction  v(vh);

    BilinearForm stiff(u, v);
    stiff = Integral(Grad(u), Grad(v));
    stiff.assemble();

    const auto& A = stiff.getOperator();
    EXPECT_EQ(A.rows(), static_cast<Eigen::Index>(vh.getSize()));
    EXPECT_EQ(A.cols(), static_cast<Eigen::Index>(vh.getSize()));
  }

  /**
   * @brief 1-D P1 mass matrix: (M * ones).sum() equals the domain length.
   *
   * By partition of unity, ∑_i φ_i(x) = 1 for all x, so
   * (M * ones)_i = ∑_j M_ij = ∫_Ω φ_i dΩ, and summing over i gives
   * ∑_i ∫_Ω φ_i dΩ = ∫_Ω 1 dΩ = |Ω| = 1 (after scaling to [0,1]).
   */
  TEST_F(Assembly_Segment_Test, P1MassMatrix_SumEqualsLength)
  {
    const auto& mesh = getMesh();
    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction  v(vh);

    BilinearForm mass(u, v);
    mass = Integral(u, v);
    mass.assemble();

    const auto& M = mass.getOperator();
    Math::Vector<Real> ones(M.cols());
    ones.setOnes();
    EXPECT_NEAR((M * ones).sum(), 1.0, 1e-12);
  }

  /**
   * @brief 1-D preassembled BilinearForm matches inline integral.
   */
  TEST_F(Assembly_Segment_Test, PreassembledBF_MatchesIntegral_P1)
  {
    const auto& mesh = getMesh();
    P1 vh(mesh);

    // Inline
    TrialFunction u1(vh);
    TestFunction  v1(vh);
    BilinearForm bfInline(u1, v1);
    bfInline = Integral(Grad(u1), Grad(v1));
    bfInline.assemble();

    // Preassembled
    TrialFunction u2(vh);
    TestFunction  v2(vh);
    BilinearForm bfPre(u2, v2);
    bfPre = Integral(Grad(u2), Grad(v2));
    bfPre.assemble();

    Problem prob(u2, v2);
    prob = bfPre;
    prob.assemble();

    const auto& A1 = bfInline.getOperator();
    const auto& A2 = prob.getLinearSystem().getOperator();

    ASSERT_EQ(A1.rows(), A2.rows());
    const Math::SparseMatrix<Real> diff = A1 - A2;
    EXPECT_NEAR(diff.norm(), 0.0, 1e-12);
  }

  // =========================================================================
  // 3-D manufactured tests — Hexahedron
  // =========================================================================

  /**
   * @brief 3-D fixture on a small Hexahedron mesh (4×4×4 grid on [0,1]³).
   */
  class Assembly_Hex_Test : public ::testing::Test
  {
    protected:
      void SetUp() override
      {
        m_mesh = Mesh<Context::Local>::UniformGrid(
          Polytope::Type::Hexahedron, { 4, 4, 4 });
        m_mesh.scale(1.0 / 3.0);
        m_mesh.getConnectivity().compute(2, 3);
        m_mesh.getConnectivity().compute(3, 0);
      }

      const Mesh<Context::Local>& getMesh() const { return m_mesh; }

    private:
      Mesh<Context::Local> m_mesh;
  };

  /**
   * @brief 3-D Hexahedron: P1 stiffness matrix has correct square dimensions.
   */
  TEST_F(Assembly_Hex_Test, StiffnessMatrix_CorrectDimensions)
  {
    const auto& mesh = getMesh();
    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction  v(vh);

    BilinearForm stiff(u, v);
    stiff = Integral(Grad(u), Grad(v));
    stiff.assemble();

    const auto& A = stiff.getOperator();
    EXPECT_EQ(A.rows(), static_cast<Eigen::Index>(vh.getSize()));
    EXPECT_EQ(A.cols(), static_cast<Eigen::Index>(vh.getSize()));
  }

  /**
   * @brief 3-D Hexahedron: P1 stiffness matrix is symmetric.
   */
  TEST_F(Assembly_Hex_Test, StiffnessMatrix_IsSymmetric)
  {
    const auto& mesh = getMesh();
    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction  v(vh);

    BilinearForm stiff(u, v);
    stiff = Integral(Grad(u), Grad(v));
    stiff.assemble();

    const auto& A = stiff.getOperator();
    const Math::SparseMatrix<Real> diff = A - Math::SparseMatrix<Real>(A.transpose());
    EXPECT_NEAR(diff.norm(), 0.0, 1e-12);
  }

  /**
   * @brief 3-D Hexahedron: P1 load vector size equals DOF count.
   */
  TEST_F(Assembly_Hex_Test, LoadVector_CorrectSize)
  {
    const auto& mesh = getMesh();
    P1 vh(mesh);
    TestFunction v(vh);

    LinearForm load(v);
    load = Integral(RealFunction(1.0), v);
    load.assemble();

    EXPECT_EQ(load.getVector().size(), static_cast<Eigen::Index>(vh.getSize()));
  }

  /**
   * @brief 3-D Hexahedron Poisson with affine P1-exact solution (f = 0).
   */
  TEST_F(Assembly_Hex_Test, Poisson3D_P1ExactSolution_ZeroError)
  {
    const auto& mesh = getMesh();
    P1 vh(mesh);

    const auto solution = F::x + RealFunction(2.0) * F::y
                        + RealFunction(3.0) * F::z + RealFunction(1.0);
    const auto f = RealFunction(0.0);

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, solution);

    CG(poisson).solve();

    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(u.getSolution() - solution, 2);
    EXPECT_NEAR(Integral(diff).compute(), 0.0, 1e-12);
  }
}
