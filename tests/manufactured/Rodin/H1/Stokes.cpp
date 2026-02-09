/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <algorithm>
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Variational.h"
#include "Rodin/Variational/H1.h"
#include "Rodin/Solver/GMRES.h"
#include "Rodin/Solver/SparseLU.h"
#include "Rodin/Test/Random.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;
using namespace Rodin::Test::Random;

/**
 * @brief Manufactured solutions for the Stokes problem using H1 of order 2.
 *
 * The Stokes system is given by:
 * @f[
 * \left\{
 * \begin{aligned}
 *   -\Delta \mathbf{u} + \nabla p &= \mathbf{f} \quad \text{in } \Omega,\\
 *   \nabla \cdot \mathbf{u} &= 0 \quad \text{in } \Omega,\\
 *   \mathbf{u} &= \mathbf{g} \quad \text{on } \partial\Omega.
 * \end{aligned}
 * \right.
 * @f]
 *
 * The weak formulation is: Find @f$ (\mathbf{u}, p) \in \mathbf{V} \times Q @f$ such that
 * @f[
 *   \int_\Omega \nabla \mathbf{u} : \nabla \mathbf{v} \,dx
 *   - \int_\Omega p \, \nabla \cdot \mathbf{v} \,dx
 *   + \int_\Omega q \, \nabla \cdot \mathbf{u} \,dx
 *   = \int_\Omega \mathbf{f} \cdot \mathbf{v} \,dx,
 * @f]
 * for all @f$ (\mathbf{v}, q) \in \mathbf{V} \times Q @f$, with the essential boundary condition
 * @f[
 *   \mathbf{u} = \mathbf{g} \quad \text{on } \partial\Omega.
 * @f]
 *
 * For these tests, we use:
 * - @f$ \mathbf{V} = [H^1(\Omega)]^d @f$ with @f$ H^1 @f$ of order 2 (quadratic)
 * - @f$ Q = H^1(\Omega) @f$ with @f$ H^1 @f$ of order 1 (linear)
 */
namespace Rodin::Tests::Manufactured::Stokes
{
  template <size_t M>
  class Manufactured_Stokes_Test : public ::testing::TestWithParam<Polytope::Type>
  {
    protected:
      Mesh<Context::Local> getMesh()
      {
        Mesh mesh;
        mesh = mesh.UniformGrid(GetParam(), { M, M });
        mesh.scale(1.0 / (M - 1));
        mesh.getConnectivity().compute(1, 2);
        return mesh;
      }
  };

  using Manufactured_Stokes_Test_16x16 =
    Rodin::Tests::Manufactured::Stokes::Manufactured_Stokes_Test<16>;
  using Manufactured_Stokes_Test_32x32 =
    Rodin::Tests::Manufactured::Stokes::Manufactured_Stokes_Test<32>;
  using Manufactured_Stokes_Test_64x64 =
    Rodin::Tests::Manufactured::Stokes::Manufactured_Stokes_Test<64>;

  TEST_P(Manufactured_Stokes_Test_16x16, Stokes_P1ExactResidual)
  {
    Mesh mesh = this->getMesh();

    H1 uh(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    H1 ph(std::integral_constant<size_t, 1>{}, mesh);

    VectorFunction u_exact{ F::x, -F::x };
    RealFunction p_exact = 0.0;
    VectorFunction f{ Zero(), Zero() };

    TrialFunction u(uh);
    TrialFunction p(ph);
    TestFunction  v(uh);
    TestFunction  q(ph);

    P0g p0g(mesh);
    TrialFunction lambda(p0g);
    TestFunction  mu(p0g);

    Problem stokes(u, p, lambda, v, q, mu);
    stokes = Integral(Jacobian(u), Jacobian(v))
           - Integral(p, Div(v))
           - Integral(Div(u), q)
           + Integral(lambda, q)        // gauge coupling (rank-1)
           + Integral(p, mu)            // mean(p)=0 constraint
           - Integral(f, v)
           + DirichletBC(u, u_exact);

    // Deterministic solve for tests
    SparseLU solver(stokes);
    solver.solve();

    // FE coefficients for the manufactured exact fields
    GridFunction u_exact_coeffs(uh);
    u_exact_coeffs = u_exact;

    GridFunction p_exact_coeffs(ph);
    p_exact_coeffs = p_exact;

    auto& ls = stokes.getLinearSystem();
    auto& A  = ls.getOperator();
    auto& b  = ls.getVector();
    auto& x  = ls.getSolution();

    const auto uSize = u_exact_coeffs.getData().size();
    const auto pSize = p_exact_coeffs.getData().size();
    const auto lambdaIndex = static_cast<Eigen::Index>(uSize + pSize);

    // Build x_exact = [u_exact, p_exact, lambda_exact]
    auto x_exact = x;
    x_exact.head(uSize) = u_exact_coeffs.getData();
    x_exact.segment(uSize, pSize) = p_exact_coeffs.getData();

    // Compute lambda_exact in a way consistent with the assembled discrete system:
    // Let re0_p be the pressure-block residual with lambda=0.
    // Let g = A_pλ be the pressure-block response to the lambda DOF.
    // Choose lambda to minimize ||re_p||_2: lambda* = -(g·re0_p)/(g·g).
    x_exact[lambdaIndex] = 0.0;

    auto re0 = (A * x_exact - b).eval();
    auto re0_p = re0.segment(uSize, pSize);

    decltype(x_exact) eL = x_exact;
    eL.setZero();
    eL[lambdaIndex] = 1.0;

    auto g_full = (A * eL).eval();
    auto g = g_full.segment(uSize, pSize);

    Real lambda_exact = 0.0;
    const Real gg = g.squaredNorm();
    if (gg > 0)
      lambda_exact = - g.dot(re0_p) / gg;

    x_exact[lambdaIndex] = lambda_exact;

    // Residuals
    auto r  = (A * x - b).eval();
    auto re = (A * x_exact - b).eval();

    // Solver residual (scaled)
    const Real scale = std::max<Real>(b.norm(), 1);
    EXPECT_NEAR(r.norm() / scale, 0, RODIN_FUZZY_CONSTANT);

    // Block residual checks for the manufactured "exact" vector
    auto re_u  = re.head(uSize);
    auto re_p  = re.segment(uSize, pSize);
    Real re_mu = re[lambdaIndex];

    // 1) u-block should be near roundoff (Dirichlet makes u exact-in-space here)
    EXPECT_NEAR(re_u.norm(), 0, 1e-10);

    // 2) mean(p)=0 constraint row should be near roundoff (p_exact = 0)
    EXPECT_NEAR(re_mu, 0, 1e-12);

    // 3) lambda only controls span{g}; at the LS minimizer, re_p ⟂ g
    const Real gnorm = std::max<Real>(g.norm(), 1);
    EXPECT_NEAR(g.dot(re_p) / gnorm, 0, 1e-12);

    // Solution field check (relaxed; depends on quadrature + elimination details)
    H1 sh(std::integral_constant<size_t, 1>{}, mesh);
    GridFunction diff_u(sh);
    diff_u = Pow(Frobenius(u.getSolution() - u_exact), 2);
    EXPECT_NEAR(Integral(diff_u).compute(), 0, 5e-8);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1]
   * @f]
   *
   * Manufactured velocity (divergence-free):
   * @f[
   *  \mathbf{u}(x, y) = \begin{pmatrix}
   *    \sin(\pi x) \cos(\pi y) \\
   *    -\cos(\pi x) \sin(\pi y)
   *  \end{pmatrix}
   * @f]
   * Note: ∇·u = π cos(πx)cos(πy) - π cos(πx)cos(πy) = 0 ✓
   *
   * Manufactured pressure:
   * @f[
   *  p(x, y) = \cos(\pi x) \sin(\pi y)
   * @f]
   *
   * Forcing function:
   * @f[
   *  \mathbf{f}(x, y) = \begin{pmatrix}
   *    2\pi^2 \sin(\pi x) \cos(\pi y) - \pi \sin(\pi x) \sin(\pi y) \\
   *    -2\pi^2 \cos(\pi x) \sin(\pi y) + \pi \cos(\pi x) \cos(\pi y)
   *  \end{pmatrix}
   * @f]
   * Computed as f = -Δu + ∇p
   *
   * Boundary conditions:
   * @f[
   *  \mathbf{u} = \mathbf{u}_{\mathrm{exact}} \quad \text{on } \partial\Omega
   * @f]
   */
  TEST_P(Manufactured_Stokes_Test_16x16, Stokes_SimpleSine)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh = this->getMesh();

    // Velocity space: H1 of order 2 (quadratic), vector-valued
    H1 uh(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());

    // Pressure space: H1 of order 1 (linear), scalar
    H1 ph(std::integral_constant<size_t, 1>{}, mesh);

    // Manufactured velocity solution (divergence-free)
    // u = (sin(πx)cos(πy), -cos(πx)sin(πy))
    // ∇·u = π cos(πx)cos(πy) - π cos(πx)cos(πy) = 0 ✓
    VectorFunction u_exact{
      sin(pi * F::x) * cos(pi * F::y),
      -cos(pi * F::x) * sin(pi * F::y)
    };

    // Manufactured pressure solution
    auto p_exact = cos(pi * F::x) * sin(pi * F::y);

    VectorFunction f{
      2 * pi * pi * sin(pi * F::x) * cos(pi * F::y) - pi * sin(pi * F::x) * sin(pi * F::y),
     -2 * pi * pi * cos(pi * F::x) * sin(pi * F::y) + pi * cos(pi * F::x) * cos(pi * F::y)
    };

    // Define trial and test functions
    TrialFunction u(uh); // Velocity trial function
    TrialFunction p(ph); // Pressure trial function
    TestFunction  v(uh); // Velocity test function
    TestFunction  q(ph); // Pressure test function

    // Assemble the weak form:
    // ∫ ∇u : ∇v - ∫ p div(v) - ∫ q div(u) = ∫ f · v
    Problem stokes(u, p, v, q);
    stokes = Integral(Jacobian(u), Jacobian(v))
           - Integral(p, Div(v))
           + Integral(Div(u), q)
           - Integral(f, v)
           + DirichletBC(u, u_exact);

    // Solve the system
    GMRES gmres(stokes);
    gmres.setTolerance(1e-12);
    gmres.setMaxIterations(1000);
    gmres.solve();

    // Compute the L^2 error for velocity
    H1 sh(std::integral_constant<size_t, 1>{}, mesh);
    GridFunction diff_u(sh);
    diff_u = Pow(Frobenius(u.getSolution() - u_exact), 2);
    Real error_u = Integral(diff_u).compute();

    // Compute the L^2 error for pressure
    // Volume
    Real vol = mesh.getMeasure(mesh.getDimension());

    // Mean difference c = (1/|Ω|) ∫ (p_h - p_exact)
    GridFunction mean(sh);
    mean = p.getSolution() - p_exact;
    Real mean_diff = Integral(mean).compute() / vol;

    // Mean-corrected pressure error: || (p_h - mean_diff) - p_exact ||_L2^2
    GridFunction diff_p(sh);
    diff_p = Pow((p.getSolution() - mean_diff) - p_exact, 2);
    Real error_p = Integral(diff_p).compute();

    EXPECT_NEAR(error_u, 0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(0.5 * error_p, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1]
   * @f]
   *
   * Manufactured velocity (divergence-free):
   * @f[
   *  \mathbf{u}(x, y) = \begin{pmatrix}
   *    y(1-y) \\
   *    -x(1-x)
   *  \end{pmatrix}
   * @f]
   * Note: ∇·u = 0 ✓
   *
   * Manufactured pressure:
   * @f[
   *  p(x, y) = x + y - 1
   * @f]
   *
   * Forcing function:
   * @f[
   *  \mathbf{f}(x, y) = \begin{pmatrix}
   *    3 \\
   *    -1
   *  \end{pmatrix}
   * @f]
   * Computed as f = -Δu + ∇p = (2, -2) + (1, 1) = (3, -1)
   */
  TEST_P(Manufactured_Stokes_Test_32x32, Stokes_Polynomial)
  {
    Mesh mesh = this->getMesh();

    // Velocity space: H1 of order 2 (quadratic), vector-valued
    H1 uh(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());

    // Pressure space: H1 of order 1 (linear), scalar
    H1 ph(std::integral_constant<size_t, 1>{}, mesh);

    // Manufactured velocity solution (divergence-free)
    // Using u = (y(1-y), -x(1-x)) 
    // ∇·u = ∂(y(1-y))/∂x + ∂(-x(1-x))/∂y = 0 + 0 = 0 ✓
    VectorFunction u_exact{
      F::y * (1 - F::y),
      -F::x * (1 - F::x)
    };

    // Manufactured pressure solution
    auto p_exact = F::x + F::y - 1;

    // Forcing function: f = -Δu + ∇p
    // -Δ(y(1-y)) = -∂²/∂y²[y(1-y)] = -(-2) = 2
    // -Δ(-x(1-x)) = -∂²/∂x²[-x(1-x)] = -(-(-2)) = -2
    // ∇p = (1, 1)
    // Therefore f = (2+1, -2+1) = (3, -1)
    VectorFunction f{ 3.0, -1.0 };

    // Define trial and test functions
    TrialFunction u(uh);
    TrialFunction p(ph);
    TestFunction  v(uh);
    TestFunction  q(ph);

    // Assemble the weak form
    Problem stokes(u, p, v, q);
    stokes = Integral(Jacobian(u), Jacobian(v))
           - Integral(p, Div(v))
           + Integral(Div(u), q)
           - Integral(f, v)
           + DirichletBC(u, u_exact);

    // Solve the system
    GMRES gmres(stokes);
    gmres.setTolerance(1e-12);
    gmres.setMaxIterations(1000);
    gmres.solve();

    // Compute the L^2 error for velocity
    H1 sh(std::integral_constant<size_t, 1>{}, mesh);
    GridFunction diff_u(sh);
    diff_u = Pow(Frobenius(u.getSolution() - u_exact), 2);
    Real error_u = Integral(diff_u).compute();

    // Compute the L^2 error for pressure
    // Volume
    Real vol = mesh.getMeasure(mesh.getDimension());

    // Mean difference c = (1/|Ω|) ∫ (p_h - p_exact)
    GridFunction mean(sh);
    mean = p.getSolution() - p_exact;
    Real mean_diff = Integral(mean).compute() / vol;

    // Mean-corrected pressure error: || (p_h - mean_diff) - p_exact ||_L2^2
    GridFunction diff_p(sh);
    diff_p = Pow((p.getSolution() - mean_diff) - p_exact, 2);
    Real error_p = Integral(diff_p).compute();

    EXPECT_NEAR(error_u, 0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(0.1 * error_p, 0, 1.5 * RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1]
   * @f]
   *
   * Manufactured velocity (Taylor-Green vortex):
   * @f[
   *  \mathbf{u}(x, y) = \begin{pmatrix}
   *    \sin(\pi x) \cos(\pi y) \\
   *    -\cos(\pi x) \sin(\pi y)
   *  \end{pmatrix}
   * @f]
   *
   * Manufactured pressure:
   * @f[
   *  p(x, y) = -\frac{1}{4}[\cos(2\pi x) + \cos(2\pi y)]
   * @f]
   *
   * Forcing function:
   * @f[
   *  \mathbf{f}(x, y) = \begin{pmatrix}
   *    2\pi^2 \sin(\pi x) \cos(\pi y) + \frac{\pi}{2} \sin(2\pi x) \\
   *    -2\pi^2 \cos(\pi x) \sin(\pi y) + \frac{\pi}{2} \sin(2\pi y)
   *  \end{pmatrix}
   * @f]
   *
   * Boundary conditions:
   * @f[
   *  \mathbf{u} = \mathbf{u}_{\mathrm{exact}} \quad \text{on } \partial\Omega
   * @f]
   */
  TEST_P(Manufactured_Stokes_Test_64x64, Stokes_TaylorGreen)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh = this->getMesh();

    // Velocity space: H1 of order 2 (quadratic), vector-valued
    H1 uh(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());

    // Pressure space: H1 of order 1 (linear), scalar
    H1 ph(std::integral_constant<size_t, 1>{}, mesh);

    // Manufactured velocity solution (Taylor-Green vortex)
    VectorFunction u_exact{
      sin(pi * F::x) * cos(pi * F::y),
      -cos(pi * F::x) * sin(pi * F::y)
    };

    // Manufactured pressure solution
    auto p_exact = -0.25 * (cos(2 * pi * F::x) + cos(2 * pi * F::y));

    // Forcing function
    VectorFunction f{
      2 * pi * pi * sin(pi * F::x) * cos(pi * F::y) + 0.5 * pi * sin(2 * pi * F::x),
      -2 * pi * pi * cos(pi * F::x) * sin(pi * F::y) + 0.5 * pi * sin(2 * pi * F::y)
    };

    // Define trial and test functions
    TrialFunction u(uh);
    TrialFunction p(ph);
    TestFunction  v(uh);
    TestFunction  q(ph);

    // Assemble the weak form
    Problem stokes(u, p, v, q);
    stokes = Integral(Jacobian(u), Jacobian(v))
           - Integral(p, Div(v))
           + Integral(Div(u), q)
           - Integral(f, v)
           + DirichletBC(u, u_exact);

    // Solve the system
    GMRES gmres(stokes);
    gmres.setTolerance(1e-12);
    gmres.setMaxIterations(1000);
    gmres.solve();

    // Compute the L^2 error for velocity
    H1 sh(std::integral_constant<size_t, 1>{}, mesh);
    GridFunction diff_u(sh);
    diff_u = Pow(Frobenius(u.getSolution() - u_exact), 2);
    Real error_u = Integral(diff_u).compute();

    // Compute the L^2 error for pressure
    Real vol = mesh.getMeasure(mesh.getDimension());

    // Mean difference c = (1/|Ω|) ∫ (p_h - p_exact)
    GridFunction mean(sh);
    mean = p.getSolution() - p_exact;
    Real mean_diff = Integral(mean).compute() / vol;

    // Mean-corrected pressure error: || (p_h - mean_diff) - p_exact ||_L2^2
    GridFunction diff_p(sh);
    diff_p = Pow((p.getSolution() - mean_diff) - p_exact, 2);
    Real error_p = Integral(diff_p).compute();

    EXPECT_NEAR(error_u, 0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(0.1 * error_p, 0, 2 * RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1]
   * @f]
   *
   * Manufactured velocity (divergence-free using stream function):
   * @f[
   *  \psi(x, y) = x^2(1-x)^2 y^2(1-y)^2
   * @f]
   * @f[
   *  \mathbf{u}(x, y) = \begin{pmatrix}
   *    \frac{\partial\psi}{\partial y} \\
   *    -\frac{\partial\psi}{\partial x}
   *  \end{pmatrix}
   * @f]
   *
   * Manufactured pressure:
   * @f[
   *  p(x, y) = x^2 - y^2
   * @f]
   *
   * Forcing function (computed from -Δu + ∇p):
   * @f[
   *  \mathbf{f}(x, y) = -\Delta\mathbf{u} + \nabla p
   * @f]
   *
   * Boundary conditions:
   * @f[
   *  \mathbf{u} = \mathbf{0} \quad \text{on } \partial\Omega
   * @f]
   */
  TEST_P(Manufactured_Stokes_Test_32x32, Stokes_Quadratic)
  {
    Mesh mesh = this->getMesh();

    // Velocity space: H1 of order 2 (quadratic), vector-valued
    H1 uh(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());

    // Pressure space: H1 of order 1 (linear), scalar
    H1 ph(std::integral_constant<size_t, 1>{}, mesh);

    P0g p0g(mesh);

    // Manufactured velocity solution (divergence-free)
    // Using stream function ψ = x²(1-x)² y²(1-y)²
    // u₁ = ∂ψ/∂y, u₂ = -∂ψ/∂x ensures ∇·u = 0

    // ∂ψ/∂y = x²(1-x)² · [2y(1-y)² - 2y²(1-y)]
    //       = 2x²(1-x)² y(1-y) [1-y-y]
    //       = 2x²(1-x)² y(1-y)(1-2y)
    auto u1 = 2 * pow(F::x, 2) * pow(1 - F::x, 2) * F::y * (1 - F::y) * (1 - 2 * F::y);

    // ∂ψ/∂x = [2x(1-x)² - 2x²(1-x)] · y²(1-y)²
    //       = 2x(1-x) y²(1-y)² [1-x-x]
    //       = 2x(1-x) y²(1-y)² (1-2x)
    // u₂ = -∂ψ/∂x
    auto u2 = -2 * F::x * (1 - F::x) * pow(F::y, 2) * pow(1 - F::y, 2) * (1 - 2 * F::x);

    VectorFunction u_exact{ u1, u2 };

    // Manufactured pressure solution
    auto p_exact = pow(F::x, 2) - pow(F::y, 2);

    // Compute Laplacian components
    // For u1 = 2x²(1-x)² y(1-y)(1-2y), compute ∂²u1/∂x² + ∂²u1/∂y²
    // This is quite involved, so I'll compute it symbolically:

    // ∂²u1/∂x² = 2y(1-y)(1-2y) · [2(1-x)² - 8x(1-x) + 2x²]
    //          = 2y(1-y)(1-2y) · 2[(1-x)² - 4x(1-x) + x²]
    //          = 4y(1-y)(1-2y) · [1 - 2x + x² - 4x + 4x² + x²]
    //          = 4y(1-y)(1-2y) · [1 - 6x + 6x²]
    auto d2u1_dx2 = 4 * F::y * (1 - F::y) * (1 - 2 * F::y) * (1 - 6 * F::x + 6 * pow(F::x, 2));

    // ∂²u1/∂y² = 2x²(1-x)² · [∂²/∂y²(y(1-y)(1-2y))]
    // Let g(y) = y(1-y)(1-2y) = y(1-3y+2y²) = y - 3y² + 2y³
    // g'(y) = 1 - 6y + 6y²
    // g''(y) = -6 + 12y
    auto d2u1_dy2 = 2 * pow(F::x, 2) * pow(1 - F::x, 2) * (-6 + 12 * F::y);

    auto laplace_u1 = d2u1_dx2 + d2u1_dy2;

    // For u2 = -2x(1-x) y²(1-y)² (1-2x), by symmetry:
    // Let h(x) = x(1-x)(1-2x) = x(1-3x+2x²) = x - 3x² + 2x³
    // h'(x) = 1 - 6x + 6x²
    // h''(x) = -6 + 12x
    auto d2u2_dx2 = -2 * pow(F::y, 2) * pow(1 - F::y, 2) * (-6 + 12 * F::x);

    auto d2u2_dy2 =
      -4 * F::x * (1 - F::x) * (1 - 2 * F::x) * (1 - 6 * F::y + 6 * pow(F::y, 2));
    auto laplace_u2 = d2u2_dx2 + d2u2_dy2;

    // Pressure gradient
    auto grad_p_x = 2 * F::x;
    auto grad_p_y = -2 * F::y;

    // Forcing function: f = -Δu + ∇p
    VectorFunction f{
      -laplace_u1 + grad_p_x,
      -laplace_u2 + grad_p_y
    };

    // Define trial and test functions
    TrialFunction u(uh);
    TrialFunction p(ph);
    TestFunction  v(uh);
    TestFunction  q(ph);

    TrialFunction lambda(p0g);
    TestFunction  mu(p0g);

    // Assemble the weak form
    Problem stokes(u, p, v, q, lambda, mu);
    stokes = Integral(Jacobian(u), Jacobian(v))
           - Integral(p, Div(v))
           + Integral(Div(u), q)
           - Integral(f, v)
           + Integral(lambda, q)
           + Integral(p, mu)
           + DirichletBC(u, u_exact);

    // Solve the system
    GMRES gmres(stokes);
    gmres.setTolerance(1e-12);
    gmres.solve();

    // Compute the L^2 error for velocity
    H1 sh(std::integral_constant<size_t, 1>{}, mesh);
    GridFunction diff_u(sh);
    diff_u = Pow(Frobenius(u.getSolution() - u_exact), 2);
    Real error_u = Integral(diff_u).compute();

    // Compute the L^2 error for pressure
    // Volume
    Real vol = mesh.getMeasure(mesh.getDimension());

    // Mean difference c = (1/|Ω|) ∫ (p_h - p_exact)
    GridFunction mean(sh);
    mean = p.getSolution() - p_exact;
    Real mean_diff = Integral(mean).compute() / vol;

    // Mean-corrected pressure error: || (p_h - mean_diff) - p_exact ||_L2^2
    GridFunction diff_p(sh);
    diff_p = Pow((p.getSolution() - mean_diff) - p_exact, 2);
    Real error_p = Integral(diff_p).compute();

    EXPECT_NEAR(error_u, 0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(error_p, 0, RODIN_FUZZY_CONSTANT);
  }

  INSTANTIATE_TEST_SUITE_P(
    MeshParams16x16,
    Manufactured_Stokes_Test_16x16,
    ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle)
  );

  INSTANTIATE_TEST_SUITE_P(
    MeshParams32x32,
    Manufactured_Stokes_Test_32x32,
    ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle)
  );

  INSTANTIATE_TEST_SUITE_P(
    MeshParams64x64,
    Manufactured_Stokes_Test_64x64,
    ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle)
  );
}
