/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 *
 * Full verification harness for Stokes(u,p,lambda) in PETSc:
 *  - solves (optionally)
 *  - checks x and b for NaN/Inf
 *  - computes ||Ax-b||/||b||
 *  - builds x_exact = [u_exact, p_exact, lambda_exact] and computes ||A x_exact - b||/||b||
 *  - prints blockwise residual norms (u,p,lambda)
 *  - checks the coupling block adjointness with the CORRECT SIGN depending on your weak form
 *
 * IMPORTANT: choose the sign convention below by setting COUPLING_SIGN.
 */

#include "Rodin/Geometry/Polytope.h"
#include <Rodin/PETSc.h>

#include <Rodin/Types.h>
#include <Rodin/Solver.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

#include <petscksp.h>
#include <petscmat.h>

#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <string>

// ------------------------------------------------------------
// Choose the expected relation between Aup and Apu^T.
//
// If your weak form contains:
//   - Integral(p, Div(v))   + Integral(Div(u), q)
// then: Aup = -Apu^T  (expected anti-adjointness)
//
// If your weak form contains:
//   + Integral(p, Div(v))   + Integral(Div(u), q)
// then: Aup = +Apu^T  (expected adjointness)
//
// Set COUPLING_SIGN = -1 for Aup + Apu^T ≈ 0
// Set COUPLING_SIGN = +1 for Aup - Apu^T ≈ 0
// ------------------------------------------------------------
static constexpr int COUPLING_SIGN = -1; // <-- set to -1 for standard Stokes: -p div(v) + div(u) q

using namespace Rodin;
using namespace Rodin::Math;
using namespace Rodin::Solver;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

// ---------------------------- Utilities ----------------------------

static void CheckPetsc(PetscErrorCode ierr)
{
  if (ierr != PETSC_SUCCESS)
  {
    PetscError(PETSC_COMM_SELF, __LINE__, "CheckPetsc", __FILE__, ierr, PETSC_ERROR_REPEAT, "PETSc call failed");
    std::abort();
  }
}

// return true if Vec contains NaN/Inf (cheap-ish: check norm/min/max)
static bool VecHasNaNInf(Vec v)
{
  PetscReal n2 = 0, vmin = 0, vmax = 0;
  PetscInt  imin = 0, imax = 0;

  CheckPetsc(VecNorm(v, NORM_2, &n2));
  CheckPetsc(VecMin(v, &imin, &vmin));
  CheckPetsc(VecMax(v, &imax, &vmax));

  auto bad = [](double x)
  {
    return std::isnan(x) || std::isinf(x);
  };

  return bad((double)n2) || bad((double)vmin) || bad((double)vmax);
}

static double RelativeResidual(Mat A, Vec x, Vec b)
{
  Vec r = nullptr;
  CheckPetsc(VecDuplicate(b, &r));
  CheckPetsc(MatMult(A, x, r));
  CheckPetsc(VecAXPY(r, -1.0, b));

  PetscReal nr = 0, nb = 0;
  CheckPetsc(VecNorm(r, NORM_2, &nr));
  CheckPetsc(VecNorm(b, NORM_2, &nb));

  CheckPetsc(VecDestroy(&r));
  if (nb > 0) return (double)(nr / nb);
  return (double)nr;
}

static void PrintVecStats(const char* name, Vec v)
{
  PetscReal n2 = 0, vmin = 0, vmax = 0;
  PetscInt  imin = 0, imax = 0;
  CheckPetsc(VecNorm(v, NORM_2, &n2));
  CheckPetsc(VecMin(v, &imin, &vmin));
  CheckPetsc(VecMax(v, &imax, &vmax));

  std::cout << name << ": ||.||2=" << std::scientific << (double)n2
            << "  min=" << (double)vmin << "  max=" << (double)vmax
            << "  (imin=" << imin << ", imax=" << imax << ")\n";
}

static double FrobeniusNorm(Mat M)
{
  PetscReal nF = 0;
  CheckPetsc(MatNorm(M, NORM_FROBENIUS, &nF));
  return (double)nF;
}

// Extract submatrix A[rows, cols] into a NEW matrix (caller destroys)
static Mat SubMat(Mat A, IS rows, IS cols)
{
  Mat S = nullptr;
  CheckPetsc(MatCreateSubMatrix(A, rows, cols, MAT_INITIAL_MATRIX, &S));
  return S;
}

static void CouplingAdjointnessCheck(Mat A, PetscInt nu, PetscInt np)
{
  // We assume global ordering [u][p][lambda]
  IS is_u = nullptr, is_p = nullptr;
  CheckPetsc(ISCreateStride(PETSC_COMM_SELF, nu, 0, 1, &is_u));
  CheckPetsc(ISCreateStride(PETSC_COMM_SELF, np, nu, 1, &is_p));

  // Aup: rows=u, cols=p
  Mat Aup = SubMat(A, is_u, is_p);
  // Apu: rows=p, cols=u
  Mat Apu = SubMat(A, is_p, is_u);

  // Compute E:
  //   if COUPLING_SIGN == -1: E = Aup + Apu^T   (expected ~0 for standard Stokes)
  //   if COUPLING_SIGN == +1: E = Aup - Apu^T   (expected ~0 for same-sign form)
  Mat At = nullptr;
  CheckPetsc(MatTranspose(Apu, MAT_INITIAL_MATRIX, &At));

  Mat E = nullptr;
  CheckPetsc(MatDuplicate(Aup, MAT_COPY_VALUES, &E));
  if constexpr (COUPLING_SIGN == -1)
  {
    // E = Aup + At
    CheckPetsc(MatAXPY(E, 1.0, At, DIFFERENT_NONZERO_PATTERN));
  }
  else
  {
    // E = Aup - At
    CheckPetsc(MatAXPY(E, -1.0, At, DIFFERENT_NONZERO_PATTERN));
  }

  const double nAup = FrobeniusNorm(Aup);
  const double nAt  = FrobeniusNorm(At);
  const double nE   = FrobeniusNorm(E);

  const double denom = nAup + nAt;
  const double rel = (denom > 0) ? (nE / denom) : nE;

  std::cout << "||Apu||_F=" << std::scientific << FrobeniusNorm(Apu)
            << "  ||Aup||_F=" << nAup << "\n";

  if constexpr (COUPLING_SIGN == -1)
  {
    std::cout << "Coupling check (standard Stokes expects Aup = -Apu^T): "
              << "||Aup + Apu^T||_F / (||Aup||_F+||Apu^T||_F) = "
              << rel << "\n";
  }
  else
  {
    std::cout << "Coupling check (same-sign expects Aup = +Apu^T): "
              << "||Aup - Apu^T||_F / (||Aup||_F+||Apu^T||_F) = "
              << rel << "\n";
  }

  CheckPetsc(MatDestroy(&E));
  CheckPetsc(MatDestroy(&At));
  CheckPetsc(MatDestroy(&Aup));
  CheckPetsc(MatDestroy(&Apu));
  CheckPetsc(ISDestroy(&is_u));
  CheckPetsc(ISDestroy(&is_p));
}

static void BlockResidualNorms(Mat A, Vec x, Vec b, PetscInt nu, PetscInt np, PetscInt nl)
{
  Vec r = nullptr;
  CheckPetsc(VecDuplicate(b, &r));
  CheckPetsc(MatMult(A, x, r));
  CheckPetsc(VecAXPY(r, -1.0, b));

  IS is_u = nullptr, is_p = nullptr, is_l = nullptr;
  CheckPetsc(ISCreateStride(PETSC_COMM_SELF, nu, 0, 1, &is_u));
  CheckPetsc(ISCreateStride(PETSC_COMM_SELF, np, nu, 1, &is_p));
  CheckPetsc(ISCreateStride(PETSC_COMM_SELF, nl, nu + np, 1, &is_l));

  auto print_block = [&](IS is, const char* name)
  {
    Vec sub = nullptr;
    CheckPetsc(VecGetSubVector(r, is, &sub));
    PetscReal n2 = 0;
    CheckPetsc(VecNorm(sub, NORM_2, &n2));
    CheckPetsc(VecRestoreSubVector(r, is, &sub));
    std::cout << "||r_" << name << "||_2 = " << std::scientific << (double)n2 << "\n";
  };

  print_block(is_u, "u");
  print_block(is_p, "p");
  print_block(is_l, "lambda");

  CheckPetsc(ISDestroy(&is_u));
  CheckPetsc(ISDestroy(&is_p));
  CheckPetsc(ISDestroy(&is_l));
  CheckPetsc(VecDestroy(&r));
}

static Vec BuildExactVector(
  Vec template_x,
  PetscInt nu, PetscInt np, PetscInt nl,
  Vec u_data, Vec p_data,
  PetscScalar lambda_value)
{
  Vec xe = nullptr;
  CheckPetsc(VecDuplicate(template_x, &xe));
  CheckPetsc(VecSet(xe, 0.0));

  IS is_u = nullptr, is_p = nullptr, is_l = nullptr;
  CheckPetsc(ISCreateStride(PETSC_COMM_SELF, nu, 0, 1, &is_u));
  CheckPetsc(ISCreateStride(PETSC_COMM_SELF, np, nu, 1, &is_p));
  CheckPetsc(ISCreateStride(PETSC_COMM_SELF, nl, nu + np, 1, &is_l));

  Vec xe_u = nullptr, xe_p = nullptr, xe_l = nullptr;
  CheckPetsc(VecGetSubVector(xe, is_u, &xe_u));
  CheckPetsc(VecGetSubVector(xe, is_p, &xe_p));
  CheckPetsc(VecGetSubVector(xe, is_l, &xe_l));

  CheckPetsc(VecCopy(u_data, xe_u));
  CheckPetsc(VecCopy(p_data, xe_p));
  CheckPetsc(VecSet(xe_l, lambda_value));

  CheckPetsc(VecRestoreSubVector(xe, is_u, &xe_u));
  CheckPetsc(VecRestoreSubVector(xe, is_p, &xe_p));
  CheckPetsc(VecRestoreSubVector(xe, is_l, &xe_l));

  CheckPetsc(ISDestroy(&is_u));
  CheckPetsc(ISDestroy(&is_p));
  CheckPetsc(ISDestroy(&is_l));

  return xe;
}

// ---------------------------- Main ----------------------------

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {12, 12, 12});
  mesh.scale(1.0 / (12 - 1));
  mesh.getConnectivity().compute(2, 3);
  mesh.getConnectivity().compute(3, 2);
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);

  H1 uh(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
  H1 ph(std::integral_constant<size_t, 1>{}, mesh);
  P0g p0g(mesh);

  std::cout << "Vector FES size: " << uh.getSize() << "\n";
  std::cout << "Scalar FES size: " << ph.getSize() << "\n";
  std::cout << "P0g   FES size: " << p0g.getSize() << "\n";
  std::cout << "Assembling...\n";

  {
    auto pi = Rodin::Math::Constants::pi();

    // Your original trig manufactured solution (NOT in FE space)
    VectorFunction u_exact{
      Sin(pi * F::x) * Cos(pi * F::y) * Cos(pi * F::z),
      -Cos(pi * F::x) * Sin(pi * F::y) * Cos(pi * F::z),
      Zero()
    };
    auto p_exact = Cos(2 * pi * F::x) * Cos(2 * pi * F::y) * Cos(2 * pi * F::z);

    VectorFunction f{
      3 * pi * pi * Sin(pi * F::x) * Cos(pi * F::y) * Cos(pi * F::z)
      - 2 * pi * Sin(2 * pi * F::x) * Cos(2 * pi * F::y) * Cos(2 * pi * F::z),

      -3 * pi * pi * Cos(pi * F::x) * Sin(pi * F::y) * Cos(pi * F::z)
      - 2 * pi * Cos(2 * pi * F::x) * Sin(2 * pi * F::y) * Cos(2 * pi * F::z),

      -2 * pi * Cos(2 * pi * F::x) * Cos(2 * pi * F::y) * Sin(2 * pi * F::z)
    };

    // Unknowns/tests
    PETSc::Variational::TrialFunction u(uh);   u.setName("u");
    PETSc::Variational::TrialFunction p(ph);   p.setName("p");
    PETSc::Variational::TrialFunction lambda(p0g); lambda.setName("lambda");

    PETSc::Variational::TestFunction v(uh);
    PETSc::Variational::TestFunction q(ph);
    PETSc::Variational::TestFunction mu(p0g);

    Problem stokes(u, p, lambda, v, q, mu);

    // Standard Stokes + mean-pressure constraint (LM) + Dirichlet on u
    stokes =
        Integral(Jacobian(u), Jacobian(v))
      - Integral(p, Div(v))
      + Integral(Div(u), q)
      + Integral(lambda, q)   // λ ∫ q
      + Integral(p, mu)       // μ ∫ p
      - Integral(f, v)
      + DirichletBC(u, u_exact);

    auto t0 = std::chrono::high_resolution_clock::now();
    stokes.assemble();
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Assembly time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
              << " ms\n";

    // Access linear system
    ::Mat& A = stokes.getLinearSystem().getOperator();
    ::Vec& x = stokes.getLinearSystem().getSolution();
    ::Vec& b = stokes.getLinearSystem().getVector(); // RHS (critical)

    // Ensure assembled (usually already done)
    CheckPetsc(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
    CheckPetsc(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
    CheckPetsc(VecAssemblyBegin(b));
    CheckPetsc(VecAssemblyEnd(b));

    // Solve using command line options (works with -ksp_type/-pc_type etc.)
    // If you want Rodin's GMRES wrapper, keep it; here is the raw PETSc path:
    {
      Solver::KSP(stokes).solve();

      // ::KSP ksp = nullptr;
      // CheckPetsc(KSPCreate(PETSC_COMM_SELF, &ksp));
      // CheckPetsc(KSPSetOperators(ksp, A, A));
      // CheckPetsc(KSPSetFromOptions(ksp));
      // CheckPetsc(KSPSolve(ksp, b, x));

      // KSPConvergedReason reason;
      // CheckPetsc(KSPGetConvergedReason(ksp, &reason));
      // std::cout << "KSP reason = " << reason << "\n";

      // CheckPetsc(KSPDestroy(&ksp));
    }

    // Basic sanity on solution/RHS
    PrintVecStats("x", x);
    PrintVecStats("b", b);

    const bool x_bad = VecHasNaNInf(x);
    const bool b_bad = VecHasNaNInf(b);
    if (x_bad || b_bad)
    {
      std::cout << "WARNING: "
                << (x_bad ? "x has NaN/Inf " : "")
                << (b_bad ? "b has NaN/Inf " : "")
                << "\n";
    }

    // Block sizes: [u][p][lambda]
    const PetscInt nu = (PetscInt)uh.getSize();
    const PetscInt np = (PetscInt)ph.getSize();
    const PetscInt nl = 1;

    // 1) Coupling sign check (correct sign controlled by COUPLING_SIGN)
    CouplingAdjointnessCheck(A, nu, np);

    // 2) Residual after solve
    const double rel_res_solve = RelativeResidual(A, x, b);
    std::cout << "||A x - b||/||b|| = " << std::scientific << rel_res_solve << "\n";
    BlockResidualNorms(A, x, b, nu, np, nl);

    // 3) Build x_exact = [u_exact, p_exact, 0] and check assembly consistency
    PETSc::Variational::GridFunction u_e(uh);
    u_e = u_exact;

    PETSc::Variational::GridFunction p_e(ph);
    p_e = p_exact;

    // IMPORTANT:
    // - With LM enforcing mean(p)=0, you generally want p_e to have zero mean too.
    // - For the trig p_exact above, mean is not guaranteed to be exactly 0 on your discrete mesh.
    //   This can contribute to a nonzero residual in the (p,mu)/(lambda,q) parts.
    // - Still, your earlier diagnostics showed the u-block residual dominates, so focus there.


    Vec xe = BuildExactVector(
      x, nu, np, nl,
      u_e.getData(),
      p_e.getData(),
      (PetscScalar)0.0
    );

    mesh.save("Stokes.mesh");
    u.getSolution().save("Stokes_velocity.gf");
    p.getSolution().save("Stokes_pressure.gf");



    PrintVecStats("x_exact", xe);
    const double rel_res_exact = RelativeResidual(A, xe, b);
    std::cout << "||A x_exact - b||/||b|| = " << std::scientific << rel_res_exact << "\n";
    BlockResidualNorms(A, xe, b, nu, np, nl);

    CheckPetsc(VecDestroy(&xe));
  }

  PetscFinalize();
  return 0;
}
