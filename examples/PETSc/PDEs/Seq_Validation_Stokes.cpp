/*
 * Standard Stokes validation harness (P2–P1), PETSc-friendly, with BC-masked assembly check.
 *
 * Goals:
 *  1) Use an MMS that is exactly representable by P2–P1:
 *        u_exact: quadratic polynomial  (P2)
 *        p_exact: linear polynomial     (P1)
 *     so the discrete system should be satisfied to roundoff ON FREE DOFS
 *     (Dirichlet rows are masked out for the assembly check).
 *
 *  2) No extra LM variable. Pressure is handled via PETSc nullspace (constant mode).
 *     We also remove the nullspace component from b: MatNullSpaceRemove(nsp,b).
 *
 *  3) Validate assembly:
 *        ||A x_exact - b|| / ||b||  on FREE DOFS (masked Dirichlet velocity dofs)
 *     Validate solve:
 *        ||A x_h - b|| / ||b||      on FREE DOFS (same mask)
 *        and print KSP reason via -ksp_converged_reason
 *
 * Notes:
 *  - This harness assumes strong enforcement of Dirichlet BC overwrites rows; hence masking.
 *  - Dirichlet dofs are obtained from DirichletBC::getDOFs().
 *
 * Run examples:
 *  - Baseline:  -ksp_type gmres -pc_type none -ksp_rtol 1e-10 -ksp_monitor_true_residual -ksp_converged_reason
 *  - Better:    -ksp_type gmres -pc_type gamg -ksp_rtol 1e-10 -ksp_monitor_true_residual -ksp_converged_reason
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

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using namespace Rodin;
using namespace Rodin::Math;
using namespace Rodin::Solver;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

static void CheckPetsc(PetscErrorCode ierr)
{
  if (ierr != PETSC_SUCCESS)
  {
    PetscError(PETSC_COMM_SELF, __LINE__, "CheckPetsc", __FILE__, ierr,
               PETSC_ERROR_REPEAT, "PETSc call failed");
    std::abort();
  }
}

static bool VecHasNaNInf(Vec v)
{
  PetscReal n2 = 0, vmin = 0, vmax = 0;
  PetscInt  imin = 0, imax = 0;

  CheckPetsc(VecNorm(v, NORM_2, &n2));
  CheckPetsc(VecMin(v, &imin, &vmin));
  CheckPetsc(VecMax(v, &imax, &vmax));

  auto bad = [](double x) { return std::isnan(x) || std::isinf(x); };
  return bad((double)n2) || bad((double)vmin) || bad((double)vmax);
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

static void ZeroEntries(Vec v, const std::vector<PetscInt>& idx)
{
  for (PetscInt i : idx)
    CheckPetsc(VecSetValue(v, i, 0.0, INSERT_VALUES));
  CheckPetsc(VecAssemblyBegin(v));
  CheckPetsc(VecAssemblyEnd(v));
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
  return (nb > 0) ? (double)(nr / nb) : (double)nr;
}

static double RelativeResidualMasked(Mat A, Vec x, Vec b, const std::vector<PetscInt>& mask)
{
  Vec r = nullptr, bb = nullptr;
  CheckPetsc(VecDuplicate(b, &r));
  CheckPetsc(VecDuplicate(b, &bb));

  // r = A x - b
  CheckPetsc(MatMult(A, x, r));
  CheckPetsc(VecAXPY(r, -1.0, b));

  // bb = b
  CheckPetsc(VecCopy(b, bb));

  // Zero masked entries in both vectors, so the norm checks only free dofs
  ZeroEntries(r, mask);
  ZeroEntries(bb, mask);

  PetscReal nr = 0, nb = 0;
  CheckPetsc(VecNorm(r, NORM_2, &nr));
  CheckPetsc(VecNorm(bb, NORM_2, &nb));

  CheckPetsc(VecDestroy(&r));
  CheckPetsc(VecDestroy(&bb));
  return (nb > 0) ? (double)(nr / nb) : (double)nr;
}

static void BlockResidualNormsMasked_Up(
  Mat A, Vec x, Vec b,
  PetscInt nu, PetscInt np,
  const std::vector<PetscInt>& mask_u_global)
{
  Vec r = nullptr;
  CheckPetsc(VecDuplicate(b, &r));
  CheckPetsc(MatMult(A, x, r));
  CheckPetsc(VecAXPY(r, -1.0, b));

  // mask constrained velocity dofs (global indices)
  ZeroEntries(r, mask_u_global);

  IS is_u = nullptr, is_p = nullptr;
  CheckPetsc(ISCreateStride(PETSC_COMM_SELF, nu, 0, 1, &is_u));
  CheckPetsc(ISCreateStride(PETSC_COMM_SELF, np, nu, 1, &is_p));

  auto print_block = [&](IS is, const char* name)
  {
    Vec sub = nullptr;
    CheckPetsc(VecGetSubVector(r, is, &sub));
    PetscReal n2 = 0;
    CheckPetsc(VecNorm(sub, NORM_2, &n2));
    CheckPetsc(VecRestoreSubVector(r, is, &sub));
    std::cout << "||r_" << name << "||_2 (masked) = " << std::scientific << (double)n2 << "\n";
  };

  print_block(is_u, "u");
  print_block(is_p, "p");

  CheckPetsc(ISDestroy(&is_u));
  CheckPetsc(ISDestroy(&is_p));
  CheckPetsc(VecDestroy(&r));
}

static Vec BuildExactVector_Up(Vec template_x, PetscInt nu, PetscInt np, Vec u_data, Vec p_data)
{
  Vec xe = nullptr;
  CheckPetsc(VecDuplicate(template_x, &xe));
  CheckPetsc(VecSet(xe, 0.0));

  IS is_u = nullptr, is_p = nullptr;
  CheckPetsc(ISCreateStride(PETSC_COMM_SELF, nu, 0, 1, &is_u));
  CheckPetsc(ISCreateStride(PETSC_COMM_SELF, np, nu, 1, &is_p));

  Vec xe_u = nullptr, xe_p = nullptr;
  CheckPetsc(VecGetSubVector(xe, is_u, &xe_u));
  CheckPetsc(VecGetSubVector(xe, is_p, &xe_p));

  CheckPetsc(VecCopy(u_data, xe_u));
  CheckPetsc(VecCopy(p_data, xe_p));

  CheckPetsc(VecRestoreSubVector(xe, is_u, &xe_u));
  CheckPetsc(VecRestoreSubVector(xe, is_p, &xe_p));

  CheckPetsc(ISDestroy(&is_u));
  CheckPetsc(ISDestroy(&is_p));

  return xe;
}

static MatNullSpace AttachPressureNullspace(Mat A, PetscInt nu, PetscInt np)
{
  // Nullspace: constant pressure mode (p += c), velocity part is zero.
  Vec ns = nullptr;
  CheckPetsc(MatCreateVecs(A, &ns, nullptr));
  CheckPetsc(VecSet(ns, 0.0));

  for (PetscInt i = 0; i < np; ++i)
    CheckPetsc(VecSetValue(ns, nu + i, 1.0, INSERT_VALUES));

  CheckPetsc(VecAssemblyBegin(ns));
  CheckPetsc(VecAssemblyEnd(ns));
  CheckPetsc(VecNormalize(ns, nullptr));

  MatNullSpace nsp = nullptr;
  CheckPetsc(MatNullSpaceCreate(PETSC_COMM_SELF, PETSC_FALSE, 1, &ns, &nsp));
  CheckPetsc(MatSetNullSpace(A, nsp));

  CheckPetsc(VecDestroy(&ns));
  return nsp; // caller destroys
}

template <class ScalarMap>
static std::vector<PetscInt> ExtractConstrainedIndices(const ScalarMap& dofs)
{
  // Assumes IndexMap<Scalar> is iterable over (index,value)-like pairs.
  // If your IndexMap API differs, adapt only this function.
  std::vector<PetscInt> idx;
  idx.reserve(dofs.size());
  for (const auto& kv : dofs)
  {
    const auto& i = kv.first; // global dof index
    idx.push_back((PetscInt)i);
  }
  std::sort(idx.begin(), idx.end());
  idx.erase(std::unique(idx.begin(), idx.end()), idx.end());
  return idx;
}

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

  H1 uh(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension()); // P2 vector
  H1 ph(std::integral_constant<size_t, 1>{}, mesh);                          // P1 scalar

  std::cout << "Vector FES size: " << uh.getSize() << "\n";
  std::cout << "Scalar FES size: " << ph.getSize() << "\n";
  std::cout << "Assembling...\n";

  {
    // FE-exact MMS:
    // u quadratic, p linear, div(u)=0.
    // p = x + y + z  => grad p = (1,1,1)
    // u = [
    //   2x^2 + y^2 + z^2,
    //   2x^2 - 2xy,
    //   2x^2 - 2xz
    // ]
    // -Δu = (-8,-4,-4) so f = -Δu + ∇p = (-7,-3,-3)
    VectorFunction u_exact{
      2 * F::x * F::x + F::y * F::y + F::z * F::z,
      2 * F::x * F::x - 2 * F::x * F::y,
      2 * F::x * F::x - 2 * F::x * F::z
    };
    auto p_exact = F::x + F::y + F::z;
    VectorFunction f{ -7.0, -3.0, -3.0 };

    PETSc::Variational::TrialFunction u(uh); u.setName("u");
    PETSc::Variational::TrialFunction p(ph); p.setName("p");
    PETSc::Variational::TestFunction  v(uh);
    PETSc::Variational::TestFunction  q(ph);

    Problem stokes(u, p, v, q);

    // Keep a handle to the BC so we can access constrained DOFs.
    auto dbc = DirichletBC(u, u_exact); // optionally: .on(attr)
    stokes =
        Integral(Jacobian(u), Jacobian(v))
      - Integral(p, Div(v))
      + Integral(Div(u), q)
      - Integral(f, v)
      + dbc;

    auto t0 = std::chrono::high_resolution_clock::now();
    stokes.assemble();
    stokes.setFieldSplits(); // optional, but recommended for block
                             // preconditioners

    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Assembly time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
              << " ms\n";

    ::Mat& A = stokes.getLinearSystem().getOperator();
    ::Vec& x = stokes.getLinearSystem().getSolution();
    ::Vec& b = stokes.getLinearSystem().getVector();

    CheckPetsc(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
    CheckPetsc(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
    CheckPetsc(VecAssemblyBegin(b));
    CheckPetsc(VecAssemblyEnd(b));

    const PetscInt nu = (PetscInt)uh.getSize();
    const PetscInt np = (PetscInt)ph.getSize();

    // Assemble BC dofs and build mask for velocity constrained dofs.
    dbc.assemble();
    const auto& bc_dofs = dbc.getDOFs();
    std::vector<PetscInt> bc_u_idx = ExtractConstrainedIndices(bc_dofs);
    std::cout << "Dirichlet constrained dofs (u): " << bc_u_idx.size() << "\n";

    // Attach nullspace and remove it from RHS (recommended for singular operators).
    MatNullSpace nsp = AttachPressureNullspace(A, nu, np);
    CheckPetsc(MatNullSpaceRemove(nsp, b));

    // Solve (configured from command line via KSPSetFromOptions in your wrapper).
    Solver::KSP(stokes).solve();

    // Sanity
    PrintVecStats("x", x);
    PrintVecStats("b", b);
    if (VecHasNaNInf(x) || VecHasNaNInf(b))
      std::cout << "WARNING: NaN/Inf detected in x or b\n";

    // Solve residuals: full + masked
    const double rel_res_full = RelativeResidual(A, x, b);
    const double rel_res_free = RelativeResidualMasked(A, x, b, bc_u_idx);
    std::cout << "||A x - b||/||b|| (full)      = " << std::scientific << rel_res_full << "\n";
    std::cout << "||A x - b||/||b|| (FREE dofs) = " << std::scientific << rel_res_free << "\n";
    BlockResidualNormsMasked_Up(A, x, b, nu, np, bc_u_idx);

    // Assembly validation with FE-exact x_exact = [I_h u_exact, I_h p_exact]
    PETSc::Variational::GridFunction u_e(uh); u_e = u_exact;
    PETSc::Variational::GridFunction p_e(ph); p_e = p_exact;
    Vec xe = BuildExactVector_Up(x, nu, np, u_e.getData(), p_e.getData());

    PrintVecStats("x_exact", xe);

    const double rel_res_exact_full = RelativeResidual(A, xe, b);
    const double rel_res_exact_free = RelativeResidualMasked(A, xe, b, bc_u_idx);
    std::cout << "||A x_exact - b||/||b|| (full)      = " << std::scientific << rel_res_exact_full << "\n";
    std::cout << "||A x_exact - b||/||b|| (FREE dofs) = " << std::scientific << rel_res_exact_free << "\n";
    BlockResidualNormsMasked_Up(A, xe, b, nu, np, bc_u_idx);

    CheckPetsc(VecDestroy(&xe));
    CheckPetsc(MatNullSpaceDestroy(&nsp));

    mesh.save("Stokes.mesh");
    u.getSolution().save("Stokes_velocity.gf");
    p.getSolution().save("Stokes_pressure.gf");
  }

  PetscFinalize();
  return 0;
}
