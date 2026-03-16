/*
 * Stokes P2–P1 convergence validation with mesh refinement.
 *
 * Problem:
 *   -Δu + ∇p = f   in Ω = (0,1)^3
 *        div u = 0 in Ω
 *             u = u_exact on ∂Ω
 *
 * This script:
 *   1) solves the problem on a sequence of refined tetrahedral meshes,
 *   2) uses a smooth divergence-free manufactured solution that is NOT exactly
 *      representable by P2–P1,
 *   3) handles the pressure nullspace with PETSc MatNullSpace,
 *   4) reports:
 *        - h
 *        - numbers of velocity / pressure dofs
 *        - assembly / solve times
 *        - algebraic residual on free dofs
 *        - ||u_h - u_exact||_{L2}
 *        - |u_h - u_exact|_{H1}
 *        - ||p_h - p_exact - c_h||_{L2}
 *        - ||div u_h||_{L2}
 *        - observed convergence rates
 *
 * Expected asymptotic rates for Taylor–Hood P2–P1:
 *   velocity L2  ~ O(h^3)
 *   velocity H1  ~ O(h^2)
 *   pressure L2  ~ O(h^2)
 *
 * Examples:
 *   OMP_NUM_THREADS=8 ./PETSc_Seq_Validation_Stokes_Convergence \
 *     -ksp_type gmres -pc_type gamg -ksp_rtol 1e-12 -ksp_converged_reason
 *
 * Optional parameters:
 *   -cells0      initial number of cells per direction   (default: 2)
 *   -levels      number of refinement levels             (default: 5)
 *   -ref_factor  refinement factor per level             (default: 2)
 *   -save_last   save last-level mesh/solutions          (default: false)
 *
 * Important API note:
 *   All post-processing integrals are performed as Integral(gridfunction).compute().
 *   No Integral(expr) is used here.
 */

#include "Rodin/Geometry/Polytope.h"
#include "Rodin/PETSc/Variational/TestFunction.h"
#include <Rodin/PETSc.h>

#include <Rodin/Types.h>
#include <Rodin/Solver.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

#include <petscksp.h>
#include <petscmat.h>
#include <petscvec.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
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

static void ZeroEntries(Vec v, const std::vector<PetscInt>& idx)
{
  for (PetscInt i : idx)
    CheckPetsc(VecSetValue(v, i, 0.0, INSERT_VALUES));
  CheckPetsc(VecAssemblyBegin(v));
  CheckPetsc(VecAssemblyEnd(v));
}

static double RelativeResidualMasked(Mat A, Vec x, Vec b, const std::vector<PetscInt>& mask)
{
  Vec r = nullptr, bb = nullptr;
  CheckPetsc(VecDuplicate(b, &r));
  CheckPetsc(VecDuplicate(b, &bb));

  CheckPetsc(MatMult(A, x, r));
  CheckPetsc(VecAXPY(r, -1.0, b));
  CheckPetsc(VecCopy(b, bb));

  ZeroEntries(r, mask);
  ZeroEntries(bb, mask);

  PetscReal nr = 0, nb = 0;
  CheckPetsc(VecNorm(r, NORM_2, &nr));
  CheckPetsc(VecNorm(bb, NORM_2, &nb));

  CheckPetsc(VecDestroy(&r));
  CheckPetsc(VecDestroy(&bb));
  return (nb > 0) ? static_cast<double>(nr / nb) : static_cast<double>(nr);
}

static MatNullSpace AttachPressureNullspace(Mat A, PetscInt nu, PetscInt np)
{
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
  return nsp;
}

template <class ScalarMap>
static std::vector<PetscInt> ExtractConstrainedIndices(const ScalarMap& dofs)
{
  std::vector<PetscInt> idx;
  idx.reserve(dofs.size());
  for (const auto& kv : dofs)
    idx.push_back(static_cast<PetscInt>(kv.first));

  std::sort(idx.begin(), idx.end());
  idx.erase(std::unique(idx.begin(), idx.end()), idx.end());
  return idx;
}

static double Rate(double e_old, double e_new, double h_old, double h_new)
{
  if (!(e_old > 0.0) || !(e_new > 0.0) || !(h_old > 0.0) || !(h_new > 0.0))
    return std::numeric_limits<double>::quiet_NaN();
  return std::log(e_old / e_new) / std::log(h_old / h_new);
}

/* ---------- Post-processing helpers ----------
 *
 * Every integral is performed on an actual GridFunction:
 *   Integral(gridfunction).compute()
 */

template <class FES, class ScalarExpr>
static double IntegrateScalarGridFunction(const FES& scalarSpace, const ScalarExpr& expr)
{
  PETSc::Variational::TestFunction v(scalarSpace);
  LinearForm lf(v);
  lf = Integral(expr, v);
  lf.assemble();
  PETSc::Variational::GridFunction one(scalarSpace);
  one = 1.0;
  return lf(one);
}

template <class FES, class VecExpr>
static double ComputeVectorL2Error(const FES& scalarSpace, const VecExpr& err)
{
  PETSc::Variational::GridFunction integrand(scalarSpace);
  integrand = Dot(err, err);
  const double e2 = static_cast<double>(Integral(integrand).compute());
  return std::sqrt(std::max(0.0, e2));
}

template <class FES, class VecExpr>
static double ComputeVectorH1SemiError(const FES& scalarSpace, const VecExpr& err)
{
  PETSc::Variational::GridFunction integrand(scalarSpace);
  integrand = Dot(Jacobian(err), Jacobian(err));
  const double e2 = static_cast<double>(Integral(integrand).compute());
  return std::sqrt(std::max(0.0, e2));
}

template <class FES, class ScalarExpr>
static double ComputeScalarL2Error(const FES& scalarSpace, const ScalarExpr& err)
{
  PETSc::Variational::GridFunction integrand(scalarSpace);
  integrand = err * err;
  const double e2 = static_cast<double>(Integral(integrand).compute());
  return std::sqrt(std::max(0.0, e2));
}

template <class FES, class VecExpr>
static double ComputeDivergenceL2(const FES& scalarSpace, const VecExpr& u)
{
  PETSc::Variational::GridFunction integrand(scalarSpace);
  integrand = Div(u) * Div(u);
  const double e2 = static_cast<double>(Integral(integrand).compute());
  return std::sqrt(std::max(0.0, e2));
}

struct LevelResult
{
  int level = 0;
  size_t cells = 0;
  size_t points = 0;
  double h = 0.0;

  size_t nu = 0;
  size_t np = 0;

  long long assembly_ms = 0;
  long long solve_ms = 0;

  double rel_res_free = 0.0;
  double u_l2 = 0.0;
  double u_h1 = 0.0;
  double p_l2 = 0.0;
  double div_l2 = 0.0;
};

static void PrintSeparator(int n = 154)
{
  for (int i = 0; i < n; ++i)
    std::cout << '-';
  std::cout << '\n';
}

static void PrintSummaryTable(const std::vector<LevelResult>& rows)
{
  std::cout << "\n";
  PrintSeparator();
  std::cout << "Stokes P2-P1 convergence summary\n";
  PrintSeparator();

  std::cout
    << std::setw(5)  << "lev"
    << std::setw(8)  << "cells"
    << std::setw(8)  << "pts"
    << std::setw(12) << "h"
    << std::setw(12) << "nu"
    << std::setw(10) << "np"
    << std::setw(12) << "asm(ms)"
    << std::setw(12) << "sol(ms)"
    << std::setw(16) << "relres(free)"
    << std::setw(16) << "||u-ue||L2"
    << std::setw(10) << "rate"
    << std::setw(16) << "|u-ue|H1"
    << std::setw(10) << "rate"
    << std::setw(16) << "||p-pe||L2"
    << std::setw(10) << "rate"
    << std::setw(16) << "||div uh||L2"
    << '\n';

  PrintSeparator();

  for (size_t i = 0; i < rows.size(); ++i)
  {
    const auto& r = rows[i];

    auto fmt_rate = [](double x) -> std::string
    {
      if (std::isnan(x))
        return "   -   ";
      std::ostringstream os;
      os << std::fixed << std::setprecision(3) << x;
      return os.str();
    };

    const double ru  = (i == 0) ? std::numeric_limits<double>::quiet_NaN()
                                : Rate(rows[i - 1].u_l2, rows[i].u_l2, rows[i - 1].h, rows[i].h);
    const double rh1 = (i == 0) ? std::numeric_limits<double>::quiet_NaN()
                                : Rate(rows[i - 1].u_h1, rows[i].u_h1, rows[i - 1].h, rows[i].h);
    const double rp  = (i == 0) ? std::numeric_limits<double>::quiet_NaN()
                                : Rate(rows[i - 1].p_l2, rows[i].p_l2, rows[i - 1].h, rows[i].h);

    std::cout
      << std::setw(5)  << r.level
      << std::setw(8)  << r.cells
      << std::setw(8)  << r.points
      << std::setw(12) << std::scientific << std::setprecision(3) << r.h
      << std::setw(12) << std::defaultfloat << r.nu
      << std::setw(10) << r.np
      << std::setw(12) << r.assembly_ms
      << std::setw(12) << r.solve_ms
      << std::setw(16) << std::scientific << std::setprecision(3) << r.rel_res_free
      << std::setw(16) << std::scientific << std::setprecision(3) << r.u_l2
      << std::setw(10) << fmt_rate(ru)
      << std::setw(16) << std::scientific << std::setprecision(3) << r.u_h1
      << std::setw(10) << fmt_rate(rh1)
      << std::setw(16) << std::scientific << std::setprecision(3) << r.p_l2
      << std::setw(10) << fmt_rate(rp)
      << std::setw(16) << std::scientific << std::setprecision(3) << r.div_l2
      << '\n';
  }

  PrintSeparator();
  std::cout << "Expected asymptotic rates: velocity L2 ~ 3, velocity H1 ~ 2, pressure L2 ~ 2\n";
  PrintSeparator();
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  PetscInt cells0 = 2;
  PetscInt levels = 5;
  PetscInt ref_factor = 2;
  PetscBool save_last = PETSC_FALSE;

  CheckPetsc(PetscOptionsGetInt(nullptr, nullptr, "-cells0", &cells0, nullptr));
  CheckPetsc(PetscOptionsGetInt(nullptr, nullptr, "-levels", &levels, nullptr));
  CheckPetsc(PetscOptionsGetInt(nullptr, nullptr, "-ref_factor", &ref_factor, nullptr));
  CheckPetsc(PetscOptionsGetBool(nullptr, nullptr, "-save_last", &save_last, nullptr));

  if (cells0 < 1)     cells0 = 1;
  if (levels < 1)     levels = 1;
  if (ref_factor < 2) ref_factor = 2;

  std::vector<LevelResult> results;
  results.reserve(static_cast<size_t>(levels));

  for (PetscInt lev = 0; lev < levels; ++lev)
  {
    const size_t cells = static_cast<size_t>(cells0 * std::pow(ref_factor, lev));
    const size_t pts   = cells + 1;
    const double h     = 1.0 / static_cast<double>(cells);

    std::cout << "\n";
    PrintSeparator();
    std::cout << "Level " << lev
              << "  |  cells per direction = " << cells
              << "  |  points per direction = " << pts
              << "  |  h = " << std::scientific << std::setprecision(6) << h << "\n";
    PrintSeparator();

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, { pts, pts, pts });
    mesh.scale(1.0 / cells);

    mesh.getConnectivity().compute(2, 3);
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 uh(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    H1 ph(std::integral_constant<size_t, 1>{}, mesh);

    /*
     * Scalar post-processing space.
     * Chosen richer than P1 to reduce pollution when representing squared errors.
     */
    H1 diag(std::integral_constant<size_t, 4>{}, mesh);

    std::cout << "Velocity dofs: " << uh.getSize() << "\n";
    std::cout << "Pressure dofs: " << ph.getSize() << "\n";

    const auto x = F::x;
    const auto y = F::y;
    const auto z = F::z;

    const auto ax   = x * x * (1.0 - x) * (1.0 - x);
    const auto ay   = y * y * (1.0 - y) * (1.0 - y);

    const auto dax  = 2.0 * x - 6.0 * x * x + 4.0 * x * x * x;
    const auto day  = 2.0 * y - 6.0 * y * y + 4.0 * y * y * y;

    const auto d2ax = 2.0 - 12.0 * x + 12.0 * x * x;
    const auto d2ay = 2.0 - 12.0 * y + 12.0 * y * y;

    const auto d3ax = -12.0 + 24.0 * x;
    const auto d3ay = -12.0 + 24.0 * y;

    // Divergence-free manufactured velocity:
    //   u = ( a(x) b'(y), -a'(x) b(y), 0 )
    // with a(t)=t^2(1-t)^2, b(t)=t^2(1-t)^2.
    VectorFunction u_exact{
      ax * day,
      -dax * ay,
      0.0
    };

    // Non-P1 pressure so that pressure convergence is nontrivial.
    auto p_exact = x * x + y * y + z * z;

    // f = -Δu + ∇p
    VectorFunction f{
      -(d2ax * day + ax * d3ay) + 2.0 * x,
       (d3ax * ay + dax * d2ay) + 2.0 * y,
       2.0 * z
    };

    PETSc::Variational::TrialFunction u(uh); u.setName("u");
    PETSc::Variational::TrialFunction p(ph); p.setName("p");
    PETSc::Variational::TestFunction  v(uh);
    PETSc::Variational::TestFunction  q(ph);

    Problem stokes(u, p, v, q);
    auto dbc = DirichletBC(u, u_exact);

    stokes =
        Integral(Jacobian(u), Jacobian(v))
      - Integral(p, Div(v))
      + Integral(Div(u), q)
      - Integral(f, v)
      + dbc;

    const auto ta0 = std::chrono::high_resolution_clock::now();
    stokes.assemble();
    stokes.setFieldSplits();
    const auto ta1 = std::chrono::high_resolution_clock::now();

    ::Mat& A = stokes.getLinearSystem().getOperator();
    ::Vec& xh = stokes.getLinearSystem().getSolution();
    ::Vec& b  = stokes.getLinearSystem().getVector();

    CheckPetsc(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
    CheckPetsc(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
    CheckPetsc(VecAssemblyBegin(b));
    CheckPetsc(VecAssemblyEnd(b));

    const PetscInt nu = static_cast<PetscInt>(uh.getSize());
    const PetscInt np = static_cast<PetscInt>(ph.getSize());

    dbc.assemble();
    const auto& bc_dofs = dbc.getDOFs();
    const std::vector<PetscInt> bc_u_idx = ExtractConstrainedIndices(bc_dofs);

    MatNullSpace nsp = AttachPressureNullspace(A, nu, np);
    CheckPetsc(MatNullSpaceRemove(nsp, b));

    const auto ts0 = std::chrono::high_resolution_clock::now();
    Solver::KSP(stokes).solve();
    const auto ts1 = std::chrono::high_resolution_clock::now();

    const long long asm_ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(ta1 - ta0).count();
    const long long sol_ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(ts1 - ts0).count();

    const double rel_res_free = RelativeResidualMasked(A, xh, b, bc_u_idx);

    auto& uh_sol = u.getSolution();
    auto& ph_sol = p.getSolution();

    /*
     * Pressure gauge correction:
     *   c_h = mean(p_h - p_exact) over Ω.
     * Since Ω = (0,1)^3, |Ω| = 1, so mean = integral.
     */
    const double p_mean_shift = IntegrateScalarGridFunction(diag, ph_sol - p_exact);

    const auto err_u = uh_sol - u_exact;
    const auto err_p = ph_sol - p_exact - p_mean_shift;

    GridFunction err_u_gf(uh); err_u_gf = err_u;
    const double u_l2   = ComputeVectorL2Error(diag, err_u);
    const double u_h1   = ComputeVectorH1SemiError(diag, err_u_gf);
    const double p_l2   = ComputeScalarL2Error(diag, err_p);
    const double div_l2 = ComputeDivergenceL2(diag, uh_sol);

    std::cout << "Assembly time [ms]       : " << asm_ms << "\n";
    std::cout << "Solve time [ms]          : " << sol_ms << "\n";
    std::cout << "Relative residual (free) : " << std::scientific << std::setprecision(6) << rel_res_free << "\n";
    std::cout << "||u_h - u_exact||_L2     : " << std::scientific << std::setprecision(6) << u_l2 << "\n";
    std::cout << "|u_h - u_exact|_H1       : " << std::scientific << std::setprecision(6) << u_h1 << "\n";
    std::cout << "||p_h - p_exact - c||_L2 : " << std::scientific << std::setprecision(6) << p_l2 << "\n";
    std::cout << "||div u_h||_L2           : " << std::scientific << std::setprecision(6) << div_l2 << "\n";

    results.push_back(LevelResult{
      static_cast<int>(lev),
      cells,
      pts,
      h,
      static_cast<size_t>(uh.getSize()),
      static_cast<size_t>(ph.getSize()),
      asm_ms,
      sol_ms,
      rel_res_free,
      u_l2,
      u_h1,
      p_l2,
      div_l2
    });

    if (save_last && lev == levels - 1)
    {
      mesh.save("StokesConvergence.mesh", IO::FileFormat::MFEM);
      uh_sol.save("StokesConvergence_velocity.gf", IO::FileFormat::MFEM);
      ph_sol.save("StokesConvergence_pressure.gf", IO::FileFormat::MFEM);

      PETSc::Variational::GridFunction u_e(uh); u_e = u_exact;
      PETSc::Variational::GridFunction p_e(ph); p_e = p_exact;
      u_e.save("StokesConvergence_velocity_exact.gf", IO::FileFormat::MFEM);
      p_e.save("StokesConvergence_pressure_exact.gf", IO::FileFormat::MFEM);
    }

    CheckPetsc(MatNullSpaceDestroy(&nsp));
  }

  PrintSummaryTable(results);

  PetscFinalize();
  return 0;
}
