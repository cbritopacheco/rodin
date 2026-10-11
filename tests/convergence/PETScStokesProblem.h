/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef RODIN_TESTS_CONVERGENCE_PETSC_STOKES_PROBLEM_H
#define RODIN_TESTS_CONVERGENCE_PETSC_STOKES_PROBLEM_H

#include <functional>
#include <utility>
#include <gtest/gtest.h>

#include "Stokes.h"
#include "Rodin/Assembly.h"
#include "Rodin/PETSc.h"

#ifdef RODIN_USE_MPI
#include "MPIConvergence.h"
#include "Rodin/MPI/Variational/H1/H1.h"
#include "Rodin/MPI/Variational/P0g.h"
#endif

namespace Rodin::Tests::Convergence
{
  /** @brief PETSc mixed Stokes workload shared by h, p, and hp studies.
   * @par Architecture
   * One mesh and continuous manufactured data define the workload. Each
   * solve creates fresh velocity degree @f$K@f$, pressure degree @f$K-1@f$,
   * and global mean-multiplier spaces, along with a fresh PETSc system.
   * Local and MPI meshes consume the same variational form without DOF-layout
   * assumptions or matrix resizing. PREONLY/LU/MUMPS handles the saddle-point
   * operator, with a 100-percent factor-workspace margin for delayed pivots.
   * An explicit positive refinement count instead selects unit-damped
   * Richardson iterations with the retained LU factors: one initial solve
   * and the requested residual corrections. This global linear-system
   * operation retains distributed vectors; MUMPS internal refinement is not
   * used because its distributed-solution path disables it. The actual
   * iteration count is checked independently of the coefficient residual.
   * MUMPS factorization status, solver status, an independent coefficient residual,
   * and the pressure mean are checked before physical-cell error integration.
   * MPI errors sum owned-cell contributions globally, excluding halo copies.
   */
  template <class MeshType>
  class PETScStokesProblem
  {
    public:
      PETScStokesProblem(const MeshType& mesh, const StokesData& data,
        size_t quadratureOrder = 16, PetscInt refinementSteps = 0)
        : m_mesh(mesh),
          m_data(data),
          m_order(quadratureOrder),
          m_refinementSteps(refinementSteps)
      {
        assert(refinementSteps >= 0 && refinementSteps < PETSC_MAX_INT);
      }

      template <size_t K>
      StokesErrors solve(Real viscosity = 1) const
      {
        return solve<K>(viscosity, 0, [](const auto&, const auto&, const auto&) {});
      }

      /** @brief Observe live converged fields with independent norm quadrature. */
      template <size_t K, class Observer>
      StokesErrors solve(Real viscosity, size_t normOrder, Observer&& observe) const
      {
        static_assert(K >= 2);
        SCOPED_TRACE(::testing::Message()
          << "velocity degree=" << K << " quadrature order=" << m_order);
        using namespace Variational;
        const auto& mesh = m_mesh.get();
        H1<K, Math::SpatialVector<Real>, MeshType> velocitySpace(
          std::integral_constant<size_t, K>{}, mesh, mesh.getDimension());
        H1<K - 1, Real, MeshType> pressureSpace(
          std::integral_constant<size_t, K - 1>{}, mesh);
        P0g<Real, MeshType> meanSpace(mesh);
        PETSc::Variational::TrialFunction u(velocitySpace);
        PETSc::Variational::TrialFunction p(pressureSpace);
        PETSc::Variational::TrialFunction lambda(meanSpace);
        PETSc::Variational::TestFunction v(velocitySpace);
        PETSc::Variational::TestFunction q(pressureSpace);
        PETSc::Variational::TestFunction mu(meanSpace);
        auto diffusion = Integral(viscosity * Jacobian(u), Jacobian(v));
        auto pressureVelocity = Integral(p, Div(v));
        auto incompressibility = Integral(Div(u), q);
        auto gaugePressure = Integral(lambda, q);
        auto gaugeMean = Integral(p, mu);
        auto body = Integral(m_data.getForcing(), v);
        diffusion.setOrder(m_order);
        pressureVelocity.setOrder(m_order);
        incompressibility.setOrder(m_order);
        gaugePressure.setOrder(m_order);
        gaugeMean.setOrder(m_order);
        body.setOrder(m_order);
        Problem problem(u, p, lambda, v, q, mu);
        problem = diffusion - pressureVelocity + incompressibility + gaugePressure +
          gaugeMean - body + DirichletBC(u, m_data.getVelocity());
        // Delayed pivots in the indefinite system can exceed MUMPS's default
        // factor-workspace estimate. Reserve space without changing the matrix.
        // Explicit runtime options retain precedence over these test defaults.
        for (const auto& [name, value] : {std::pair{"-mat_mumps_icntl_14", "100"},
               std::pair{"-mat_mumps_icntl_20", "0"}})
        {
          PetscBool set = PETSC_FALSE;
          EXPECT_EQ(PetscOptionsHasName(nullptr, nullptr, name, &set), PETSC_SUCCESS);
          if (!set)
          {
            EXPECT_EQ(PetscOptionsSetValue(nullptr, name, value), PETSC_SUCCESS);
          }
        }
        PETSc::Solver::KSP solver(problem);
        if (m_refinementSteps == 0)
          solver.setType(KSPPREONLY);
        else
        {
          solver.setType(KSPRICHARDSON);
          solver.setTolerances(PETSc::Solver::KSP::DEFAULT_RTOL,
            PETSc::Solver::KSP::DEFAULT_ABSTOL, PETSc::Solver::KSP::DEFAULT_DTOL,
            m_refinementSteps + 1);
          EXPECT_EQ(KSPSetNormType(solver.getHandle(), KSP_NORM_NONE), PETSC_SUCCESS);
          EXPECT_EQ(KSPRichardsonSetScale(solver.getHandle(), 1), PETSC_SUCCESS);
        }
        PC pc = nullptr;
        EXPECT_EQ(KSPGetPC(solver.getHandle(), &pc), PETSC_SUCCESS);
        EXPECT_EQ(PCSetType(pc, PCLU), PETSC_SUCCESS);
        EXPECT_EQ(PCFactorSetMatSolverType(pc, MATSOLVERMUMPS), PETSC_SUCCESS);
        solver.solve();
        PCFailedReason failure = PC_NOERROR;
        EXPECT_EQ(PCGetFailedReason(pc, &failure), PETSC_SUCCESS);
        Mat factor = nullptr;
        EXPECT_EQ(PCFactorGetMatrix(pc, &factor), PETSC_SUCCESS);
        PetscInt info = 0, detail = 0;
        EXPECT_EQ(MatMumpsGetInfog(factor, 1, &info), PETSC_SUCCESS);
        EXPECT_EQ(MatMumpsGetInfog(factor, 2, &detail), PETSC_SUCCESS);
        EXPECT_EQ(failure, PC_NOERROR)
          << "MUMPS INFOG(1)=" << info << " INFOG(2)=" << detail;
        EXPECT_EQ(info, 0) << "MUMPS INFOG(2)=" << detail;
        KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
        EXPECT_EQ(KSPGetConvergedReason(solver.getHandle(), &reason), PETSC_SUCCESS);
        EXPECT_GT(reason, 0);
        if (m_refinementSteps > 0)
        {
          const char* type = nullptr;
          EXPECT_EQ(KSPGetType(solver.getHandle(), &type), PETSC_SUCCESS);
          EXPECT_STREQ(type, KSPRICHARDSON);
          KSPNormType normType = KSP_NORM_DEFAULT;
          EXPECT_EQ(KSPGetNormType(solver.getHandle(), &normType), PETSC_SUCCESS);
          EXPECT_EQ(normType, KSP_NORM_NONE);
          EXPECT_EQ(
            solver.getIterationNumber(), static_cast<size_t>(m_refinementSteps + 1));
        }
        auto& system = problem.getLinearSystem();
        Vec residual = nullptr;
        EXPECT_EQ(VecDuplicate(system.getVector(), &residual), PETSC_SUCCESS);
        EXPECT_EQ(
          MatMult(system.getOperator(), system.getSolution(), residual), PETSC_SUCCESS);
        EXPECT_EQ(VecAXPY(residual, -1, system.getVector()), PETSC_SUCCESS);
        PetscReal norm = 0, rhsNorm = 0;
        EXPECT_EQ(VecNorm(residual, NORM_2, &norm), PETSC_SUCCESS);
        EXPECT_EQ(VecNorm(system.getVector(), NORM_2, &rhsNorm), PETSC_SUCCESS);
        EXPECT_TRUE(std::isfinite(norm));
        EXPECT_TRUE(std::isfinite(rhsNorm));
        EXPECT_LT(norm / std::max(Real(1), rhsNorm), 1e-11);
        EXPECT_EQ(VecDestroy(&residual), PETSC_SUCCESS);
        auto pressureMean = Integral(p.getSolution());
        pressureMean.setOrder(m_order + 2);
        EXPECT_LT(std::abs(pressureMean.compute()), 1e-10);
        const auto& velocity = u.getSolution();
        const auto& pressure = p.getSolution();
        observe(velocity, pressure, m_data);
        const size_t integrationOrder = normOrder == 0 ? m_order : normOrder;
        return {ErrorNorm::computeVector(mesh, velocity, m_data.getVelocity(),
                  m_data.getVelocityJacobian(), integrationOrder),
          ErrorNorm::compute(mesh, pressure, m_data.getPressure(),
            m_data.getPressureGradient(), integrationOrder),
          ErrorNorm::computeDivergenceL2(mesh, velocity, integrationOrder)};
      }

    private:
      std::reference_wrapper<const MeshType> m_mesh;
      StokesData m_data;
      size_t m_order;
      PetscInt m_refinementSteps;
  };
}

#endif
