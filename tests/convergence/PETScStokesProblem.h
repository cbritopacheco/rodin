/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef RODIN_TESTS_CONVERGENCE_PETSC_STOKES_PROBLEM_H
#define RODIN_TESTS_CONVERGENCE_PETSC_STOKES_PROBLEM_H

#include <functional>
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
   * operator. Solver status, an independently recomputed coefficient residual,
   * and the pressure mean are checked before physical-cell error integration.
   * MPI errors sum owned-cell contributions globally, excluding halo copies.
   */
  template <class MeshType>
  class PETScStokesProblem
  {
    public:
      PETScStokesProblem(
        const MeshType& mesh, const StokesData& data, size_t quadratureOrder = 16)
        : m_mesh(mesh),
          m_data(data),
          m_order(quadratureOrder)
      {}

      template <size_t K>
      StokesErrors solve(Real viscosity = 1) const
      {
        static_assert(K >= 2);
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
        PETSc::Solver::KSP solver(problem);
        solver.setType(KSPPREONLY);
        PC pc = nullptr;
        EXPECT_EQ(KSPGetPC(solver.getHandle(), &pc), PETSC_SUCCESS);
        EXPECT_EQ(PCSetType(pc, PCLU), PETSC_SUCCESS);
        EXPECT_EQ(PCFactorSetMatSolverType(pc, MATSOLVERMUMPS), PETSC_SUCCESS);
        solver.solve();
        KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
        EXPECT_EQ(KSPGetConvergedReason(solver.getHandle(), &reason), PETSC_SUCCESS);
        EXPECT_GT(reason, 0);
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
        return {ErrorNorm::computeVector(mesh, u.getSolution(), m_data.getVelocity(),
                  m_data.getVelocityJacobian(), m_order),
          ErrorNorm::compute(mesh, p.getSolution(), m_data.getPressure(),
            m_data.getPressureGradient(), m_order),
          ErrorNorm::computeDivergenceL2(mesh, u.getSolution(), m_order)};
      }

    private:
      std::reference_wrapper<const MeshType> m_mesh;
      StokesData m_data;
      size_t m_order;
  };
}

#endif
