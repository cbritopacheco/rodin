/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef RODIN_TESTS_CONVERGENCE_PETSC_DIFFUSIONPROBLEM_H
#define RODIN_TESTS_CONVERGENCE_PETSC_DIFFUSIONPROBLEM_H

#include <functional>
#include <gtest/gtest.h>

#include "Conductivity.h"
#include "Rodin/PETSc.h"

namespace Rodin::Tests::Convergence
{
  /** @brief Shared real-PETSc scalar diffusion verification workload.
   * @par Mathematical formulation
   * @f$-\nabla\cdot(\gamma\nabla u)=f@f$, with @f$\gamma=1@f$ for
   * Poisson or @f$\gamma=1+\sum_jx_j@f$ for variable conductivity.
   * Sources and Dirichlet traces retain the selected correct coefficient
   * when a deliberately wrong operator is requested.
   * @par Architecture
   * Each measurement owns a fresh fixed-layout space and linear system.
   * Independent norm quadrature and a solve-scoped const observer support
   * flat and mapped domains without duplicating field evaluation or MPI norms.
   */
  template <size_t K, class MeshType>
  class PETScDiffusionProblem
  {
    public:
      PETScDiffusionProblem(const MeshType& mesh, bool poisson,
        ConductivityData::Field field, size_t assemblyOrder = 16)
        : m_mesh(mesh),
          m_data(mesh.getDimension(), field),
          m_poisson(poisson),
          m_order(assemblyOrder)
      {
        assert(mesh.getDimension() > 0);
        assert(mesh.getDimension() == mesh.getSpaceDimension());
      }

      ErrorNorms solve(
        bool wrong = false, Real tolerance = SolverTolerance, size_t normOrder = 18) const
      {
        return solve(wrong, tolerance, normOrder, [](const auto&, const auto&) {});
      }

      template <class Observer>
      ErrorNorms solve(
        bool wrong, Real tolerance, size_t normOrder, Observer&& observe) const
      {
        using namespace Variational;
        const auto& mesh = m_mesh.get();
        const auto exact = m_data.getSolution();
        const auto coefficient = m_data.getCoefficient(m_poisson || wrong);
        const Real scale = m_poisson && wrong ? Real(2) : Real(1);
        H1<K, Real, MeshType> space(std::integral_constant<size_t, K>{}, mesh);
        PETSc::Variational::TrialFunction u(space);
        PETSc::Variational::TestFunction v(space);
        auto stiffness = Integral(scale * coefficient * Grad(u), Grad(v));
        auto load = Integral(m_data.getSource(m_poisson), v);
        stiffness.setOrder(m_order);
        load.setOrder(m_order);
        Problem problem(u, v);
        problem = stiffness - load + DirichletBC(u, exact);
        PETSc::Solver::CG solver(problem);
        solver.setTolerances(
          tolerance, AbsoluteTolerance, DivergenceTolerance, MaxIterations);
        solver.solve();
        KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
        EXPECT_EQ(KSPGetConvergedReason(solver.getHandle(), &reason), PETSC_SUCCESS);
        EXPECT_GT(reason, 0);
        EXPECT_TRUE(std::isfinite(solver.getError()));
        EXPECT_LT(solver.getError(), ReportedResidualTolerance);

        const auto& system = problem.getLinearSystem();
        Vec residual = nullptr;
        EXPECT_EQ(VecDuplicate(system.getVector(), &residual), PETSC_SUCCESS);
        EXPECT_EQ(
          MatMult(system.getOperator(), system.getSolution(), residual), PETSC_SUCCESS);
        EXPECT_EQ(VecAXPY(residual, -1, system.getVector()), PETSC_SUCCESS);
        PetscReal norm = 0, rhsNorm = 0;
        EXPECT_EQ(VecNorm(residual, NORM_2, &norm), PETSC_SUCCESS);
        EXPECT_EQ(VecNorm(system.getVector(), NORM_2, &rhsNorm), PETSC_SUCCESS);
        const Real relative = norm / std::max(Real(1), rhsNorm);
        EXPECT_EQ(VecDestroy(&residual), PETSC_SUCCESS);
        EXPECT_TRUE(std::isfinite(relative));
        EXPECT_LT(relative, IndependentResidualTolerance);
        observe(u.getSolution(), m_data);
        return ErrorNorm::compute(
          mesh, u.getSolution(), exact, m_data.getGradient(), normOrder);
      }

    private:
      static constexpr Real SolverTolerance = 1e-13;
      static constexpr Real AbsoluteTolerance = SolverTolerance / 10;
      static constexpr Real ReportedResidualTolerance = 1e-8;
      static constexpr Real IndependentResidualTolerance = 1e-11;
      static constexpr Real DivergenceTolerance = 1e5;
      static constexpr size_t MaxIterations = 50000;
      std::reference_wrapper<const MeshType> m_mesh;
      ConductivityData m_data;
      bool m_poisson;
      size_t m_order;
  };
}

#endif
