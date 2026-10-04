/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef RODIN_TESTS_CONVERGENCE_PETSC_HELMHOLTZPROBLEM_H
#define RODIN_TESTS_CONVERGENCE_PETSC_HELMHOLTZPROBLEM_H

#include <functional>
#include <gtest/gtest.h>

#include "Helmholtz.h"
#include "Rodin/PETSc.h"

namespace Rodin::Tests::Convergence
{
  /** @brief Shared complex-PETSc Helmholtz verification workload.
   * @par Mathematical formulation
   * @f$-\Delta u-\tfrac14u=f@f$ with manufactured complex Dirichlet data.
   * Sources and traces retain the original mass term when the solved
   * operator deliberately omits it. The full-Dirichlet unit-box problem
   * is coercive since @f$1/4<\lambda_1=d\pi^2@f$.
   * @par Architecture
   * Each measurement owns a fresh fixed-layout space and linear system.
   * Independent norm quadrature and a solve-scoped const observer support
   * flat and mapped domains without duplicating field evaluation or MPI norms.
   */
  template <size_t K, class MeshType>
  class PETScHelmholtzProblem
  {
    public:
      PETScHelmholtzProblem(
        const MeshType& mesh, HelmholtzData::Field field, size_t assemblyOrder = 16)
        : m_mesh(mesh),
          m_data(mesh.getDimension(), field),
          m_order(assemblyOrder)
      {
        assert(mesh.getDimension() > 0);
        assert(mesh.getDimension() == mesh.getSpaceDimension());
      }

      ErrorNorms solve(bool omitMass = false, Real tolerance = SolverTolerance,
        size_t normOrder = 18) const
      {
        return solve(omitMass, tolerance, normOrder, [](const auto&, const auto&) {});
      }

      template <class Observer>
      ErrorNorms solve(
        bool omitMass, Real tolerance, size_t normOrder, Observer&& observe) const
      {
        using namespace Variational;
        const auto& mesh = m_mesh.get();
        const auto exact = m_data.getSolution();
        H1<K, Complex, MeshType> space(std::integral_constant<size_t, K>{}, mesh);
        PETSc::Variational::TrialFunction u(space);
        PETSc::Variational::TestFunction v(space);
        auto stiffness = Integral(Grad(u), Grad(v));
        auto mass =
          Integral((omitMass ? Real(0) : HelmholtzData::WaveNumberSquared) * u, v);
        auto load = Integral(m_data.getSource(), v);
        stiffness.setOrder(m_order);
        load.setOrder(m_order);
        mass.setOrder(m_order);
        Problem problem(u, v);
        problem = stiffness - mass - load + DirichletBC(u, exact);
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
      static constexpr Real ReportedResidualTolerance = 1e-9;
      static constexpr Real IndependentResidualTolerance = 1e-11;
      static constexpr Real DivergenceTolerance = 1e5;
      static constexpr size_t MaxIterations = 50000;
      std::reference_wrapper<const MeshType> m_mesh;
      HelmholtzData m_data;
      size_t m_order;
  };
}

#endif
