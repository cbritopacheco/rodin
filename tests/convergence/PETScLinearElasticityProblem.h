/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef RODIN_TESTS_CONVERGENCE_PETSC_LINEARELASTICITYPROBLEM_H
#define RODIN_TESTS_CONVERGENCE_PETSC_LINEARELASTICITYPROBLEM_H

#include <functional>
#include <gtest/gtest.h>

#include "LinearElasticity.h"
#include "Rodin/PETSc.h"

namespace Rodin::Tests::Convergence
{
  /** @brief Shared vector linear-elasticity PETSc verification workload.
   * @par Mathematical formulation
   * @f$-\nabla\cdot\sigma(u)=f@f$ with
   * @f$\sigma=\lambda(\nabla\cdot u)I+2\mu\varepsilon(u)@f$,
   * @f$\lambda=1.5@f$, @f$\mu=0.5@f$ and full manufactured traces.
   * Removing the volumetric stiffness retains the original source and trace.
   * @par Architecture
   * Each solve creates a fresh fixed-layout vector space and system. Shared
   * physical data and vector norms are independent of assembly; a const
   * solve-scoped observer supports tensor and mapped-domain verification.
   */
  template <size_t K, class MeshType>
  class PETScLinearElasticityProblem
  {
    public:
      using Data = LinearElasticity::ManufacturedSolution;

      PETScLinearElasticityProblem(
        const MeshType& mesh, Data::Field field, size_t assemblyOrder = 16)
        : m_mesh(mesh),
          m_data(mesh.getDimension(), Lambda, Mu, field),
          m_order(assemblyOrder)
      {
        assert(mesh.getDimension() > 0);
        assert(mesh.getDimension() == mesh.getSpaceDimension());
      }

      ErrorNorms solve(bool omitVolumetric = false, Real tolerance = SolverTolerance,
        size_t normOrder = 18) const
      {
        return solve(
          omitVolumetric, tolerance, normOrder, [](const auto&, const auto&) {});
      }

      template <class Observer>
      ErrorNorms solve(
        bool omitVolumetric, Real tolerance, size_t normOrder, Observer&& observe) const
      {
        using namespace Variational;
        const auto& mesh = m_mesh.get();
        H1<K, Math::SpatialVector<Real>, MeshType> space(
          std::integral_constant<size_t, K>{}, mesh, mesh.getDimension());
        PETSc::Variational::TrialFunction u(space);
        PETSc::Variational::TestFunction v(space);
        auto volumetric =
          Integral((omitVolumetric ? Real(0) : m_data.getLambda()) * Div(u), Div(v));
        auto shear = Integral(m_data.getMu() * (Jacobian(u) + Jacobian(u).T()),
          Real(0.5) * (Jacobian(v) + Jacobian(v).T()));
        auto load = Integral(m_data.forcing, v);
        volumetric.setOrder(m_order);
        shear.setOrder(m_order);
        load.setOrder(m_order);
        Problem problem(u, v);
        problem = volumetric + shear - load + DirichletBC(u, m_data.exact);
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
        return ErrorNorm::computeVector(
          mesh, u.getSolution(), m_data.exact,
          [this](const Geometry::Point& p) { return m_data.jacobian(p); }, normOrder);
      }

    private:
      static constexpr Real Lambda = 1.5, Mu = 0.5;
      static constexpr Real SolverTolerance = 1e-13;
      static constexpr Real AbsoluteTolerance = SolverTolerance / 10;
      static constexpr Real ReportedResidualTolerance = 1e-8;
      static constexpr Real IndependentResidualTolerance = 1e-11;
      static constexpr Real DivergenceTolerance = 1e5;
      static constexpr size_t MaxIterations = 50000;
      std::reference_wrapper<const MeshType> m_mesh;
      Data m_data;
      size_t m_order;
  };
}

#endif
