/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef RODIN_TESTS_CONVERGENCE_PETSC_REACTIONDIFFUSIONPROBLEM_H
#define RODIN_TESTS_CONVERGENCE_PETSC_REACTIONDIFFUSIONPROBLEM_H

#include <array>
#include <cmath>
#include <functional>
#include <gtest/gtest.h>

#include "ReactionDiffusion.h"
#include "Rodin/PETSc.h"

namespace Rodin::Tests::Convergence
{
  /** @brief Coupled PETSc verification workload with shared manufactured data.
   * @par Mathematical formulation
   * @f$-\kappa_i\Delta u_i+u_i+\alpha u_{1-i}=f_i@f$ with
   * @f$\kappa=(1,2)@f$ and @f$\alpha=0.2@f$, on a full-dimensional domain.
   * Both fields have their manufactured Dirichlet trace by default. Mixed
   * data prescribe this trace on one coordinate face and flux or Robin data
   * elsewhere. Pure Neumann data need no mean constraint because the reaction
   * matrix has eigenvalues 0.8 and 1.2. Omitting coupling
   * changes only the solved operator; sources retain the original reaction.
   * @par Architecture
   * Each solve constructs a fresh space and two-field fixed-layout system.
   * The common data supplies physical sources and analytic gradients. After
   * CG convergence, an independent coefficient residual is evaluated; const
   * solution views are available only during a solve-scoped observer. Field
   * norms use independent quadrature and owned-cell MPI reduction.
   */
  template <size_t K, class MeshType>
  class PETScReactionDiffusionProblem
  {
    public:
      enum class Boundary
      {
        Dirichlet,
        MixedNeumann,
        Robin,
        PureNeumann
      };
      static constexpr Geometry::Attribute DirichletAttribute = 501;
      static constexpr Geometry::Attribute NaturalAttribute = 502;

      explicit PETScReactionDiffusionProblem(const MeshType& mesh,
        ReactionDiffusionData::Field field = ReactionDiffusionData::Field::Smooth,
        size_t assemblyOrder = 16, Boundary boundary = Boundary::Dirichlet)
        : m_mesh(mesh),
          m_data(mesh.getDimension(), field),
          m_order(assemblyOrder),
          m_boundary(boundary)
      {
        assert(mesh.getDimension() > 0);
        assert(mesh.getDimension() == mesh.getSpaceDimension());
      }

      /** @brief Assemble, solve, and measure the global two-field problem.
       * @note All ranks of an MPI mesh communicator must participate in
       * distributed assembly, solving, residuals, and error norms. Pointwise
       * manufactured field and boundary evaluations are noncollective.
       */
      std::array<ErrorNorms, 2> solve(bool omitCoupling = false,
        Real tolerance = SolverTolerance, size_t normOrder = 0,
        bool omitFlux = false) const
      {
        return solve(
          omitCoupling, tolerance, normOrder,
          [](const auto&, const auto&, const auto&) {}, omitFlux);
      }

      /** @brief Solve with a rank-local solution observer before global norms.
       * @note The same communicator participation as ordinary solve applies.
       */
      template <class Observer>
      std::array<ErrorNorms, 2> solve(bool omitCoupling, Real tolerance, size_t normOrder,
        Observer&& observe, bool omitFlux = false) const
      {
        using namespace Variational;
        const auto& mesh = m_mesh.get();
        const auto exactU = m_data.getSolution(0), exactW = m_data.getSolution(1);
        H1<K, Real, MeshType> space(std::integral_constant<size_t, K>{}, mesh);
        PETSc::Variational::TrialFunction u(space), w(space);
        PETSc::Variational::TestFunction v(space), z(space);
        const Real alpha = omitCoupling ? Real(0) : ReactionDiffusionData::Coupling;
        auto aUU = Integral(Grad(u), Grad(v));
        auto aWW = Integral(Real(2) * Grad(w), Grad(z));
        auto rUU = Integral(u, v), rWW = Integral(w, z);
        auto rUW = Integral(alpha * w, v), rWU = Integral(alpha * u, z);
        auto bU = Integral(m_data.getSource(0), v), bW = Integral(m_data.getSource(1), z);
        aUU.setOrder(m_order);
        aWW.setOrder(m_order);
        rUU.setOrder(m_order);
        rWW.setOrder(m_order);
        rUW.setOrder(m_order);
        rWU.setOrder(m_order);
        bU.setOrder(m_order);
        bW.setOrder(m_order);
        Problem problem(u, w, v, z);
        if (m_boundary == Boundary::Dirichlet)
          problem = aUU + rUU + rUW - bU + aWW + rWW + rWU - bW + DirichletBC(u, exactU) +
            DirichletBC(w, exactW);
        else
        {
          const BoundaryNormal normal(mesh);
          const size_t dim = mesh.getDimension();
          const auto flux = [&](size_t component) {
            return RealFunction(
              [gradient = m_data.getGradient(component), normal, dim,
                coefficient = Real(component + 1)](const Geometry::Point& point) {
                const auto derivative = gradient(point);
                if (dim == 1)
                  return coefficient * (point(0) < Real(0.5) ? Real(-1) : Real(1)) *
                    derivative(0);
                const auto outward = normal(point);
                Real value = 0;
                for (size_t j = 0; j < dim; ++j)
                  value += coefficient * derivative(j) * outward(j);
                return value;
              });
          };
          const Real beta = m_boundary == Boundary::Robin ? Real(1) : Real(0);
          const Real includeFlux = omitFlux ? Real(0) : Real(1);
          auto robinU = BoundaryIntegral(beta * u, v);
          auto robinW = BoundaryIntegral(beta * w, z);
          auto fluxU = BoundaryIntegral(includeFlux * flux(0) + beta * exactU, v);
          auto fluxW = BoundaryIntegral(includeFlux * flux(1) + beta * exactW, z);
          robinU.setOrder(m_order);
          robinW.setOrder(m_order);
          fluxU.setOrder(m_order);
          fluxW.setOrder(m_order);
          if (m_boundary == Boundary::PureNeumann)
            problem = aUU + rUU + rUW - bU + aWW + rWW + rWU - bW - fluxU - fluxW;
          else
            problem = aUU + rUU + rUW - bU + aWW + rWW + rWU - bW +
              robinU.over(NaturalAttribute) + robinW.over(NaturalAttribute) -
              fluxU.over(NaturalAttribute) - fluxW.over(NaturalAttribute) +
              DirichletBC(u, exactU).on(DirichletAttribute) +
              DirichletBC(w, exactW).on(DirichletAttribute);
        }
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

        const auto& solutionU = u.getSolution();
        const auto& solutionW = w.getSolution();
        observe(solutionU, solutionW, m_data);
        const size_t integrationOrder = normOrder == 0 ? m_order : normOrder;
        return {ErrorNorm::compute(
                  mesh, solutionU, exactU, m_data.getGradient(0), integrationOrder),
          ErrorNorm::compute(
            mesh, solutionW, exactW, m_data.getGradient(1), integrationOrder)};
      }

    private:
      // Relative stopping tolerance and dimensionless independent residual budgets.
      static constexpr Real SolverTolerance = 1e-13;
      static constexpr Real AbsoluteTolerance = SolverTolerance / 10;
      static constexpr Real ReportedResidualTolerance = 1e-8;
      static constexpr Real IndependentResidualTolerance = 1e-11;
      static constexpr Real DivergenceTolerance =
        1e5; // Remote divergence safety threshold.
      static constexpr size_t MaxIterations =
        50000; // Safety cap, not an acceptance condition.
      std::reference_wrapper<const MeshType> m_mesh;
      ReactionDiffusionData m_data;
      size_t m_order;
      Boundary m_boundary;
  };
}

#endif
