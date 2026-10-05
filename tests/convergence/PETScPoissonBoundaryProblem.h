/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef RODIN_TESTS_CONVERGENCE_PETSC_POISSONBOUNDARYPROBLEM_H
#define RODIN_TESTS_CONVERGENCE_PETSC_POISSONBOUNDARYPROBLEM_H

#include <gtest/gtest.h>

#include "Conductivity.h"
#include "Rodin/PETSc.h"

namespace Rodin::Tests::Convergence
{
  /** @brief Poisson with manufactured essential, natural and Robin data.
   * @par Mathematical formulation
   * @f$-\Delta u=f@f$, @f$u=e^{\sum_jx_j}@f$, with
   * @f$g=\nabla u\cdot n@f$ or @f$r=g+2u@f$ on the natural partition.
   * Pure Neumann data use @f$u=e^{\sum_jx_j}-(e-1)^d@f$ and a P0g
   * multiplier imposing @f$\int_\Omega u=0@f$.
   * @par Architecture
   * Each solve owns its fixed-layout space and system. Boundary attributes
   * are prescribed before partitioning; no entity is matched by coordinates.
   * ErrorNorm integrates the field independently. Solver status, coefficient
   * residual, mean and compatibility multiplier have separate acceptance.
   */
  template <size_t K, class MeshType>
  class PETScPoissonBoundaryProblem
  {
    public:
      enum class Condition
      {
        Neumann,
        Robin,
        PureNeumann
      };
      static constexpr Geometry::Attribute DirichletAttribute = 101;
      static constexpr Geometry::Attribute NaturalAttribute = 102;

      PETScPoissonBoundaryProblem(
        const MeshType& mesh, Condition condition, size_t assemblyOrder = 16)
        : m_mesh(mesh),
          m_condition(condition),
          m_order(assemblyOrder)
      {}

      ErrorNorms solve(bool wrongFlux = false, size_t normOrder = 18,
        Real tolerance = SolverTolerance) const
      {
        using namespace Variational;
        const auto& mesh = m_mesh.get();
        const size_t dim = mesh.getDimension();
        const bool pure = m_condition == Condition::PureNeumann;
        const ConductivityData data(dim, ConductivityData::Field::Exponential);
        const auto exponential = data.getSolution();
        const Real mean = pure ? std::pow(std::expm1(Real(1)), Real(dim)) : 0;
        const RealFunction exact([exponential, mean](const Geometry::Point& p) {
          return exponential(p) - mean;
        });
        const auto gradient = data.getGradient();
        const BoundaryNormal normal(mesh);
        const RealFunction flux([dim, gradient, normal](const Geometry::Point& p) {
          const auto derivative = gradient(p);
          // Unit-box endpoint classification, not a correspondence between meshes.
          if (dim == 1)
            return (p(0) < Real(0.5) ? -1 : 1) * derivative(0);
          return derivative.dot(normal(p));
        });
        const Real scale = wrongFlux ? 0 : 1;
        H1<K, Real, MeshType> space(std::integral_constant<size_t, K>{}, mesh);
        PETSc::Variational::TrialFunction u(space);
        PETSc::Variational::TestFunction v(space);
        auto stiffness = Integral(Grad(u), Grad(v));
        auto load = Integral(data.getSource(true), v);
        auto natural = BoundaryIntegral(scale * flux, v);
        stiffness.setOrder(m_order);
        load.setOrder(m_order);
        natural.setOrder(m_order);
        const auto measure = [&] {
          return ErrorNorm::compute(mesh, u.getSolution(), exact, gradient, normOrder);
        };
        if (!pure)
        {
          Problem problem(u, v);
          if (m_condition == Condition::Neumann)
            problem = stiffness - load - natural.over(NaturalAttribute) +
              DirichletBC(u, exact).on(DirichletAttribute);
          else
          {
            auto mass = BoundaryIntegral(u, v);
            auto robin = BoundaryIntegral(scale * flux + RobinCoefficient * exact, v);
            mass.setOrder(m_order);
            robin.setOrder(m_order);
            problem = stiffness + RobinCoefficient * mass.over(NaturalAttribute) - load -
              robin.over(NaturalAttribute) + DirichletBC(u, exact).on(DirichletAttribute);
          }
          PETSc::Solver::CG solver(problem);
          solver.setTolerances(
            tolerance, tolerance / 10, DivergenceTolerance, MaxIterations);
          solver.solve();
          checkSystem(problem, solver);
          return measure();
        }
#ifdef RODIN_POISSON_BOUNDARY_MUMPS
        P0g<Real, MeshType> constants(mesh);
        PETSc::Variational::TrialFunction lambda(constants);
        PETSc::Variational::TestFunction mu(constants);
        auto gauge = Integral(lambda, v);
        auto constraint = Integral(u, mu);
        gauge.setOrder(m_order);
        constraint.setOrder(m_order);
        Problem problem(u, lambda, v, mu);
        problem = stiffness + gauge + constraint - load - natural.over(NaturalAttribute);
        // Retain the already verified Stokes MUMPS policy. In particular,
        // centralized RHS avoids the local distributed-RHS scatter fault.
        for (const auto& [name, value] : {std::pair{"-mat_mumps_icntl_14", "100"},
               std::pair{"-mat_mumps_icntl_20", "0"}})
        {
          PetscBool set = PETSC_FALSE;
          EXPECT_EQ(PetscOptionsHasName(nullptr, nullptr, name, &set), PETSC_SUCCESS);
          if (!set)
            EXPECT_EQ(PetscOptionsSetValue(nullptr, name, value), PETSC_SUCCESS);
        }
        PETSc::Solver::KSP solver(problem);
        solver.setType(KSPPREONLY);
        PC pc = nullptr;
        EXPECT_EQ(KSPGetPC(solver.getHandle(), &pc), PETSC_SUCCESS);
        EXPECT_EQ(PCSetType(pc, PCLU), PETSC_SUCCESS);
        EXPECT_EQ(PCFactorSetMatSolverType(pc, MATSOLVERMUMPS), PETSC_SUCCESS);
        solver.solve();
        PCFailedReason failure = PC_NOERROR;
        EXPECT_EQ(PCGetFailedReason(pc, &failure), PETSC_SUCCESS);
        EXPECT_EQ(failure, PC_NOERROR);
        checkSystem(problem, solver);
        auto solutionMean = Integral(u.getSolution());
        auto compatibility = Integral(lambda.getSolution());
        solutionMean.setOrder(normOrder);
        compatibility.setOrder(normOrder);
        EXPECT_LT(std::abs(solutionMean.compute()), GaugeTolerance);
        // Integrating -Delta u + lambda = f gives lambda = integral(f+g)
        // on the unit-volume box. Missing flux must violate compatibility.
        const Real defect = std::abs(compatibility.compute());
        if (wrongFlux)
          EXPECT_GT(defect, WrongFluxFloor);
        else
          EXPECT_LT(defect, GaugeTolerance);
        return measure();
#else
        ADD_FAILURE() << "Pure Neumann certification requires PETSc MUMPS";
        return {};
#endif
      }

    private:
      template <class ProblemType, class SolverType>
      void checkSystem(const ProblemType& problem, const SolverType& solver) const
      {
        KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
        EXPECT_EQ(KSPGetConvergedReason(solver.getHandle(), &reason), PETSC_SUCCESS);
        EXPECT_GT(reason, 0);
        const auto& system = problem.getLinearSystem();
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
        EXPECT_LT(norm / std::max(Real(1), rhsNorm), ResidualTolerance);
        EXPECT_EQ(VecDestroy(&residual), PETSC_SUCCESS);
      }

      static constexpr Real RobinCoefficient = 2;
      // Algebraic and gauge budgets are below the finest measured field errors.
      static constexpr Real SolverTolerance = 1e-13;
      static constexpr Real ResidualTolerance = 1e-11;
      static constexpr Real GaugeTolerance = 1e-10;
      static constexpr Real WrongFluxFloor = 1e-2;
      static constexpr Real DivergenceTolerance = 1e5;
      static constexpr size_t MaxIterations = 50000;
      std::reference_wrapper<const MeshType> m_mesh;
      Condition m_condition;
      size_t m_order;
  };
}

#endif
