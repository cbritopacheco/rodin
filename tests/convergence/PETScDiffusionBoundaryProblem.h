/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef RODIN_TESTS_CONVERGENCE_PETSC_DIFFUSIONBOUNDARYPROBLEM_H
#define RODIN_TESTS_CONVERGENCE_PETSC_DIFFUSIONBOUNDARYPROBLEM_H

#include <gtest/gtest.h>

#include "Conductivity.h"
#include "Rodin/PETSc.h"

namespace Rodin::Tests::Convergence
{
  /** @brief Scalar diffusion with essential, natural and Robin data.
   * @par Mathematical formulation
   * On the unit box, @f$-\nabla\cdot(\gamma\nabla u)=f@f$ and
   * @f$u=e^{\sum_jx_j}@f$, with @f$\gamma=1@f$ for Poisson or
   * @f$\gamma=1+\sum_jx_j\ge1@f$ for variable conductivity.
   * Natural data are @f$g=\gamma\nabla u\cdot n@f$ or
   * @f$r=g+2u@f$ on the natural partition.
   * Pure Neumann data use @f$u=e^{\sum_jx_j}-(e-1)^d@f$ and a P0g
   * multiplier imposing @f$\int_\Omega u=0@f$.
   * Affine and quadratic patches use their analytic unit-box means instead.
   * @par Architecture
   * Each solve owns its fixed-layout space and system. Boundary attributes
   * are prescribed before partitioning; no entity is matched by coordinates.
   * ErrorNorm integrates the field independently. Solver status, coefficient
   * residual, mean and compatibility multiplier have separate acceptance.
   * Distributed assembly, solution and global error/residual integration
   * have collective semantics; prescribing the reference mean does not.
   * An explicit reference mean permits other prescribed domains without
   * changing the default analytic unit-box path.
   */
  template <size_t K, class MeshType, bool ConstantCoefficient = true>
  class PETScDiffusionBoundaryProblem
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

      /** @brief Bind the mesh and prescribed manufactured boundary data.
       * @param referenceMean Optional physical-domain mean for pure Neumann
       * data; omission retains the analytic unit-box mean. Supplying this
       * datum neither integrates the mesh nor communicates between ranks.
       */
      PETScDiffusionBoundaryProblem(const MeshType& mesh, Condition condition,
        size_t assemblyOrder = 16,
        ConductivityData::Field field = ConductivityData::Field::Exponential,
        Optional<Real> referenceMean = std::nullopt)
        : m_mesh(mesh),
          m_condition(condition),
          m_order(assemblyOrder),
          m_field(field),
          m_referenceMean(referenceMean)
      {}

      /** @brief Assemble, solve, and measure the global boundary problem.
       * @note For an MPI mesh, all ranks of its communicator must participate
       * in assembly, the solve, and the global residual and error norms.
       * Pointwise coefficients, sources, and the analytic mean are local data.
       */
      ErrorNorms solve(bool wrongFlux = false, size_t normOrder = 18,
        Real tolerance = SolverTolerance) const
      {
        using namespace Variational;
        const auto& mesh = m_mesh.get();
        const size_t dim = mesh.getDimension();
        const bool pure = m_condition == Condition::PureNeumann;
        const ConductivityData data(dim, m_field);
        const auto solution = data.getSolution();
        const Real mean = pure ? m_referenceMean.value_or(getUnitBoxMean(dim)) : 0;
        const RealFunction exact(
          [solution, mean](const Geometry::Point& p) { return solution(p) - mean; });
        const auto gradient = data.getGradient();
        const auto coefficient = data.getCoefficient(ConstantCoefficient);
        const BoundaryNormal normal(mesh);
        const RealFunction flux(
          [dim, gradient, coefficient, normal](const Geometry::Point& p) {
            const auto derivative = gradient(p);
            // Unit-box endpoint classification, not a correspondence between meshes.
            const Real normalDerivative = dim == 1
              ? (p(0) < Real(0.5) ? -1 : 1) * derivative(0)
              : derivative.dot(normal(p));
            if constexpr (ConstantCoefficient)
              return normalDerivative;
            else
              return coefficient(p) * normalDerivative;
          });
        const Real scale = wrongFlux ? 0 : 1;
        H1<K, Real, MeshType> space(std::integral_constant<size_t, K>{}, mesh);
        PETSc::Variational::TrialFunction u(space);
        PETSc::Variational::TestFunction v(space);
        auto stiffness = [&] {
          if constexpr (ConstantCoefficient)
            return Integral(Grad(u), Grad(v));
          else
            return Integral(coefficient * Grad(u), Grad(v));
        }();
        auto load = Integral(data.getSource(ConstantCoefficient), v);
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
#if defined(RODIN_POISSON_BOUNDARY_MUMPS) || defined(RODIN_DIFFUSION_BOUNDARY_MUMPS)
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
          {
            EXPECT_EQ(PetscOptionsSetValue(nullptr, name, value), PETSC_SUCCESS);
          }
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
        // Integrating -div(gamma grad u) + lambda = f gives
        // integral(lambda) = integral(f) + boundary_integral(g), for any volume.
        // Missing flux must violate compatibility.
        const Real defect = std::abs(compatibility.compute());
        if (wrongFlux)
        {
          EXPECT_GT(defect, WrongFluxFloor);
        }
        else
          EXPECT_LT(defect, GaugeTolerance);
        return measure();
#else
        ADD_FAILURE() << "Pure Neumann certification requires PETSc MUMPS";
        return {};
#endif
      }

    private:
      // Analytic unit-box integrals. No distributed numerical reduction is
      // needed to prescribe the manufactured field's additive constant.
      Real getUnitBoxMean(size_t dim) const
      {
        switch (m_field)
        {
          case ConductivityData::Field::Constant:
            return 1;
          case ConductivityData::Field::Affine:
            return 1 + Real(dim) / 2;
          case ConductivityData::Field::Quadratic:
            return 1 + Real(dim) / 3;
          case ConductivityData::Field::Smooth:
            return 1 + std::pow(2 / Math::Constants::pi(), Real(dim));
          case ConductivityData::Field::Exponential:
            return std::pow(std::expm1(Real(1)), Real(dim));
        }
        assert(false);
        return 0;
      }

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
      ConductivityData::Field m_field;
      Optional<Real> m_referenceMean;
  };
}

#endif
