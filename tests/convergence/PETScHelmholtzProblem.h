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
   * Mixed unit-box data prescribe the trace on @f$x_0=0@f$ and
   * @f$\partial_nu+i\beta u@f$ elsewhere, with @f$\beta=0@f$ or 1.
   * The mixed Poincare bound has @f$\lambda_1=\pi^2/4>1/4@f$.
   * Complex impedance requires GMRES rather than Hermitian CG.
   * @par Architecture
   * Each measurement owns a fresh fixed-layout space and linear system.
   * Independent norm quadrature and a solve-scoped const observer support
   * flat and mapped domains without duplicating field evaluation or MPI norms.
   */
  template <size_t K, class MeshType>
  class PETScHelmholtzProblem
  {
    public:
      enum class Boundary
      {
        Dirichlet,
        MixedNeumann,
        Impedance
      };
      static constexpr Geometry::Attribute DirichletAttribute = 301;
      static constexpr Geometry::Attribute NaturalAttribute = 302;

      PETScHelmholtzProblem(const MeshType& mesh, HelmholtzData::Field field,
        size_t assemblyOrder = 16, Boundary boundary = Boundary::Dirichlet)
        : m_mesh(mesh),
          m_data(mesh.getDimension(), field),
          m_order(assemblyOrder),
          m_boundary(boundary)
      {
        assert(mesh.getDimension() > 0);
        assert(mesh.getDimension() == mesh.getSpaceDimension());
      }

      /** @brief Assemble, solve, and measure the global manufactured problem.
       * @note For an MPI mesh, all ranks of its communicator must participate:
       * distributed assembly, the linear solve, residual norms, and error
       * norms are global operations. Manufactured pointwise data remain local.
       */
      ErrorNorms solve(bool omitMass = false, Real tolerance = SolverTolerance,
        size_t normOrder = 18, bool omitFlux = false) const
      {
        return solve(
          omitMass, tolerance, normOrder, [](const auto&, const auto&) {}, omitFlux);
      }

      /** @brief Solve with a solve-scoped observer before error integration.
       * @note The same collective contract as the ordinary solve applies.
       * The observer is invoked on each rank with that rank's solution view;
       * it must not assume an independently replicated global field.
       */
      template <class Observer>
      ErrorNorms solve(bool omitMass, Real tolerance, size_t normOrder,
        Observer&& observe, bool omitFlux = false) const
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
        if (m_boundary == Boundary::Dirichlet)
          problem = stiffness - mass - load + DirichletBC(u, exact);
        else
        {
          const auto gradient = m_data.getGradient();
          const Variational::BoundaryNormal normal(mesh);
          const size_t dim = mesh.getDimension();
          const ComplexFunction flux([gradient, normal, dim](const Geometry::Point& p) {
            const auto derivative = gradient(p);
            if (dim == 1)
              return (p(0) < Real(0.5) ? Real(-1) : Real(1)) * derivative(0);
            const auto outward = normal(p);
            Complex result = 0;
            // Physical normal contraction is bilinear, not a Hermitian dot
            // product: conjugating the complex gradient would change the PDE.
            for (size_t j = 0; j < dim; ++j)
              result += derivative(j) * outward(j);
            return result;
          });
          const Complex impedance =
            m_boundary == Boundary::Impedance ? Complex(0, 1) : Complex(0);
          auto boundaryMass = BoundaryIntegral(impedance * u, v);
          auto boundaryLoad = BoundaryIntegral(
            (omitFlux ? Real(0) : Real(1)) * flux + impedance * exact, v);
          boundaryMass.setOrder(m_order);
          boundaryLoad.setOrder(m_order);
          problem = stiffness - mass + boundaryMass.over(NaturalAttribute) - load -
            boundaryLoad.over(NaturalAttribute) +
            DirichletBC(u, exact).on(DirichletAttribute);
        }
        PETSc::Solver::KSP solver(problem);
        solver.setType(m_boundary == Boundary::Impedance ? KSPGMRES : KSPCG);
        if (m_boundary == Boundary::Impedance)
        {
          PC pc = nullptr;
          EXPECT_EQ(KSPGetPC(solver.getHandle(), &pc), PETSC_SUCCESS);
          EXPECT_EQ(PCSetType(pc, PCJACOBI), PETSC_SUCCESS);
        }
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
      Boundary m_boundary;
  };
}

#endif
