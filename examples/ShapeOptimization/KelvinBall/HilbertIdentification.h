/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef KELVIN_BALL_HILBERT_IDENTIFICATION_H
#define KELVIN_BALL_HILBERT_IDENTIFICATION_H

#include "RotatedNitscheIntegrator.h"

namespace KelvinBall
{
  /**
   * @brief Riesz identification with one operator per design.
   *
   * The metric is @f$ a(g,v)=\ell^2\int_D\nabla g:\nabla v+
   * \int_D g\cdot v+j_\Sigma(g,v) @f$, with zero values on the outer wall.
   * Architecture: assemble the metric once; use Rodin RHS-only assembly for
   * successive differentials; retain the MUMPS factorization when available.
   * The object and its finite element space expire before reconstruction.
   */
  class HilbertIdentification
  {
    public:
      struct Diagnostics
      {
          Real residual;
          Real jump;
      };

      using Field = GridFunction<VelocitySpace, Math::Vector<Real>>;
      using Trial = TrialFunction<Field, VelocitySpace>;
      using Test = TestFunction<VelocitySpace>;
      using System = Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;
      using ProblemType =
        decltype(Problem(std::declval<Trial&>(), std::declval<Test&>()));
      using Trace = decltype(std::declval<const RotatedNitscheIntegrator&>().vectorTrace(
        std::declval<const Trial&>(), std::declval<const Test&>(), Real(1), Real(1)));
#ifdef RODIN_USE_MUMPS
      using DirectSolver = Solver::MUMPS<System>;
#elif defined(RODIN_USE_UMFPACK)
      using DirectSolver = Solver::UMFPack<System>;
#else
      using DirectSolver = Solver::SparseLU<System>;
#endif

      HilbertIdentification(const VelocitySpace& space,
        const RotatedNitscheIntegrator& coupling, Real length, Real penalty)
        : m_coupling(coupling),
          m_length(length),
          m_trial(space),
          m_test(space),
          m_trace(coupling.vectorTrace(m_trial, m_test, length * length, penalty)),
          m_problem(m_trial, m_test),
          m_solver(m_problem)
      {
        m_problem = Integral(length * length * Jacobian(m_trial), Jacobian(m_test)) +
          m_trace + Integral(m_trial, m_test) +
          DirichletBC(m_trial, VectorFunction{0, 0, 0}).on(Outer);
        m_problem.assemble();
#ifdef RODIN_USE_MUMPS
        m_solver.setSymmetric(DirectSolver::Symmetry::General);
#endif
      }

      template <class Density>
      Diagnostics identify(const Density& density, Field& gradient,
        const Math::Vector<Real>* additionalLoad = nullptr)
      {
        auto normal = FaceNormal(m_trial.getFiniteElementSpace().getMesh());
        normal.traceOf(Fluid);
        m_problem = Integral(m_length * m_length * Jacobian(m_trial), Jacobian(m_test)) +
          m_trace + Integral(m_trial, m_test) -
          FaceIntegral(density, Dot(normal, m_test)).over(Gamma) +
          DirichletBC(m_trial, VectorFunction{0, 0, 0}).on(Outer);
        m_problem.assemble(AssemblyTarget::RHS);
        auto& system = m_problem.getLinearSystem();
        if (additionalLoad)
          system.getVector() += *additionalLoad;
        m_problem.solve(m_solver);
        if (!m_solver.success())
          throw std::runtime_error(
            "Hilbert identification factorization or solve failed.");
        const Real residual =
          (system.getOperator() * system.getSolution() - system.getVector()).norm() /
          std::max(system.getVector().norm(), Real(1));
        if (!std::isfinite(residual) || residual > LinearResidualTolerance)
          throw std::runtime_error(
            "The gradient identification problem did not converge.");
        gradient = m_trial.getSolution();
        return {residual, m_coupling.get().vectorJump(gradient)};
      }

    private:
      std::reference_wrapper<const RotatedNitscheIntegrator> m_coupling;
      Real m_length;
      Trial m_trial;
      Test m_test;
      Trace m_trace;
      ProblemType m_problem;
      DirectSolver m_solver;
  };
}

#endif
