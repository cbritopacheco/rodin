/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef RODIN_TESTS_CONVERGENCE_STOKES_PROBLEM_H
#define RODIN_TESTS_CONVERGENCE_STOKES_PROBLEM_H

#include <functional>
#include <gtest/gtest.h>

#include "Stokes.h"
#include "Rodin/Assembly.h"
#include "Rodin/Solver/SparseLU.h"

namespace Rodin::Tests::Convergence
{
  /** @brief Native mixed Stokes workload shared by h, p, and hp studies.
   * @par Architecture
   * One immutable mesh and continuous data set define the workload. Each
   * solve constructs velocity degree @f$K@f$, pressure degree @f$K-1@f$,
   * and a global pressure-mean multiplier in fresh spaces and a fresh linear
   * system. The same form is used at every degree; no DOF layout is assumed.
   * SparseLU solves the saddle-point problem. Algebraic residual and gauge
   * diagnostics precede independent physical-cell field-error integration.
   */
  class StokesProblem
  {
    public:
      StokesProblem(const Geometry::LocalMesh& mesh, const StokesData& data,
        size_t quadratureOrder = 16)
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
        H1 velocitySpace(std::integral_constant<size_t, K>{}, mesh, mesh.getDimension());
        H1 pressureSpace(std::integral_constant<size_t, K - 1>{}, mesh);
        P0g meanSpace(mesh);
        TrialFunction u(velocitySpace);
        TrialFunction p(pressureSpace);
        TrialFunction lambda(meanSpace);
        TestFunction v(velocitySpace);
        TestFunction q(pressureSpace);
        TestFunction mu(meanSpace);
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
        Solver::SparseLU solver(problem);
        solver.solve();
        EXPECT_TRUE(solver.success());
        const auto& system = problem.getLinearSystem();
        const Real residual =
          (system.getOperator() * system.getSolution() - system.getVector()).norm() /
          std::max(Real(1), system.getVector().norm());
        EXPECT_TRUE(std::isfinite(residual));
        EXPECT_LT(residual, 1e-11);
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
      std::reference_wrapper<const Geometry::LocalMesh> m_mesh;
      StokesData m_data;
      size_t m_order;
  };
}

#endif
