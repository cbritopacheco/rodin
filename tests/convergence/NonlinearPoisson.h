/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef RODIN_TESTS_CONVERGENCE_NONLINEAR_POISSON_H
#define RODIN_TESTS_CONVERGENCE_NONLINEAR_POISSON_H

#include <functional>
#include <gtest/gtest.h>

#include "Convergence.h"
#include "Rodin/Assembly.h"
#include "Rodin/Solver/NewtonSolver.h"
#include "Rodin/Solver/SparseLU.h"

namespace Rodin::Tests::Convergence
{
  /** @brief Backend-independent continuous data for
   * @f$-\Delta u+u+u^3=f@f$ on the unit box.
   * The amplitude is prescribed independently of the discrete solution.
   */
  class NonlinearPoissonData
  {
    public:
      explicit NonlinearPoissonData(size_t dimension, Real amplitude = 1)
        : m_dimension(dimension),
          m_amplitude(amplitude)
      {}

      auto getSolution() const
      {
        return Variational::RealFunction(
          [dim = m_dimension, amplitude = m_amplitude](const Geometry::Point& x) {
            Real value = amplitude;
            for (size_t j = 0; j < dim; ++j)
              value *= std::sin(Math::Constants::pi() * x(j));
            return value;
          });
      }

      auto getGradient() const
      {
        return Variational::VectorFunction(m_dimension,
          [dim = m_dimension, amplitude = m_amplitude](const Geometry::Point& x) {
            const Real pi = Math::Constants::pi();
            Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
            for (size_t j = 0; j < dim; ++j)
            {
              value(j) = amplitude * pi * std::cos(pi * x(j));
              for (size_t k = 0; k < dim; ++k)
                if (k != j)
                  value(j) *= std::sin(pi * x(k));
            }
            return value;
          });
      }

      auto getSource() const
      {
        return Variational::RealFunction(
          [dim = m_dimension, exact = getSolution()](const Geometry::Point& x) {
            const Real u = exact(x), pi = Math::Constants::pi();
            return (Real(dim) * pi * pi + 1) * u + u * u * u;
          });
      }

    private:
      size_t m_dimension;
      Real m_amplitude;
  };

  /** @brief Native semilinear Poisson verification problem.
   * The common reference is @f$u=\prod_j\sin(\pi x_j)@f$ with
   * @f$f=(d\pi^2+1)u+u^3@f$ and zero trace. Newton iteration starts at
   * zero and uses @f$J(u)=-\Delta+1+3u^2@f$.
   *
   * @par Architecture
   * The mesh and quadrature order define one immutable workload. Analytic
   * callbacks provide continuous data independently of the discrete state.
   * Each solve constructs fresh spaces, fields, a tangential problem, and
   * Newton/SparseLU solvers; no linear system is resized across refinements.
   * The derivative check separately evaluates residuals at two perturbed
   * states and compares their difference with the baseline tangent action.
   */
  class NonlinearPoissonProblem
  {
    public:
      explicit NonlinearPoissonProblem(
        const Geometry::LocalMesh& mesh, size_t quadratureOrder = 16, Real amplitude = 1)
        : m_mesh(mesh),
          m_order(quadratureOrder),
          m_data(mesh.getDimension(), amplitude)
      {}

      auto getSolution() const
      {
        return m_data.getSolution();
      }

      auto getGradient() const
      {
        return m_data.getGradient();
      }

      auto getSource() const
      {
        return m_data.getSource();
      }

      template <size_t K>
      ErrorNorms solve(bool omitCubic = false, Real tolerance = 1e-11) const
      {
        using namespace Variational;
        H1<K, Real> space(std::integral_constant<size_t, K>{}, m_mesh.get());
        GridFunction current(space);
        current = Zero();
        TrialFunction du(space);
        TestFunction v(space);
        const Real gamma = omitCubic ? 0 : 1;
        auto a = Integral(Grad(du), Grad(v));
        auto c = Integral((1 + 3 * gamma * current * current) * du, v);
        auto r = Integral(Grad(current), Grad(v));
        auto s = Integral(current + gamma * current * current * current, v);
        auto b = Integral(getSource(), v);
        a.setOrder(m_order);
        c.setOrder(m_order);
        r.setOrder(m_order);
        s.setOrder(m_order);
        b.setOrder(m_order);
        Problem problem(du, v);
        problem = a + c + r + s - b + DirichletBC(du, Zero());
        Solver::SparseLU linear(problem);
        Solver::NewtonSolver newton(linear);
        newton.setMaxIterations(20).setAbsoluteTolerance(tolerance).setRelativeTolerance(
          tolerance);
        newton.solve(current);
        EXPECT_TRUE(newton.converged());
        EXPECT_TRUE(linear.success());
        EXPECT_TRUE(std::isfinite(newton.getReport().finalResidual));
        return ErrorNorm::compute(
          m_mesh.get(), current, getSolution(), getGradient(), m_order);
      }

      /** @brief Central difference of the assembled residual, compared with
       * the tangent action at a nonzero state. The deliberately incorrect
       * tangent uses @f$1+u^2@f$ instead of @f$1+3u^2@f$.
       */
      template <size_t K>
      Real tangentDefect(bool wrongTangent = false) const
      {
        using namespace Variational;
        H1<K, Real> space(std::integral_constant<size_t, K>{}, m_mesh.get());
        GridFunction state(space), direction(space);
        state = 0.5 * getSolution();
        direction = 0.25 * getSolution();
        const auto baseline = state.getData();
        TrialFunction du(space);
        TestFunction v(space);
        auto a = Integral(Grad(du), Grad(v));
        auto c = Integral((1 + (wrongTangent ? 1 : 3) * state * state) * du, v);
        auto r = Integral(Grad(state), Grad(v));
        auto s = Integral(state + state * state * state, v);
        auto b = Integral(getSource(), v);
        a.setOrder(m_order);
        c.setOrder(m_order);
        r.setOrder(m_order);
        s.setOrder(m_order);
        b.setOrder(m_order);
        Problem problem(du, v);
        problem = a + c + r + s - b + DirichletBC(du, Zero());
        problem.assemble();
        const Math::Vector<Real> action =
          problem.getLinearSystem().getOperator() * direction.getData();
        constexpr Real epsilon = 1e-5;
        state.getData() = baseline + epsilon * direction.getData();
        problem.assemble();
        const Math::Vector<Real> plus = -problem.getLinearSystem().getVector();
        state.getData() = baseline - epsilon * direction.getData();
        problem.assemble();
        const Math::Vector<Real> minus = -problem.getLinearSystem().getVector();
        const Math::Vector<Real> difference = (plus - minus) / (2 * epsilon);
        EXPECT_TRUE(std::isfinite(difference.norm()));
        EXPECT_GT(difference.norm(), 0);
        return (action - difference).norm() / difference.norm();
      }

    private:
      std::reference_wrapper<const Geometry::LocalMesh> m_mesh;
      size_t m_order;
      NonlinearPoissonData m_data;
  };
}

#endif
