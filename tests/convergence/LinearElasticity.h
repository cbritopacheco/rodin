/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Manufactured vector elasticity data shared by p and hp studies. */

#ifndef RODIN_TESTS_CONVERGENCE_LINEARELASTICITY_H
#define RODIN_TESTS_CONVERGENCE_LINEARELASTICITY_H

#include <cmath>
#include <cstdint>
#include <functional>

#include <gtest/gtest.h>

#include "Convergence.h"
#include "Rodin/Assembly.h"
#include "Rodin/Solver/CG.h"

namespace Rodin::Tests::Convergence::LinearElasticity
{
  /**
   * @brief Exact isotropic-elasticity field on the unit box.
   *
   * With @f$u_i=(i+1)e^{\sum_jx_j}@f$ and
   * @f$\sigma=\lambda(\nabla\cdot u)I+\mu(\nabla u+\nabla u^T)@f$,
   * the source is @f$f_i=-e^{\sum_jx_j}
   * [\mu d(i+1)+(\lambda+\mu)\sum_j(j+1)]@f$.
   */
  class ManufacturedSolution
  {
    private:
      size_t m_dim;
      Real m_lambda;
      Real m_mu;

    public:
      ManufacturedSolution(size_t dim, Real lambda, Real mu)
        : m_dim(dim),
          m_lambda(lambda),
          m_mu(mu),
          exact(dim,
            [dim](const Geometry::Point& p) {
              Real exponent = 0;
              for (size_t j = 0; j < dim; ++j)
                exponent += p(j);
              Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
              for (size_t i = 0; i < dim; ++i)
                value(i) = Real(i + 1) * std::exp(exponent);
              return value;
            }),
          forcing(dim, [dim, lambda, mu](const Geometry::Point& p) {
            Real exponent = 0;
            Real coefficientSum = 0;
            for (size_t j = 0; j < dim; ++j)
            {
              exponent += p(j);
              coefficientSum += Real(j + 1);
            }
            Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
            for (size_t i = 0; i < dim; ++i)
              value(i) = -std::exp(exponent) *
                (mu * Real(dim) * Real(i + 1) + (lambda + mu) * coefficientSum);
            return value;
          })
      {}

      Math::SpatialMatrix<Real> jacobian(const Geometry::Point& p) const
      {
        Real exponent = 0;
        for (size_t j = 0; j < m_dim; ++j)
          exponent += p(j);
        Math::SpatialMatrix<Real> value(
          static_cast<std::uint8_t>(m_dim), static_cast<std::uint8_t>(m_dim));
        for (size_t i = 0; i < m_dim; ++i)
          for (size_t j = 0; j < m_dim; ++j)
            value(i, j) = Real(i + 1) * std::exp(exponent);
        return value;
      }

      size_t getDimension() const
      {
        return m_dim;
      }
      Real getLambda() const
      {
        return m_lambda;
      }
      Real getMu() const
      {
        return m_mu;
      }

      Variational::VectorFunction<
        std::function<Math::SpatialVector<Real>(const Geometry::Point&)>>
        exact;
      Variational::VectorFunction<
        std::function<Math::SpatialVector<Real>(const Geometry::Point&)>>
        forcing;
  };

  template <size_t K>
  ErrorNorms solve(const Geometry::LocalMesh& mesh, const ManufacturedSolution& data)
  {
    using namespace Variational;
    const auto evaluate = [&](auto& space) {
      TrialFunction u(space);
      TestFunction v(space);
      auto volumetric = Integral(data.getLambda() * Div(u), Div(v));
      auto shear = Integral(data.getMu() * (Jacobian(u) + Jacobian(u).T()),
        0.5 * (Jacobian(v) + Jacobian(v).T()));
      auto body = Integral(data.forcing, v);
      volumetric.setOrder(12);
      shear.setOrder(12);
      body.setOrder(12);
      Problem problem(u, v);
      problem = volumetric + shear - body + DirichletBC(u, data.exact);
      Solver::CG solver(problem);
      solver.setTolerance(1e-13).setMaxIterations(50000).solve();
      EXPECT_TRUE(solver.success());
      return ErrorNorm::computeVector(
        mesh, u.getSolution(), data.exact,
        [&data](const Geometry::Point& p) { return data.jacobian(p); }, 12);
    };
    if constexpr (K == 1)
    {
      P1 space(mesh, data.getDimension());
      return evaluate(space);
    }
    else
    {
      H1 space(std::integral_constant<size_t, K>{}, mesh, data.getDimension());
      return evaluate(space);
    }
  }
}

#endif
