/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Manufactured vector elasticity data shared by convergence studies. */

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
   * The constant field @f$u_i=1@f$ has zero source and Jacobian.
   * The default exponential field is retained for the existing p/hp studies.
   * Optional affine and quadratic fields use
   * @f$u_i=1+(i+1)s@f$ and @f$u_i=1+(i+1)s^2@f$, respectively,
   * where @f$s=\sum_jx_j@f$. Their sources are zero and
   * @f$-2[\mu d(i+1)+(\lambda+\mu)\sum_j(j+1)]@f$.
   * The asymmetric affine oracle uses @f$u=\mathbf{1}+Ax@f$ with
   * @f$A_{ij}=(i+1)(j+1)+\delta_{i0}\delta_{j,d-1}@f$ and zero source.
   * Coordinate overloads permit exact-domain lift evaluation without constructing
   * an artificial Geometry::Point; Point-based field factories delegate to them.
   * For @f$d\ge2@f$, the divergence-free shear field is
   * @f$u=(\sin(\pi x_1),0,\ldots,0)@f$ with
   * @f$f=\mu\pi^2u@f$. Its exact field and stress are independent of lambda;
   * it supports a case-specific nearly incompressible convergence study.
   */
  class ManufacturedSolution
  {
    public:
      enum class Field
      {
        Constant,
        Affine,
        Quadratic,
        Exponential,
        AsymmetricAffine,
        DivergenceFree
      };

    private:
      size_t m_dim;
      Real m_lambda;
      Real m_mu;
      Field m_field;

      static Math::SpatialVector<Real> solution(
        size_t dim, Field field, const Math::SpatialPoint& x)
      {
        Real exponent = 0;
        for (size_t j = 0; j < dim; ++j)
          exponent += x(j);
        Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
        if (field == Field::DivergenceFree)
        {
          assert(dim > 1);
          value.setZero();
          value(0) = std::sin(Math::Constants::pi() * x(1));
          return value;
        }
        for (size_t i = 0; i < dim; ++i)
        {
          if (field == Field::AsymmetricAffine)
          {
            value(i) = 1;
            for (size_t j = 0; j < dim; ++j)
              value(i) += (Real((i + 1) * (j + 1)) + (i == 0 && j == dim - 1)) * x(j);
          }
          else
            value(i) = field == Field::Constant ? Real(1)
              : field == Field::Exponential     ? Real(i + 1) * std::exp(exponent)
                                                : 1 +
                Real(i + 1) * (field == Field::Affine ? exponent : exponent * exponent);
        }
        return value;
      }

    public:
      ManufacturedSolution(
        size_t dim, Real lambda, Real mu, Field field = Field::Exponential)
        : m_dim(dim),
          m_lambda(lambda),
          m_mu(mu),
          m_field(field),
          exact(dim,
            [dim, field](const Geometry::Point& p) {
              return solution(dim, field, p.getPhysicalCoordinates());
            }),
          forcing(dim, [dim, lambda, mu, field](const Geometry::Point& p) {
            Real exponent = 0;
            Real coefficientSum = 0;
            for (size_t j = 0; j < dim; ++j)
            {
              exponent += p(j);
              coefficientSum += Real(j + 1);
            }
            Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
            for (size_t i = 0; i < dim; ++i)
              value(i) = field == Field::DivergenceFree
                ? (i == 0 ? mu * Math::Constants::pi() * Math::Constants::pi() *
                        std::sin(Math::Constants::pi() * p(1))
                          : Real(0))
                : -(field == Field::Exponential   ? std::exp(exponent)
                      : field == Field::Quadratic ? Real(2)
                                                  : Real(0)) *
                  (mu * Real(dim) * Real(i + 1) + (lambda + mu) * coefficientSum);
            return value;
          })
      {
        assert(field != Field::DivergenceFree || dim > 1);
      }

      Math::SpatialMatrix<Real> jacobian(const Geometry::Point& p) const
      {
        return getJacobian(p.getPhysicalCoordinates());
      }

      Math::SpatialVector<Real> getSolution(const Math::SpatialPoint& x) const
      {
        return solution(m_dim, m_field, x);
      }

      Math::SpatialMatrix<Real> getJacobian(const Math::SpatialPoint& x) const
      {
        Real exponent = 0;
        for (size_t j = 0; j < m_dim; ++j)
          exponent += x(j);
        Math::SpatialMatrix<Real> value(
          static_cast<std::uint8_t>(m_dim), static_cast<std::uint8_t>(m_dim));
        if (m_field == Field::DivergenceFree)
        {
          assert(m_dim > 1);
          value.setZero();
          value(0, 1) = Math::Constants::pi() * std::cos(Math::Constants::pi() * x(1));
          return value;
        }
        for (size_t i = 0; i < m_dim; ++i)
          for (size_t j = 0; j < m_dim; ++j)
            value(i, j) = m_field == Field::AsymmetricAffine
              ? Real((i + 1) * (j + 1)) + (i == 0 && j == m_dim - 1)
              : Real(i + 1) *
                (m_field == Field::Exponential    ? std::exp(exponent)
                    : m_field == Field::Quadratic ? 2 * exponent
                    : m_field == Field::Constant  ? Real(0)
                                                  : Real(1));
        return value;
      }

      size_t getDimension() const
      {
        return m_dim;
      }

      /** @brief Analytic symmetric gradient, independent of discrete operators. */
      Math::SpatialMatrix<Real> strain(const Geometry::Point& p) const
      {
        if (m_field == Field::DivergenceFree)
        {
          Math::SpatialMatrix<Real> value(m_dim, m_dim);
          value.setZero();
          value(0, 1) = value(1, 0) =
            Math::Constants::pi() * std::cos(Math::Constants::pi() * p(1)) / 2;
          return value;
        }
        Real sum = 0;
        for (size_t j = 0; j < m_dim; ++j)
          sum += p(j);
        const Real factor = m_field == Field::Exponential ? std::exp(sum)
          : m_field == Field::Quadratic                   ? 2 * sum
          : m_field == Field::Constant                    ? Real(0)
                                                          : Real(1);
        Math::SpatialMatrix<Real> value(
          static_cast<std::uint8_t>(m_dim), static_cast<std::uint8_t>(m_dim));
        for (size_t i = 0; i < m_dim; ++i)
          for (size_t j = 0; j < m_dim; ++j)
            value(i, j) = m_field == Field::AsymmetricAffine
              ? Real((i + 1) * (j + 1)) +
                Real((i == 0 && j == m_dim - 1) + (j == 0 && i == m_dim - 1)) / 2
              : Real(i + j + 2) * factor / 2;
        return value;
      }

      /** @brief Analytic Cauchy stress for the stated isotropic law. */
      Math::SpatialMatrix<Real> stress(const Geometry::Point& p) const
      {
        if (m_field == Field::DivergenceFree)
        {
          Math::SpatialMatrix<Real> value(m_dim, m_dim);
          value.setZero();
          value(0, 1) = value(1, 0) =
            m_mu * Math::Constants::pi() * std::cos(Math::Constants::pi() * p(1));
          return value;
        }
        Real sum = 0;
        for (size_t j = 0; j < m_dim; ++j)
          sum += p(j);
        const Real factor = m_field == Field::Exponential ? std::exp(sum)
          : m_field == Field::Quadratic                   ? 2 * sum
          : m_field == Field::Constant                    ? Real(0)
                                                          : Real(1);
        const Real coefficientSum = Real(m_dim * (m_dim + 1)) / 2;
        Math::SpatialMatrix<Real> value(
          static_cast<std::uint8_t>(m_dim), static_cast<std::uint8_t>(m_dim));
        for (size_t i = 0; i < m_dim; ++i)
          for (size_t j = 0; j < m_dim; ++j)
            value(i, j) = m_field == Field::AsymmetricAffine
              ? 2 * m_mu * Real((i + 1) * (j + 1)) +
                m_mu * ((i == 0 && j == m_dim - 1) + (j == 0 && i == m_dim - 1)) +
                m_lambda *
                  (Real(m_dim * (m_dim + 1) * (2 * m_dim + 1)) / 6 + (m_dim == 1)) *
                  (i == j)
              : factor * (m_mu * Real(i + j + 2) + m_lambda * coefficientSum * (i == j));
        return value;
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
