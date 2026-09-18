/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Equation- and refinement-independent convergence test utilities.
 */

#ifndef RODIN_TESTS_CONVERGENCE_CONVERGENCE_H
#define RODIN_TESTS_CONVERGENCE_CONVERGENCE_H

#include <cassert>
#include <cmath>
#include <string_view>
#include <type_traits>
#include <vector>

#include "Rodin/Geometry.h"
#include "Rodin/QF.h"
#include "Rodin/Variational.h"

namespace Rodin::Tests::Convergence
{
  class ErrorNorms
  {
    public:
      ErrorNorms(Real l2, Real h1Semi)
        : m_l2(l2), m_h1Semi(h1Semi)
      {}

      Real getL2() const { return m_l2; }
      Real getH1Seminorm() const { return m_h1Semi; }

      bool isFinite() const
      {
        return std::isfinite(m_l2) && std::isfinite(m_h1Semi);
      }

    private:
      Real m_l2;
      Real m_h1Semi;
  };

  class Rates
  {
    public:
      Rates(Real l2, Real h1Semi)
        : m_l2(l2), m_h1Semi(h1Semi)
      {}

      Real getL2() const { return m_l2; }
      Real getH1Seminorm() const { return m_h1Semi; }

    private:
      Real m_l2;
      Real m_h1Semi;
  };

  /**
   * @brief Error samples for algebraic or exponential convergence studies.
   *
   * `getAlgebraicRates()` evaluates
   * @f$\log(e_c/e_f)/\log(s_c/s_f)@f$ for a decreasing discretization scale
   * such as mesh size. `getExponentialRates()` evaluates
   * @f$\log(e_c/e_f)/(p_f-p_c)@f$ for an increasing parameter such as degree.
   */
  class ErrorHistory
  {
    public:
      struct Sample
      {
        Real parameter;
        ErrorNorms error;
      };

      ErrorHistory& append(Real parameter, const ErrorNorms& error)
      {
        m_samples.push_back({parameter, error});
        return *this;
      }

      size_t getSize() const { return m_samples.size(); }

      const Sample& getSample(size_t i) const
      {
        return m_samples.at(i);
      }

      Rates getAlgebraicRates(size_t fineSample) const
      {
        const auto& coarse = getCoarse(fineSample);
        const auto& fine = getFine(fineSample);
        assert(coarse.parameter > fine.parameter);
        return getRates(
          coarse, fine, std::log(coarse.parameter / fine.parameter));
      }

      Rates getExponentialRates(size_t fineSample) const
      {
        const auto& coarse = getCoarse(fineSample);
        const auto& fine = getFine(fineSample);
        assert(fine.parameter > coarse.parameter);
        return getRates(coarse, fine, fine.parameter - coarse.parameter);
      }

    private:
      const Sample& getCoarse(size_t fineSample) const
      {
        assert(fineSample > 0);
        assert(fineSample < m_samples.size());
        return m_samples.at(fineSample - 1);
      }

      const Sample& getFine(size_t fineSample) const
      {
        assert(fineSample > 0);
        assert(fineSample < m_samples.size());
        return m_samples.at(fineSample);
      }

      static Rates getRates(
        const Sample& coarse, const Sample& fine, Real denominator)
      {
        assert(denominator > 0);
        assert(coarse.error.getL2() > 0);
        assert(fine.error.getL2() > 0);
        assert(coarse.error.getH1Seminorm() > 0);
        assert(fine.error.getH1Seminorm() > 0);
        return {
          std::log(coarse.error.getL2() / fine.error.getL2()) / denominator,
          std::log(
            coarse.error.getH1Seminorm() / fine.error.getH1Seminorm())
            / denominator
        };
      }

      std::vector<Sample> m_samples;
  };

  /**
   * @brief Unit-box UniformGrid factory shared by refinement strategies.
   */
  class UniformGrid
  {
    public:
      explicit UniformGrid(Geometry::Polytope::Type geometry)
        : m_geometry(geometry)
      {
        assert(getDimension() >= 1);
        assert(getDimension() <= 3);
      }

      size_t getDimension() const
      {
        return Geometry::Polytope::Traits(m_geometry).getDimension();
      }

      Geometry::Polytope::Type getGeometry() const { return m_geometry; }

      Geometry::LocalMesh makeMesh(size_t pointsPerAxis) const
      {
        assert(pointsPerAxis >= 2);
        const size_t dim = getDimension();
        Geometry::LocalMesh mesh;
        switch (dim)
        {
          case 1:
            mesh = Geometry::LocalMesh::UniformGrid(
              m_geometry, {pointsPerAxis});
            break;
          case 2:
            mesh = Geometry::LocalMesh::UniformGrid(
              m_geometry, {pointsPerAxis, pointsPerAxis});
            break;
          case 3:
            mesh = Geometry::LocalMesh::UniformGrid(
              m_geometry,
              {pointsPerAxis, pointsPerAxis, pointsPerAxis});
            break;
          default:
            assert(false);
            return mesh;
        }

        mesh.scale(Real(1) / Real(pointsPerAxis - 1));
        if (dim == 3)
        {
          mesh.getConnectivity().compute(2, 3);
          mesh.getConnectivity().compute(3, 2);
          mesh.getConnectivity().compute(2, 1);
        }
        else if (dim == 2)
        {
          mesh.getConnectivity().compute(1, 2);
          mesh.getConnectivity().compute(2, 1);
        }
        else
        {
          mesh.getConnectivity().compute(0, 1);
          mesh.getConnectivity().compute(1, 0);
        }
        if (dim > 1)
          mesh.getConnectivity().compute(1, 0);
        return mesh;
      }

      static std::string_view getGeometryName(
        Geometry::Polytope::Type geometry)
      {
        using Type = Geometry::Polytope::Type;
        switch (geometry)
        {
          case Type::Segment:       return "Segment";
          case Type::Triangle:      return "Triangle";
          case Type::Quadrilateral: return "Quadrilateral";
          case Type::Tetrahedron:   return "Tetrahedron";
          case Type::Pyramid:       return "Pyramid";
          case Type::Hexahedron:    return "Hexahedron";
          case Type::Wedge:         return "Wedge";
          default:                  return "Unsupported";
        }
      }

    private:
      Geometry::Polytope::Type m_geometry;
  };

  /**
   * @brief Attribute partition of the boundary of a unit box.
   *
   * The lower and upper facets perpendicular to a selected coordinate axis
   * receive their own attributes. Every other exterior facet receives the
   * remainder attribute. This gives convergence suites a geometry-independent
   * way to define mixed boundary conditions on every UniformGrid cell type.
   */
  class UnitBoxBoundary
  {
    public:
      static void labelCoordinatePartition(
        Geometry::LocalMesh& mesh,
        size_t axis,
        Geometry::Attribute lowerAttribute,
        Geometry::Attribute upperAttribute,
        Geometry::Attribute remainderAttribute)
      {
        assert(axis < mesh.getSpaceDimension());
        constexpr Real tolerance = 1e-12;
        for (auto boundary = mesh.getBoundary(); boundary; ++boundary)
        {
          Real coordinate = 0;
          size_t vertexCount = 0;
          for (const auto vertex : boundary->getVertices())
          {
            coordinate += mesh.getVertexCoordinates(vertex)(axis);
            ++vertexCount;
          }
          assert(vertexCount > 0);
          coordinate /= Real(vertexCount);

          if (std::abs(coordinate) < tolerance)
            mesh.setAttribute(boundary.key(), lowerAttribute);
          else if (std::abs(coordinate - 1) < tolerance)
            mesh.setAttribute(boundary.key(), upperAttribute);
          else
            mesh.setAttribute(boundary.key(), remainderAttribute);
        }
      }
  };

  class ErrorNorm
  {
    public:
      template <class GF, class Exact, class ExactGradient>
      static ErrorNorms compute(
        const Geometry::LocalMesh& mesh,
        const GF& uh,
        const Exact& exact,
        const ExactGradient& exactGradient,
        size_t quadratureOrder = 8)
      {
        const auto gradient = Variational::Grad(uh);
        return computeWithDerivative(
          mesh, uh, gradient, exact, exactGradient, quadratureOrder);
      }

      /**
       * @brief Computes vector L2 and H1-seminorm errors.
       *
       * The value error uses the Euclidean norm and the derivative error uses
       * the Frobenius norm of the displacement Jacobian.
       */
      template <class GF, class Exact, class ExactJacobian>
      static ErrorNorms computeVector(
        const Geometry::LocalMesh& mesh,
        const GF& uh,
        const Exact& exact,
        const ExactJacobian& exactJacobian,
        size_t quadratureOrder = 8)
      {
        const auto jacobian = Variational::Jacobian(uh);
        return computeWithDerivative(
          mesh, uh, jacobian, exact, exactJacobian, quadratureOrder);
      }

      /**
       * @brief Computes the L2 norm of the divergence of a vector field.
       *
       * The Jacobian is evaluated directly because this diagnostic concerns
       * the trace of a vector-field derivative rather than a scalar H1 norm.
       */
      template <class GF>
      static Real computeDivergenceL2(
        const Geometry::LocalMesh& mesh,
        const GF& uh,
        size_t quadratureOrder = 8)
      {
        const auto jacobian = Variational::Jacobian(uh);
        Real squared = 0;
        for (auto cell = mesh.getCell(); cell; ++cell)
        {
          const auto& qf = QF::PolytopeQuadratureFormula::get(
            quadratureOrder, cell->getGeometry());
          const auto& quadrature = cell->getQuadrature(qf);
          for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
          {
            const auto& p = quadrature.getPoint(qp);
            const Variational::IntegrationPoint ip(p, &qf, qp);
            const auto value = jacobian(ip);
            Real divergence = 0;
            const size_t dimension = std::min(
              static_cast<size_t>(value.rows()), static_cast<size_t>(value.cols()));
            for (size_t i = 0; i < dimension; ++i)
              divergence += value(i, i);
            squared += qf.getWeight(qp) * p.getDistortion()
                     * divergence * divergence;
          }
        }
        return std::sqrt(squared);
      }

    private:
      template <class GF, class Derivative, class Exact,
        class ExactDerivative>
      static ErrorNorms computeWithDerivative(
        const Geometry::LocalMesh& mesh,
        const GF& uh,
        const Derivative& derivative,
        const Exact& exact,
        const ExactDerivative& exactDerivative,
        size_t quadratureOrder)
      {
        Real l2Squared = 0;
        Real h1SemiSquared = 0;
        for (auto cell = mesh.getCell(); cell; ++cell)
        {
          const auto& qf = QF::PolytopeQuadratureFormula::get(
            quadratureOrder, cell->getGeometry());
          const auto& quadrature = cell->getQuadrature(qf);
          for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
          {
            const auto& p = quadrature.getPoint(qp);
            const Variational::IntegrationPoint ip(p, &qf, qp);
            const Real weight = qf.getWeight(qp) * p.getDistortion();
            const auto ev = uh(ip) - evaluate(exact, ip);
            const auto eg = derivative(ip) - evaluate(exactDerivative, ip);
            l2Squared += weight * squaredMagnitude(ev);
            h1SemiSquared += weight * squaredMagnitude(eg);
          }
        }
        return {std::sqrt(l2Squared), std::sqrt(h1SemiSquared)};
      }

      template <class F>
      static decltype(auto) evaluate(
        const F& f, const Variational::IntegrationPoint& ip)
      {
        if constexpr (
          std::is_invocable_v<const F&, const Variational::IntegrationPoint&>)
          return f(ip);
        else
          return f(ip.getPoint());
      }

      template <class Scalar>
      static Real squaredMagnitude(const Scalar& value)
      {
        if constexpr (requires { value.rows(); value.cols(); value(0, 0); })
        {
          Real result = 0;
          for (size_t i = 0; i < static_cast<size_t>(value.rows()); ++i)
            for (size_t j = 0; j < static_cast<size_t>(value.cols()); ++j)
              result += squaredMagnitude(value(i, j));
          return result;
        }
        else if constexpr (requires { value.size(); value(0); })
        {
          Real result = 0;
          for (size_t i = 0; i < static_cast<size_t>(value.size()); ++i)
            result += squaredMagnitude(value(i));
          return result;
        }
        else
        {
          const Real magnitude = std::abs(value);
          return magnitude * magnitude;
        }
      }
  };
}

#endif
