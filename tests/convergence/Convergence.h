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
            const auto eg = gradient(ip) - evaluate(exactGradient, ip);
            l2Squared += weight * squaredMagnitude(ev);
            for (size_t i = 0; i < eg.size(); ++i)
              h1SemiSquared += weight * squaredMagnitude(eg(i));
          }
        }
        return {std::sqrt(l2Squared), std::sqrt(h1SemiSquared)};
      }

    private:
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
        const Real magnitude = std::abs(value);
        return magnitude * magnitude;
      }
  };
}

#endif
