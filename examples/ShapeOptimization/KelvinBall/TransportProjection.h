/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef KELVIN_BALL_TRANSPORT_PROJECTION_H
#define KELVIN_BALL_TRANSPORT_PROJECTION_H

#include <Eigen/LU>
#include <Rodin/Variational.h>

#include "Common.h"

namespace KelvinBall
{
  /**
   * Chamber-local P1 transport projection with matched quadrature omission.
   *
   * Architecture: trace each sample once, form the retained local Gram matrix,
   * and expose matching coefficients for the mass and load integrals. A
   * rank-deficient cell retains its old distance through the full local mass
   * form. The global projection and trace coupling remain continuous, so this
   * fallback supplies old-distance data rather than fixing nodal values.
   */
  class TransportProjection
  {
    public:
      template <class Mesh, class Distance, class Flow>
      TransportProjection(
        const Mesh& mesh, const Distance& distance, const Flow& flow, size_t order)
        : m_samples(mesh.getCellCount())
      {
        for (auto cell = mesh.getPolytope(mesh.getDimension()); cell; ++cell)
        {
          const auto& formula =
            QF::PolytopeQuadratureFormula::get(order, cell->getGeometry());
          const auto& quadrature = cell->getQuadrature(formula);
          const Variational::H1Element<1, Real> element(cell->getGeometry());
          Math::Matrix<Real> gram(element.getCount(), element.getCount());
          gram.setZero();
          auto& samples = m_samples[cell->getIndex()];
          for (size_t q = 0; q < quadrature.getSize(); ++q)
          {
            const auto& point = quadrature.getPoint(q);
            const auto trace = flow.trace(point);
            const Real weight = formula.getWeight(q) * point.getDistortion();
            const bool retained = !trace.exited();
            samples.push_back(
              {point.getReferenceCoordinates(), retained ? Real(1) : Real(0),
                retained ? distance.getValue(trace.getPoint()) + trace.getCorrection()
                         : Real(0)});
            ++m_attempted;
            m_totalWeight += weight;
            if (!retained)
            {
              ++m_omitted;
              m_omittedWeight += weight;
              continue;
            }
            Math::Vector<Real> basis(element.getCount());
            for (size_t i = 0; i < element.getCount(); ++i)
              basis(i) = element.getBasis(i)(point.getReferenceCoordinates());
            gram.noalias() += weight * basis * basis.transpose();
          }
          if (gram.fullPivLu().rank() != gram.rows())
          {
            ++m_fallbackCells;
            for (size_t q = 0; q < samples.size(); ++q)
            {
              samples[q].mask = 1;
              samples[q].value = distance.getValue(quadrature.getPoint(q));
            }
          }
        }
      }

      Real mask(const Geometry::Point& point) const
      {
        return sample(point).mask;
      }
      Real value(const Geometry::Point& point) const
      {
        return sample(point).value;
      }
      size_t getAttemptedCount() const
      {
        return m_attempted;
      }
      size_t getOmittedCount() const
      {
        return m_omitted;
      }
      size_t getFallbackCellCount() const
      {
        return m_fallbackCells;
      }
      Real getOmittedWeightFraction() const
      {
        return m_totalWeight > 0 ? m_omittedWeight / m_totalWeight : Real(0);
      }

    private:
      struct Sample
      {
          Math::SpatialPoint reference;
          Real mask;
          Real value;
      };

      const Sample& sample(const Geometry::Point& point) const
      {
        // Both forms use exactly the sampling quadrature order. This lookup
        // avoids tracing twice and is read-only during parallel assembly.
        for (const auto& entry : m_samples.at(point.getPolytope().getIndex()))
          if ((entry.reference - point.getReferenceCoordinates()).squaredNorm() == 0)
            return entry;
        throw std::runtime_error(
          "Transport projection quadrature does not match its samples.");
      }

      std::vector<std::vector<Sample>> m_samples;
      size_t m_attempted = 0;
      size_t m_omitted = 0;
      size_t m_fallbackCells = 0;
      Real m_totalWeight = 0;
      Real m_omittedWeight = 0;
  };
}

#endif
