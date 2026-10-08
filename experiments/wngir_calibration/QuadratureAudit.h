/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXPERIMENTS_QUADRATUREAUDIT_H
#define RODIN_EXPERIMENTS_QUADRATUREAUDIT_H

#include <array>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "Rodin/Adaptation.h"
#include "Rodin/Assembly.h"
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"

namespace Rodin::Experiments
{
  // Fixed-state diagnostics only: no alternative solver or production policy.
  template <class Displacement>
  class QuadratureAudit
  {
    public:
      QuadratureAudit(const Displacement& current,
        const Adaptation::WNGIR::Parameters& parameters)
        : m_current(current), m_parameters(parameters)
      {}

      void volume() const
      {
        using namespace Variational;
        using namespace Geometry;
        const auto& fes = m_current.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const size_t dimension = mesh.getDimension();
        std::array<Displacement, 3> probes{
          Displacement(fes), Displacement(fes), Displacement(fes)};
        for (size_t i = 0; i < probes.size(); ++i)
          probes[i] = Adaptation::AnalyticVectorFunction([=](const Point& point) {
            Math::SpatialVector<Real> value = Math::SpatialVector<Real>::Zero(dimension);
            if (i < 2)
              value(i) = point.getCoordinates()(i) * point.getCoordinates()(i);
            else
              for (size_t axis = 0; axis < dimension; ++axis)
                value(axis) = point.x() * point.y();
            return value;
          }, dimension);
        auto gradient = Jacobian(m_current);
        std::array probeGradients{Jacobian(probes[0]), Jacobian(probes[1]), Jacobian(probes[2])};
        const bool affineP1 = [&] {
          for (auto cell = mesh.getCell(); cell; ++cell)
            if (fes.getFiniteElement(dimension, cell->getIndex()).getOrder() != 1 ||
                cell->getTransformation().getOrder() != 1)
              return false;
          return true;
        }();
        Math::Matrix<Real> reference;
        Math::Matrix<Real> referenceHinge;
        Real referenceHingeEnergy = 0;
        std::set<Index> critical;
        Real referenceMinJ = std::numeric_limits<Real>::infinity(), referenceMaxQ = 0;
        for (const size_t order : Orders)
        {
          const auto start = std::chrono::steady_clock::now();
          Math::Matrix<Real> action = Math::Matrix<Real>::Zero(3, 3);
          std::array<Math::SpatialMatrix<Real>, 3> means;
          for (auto& mean : means)
          {
            mean = Math::SpatialMatrix<Real>(dimension, dimension);
            mean.setZero();
          }
          Math::Vector<Real> traces = Math::Vector<Real>::Zero(3);
          Math::Matrix<Real> hinge = Math::Matrix<Real>::Zero(3, 3);
          Real hingeEnergy = 0;
          Real volume = 0, minJ = std::numeric_limits<Real>::infinity(), maxQ = 0;
          Index minCell = 0, maxCell = 0;
          Adaptation::CellDeformation deformation(dimension);
          const auto quality = [&](const auto& ip, Index index) {
            deformation.setDisplacementGradient(gradient.getValue(ip));
            const Real j = deformation.getJacobian();
            const Real q = deformation.isAdmissible()
              ? deformation.getRelativeDistortion() : std::numeric_limits<Real>::infinity();
            if (j < minJ) { minJ = j; minCell = index; }
            if (q > maxQ) { maxQ = q; maxCell = index; }
          };
          for (auto cell = mesh.getCell(); cell; ++cell)
          {
            const auto& qf = QF::PolytopeQuadratureFormula::get(
              affineP1 ? 2 : order, cell->getGeometry());
            const auto& quadrature = cell->getQuadrature(qf);
            for (size_t q = 0; q < quadrature.getSize(); ++q)
            {
              const IntegrationPoint ip(quadrature.getPoint(q), &qf, q);
              quality(ip, cell->getIndex());
              if (!deformation.isAdmissible())
                throw std::runtime_error("Fixed state is inverted; distribution is undefined.");
              const Real weight = qf.getWeight(q) * ip.getPoint().getDistortion() *
                deformation.getJacobian();
              // Frozen affine hinge, with a compressive trial that activates
              // the guard in the stressed state. Coefficient is fixed at one.
              Math::SpatialMatrix<Real> inner(dimension, dimension);
              inner.setZero();
              inner(0, 0) = -Real(0.2);
              const Adaptation::WNGIR::HingeState state(deformation, inner, m_parameters, 1);
              const Real referenceWeight = qf.getWeight(q) * ip.getPoint().getDistortion();
              hingeEnergy += referenceWeight * state.getEnergy(m_parameters, 1);
              for (size_t i = 0; i < 3; ++i)
                for (size_t j = 0; j < 3; ++j)
                {
                  const auto gi = probeGradients[i].getValue(ip);
                  const auto gj = probeGradients[j].getValue(ip);
                  hinge(i, j) += referenceWeight *
                    (state.getJacobianHessian() * state.getJacobianRow(gi) * state.getJacobianRow(gj) +
                     state.getDistortionHessian() * state.getDistortionRow(gi) * state.getDistortionRow(gj));
                }
              volume += weight;
              const Math::SpatialMatrix<Real> inverse =
                deformation.getInverseTranspose().transpose();
              std::array<Math::SpatialMatrix<Real>, 3> strains;
              Math::Vector<Real> trace(3);
              for (size_t i = 0; i < 3; ++i)
              {
                const Math::SpatialMatrix<Real> h = probeGradients[i].getValue(ip) * inverse;
                strains[i] = Real(0.5) * (h + h.transpose());
                trace(i) = strains[i].trace();
                for (size_t axis = 0; axis < dimension; ++axis)
                  strains[i](axis, axis) -= trace(i) / Real(dimension);
                means[i] += weight * strains[i];
                traces(i) += weight * trace(i);
              }
              for (size_t i = 0; i < 3; ++i)
                for (size_t j = 0; j < 3; ++j)
                  action(i, j) += weight *
                    (m_parameters.model.distribution.deviatoric *
                      (strains[i].transpose() * strains[j]).trace() +
                     m_parameters.model.distribution.divergence / Real(dimension) *
                       trace(i) * trace(j));
            }
            const Polytope::Traits traits(cell->getGeometry());
            for (size_t vertex = 0; vertex < traits.getVertexCount(); ++vertex)
              quality(Point(*cell, traits.getVertex(vertex)), cell->getIndex());
          }
          for (size_t i = 0; i < 3; ++i)
            for (size_t j = 0; j < 3; ++j)
              action(i, j) -= (m_parameters.model.distribution.deviatoric *
                (means[i].transpose() * means[j]).trace() +
                m_parameters.model.distribution.divergence / Real(dimension) *
                  traces(i) * traces(j)) / volume;
          action *= m_parameters.model.h;
          if (reference.size() == 0)
          {
            reference = action;
            referenceHinge = hinge;
            referenceHingeEnergy = hingeEnergy;
          }
          critical.insert(minCell);
          critical.insert(maxCell);
          referenceMinJ = std::min(referenceMinJ, minJ);
          referenceMaxQ = std::max(referenceMaxQ, maxQ);
          std::cout << "audit volume order=" << order
                    << " evaluated_order=" << (affineP1 ? 2 : order)
                    << " action_error=" << (action - reference).norm() /
                         std::max(reference.norm(), NumericalFloor)
                    << " hinge_error=" << (hinge - referenceHinge).norm() /
                         std::max(referenceHinge.norm(), NumericalFloor)
                    << " hinge_energy=" << hingeEnergy
                    << " hinge_energy_error=" << std::abs(hingeEnergy - referenceHingeEnergy) /
                         std::max(std::abs(referenceHingeEnergy), NumericalFloor)
                    << " min_j=" << minJ << " max_q=" << maxQ
                    << " seconds=" << elapsed(start) << '\n';
        }
        // All cells get a boundary-inclusive lattice at density 8; cells
        // limiting any rule also get nested 16/32/64 refinement. Not a proof.
        for (const size_t density : {8u, 16u, 32u, 64u})
        {
          const auto start = std::chrono::steady_clock::now();
          if (!affineP1)
            for (auto cell = mesh.getCell(); cell; ++cell)
              if (density == 8 || critical.contains(cell->getIndex()))
                lattice(*cell, density, [&](const Point& point) {
                  Adaptation::CellDeformation deformation(dimension);
                  deformation.setDisplacementGradient(gradient.getValue(point));
                  const Real j = deformation.getJacobian();
                  const Real q = deformation.isAdmissible()
                    ? deformation.getRelativeDistortion() : std::numeric_limits<Real>::infinity();
                  if (j < referenceMinJ || q > referenceMaxQ)
                    critical.insert(cell->getIndex());
                  referenceMinJ = std::min(referenceMinJ, j);
                  referenceMaxQ = std::max(referenceMaxQ, q);
                });
          std::cout << "audit quality density=" << density << " min_j=" << referenceMinJ
                    << " max_q=" << referenceMaxQ << " critical_cells=" << critical.size()
                    << " affine_p1=" << affineP1 << " seconds=" << elapsed(start) << '\n';
        }
      }

      template <class Phi, class Grad>
      void surface(const Phi& phi, const Grad& grad, Real sigma, Real normalization) const
      {
        using namespace Variational;
        const auto& fes = m_current.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const size_t dimension = mesh.getDimension();
        const Location::AABB<std::remove_cvref_t<decltype(mesh)>> locator(mesh);
        TrialFunction trial(fes);
        TestFunction test(fes);
        const Adaptation::WNGIR::Loss loss(sigma);
        const Adaptation::WNGIR::FittingTensor tensor(
          grad, m_current, locator, m_parameters, normalization, dimension);
        const Adaptation::WNGIR::FittingForce force(
          phi, grad, m_current, locator, loss, normalization, dimension);
        Math::SparseMatrix<Real> referenceMetric;
        Math::Vector<Real> referenceForce;
        Real referenceEnergy = 0, maximum = 0;
        for (const size_t order : SurfaceOrders)
        {
          const auto start = std::chrono::steady_clock::now();
          auto metricIntegral = FaceIntegral(Dot(tensor * trial, test));
          metricIntegral.setOrder(order);
          metricIntegral.over(*m_parameters.interfaceAttribute);
          auto forceIntegral = FaceIntegral(force, test);
          forceIntegral.setOrder(order);
          forceIntegral.over(*m_parameters.interfaceAttribute);
          BilinearForm metric(trial, test);
          LinearForm rhs(test);
          metric = metricIntegral;
          rhs = forceIntegral;
          metric.assemble();
          rhs.assemble();
          Adaptation::DeformationMap deformation(m_current, locator);
          Real energy = 0, sup = 0;
          const auto distance = [&](const Geometry::Point& point) {
            const auto& moved = deformation.getMovedPoint(IntegrationPoint(point));
            const Real value = std::abs(phi.getValue(moved)) / grad.getValue(moved).norm();
            sup = std::max(sup, value);
          };
          for (auto face = mesh.getFace(); face; ++face)
            if (face->getAttribute() == *m_parameters.interfaceAttribute)
            {
              const auto& qf = QF::PolytopeQuadratureFormula::get(order, face->getGeometry());
              const auto& quadrature = face->getQuadrature(qf);
              for (size_t q = 0; q < quadrature.getSize(); ++q)
              {
                const auto& point = quadrature.getPoint(q);
                const IntegrationPoint ip(point, &qf, q);
                const auto& moved = deformation.getMovedPoint(ip);
                energy += normalization * qf.getWeight(q) * point.getDistortion() *
                  loss.getValue(phi.getValue(moved));
                distance(point);
              }
              const Geometry::Polytope::Traits traits(face->getGeometry());
              for (size_t vertex = 0; vertex < traits.getVertexCount(); ++vertex)
                distance(Geometry::Point(*face, traits.getVertex(vertex)));
            }
          if (referenceForce.size() == 0)
          {
            referenceMetric = metric.getOperator();
            referenceForce = rhs.getVector();
            referenceEnergy = energy;
          }
          maximum = std::max(maximum, sup);
          std::cout << "audit surface order=" << order << " energy=" << energy
                    << " energy_error=" << std::abs(energy - referenceEnergy) /
                         std::max(std::abs(referenceEnergy), NumericalFloor)
                    << " force_error=" << (rhs.getVector() - referenceForce).norm() /
                         std::max(referenceForce.norm(), NumericalFloor)
                    << " metric_error=" << (metric.getOperator() - referenceMetric).norm() /
                         std::max(referenceMetric.norm(), NumericalFloor)
                    << " geom_sup=" << sup << " seconds=" << elapsed(start) << '\n';
        }
        const Adaptation::DeformationMap deformation(m_current, locator);
        for (const size_t density : {8u, 16u, 32u, 64u})
        {
          for (auto face = mesh.getFace(); face; ++face)
            if (face->getAttribute() == *m_parameters.interfaceAttribute)
              lattice(*face, density, [&](const Geometry::Point& point) {
                const auto& moved = deformation.getMovedPoint(IntegrationPoint(point));
                maximum = std::max(maximum,
                  std::abs(phi.getValue(moved)) / grad.getValue(moved).norm());
              });
          std::cout << "audit geometry density=" << density << " geom_sup=" << maximum << '\n';
        }
      }

    private:
      template <class Evaluate>
      static void lattice(const Geometry::Polytope& cell, size_t density, Evaluate evaluate)
      {
        const Geometry::Polytope::Traits traits(cell.getGeometry());
        if (traits.getVertexCount() != cell.getDimension() + 1)
          throw std::runtime_error("Calibration lattice requires a simplex.");
        auto reference = Math::SpatialPoint::Zero(cell.getDimension());
        const auto visit = [&](auto&& self, size_t axis, size_t remaining) -> void {
          if (axis == cell.getDimension())
          {
            evaluate(Geometry::Point(cell, reference));
            return;
          }
          for (size_t i = 0; i <= remaining; ++i)
          {
            reference(axis) = Real(i) / Real(density);
            self(self, axis + 1, remaining - i);
          }
        };
        visit(visit, 0, density);
      }

      static Real elapsed(std::chrono::steady_clock::time_point start)
      {
        return std::chrono::duration<Real>(std::chrono::steady_clock::now() - start).count();
      }
      static constexpr Real NumericalFloor = Real(1e-14);
      static constexpr std::array<size_t, 8> Orders{32, 24, 2, 4, 6, 8, 12, 16};
      static constexpr std::array<size_t, 9> SurfaceOrders{64, 32, 2, 4, 6, 8, 12, 16, 24};
      const Displacement& m_current;
      const Adaptation::WNGIR::Parameters& m_parameters;
  };
}

#endif
