/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_WNGIRREGULARITYMETRIC_H
#define RODIN_ADAPTATION_WNGIRREGULARITYMETRIC_H

#include <cmath>
#include <limits>
#include <vector>
#include "Rodin/Assembly.h"
#include "Rodin/Variational.h"
#include "CellDeformation.h"

namespace Rodin::Adaptation::Detail
{
  /**
   * @brief Frozen current-strain bilinear form, excluding the global dilation correction.
   *
   * Binds cached geometry and basis Jacobians, evaluates j and F^-1 once per
   * quadrature point, then tabulates sym(grad v F^-1) before contracting pairs.
   * This is the same integral as wngirCurrentStrain, without repeated coefficient
   * evaluation inside the basis-pair loop. Assembly clones isolate bind state.
   */
  template <class TrialFunction, class TestFunction, class Displacement>
  class WNGIRCurrentStrainMetric final
    : public Variational::LocalBilinearFormIntegratorBase<typename TrialFunction::ScalarType>
  {
    public:
      using ScalarType = typename TrialFunction::ScalarType;
      using Parent = Variational::LocalBilinearFormIntegratorBase<ScalarType>;

      WNGIRCurrentStrainMetric(const TrialFunction& trial, const TestFunction& test,
        const Displacement& current, Real coefficient, size_t order)
        : Parent(trial.getLeaf(), test.getLeaf()), m_trial(trial), m_test(test),
          m_current(current), m_coefficient(coefficient), m_order(order)
      {}

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      WNGIRCurrentStrainMetric& setPolytope(const Geometry::Polytope& cell) final override
      {
        m_polytope = &cell;
        const auto d = cell.getDimension();
        const auto& trialFES = m_trial.get().getFiniteElementSpace();
        const auto& testFES = m_test.get().getFiniteElementSpace();
        const auto nTrial = trialFES.getFiniteElement(d, cell.getIndex()).getCount();
        const auto nTest = testFES.getFiniteElement(d, cell.getIndex()).getCount();
        const auto& qf = QF::PolytopeQuadratureFormula::get(m_order, cell.getGeometry());
        const auto& quadrature = cell.getQuadrature(qf);
        m_matrix = Math::Matrix<Real>::Zero(nTest, nTrial);
        m_trialStrains.resize(nTrial);
        m_testStrains.resize(nTest);
        auto currentJacobian = Variational::Jacobian(m_current.get());
        auto trialJacobian = Variational::Jacobian(m_trial.get());
        auto testJacobian = Variational::Jacobian(m_test.get());
        for (size_t q = 0; q < quadrature.getSize(); ++q)
        {
          const auto& point = quadrature.getPoint(q);
          const Variational::IntegrationPoint ip(point, &qf, q);
          CellDeformation deformation(d);
          deformation.setDisplacementGradient(currentJacobian.getValue(ip));
          const Math::SpatialMatrix<Real> inverse(deformation.getInverseTranspose().transpose());
          trialJacobian.setIntegrationPoint(ip);
          for (size_t local = 0; local < nTrial; ++local)
          {
            const Math::SpatialMatrix<Real> L(trialJacobian.getBasis(local) * inverse);
            m_trialStrains[local] = Real(0.5) * (L + L.transpose());
          }
          if (trialFES == testFES)
            m_testStrains = m_trialStrains;
          else
          {
            testJacobian.setIntegrationPoint(ip);
            for (size_t local = 0; local < nTest; ++local)
            {
              const Math::SpatialMatrix<Real> L(testJacobian.getBasis(local) * inverse);
              m_testStrains[local] = Real(0.5) * (L + L.transpose());
            }
          }
          const Real weight = m_coefficient * qf.getWeight(q) * point.getDistortion() *
            deformation.getJacobian();
          for (size_t test = 0; test < nTest; ++test)
            for (size_t trial = 0; trial < nTrial; ++trial)
              m_matrix(test, trial) += weight *
                Math::dot(m_trialStrains[trial], m_testStrains[test]);
        }
        return *this;
      }

      ScalarType integrate(size_t trial, size_t test) final override
      {
        return m_matrix(test, trial);
      }
      Geometry::Region getRegion() const final override { return Geometry::Region::Cells; }
      WNGIRCurrentStrainMetric* copy() const noexcept final override
      {
        return new WNGIRCurrentStrainMetric(*this);
      }

    private:
      std::reference_wrapper<const TrialFunction> m_trial;
      std::reference_wrapper<const TestFunction> m_test;
      std::reference_wrapper<const Displacement> m_current;
      Real m_coefficient;
      size_t m_order;
      const Geometry::Polytope* m_polytope = nullptr;
      Math::Matrix<Real> m_matrix;
      std::vector<Math::SpatialMatrix<Real>> m_trialStrains, m_testStrains;
  };

  /// Frozen linear trace coupling, using one coefficient evaluation per quadrature point.
  template <class TestFunction, class Displacement>
  class WNGIRCurrentStrainTrace final
    : public Variational::LinearFormIntegratorBase<typename TestFunction::ScalarType>
  {
    public:
      using ScalarType = typename TestFunction::ScalarType;
      using Parent = Variational::LinearFormIntegratorBase<ScalarType>;
      WNGIRCurrentStrainTrace(const TestFunction& test, const Displacement& current, size_t order)
        : Parent(test.getLeaf()), m_test(test), m_current(current), m_order(order)
      {}
      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }
      WNGIRCurrentStrainTrace& setPolytope(const Geometry::Polytope& cell) final override
      {
        m_polytope = &cell;
        const auto d = cell.getDimension();
        const auto n = m_test.get().getFiniteElementSpace().getFiniteElement(d, cell.getIndex()).getCount();
        const auto& qf = QF::PolytopeQuadratureFormula::get(m_order, cell.getGeometry());
        const auto& quadrature = cell.getQuadrature(qf);
        m_vector = Math::Vector<Real>::Zero(n);
        auto currentJacobian = Variational::Jacobian(m_current.get());
        auto testJacobian = Variational::Jacobian(m_test.get());
        for (size_t q = 0; q < quadrature.getSize(); ++q)
        {
          const auto& point = quadrature.getPoint(q);
          const Variational::IntegrationPoint ip(point, &qf, q);
          CellDeformation deformation(d);
          deformation.setDisplacementGradient(currentJacobian.getValue(ip));
          const Math::SpatialMatrix<Real> inverse(deformation.getInverseTranspose().transpose());
          const Real weight = qf.getWeight(q) * point.getDistortion() *
            deformation.getJacobian() / std::sqrt(Real(d));
          testJacobian.setIntegrationPoint(ip);
          for (size_t local = 0; local < n; ++local)
            m_vector(local) += weight * (testJacobian.getBasis(local) * inverse).trace();
        }
        return *this;
      }
      ScalarType integrate(size_t local) final override { return m_vector(local); }
      Geometry::Region getRegion() const final override { return Geometry::Region::Cells; }
      WNGIRCurrentStrainTrace* copy() const noexcept final override
      {
        return new WNGIRCurrentStrainTrace(*this);
      }
    private:
      std::reference_wrapper<const TestFunction> m_test;
      std::reference_wrapper<const Displacement> m_current;
      size_t m_order;
      const Geometry::Polytope* m_polytope = nullptr;
      Math::Vector<Real> m_vector;
  };
  /// Frozen inverse deformation gradient, expressed as a form-language coefficient.
  template <class Displacement>
  class WNGIRCurrentInverse final
    : public Variational::MatrixFunctionBase<Real, WNGIRCurrentInverse<Displacement>>
  {
    public:
      WNGIRCurrentInverse(const Displacement& current, size_t dimension)
        : m_gradient(Variational::Jacobian(current)),
          m_dimension(dimension)
      {}
      template <class Point>
      Math::SpatialMatrix<Real> getValue(const Point& point) const
      {
        CellDeformation deformation(m_dimension);
        deformation.setDisplacementGradient(m_gradient.getValue(point));
        return Math::SpatialMatrix<Real>(deformation.getInverseTranspose().transpose());
      }
      size_t getRows() const
      {
        return m_dimension;
      }
      size_t getColumns() const
      {
        return m_dimension;
      }
      Optional<size_t> getOrder(const Geometry::Polytope&) const
      {
        return std::nullopt;
      }
      WNGIRCurrentInverse* copy() const noexcept override
      {
        return new WNGIRCurrentInverse(*this);
      }

    private:
      decltype(Variational::Jacobian(std::declval<const Displacement&>())) m_gradient;
      size_t m_dimension;
  };

  template <class Displacement>
  auto wngirCurrentVolumeWeight(const Displacement& current, size_t dimension)
  {
    return Variational::RealFunction(
      [gradient = Variational::Jacobian(current), dimension](const auto& point) -> Real {
        CellDeformation deformation(dimension);
        deformation.setDisplacementGradient(gradient.getValue(point));
        return deformation.getJacobian();
      });
  }

  /// Pullback of sym(grad_y v), with y=x+u_k(x).
  template <class Function, class Displacement>
  auto wngirCurrentStrain(
    const Function& function, const Displacement& current, size_t dimension)
  {
    const auto gradient =
      Variational::Jacobian(function) * WNGIRCurrentInverse(current, dimension);
    return Real(0.5) * (gradient + Variational::Transpose(gradient));
  }

  /// Distribution core = integral j*eps(v):eps(z) - b(v)b(z)/(d*integral j).
  /// Only uniform current-configuration dilation is projected out.
  template <class TestFunction, class Displacement>
  std::vector<Math::Vector<Real>> wngirCurrentStrainCouplings(const TestFunction& test,
    const Displacement& current, size_t dimension, Real coefficient, size_t order)
  {
    if (coefficient == Real(0))
      return {};
    if (!(coefficient > Real(0)))
      Alert::Exception() << "Current strain regularity requires positive coefficient."
                         << Alert::Raise;
    Variational::LinearForm form(test);
    form = WNGIRCurrentStrainTrace(test, current, order);
    form.assemble();
    Variational::GridFunction dilation(test.getFiniteElementSpace());
    dilation = Variational::VectorFunction(dimension, [](const Geometry::Point& point) {
      return Math::SpatialVector<Real>(point.getCoordinates());
    });
    dilation += current;
    dilation *= Real(1) / std::sqrt(Real(dimension));
    const Real measure = form.getVector().dot(dilation.getData());
    if (!(measure > Real(0)) || !std::isfinite(measure))
      Alert::Exception()
        << "Current strain regularity requires positive finite current volume."
        << Alert::Raise;
    return {Math::Vector<Real>(std::sqrt(coefficient / measure) * form.getVector())};
  }

}
#endif
