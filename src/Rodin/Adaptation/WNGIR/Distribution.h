/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_WNGIR_DISTRIBUTION_H
#define RODIN_ADAPTATION_WNGIR_DISTRIBUTION_H

#include <vector>
#include "Rodin/Assembly.h"
#include "Rodin/Variational.h"
#include "../CellDeformation.h"

namespace Rodin::Adaptation::WNGIR
{
  template <class Displacement>
  class CurrentInverse;

  /**
   * @brief Sparse core of the frozen centered current-strain bilinear form.
   *
   * Binds cached geometry and basis Jacobians, evaluates j and F^-1 once per
   * quadrature point, then tabulates sym(grad v F^-1) before contracting pairs.
   * Computes the local weighted deviatoric and divergence products. Subtracting
   * the global mean-strain couplings gives the centered form; no independent
   * element means are subtracted. Assembly clones isolate bind state.
   */
  template <class TrialFunction, class TestFunction, class Displacement>
  class Distribution final : public Variational::LocalBilinearFormIntegratorBase<
                                    typename TrialFunction::ScalarType>
  {
    public:
      using ScalarType = typename TrialFunction::ScalarType;
      using Parent = Variational::LocalBilinearFormIntegratorBase<ScalarType>;

      Distribution(const TrialFunction& trial, const TestFunction& test,
        const Displacement& current, Real deviatoric, Real divergence, size_t order)
        : Parent(trial.getLeaf(), test.getLeaf()),
          m_trial(trial),
          m_test(test),
          m_current(current),
          m_deviatoric(deviatoric),
          m_divergence(divergence),
          m_order(order)
      {}

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      Distribution& setPolytope(const Geometry::Polytope& cell) final override
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
          const Math::SpatialMatrix<Real> inverse(
            deformation.getInverseTranspose().transpose());
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
          const Real weight = qf.getWeight(q) * point.getDistortion() *
            deformation.getJacobian();
          for (size_t test = 0; test < nTest; ++test)
            for (size_t trial = 0; trial < nTrial; ++trial)
              m_matrix(test, trial) +=
                weight * (m_deviatoric *
                  Math::dot(m_trialStrains[trial], m_testStrains[test]) +
                  (m_divergence - m_deviatoric) / Real(d) *
                    m_trialStrains[trial].trace() * m_testStrains[test].trace());
        }
        return *this;
      }

      ScalarType integrate(size_t trial, size_t test) final override
      {
        return m_matrix(test, trial);
      }
      Geometry::Region getRegion() const final override
      {
        return Geometry::Region::Cells;
      }
      Distribution* copy() const noexcept final override
      {
        return new Distribution(*this);
      }

      /**
       * @brief Global mean-strain subtraction represented by @f$UU^\top@f$.
       * Native linear forms integrate the symmetric strain in an orthonormal
       * tensor basis. The volume and moments use the sparse core's quadrature,
       * preserving the variance identity for the centered form.
       * @returns Mean-strain coupling matrix @f$U@f$ for the current fields and weights.
       */
      Math::Matrix<Real> getCentering() const
      {
        const auto& test = m_test.get();
        const auto& current = m_current.get();
        const auto dimension = test.getFiniteElementSpace().getMesh().getDimension();
        const auto weight = Variational::RealFunction(
          [gradient = Variational::Jacobian(current), dimension](const auto& point) -> Real {
            CellDeformation deformation(dimension);
            deformation.setDisplacementGradient(gradient.getValue(point));
            return deformation.getJacobian();
          });
        const auto gradient = Variational::Jacobian(test) *
          CurrentInverse<Displacement>(current, dimension);
        const auto strain = Real(0.5) * (gradient + Variational::Transpose(gradient));
        auto traceIntegral = Variational::Integral(weight * Variational::Trace(strain));
        traceIntegral.setOrder(m_order);
        Variational::LinearForm trace(test);
        trace = traceIntegral;
        trace.assemble();
        Variational::GridFunction position(test.getFiniteElementSpace());
        position = Variational::VectorFunction(dimension, [](const Geometry::Point& point) {
          return Math::SpatialVector<Real>(point.getCoordinates());
        });
        position += current;
        const Real volume = trace.getVector().dot(position.getData()) / Real(dimension);
        assert(volume > Real(0));
        Math::Matrix<Real> centering(trace.getVector().size(), dimension * (dimension + 1) / 2);
        size_t column = 0;
        for (size_t a = 0; a < dimension; ++a)
          for (size_t b = a; b < dimension; ++b)
          {
            Math::Matrix<Real> tensor = Math::Matrix<Real>::Zero(dimension, dimension);
            if (a == b)
              tensor(a, b) = Real(1);
            else
              tensor(a, b) = tensor(b, a) = Real(1) / std::sqrt(Real(2));
            const Real mean = tensor.trace() / Real(dimension);
            tensor *= std::sqrt(m_deviatoric);
            for (size_t axis = 0; axis < dimension; ++axis)
              tensor(axis, axis) +=
                (std::sqrt(m_divergence) - std::sqrt(m_deviatoric)) * mean;
            auto integral = Variational::Integral(weight *
              Variational::Dot(Variational::MatrixFunction(tensor), strain));
            integral.setOrder(m_order);
            Variational::LinearForm form(test);
            form = integral;
            form.assemble();
            centering.col(column++) = form.getVector() / std::sqrt(volume);
          }
        return centering;
      }

    private:
      std::reference_wrapper<const TrialFunction> m_trial;
      std::reference_wrapper<const TestFunction> m_test;
      std::reference_wrapper<const Displacement> m_current;
      Real m_deviatoric;
      Real m_divergence;
      size_t m_order;
      const Geometry::Polytope* m_polytope = nullptr;
      Math::Matrix<Real> m_matrix;
      std::vector<Math::SpatialMatrix<Real>> m_trialStrains, m_testStrains;
  };

  /// Frozen inverse deformation gradient, expressed as a form-language coefficient.
  template <class Displacement>
  class CurrentInverse final
    : public Variational::MatrixFunctionBase<Real, CurrentInverse<Displacement>>
  {
    public:
      CurrentInverse(const Displacement& current, size_t dimension)
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
      CurrentInverse* copy() const noexcept override
      {
        return new CurrentInverse(*this);
      }

    private:
      decltype(Variational::Jacobian(std::declval<const Displacement&>())) m_gradient;
      size_t m_dimension;
  };

}
#endif
