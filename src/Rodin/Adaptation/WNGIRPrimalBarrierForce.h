/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_WNGIRPRIMALBARRIERFORCE_H
#define RODIN_ADAPTATION_WNGIRPRIMALBARRIERFORCE_H

#include "CellDeformation.h"
#include "WNGIRPrimalBarrierState.h"
#include "WNGIRValidationWeights.h"

namespace Rodin::Adaptation::Detail
{
  /// @brief Barrier contribution to the next primal Newton iterate.
  template <class TestFunction, class Displacement>
  class WNGIRPrimalBarrierForce final
    : public Variational::LinearFormIntegratorBase<typename TestFunction::ScalarType>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = typename TestFunction::ScalarType;
      /// @brief Parent class type.
      using Parent = Variational::LinearFormIntegratorBase<ScalarType>;

      /// @brief Constructs the w n g i r primal barrier force.
      WNGIRPrimalBarrierForce(const TestFunction& z, const Displacement& current,
        const Displacement& inner, const WNGIRParameters& parameters,
        Real barrierCoefficient)
        : Parent(z.getLeaf()),
          m_z(z),
          m_current(current),
          m_inner(inner),
          m_parameters(parameters),
          m_barrierCoefficient(barrierCoefficient)
      {}

      /// @brief Copy constructor.
      WNGIRPrimalBarrierForce(const WNGIRPrimalBarrierForce&) = default;

      /// @brief Returns the current polytope.
      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      /// @brief Binds to a polytope and assembles the local system.
      WNGIRPrimalBarrierForce& setPolytope(
        const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;
        const std::size_t dim = polytope.getDimension();
        const Index index = polytope.getIndex();
        const auto& fes = m_z.getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(dim, index);
        const auto& parameters = m_parameters.get();
        const std::size_t order = parameters.quadratureOrder > 0
          ? parameters.quadratureOrder
          : std::max<std::size_t>(2, 2 * fe.getOrder());
        const auto& qf =
          QF::PolytopeQuadratureFormula::get(order, polytope.getGeometry());
        const auto& quadrature = polytope.getQuadrature(qf);
        const WNGIRValidationWeights validationWeights(qf);

        m_vector.resize(static_cast<Eigen::Index>(fe.getCount()));
        m_vector.setZero();
        auto testJacobian = Variational::Jacobian(m_z);
        auto currentJacobian = Variational::Jacobian(m_current.get());
        auto innerJacobian = Variational::Jacobian(m_inner.get());
        CellDeformation deformation(dim);

        for (std::size_t q = 0; q < quadrature.getSize(); ++q)
        {
          const auto& point = quadrature.getPoint(q);
          const Variational::IntegrationPoint ip(point, &qf, q);
          deformation.setDisplacementGradient(currentJacobian.getValue(ip));
          if (!deformation.isAdmissible())
            continue;
          const WNGIRPrimalBarrierState state(
            deformation, innerJacobian.getValue(ip), parameters, m_barrierCoefficient);
          assert(state.isFeasible());
          if (!state.isFeasible())
            continue;

          const Real coefficientJ = state.getJacobianForce();
          const Real coefficientQ = state.getDistortionForce();
          const Real weight = validationWeights.getWeight(q) * point.getDistortion();
          testJacobian.setIntegrationPoint(ip);
          for (std::size_t local = 0; local < fe.getCount(); ++local)
          {
            const auto gradient = testJacobian.getBasis(local);
            const Real rowJ = -deformation.getJacobianAction(gradient);
            const Real rowQ = deformation.getRelativeDistortionAction(gradient);
            m_vector(static_cast<Eigen::Index>(local)) +=
              weight * (coefficientJ * rowJ + coefficientQ * rowQ);
          }
        }
        return *this;
      }

      /// @brief Returns an entry of the assembled local system.
      ScalarType integrate(std::size_t local) final override
      {
        return m_vector(static_cast<Eigen::Index>(local));
      }

      /// @brief Returns the integration region.
      Geometry::Region getRegion() const final override
      {
        return Geometry::Region::Cells;
      }

      /// @brief Clones this object.
      WNGIRPrimalBarrierForce* copy() const noexcept final override
      {
        return new WNGIRPrimalBarrierForce(*this);
      }

    private:
      TestFunction m_z;
      std::reference_wrapper<const Displacement> m_current;
      std::reference_wrapper<const Displacement> m_inner;
      std::reference_wrapper<const WNGIRParameters> m_parameters;
      Real m_barrierCoefficient;
      const Geometry::Polytope* m_polytope = nullptr;
      Math::Vector<Real> m_vector;
  };
}

#endif
