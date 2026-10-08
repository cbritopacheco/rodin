/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_SWIFT_HINGEFORCE_H
#define RODIN_ADAPTATION_SWIFT_HINGEFORCE_H

#include "../CellDeformation.h"
#include "Hinge.h"

namespace Rodin::Adaptation::SWIFT
{
  /// @brief Hinge contribution to the next primal Newton iterate.
  template <class TestFunction, class Displacement>
  class HingeForce final
    : public Variational::LinearFormIntegratorBase<typename TestFunction::ScalarType>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = typename TestFunction::ScalarType;
      /// @brief Parent class type.
      using Parent = Variational::LinearFormIntegratorBase<ScalarType>;

      /**
       * @brief Constructs the affine quadratic hinge load.
       * @param z Test function for the displacement increment.
       * @param current Frozen outer displacement.
       * @param inner Current inner displacement increment.
       * @param parameters Quality budgets and hinge activation weights.
       * @param hingeCoefficient Effective hinge coefficient for this outer model.
       */
      HingeForce(const TestFunction& z, const Displacement& current,
        const Displacement& inner, const Parameters& parameters, Real hingeCoefficient)
        : Parent(z.getLeaf()),
          m_z(z),
          m_current(current),
          m_inner(inner),
          m_parameters(parameters),
          m_hingeCoefficient(hingeCoefficient)
      {}

      /**
       * @brief Copy constructor.
       * @param other Integrator to copy, retaining its field references.
       */
      HingeForce(const HingeForce& other) = default;

      /**
       * @brief Returns the current polytope.
       * @returns The current polytope.
       */
      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      /**
       * @brief Binds to a polytope and assembles the local load.
       * @param polytope Cell to integrate.
       * @returns This integrator after binding and assembly.
       */
      HingeForce& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;
        const std::size_t dim = polytope.getDimension();
        const Index index = polytope.getIndex();
        const auto& fes = m_z.getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(dim, index);
        const auto& parameters = m_parameters.get();
        const std::size_t order = parameters.quadrature.getVolumeOrder(fe.getOrder(),
          polytope.getTransformation().getOrder(),
          Geometry::Polytope::Traits(polytope.getGeometry()).getVertexCount() == dim + 1);
        const auto& qf =
          QF::PolytopeQuadratureFormula::get(order, polytope.getGeometry());
        const auto& quadrature = polytope.getQuadrature(qf);

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
          const HingeState state(
            deformation, innerJacobian.getValue(ip), parameters, m_hingeCoefficient);
          const Real coefficientJ = state.getJacobianForce();
          const Real coefficientQ = state.getDistortionForce();
          const Real weight = qf.getWeight(q) * point.getDistortion();
          testJacobian.setIntegrationPoint(ip);
          for (std::size_t local = 0; local < fe.getCount(); ++local)
          {
            const auto gradient = testJacobian.getBasis(local);
            const Real rowJ = state.getJacobianRow(gradient);
            const Real rowQ = state.getDistortionRow(gradient);
            m_vector(static_cast<Eigen::Index>(local)) +=
              weight * (coefficientJ * rowJ + coefficientQ * rowQ);
          }
        }
        return *this;
      }

      /**
       * @brief Returns an entry of the assembled local system.
       * @param local Index in the local numbering.
       * @returns Integral computed by the quadrature rule.
       */
      ScalarType integrate(std::size_t local) final override
      {
        return m_vector(static_cast<Eigen::Index>(local));
      }

      /**
       * @brief Returns the integration region.
       * @returns The integration region.
       */
      Geometry::Region getRegion() const final override
      {
        return Geometry::Region::Cells;
      }

      /**
       * @brief Clones this integrator.
       * @returns Newly allocated copy owned by the caller.
       */
      HingeForce* copy() const noexcept final override
      {
        return new HingeForce(*this);
      }

    private:
      TestFunction m_z;
      std::reference_wrapper<const Displacement> m_current;
      std::reference_wrapper<const Displacement> m_inner;
      std::reference_wrapper<const Parameters> m_parameters;
      Real m_hingeCoefficient;
      const Geometry::Polytope* m_polytope = nullptr;
      Math::Vector<Real> m_vector;
  };
}

#endif
