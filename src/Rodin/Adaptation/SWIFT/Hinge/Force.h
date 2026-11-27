/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_SWIFT_HINGE_FORCE_H
#define RODIN_ADAPTATION_SWIFT_HINGE_FORCE_H

#include "Rodin/Variational/Integrator.h"
#include "Model.h"

namespace Rodin::Adaptation::SWIFT::Hinge
{
  /// @brief Hinge contribution to the next primal Newton iterate.
  template <class TestFunction, class Displacement>
  class Force final
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
       * @param model Frozen outer hinge model, shared with merit evaluation.
       * @param inner Increment at which to evaluate the hinges.
       */
      Force(const TestFunction& z, const Model<Displacement>& model,
        const Displacement& inner)
        : Parent(z.getLeaf()),
          m_z(z),
          m_model(model),
          m_inner(inner)
      {}

      /**
       * @brief Copy constructor.
       * @param other Integrator to copy, retaining its field references.
       */
      Force(const Force& other) = default;

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
      Force& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;
        const std::size_t dim = polytope.getDimension();
        const Index index = polytope.getIndex();
        const auto& fes = m_z.getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(dim, index);
        const auto& parameters = m_model.getParameters();
        const auto hinge =
          m_model.on(polytope, fe.getDegree().value_or(2 * fe.getOrder()));

        m_vector.resize(static_cast<Eigen::Index>(fe.getCount()));
        m_vector.setZero();
        auto testJacobian = Variational::Jacobian(m_z);
        testJacobian.setPolytope(polytope);
        auto currentJacobian = Variational::Jacobian(m_model.getDisplacement());
        auto innerJacobian = Variational::Jacobian(m_inner.get());
        CellDeformation deformation(dim);

        const auto& witnesses = hinge.getWitnesses();
        const auto& weights = hinge.getWeights();
        for (size_t q = 0; q < witnesses.getSize(); ++q)
        {
          const auto ip = witnesses.getPoint(q);
          deformation.setDisplacementGradient(currentJacobian.getValue(ip));
          if (!deformation.isAdmissible())
            continue;
          const State state(deformation, innerJacobian.getValue(ip), parameters,
            m_model.getCoefficient());
          const Real coefficientJ = state.getJacobianForce();
          const Real coefficientQ = state.getDistortionForce();
          if (coefficientJ == Real(0) && coefficientQ == Real(0))
            continue;
          const auto [weightJ, weightQ] = weights[q];
          testJacobian.setIntegrationPoint(ip);
          for (std::size_t local = 0; local < fe.getCount(); ++local)
          {
            const auto gradient = testJacobian.getBasis(local);
            const Real rowJ = state.getJacobianRow(gradient);
            const Real rowQ = state.getDistortionRow(gradient);
            m_vector(static_cast<Eigen::Index>(local)) +=
              weightJ * coefficientJ * rowJ + weightQ * coefficientQ * rowQ;
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
      Force* copy() const noexcept final override
      {
        return new Force(*this);
      }

    private:
      TestFunction m_z;
      Model<Displacement> m_model;
      std::reference_wrapper<const Displacement> m_inner;
      const Geometry::Polytope* m_polytope = nullptr;
      Math::Vector<Real> m_vector;
  };
}

#endif
