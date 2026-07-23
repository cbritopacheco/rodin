/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_WNGIRMARGINFORCE_H
#define RODIN_ADAPTATION_WNGIRMARGINFORCE_H

#include "CellDeformation.h"
#include "WNGIRConstraintState.h"
#include "WNGIRParameters.h"

namespace Rodin::Adaptation::Detail
{
  /**
   * @brief Affine load of the fractional-margin metric model.
   *
   * For a linearized row @f$C v\leq b@f$, the local penalty
   * @f$\frac{\gamma}{2}(Cv-\theta b)^2@f$ contributes
   * @f$C^\top\gamma C@f$ to the metric and
   * @f$C^\top\gamma\theta b@f$ to its right-hand side.
   */
  template <class TestFunction, class Displacement>
  class WNGIRMarginForce final
    : public Variational::LinearFormIntegratorBase<typename TestFunction::ScalarType>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = typename TestFunction::ScalarType;
      /// @brief Parent class type.
      using Parent = Variational::LinearFormIntegratorBase<ScalarType>;

      /// @brief Constructs the w n g i r margin force.
      WNGIRMarginForce(const TestFunction& v, const Displacement& current,
        const Displacement& predictor, const WNGIRParameters& parameters,
        Real constraintMultiplier = Real(1))
        : Parent(v.getLeaf()),
          m_v(v),
          m_current(current),
          m_predictor(predictor),
          m_parameters(parameters),
          m_constraintMultiplier(constraintMultiplier)
      {}

      /// @brief Copy constructor.
      WNGIRMarginForce(const WNGIRMarginForce&) = default;

      /// @brief Returns the current polytope.
      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      /// @brief Binds to a polytope and assembles the local system.
      WNGIRMarginForce& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;
        const std::size_t dim = polytope.getDimension();
        const auto idx = polytope.getIndex();
        const auto& fes = m_v.getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(dim, idx);
        const auto& params = m_parameters.get();
        const std::size_t qOrder = params.quadratureOrder > 0
          ? params.quadratureOrder
          : std::max<std::size_t>(2, 2 * fe.getOrder());
        const auto& qf =
          QF::PolytopeQuadratureFormula::get(qOrder, polytope.getGeometry());
        const auto& quad = polytope.getQuadrature(qf);

        m_vector.resize(static_cast<Eigen::Index>(fe.getCount()));
        m_vector.setZero();
        auto testJacobian = Variational::Jacobian(m_v);
        auto currentJacobian = Variational::Jacobian(m_current.get());
        auto predictorJacobian = Variational::Jacobian(m_predictor.get());
        CellDeformation deformation(dim);

        for (std::size_t qp = 0; qp < quad.getSize(); ++qp)
        {
          const auto& pt = quad.getPoint(qp);
          const Variational::IntegrationPoint ip(pt, &qf, qp);
          deformation.setDisplacementGradient(currentJacobian.getValue(ip));
          if (!deformation.isInvertible())
            continue;
          const auto predictorGradient = predictorJacobian.getValue(ip);
          const WNGIRConstraintState state(deformation, params, &predictorGradient);
          const Real jWeight = state.getJacobianWeight();
          const Real qWeight = state.getDistortionWeight();
          if (jWeight <= Real(0) && qWeight <= Real(0))
            continue;
          const Real w = qf.getWeight(qp) * pt.getDistortion();
          testJacobian.setIntegrationPoint(ip);
          for (std::size_t local = 0; local < fe.getCount(); ++local)
          {
            const auto basisGradient = testJacobian.getBasis(local);
            const Real cJ = -deformation.getJacobianAction(basisGradient);
            const Real cQ = qWeight > Real(0)
              ? deformation.getRelativeDistortionAction(basisGradient)
              : Real(0);
            m_vector(static_cast<Eigen::Index>(local)) += w * m_constraintMultiplier *
              (params.gammaJ * jWeight * state.getJacobianTarget() * cJ +
                params.gammaQ * qWeight * state.getDistortionTarget() * cQ);
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
      WNGIRMarginForce* copy() const noexcept final override
      {
        return new WNGIRMarginForce(*this);
      }

    private:
      TestFunction m_v;
      std::reference_wrapper<const Displacement> m_current;
      std::reference_wrapper<const Displacement> m_predictor;
      std::reference_wrapper<const WNGIRParameters> m_parameters;
      Real m_constraintMultiplier;
      const Geometry::Polytope* m_polytope = nullptr;
      Math::Vector<Real> m_vector;
  };
}

#endif
