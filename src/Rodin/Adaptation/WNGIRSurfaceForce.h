/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file WNGIRSurfaceForce.h
 * @brief Linear-form integrator assembling the Welsch first-variation
 * surface force on interface facets (the WNGIR right-hand side).
 */
#ifndef RODIN_ADAPTATION_WNGIRSURFACEFORCE_H
#define RODIN_ADAPTATION_WNGIRSURFACEFORCE_H

#include "WNGIRDetail.h"
#include "WNGIRParameters.h"

namespace Rodin::Adaptation::Detail
{
    template <class PhiDerived, class GradDerived,
              class TestFunction, class Displacement>
    class WNGIRSurfaceForce final
      : public Variational::LinearFormIntegratorBase<
          typename TestFunction::ScalarType>
    {
      public:
        using ScalarType = typename TestFunction::ScalarType;
        using Parent = Variational::LinearFormIntegratorBase<ScalarType>;
        using PhiType = Variational::RealFunctionBase<PhiDerived>;
        using GradType = Variational::VectorFunctionBase<Real, GradDerived>;

        WNGIRSurfaceForce(
            const PhiType& phi,
            const GradType& grad,
            const TestFunction& v,
            const Displacement& current,
            const WNGIRParameters& parameters,
            Real sigma2)
          : Parent(v.getLeaf()),
            m_phi(phi.copy()),
            m_grad(grad.copy()),
            m_v(v),
            m_current(current),
            m_parameters(parameters),
            m_sigma2(sigma2)
        {}

        WNGIRSurfaceForce(const WNGIRSurfaceForce& other)
          : Parent(other),
            m_phi(other.m_phi ? other.m_phi->copy() : nullptr),
            m_grad(other.m_grad ? other.m_grad->copy() : nullptr),
            m_v(other.m_v),
            m_current(other.m_current),
            m_parameters(other.m_parameters),
            m_sigma2(other.m_sigma2),
            m_polytope(other.m_polytope)
        {}

        const Geometry::Polytope& getPolytope() const final override
        {
          assert(m_polytope);
          return *m_polytope;
        }

        WNGIRSurfaceForce& setPolytope(
            const Geometry::Polytope& polytope) final override
        {
          m_polytope = &polytope;

          const std::size_t faceDim = polytope.getDimension();
          const Index faceIdx = polytope.getIndex();
          const auto& testFES = m_v.getFiniteElementSpace();
          const auto& testFE = testFES.getFiniteElement(faceDim, faceIdx);
          const auto& params = m_parameters.get();
          const std::size_t qOrder = params.quadratureOrder > 0
            ? params.quadratureOrder
            : std::max<std::size_t>(2, 2 * testFE.getOrder());
          const auto& qf =
            QF::PolytopeQuadratureFormula::get(qOrder, polytope.getGeometry());
          const auto& quad = polytope.getQuadrature(qf);

          const std::size_t nte = testFE.getCount();
          m_vector.resize(static_cast<Eigen::Index>(nte));
          m_vector.setZero();

          for (std::size_t qp = 0; qp < quad.getSize(); ++qp)
          {
            const auto& pt = quad.getPoint(qp);
            const Variational::IntegrationPoint ip(pt, &qf, qp);
            const Real w = qf.getWeight(qp) * pt.getDistortion();
            const auto uq = m_current.get().getValue(ip);
            Math::SpatialVector<Real> displacement(polytope.getMesh().getDimension());
            displacement.setZero();
            for (std::size_t r = 0; r < static_cast<std::size_t>(displacement.size()); ++r)
              displacement(static_cast<Eigen::Index>(r)) = uq(r);

            const auto moved =
              makeTranslatedPoint(pt, pt.getPhysicalCoordinates() + displacement,
                                  params.pointLocationTolerance);
            const Real r = m_phi->getValue(moved);
            const auto g = m_grad->getValue(moved);
            const Real omega = std::exp(-r * r / m_sigma2);
            const auto dGamma = (-omega * r) * g;

            const auto& rc = pt.getReferenceCoordinates();
            for (std::size_t te = 0; te < nte; ++te)
            {
              const auto testValue = testFE.getBasis(te)(rc);
              m_vector(static_cast<Eigen::Index>(te))
                += w * dGamma.dot(testValue);
            }
          }
          return *this;
        }

        ScalarType integrate(std::size_t local) final override
        {
          return m_vector(static_cast<Eigen::Index>(local));
        }

        Geometry::Region getRegion() const final override
        {
          return Geometry::Region::Faces;
        }

        WNGIRSurfaceForce* copy() const noexcept final override
        {
          return new WNGIRSurfaceForce(*this);
        }

      private:
        std::unique_ptr<PhiType> m_phi;
        std::unique_ptr<GradType> m_grad;
        TestFunction m_v;
        std::reference_wrapper<const Displacement> m_current;
        std::reference_wrapper<const WNGIRParameters> m_parameters;
        Real m_sigma2;
        const Geometry::Polytope* m_polytope = nullptr;
        Math::Vector<Real> m_vector;
    };
}

#endif
