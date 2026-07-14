/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file WNGIRNormalOffsetInitializer.h
 * @brief Normal-offset initializer for WNGIR interface fitting.
 */
#ifndef RODIN_ADAPTATION_WNGIRNORMALOFFSETINITIALIZER_H
#define RODIN_ADAPTATION_WNGIRNORMALOFFSETINITIALIZER_H

#include "WNGIRDetail.h"
#include "WNGIRParameters.h"

namespace Rodin::Adaptation::Detail
{
  /// @cond RODIN_INTERNAL
  template <class PhiDerived, class GradDerived, class TrialFunction, class TestFunction>
  class WNGIRNormalOffsetMetric final
    : public Variational::LocalBilinearFormIntegratorBase<
        typename TrialFunction::ScalarType>
  {
    public:
      using ScalarType = typename TrialFunction::ScalarType;
      using Parent = Variational::LocalBilinearFormIntegratorBase<ScalarType>;
      using PhiType = Variational::RealFunctionBase<PhiDerived>;
      using GradType = Variational::VectorFunctionBase<Real, GradDerived>;

      WNGIRNormalOffsetMetric(const PhiType& phi, const GradType& grad,
        const TrialFunction& du, const TestFunction& v, const WNGIRParameters& parameters,
        Real sigma2)
        : Parent(du.getLeaf(), v.getLeaf()),
          m_phi(phi.copy()),
          m_grad(grad.copy()),
          m_du(std::cref(du)),
          m_v(std::cref(v)),
          m_parameters(parameters),
          m_sigma2(sigma2)
      {}

      WNGIRNormalOffsetMetric(const WNGIRNormalOffsetMetric& other)
        : Parent(other),
          m_phi(other.m_phi ? other.m_phi->copy() : nullptr),
          m_grad(other.m_grad ? other.m_grad->copy() : nullptr),
          m_du(other.m_du),
          m_v(other.m_v),
          m_parameters(other.m_parameters),
          m_sigma2(other.m_sigma2),
          m_polytope(other.m_polytope)
      {}

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      WNGIRNormalOffsetMetric& setPolytope(
        const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const std::size_t faceDim = polytope.getDimension();
        const Index faceIdx = polytope.getIndex();
        const auto& trialFES = m_du.get().getFiniteElementSpace();
        const auto& testFES = m_v.get().getFiniteElementSpace();
        const auto& trialFE = trialFES.getFiniteElement(faceDim, faceIdx);
        const auto& testFE = testFES.getFiniteElement(faceDim, faceIdx);
        const auto& params = m_parameters.get();
        const std::size_t qOrder = params.quadratureOrder > 0
          ? params.quadratureOrder
          : std::max<std::size_t>(2, 2 * std::max(trialFE.getOrder(), testFE.getOrder()));
        const auto& qf =
          QF::PolytopeQuadratureFormula::get(qOrder, polytope.getGeometry());
        const auto& quad = polytope.getQuadrature(qf);

        const std::size_t ntr = trialFE.getCount();
        const std::size_t nte = testFE.getCount();
        m_matrix.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_matrix.setZero();

        constexpr Real epsG = Real(1e-12);
        for (std::size_t qp = 0; qp < quad.getSize(); ++qp)
        {
          const auto& pt = quad.getPoint(qp);
          const Real w = qf.getWeight(qp) * pt.getDistortion();
          const Real r = m_phi->getValue(pt);
          const auto g = m_grad->getValue(pt);
          const Real gNorm = std::sqrt(std::max(g.dot(g), Real(0)));
          if (gNorm <= Real(1e-14))
            continue;
          const auto n = g / (gNorm + epsG);
          Real s = -r / (gNorm + epsG);
          if (params.initialGuessCapH > Real(0) && params.h > Real(0))
          {
            const Real cap = params.initialGuessCapH * params.h;
            s = std::clamp(s, -cap, cap);
          }
          const Real chi = std::exp(-s * s / m_sigma2);
          const auto& rc = pt.getReferenceCoordinates();
          for (std::size_t te = 0; te < nte; ++te)
          {
            const auto testValue = testFE.getBasis(te)(rc);
            const Real vn = n.dot(testValue);
            for (std::size_t tr = 0; tr < ntr; ++tr)
            {
              const auto trialValue = trialFE.getBasis(tr)(rc);
              m_matrix(static_cast<Eigen::Index>(te), static_cast<Eigen::Index>(tr)) +=
                w * params.initialGuessGamma * chi * n.dot(trialValue) * vn;
            }
          }
        }
        return *this;
      }

      ScalarType integrate(std::size_t tr, std::size_t te) final override
      {
        return m_matrix(static_cast<Eigen::Index>(te), static_cast<Eigen::Index>(tr));
      }

      Geometry::Region getRegion() const final override
      {
        return Geometry::Region::Faces;
      }

      WNGIRNormalOffsetMetric* copy() const noexcept final override
      {
        return new WNGIRNormalOffsetMetric(*this);
      }

    private:
      std::unique_ptr<PhiType> m_phi;
      std::unique_ptr<GradType> m_grad;
      std::reference_wrapper<const TrialFunction> m_du;
      std::reference_wrapper<const TestFunction> m_v;
      std::reference_wrapper<const WNGIRParameters> m_parameters;
      Real m_sigma2;
      const Geometry::Polytope* m_polytope = nullptr;
      Math::Matrix<Real> m_matrix;
  };

  template <class PhiDerived, class GradDerived, class TestFunction>
  class WNGIRNormalOffsetForce final
    : public Variational::LinearFormIntegratorBase<typename TestFunction::ScalarType>
  {
    public:
      using ScalarType = typename TestFunction::ScalarType;
      using Parent = Variational::LinearFormIntegratorBase<ScalarType>;
      using PhiType = Variational::RealFunctionBase<PhiDerived>;
      using GradType = Variational::VectorFunctionBase<Real, GradDerived>;

      WNGIRNormalOffsetForce(const PhiType& phi, const GradType& grad,
        const TestFunction& v, const WNGIRParameters& parameters, Real sigma2)
        : Parent(v.getLeaf()),
          m_phi(phi.copy()),
          m_grad(grad.copy()),
          m_v(std::cref(v)),
          m_parameters(parameters),
          m_sigma2(sigma2)
      {}

      WNGIRNormalOffsetForce(const WNGIRNormalOffsetForce& other)
        : Parent(other),
          m_phi(other.m_phi ? other.m_phi->copy() : nullptr),
          m_grad(other.m_grad ? other.m_grad->copy() : nullptr),
          m_v(other.m_v),
          m_parameters(other.m_parameters),
          m_sigma2(other.m_sigma2),
          m_polytope(other.m_polytope)
      {}

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      WNGIRNormalOffsetForce& setPolytope(
        const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const std::size_t faceDim = polytope.getDimension();
        const Index faceIdx = polytope.getIndex();
        const auto& testFES = m_v.get().getFiniteElementSpace();
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

        constexpr Real epsG = Real(1e-12);
        for (std::size_t qp = 0; qp < quad.getSize(); ++qp)
        {
          const auto& pt = quad.getPoint(qp);
          const Real w = qf.getWeight(qp) * pt.getDistortion();
          const Real r = m_phi->getValue(pt);
          const auto g = m_grad->getValue(pt);
          const Real gNorm = std::sqrt(std::max(g.dot(g), Real(0)));
          if (gNorm <= Real(1e-14))
            continue;
          const auto n = g / (gNorm + epsG);
          Real s = -r / (gNorm + epsG);
          if (params.initialGuessCapH > Real(0) && params.h > Real(0))
          {
            const Real cap = params.initialGuessCapH * params.h;
            s = std::clamp(s, -cap, cap);
          }
          const Real chi = std::exp(-s * s / m_sigma2);
          const auto& rc = pt.getReferenceCoordinates();
          for (std::size_t te = 0; te < nte; ++te)
          {
            const auto testValue = testFE.getBasis(te)(rc);
            m_vector(static_cast<Eigen::Index>(te)) +=
              w * params.initialGuessGamma * chi * s * n.dot(testValue);
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

      WNGIRNormalOffsetForce* copy() const noexcept final override
      {
        return new WNGIRNormalOffsetForce(*this);
      }

    private:
      std::unique_ptr<PhiType> m_phi;
      std::unique_ptr<GradType> m_grad;
      std::reference_wrapper<const TestFunction> m_v;
      std::reference_wrapper<const WNGIRParameters> m_parameters;
      Real m_sigma2;
      const Geometry::Polytope* m_polytope = nullptr;
      Math::Vector<Real> m_vector;
  };
  /// @endcond
}

#endif
