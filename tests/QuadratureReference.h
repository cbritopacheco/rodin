/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TESTS_QUADRATUREREFERENCE_H
#define RODIN_TESTS_QUADRATUREREFERENCE_H

#include "Rodin/Variational.h"

namespace Rodin::Tests
{
  /**
   * @brief Test-only pre-reuse bilinear quadrature algorithm.
   * Retains the expression and matrix allocation across binds, as the original
   * production rule did. There is no alternate production dispatch path.
   * Every test expression is evaluated inside the trial/test pair loop.
   */
  template <class Integral>
  class QuadratureReference
  {
    public:
      using Integrand = typename Integral::IntegrandType;
      using Scalar = typename FormLanguage::Traits<Integrand>::ScalarType;

      explicit QuadratureReference(const Integral& integral)
        : m_integrand(integral.getIntegrand()),
          m_matrix(),
          m_qf(nullptr),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      void assemble(const Geometry::Polytope& cell, size_t order)
      {
        auto& trial = m_integrand.getLHS();
        auto& test = m_integrand.getRHS();
        const size_t nr = trial.getDOFs(cell);
        const size_t nt = test.getDOFs(cell);
        m_matrix.resize(nt, nr);
        m_matrix.setZero();
        if (!m_qf || m_order != order || m_geometry != cell.getGeometry())
        {
          m_order = order;
          m_geometry = cell.getGeometry();
          m_qf = &QF::PolytopeQuadratureFormula::get(order, m_geometry);
        }
        const auto& quadrature = cell.getQuadrature(*m_qf);
        for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
        {
          const auto& p = quadrature.getPoint(qp);
          const Scalar weight = static_cast<Scalar>(m_qf->getWeight(qp)) *
            static_cast<Scalar>(p.getDistortion());
          const Variational::IntegrationPoint ip(p, m_qf, qp);
          m_integrand.setIntegrationPoint(ip);
          for (size_t tr = 0; tr < nr; ++tr)
          {
            const auto& lhs = trial.getBasis(tr);
            for (size_t te = 0; te < nt; ++te)
              m_matrix(te, tr) += weight * Math::dot(lhs, test.getBasis(te));
          }
        }
      }

      const Math::Matrix<Scalar>& getOperator() const
      {
        return m_matrix;
      }

    private:
      Integrand m_integrand;
      Math::Matrix<Scalar> m_matrix;
      const QF::QuadratureFormulaBase* m_qf;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;
  };

  /** @brief Original-loop integrator for backend-equivalence tests only. */
  template <class Integral>
  class ReferenceIntegral final : public Variational::LocalBilinearFormIntegratorBase<
                                    typename QuadratureReference<Integral>::Scalar>
  {
    public:
      using Scalar = typename QuadratureReference<Integral>::Scalar;
      using Parent = Variational::LocalBilinearFormIntegratorBase<Scalar>;

      explicit ReferenceIntegral(const Integral& integral)
        : Parent(integral),
          m_integral(integral),
          m_reference(integral),
          m_polytope(nullptr)
      {}

      ReferenceIntegral(const ReferenceIntegral& other)
        : Parent(other),
          m_integral(other.m_integral),
          m_reference(m_integral),
          m_polytope(nullptr)
      {}

      ReferenceIntegral& setPolytope(const Geometry::Polytope& cell) override
      {
        m_polytope = &cell;
        m_reference.assemble(
          cell, this->getOrder(cell).value_or(m_integral.getOrder(cell).value_or(6)));
        return *this;
      }

      const Geometry::Polytope& getPolytope() const override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      Scalar integrate(size_t tr, size_t te) override
      {
        return m_reference.getOperator()(te, tr);
      }

      Geometry::Region getRegion() const override
      {
        return m_integral.getRegion();
      }

      ReferenceIntegral* copy() const noexcept override
      {
        return new ReferenceIntegral(*this);
      }

    private:
      Integral m_integral;
      QuadratureReference<Integral> m_reference;
      const Geometry::Polytope* m_polytope;
  };
}
#endif
