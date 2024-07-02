#ifndef RODIN_VARIATIONAL_QUADRATURERULE_H
#define RODIN_VARIATIONAL_QUADRATURERULE_H

#include "ForwardDecls.h"

#include "Rodin/QF/GaussLegendre.h"

#include "Dot.h"
#include "Sum.h"
#include "ShapeFunction.h"
#include "LinearFormIntegrator.h"
#include "BilinearFormIntegrator.h"

namespace Rodin::Variational
{
  /**
   * @defgroup QuadratureRuleSpecializations QuadratureRule Template Specializations
   * @brief Template specializations of the QuadratureRule class.
   *
   * @see QuadratureRule
   * @see RodinQuadrature
   */

  /**
   * @brief Quadrature rule on polytope for any function defined on the mesh.
   */
  template <class FunctionDerived>
  class QuadratureRule<FunctionBase<FunctionDerived>> final
    : public FormLanguage::Base
  {
    public:
      using IntegrandType = FunctionBase<FunctionDerived>;

      using IntegrandRangeType = typename FormLanguage::Traits<IntegrandType>::RangeType;

      using ScalarType = typename FormLanguage::Traits<IntegrandRangeType>::ScalarType;

      using Parent = FormLanguage::Base;

      QuadratureRule(const IntegrandType& f)
        : m_integrand(f.copy())
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_polytope(other.m_polytope),
          m_integrand(other.m_integrand->copy()),
          m_qfgg(other.m_qfgg),
          m_qf(other.m_qf)
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_polytope(std::move(other.m_polytope)),
          m_integrand(std::move(other.m_integrand)),
          m_qfgg(std::move(other.m_qfgg)),
          m_qf(std::move(other.m_qf))
      {}

      const Geometry::Polytope& getPolytope() const
      {
        return m_polytope.value().get();
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope)
      {
        m_polytope = polytope;
        if (!m_qf)
        {
          m_qfgg.emplace(polytope.getGeometry());
          m_qf = m_qfgg.value();
        }
        const auto& trans = polytope.getTransformation();
        const auto& qf = m_qf.value().get();
        for (size_t i = 0; i < qf.getSize(); i++)
          m_ps.emplace_back(polytope, trans, std::cref(qf.getPoint(i)));
        return *this;
      }

      ScalarType compute()
      {
        auto& res = m_value.emplace(0);
        const auto& qf = getQuadratureFormula();
        const auto& f = getIntegrand();
        assert(m_ps.size() == qf.getSize());
        for (size_t i = 0; i < m_ps.size(); i++)
          res += qf.getWeight(i) * m_ps[i].getDistortion() * f(m_ps[i]);
        return res;
      }

      inline
      const std::optional<ScalarType>& getValue() const
      {
        return m_value;
      }

      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const QF::QuadratureFormulaBase& getQuadratureFormula() const
      {
        assert(m_qf);
        return m_qf.value().get();
      }

      QuadratureRule& setQuadratureFormula(const QF::QuadratureFormulaBase& qf)
      {
        m_qf = qf;
        return *this;
      }

    private:
      std::optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      std::unique_ptr<IntegrandType> m_integrand;
      std::optional<const QF::GenericPolytopeQuadrature> m_qfgg;
      std::optional<std::reference_wrapper<const QF::QuadratureFormulaBase>> m_qf;
      std::optional<ScalarType> m_value;

      std::vector<Geometry::Point> m_ps;
  };

  /**
   * @ingroup IntegralSpecializations
   * @brief Integration of a GridFunction object.
   */
  template <class FES>
  class QuadratureRule<GridFunction<FES>> : public Integrator
  {
    public:
      using FESType = FES;

      using IntegrandType = GridFunction<FESType>;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using Parent = Integrator;

      /**
       * @brief Constructs the integral object from the given integrand.
       */
      QuadratureRule(const IntegrandType& u)
        : m_u(u),
          m_v(u.getFiniteElementSpace()),
          m_lf(m_v)
      {
        assert(u.getFiniteElementSpace().getVectorDimension() == 1);
      }

      /**
       * @brief Copy constructor.
       */
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_u(other.m_u),
          m_v(other.m_u.get().getFiniteElementSpace()),
          m_lf(m_v)
      {}

      /**
       * @brief Move constructor.
       */
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_v(std::move(other.m_v)),
          m_lf(std::move(other.m_lf))
      {}

      /**
       * @brief Integrates the expression and returns the value.
       *
       * Compute the value of the integral, caches it and returns it.
       *
       * @returns Value of integral
       */
      inline
      ScalarType compute()
      {
        switch (getRegion())
        {
          case Region::Cells:
          {
            auto lfi = Variational::Integral(m_v);
            if (m_attrs.size() > 0)
              lfi.over(m_attrs);
            m_lf.from(lfi).assemble();
            return m_value.emplace(m_lf(m_u));
          }
          case Region::Boundary:
          {
            auto lfi = Variational::BoundaryIntegral(m_v);
            if (m_attrs.size() > 0)
              lfi.over(m_attrs);
            m_lf.from(lfi).assemble();
            return m_value.emplace(m_lf(m_u));
          }
          case Region::Faces:
          {
            auto lfi = Variational::FaceIntegral(m_v);
            if (m_attrs.size() > 0)
              lfi.over(m_attrs);
            m_lf.from(lfi).assemble();
            return m_value.emplace(m_lf(m_u));
          }
          case Region::Interface:
          {
            auto lfi = Variational::InterfaceIntegral(m_v);
            if (m_attrs.size() > 0)
              lfi.over(m_attrs);
            m_lf.from(lfi).assemble();
            return m_value.emplace(m_lf(m_u));
          }
        }
        assert(false);
        return NAN;
      }

      /**
       * @brief Returns the value of the integral, computing it if necessary.
       *
       * If compute() has been called before, returns the value of the cached
       * value. Otherwise, it will call compute() and return the newly computed
       * value.
       *
       * @returns Value of integral
       */
      inline
      operator ScalarType()
      {
        if (!m_value.has_value())
          return compute();
        else
          return m_value.value();
      }

      inline
      QuadratureRule& over(Geometry::Attribute attr)
      {
        return over(FlatSet<Geometry::Attribute>{attr});
      }

      inline
      QuadratureRule& over(const FlatSet<Geometry::Attribute>& attrs)
      {
        m_attrs = attrs;
        return *this;
      }

      inline
      const std::optional<ScalarType>& getValue() const
      {
        return m_value;
      }

      inline
      Type getType() const override
      {
        return Integrator::Type::Linear;
      }

      virtual Region getRegion() const = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::reference_wrapper<const GridFunction<FES>>   m_u;
      TestFunction<FES>                                 m_v;

      FlatSet<Geometry::Attribute>      m_attrs;
      LinearForm<FES, Math::Vector<ScalarType>>     m_lf;

      std::optional<ScalarType> m_value;
  };

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Approximation of the integral of the the dot product between a
   * trial shape function and a test shape function.
   *
   * Represents the quadrature rule approximation of an integral:
   * @f[
   *  \int_{\mathcal{R}_h} \mathrm{IntegrandType} \ dx \approx \sum_{i = 1}^{n}
   *  w_i \ \mathrm{IntegrandType} (x_i)
   * @f]
   * where @f$ \mathcal{R}_h @f$ is some region of the mesh @f$ \mathcal{T}_h
   * @f$, the quadrature point @f$ x_i @f$ has an associated weight @f$ w_i @f$
   * and @f$ \mathrm{IntegrandType}(x_i) @f$ is the value of the integrand at the
   * quadrature point.
   */
  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>,
      ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<
          Dot<
            ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>,
            ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>::ScalarType>
  {
    public:
      using LHSType = ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>;

      using RHSType = ShapeFunctionBase<RHSDerived, TestFES, TestSpace>;

      using IntegrandType = Dot<LHSType, RHSType>;

      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      QuadratureRule(const LHSType& lhs, const RHSType& rhs)
        : QuadratureRule(Dot(lhs, rhs))
      {}

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy())
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy())
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand))
      {}

      inline
      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      inline
      const Geometry::Polytope& getPolytope() const override
      {
        return m_polytope.value().get();
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) override
      {
        m_polytope = polytope;
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& integrand = *m_integrand;
        const auto& trial = integrand.getLHS();
        const auto& test = integrand.getRHS();
        const auto& trans = polytope.getTransformation();
        const auto& trialfes = trial.getFiniteElementSpace();
        const auto& testfes = test.getFiniteElementSpace();
        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe = testfes.getFiniteElement(d, idx);
        const size_t order = std::max(trialfe.getOrder(), testfe.getOrder());
        m_qf.reset(new QF::GenericPolytopeQuadrature(order, polytope.getGeometry()));
        m_ps.clear();
        m_ps.reserve(m_qf->getSize());
        for (size_t i = 0; i < m_qf->getSize(); i++)
          m_ps.emplace_back(polytope, trans, std::cref(m_qf->getPoint(i)));
        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        ScalarType res = 0;
        auto& integrand = *m_integrand;
        for (size_t i = 0; i < m_ps.size(); i++)
        {
          integrand.setPoint(m_ps[i]);
          res += m_qf->getWeight(i) * m_ps[i].getDistortion() * integrand(tr, te);
        }
        return res;
      }

      virtual Integrator::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      std::optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      std::unique_ptr<QF::QuadratureFormulaBase> m_qf;
      std::vector<Geometry::Point> m_ps;
  };

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Approximation of the integral of a test shape function.
   */
  template <class NestedDerived, class FES>
  class QuadratureRule<ShapeFunctionBase<NestedDerived, FES, TestSpace>>
    : public LinearFormIntegratorBase<
        typename FormLanguage::Traits<ShapeFunctionBase<NestedDerived, FES, TestSpace>>::ScalarType>
  {
    public:
      using FESType = FES;

      using IntegrandType = ShapeFunctionBase<NestedDerived, FESType, TestSpace>;

      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent = LinearFormIntegratorBase<ScalarType>;

      template <class LHSDerived, class RHSDerived>
      constexpr
      QuadratureRule(const FunctionBase<LHSDerived>& lhs, const ShapeFunctionBase<RHSDerived, FES, TestSpace>& rhs)
        : QuadratureRule(Dot(lhs, rhs))
      {}

      constexpr
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLeaf()),
          m_integrand(integrand.copy())
      {}

      constexpr
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy())
      {}

      constexpr
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand))
      {}

      inline
      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const Geometry::Polytope& getPolytope() const override
      {
        return m_polytope.value().get();
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& trans = polytope.getTransformation();
        const auto& integrand = getIntegrand();
        const auto& fes = integrand.getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, idx);
        const size_t order = fe.getOrder();
        m_qf.reset(new QF::GenericPolytopeQuadrature(order, polytope.getGeometry()));
        m_ps.clear();
        m_ps.reserve(m_qf->getSize());
        for (size_t i = 0; i < m_qf->getSize(); i++)
          m_ps.emplace_back(polytope, trans, std::cref(m_qf->getPoint(i)));
        return *this;
      }

      ScalarType integrate(size_t local) final override
      {
        ScalarType res = 0;
        auto& integrand = *m_integrand;
        for (size_t i = 0; i < m_ps.size(); i++)
        {
          integrand.setPoint(m_ps[i]);
          res += m_qf->getWeight(i) * m_ps[i].getDistortion() * integrand.getBasis(local);
        }
        return res;
      }

      virtual Integrator::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      std::optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      std::unique_ptr<QF::QuadratureFormulaBase> m_qf;
      std::vector<Geometry::Point> m_ps;
  };
}

#endif
