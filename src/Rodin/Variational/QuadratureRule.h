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

      using NumberType = typename FormLanguage::Traits<IntegrandRangeType>::NumberType;

      using Parent = FormLanguage::Base;

      QuadratureRule(
          std::reference_wrapper<const Geometry::Polytope> polytope, const IntegrandType& f)
        : m_polytope(polytope),
          m_integrand(f.copy()),
          m_qfgg(polytope.get().getGeometry()),
          m_qf(m_qfgg)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_polytope(other.m_polytope.get()),
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


      NumberType compute()
      {
        auto& res = m_value.emplace(0);
        const auto& qf = m_qf.get();
        const auto& f = getIntegrand();
        const auto& polytope = m_polytope.get();
        const auto& trans = polytope.getTransformation();
        for (size_t i = 0; i < qf.getSize(); i++)
        {
          const Geometry::Point p(polytope, trans, std::cref(qf.getPoint(i)));
          res += qf.getWeight(i) * p.getDistortion() * f(p);
        }
        return res;
      }

      inline
      const std::optional<NumberType>& getValue() const
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
        return m_qf.get();
      }

      QuadratureRule& setQuadratureFormula(const QF::QuadratureFormulaBase& qf)
      {
        m_qf = qf;
        return *this;
      }

    private:
      std::reference_wrapper<const Geometry::Polytope> m_polytope;
      std::unique_ptr<IntegrandType> m_integrand;
      const QF::GenericPolytopeQuadrature m_qfgg;
      std::reference_wrapper<const QF::QuadratureFormulaBase> m_qf;
      std::optional<NumberType> m_value;
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

      using NumberType = typename FormLanguage::Traits<FESType>::NumberType;

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
      NumberType compute()
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
      operator NumberType()
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
      const std::optional<NumberType>& getValue() const
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
      LinearForm<FES, Math::Vector<NumberType>>     m_lf;

      std::optional<NumberType> m_value;
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
    : public LocalBilinearFormIntegratorBase
  {
    public:
      using LHSType = ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>;

      using RHSType = ShapeFunctionBase<RHSDerived, TestFES, TestSpace>;

      using IntegrandType = Dot<LHSType, RHSType>;

      using IntegrandRangeType = typename FormLanguage::Traits<IntegrandType>::RangeType;

      using NumberType = typename FormLanguage::Traits<IntegrandRangeType>::NumberType;

      using Parent = LocalBilinearFormIntegratorBase;

      QuadratureRule(const LHSType& lhs, const RHSType& rhs)
        : QuadratureRule(Dot(lhs, rhs))
      {}

      QuadratureRule(const IntegrandType& prod)
        : LocalBilinearFormIntegratorBase(prod.getLHS().getLeaf(), prod.getRHS().getLeaf()),
          m_prod(prod.copy())
      {}

      QuadratureRule(const QuadratureRule& other)
        : LocalBilinearFormIntegratorBase(other),
          m_prod(other.m_prod->copy())
      {}

      QuadratureRule(QuadratureRule&& other)
        : LocalBilinearFormIntegratorBase(std::move(other)),
          m_prod(std::move(other.m_prod))
      {}

      inline
      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_prod);
        return *m_prod;
      }

      void assemble(const Geometry::Polytope& polytope) final override
      {
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        auto& integrand = *m_prod;
        const auto& trial = integrand.getLHS();
        const auto& test = integrand.getRHS();
        const auto& trans = polytope.getTransformation();
        const auto& trialfes = trial.getFiniteElementSpace();
        const auto& testfes = test.getFiniteElementSpace();
        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe = testfes.getFiniteElement(d, idx);
        const size_t order = std::max(trialfe.getOrder(), testfe.getOrder());
        const QF::GenericPolytopeQuadrature qf(order, polytope.getGeometry());
        auto& res = getMatrix();
        res.resize(test.getDOFs(polytope), trial.getDOFs(polytope));
        res.setZero();
        for (size_t i = 0; i < qf.getSize(); i++)
        {
          const Geometry::Point p(polytope, trans, std::cref(qf.getPoint(i)));
          integrand.assemble(p);
          res.noalias() += qf.getWeight(i) * p.getDistortion() * integrand.getMatrix();
        }
      }

      virtual Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_prod;
  };

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Approximation of the integral of a test shape function.
   */
  template <class NestedDerived, class FES>
  class QuadratureRule<ShapeFunctionBase<NestedDerived, FES, TestSpace>>
    : public LinearFormIntegratorBase<typename FormLanguage::Traits<FES>::NumberType>
  {
    public:
      using FESType = FES;

      using NumberType = typename FormLanguage::Traits<FESType>::NumberType;

      using IntegrandType = ShapeFunctionBase<NestedDerived, FESType, TestSpace>;

      using Parent = LinearFormIntegratorBase<NumberType>;

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

      void assemble(const Geometry::Polytope& polytope) final override
      {
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& trans = polytope.getTransformation();
        const auto& integrand = getIntegrand();
        const auto& fes = integrand.getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, idx);
        const size_t order = fe.getOrder();
        const QF::GenericPolytopeQuadrature qf(order, polytope.getGeometry());
        auto& res = this->getVector();
        res.resize(integrand.getDOFs(polytope));
        res.setZero();
        for (size_t i = 0; i < qf.getSize(); i++)
        {
          const Geometry::Point p(polytope, trans, std::cref(qf.getPoint(i)));
          auto basis = integrand.getTensorBasis(p);
          const Scalar distortion = p.getDistortion();
          for (size_t local = 0; local < basis.getDOFs(); local++)
            res.coeffRef(local) += qf.getWeight(i) * distortion * basis(local);
        }
      }

      virtual Integrator::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;
  };
}

#endif
