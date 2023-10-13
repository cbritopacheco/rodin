#ifndef RODIN_VARIATIONAL_QUADRATURERULE_H
#define RODIN_VARIATIONAL_QUADRATURERULE_H

#include "ForwardDecls.h"

#include "Rodin/QF/QFGG.h"

#include "Dot.h"
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
   */

  template <class FunctionDerived>
  class QuadratureRule<FunctionBase<FunctionDerived>> final
    : public FormLanguage::Base
  {
    public:
      using Parent = FormLanguage::Base;

      using Integrand = FunctionBase<FunctionDerived>;

      QuadratureRule(
          std::reference_wrapper<const Geometry::Polytope> polytope, const Integrand& f)
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


      Scalar compute()
      {
        auto& res = m_value.emplace(0);
        const auto& qf = m_qf.get();
        const auto& f = getIntegrand();
        const auto& polytope = m_polytope.get();
        const auto& trans = polytope.getTransformation();
        for (size_t i = 0; i < qf.getSize(); i++)
        {
          Geometry::Point p(polytope, trans, std::ref(qf.getPoint(i)));
          res += qf.getWeight(i) * p.getDistortion() * f(p);
        }
        return res;
      }

      inline
      const std::optional<Scalar>& getValue() const
      {
        return m_value;
      }

      const Integrand& getIntegrand() const
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
      std::unique_ptr<Integrand> m_integrand;
      const QF::QFGG m_qfgg;
      std::reference_wrapper<const QF::QuadratureFormulaBase> m_qf;
      std::optional<Scalar> m_value;
  };

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Approximation of the integral of the the dot product between a
   * trial shape function and a test shape function.
   *
   * Represents the quadrature rule approximation of an integral:
   * @f[
   *  \int_{\mathcal{R}_h} \mathrm{Integrand} \ dx \approx \sum_{i = 1}^{n}
   *  w_i \ \mathrm{Integrand} (x_i)
   * @f]
   * where @f$ \mathcal{R}_h @f$ is some region of the mesh @f$ \mathcal{T}_h
   * @f$, the quadrature point @f$ x_i @f$ has an associated weight @f$ w_i @f$
   * and @f$ \mathrm{Integrand}(x_i) @f$ is the value of the integrand at the
   * quadrature point.
   */
  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>,
      ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>
        : public BilinearFormIntegratorBase
  {
    public:
      using LHS = ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>;

      using RHS = ShapeFunctionBase<RHSDerived, TestFES, TestSpace>;

      using Integrand = Dot<LHS, RHS>;

      using Parent = BilinearFormIntegratorBase;

      QuadratureRule(const LHS& lhs, const RHS& rhs)
        : QuadratureRule(Dot(lhs, rhs))
      {}

      QuadratureRule(const Integrand& prod)
        : BilinearFormIntegratorBase(prod.getLHS().getLeaf(), prod.getRHS().getLeaf()),
          m_prod(prod.copy())
      {}

      QuadratureRule(const QuadratureRule& other)
        : BilinearFormIntegratorBase(other),
          m_prod(other.m_prod->copy()),
          m_mat(other.m_mat)
      {}

      QuadratureRule(QuadratureRule&& other)
        : BilinearFormIntegratorBase(std::move(other)),
          m_prod(std::move(other.m_prod)),
          m_mat(std::move(other.m_mat))
      {}

      inline
      constexpr
      const Integrand& getIntegrand() const
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
        const size_t vc = Geometry::Polytope::getVertexCount(polytope.getGeometry());
        const size_t order = trialfe.getCount() + testfe.getCount() + vc;
        const QF::QFGG qf(order, polytope.getGeometry());
        auto& res = getMatrix();
        res.resize(test.getDOFs(polytope), trial.getDOFs(polytope));
        res.setZero();
        for (size_t i = 0; i < qf.getSize(); i++)
        {
          Geometry::Point p(polytope, trans, std::ref(qf.getPoint(i)));
          integrand.assemble(p);
          res.noalias() += qf.getWeight(i) * p.getDistortion() * integrand.getMatrix();
        }
      }

      virtual Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<Integrand> m_prod;

      Math::Matrix m_mat;
  };

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Approximation of the integral of a test shape function.
   */
  template <class NestedDerived, class FES>
  class QuadratureRule<ShapeFunctionBase<NestedDerived, FES, TestSpace>>
    : public LinearFormIntegratorBase
  {
    public:
      using Integrand = ShapeFunctionBase<NestedDerived, FES, TestSpace>;
      using Parent = LinearFormIntegratorBase;

      template <class LHSDerived, class RHSDerived>
      constexpr
      QuadratureRule(const FunctionBase<LHSDerived>& lhs, const ShapeFunctionBase<RHSDerived, FES, TestSpace>& rhs)
        : QuadratureRule(Dot(lhs, rhs))
      {}

      constexpr
      QuadratureRule(const Integrand& integrand)
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
      const Integrand& getIntegrand() const
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
        const size_t vc = Geometry::Polytope::getVertexCount(polytope.getGeometry());
        assert(integrand.getRangeType() == RangeType::Scalar);
        const size_t order = fe.getCount() + vc;
        const QF::QFGG qf(order, polytope.getGeometry());
        auto& res = getVector();
        res = Math::Vector::Zero(integrand.getDOFs(polytope));
        for (size_t i = 0; i < qf.getSize(); i++)
        {
          Geometry::Point p(polytope, trans, std::ref(qf.getPoint(i)));
          auto basis = integrand.getTensorBasis(p);
          for (size_t local = 0; local < basis.getDOFs(); local++)
            res.coeffRef(local) += qf.getWeight(i) * p.getDistortion() * basis(local);
        }
      }

      virtual Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<Integrand> m_integrand;
  };
}

#endif
