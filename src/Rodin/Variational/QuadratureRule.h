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
    Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>
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
          m_prod(other.m_prod->copy())
      {}

      QuadratureRule(QuadratureRule&& other)
        : BilinearFormIntegratorBase(std::move(other)),
          m_prod(std::move(other.m_prod))
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
        res = Math::Matrix::Zero(test.getDOFs(polytope), trial.getDOFs(polytope));
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
          res.noalias() += qf.getWeight(i) * p.getDistortion() * integrand.getTensorBasis(p).getVector();
        }
      }

      virtual Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<Integrand> m_integrand;
  };

  /* ||-- OPTIMIZATIONS -----------------------------------------------------
   * QuadratureRule<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
   * ---------------------------------------------------------------------->>
   */

  /**
   * @ingroup QuadratureRuleSpecializations
   *
   * @f[
   * \int \nabla u \cdot \nabla v \ dx
   * @f]
   */
  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  class QuadratureRule<Dot<
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, TrialFES, TrialSpace>>, TrialFES, TrialSpace>,
        ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, TestFES, TestSpace>>, TestFES, TestSpace>>>
    : public BilinearFormIntegratorBase
  {
    public:
      using Parent = BilinearFormIntegratorBase;
      using LHS = ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, TrialFES, TrialSpace>>, TrialFES, TrialSpace>;
      using RHS = ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, TestFES, TestSpace>>, TestFES, TestSpace>;
      using Integrand = Dot<
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, TrialFES, TrialSpace>>, TrialFES, TrialSpace>,
        ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, TestFES, TestSpace>>, TestFES, TestSpace>>;

      constexpr
      QuadratureRule(const Integrand& integrand)
        : BilinearFormIntegratorBase(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
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

      void assemble(const Geometry::Polytope& polytope) override
      {
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& integrand = getIntegrand();
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
        res = Math::Matrix::Zero(test.getDOFs(polytope), trial.getDOFs(polytope));
        if (&trialfe == &testfe)
        {
          m_trgrad.resize(d, trialfe.getCount());
          for (size_t i = 0; i < qf.getSize(); i++)
          {
            Geometry::Point p(polytope, trans, std::ref(qf.getPoint(i)));
            const Math::Vector& rc = p.getCoordinates(Geometry::Point::Coordinates::Reference);
            for (size_t local = 0; local < static_cast<size_t>(m_trgrad.cols()); local++)
              m_trgrad.col(local) = trialfe.getGradient(local)(CacheResult, rc);
            const auto jact = p.getJacobianInverse().transpose();
            res.noalias() += qf.getWeight(i) * p.getDistortion() * (jact * m_trgrad).transpose() * (jact * m_trgrad);
          }
        }
        else
        {
          m_trgrad.resize(d, trialfe.getCount());
          m_tegrad.resize(d, testfe.getCount());
          for (size_t i = 0; i < qf.getSize(); i++)
          {
            Geometry::Point p(polytope, trans, std::ref(qf.getPoint(i)));
            const Math::Vector& rc = p.getCoordinates(Geometry::Point::Coordinates::Reference);
            for (size_t local = 0; local < static_cast<size_t>(m_trgrad.cols()); local++)
              m_trgrad.col(local) = trialfe.getGradient(local)(CacheResult, rc);
            for (size_t local = 0; local < static_cast<size_t>(m_tegrad.cols()); local++)
              m_tegrad.col(local) = testfe.getGradient(local)(CacheResult, rc);
            const auto jact = p.getJacobianInverse().transpose();
            res.noalias() += qf.getWeight(i) * p.getDistortion() * (jact * m_tegrad).transpose() * (jact * m_trgrad);
          }
        }
      }

      virtual Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      Math::Matrix m_trgrad;
      Math::Matrix m_tegrad;
      std::unique_ptr<Integrand> m_integrand;
  };

  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  QuadratureRule(const Dot<
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, TrialFES, TrialSpace>>, TrialFES, TrialSpace>,
        ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, TestFES, TestSpace>>, TestFES, TestSpace>>&)
    -> QuadratureRule<Dot<
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, TrialFES, TrialSpace>>, TrialFES, TrialSpace>,
        ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, TestFES, TestSpace>>, TestFES, TestSpace>>>;

  /* <<-- OPTIMIZATIONS -----------------------------------------------------
   * Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
   * ----------------------------------------------------------------------||
   */
}

#endif
