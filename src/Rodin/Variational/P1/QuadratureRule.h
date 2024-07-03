#ifndef RODIN_VARIATIONAL_P1_QUADRATURERULE_H
#define RODIN_VARIATIONAL_P1_QUADRATURERULE_H

#include "Rodin/Variational/QuadratureRule.h"
#include "Rodin/QF/QF1P1.h"
#include "Rodin/QF/GrundmannMoller.h"

#include "P1.h"
#include "P1Element.h"

namespace Rodin::Variational
{
  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Integration of a P1 ShapeFunction.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int v \ dx \: ,
   * @f]
   * where @f$ v \in \mathbb{P}_1 @f$.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int v \ dx :
   * \texttt{QuadratureRule}}
   * {\vdash v : \mathbb{P}_1}
   * @f]
   */
  template <class NestedDerived, class Range, class Mesh>
  class QuadratureRule<
    ShapeFunctionBase<
      ShapeFunction<NestedDerived, P1<Range, Mesh>, TestSpace>, P1<Range, Mesh>, TestSpace>>
    : public LinearFormIntegratorBase<
        typename FormLanguage::Traits<
          ShapeFunctionBase<ShapeFunction<NestedDerived, P1<Range, Mesh>, TestSpace>, P1<Range, Mesh>, TestSpace>>
        ::ScalarType>
  {
    public:
      using FESType = P1<Range, Mesh>;

      using IntegrandType =
        ShapeFunctionBase<ShapeFunction<NestedDerived, FESType, TestSpace>>;

      using IntegrandRangeType = typename FormLanguage::Traits<FESType>::RangeType;

      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent = LinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<IntegrandRangeType, ScalarType>);

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

      inline
      const Geometry::Polytope& getPolytope() const final override
      {
        return m_polytope.value().get();
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;
        const auto& trans = polytope.getTransformation();
        m_qf.emplace(polytope.getGeometry());
        assert(m_qf->getSize() == 1);
        m_p.emplace(polytope, trans, std::cref(m_qf->getPoint(0)));
        return *this;
      }

      ScalarType integrate(size_t local) final override
      {
        const auto& polytope = getPolytope();
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& integrand = getIntegrand().getDerived();
        const auto& fes = integrand.getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, idx);
        const auto& qf = *m_qf;
        assert(qf.getSize() == 1);
        const auto& p = m_p.value();
        return qf.getWeight(0) * p.getDistortion() * fe.getBasis(local)(qf.getPoint(0));
      }

      virtual Integrator::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      std::optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      std::optional<QF::QF1P1> m_qf;
      std::optional<Geometry::Point> m_p;
  };

  /**
   * @ingroup RodinCTAD
   */
  template <class NestedDerived, class Range, class Mesh>
  QuadratureRule(const ShapeFunctionBase<ShapeFunction<NestedDerived, P1<Range, Mesh>, TestSpace>>&)
    -> QuadratureRule<ShapeFunctionBase<ShapeFunction<NestedDerived, P1<Range, Mesh>, TestSpace>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Integration of the Dot product of some coefficient function and a
   * P1 ShapeFunction.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int f \cdot v \ dx \: ,
   * @f]
   * where @f$ v \in \mathbb{P}_1 @f$.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int f \cdot v \ dx :
   * \texttt{QuadratureRule}}
   * {\vdash v : \mathbb{P}_1}
   * @f]
   */
  template <class LHSDerived, class RHSDerived, class Range, class Mesh>
  class QuadratureRule<
    ShapeFunctionBase<
      Dot<
        FunctionBase<LHSDerived>,
        ShapeFunctionBase<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>, P1<Range, Mesh>, TestSpace>>,
          P1<Range, Mesh>, TestSpace>>
    : public LinearFormIntegratorBase<
        typename FormLanguage::Traits<
          ShapeFunctionBase<
            Dot<
              FunctionBase<LHSDerived>,
              ShapeFunctionBase<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>>
        ::ScalarType>
  {
    public:
      using FESType = P1<Range, Mesh>;

      using LHSType = FunctionBase<LHSDerived>;

      using RHSType =
        ShapeFunctionBase<ShapeFunction<RHSDerived, FESType, TestSpace>, FESType, TestSpace>;

      using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;

      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

      using IntegrandType = ShapeFunctionBase<Dot<LHSType, RHSType>>;

      using IntegrandRangeType = typename FormLanguage::Traits<IntegrandType>::RangeType;

      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent = LinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);

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

      inline
      const Geometry::Polytope& getPolytope() const final override
      {
        return m_polytope.value().get();
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;
        const auto& trans = polytope.getTransformation();
        m_qf.emplace(polytope.getGeometry());
        assert(m_qf->getSize() == 1);
        m_p.emplace(polytope, trans, std::cref(m_qf->getPoint(0)));
        return *this;
      }

      ScalarType integrate(size_t local) final override
      {
        const auto& polytope = getPolytope();
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& integrand = getIntegrand().getDerived();
        const auto& f = integrand.getLHS();
        const auto& fes = integrand.getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, idx);
        const auto& qf = *m_qf;
        assert(qf.getSize() == 1);
        const auto& p = m_p.value();
        const auto& w = qf.getWeight(0);
        const auto& rc = qf.getPoint(0);
        if constexpr (std::is_same_v<ScalarType, LHSRangeType>)
        {
          return w * p.getDistortion() * f(p) * fe.getBasis(local)(rc);
        }
        else if constexpr (std::is_same_v<Math::Vector<ScalarType>, LHSRangeType>)
        {
          return w * p.getDistortion() * f(p).dot(fe.getBasis(local)(rc));
        }
        else if constexpr (std::is_same_v<Math::Matrix<ScalarType>, LHSRangeType>)
        {
          return w * p.getDistortion() * (f(p).array() *
              fe.getBasis(local)(rc).array()).rowwise().sum().colwise().sum().value();
        }
        else
        {
          assert(false);
          return NAN;
        }
      }

      virtual Integrator::Region getRegion() const override = 0;


      virtual QuadratureRule* copy() const noexcept override = 0;

    private:

      std::unique_ptr<IntegrandType> m_integrand;

      std::optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      std::optional<QF::QF1P1> m_qf;
      std::optional<Geometry::Point> m_p;
  };

  /**
   * @ingroup RodinCTAD
   */
  template <class LHSDerived, class RHSDerived, class Range, class Mesh>
  QuadratureRule(
      const ShapeFunctionBase<
        Dot<FunctionBase<LHSDerived>,
        ShapeFunctionBase<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>&)
    -> QuadratureRule<ShapeFunctionBase<
        Dot<
          FunctionBase<LHSDerived>,
          ShapeFunctionBase<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Integration of the anisotropic Dot product of two instances of the
   * P1 ShapeFunction.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int (A u) \cdot v \ dx \: ,
   * @f]
   * where @f$ u \in \mathbb{P}_1 @f$ and @f$ v \in \mathbb{P}_1 @f$, and @f$
   * A @f$ is a coefficient function.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int (A u) \cdot v \ dx :
   * \texttt{QuadratureRule}}
   * {\vdash u, v : \mathbb{P}_1}
   * @f]
   */
  template <
    class CoefficientDerived, class LHSDerived, class RHSDerived,
    class LHSRange, class RHSRange, class LHSMesh, class RHSMesh>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Mult<
          FunctionBase<CoefficientDerived>,
          ShapeFunctionBase<ShapeFunction<LHSDerived, P1<LHSRange, LHSMesh>, TrialSpace>,
            P1<LHSRange, LHSMesh>, TrialSpace>>, P1<LHSRange, LHSMesh>, TrialSpace>,
      ShapeFunctionBase<
        ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>, P1<RHSRange, RHSMesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<
          Dot<
            ShapeFunctionBase<
              Mult<
                FunctionBase<CoefficientDerived>,
                ShapeFunctionBase<ShapeFunction<LHSDerived, P1<LHSRange, LHSMesh>, TrialSpace>>>>,
            ShapeFunctionBase<
              ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>>>>
        ::ScalarType>
  {
    public:
      using LHSFESType = P1<LHSRange, LHSMesh>;

      using RHSFESType = P1<RHSRange, RHSMesh>;

      using CoefficientType = FunctionBase<CoefficientDerived>;

      using MultiplicandType =
        ShapeFunctionBase<ShapeFunction<LHSDerived, LHSFESType, TrialSpace>>;

      using LHSType = ShapeFunctionBase<
        Mult<CoefficientType, MultiplicandType>>;

      using RHSType = ShapeFunctionBase<ShapeFunction<RHSDerived, RHSFESType, TestSpace>>;

      using IntegrandType = Dot<LHSType, RHSType>;

      using CoefficientRangeType = typename FormLanguage::Traits<CoefficientType>::RangeType;

      using MultiplicandRangeType = typename FormLanguage::Traits<MultiplicandType>::RangeType;

      using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;

      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);

      constexpr
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
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

      inline
      const Geometry::Polytope& getPolytope() const final override
      {
        return m_polytope.value().get();
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;
        const auto& trans = polytope.getTransformation();
        m_qf.emplace(polytope.getGeometry());
        assert(m_qf->getSize() == 1);
        m_p.emplace(polytope, trans, std::cref(m_qf->getPoint(0)));
        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        const auto& polytope = getPolytope();
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS().getDerived();
        const auto& coeff = lhs.getLHS();
        const auto& multiplicand = lhs.getRHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = multiplicand.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();
        const auto& qf = *m_qf;
        assert(qf.getSize() == 1);
        const auto& p = m_p.value();
        const auto& w = qf.getWeight(0);
        const auto& rc = qf.getPoint(0);
        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe = testfes.getFiniteElement(d, idx);
        if constexpr (std::is_same_v<MultiplicandRangeType, ScalarType>)
        {
          if constexpr (std::is_same_v<RHSType, ScalarType>)
          {
            static_assert(std::is_same_v<CoefficientType, ScalarType>);
            const ScalarType s = coeff.getValue(p);
            return w * p.getDistortion() * s * testfe.getBasis(te)(rc) * trialfe.getBasis(tr)(rc);
          }
          else
          {
            assert(false);
            return NAN;
          }
        }
        else
        {
          assert(false); // Not handled yet
          return NAN;
        }
      }

      virtual Integrator::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      std::optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      std::optional<QF::QF1P1> m_qf;
      std::optional<Geometry::Point> m_p;
  };

  template <class CoefficientDerived, class LHSDerived, class RHSDerived, class Number, class Mesh>
  QuadratureRule(const
    Dot<
      ShapeFunctionBase<
        Mult<
          FunctionBase<CoefficientDerived>,
          ShapeFunctionBase<ShapeFunction<LHSDerived, P1<Number, Mesh>, TrialSpace>>>>,
      ShapeFunctionBase<
        ShapeFunction<RHSDerived, P1<Number, Mesh>, TestSpace>>>&)
  ->
  QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Mult<
          FunctionBase<CoefficientDerived>,
          ShapeFunctionBase<ShapeFunction<LHSDerived, P1<Number, Mesh>, TrialSpace>>>>,
      ShapeFunctionBase<

        ShapeFunction<RHSDerived, P1<Number, Mesh>, TestSpace>>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Integration of the isotropic Dot product of two instances of the P1
   * Grad of ShapeFunction.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int \nabla u \cdot \nabla v \ dx \: ,
   * @f]
   * where @f$ u \in \mathbb{P}_1 @f$ and @f$ v \in \mathbb{P}_1 @f$, and @f$ A
   * @f$ is a coefficient function.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int \nabla u \cdot \nabla v \ dx :
   * \texttt{QuadratureRule}}
   * {\vdash u, v : \mathbb{P}_1}
   * @f]
   */
  template <class LHSDerived, class RHSDerived, class LHSRange, class RHSRange, class LHSMesh, class RHSMesh>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Grad<ShapeFunction<LHSDerived, P1<LHSRange, LHSMesh>, TrialSpace>>, P1<LHSRange, LHSMesh>, TrialSpace>,
      ShapeFunctionBase<
        Grad<ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>>, P1<RHSRange, RHSMesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<
          Dot<
            ShapeFunctionBase<
              Grad<ShapeFunction<LHSDerived, P1<LHSRange, LHSMesh>, TrialSpace>>, P1<LHSRange, LHSMesh>, TrialSpace>,
            ShapeFunctionBase<
              Grad<ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>>, P1<RHSRange, RHSMesh>, TestSpace>>>
        ::ScalarType>
  {
    public:
      using LHSFESType = P1<LHSRange, LHSMesh>;

      using RHSFESType = P1<RHSRange, RHSMesh>;

      using LHSType =
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, LHSFESType, TrialSpace>>>;

      using LHSOperandType =
        ShapeFunction<LHSDerived, LHSFESType, TrialSpace>;

      using LHSOperandRangeType =
        typename FormLanguage::Traits<LHSOperandType>::RangeType;

      using RHSType =
        ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, RHSFESType, TestSpace>>>;

      using RHSOperandType =
        ShapeFunction<RHSDerived, RHSFESType, TestSpace>;

      using RHSOperandRangeType =
        typename FormLanguage::Traits<RHSOperandType>::RangeType;

      using IntegrandType = Dot<LHSType, RHSType>;

      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<LHSOperandRangeType, ScalarType>);

      static_assert(std::is_same_v<RHSOperandRangeType, ScalarType>);

      constexpr
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
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

      inline
      const Geometry::Polytope& getPolytope() const final override
      {
        return m_polytope.value().get();
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;
        const auto& trans = polytope.getTransformation();
        m_qf.emplace(polytope.getGeometry());
        assert(m_qf->getSize() == 1);
        m_p.emplace(polytope, trans, std::cref(m_qf->getPoint(0)));
        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        const auto& polytope = getPolytope();
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();
        const auto& qf = *m_qf;
        const auto& p = m_p.value();
        const auto& w = qf.getWeight(0);
        const auto& rc = qf.getPoint(0);
        if (trialfes == testfes)
        {
          const auto& fe = trialfes.getFiniteElement(d, idx);
          if (tr == te)
          {
            fe.getGradient(tr)(m_grad1, rc);
            return w * p.getDistortion() * (p.getJacobianInverse().transpose() * m_grad1).squaredNorm();
          }
          else
          {
            fe.getGradient(tr)(m_grad1, rc);
            fe.getGradient(te)(m_grad2, rc);
            return w * p.getDistortion() * (
                p.getJacobianInverse().transpose() * m_grad2).dot(p.getJacobianInverse().transpose() * m_grad1);
          }
        }
        else
        {
          const auto& trialfe = lhs.getFiniteElementSpace().getFiniteElement(d, idx);
          const auto& testfe = rhs.getFiniteElementSpace().getFiniteElement(d, idx);
          assert(qf.getSize() == 1);
          trialfe.getGradient(tr)(m_grad1, rc);
          testfe.getGradient(te)(m_grad2, rc);
          return w * p.getDistortion() * (
              p.getJacobianInverse().transpose() * m_grad2).dot(p.getJacobianInverse().transpose() * m_grad1);
        }
      }

      virtual Integrator::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      std::optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      std::optional<QF::QF1P1> m_qf;
      std::optional<Geometry::Point> m_p;

      Math::SpatialVector<ScalarType> m_grad1, m_grad2;
  };

  template <class LHSDerived, class RHSDerived, class Range, class Mesh>
  QuadratureRule(
      const Dot<
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>>>,
        ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>&)
    -> QuadratureRule<Dot<
          ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>>>,
          ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Integration of the anisotropic Dot product of two instances of the
   * P1 Grad of ShapeFunction.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int (A \nabla u) \cdot \nabla v \ dx \: ,
   * @f]
   * where @f$ u \in \mathbb{P}_1 @f$ and @f$ v \in \mathbb{P}_1 @f$, and @f$ A
   * @f$ is a coefficient function.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
      static_assert(std::is_same_v<LHSOperandRangeType, ScalarType>);
   * \dfrac
   * {\vdash \int (A \nabla u) \cdot \nabla v \ dx :
   * \texttt{QuadratureRule}}
   * {\vdash u, v : \mathbb{P}_1}
   * @f]
   */
  template <
    class CoefficientDerived, class LHSDerived, class RHSDerived,
    class LHSRange, class RHSRange, class LHSMesh, class RHSMesh>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Mult<
          FunctionBase<CoefficientDerived>,
          ShapeFunctionBase<

            Grad<ShapeFunction<LHSDerived, P1<LHSRange, LHSMesh>, TrialSpace>>, P1<LHSRange, LHSMesh>, TrialSpace>>,
          P1<LHSRange, LHSMesh>, TrialSpace>,
      ShapeFunctionBase<
        Grad<ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>>, P1<RHSRange, RHSMesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<
          Dot<
            ShapeFunctionBase<
              Mult<
                FunctionBase<CoefficientDerived>,
                ShapeFunctionBase<
                  Grad<ShapeFunction<LHSDerived, P1<LHSRange, LHSMesh>, TrialSpace>>, P1<LHSRange, LHSMesh>, TrialSpace>>,
                P1<LHSRange, LHSMesh>, TrialSpace>,
            ShapeFunctionBase<
              Grad<ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>>, P1<RHSRange, RHSMesh>, TestSpace>>>
        ::ScalarType>
  {
    public:
      using LHSFESType = P1<LHSRange, LHSMesh>;

      using RHSFESType = P1<RHSRange, RHSMesh>;

      using CoefficientType = FunctionBase<CoefficientDerived>;

      using MultiplicandType =
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, LHSFESType, TrialSpace>>>;

      using MultiplicandOperandType =
        ShapeFunction<LHSDerived, LHSFESType, TrialSpace>;

      using CoefficientRangeType = typename FormLanguage::Traits<CoefficientType>::RangeType;

      using MultiplicandRangeType = typename FormLanguage::Traits<MultiplicandType>::RangeType;

      using LHSType = ShapeFunctionBase<Mult<CoefficientType, MultiplicandType>>;

      using MultiplicandOperandRangeType =
        typename FormLanguage::Traits<MultiplicandOperandType>::RangeType;

      using RHSType =
        ShapeFunctionBase<
          Grad<ShapeFunction<RHSDerived, RHSFESType, TestSpace>>>;

      using RHSOperandType =
        ShapeFunction<RHSDerived, RHSFESType, TestSpace>;

      using RHSOperandRangeType =
        typename FormLanguage::Traits<RHSOperandType>::RangeType;

      using IntegrandType = Dot<LHSType, RHSType>;

      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<MultiplicandOperandRangeType, ScalarType>);

      static_assert(std::is_same_v<RHSOperandRangeType, ScalarType>);

      constexpr
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
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

      inline
      const Geometry::Polytope& getPolytope() const final override
      {
        return m_polytope.value().get();
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;
        const auto& trans = polytope.getTransformation();
        m_qf.emplace(polytope.getGeometry());
        assert(m_qf->getSize() == 1);
        m_p.emplace(polytope, trans, std::cref(m_qf->getPoint(0)));
        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        const auto& polytope = getPolytope();
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS().getDerived();
        const auto& coeff = lhs.getLHS();
        const auto& multiplicand = lhs.getRHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = rhs.getFiniteElementSpace();
        const auto& testfes = multiplicand.getFiniteElementSpace();
        const auto& qf = *m_qf;
        assert(qf.getSize() == 1);
        const auto& p = m_p.value();
        const auto& w = qf.getWeight(0);
        const auto& rc = qf.getPoint(0);
        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe = testfes.getFiniteElement(d, idx);
        if constexpr (std::is_same_v<CoefficientRangeType, ScalarType>)
        {
          const ScalarType s = coeff.getValue(p);
          trialfe.getGradient(tr)(m_grad1, rc);
          testfe.getGradient(te)(m_grad2, rc);
          return w * p.getDistortion() * s * (
              p.getJacobianInverse().transpose() * m_grad2).dot(
                p.getJacobianInverse().transpose() * m_grad1);
        }
        else if constexpr (std::is_same_v<CoefficientRangeType, Math::Matrix<ScalarType>>)
        {
          coeff.getValue(m_mv, p);
          trialfe.getGradient(tr)(m_grad1, rc);
          testfe.getGradient(te)(m_grad2, rc);
          return w * p.getDistortion() * (
              p.getJacobianInverse().transpose() * m_grad2).dot(
                m_mv * p.getJacobianInverse().transpose() * m_grad1);
        }
        else
        {
          assert(false);
          return NAN;
        }
      }

      virtual Integrator::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      std::optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      std::optional<QF::QF1P1> m_qf;
      std::optional<Geometry::Point> m_p;

      Math::Matrix<ScalarType> m_mv;
      Math::SpatialVector<ScalarType> m_grad1, m_grad2;
  };

  /**
   * @ingroup RodinCTAD
   */
  template <class LHSFunctionDerived, class LHSDerived, class RHSDerived, class Range, class Mesh>
  QuadratureRule(
      const Dot<
        ShapeFunctionBase<
          Mult<
            FunctionBase<LHSFunctionDerived>,
            ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>>>>>,
        ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>&)
    -> QuadratureRule<
          Dot<ShapeFunctionBase<
            Mult<
              FunctionBase<LHSFunctionDerived>,
              ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>>>>>,
          ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Integration of the isotropic Frobenius inner product two instances
   * of the P1 Jacobian of ShapeFunction.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int \mathbf{J} \: u : \mathbf{J} \: v \ dx \: ,
   * @f]
   * where @f$ u \in \mathbb{P}_1 @f$ and @f$ v \in
   * \mathbb{P}_1 @f$, and @f$ A @f$ is coefficient function.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int \mathbf{J} \: u : \mathbf{J} \: v \ dx :
   * \texttt{QuadratureRule}}
   * {\vdash u, v : \mathbb{P}_1}
   * @f]
   */
  template <class LHSDerived, class RHSDerived, class LHSRange, class RHSRange, class LHSMesh, class RHSMesh>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Jacobian<ShapeFunction<LHSDerived, P1<LHSRange, LHSMesh>, TrialSpace>>,
          P1<LHSRange, LHSMesh>, TrialSpace>,
      ShapeFunctionBase<
        Jacobian<ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>>,
          P1<RHSRange, RHSMesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<Dot<
          ShapeFunctionBase<
            Jacobian<ShapeFunction<LHSDerived, P1<LHSRange, LHSMesh>, TrialSpace>>,
              P1<LHSRange, LHSMesh>, TrialSpace>,
          ShapeFunctionBase<
            Jacobian<ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>>,
              P1<RHSRange, RHSMesh>, TestSpace>>>
        ::ScalarType>
  {
    public:
      using LHSFESType = P1<LHSRange, LHSMesh>;

      using RHSFESType = P1<RHSRange, RHSMesh>;

      using LHSType =
        ShapeFunctionBase<
          Jacobian<ShapeFunction<LHSDerived, LHSFESType, TrialSpace>>>;

      using LHSOperandType =
        ShapeFunction<LHSDerived, LHSFESType, TrialSpace>;

      using LHSOperandRangeType =
        typename FormLanguage::Traits<LHSOperandType>::RangeType;

      using RHSType =
        ShapeFunctionBase<
          Jacobian<ShapeFunction<RHSDerived, RHSFESType, TestSpace>>>;

      using RHSOperandType =
        ShapeFunction<RHSDerived, RHSFESType, TestSpace>;

      using RHSOperandRangeType =
        typename FormLanguage::Traits<RHSOperandType>::RangeType;

      using IntegrandType = Dot<LHSType, RHSType>;

      using IntegrandRangeType = typename FormLanguage::Traits<IntegrandType>::RangeType;

      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<LHSOperandRangeType, Math::Vector<ScalarType>>);

      static_assert(std::is_same_v<RHSOperandRangeType, Math::Vector<ScalarType>>);

      constexpr
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
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

      inline
      const Geometry::Polytope& getPolytope() const final override
      {
        return m_polytope.value().get();
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;
        const auto& trans = polytope.getTransformation();
        m_qf.emplace(polytope.getGeometry());
        assert(m_qf->getSize() == 1);
        m_p.emplace(polytope, trans, std::cref(m_qf->getPoint(0)));
        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        const auto& polytope = getPolytope();
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();
        const auto& qf = *m_qf;
        const auto& p = m_p.value();
        const auto& w = qf.getWeight(0);
        const auto& rc = qf.getPoint(0);
        if (trialfes == testfes)
        {
          const auto& fe = trialfes.getFiniteElement(d, idx);
          if (tr == te)
          {
            fe.getJacobian(tr)(m_jac1, rc);
            return w * p.getDistortion() * (m_jac1 * p.getJacobianInverse()).squaredNorm();
          }
          else
          {
            fe.getJacobian(tr)(m_jac1, rc);
            fe.getJacobian(te)(m_jac2, rc);
            return w * p.getDistortion() * (
                (m_jac2 * p.getJacobianInverse()).array() * (
                  m_jac1 * p.getJacobianInverse()).array()).rowwise().sum().colwise().sum().value();
          }
        }
        else
        {
          const auto& trialfe = trialfes.getFiniteElement(d, idx);
          const auto& testfe = testfes.getFiniteElement(d, idx);
          trialfe.getJacobian(tr)(m_jac1, rc);
          testfe.getJacobian(te)(m_jac2, rc);
          return w * p.getDistortion() * (
              (m_jac2 * p.getJacobianInverse()).array() * (
                m_jac1 * p.getJacobianInverse()).array()).rowwise().sum().colwise().sum().value();
        }
      }

      virtual Integrator::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      std::optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      std::optional<QF::QF1P1> m_qf;
      std::optional<Geometry::Point> m_p;

      Math::SpatialMatrix<ScalarType> m_jac1, m_jac2;
  };

  /**
   * @ingroup RodinCTAD
   */
  template <class LHSDerived, class RHSDerived, class Range, class Mesh>
  QuadratureRule(const
    Dot<
      ShapeFunctionBase<
        Jacobian<ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>>>,
      ShapeFunctionBase<
        Jacobian<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>&)
    ->
      QuadratureRule<
        Dot<
          ShapeFunctionBase<
            Jacobian<ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>>>,
          ShapeFunctionBase<
            Jacobian<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Integration of the anisotropic Frobenius inner product two
   * instances of the P1 Jacobian of ShapeFunction.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int (A \: \mathbf{J} \: u) : \mathbf{J} \: v \ dx \: ,
   * @f]
   * where @f$ u \in \mathbb{P}_1 @f$ and @f$ v \in
   * \mathbb{P}_1 @f$, and @f$ A @f$ is a coefficient function.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int (A \: \mathbf{J} \: u) : \mathbf{J} \: v \ dx :
   * \texttt{QuadratureRule}}
   * {\vdash u, v : \mathbb{P}_1}
   * @f]
   */
  template <
    class CoefficientDerived, class LHSDerived, class RHSDerived,
    class LHSRange, class RHSRange, class LHSMesh, class RHSMesh>
  class QuadratureRule<
  Dot<
    ShapeFunctionBase<
      Mult<
        FunctionBase<CoefficientDerived>,
        ShapeFunctionBase<
          Jacobian<ShapeFunction<LHSDerived, P1<LHSRange, LHSMesh>, TrialSpace>>,
        P1<LHSRange, LHSMesh>, TrialSpace>>,
      P1<LHSRange, LHSMesh>, TrialSpace>,
    ShapeFunctionBase<
      Jacobian<ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>>,
      P1<RHSRange, RHSMesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<
          Dot<ShapeFunctionBase<
            Mult<
              FunctionBase<CoefficientDerived>,
              ShapeFunctionBase<
                Jacobian<ShapeFunction<LHSDerived, P1<LHSRange, LHSMesh>, TrialSpace>>,
              P1<LHSRange, LHSMesh>, TrialSpace>>,
            P1<LHSRange, LHSMesh>, TrialSpace>,
          ShapeFunctionBase<
            Jacobian<ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>>,
            P1<RHSRange, RHSMesh>, TestSpace>>>
        ::ScalarType>
  {
    public:
      using LHSFESType = P1<LHSRange, LHSMesh>;

      using RHSFESType = P1<RHSRange, RHSMesh>;

      using CoefficientType = FunctionBase<CoefficientDerived>;

      using CoefficientRangeType =
        typename FormLanguage::Traits<CoefficientType>::RangeType;

      using MultiplicandType = ShapeFunctionBase<
        Jacobian<ShapeFunction<LHSDerived, LHSFESType, TrialSpace>>>;

      using LHSType = ShapeFunctionBase<Mult<CoefficientType, MultiplicandType>>;

      using LHSOperandType =
        ShapeFunction<LHSDerived, LHSFESType, TrialSpace>;

      using LHSOperandRangeType =
        typename FormLanguage::Traits<LHSOperandType>::RangeType;

      using RHSType =
        ShapeFunctionBase<
          Jacobian<ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>>>;

      using RHSOperandType =
        ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>;

      using RHSOperandRangeType =
        typename FormLanguage::Traits<RHSOperandType>::RangeType;

      using IntegrandType = Dot<LHSType, RHSType>;

      using IntegrandRangeType = typename FormLanguage::Traits<IntegrandType>::RangeType;

      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<LHSOperandRangeType, Math::Vector<ScalarType>>);

      static_assert(std::is_same_v<RHSOperandRangeType, Math::Vector<ScalarType>>);

      constexpr
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
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

      /**
       * @brief Gets the integrand.
       */
      inline
      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      inline
      const Geometry::Polytope& getPolytope() const final override
      {
        return m_polytope.value().get();
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;
        const auto& trans = polytope.getTransformation();
        m_qf.emplace(polytope.getGeometry());
        assert(m_qf->getSize() == 1);
        m_p.emplace(polytope, trans, std::cref(m_qf->getPoint(0)));
        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        const auto& polytope = getPolytope();
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS().getDerived();
        const auto& coeff = lhs.getLHS();
        const auto& multiplicand = lhs.getRHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = multiplicand.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();
        const auto& qf = *m_qf;
        assert(qf.getSize() == 1);
        const auto& p = m_p.value();
        const auto& w = qf.getWeight(0);
        const auto& rc = qf.getPoint(0);
        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe = testfes.getFiniteElement(d, idx);
        if constexpr (std::is_same_v<CoefficientRangeType, ScalarType>)
        {
          const ScalarType s = coeff.getValue(p);
          trialfe.getJacobian(tr)(m_jac1, rc);
          testfe.getJacobian(te)(m_jac2, rc);
          return w * p.getDistortion() * s * (
              (m_jac2 * p.getJacobianInverse()).array() * (
                m_jac1 * p.getJacobianInverse()).array()).rowwise().sum().colwise().sum().value();
        }
        else if constexpr (std::is_same_v<CoefficientRangeType, Math::Matrix<ScalarType>>)
        {
          coeff.getValue(m_mv, p);
          trialfe.getJacobian(tr)(m_jac1, rc);
          testfe.getJacobian(te)(m_jac2, rc);
          return w * p.getDistortion() * (
              (m_jac2 * p.getJacobianInverse()).array() * (
                m_mv * m_jac1 * p.getJacobianInverse()).array()).rowwise().sum().colwise().sum().value();
        }
        else
        {
          assert(false);
          return NAN;
        }
      }

      virtual Integrator::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      std::optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      std::optional<QF::QF1P1> m_qf;
      std::optional<Geometry::Point> m_p;

      Math::SpatialMatrix<ScalarType> m_jac1, m_jac2;
      Math::Matrix<ScalarType> m_mv;
  };

  /**
   * @ingroup RodinCTAD
   */
  template <class LHSFunctionDerived, class LHSDerived, class RHSDerived, class Mesh>
  QuadratureRule(const
    Dot<
      ShapeFunctionBase<
        Mult<
          FunctionBase<LHSFunctionDerived>,
          ShapeFunctionBase<Jacobian<ShapeFunction<LHSDerived, P1<Math::Vector<Real>, Mesh>, TrialSpace>>>>>,
      ShapeFunctionBase<
        Jacobian<ShapeFunction<RHSDerived, P1<Math::Vector<Real>, Mesh>, TestSpace>>>>&)
  ->
  QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Mult<
          FunctionBase<LHSFunctionDerived>,
          ShapeFunctionBase<
            Jacobian<ShapeFunction<LHSDerived, P1<Math::Vector<Real>, Mesh>, TrialSpace>>>>>,
      ShapeFunctionBase<
        Jacobian<ShapeFunction<RHSDerived, P1<Math::Vector<Real>, Mesh>, TestSpace>>>>>;
}

#endif

