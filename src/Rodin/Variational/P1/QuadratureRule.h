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

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

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
        m_weight = m_qf->getWeight(0);
        m_distortion = m_p->getDistortion();
        const auto& rc = m_qf->getPoint(0);
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& integrand = getIntegrand().getDerived();
        const auto& fes = integrand.getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, idx);
        m_basis.resize(fe.getCount());
        for (size_t i = 0; i < fe.getCount(); i++)
          m_basis[i] = fe.getBasis(i)(rc);
        return *this;
      }

      ScalarType integrate(size_t local) final override
      {
        return m_weight * m_distortion * m_basis[local];
      }

      virtual Integrator::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      std::optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      std::optional<QF::QF1P1> m_qf;
      std::optional<Geometry::Point> m_p;

      Real m_distortion;
      Real m_weight;
      std::vector<ScalarType> m_basis;
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

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

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
        m_weight = m_qf->getWeight(0);
        m_distortion = m_p->getDistortion();
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& p = *m_p;
        const auto& rc = m_qf->getPoint(0);
        const auto& integrand = getIntegrand().getDerived();
        const auto& f = integrand.getLHS();
        const auto& fes = integrand.getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, idx);
        m_basis.resize(fe.getCount());
        if constexpr (std::is_same_v<ScalarType, LHSRangeType>)
        {
          const ScalarType sv = f.getValue(p);
          for (size_t i = 0; i < fe.getCount(); i++)
            m_basis[i] = Math::dot(sv, fe.getBasis(i)(rc));
        }
        else if constexpr (std::is_same_v<Math::Vector<ScalarType>, LHSRangeType>)
        {
          f.getDerived().getValue(m_vv, p);
          for (size_t i = 0; i < fe.getCount(); i++)
          {
            fe.getBasis(i)(m_vb, rc);
            m_basis[i] = Math::dot(m_vv, m_vb);
          }
        }
        else if constexpr (std::is_same_v<Math::Matrix<ScalarType>, LHSRangeType>)
        {
          f.getValue(m_mv, p);
          for (size_t i = 0; i < fe.getCount(); i++)
          {
            fe.getBasis(i)(m_mb, rc);
            m_basis[i] = Math::dot(m_mv, m_mb);
          }
        }
        else
        {
          assert(false);
        }

        return *this;
      }

      ScalarType integrate(size_t local) final override
      {
        return m_weight * m_distortion * m_basis[local];
      }

      virtual Integrator::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:

      std::unique_ptr<IntegrandType> m_integrand;

      std::optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      std::optional<QF::QF1P1> m_qf;
      std::optional<Geometry::Point> m_p;

      Real m_weight;
      Real m_distortion;

      Math::Vector<ScalarType> m_vv;
      Math::Matrix<ScalarType> m_mv;

      Math::Vector<ScalarType> m_vb;
      Math::Matrix<ScalarType> m_mb;

      std::vector<ScalarType> m_basis;
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

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

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
        m_weight = m_qf->getWeight(0);
        m_distortion = m_p->getDistortion();
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS().getDerived();
        const auto& coeff = lhs.getLHS();
        const auto& multiplicand = lhs.getRHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();
        const auto& rc = m_qf->getPoint(0);
        const auto& p = *m_p;
        if (trialfes == testfes)
        {
          const auto& fes = trialfes;
          const auto& fe = fes.getFiniteElement(d, idx);
          m_matrix.resize(fe.getCount(), fe.getCount());
          if constexpr (std::is_same_v<CoefficientRangeType, ScalarType>)
          {
            const ScalarType csv = coeff.getValue(p);
            if constexpr (std::is_same_v<MultiplicandRangeType, ScalarType>)
            {
              m_sb1.resize(fe.getCount());
              for (size_t i = 0; i < fe.getCount(); i++)
                m_sb1[i] = fe.getBasis(i)(rc);
              for (size_t i = 0; i < fe.getCount(); i++)
                  m_matrix(i, i) = Math::dot(csv * m_sb1[i], m_sb1[i]);
              for (size_t i = 0; i < fe.getCount(); i++)
                for (size_t j = 0; j < i; j++)
                  m_matrix(i, j) = Math::dot(csv * m_sb1[j], m_sb1[i]);
              m_matrix.template triangularView<Eigen::Upper>() = m_matrix.adjoint();
            }
            else if constexpr (std::is_same_v<MultiplicandRangeType, Math::Vector<ScalarType>>)
            {
              m_vb1.resize(fe.getCount());
              for (size_t i = 0; i < fe.getCount(); i++)
                fe.getBasis(i)(m_vb1[i], rc);
              for (size_t i = 0; i < fe.getCount(); i++)
                  m_matrix(i, i) = Math::dot(csv * m_vb1[i], m_vb1[i]);
              for (size_t i = 0; i < fe.getCount(); i++)
                for (size_t j = 0; j < i; j++)
                  m_matrix(i, j) = Math::dot(csv * m_vb1[j], m_vb1[i]);
              m_matrix.template triangularView<Eigen::Upper>() = m_matrix.adjoint();
            }
            else if constexpr (std::is_same_v<MultiplicandRangeType, Math::Matrix<ScalarType>>)
            {
              m_mb1.resize(fe.getCount());
              for (size_t i = 0; i < fe.getCount(); i++)
                fe.getBasis(i)(m_mb1[i], rc);
              for (size_t i = 0; i < fe.getCount(); i++)
                  m_matrix(i, i) = Math::dot(csv * m_mb1[i], m_mb1[i]);
              for (size_t i = 0; i < fe.getCount(); i++)
                for (size_t j = 0; j < i; j++)
                  m_matrix(i, j) = Math::dot(csv * m_mb1[j], m_mb1[i]);
              m_matrix.template triangularView<Eigen::Upper>() = m_matrix.adjoint();
            }
            else
            {
              assert(false);
            }
          }
          else if constexpr (std::is_same_v<CoefficientRangeType, Math::Matrix<ScalarType>>)
          {
            coeff.getValue(m_cmv, p);
            static_assert(std::is_same_v<MultiplicandRangeType, Math::Vector<ScalarType>>);
            m_vb1.resize(fe.getCount());
            for (size_t i = 0; i < fe.getCount(); i++)
              fe.getBasis(m_vb1[i], rc);
            for (size_t i = 0; i < fe.getCount(); i++)
                m_matrix(i, i) = Math::dot(m_cmv * m_vb1[i], m_vb1[i]);
            for (size_t i = 0; i < fe.getCount(); i++)
              for (size_t j = 0; j < i; j++)
                m_matrix(i, j) = Math::dot(m_cmv * m_vb1[j], m_vb1[i]);
            m_matrix.template triangularView<Eigen::Upper>() = m_matrix.adjoint();
          }
          else
          {
            assert(false);
          }
        }
        else
        {
          const auto& trialfe = trialfes.getFiniteElement(d, idx);
          const auto& testfe = testfes.getFiniteElement(d, idx);
          m_matrix.resize(testfe.getCount(), trialfe.getCount());
          assert(false);
        }
        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_weight * m_distortion * m_matrix(te, tr);
      }

      virtual Integrator::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      std::optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      std::optional<QF::QF1P1> m_qf;
      std::optional<Geometry::Point> m_p;

      Real m_distortion;
      Real m_weight;

      Math::Vector<ScalarType> m_vv;
      Math::Matrix<ScalarType> m_cmv;

      std::vector<ScalarType> m_sb1, m_sb2;
      std::vector<Math::Vector<ScalarType>> m_vb1, m_vb2;
      std::vector<Math::Matrix<ScalarType>> m_mb1, m_mb2;

      Math::Matrix<ScalarType> m_matrix;
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

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

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
        m_weight = m_qf->getWeight(0);
        m_distortion = m_p->getDistortion();
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();
        const auto& rc = m_qf->getPoint(0);
        const auto& p = *m_p;
        if (trialfes == testfes)
        {
          const auto& fe = trialfes.getFiniteElement(d, idx);
          m_matrix.resize(fe.getCount(), fe.getCount());

          m_grad1.resize(fe.getCount());
          for (size_t i = 0; i < fe.getCount(); i++)
            fe.getGradient(i)(m_grad1[i], rc);

          for (size_t i = 0; i < fe.getCount(); i++)
            m_matrix(i, i) = (p.getJacobianInverse().transpose() * m_grad1[i]).squaredNorm();

          for (size_t i = 0; i < fe.getCount(); i++)
          {
            for (size_t j = 0; j < i; j++)
            {
              m_matrix(i, j) = Math::dot(
                  p.getJacobianInverse().transpose() * m_grad1[j],
                  p.getJacobianInverse().transpose() * m_grad1[i]);
            }
          }
          m_matrix.template triangularView<Eigen::Upper>() = m_matrix.adjoint();
        }
        else
        {
          const auto& trialfe = lhs.getFiniteElementSpace().getFiniteElement(d, idx);
          m_grad1.resize(trialfe.getCount());
          for (size_t i = 0; i < trialfe.getCount(); i++)
            trialfe.getGradient(i)(m_grad1[i], rc);

          const auto& testfe = rhs.getFiniteElementSpace().getFiniteElement(d, idx);
          m_grad2.resize(testfe.getCount());
          for (size_t i = 0; i < testfe.getCount(); i++)
            testfe.getGradient(i)(m_grad2[i], rc);

          for (size_t i = 0; i < testfe.getCount(); i++)
          {
            for (size_t j = 0; j < trialfe.getCount(); j++)
            {
              m_matrix(i, j) = Math::dot(
                  p.getJacobianInverse().transpose() * m_grad1[j],
                  p.getJacobianInverse().transpose() * m_grad2[i]);
            }
          }
        }
        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_weight * m_distortion * m_matrix(te, tr);
      }

      virtual Integrator::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      std::optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      std::optional<QF::QF1P1> m_qf;
      std::optional<Geometry::Point> m_p;

      Real m_weight;
      Real m_distortion;
      std::vector<Math::SpatialVector<ScalarType>> m_grad1, m_grad2;

      Math::Matrix<ScalarType> m_matrix;
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

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

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
        m_weight = m_qf->getWeight(0);
        m_distortion = m_p->getDistortion();
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS().getDerived();
        const auto& rhs = integrand.getRHS();
        const auto& coeff = lhs.getLHS();
        const auto& multiplicand = lhs.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();
        const auto& rc = m_qf->getPoint(0);
        const auto& p = *m_p;
        if (trialfes == testfes)
        {
          const auto& fe = trialfes.getFiniteElement(d, idx);
          m_matrix.resize(fe.getCount(), fe.getCount());
          m_grad1.resize(fe.getCount());
          for (size_t i = 0; i < fe.getCount(); i++)
            fe.getGradient(i)(m_grad1[i], rc);
          if constexpr (std::is_same_v<CoefficientRangeType, ScalarType>)
          {
            const ScalarType csv = coeff.getValue(p);
            for (size_t i = 0; i < fe.getCount(); i++)
            {
              m_matrix(i, i) =
                Math::conj(csv) * (p.getJacobianInverse().transpose() * m_grad1[i]).squaredNorm();
            }
            for (size_t i = 0; i < fe.getCount(); i++)
            {
              for (size_t j = 0; j < i; j++)
              {
                m_matrix(i, j) = Math::dot(
                    csv * p.getJacobianInverse().transpose() * m_grad1[j],
                    p.getJacobianInverse().transpose() * m_grad1[i]);
              }
            }
            m_matrix.template triangularView<Eigen::Upper>() = m_matrix.adjoint();
          }
          else if constexpr (std::is_same_v<CoefficientRangeType, Math::Matrix<ScalarType>>)
          {
            coeff.getValue(m_cmv, p);
            for (size_t i = 0; i < fe.getCount(); i++)
            {
              m_matrix(i, i) = Math::dot(
                  m_cmv * p.getJacobianInverse().transpose() * m_grad1[i],
                  p.getJacobianInverse().transpose() * m_grad1[i]);
            }
            for (size_t i = 0; i < fe.getCount(); i++)
            {
              for (size_t j = 0; j < i; j++)
              {
                m_matrix(i, j) = Math::dot(
                    m_cmv * p.getJacobianInverse().transpose() * m_grad1[j],
                    p.getJacobianInverse().transpose() * m_grad1[i]);
              }
            }
            m_matrix.template triangularView<Eigen::Upper>() = m_matrix.adjoint();
          }
          else
          {
            assert(false);
          }
        }
        else
        {
          assert(false);
        }
        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_weight * m_distortion * m_matrix(te, tr);
      }

      virtual Integrator::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      std::optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      std::optional<QF::QF1P1> m_qf;
      std::optional<Geometry::Point> m_p;

      Real m_weight;
      Real m_distortion;

      Math::Matrix<ScalarType> m_cmv;

      std::vector<Math::SpatialVector<ScalarType>> m_grad1, m_grad2;

      Math::Matrix<ScalarType> m_matrix;
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

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

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
        m_weight = m_qf->getWeight(0);
        m_distortion = m_p->getDistortion();
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();
        const auto& rc = m_qf->getPoint(0);
        const auto& p = *m_p;
        if (trialfes == testfes)
        {
          const auto& fe = trialfes.getFiniteElement(d, idx);
          m_matrix.resize(fe.getCount(), fe.getCount());

          m_jac1.resize(fe.getCount());
          for (size_t i = 0; i < fe.getCount(); i++)
            fe.getJacobian(i)(m_jac1[i], rc);

          for (size_t i = 0; i < fe.getCount(); i++)
            m_matrix(i, i) = (m_jac1[i] * p.getJacobianInverse()).squaredNorm();

          for (size_t i = 0; i < fe.getCount(); i++)
          {
            for (size_t j = 0; j < fe.getCount(); j++)
            {
              m_matrix(i, j) = Math::dot(
                  m_jac1[j] * p.getJacobianInverse(), m_jac1[i] * p.getJacobianInverse());
            }
          }
          m_matrix.template triangularView<Eigen::Upper>() = m_matrix.adjoint();
        }
        else
        {
          const auto& trialfe = lhs.getFiniteElementSpace().getFiniteElement(d, idx);
          m_jac1.resize(trialfe.getCount());
          for (size_t i = 0; i < trialfe.getCount(); i++)
            trialfe.getJacobian(i)(m_jac1[i], rc);

          const auto& testfe = rhs.getFiniteElementSpace().getFiniteElement(d, idx);
          m_jac2.resize(testfe.getCount());
          for (size_t i = 0; i < testfe.getCount(); i++)
            testfe.getJacobian(i)(m_jac2[i], rc);

          for (size_t i = 0; i < testfe.getCount(); i++)
          {
            for (size_t j = 0; j < trialfe.getCount(); j++)
            {
              m_matrix(i, j) = Math::dot(
                  m_jac1[j] * p.getJacobianInverse(), m_jac2[i] * p.getJacobianInverse());
            }
          }
        }
        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_weight * m_distortion * m_matrix(te, tr);
      }

      virtual Integrator::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      std::optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      std::optional<QF::QF1P1> m_qf;
      std::optional<Geometry::Point> m_p;

      Real m_distortion;
      Real m_weight;

      std::vector<Math::SpatialMatrix<ScalarType>> m_jac1, m_jac2;

      Math::Matrix<ScalarType> m_matrix;
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
      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

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
        m_weight = m_qf->getWeight(0);
        m_distortion = m_p->getDistortion();
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS().getDerived();
        const auto& rhs = integrand.getRHS();
        const auto& coeff = lhs.getLHS();
        const auto& multiplicand = lhs.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();
        const auto& rc = m_qf->getPoint(0);
        const auto& p = *m_p;
        if (trialfes == testfes)
        {
          const auto& fe = trialfes.getFiniteElement(d, idx);
          m_matrix.resize(fe.getCount(), fe.getCount());
          m_jac1.resize(fe.getCount());
          for (size_t i = 0; i < fe.getCount(); i++)
            fe.getJacobian(i)(m_jac1[i], rc);
          if constexpr (std::is_same_v<CoefficientRangeType, ScalarType>)
          {
            const ScalarType csv = coeff.getValue(p);
            for (size_t i = 0; i < fe.getCount(); i++)
              m_matrix(i, i) = Math::conj(csv) * (m_jac1[i] * p.getJacobianInverse()).squaredNorm();
            for (size_t i = 0; i < fe.getCount(); i++)
            {
              for (size_t j = 0; j < i; j++)
              {
                m_matrix(i, j) = Math::dot(
                    csv * m_jac1[j] * p.getJacobianInverse(),
                    m_jac1[i] * p.getJacobianInverse());
              }
            }
            m_matrix.template triangularView<Eigen::Upper>() = m_matrix.adjoint();
          }
          else if constexpr (std::is_same_v<CoefficientRangeType, Math::Matrix<ScalarType>>)
          {
            coeff.getValue(m_cmv, p);
            for (size_t i = 0; i < fe.getCount(); i++)
            {
              m_matrix(i, i) = Math::dot(
                  m_cmv * m_jac1[i] * p.getJacobianInverse(),
                  m_jac2[i] * p.getJacobianInverse());
            }
            for (size_t i = 0; i < fe.getCount(); i++)
            {
              for (size_t j = 0; j < i; j++)
              {
                m_matrix(i, j) = Math::dot(
                    m_cmv * m_jac1[j] * p.getJacobianInverse(),
                    m_jac2[i] * p.getJacobianInverse());
              }
            }
            m_matrix.template triangularView<Eigen::Upper>() = m_matrix.adjoint();
          }
          else
          {
            assert(false);
          }
        }
        else
        {
          assert(false);
        }
        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_weight * m_distortion * m_matrix(te, tr);
      }

      virtual Integrator::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      std::optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      std::optional<QF::QF1P1> m_qf;
      std::optional<Geometry::Point> m_p;

      Real m_weight;
      Real m_distortion;

      Math::Matrix<ScalarType> m_cmv;

      std::vector<Math::SpatialMatrix<ScalarType>> m_jac1, m_jac2;

      Math::Matrix<ScalarType> m_matrix;
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

  template <class Kernel, class Range, class Mesh, class LHSDerived, class RHSDerived>
  class QuadratureRule<
    Dot<
      Potential<Kernel, ShapeFunctionBase<LHSDerived, P1<Range, Mesh>, TrialSpace>>,
      ShapeFunctionBase<RHSDerived, P1<Range, Mesh>, TestSpace>>>
        : public GlobalBilinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;

      using KernelType = Kernel;

      using TrialFESType = P1<Range, Mesh>;

      using TestFESType = P1<Range, Mesh>;

      using LHSType =
        Potential<KernelType, ShapeFunctionBase<LHSDerived, TrialFESType, TrialSpace>>;

      using RHSType =
        ShapeFunctionBase<RHSDerived, TestFESType, TestSpace>;

      using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;

      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

      using IntegrandType = Dot<LHSType, RHSType>;

      using Parent = GlobalBilinearFormIntegratorBase;

      static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);

      constexpr
      QuadratureRule(const LHSType& lhs, const RHSType& rhs)
        : QuadratureRule(Dot(lhs, rhs))
      {}

      constexpr
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getOperand().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qfs(Geometry::Polytope::Type::Segment)
      {}

      constexpr
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qfs(other.m_qfs)
      {}

      constexpr
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qfs(std::move(other.m_qfs))
      {}

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& trp, const Geometry::Polytope& tep)
      {
        m_trp = trp;
        m_tep = tep;
        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = lhs.getOperand().getFiniteElementSpace();
        const auto& trialfe = trialfes.getFiniteElement(trp.getDimension(), trp.getIndex());
        const auto& testfes = rhs.getFiniteElementSpace();
        const auto& testfe = testfes.getFiniteElement(tep.getDimension(), tep.getIndex());
        const auto& kernel = lhs.getKernel();
        if (trp == tep)
        {
          const auto& polytope = trp;
          const auto& trans = trp.getTransformation();
          switch (polytope.getGeometry())
          {
            case Geometry::Polytope::Type::Point:
            case Geometry::Polytope::Type::Segment:
            case Geometry::Polytope::Type::Tetrahedron:
            {
              assert(false);
              break;
            }
            case Geometry::Polytope::Type::Quadrilateral:
            {
              assert(false);
              break;
            }
            case Geometry::Polytope::Type::Triangle:
            {

              Math::SpatialVector<ScalarType> rx0(2), rz0(2),
                                              rx1(2), rz1(2),
                                              rx2(2), rz2(2),
                                              rx3(2), rz3(2),
                                              rx4(2), rz4(2),
                                              rx5(2), rz5(2);

              assert(m_qfs.getSize() == 1);
              const auto r = m_qfs.getPoint(0).value();
              const auto w = m_qfs.getWeight(0);

              const ScalarType xi = 1 - r;
              const ScalarType eta1 = r;
              const ScalarType eta2 = r;
              const ScalarType eta3 = r;

              rx0 << xi, xi * (1 - eta1 + eta1 * eta2);
              rz1 << xi, xi * (1 - eta1 + eta1 * eta2);
              rz2 << xi * (1 - eta1 * eta2), xi * eta1 * (1 - eta2);
              rx3 << xi * (1 - eta1 * eta2), xi * eta1 * (1 - eta2);
              rz4 << xi, xi * eta1 * (1 - eta2);
              rx5 << xi, xi * eta1 * (1 - eta2);

              rz0 << xi * (1 - eta1 * eta2 * eta3), xi * (1 - eta1);
              rx1 << xi * (1 - eta1 * eta2 * eta3), xi * (1 - eta1);
              rx2 << xi, xi * eta1 * (1 - eta2 + eta2 * eta3);
              rz3 << xi, xi * eta1 * (1 - eta2 + eta2 * eta3);
              rx4 << xi * (1 - eta1 * eta2 * eta3), xi * eta1 * (1 - eta2 * eta3);
              rz5 << xi * (1 - eta1 * eta2 * eta3), xi * eta1 * (1 - eta2 * eta3);

              const Geometry::Point x0(polytope, trans, std::cref(rx0));
              const Geometry::Point x1(polytope, trans, std::cref(rx1));
              const Geometry::Point x2(polytope, trans, std::cref(rx2));
              const Geometry::Point x3(polytope, trans, std::cref(rx3));
              const Geometry::Point x4(polytope, trans, std::cref(rx4));
              const Geometry::Point x5(polytope, trans, std::cref(rx5));

              const Geometry::Point z0(polytope, trans, std::cref(rz0));
              const Geometry::Point z1(polytope, trans, std::cref(rz1));
              const Geometry::Point z2(polytope, trans, std::cref(rz2));
              const Geometry::Point z3(polytope, trans, std::cref(rz3));
              const Geometry::Point z4(polytope, trans, std::cref(rz4));
              const Geometry::Point z5(polytope, trans, std::cref(rz5));

              m_distortion = xi * xi * xi * eta1 * eta1 * eta2;
              m_weight = w * w * w * w;
              m_matrix.resize(testfe.getCount(), trialfe.getCount());

              Real s0 = x0.getDistortion() * z0.getDistortion();
              Real s1 = x1.getDistortion() * z1.getDistortion();
              Real s2 = x2.getDistortion() * z2.getDistortion();
              Real s3 = x3.getDistortion() * z3.getDistortion();
              Real s4 = x4.getDistortion() * z4.getDistortion();
              Real s5 = x5.getDistortion() * z5.getDistortion();

              assert(std::isfinite(s0));
              assert(std::isfinite(s1));
              assert(std::isfinite(s2));
              assert(std::isfinite(s3));
              assert(std::isfinite(s4));
              assert(std::isfinite(s5));

              if constexpr (std::is_same_v<LHSRangeType, ScalarType>)
              {
                for (size_t l = 0; l < testfe.getCount(); l++)
                {
                  const auto& teb = testfe.getBasis(l);
                  for (size_t m = 0; m < trialfe.getCount(); m++)
                  {
                    const auto& trb = trialfe.getBasis(m);
                    m_matrix(l, m) = s0 * kernel(x0, z0) * trb(rx0) * teb(rz0);
                    m_matrix(l, m) += s1 * kernel(x1, z1) * trb(rx1) * teb(rz1);
                    m_matrix(l, m) += s2 * kernel(x2, z2) * trb(rx2) * teb(rz2);
                    m_matrix(l, m) += s3 * kernel(x3, z3) * trb(rx3) * teb(rz3);
                    m_matrix(l, m) += s4 * kernel(x4, z4) * trb(rx4) * teb(rz4);
                    m_matrix(l, m) += s5 * kernel(x5, z5) * trb(rx5) * teb(rz5);
                  }
                }
              }
              else if constexpr (std::is_same_v<LHSRangeType, Math::Vector<ScalarType>>)
              {
                kernel(m_k0, x0, z0);
                kernel(m_k1, x1, z1);
                kernel(m_k2, x2, z2);
                kernel(m_k3, x3, z3);
                kernel(m_k4, x4, z4);
                kernel(m_k5, x5, z5);

                for (size_t l = 0; l < testfe.getCount(); l++)
                {
                  const auto& teb = testfe.getBasis(l);
                  for (size_t m = 0; m < trialfe.getCount(); m++)
                  {
                    const auto& trb = trialfe.getBasis(m);
                    trb(m_trv, rx0);
                    teb(m_tev, rz0);
                    m_matrix(l, m) = s0 * (m_k0 * m_trv).dot(m_tev);

                    trb(m_trv, rx1);
                    teb(m_tev, rz1);
                    m_matrix(l, m) += s1 * (m_k1 * m_trv).dot(m_tev);

                    trb(m_trv, rx2);
                    teb(m_tev, rz2);
                    m_matrix(l, m) += s2 * (m_k2 * m_trv).dot(m_tev);

                    trb(m_trv, rx3);
                    teb(m_tev, rz3);
                    m_matrix(l, m) += s3 * (m_k3 * m_trv).dot(m_tev);

                    trb(m_trv, rx4);
                    teb(m_tev, rz4);
                    m_matrix(l, m) += s4 * (m_k4 * m_trv).dot(m_tev);

                    trb(m_trv, rx5);
                    teb(m_tev, rz5);
                    m_matrix(l, m) += s5 * (m_k5 * m_trv).dot(m_tev);
                  }
                }
              }
              else
              {
                assert(false);
              }
              break;
            }
          }
        }
        else
        {
          m_qftr.emplace(trp.getGeometry());
          assert(m_qftr.getSize() == 1);
          m_qfte.emplace(tep.getGeometry());
          assert(m_qfte.getSize() == 1);
          const auto& trptrans = trp.getTransformation();
          const auto& teptrans = tep.getTransformation();
          const auto& rx = m_qftr->getPoint(0);
          const auto& ry = m_qfte->getPoint(0);
          const Geometry::Point x(trp, trptrans, std::cref(rx));
          const Geometry::Point y(tep, teptrans, std::cref(ry));
          m_distortion = x.getDistortion() * y.getDistortion();
          m_weight = m_qftr->getWeight(0) * m_qfte->getWeight(0);
          m_matrix.resize(testfe.getCount(), trialfe.getCount());
          if constexpr (std::is_same_v<LHSRangeType, ScalarType>)
          {
            m_sk = kernel(x, y);
            for (size_t l = 0; l < testfe.getCount(); l++)
            {
              const ScalarType teb = testfe.getBasis(l)(ry);
              for (size_t m = 0; m < trialfe.getCount(); m++)
              {
                const ScalarType trb = trialfe.getBasis(m)(rx);
                m_matrix(l, m) = m_sk * trb * teb;
              }
            }
          }
          else if constexpr (std::is_same_v<LHSRangeType, Math::Vector<ScalarType>>)
          {
            kernel(m_mk, x, y);
            for (size_t l = 0; l < testfe.getCount(); l++)
            {
              const auto teb = testfe.getBasis(l)(ry);
              for (size_t m = 0; m < trialfe.getCount(); m++)
              {
                const auto trb = trialfe.getBasis(m)(rx);
                m_matrix(l, m) = (m_mk * trb).dot(teb);
              }
            }
          }
          else
          {
            assert(false);
          }
        }
        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) override
      {
        return m_distortion * m_weight * m_matrix(te, tr);
      }

      Region getTrialRegion() const override
      {
        return getIntegrand().getLHS().getRegion();
      }

      virtual Region getTestRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      std::optional<std::reference_wrapper<const Geometry::Polytope>> m_trp;
      std::optional<std::reference_wrapper<const Geometry::Polytope>> m_tep;

      const QF::QF1P1 m_qfs;
      std::optional<QF::QF1P1> m_qftr;
      std::optional<QF::QF1P1> m_qfte;

      Real m_weight;
      Real m_distortion;

      ScalarType m_sk;
      Math::Matrix<ScalarType> m_mk;

      Math::Vector<ScalarType> m_trv, m_tev;
      Math::Matrix<ScalarType> m_k0, m_k1, m_k2, m_k3, m_k4, m_k5;

      Math::Matrix<ScalarType> m_matrix;
  };
}

#endif

