#ifndef RODIN_VARIATIONAL_H1_QUADRATURERULE_H
#define RODIN_VARIATIONAL_H1_QUADRATURERULE_H

#include "Rodin/QF/GenericPolytopeQuadrature.h"

#include "Rodin/FormLanguage/Traits.h"
#include "Rodin/Math/Common.h"
#include "Rodin/Variational/ShapeFunction.h"

#include "H1.h"
#include "H1Element.h"

namespace Rodin::Variational
{
  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Specialization for @f$\int v \ dx@f$ with an H1 test shape function.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int v \ dx \: ,
   * @f]
   * where @f$ v @f$ is an H1 shape function of order @f$K@f$.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int v \ dx : \texttt{QuadratureRule}}
   * {\vdash v : \texttt{H1}<K>}
   * @f]
   */
  template <size_t K, class NestedDerived, class Scalar, class Mesh>
  class QuadratureRule<
    ShapeFunctionBase<
      ShapeFunction<NestedDerived, H1<K, Scalar, Mesh>, TestSpace>,
      H1<K, Scalar, Mesh>, TestSpace>>
    : public LinearFormIntegratorBase<
        typename FormLanguage::Traits<
          ShapeFunctionBase<
            ShapeFunction<NestedDerived, H1<K, Scalar, Mesh>, TestSpace>,
            H1<K, Scalar, Mesh>, TestSpace>>::ScalarType>
  {
    public:
      using FESType = H1<K, Scalar, Mesh>;
      using IntegrandType =
        ShapeFunctionBase<ShapeFunction<NestedDerived, FESType, TestSpace>, FESType, TestSpace>;
      using IntegrandRangeType = typename FormLanguage::Traits<IntegrandType>::RangeType;
      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;
      using Parent = LinearFormIntegratorBase<ScalarType>;

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_ps(std::move(other.m_ps)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_vec(std::move(other.m_vec))
      {}

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& fes = integrand.getFiniteElementSpace();
        const auto& fe  = fes.getFiniteElement(d, idx);

        const size_t order =
          integrand.getOrder(polytope).value_or(fe.getOrder());

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::GenericPolytopeQuadrature::get(order, geometry);

          m_ps.clear();
          m_ps.reserve(m_qf->getSize());
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps.emplace_back(polytope, m_qf->getPoint(qp));
        }
        else
        {
          assert(m_qf);
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps[qp].setPolytope(polytope);
        }

        const size_t nte = integrand.getDOFs(polytope);

        const H1Element<K, ScalarType> scalarFE(geometry);
        const auto& tab = scalarFE.getTabulation(*m_qf);

        m_vec.resize(static_cast<Eigen::Index>(nte));
        m_vec.setZero();

        for (size_t qp = 0; qp < m_ps.size(); ++qp)
        {
          const auto& p = m_ps[qp];
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          if constexpr (std::is_same_v<IntegrandRangeType, ScalarType>)
          {
            for (size_t local = 0; local < nte; ++local)
              m_vec(local) += wdet * tab.getBasis(qp, local);
          }
          else if constexpr (std::is_same_v<IntegrandRangeType, Math::Vector<ScalarType>>)
          {
            const size_t vdim = fes.getVectorDimension();
            assert(nte == scalarFE.getCount() * vdim);

            for (size_t local = 0; local < nte; ++local)
            {
              const size_t scalarLocal = local / vdim;
              m_vec(local) += wdet * tab.getBasis(qp, scalarLocal);
            }
          }
          else
          {
            static_assert(
              std::is_same_v<IntegrandRangeType, ScalarType>
              || std::is_same_v<IntegrandRangeType, Math::Vector<ScalarType>>,
              "Unsupported H1 Integral(v) range type. Expected scalar or vector-valued shape function.");
          }
        }

        return *this;
      }

      inline ScalarType integrate(size_t local) final override
      {
        return m_vec(local);
      }

      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      std::vector<Geometry::Point> m_ps;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      Math::Vector<ScalarType> m_vec;
  };

  template <size_t K, class NestedDerived, class Scalar, class Mesh>
  QuadratureRule(
    const ShapeFunctionBase<
      ShapeFunction<NestedDerived, H1<K, Scalar, Mesh>, TestSpace>,
      H1<K, Scalar, Mesh>, TestSpace>&)
    -> QuadratureRule<
        ShapeFunctionBase<
          ShapeFunction<NestedDerived, H1<K, Scalar, Mesh>, TestSpace>,
          H1<K, Scalar, Mesh>, TestSpace>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Specialization for @f$\int f \cdot v \ dx@f$ with an H1 test shape function.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int f \cdot v \ dx \: ,
   * @f]
   * where @f$ v @f$ is an H1 shape function of order @f$K@f$.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int f \cdot v \ dx : \texttt{QuadratureRule}}
   * {\vdash v : \texttt{H1}<K>}
   * @f]
   */
  template <size_t K, class LHSDerived, class RHSDerived, class Scalar, class Mesh>
  class QuadratureRule<
    ShapeFunctionBase<
      Dot<
        FunctionBase<LHSDerived>,
        ShapeFunctionBase<
          ShapeFunction<RHSDerived, H1<K, Scalar, Mesh>, TestSpace>,
          H1<K, Scalar, Mesh>, TestSpace>>,
      H1<K, Scalar, Mesh>, TestSpace>>
    : public LinearFormIntegratorBase<
        typename FormLanguage::Traits<
          ShapeFunctionBase<
            Dot<
              FunctionBase<LHSDerived>,
              ShapeFunctionBase<
                ShapeFunction<RHSDerived, H1<K, Scalar, Mesh>, TestSpace>,
                H1<K, Scalar, Mesh>, TestSpace>>,
            H1<K, Scalar, Mesh>, TestSpace>>::ScalarType>
  {
    public:
      using FESType = H1<K, Scalar, Mesh>;

      using LHSType = FunctionBase<LHSDerived>;

      using RHSType =
        ShapeFunctionBase<ShapeFunction<RHSDerived, FESType, TestSpace>, FESType, TestSpace>;

      using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;
      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

      using IntegrandType =
        ShapeFunctionBase<Dot<LHSType, RHSType>, FESType, TestSpace>;

      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent = LinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_ps(std::move(other.m_ps)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_vec(std::move(other.m_vec))
      {}

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& f   = integrand.getDerived().getLHS();
        const auto& fes = integrand.getFiniteElementSpace();
        const auto& fe  = fes.getFiniteElement(d, idx);

        const size_t order =
          integrand.getOrder(polytope).value_or(fe.getOrder());

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::GenericPolytopeQuadrature::get(order, geometry);

          m_ps.clear();
          m_ps.reserve(m_qf->getSize());
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps.emplace_back(polytope, m_qf->getPoint(qp));
        }
        else
        {
          assert(m_qf);
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps[qp].setPolytope(polytope);
        }

        const size_t nte = integrand.getDOFs(polytope);

        const H1Element<K, ScalarType> scalarFE(geometry);
        const auto& tab = scalarFE.getTabulation(*m_qf);

        m_vec.resize(static_cast<Eigen::Index>(nte));
        m_vec.setZero();

        for (size_t qp = 0; qp < m_ps.size(); ++qp)
        {
          const auto& p = m_ps[qp];
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const LHSRangeType fval = f(p);

          if constexpr (std::is_same_v<RHSRangeType, ScalarType>)
          {
            for (size_t local = 0; local < nte; ++local)
              m_vec(local) += wdet * fval * tab.getBasis(qp, local);
          }
          else if constexpr (std::is_same_v<RHSRangeType, Math::Vector<ScalarType>>)
          {
            const size_t vdim = fes.getVectorDimension();
            assert(nte == scalarFE.getCount() * vdim);

            for (size_t local = 0; local < nte; ++local)
            {
              const size_t scalarLocal = local / vdim;
              const size_t comp = local % vdim;
              m_vec(local) += wdet * fval.coeff(comp) * tab.getBasis(qp, scalarLocal);
            }
          }
          else
          {
            static_assert(
              std::is_same_v<RHSRangeType, ScalarType>
              || std::is_same_v<RHSRangeType, Math::Vector<ScalarType>>,
              "Unsupported H1 Integral(f,v) RHS range type. Expected scalar or vector-valued shape function.");
          }
        }

        return *this;
      }

      inline ScalarType integrate(size_t local) final override
      {
        return m_vec(local);
      }

      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      std::vector<Geometry::Point> m_ps;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      Math::Vector<ScalarType> m_vec;
  };

  template <size_t K, class LHSDerived, class RHSDerived, class Scalar, class Mesh>
  QuadratureRule(
    const ShapeFunctionBase<
      Dot<
        FunctionBase<LHSDerived>,
        ShapeFunctionBase<
          ShapeFunction<RHSDerived, H1<K, Scalar, Mesh>, TestSpace>,
          H1<K, Scalar, Mesh>, TestSpace>>,
      H1<K, Scalar, Mesh>, TestSpace>&)
    -> QuadratureRule<
        ShapeFunctionBase<
          Dot<
            FunctionBase<LHSDerived>,
            ShapeFunctionBase<
              ShapeFunction<RHSDerived, H1<K, Scalar, Mesh>, TestSpace>,
              H1<K, Scalar, Mesh>, TestSpace>>,
          H1<K, Scalar, Mesh>, TestSpace>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Specialization for @f$\int u \cdot v \ dx@f$ with H1 trial and test shape functions.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int u \cdot v \ dx \: ,
   * @f]
   * where @f$ u @f$ and @f$ v @f$ are H1 shape functions.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int u \cdot v \ dx : \texttt{QuadratureRule}}
   * {\vdash u : \texttt{H1}<K_{\mathrm{trial}}>, \ \vdash v : \texttt{H1}<K_{\mathrm{test}}>}
   * @f]
   */
  template <
    size_t KTrial, size_t KTest,
    class LHSDerived, class RHSDerived,
    class Scalar, class Mesh>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<
        ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>,
        H1<KTrial, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
        H1<KTest, Scalar, Mesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<
          Dot<
            ShapeFunctionBase<
              ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>,
              H1<KTrial, Scalar, Mesh>, TrialSpace>,
            ShapeFunctionBase<
              ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
              H1<KTest, Scalar, Mesh>, TestSpace>>>::ScalarType>
  {
    public:
      using TrialFESType = H1<KTrial, Scalar, Mesh>;
      using TestFESType  = H1<KTest, Scalar, Mesh>;

      using LHSType =
        ShapeFunctionBase<ShapeFunction<LHSDerived, TrialFESType, TrialSpace>, TrialFESType, TrialSpace>;

      using RHSType =
        ShapeFunctionBase<ShapeFunction<RHSDerived, TestFESType, TestSpace>, TestFESType, TestSpace>;

      using IntegrandType = Dot<LHSType, RHSType>;
      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_ps(std::move(other.m_ps)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_mat(std::move(other.m_mat))
      {}

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();

        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes  = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes .getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe .getOrder();
        const size_t order = k_tr + k_te;

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::GenericPolytopeQuadrature::get(order, geometry);

          m_ps.clear();
          m_ps.reserve(m_qf->getSize());
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps.emplace_back(polytope, m_qf->getPoint(qp));
        }
        else
        {
          assert(m_qf);
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps[qp].setPolytope(polytope);
        }

        const size_t ntr = lhs.getDOFs(polytope);
        const size_t nte = rhs.getDOFs(polytope);

        const H1Element<KTrial, ScalarType> trialScalarFE(geometry);
        const H1Element<KTest,  ScalarType> testScalarFE(geometry);
        const size_t scalarCountTr = trialScalarFE.getCount();
        const size_t scalarCountTe = testScalarFE.getCount();

        assert(scalarCountTr > 0 && ntr % scalarCountTr == 0);
        assert(scalarCountTe > 0 && nte % scalarCountTe == 0);
        const size_t vdim = ntr / scalarCountTr;
        assert(vdim == nte / scalarCountTe);

        const bool symmetric =
          (&trialfes.getMesh() == &testfes.getMesh()) && (ntr == nte);

        const auto& trTab = trialScalarFE.getTabulation(*m_qf);
        const auto& teTab = testScalarFE .getTabulation(*m_qf);

        m_mat.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_mat.setZero();
        ScalarType* A = m_mat.data();

        for (size_t qp = 0; qp < m_ps.size(); ++qp)
        {
          const auto& p = m_ps[qp];
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          if (symmetric)
          {
            for (size_t ib = 0; ib < scalarCountTe; ++ib)
            {
              const ScalarType phi_te = teTab.getBasis(qp, ib);

              {
                const ScalarType kii = wdet * phi_te * trTab.getBasis(qp, ib);
                for (size_t c = 0; c < vdim; ++c)
                  A[(ib * vdim + c) * ntr + (ib * vdim + c)] += kii;
              }

              for (size_t ia = 0; ia < ib; ++ia)
              {
                const ScalarType kij = wdet * phi_te * trTab.getBasis(qp, ia);
                for (size_t c = 0; c < vdim; ++c)
                  A[(ib * vdim + c) * ntr + (ia * vdim + c)] += kij;
              }
            }
          }
          else
          {
            for (size_t ib = 0; ib < scalarCountTe; ++ib)
            {
              const ScalarType phi_te = teTab.getBasis(qp, ib);
              for (size_t ia = 0; ia < scalarCountTr; ++ia)
              {
                const ScalarType kij = wdet * phi_te * trTab.getBasis(qp, ia);
                for (size_t c = 0; c < vdim; ++c)
                  A[(ib * vdim + c) * ntr + (ia * vdim + c)] += kij;
              }
            }
          }
        }

        if (symmetric)
        {
          m_mat.template triangularView<Eigen::Upper>() =
            m_mat.transpose().template triangularView<Eigen::Upper>();
        }

        return *this;
      }

      inline ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_mat(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      std::vector<Geometry::Point> m_ps;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_mat;
  };

  template <size_t KTrial, size_t KTest, class LHSDerived, class RHSDerived, class Scalar, class Mesh>
  QuadratureRule(
    const Dot<
      ShapeFunctionBase<
        ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>,
        H1<KTrial, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
        H1<KTest, Scalar, Mesh>, TestSpace>>&)
    -> QuadratureRule<
        Dot<
          ShapeFunctionBase<
            ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>,
            H1<KTrial, Scalar, Mesh>, TrialSpace>,
          ShapeFunctionBase<
            ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
            H1<KTest, Scalar, Mesh>, TestSpace>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Specialization for @f$\int (A u) \cdot v \ dx@f$ with H1 trial and test shape functions.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int (A u) \cdot v \ dx \: ,
   * @f]
   * where @f$u@f$ and @f$v@f$ are H1 shape functions and @f$A@f$ is a coefficient function.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int (A u) \cdot v \ dx : \texttt{QuadratureRule}}
   * {\vdash u : \texttt{H1}<K_{\mathrm{trial}}>, \ \vdash v : \texttt{H1}<K_{\mathrm{test}}>}
   * @f]
   */
  template <
    size_t KTrial, size_t KTest,
    class CoefficientDerived, class LHSDerived, class RHSDerived,
    class Scalar, class Mesh>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Mult<
          FunctionBase<CoefficientDerived>,
          ShapeFunctionBase<
            ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>,
            H1<KTrial, Scalar, Mesh>, TrialSpace>>,
        H1<KTrial, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
        H1<KTest, Scalar, Mesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<
          Dot<
            ShapeFunctionBase<
              Mult<
                FunctionBase<CoefficientDerived>,
                ShapeFunctionBase<
                  ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>,
                  H1<KTrial, Scalar, Mesh>, TrialSpace>>,
              H1<KTrial, Scalar, Mesh>, TrialSpace>,
            ShapeFunctionBase<
              ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
              H1<KTest, Scalar, Mesh>, TestSpace>>>::ScalarType>
  {
    public:
      using TrialFESType = H1<KTrial, Scalar, Mesh>;
      using TestFESType  = H1<KTest, Scalar, Mesh>;

      using CoefficientType = FunctionBase<CoefficientDerived>;
      using CoefficientRangeType =
        typename FormLanguage::Traits<CoefficientType>::RangeType;

      using MultiplicandType =
        ShapeFunctionBase<ShapeFunction<LHSDerived, TrialFESType, TrialSpace>, TrialFESType, TrialSpace>;

      using LHSType =
        ShapeFunctionBase<Mult<CoefficientType, MultiplicandType>, TrialFESType, TrialSpace>;

      using RHSType =
        ShapeFunctionBase<ShapeFunction<RHSDerived, TestFESType, TestSpace>, TestFESType, TestSpace>;

      using IntegrandType = Dot<LHSType, RHSType>;
      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_ps(std::move(other.m_ps)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_cmv(std::move(other.m_cmv)),
          m_mat(std::move(other.m_mat))
      {}

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& lhs   = integrand.getLHS();
        const auto& coeff = lhs.getDerived().getLHS();
        const auto& rhs   = integrand.getRHS();

        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes  = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes .getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe .getOrder();
        const size_t order = k_tr + k_te;

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::GenericPolytopeQuadrature::get(order, geometry);

          m_ps.clear();
          m_ps.reserve(m_qf->getSize());
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps.emplace_back(polytope, m_qf->getPoint(qp));
        }
        else
        {
          assert(m_qf);
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps[qp].setPolytope(polytope);
        }

        const size_t ntr = lhs.getDOFs(polytope);
        const size_t nte = rhs.getDOFs(polytope);

        const H1Element<KTrial, ScalarType> trialScalarFE(geometry);
        const H1Element<KTest,  ScalarType> testScalarFE(geometry);
        const size_t scalarCountTr = trialScalarFE.getCount();
        const size_t scalarCountTe = testScalarFE.getCount();

        assert(scalarCountTr > 0 && ntr % scalarCountTr == 0);
        assert(scalarCountTe > 0 && nte % scalarCountTe == 0);
        const size_t vdim = ntr / scalarCountTr;
        assert(vdim == nte / scalarCountTe);

        const auto& trTab = trialScalarFE.getTabulation(*m_qf);
        const auto& teTab = testScalarFE .getTabulation(*m_qf);

        m_mat.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_mat.setZero();
        ScalarType* A = m_mat.data();

        for (size_t qp = 0; qp < m_ps.size(); ++qp)
        {
          const auto& p = m_ps[qp];
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          if constexpr (std::is_same_v<CoefficientRangeType, ScalarType>)
          {
            const ScalarType csv = coeff.getValue(p);
            for (size_t ib = 0; ib < scalarCountTe; ++ib)
            {
              const ScalarType phi_te = teTab.getBasis(qp, ib);
              for (size_t ia = 0; ia < scalarCountTr; ++ia)
              {
                const ScalarType kij = wdet * csv * phi_te * trTab.getBasis(qp, ia);
                for (size_t c = 0; c < vdim; ++c)
                  A[(ib * vdim + c) * ntr + (ia * vdim + c)] += kij;
              }
            }
          }
          else if constexpr (std::is_same_v<CoefficientRangeType, Math::Matrix<ScalarType>>)
          {
            coeff.getValue(m_cmv, p);
            for (size_t ib = 0; ib < scalarCountTe; ++ib)
            {
              const ScalarType phi_te = teTab.getBasis(qp, ib);
              for (size_t ia = 0; ia < scalarCountTr; ++ia)
              {
                const ScalarType basis_prod = wdet * phi_te * trTab.getBasis(qp, ia);
                for (size_t dd = 0; dd < vdim; ++dd)
                  for (size_t cc = 0; cc < vdim; ++cc)
                    A[(ib * vdim + dd) * ntr + (ia * vdim + cc)] += basis_prod * m_cmv(dd, cc);
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

      inline ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_mat(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      std::vector<Geometry::Point> m_ps;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      Math::Matrix<ScalarType> m_cmv;
      Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_mat;
  };

  template <
    size_t KTrial, size_t KTest,
    class CoefficientDerived, class LHSDerived, class RHSDerived,
    class Scalar, class Mesh>
  QuadratureRule(
    const Dot<
      ShapeFunctionBase<
        Mult<
          FunctionBase<CoefficientDerived>,
          ShapeFunctionBase<
            ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>,
            H1<KTrial, Scalar, Mesh>, TrialSpace>>,
        H1<KTrial, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
        H1<KTest, Scalar, Mesh>, TestSpace>>&)
    -> QuadratureRule<
        Dot<
          ShapeFunctionBase<
            Mult<
              FunctionBase<CoefficientDerived>,
              ShapeFunctionBase<
                ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>,
                H1<KTrial, Scalar, Mesh>, TrialSpace>>,
            H1<KTrial, Scalar, Mesh>, TrialSpace>,
          ShapeFunctionBase<
            ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
            H1<KTest, Scalar, Mesh>, TestSpace>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Specialization for @f$\int (A \nabla u) \cdot \nabla v \ dx@f$ with H1 trial and test shape functions.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int (A \nabla u) \cdot \nabla v \ dx \: ,
   * @f]
   * where @f$u@f$ and @f$v@f$ are H1 shape functions and @f$A@f$ is a coefficient function.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int (A \nabla u) \cdot \nabla v \ dx : \texttt{QuadratureRule}}
   * {\vdash u : \texttt{H1}<K_{\mathrm{trial}}>, \ \vdash v : \texttt{H1}<K_{\mathrm{test}}>}
   * @f]
   */
  template <
    size_t KTrial, size_t KTest,
    class CoefficientDerived, class LHSDerived, class RHSDerived,
    class Scalar, class Mesh>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Mult<
          FunctionBase<CoefficientDerived>,
          ShapeFunctionBase<
            Grad<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
            H1<KTrial, Scalar, Mesh>, TrialSpace>>,
        H1<KTrial, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        Grad<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
        H1<KTest, Scalar, Mesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<
          Dot<
            ShapeFunctionBase<
              Mult<
                FunctionBase<CoefficientDerived>,
                ShapeFunctionBase<
                  Grad<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
                  H1<KTrial, Scalar, Mesh>, TrialSpace>>,
              H1<KTrial, Scalar, Mesh>, TrialSpace>,
            ShapeFunctionBase<
              Grad<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
              H1<KTest, Scalar, Mesh>, TestSpace>>>::ScalarType>
  {
    public:
      using TrialFESType = H1<KTrial, Scalar, Mesh>;
      using TestFESType  = H1<KTest, Scalar, Mesh>;

      using CoefficientType = FunctionBase<CoefficientDerived>;
      using CoefficientRangeType =
        typename FormLanguage::Traits<CoefficientType>::RangeType;

      using MultiplicandType =
        ShapeFunctionBase<
          Grad<ShapeFunction<LHSDerived, TrialFESType, TrialSpace>>,
          TrialFESType, TrialSpace>;

      using LHSType =
        ShapeFunctionBase<Mult<CoefficientType, MultiplicandType>, TrialFESType, TrialSpace>;

      using RHSType =
        ShapeFunctionBase<
          Grad<ShapeFunction<RHSDerived, TestFESType, TestSpace>>,
          TestFESType, TestSpace>;

      using IntegrandType = Dot<LHSType, RHSType>;
      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_ps(std::move(other.m_ps)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_cmv(std::move(other.m_cmv)),
          m_mat(std::move(other.m_mat))
      {}

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& lhs   = integrand.getLHS();
        const auto& coeff = lhs.getDerived().getLHS();
        const auto& rhs   = integrand.getRHS();

        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes  = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes .getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe .getOrder();
        const size_t order = (k_tr == 0 || k_te == 0) ? 0 : (k_tr + k_te - 2);

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::GenericPolytopeQuadrature::get(order, geometry);

          m_ps.clear();
          m_ps.reserve(m_qf->getSize());
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps.emplace_back(polytope, m_qf->getPoint(qp));
        }
        else
        {
          assert(m_qf);
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps[qp].setPolytope(polytope);
        }

        const size_t ntr = lhs.getDOFs(polytope);
        const size_t nte = rhs.getDOFs(polytope);

        m_mat.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_mat.setZero();
        ScalarType* A_data = m_mat.data();

        const auto& trTab = trialfe.getTabulation(*m_qf);
        const auto& teTab = testfe .getTabulation(*m_qf);

        static thread_local std::vector<Math::SpatialVector<ScalarType>> Gtr;
        static thread_local std::vector<Math::SpatialVector<ScalarType>> Gte;

        if (Gtr.size() < ntr) Gtr.resize(ntr);
        if (Gte.size() < nte) Gte.resize(nte);
        for (size_t a = 0; a < ntr; ++a) Gtr[a].resize(static_cast<std::uint8_t>(d));
        for (size_t b = 0; b < nte; ++b) Gte[b].resize(static_cast<std::uint8_t>(d));

        for (size_t qp = 0; qp < m_ps.size(); ++qp)
        {
          const auto& p = m_ps[qp];
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto Jinv = p.getJacobianInverse();

          // Build physical gradients for trial and test
          if (d == 3)
          {
            const ScalarType a00 = Jinv(0,0), a10 = Jinv(1,0), a20 = Jinv(2,0);
            const ScalarType a01 = Jinv(0,1), a11 = Jinv(1,1), a21 = Jinv(2,1);
            const ScalarType a02 = Jinv(0,2), a12 = Jinv(1,2), a22 = Jinv(2,2);

            for (size_t a = 0; a < ntr; ++a)
            {
              const auto g = trTab.getGradient(qp, a);
              const ScalarType gx = g[0], gy = g[1], gz = g[2];
              Gtr[a][0] = a00*gx + a10*gy + a20*gz;
              Gtr[a][1] = a01*gx + a11*gy + a21*gz;
              Gtr[a][2] = a02*gx + a12*gy + a22*gz;
            }
            for (size_t b = 0; b < nte; ++b)
            {
              const auto g = teTab.getGradient(qp, b);
              const ScalarType gx = g[0], gy = g[1], gz = g[2];
              Gte[b][0] = a00*gx + a10*gy + a20*gz;
              Gte[b][1] = a01*gx + a11*gy + a21*gz;
              Gte[b][2] = a02*gx + a12*gy + a22*gz;
            }
          }
          else if (d == 2)
          {
            const ScalarType a00 = Jinv(0,0), a10 = Jinv(1,0);
            const ScalarType a01 = Jinv(0,1), a11 = Jinv(1,1);

            for (size_t a = 0; a < ntr; ++a)
            {
              const auto g = trTab.getGradient(qp, a);
              const ScalarType gx = g[0], gy = g[1];
              Gtr[a][0] = a00*gx + a10*gy;
              Gtr[a][1] = a01*gx + a11*gy;
            }
            for (size_t b = 0; b < nte; ++b)
            {
              const auto g = teTab.getGradient(qp, b);
              const ScalarType gx = g[0], gy = g[1];
              Gte[b][0] = a00*gx + a10*gy;
              Gte[b][1] = a01*gx + a11*gy;
            }
          }
          else if (d == 1)
          {
            const ScalarType a00 = Jinv(0,0);
            for (size_t a = 0; a < ntr; ++a)
            {
              const auto g = trTab.getGradient(qp, a);
              Gtr[a][0] = a00 * g[0];
            }
            for (size_t b = 0; b < nte; ++b)
            {
              const auto g = teTab.getGradient(qp, b);
              Gte[b][0] = a00 * g[0];
            }
          }
          else
          {
            assert(false);
          }

          if constexpr (std::is_same_v<CoefficientRangeType, ScalarType>)
          {
            const ScalarType csv = coeff.getValue(p);
            for (size_t b = 0; b < nte; ++b)
            {
              ScalarType* row = A_data + b * ntr;
              const auto& gb = Gte[b];
              for (size_t a = 0; a < ntr; ++a)
                row[a] += wdet * csv * Math::dot(gb, Gtr[a]);
            }
          }
          else if constexpr (std::is_same_v<CoefficientRangeType, Math::Matrix<ScalarType>>)
          {
            coeff.getValue(m_cmv, p);
            for (size_t b = 0; b < nte; ++b)
            {
              ScalarType* row = A_data + b * ntr;
              const auto& gb = Gte[b];
              for (size_t a = 0; a < ntr; ++a)
              {
                Math::SpatialVector<ScalarType> AGtr_a = m_cmv * Gtr[a];
                row[a] += wdet * Math::dot(gb, AGtr_a);
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

      inline ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_mat(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      std::vector<Geometry::Point> m_ps;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      Math::Matrix<ScalarType> m_cmv;
      Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_mat;
  };

  template <
    size_t KTrial, size_t KTest,
    class CoefficientDerived, class LHSDerived, class RHSDerived,
    class Scalar, class Mesh>
  QuadratureRule(
    const Dot<
      ShapeFunctionBase<
        Mult<
          FunctionBase<CoefficientDerived>,
          ShapeFunctionBase<
            Grad<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
            H1<KTrial, Scalar, Mesh>, TrialSpace>>,
        H1<KTrial, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        Grad<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
        H1<KTest, Scalar, Mesh>, TestSpace>>&)
    -> QuadratureRule<
        Dot<
          ShapeFunctionBase<
            Mult<
              FunctionBase<CoefficientDerived>,
              ShapeFunctionBase<
                Grad<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
                H1<KTrial, Scalar, Mesh>, TrialSpace>>,
            H1<KTrial, Scalar, Mesh>, TrialSpace>,
          ShapeFunctionBase<
            Grad<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
            H1<KTest, Scalar, Mesh>, TestSpace>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Specialization for @f$\int c \, (u \cdot v) \ dx@f$ with H1 trial and test shape functions.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int c \, (u \cdot v) \ dx \: ,
   * @f]
   * where @f$u@f$ and @f$v@f$ are H1 shape functions and @f$c@f$ is a scalar coefficient function.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int c \, (u \cdot v) \ dx : \texttt{QuadratureRule}}
   * {\vdash u : \texttt{H1}<K_{\mathrm{trial}}>, \ \vdash v : \texttt{H1}<K_{\mathrm{test}}>}
   * @f]
   */
  template <
    size_t KTrial, size_t KTest,
    class CoefficientDerived, class LHSDerived, class RHSDerived,
    class Scalar, class Mesh>
  class QuadratureRule<
    Mult<
      FunctionBase<CoefficientDerived>,
      Dot<
        ShapeFunctionBase<
          ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>,
          H1<KTrial, Scalar, Mesh>, TrialSpace>,
        ShapeFunctionBase<
          ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
          H1<KTest, Scalar, Mesh>, TestSpace>>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<
          Mult<
            FunctionBase<CoefficientDerived>,
            Dot<
              ShapeFunctionBase<
                ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>,
                H1<KTrial, Scalar, Mesh>, TrialSpace>,
              ShapeFunctionBase<
                ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
                H1<KTest, Scalar, Mesh>, TestSpace>>>>::ScalarType>
  {
    public:
      using TrialFESType = H1<KTrial, Scalar, Mesh>;
      using TestFESType  = H1<KTest, Scalar, Mesh>;

      using CoefficientType = FunctionBase<CoefficientDerived>;

      using LHSType =
        ShapeFunctionBase<ShapeFunction<LHSDerived, TrialFESType, TrialSpace>, TrialFESType, TrialSpace>;

      using RHSType =
        ShapeFunctionBase<ShapeFunction<RHSDerived, TestFESType, TestSpace>, TestFESType, TestSpace>;

      using InnerIntegrandType = Dot<LHSType, RHSType>;
      using IntegrandType = Mult<CoefficientType, InnerIntegrandType>;
      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getRHS().getLHS().getLeaf(), integrand.getRHS().getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_ps(std::move(other.m_ps)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_mat(std::move(other.m_mat))
      {}

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& coeff = integrand.getLHS();
        const auto& inner = integrand.getRHS();
        const auto& lhs   = inner.getLHS();
        const auto& rhs   = inner.getRHS();

        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes  = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes .getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe .getOrder();
        const size_t order = k_tr + k_te;

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::GenericPolytopeQuadrature::get(order, geometry);

          m_ps.clear();
          m_ps.reserve(m_qf->getSize());
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps.emplace_back(polytope, m_qf->getPoint(qp));
        }
        else
        {
          assert(m_qf);
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps[qp].setPolytope(polytope);
        }

        const size_t ntr = lhs.getDOFs(polytope);
        const size_t nte = rhs.getDOFs(polytope);

        const H1Element<KTrial, ScalarType> trialScalarFE(geometry);
        const H1Element<KTest,  ScalarType> testScalarFE(geometry);
        const size_t scalarCountTr = trialScalarFE.getCount();
        const size_t scalarCountTe = testScalarFE.getCount();

        assert(scalarCountTr > 0 && ntr % scalarCountTr == 0);
        assert(scalarCountTe > 0 && nte % scalarCountTe == 0);
        const size_t vdim = ntr / scalarCountTr;
        assert(vdim == nte / scalarCountTe);

        const auto& trTab = trialScalarFE.getTabulation(*m_qf);
        const auto& teTab = testScalarFE .getTabulation(*m_qf);

        m_mat.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_mat.setZero();
        ScalarType* A = m_mat.data();

        for (size_t qp = 0; qp < m_ps.size(); ++qp)
        {
          const auto& p = m_ps[qp];
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const ScalarType csv = coeff(p);

          for (size_t ib = 0; ib < scalarCountTe; ++ib)
          {
            const ScalarType phi_te = teTab.getBasis(qp, ib);
            for (size_t ia = 0; ia < scalarCountTr; ++ia)
            {
              const ScalarType kij = wdet * csv * phi_te * trTab.getBasis(qp, ia);
              for (size_t c = 0; c < vdim; ++c)
                A[(ib * vdim + c) * ntr + (ia * vdim + c)] += kij;
            }
          }
        }

        return *this;
      }

      inline ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_mat(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      std::vector<Geometry::Point> m_ps;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_mat;
  };

  template <
    size_t KTrial, size_t KTest,
    class CoefficientDerived, class LHSDerived, class RHSDerived,
    class Scalar, class Mesh>
  QuadratureRule(
    const Mult<
      FunctionBase<CoefficientDerived>,
      Dot<
        ShapeFunctionBase<
          ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>,
          H1<KTrial, Scalar, Mesh>, TrialSpace>,
        ShapeFunctionBase<
          ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
          H1<KTest, Scalar, Mesh>, TestSpace>>>&)
    -> QuadratureRule<
        Mult<
          FunctionBase<CoefficientDerived>,
          Dot<
            ShapeFunctionBase<
              ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>,
              H1<KTrial, Scalar, Mesh>, TrialSpace>,
            ShapeFunctionBase<
              ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
              H1<KTest, Scalar, Mesh>, TestSpace>>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Specialization for @f$\int \nabla \cdot u \, q \ dx@f$ with H1 trial and test shape functions.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int \nabla \cdot u \, q \ dx \: ,
   * @f]
   * where @f$u@f$ is an H1 trial shape function and @f$q@f$ is an H1 test shape function.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int \nabla \cdot u \, q \ dx : \texttt{QuadratureRule}}
   * {\vdash u : \texttt{H1}<K_{\mathrm{trial}}>, \ \vdash q : \texttt{H1}<K_{\mathrm{test}}>}
   * @f]
   */
  template <
    size_t KTrial, size_t KTest,
    class LHSDerived, class RHSDerived,
    class Scalar, class Mesh>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Div<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
        H1<KTrial, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
        H1<KTest, Scalar, Mesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<
          Dot<
            ShapeFunctionBase<
              Div<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
              H1<KTrial, Scalar, Mesh>, TrialSpace>,
            ShapeFunctionBase<
              ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
              H1<KTest, Scalar, Mesh>, TestSpace>>>::ScalarType>
  {
    public:
      using TrialFESType = H1<KTrial, Scalar, Mesh>;
      using TestFESType  = H1<KTest, Scalar, Mesh>;

      using LHSType =
        ShapeFunctionBase<
          Div<ShapeFunction<LHSDerived, TrialFESType, TrialSpace>>,
          TrialFESType, TrialSpace>;

      using RHSType =
        ShapeFunctionBase<ShapeFunction<RHSDerived, TestFESType, TestSpace>, TestFESType, TestSpace>;

      using IntegrandType = Dot<LHSType, RHSType>;
      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_ps(std::move(other.m_ps)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_mat(std::move(other.m_mat))
      {}

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();

        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes  = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes .getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe .getOrder();
        const size_t order = (k_tr == 0) ? 0 : (k_tr + k_te - 1);

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::GenericPolytopeQuadrature::get(order, geometry);

          m_ps.clear();
          m_ps.reserve(m_qf->getSize());
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps.emplace_back(polytope, m_qf->getPoint(qp));
        }
        else
        {
          assert(m_qf);
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps[qp].setPolytope(polytope);
        }

        const size_t ntr = lhs.getDOFs(polytope);
        const size_t nte = rhs.getDOFs(polytope);

        const H1Element<KTrial, ScalarType> trialScalarFE(geometry);
        const H1Element<KTest,  ScalarType> testScalarFE(geometry);
        const size_t scalarCountTr = trialScalarFE.getCount();

        assert(scalarCountTr > 0 && ntr % scalarCountTr == 0);
        const size_t vdim = ntr / scalarCountTr;

        const auto& trTabS = trialScalarFE.getTabulation(*m_qf);
        const auto& teTab  = testScalarFE .getTabulation(*m_qf);

        m_mat.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_mat.setZero();
        ScalarType* A = m_mat.data();

        static thread_local std::vector<Math::SpatialVector<ScalarType>> GtrS;
        if (GtrS.size() < scalarCountTr) GtrS.resize(scalarCountTr);
        for (size_t a = 0; a < scalarCountTr; ++a) GtrS[a].resize(static_cast<std::uint8_t>(d));

        for (size_t qp = 0; qp < m_ps.size(); ++qp)
        {
          const auto& p = m_ps[qp];
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto Jinv = p.getJacobianInverse();

          // Physical gradients of scalar trial basis functions
          if (d == 3)
          {
            const ScalarType j00 = Jinv(0,0), j10 = Jinv(1,0), j20 = Jinv(2,0);
            const ScalarType j01 = Jinv(0,1), j11 = Jinv(1,1), j21 = Jinv(2,1);
            const ScalarType j02 = Jinv(0,2), j12 = Jinv(1,2), j22 = Jinv(2,2);

            for (size_t a = 0; a < scalarCountTr; ++a)
            {
              const auto g = trTabS.getGradient(qp, a);
              const ScalarType gx = g[0], gy = g[1], gz = g[2];
              GtrS[a][0] = j00*gx + j10*gy + j20*gz;
              GtrS[a][1] = j01*gx + j11*gy + j21*gz;
              GtrS[a][2] = j02*gx + j12*gy + j22*gz;
            }
          }
          else if (d == 2)
          {
            const ScalarType j00 = Jinv(0,0), j10 = Jinv(1,0);
            const ScalarType j01 = Jinv(0,1), j11 = Jinv(1,1);

            for (size_t a = 0; a < scalarCountTr; ++a)
            {
              const auto g = trTabS.getGradient(qp, a);
              const ScalarType gx = g[0], gy = g[1];
              GtrS[a][0] = j00*gx + j10*gy;
              GtrS[a][1] = j01*gx + j11*gy;
            }
          }
          else if (d == 1)
          {
            const ScalarType j00 = Jinv(0,0);
            for (size_t a = 0; a < scalarCountTr; ++a)
            {
              const auto g = trTabS.getGradient(qp, a);
              GtrS[a][0] = j00 * g[0];
            }
          }
          else
          {
            assert(false);
          }

          // Assemble: div(e_c * phi_s) = d(phi_s)/dx_c = GtrS[s][c]
          for (size_t b = 0; b < nte; ++b)
          {
            ScalarType* row = A + b * ntr;
            const ScalarType phi_te = teTab.getBasis(qp, b);

            for (size_t s = 0; s < scalarCountTr; ++s)
            {
              // Only components c < d contribute: physical gradient has d entries.
              // Components c >= d have zero divergence (no spatial dimension for them).
              for (size_t c = 0; c < std::min(vdim, d); ++c)
              {
                const size_t a = s * vdim + c;
                row[a] += wdet * phi_te * GtrS[s][c];
              }
            }
          }
        }

        return *this;
      }

      inline ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_mat(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      std::vector<Geometry::Point> m_ps;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_mat;
  };

  template <size_t KTrial, size_t KTest, class LHSDerived, class RHSDerived, class Scalar, class Mesh>
  QuadratureRule(
    const Dot<
      ShapeFunctionBase<
        Div<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
        H1<KTrial, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
        H1<KTest, Scalar, Mesh>, TestSpace>>&)
    -> QuadratureRule<
        Dot<
          ShapeFunctionBase<
            Div<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
            H1<KTrial, Scalar, Mesh>, TrialSpace>,
          ShapeFunctionBase<
            ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
            H1<KTest, Scalar, Mesh>, TestSpace>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Specialization for @f$\int p \, \nabla \cdot v \ dx@f$ with H1 trial and test shape functions.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int p \, \nabla \cdot v \ dx \: ,
   * @f]
   * where @f$p@f$ is an H1 trial shape function and @f$v@f$ is an H1 test shape function.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int p \, \nabla \cdot v \ dx : \texttt{QuadratureRule}}
   * {\vdash p : \texttt{H1}<K_{\mathrm{trial}}>, \ \vdash v : \texttt{H1}<K_{\mathrm{test}}>}
   * @f]
   */
  template <
    size_t KTrial, size_t KTest,
    class LHSDerived, class RHSDerived,
    class Scalar, class Mesh>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<
        ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>,
        H1<KTrial, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        Div<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
        H1<KTest, Scalar, Mesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<
          Dot<
            ShapeFunctionBase<
              ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>,
              H1<KTrial, Scalar, Mesh>, TrialSpace>,
            ShapeFunctionBase<
              Div<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
              H1<KTest, Scalar, Mesh>, TestSpace>>>::ScalarType>
  {
    public:
      using TrialFESType = H1<KTrial, Scalar, Mesh>;
      using TestFESType  = H1<KTest, Scalar, Mesh>;

      using LHSType =
        ShapeFunctionBase<ShapeFunction<LHSDerived, TrialFESType, TrialSpace>, TrialFESType, TrialSpace>;

      using RHSType =
        ShapeFunctionBase<
          Div<ShapeFunction<RHSDerived, TestFESType, TestSpace>>,
          TestFESType, TestSpace>;

      using IntegrandType = Dot<LHSType, RHSType>;
      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_ps(std::move(other.m_ps)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_mat(std::move(other.m_mat))
      {}

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();

        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes  = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes .getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe .getOrder();
        const size_t order = (k_te == 0) ? 0 : (k_tr + k_te - 1);

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::GenericPolytopeQuadrature::get(order, geometry);

          m_ps.clear();
          m_ps.reserve(m_qf->getSize());
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps.emplace_back(polytope, m_qf->getPoint(qp));
        }
        else
        {
          assert(m_qf);
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps[qp].setPolytope(polytope);
        }

        const size_t ntr = lhs.getDOFs(polytope);
        const size_t nte = rhs.getDOFs(polytope);

        const H1Element<KTrial, ScalarType> trialScalarFE(geometry);
        const H1Element<KTest,  ScalarType> testScalarFE(geometry);
        const size_t scalarCountTe = testScalarFE.getCount();

        assert(scalarCountTe > 0 && nte % scalarCountTe == 0);
        const size_t vdim = nte / scalarCountTe;

        const auto& trTab  = trialScalarFE.getTabulation(*m_qf);
        const auto& teTabS = testScalarFE .getTabulation(*m_qf);

        m_mat.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_mat.setZero();
        ScalarType* A = m_mat.data();

        static thread_local std::vector<Math::SpatialVector<ScalarType>> GteS;
        if (GteS.size() < scalarCountTe) GteS.resize(scalarCountTe);
        for (size_t b = 0; b < scalarCountTe; ++b) GteS[b].resize(static_cast<std::uint8_t>(d));

        for (size_t qp = 0; qp < m_ps.size(); ++qp)
        {
          const auto& p = m_ps[qp];
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto Jinv = p.getJacobianInverse();

          // Physical gradients of scalar test basis functions
          if (d == 3)
          {
            const ScalarType j00 = Jinv(0,0), j10 = Jinv(1,0), j20 = Jinv(2,0);
            const ScalarType j01 = Jinv(0,1), j11 = Jinv(1,1), j21 = Jinv(2,1);
            const ScalarType j02 = Jinv(0,2), j12 = Jinv(1,2), j22 = Jinv(2,2);

            for (size_t b = 0; b < scalarCountTe; ++b)
            {
              const auto g = teTabS.getGradient(qp, b);
              const ScalarType gx = g[0], gy = g[1], gz = g[2];
              GteS[b][0] = j00*gx + j10*gy + j20*gz;
              GteS[b][1] = j01*gx + j11*gy + j21*gz;
              GteS[b][2] = j02*gx + j12*gy + j22*gz;
            }
          }
          else if (d == 2)
          {
            const ScalarType j00 = Jinv(0,0), j10 = Jinv(1,0);
            const ScalarType j01 = Jinv(0,1), j11 = Jinv(1,1);

            for (size_t b = 0; b < scalarCountTe; ++b)
            {
              const auto g = teTabS.getGradient(qp, b);
              const ScalarType gx = g[0], gy = g[1];
              GteS[b][0] = j00*gx + j10*gy;
              GteS[b][1] = j01*gx + j11*gy;
            }
          }
          else if (d == 1)
          {
            const ScalarType j00 = Jinv(0,0);
            for (size_t b = 0; b < scalarCountTe; ++b)
            {
              const auto g = teTabS.getGradient(qp, b);
              GteS[b][0] = j00 * g[0];
            }
          }
          else
          {
            assert(false);
          }

          // Assemble: div(e_c * psi_s) = d(psi_s)/dx_c = GteS[s][c]
          for (size_t a = 0; a < ntr; ++a)
          {
            const ScalarType phi_tr = trTab.getBasis(qp, a);

            for (size_t s = 0; s < scalarCountTe; ++s)
            {
              // Only components c < d contribute: physical gradient has d entries.
              // Components c >= d have zero divergence (no spatial dimension for them).
              for (size_t c = 0; c < std::min(vdim, d); ++c)
              {
                const size_t b = s * vdim + c;
                A[b * ntr + a] += wdet * phi_tr * GteS[s][c];
              }
            }
          }
        }

        return *this;
      }

      inline ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_mat(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      std::vector<Geometry::Point> m_ps;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_mat;
  };

  template <size_t KTrial, size_t KTest, class LHSDerived, class RHSDerived, class Scalar, class Mesh>
  QuadratureRule(
    const Dot<
      ShapeFunctionBase<
        ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>,
        H1<KTrial, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        Div<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
        H1<KTest, Scalar, Mesh>, TestSpace>>&)
    -> QuadratureRule<
        Dot<
          ShapeFunctionBase<
            ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>,
            H1<KTrial, Scalar, Mesh>, TrialSpace>,
          ShapeFunctionBase<
            Div<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
            H1<KTest, Scalar, Mesh>, TestSpace>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Specialization for @f$\int (A J u) : J v \ dx@f$ with H1 trial and test shape functions.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int (A J u) : J v \ dx \: ,
   * @f]
   * where @f$u@f$ and @f$v@f$ are vector-valued H1 shape functions and @f$A@f$ is a coefficient function.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int (A J u) : J v \ dx : \texttt{QuadratureRule}}
   * {\vdash u : \texttt{H1}<K_{\mathrm{trial}}>, \ \vdash v : \texttt{H1}<K_{\mathrm{test}}>}
   * @f]
   */
  template <
    size_t KTrial, size_t KTest,
    class CoefficientDerived, class LHSDerived, class RHSDerived,
    class Scalar, class Mesh>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Mult<
          FunctionBase<CoefficientDerived>,
          ShapeFunctionBase<
            Jacobian<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
            H1<KTrial, Scalar, Mesh>, TrialSpace>>,
        H1<KTrial, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        Jacobian<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
        H1<KTest, Scalar, Mesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<
          Dot<
            ShapeFunctionBase<
              Mult<
                FunctionBase<CoefficientDerived>,
                ShapeFunctionBase<
                  Jacobian<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
                  H1<KTrial, Scalar, Mesh>, TrialSpace>>,
              H1<KTrial, Scalar, Mesh>, TrialSpace>,
            ShapeFunctionBase<
              Jacobian<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
              H1<KTest, Scalar, Mesh>, TestSpace>>>::ScalarType>
  {
    public:
      using TrialFESType = H1<KTrial, Scalar, Mesh>;
      using TestFESType  = H1<KTest, Scalar, Mesh>;

      using CoefficientType = FunctionBase<CoefficientDerived>;
      using CoefficientRangeType =
        typename FormLanguage::Traits<CoefficientType>::RangeType;

      using MultiplicandType =
        ShapeFunctionBase<
          Jacobian<ShapeFunction<LHSDerived, TrialFESType, TrialSpace>>,
          TrialFESType, TrialSpace>;

      using LHSType =
        ShapeFunctionBase<Mult<CoefficientType, MultiplicandType>, TrialFESType, TrialSpace>;

      using RHSType =
        ShapeFunctionBase<
          Jacobian<ShapeFunction<RHSDerived, TestFESType, TestSpace>>,
          TestFESType, TestSpace>;

      using IntegrandType = Dot<LHSType, RHSType>;
      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_ps(std::move(other.m_ps)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_cmv(std::move(other.m_cmv)),
          m_mat(std::move(other.m_mat))
      {}

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& lhs   = integrand.getLHS();
        const auto& coeff = lhs.getDerived().getLHS();
        const auto& rhs   = integrand.getRHS();

        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes  = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes .getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe .getOrder();
        const size_t order = (k_tr == 0 || k_te == 0) ? 0 : (k_tr + k_te - 2);

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::GenericPolytopeQuadrature::get(order, geometry);

          m_ps.clear();
          m_ps.reserve(m_qf->getSize());
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps.emplace_back(polytope, m_qf->getPoint(qp));
        }
        else
        {
          assert(m_qf);
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps[qp].setPolytope(polytope);
        }

        const size_t ntr = lhs.getDOFs(polytope);
        const size_t nte = rhs.getDOFs(polytope);

        const H1Element<KTrial, ScalarType> trialScalarElement(geometry);
        const H1Element<KTest,  ScalarType> testScalarElement(geometry);
        const size_t trialScalarCount = trialScalarElement.getCount();
        const size_t testScalarCount  = testScalarElement.getCount();

        assert(trialScalarCount > 0);
        assert(testScalarCount  > 0);
        assert(ntr % trialScalarCount == 0);
        assert(nte % testScalarCount  == 0);
        const size_t vdim = ntr / trialScalarCount;
        assert(vdim == nte / testScalarCount);

        const bool symmetric =
          (&trialfes.getMesh() == &testfes.getMesh()) && (ntr == nte);

        m_mat.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_mat.setZero();
        ScalarType* A_data = m_mat.data();

        const auto& trTabS = trialScalarElement.getTabulation(*m_qf);
        const auto& teTabS = testScalarElement .getTabulation(*m_qf);

        static thread_local std::vector<Math::SpatialVector<ScalarType>> GtrS;
        static thread_local std::vector<Math::SpatialVector<ScalarType>> GteS;

        if (GtrS.size() < trialScalarCount) GtrS.resize(trialScalarCount);
        if (GteS.size() < testScalarCount)  GteS.resize(testScalarCount);
        for (size_t a = 0; a < trialScalarCount; ++a) GtrS[a].resize(static_cast<std::uint8_t>(d));
        for (size_t b = 0; b < testScalarCount;  ++b) GteS[b].resize(static_cast<std::uint8_t>(d));

        for (size_t qp = 0; qp < m_ps.size(); ++qp)
        {
          const auto& p = m_ps[qp];
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto Jinv = p.getJacobianInverse();

          // Map scalar reference gradients to physical gradients
          if (d == 3)
          {
            const ScalarType a00 = Jinv(0,0), a10 = Jinv(1,0), a20 = Jinv(2,0);
            const ScalarType a01 = Jinv(0,1), a11 = Jinv(1,1), a21 = Jinv(2,1);
            const ScalarType a02 = Jinv(0,2), a12 = Jinv(1,2), a22 = Jinv(2,2);

            for (size_t a = 0; a < trialScalarCount; ++a)
            {
              const auto g = trTabS.getGradient(qp, a);
              const ScalarType gx = g[0], gy = g[1], gz = g[2];
              GtrS[a][0] = a00*gx + a10*gy + a20*gz;
              GtrS[a][1] = a01*gx + a11*gy + a21*gz;
              GtrS[a][2] = a02*gx + a12*gy + a22*gz;
            }
            for (size_t b = 0; b < testScalarCount; ++b)
            {
              const auto g = teTabS.getGradient(qp, b);
              const ScalarType gx = g[0], gy = g[1], gz = g[2];
              GteS[b][0] = a00*gx + a10*gy + a20*gz;
              GteS[b][1] = a01*gx + a11*gy + a21*gz;
              GteS[b][2] = a02*gx + a12*gy + a22*gz;
            }
          }
          else if (d == 2)
          {
            const ScalarType a00 = Jinv(0,0), a10 = Jinv(1,0);
            const ScalarType a01 = Jinv(0,1), a11 = Jinv(1,1);

            for (size_t a = 0; a < trialScalarCount; ++a)
            {
              const auto g = trTabS.getGradient(qp, a);
              const ScalarType gx = g[0], gy = g[1];
              GtrS[a][0] = a00*gx + a10*gy;
              GtrS[a][1] = a01*gx + a11*gy;
            }
            for (size_t b = 0; b < testScalarCount; ++b)
            {
              const auto g = teTabS.getGradient(qp, b);
              const ScalarType gx = g[0], gy = g[1];
              GteS[b][0] = a00*gx + a10*gy;
              GteS[b][1] = a01*gx + a11*gy;
            }
          }
          else if (d == 1)
          {
            const ScalarType a00 = Jinv(0,0);
            for (size_t a = 0; a < trialScalarCount; ++a)
            {
              const auto g = trTabS.getGradient(qp, a);
              GtrS[a][0] = a00 * g[0];
            }
            for (size_t b = 0; b < testScalarCount; ++b)
            {
              const auto g = teTabS.getGradient(qp, b);
              GteS[b][0] = a00 * g[0];
            }
          }
          else
          {
            assert(false);
          }

          if constexpr (std::is_same_v<CoefficientRangeType, ScalarType>)
          {
            const ScalarType csv = coeff.getValue(p);

            // Block-diagonal assembly with scalar coefficient
            if (symmetric)
            {
              for (size_t ib = 0; ib < testScalarCount; ++ib)
              {
                const auto& gb = GteS[ib];
                for (size_t ia = 0; ia <= ib; ++ia)
                {
                  const ScalarType kij = wdet * csv * Math::dot(gb, GtrS[ia]);
                  if (kij == ScalarType(0))
                    continue;
                  for (size_t c = 0; c < vdim; ++c)
                  {
                    const size_t row = ib * vdim + c;
                    const size_t col = ia * vdim + c;
                    A_data[row * ntr + col] += kij;
                  }
                }
              }
            }
            else
            {
              for (size_t ib = 0; ib < testScalarCount; ++ib)
              {
                const auto& gb = GteS[ib];
                for (size_t ia = 0; ia < trialScalarCount; ++ia)
                {
                  const ScalarType kij = wdet * csv * Math::dot(gb, GtrS[ia]);
                  if (kij == ScalarType(0))
                    continue;
                  for (size_t c = 0; c < vdim; ++c)
                  {
                    const size_t row = ib * vdim + c;
                    const size_t col = ia * vdim + c;
                    A_data[row * ntr + col] += kij;
                  }
                }
              }
            }
          }
          else if constexpr (std::is_same_v<CoefficientRangeType, Math::Matrix<ScalarType>>)
          {
            coeff.getValue(m_cmv, p);

            // With matrix coefficient, the Frobenius inner product becomes:
            // (A * J_trial) : J_test = sum_{c,j} A(c,c2) * grad_trial_c2[j] * grad_test_c[j]
            // This breaks block-diagonal structure.
            for (size_t ib = 0; ib < testScalarCount; ++ib)
            {
              const auto& gb = GteS[ib];
              for (size_t ia = 0; ia < trialScalarCount; ++ia)
              {
                // grad_tr . grad_te (scalar dot of physical gradients)
                const ScalarType grad_dot = Math::dot(gb, GtrS[ia]);

                for (size_t d_comp = 0; d_comp < vdim; ++d_comp)
                {
                  const size_t row = ib * vdim + d_comp;
                  for (size_t c_comp = 0; c_comp < vdim; ++c_comp)
                  {
                    const size_t col = ia * vdim + c_comp;
                    A_data[row * ntr + col] += wdet * m_cmv(d_comp, c_comp) * grad_dot;
                  }
                }
              }
            }
          }
          else
          {
            assert(false);
          }
        }

        if (symmetric)
        {
          m_mat.template triangularView<Eigen::Upper>() =
            m_mat.transpose().template triangularView<Eigen::Upper>();
        }

        return *this;
      }

      inline ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_mat(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      std::vector<Geometry::Point> m_ps;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      Math::Matrix<ScalarType> m_cmv;
      Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_mat;
  };

  template <
    size_t KTrial, size_t KTest,
    class CoefficientDerived, class LHSDerived, class RHSDerived,
    class Scalar, class Mesh>
  QuadratureRule(
    const Dot<
      ShapeFunctionBase<
        Mult<
          FunctionBase<CoefficientDerived>,
          ShapeFunctionBase<
            Jacobian<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
            H1<KTrial, Scalar, Mesh>, TrialSpace>>,
        H1<KTrial, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        Jacobian<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
        H1<KTest, Scalar, Mesh>, TestSpace>>&)
    -> QuadratureRule<
        Dot<
          ShapeFunctionBase<
            Mult<
              FunctionBase<CoefficientDerived>,
              ShapeFunctionBase<
                Jacobian<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
                H1<KTrial, Scalar, Mesh>, TrialSpace>>,
            H1<KTrial, Scalar, Mesh>, TrialSpace>,
          ShapeFunctionBase<
            Jacobian<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
            H1<KTest, Scalar, Mesh>, TestSpace>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Specialization for @f$\int \nabla u \cdot \nabla v \ dx@f$ with H1 trial and test shape functions.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int \nabla u \cdot \nabla v \ dx \: ,
   * @f]
   * where @f$u@f$ and @f$v@f$ are H1 shape functions.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int \nabla u \cdot \nabla v \ dx : \texttt{QuadratureRule}}
   * {\vdash u : \texttt{H1}<K_{\mathrm{trial}}>, \ \vdash v : \texttt{H1}<K_{\mathrm{test}}>}
   * @f]
   */
  template <
    size_t KTrial, size_t KTest,
    class LHSDerived, class RHSDerived,
    class Scalar, class Mesh>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Grad<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
        H1<KTrial, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        Grad<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
        H1<KTest, Scalar, Mesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<
          Dot<
            ShapeFunctionBase<
              Grad<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
              H1<KTrial, Scalar, Mesh>, TrialSpace>,
            ShapeFunctionBase<
              Grad<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
              H1<KTest, Scalar, Mesh>, TestSpace>>>::ScalarType>
  {
    public:
      using TrialFESType = H1<KTrial, Scalar, Mesh>;
      using TestFESType  = H1<KTest, Scalar, Mesh>;

      using LHSType =
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, TrialFESType, TrialSpace>>, TrialFESType, TrialSpace>;

      using RHSType =
        ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, TestFESType, TestSpace>>,  TestFESType, TestSpace>;

      using IntegrandType = Dot<LHSType, RHSType>;

      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point),
          m_dim(0),
          m_trCount(0),
          m_teCount(0)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point),
          m_dim(0),
          m_trCount(0),
          m_teCount(0)
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_ps(std::move(other.m_ps)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_dim(std::exchange(other.m_dim, 0)),
          m_trCount(std::exchange(other.m_trCount, 0)),
          m_teCount(std::exchange(other.m_teCount, 0)),
          m_trRefGrad(std::move(other.m_trRefGrad)),
          m_teRefGrad(std::move(other.m_teRefGrad)),
          m_mat(std::move(other.m_mat))
      {}

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();

        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes  = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes .getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe .getOrder();
        const size_t order = (k_tr == 0 || k_te == 0) ? 0 : (k_tr + k_te - 2);

        const auto geometry = polytope.getGeometry();
        const bool recompute_qf = (!m_set || m_order != order || m_geometry != geometry);

        if (recompute_qf)
        {
          m_set = true;
          m_order = order;
          m_geometry = geometry;

          m_qf = &QF::GenericPolytopeQuadrature::get(order, geometry);

          m_ps.clear();
          m_ps.reserve(m_qf->getSize());
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps.emplace_back(polytope, m_qf->getPoint(qp));
        }
        else
        {
          assert(m_qf);
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps[qp].setPolytope(polytope);
        }

        const size_t ntr = lhs.getDOFs(polytope);
        const size_t nte = rhs.getDOFs(polytope);

        // Symmetry applicability for this specialization:
        const bool symmetric =
          (&trialfes.getMesh() == &testfes.getMesh()) && (ntr == nte);
        // If you can assert stronger (same FE object / same ordering), do it here.

        // Row-major local matrix for fast row writes.
        m_mat.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_mat.setZero();
        ScalarType* A = m_mat.data();

        const auto& trTab = trialfe.getTabulation(*m_qf);
        const auto& teTab = testfe .getTabulation(*m_qf);

        // Scratch for physical gradients at this qp.
        static thread_local std::vector<Math::SpatialVector<ScalarType>> Gtr;
        static thread_local std::vector<Math::SpatialVector<ScalarType>> Gte;

        if (Gtr.size() < ntr) Gtr.resize(ntr);
        if (Gte.size() < nte) Gte.resize(nte);
        for (size_t a = 0; a < ntr; ++a) Gtr[a].resize(static_cast<std::uint8_t>(d));
        for (size_t b = 0; b < nte; ++b) Gte[b].resize(static_cast<std::uint8_t>(d));

        for (size_t qp = 0; qp < m_ps.size(); ++qp)
        {
          const auto& p = m_ps[qp];
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto Jinv = p.getJacobianInverse();

          // Build physical gradients for trial and test
          if (d == 3)
          {
            const ScalarType a00 = Jinv(0,0), a10 = Jinv(1,0), a20 = Jinv(2,0);
            const ScalarType a01 = Jinv(0,1), a11 = Jinv(1,1), a21 = Jinv(2,1);
            const ScalarType a02 = Jinv(0,2), a12 = Jinv(1,2), a22 = Jinv(2,2);

            for (size_t a = 0; a < ntr; ++a)
            {
              const auto g = trTab.getGradient(qp, a);
              const ScalarType gx = g[0], gy = g[1], gz = g[2];
              Gtr[a][0] = a00*gx + a10*gy + a20*gz;
              Gtr[a][1] = a01*gx + a11*gy + a21*gz;
              Gtr[a][2] = a02*gx + a12*gy + a22*gz;
            }
            for (size_t b = 0; b < nte; ++b)
            {
              const auto g = teTab.getGradient(qp, b);
              const ScalarType gx = g[0], gy = g[1], gz = g[2];
              Gte[b][0] = a00*gx + a10*gy + a20*gz;
              Gte[b][1] = a01*gx + a11*gy + a21*gz;
              Gte[b][2] = a02*gx + a12*gy + a22*gz;
            }
          }
          else if (d == 2)
          {
            const ScalarType a00 = Jinv(0,0), a10 = Jinv(1,0);
            const ScalarType a01 = Jinv(0,1), a11 = Jinv(1,1);

            for (size_t a = 0; a < ntr; ++a)
            {
              const auto g = trTab.getGradient(qp, a);
              const ScalarType gx = g[0], gy = g[1];
              Gtr[a][0] = a00*gx + a10*gy;
              Gtr[a][1] = a01*gx + a11*gy;
            }
            for (size_t b = 0; b < nte; ++b)
            {
              const auto g = teTab.getGradient(qp, b);
              const ScalarType gx = g[0], gy = g[1];
              Gte[b][0] = a00*gx + a10*gy;
              Gte[b][1] = a01*gx + a11*gy;
            }
          }
          else if (d == 1)
          {
            const ScalarType a00 = Jinv(0,0);
            for (size_t a = 0; a < ntr; ++a)
            {
              const auto g = trTab.getGradient(qp, a);
              Gtr[a][0] = a00 * g[0];
            }
            for (size_t b = 0; b < nte; ++b)
            {
              const auto g = teTab.getGradient(qp, b);
              Gte[b][0] = a00 * g[0];
            }
          }
          else
          {
            assert(false);
          }

          if (symmetric)
          {
            // Accumulate only lower triangle (i,j with j <= i)
            // Use trial index as column, test index as row: A[row*ntr + col]
            for (size_t i = 0; i < nte; ++i)
            {
              ScalarType* row = A + i * ntr;
              const auto& gi = Gte[i];

              // diagonal
              row[i] += wdet * Math::dot(gi, Gtr[i]);

              for (size_t j = 0; j < i; ++j)
                row[j] += wdet * Math::dot(gi, Gtr[j]);
            }
          }
          else
          {
            // General (non-symmetric) accumulation
            for (size_t b = 0; b < nte; ++b)
            {
              ScalarType* row = A + b * ntr;
              const auto& gb = Gte[b];
              for (size_t a = 0; a < ntr; ++a)
                row[a] += wdet * Math::dot(gb, Gtr[a]);
            }
          }
        }

        if (symmetric)
        {
          // Mirror lower -> upper
          // m_mat is row-major, but Eigen triangular views still work.
          // We filled (i,i) and (i,j) for j<i. Mirror to upper.
          m_mat.template triangularView<Eigen::Upper>() =
            m_mat.transpose().template triangularView<Eigen::Upper>();
        }

        return *this;
      }

      inline ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_mat(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      std::vector<Geometry::Point> m_ps;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      // sizes
      size_t m_dim;
      size_t m_qps;
      size_t m_trCount;
      size_t m_teCount;

      // reference gradients (qp-major)
      std::vector<Math::SpatialVector<ScalarType>> m_trRefGrad;
      std::vector<Math::SpatialVector<ScalarType>> m_teRefGrad;

      // local matrix (rows=test, cols=trial)
      Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_mat;
    };

  // CTAD helper
  template <size_t KTrial, size_t KTest, class LHSDerived, class RHSDerived, class Scalar, class Mesh>
  QuadratureRule(
    const Dot<
      ShapeFunctionBase<
        Grad<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
        H1<KTrial, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        Grad<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
        H1<KTest, Scalar, Mesh>, TestSpace>>&)
  -> QuadratureRule<
       Dot<
         ShapeFunctionBase<
           Grad<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
           H1<KTrial, Scalar, Mesh>, TrialSpace>,
         ShapeFunctionBase<
           Grad<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
           H1<KTest, Scalar, Mesh>, TestSpace>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Specialization for @f$\int J u : J v \ dx@f$ with H1 trial and test shape functions.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int J u : J v \ dx \: ,
   * @f]
   * where @f$u@f$ and @f$v@f$ are vector-valued H1 shape functions.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int J u : J v \ dx : \texttt{QuadratureRule}}
   * {\vdash u : \texttt{H1}<K_{\mathrm{trial}}>, \ \vdash v : \texttt{H1}<K_{\mathrm{test}}>}
   * @f]
   */
  template <
    size_t KTrial, size_t KTest,
    class LHSDerived, class RHSDerived,
    class Scalar, class Mesh>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Jacobian<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
        H1<KTrial, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        Jacobian<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
        H1<KTest, Scalar, Mesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<
          Dot<
            ShapeFunctionBase<
              Jacobian<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
              H1<KTrial, Scalar, Mesh>, TrialSpace>,
            ShapeFunctionBase<
              Jacobian<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
              H1<KTest, Scalar, Mesh>, TestSpace>>>::ScalarType>
  {
    public:
      using TrialFESType = H1<KTrial, Scalar, Mesh>;
      using TestFESType  = H1<KTest, Scalar, Mesh>;

      using LHSType =
        ShapeFunctionBase<
          Jacobian<ShapeFunction<LHSDerived, TrialFESType, TrialSpace>>,
          TrialFESType, TrialSpace>;

      using RHSType =
        ShapeFunctionBase<
          Jacobian<ShapeFunction<RHSDerived, TestFESType, TestSpace>>,
          TestFESType, TestSpace>;

      using IntegrandType = Dot<LHSType, RHSType>;
      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_ps(std::move(other.m_ps)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_mat(std::move(other.m_mat))
      {}

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();

        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes  = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes .getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe .getOrder();
        const size_t order = (k_tr == 0 || k_te == 0) ? 0 : (k_tr + k_te - 2);

        const auto geometry = polytope.getGeometry();
        const bool recompute_qf = (!m_set || m_order != order || m_geometry != geometry);

        if (recompute_qf)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::GenericPolytopeQuadrature::get(order, geometry);

          m_ps.clear();
          m_ps.reserve(m_qf->getSize());
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps.emplace_back(polytope, m_qf->getPoint(qp));
        }
        else
        {
          assert(m_qf);
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps[qp].setPolytope(polytope);
        }

        const size_t ntr = lhs.getDOFs(polytope);
        const size_t nte = rhs.getDOFs(polytope);

        // --- infer vdim from H1Element<K, Scalar> conventions ---
        // Each element's scalar count depends on its polynomial order.
        const H1Element<KTrial, ScalarType> trialScalarElement(polytope.getGeometry());
        const H1Element<KTest,  ScalarType> testScalarElement(polytope.getGeometry());
        const size_t trialScalarCount = trialScalarElement.getCount();
        const size_t testScalarCount  = testScalarElement.getCount();

        assert(trialScalarCount > 0);
        assert(testScalarCount  > 0);
        assert(ntr % trialScalarCount == 0);
        assert(nte % testScalarCount  == 0);

        const size_t vdim_tr = ntr / trialScalarCount;
        assert(vdim_tr == nte / testScalarCount);
        const size_t vdim = vdim_tr;

        const bool symmetric =
          (&trialfes.getMesh() == &testfes.getMesh()) && (ntr == nte);

        m_mat.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_mat.setZero();
        ScalarType* A = m_mat.data(); // row-major (rows=test, cols=trial)

        // Use scalar tabulations (fast, cached in H1Element<K, Scalar>::getTabulation)
        const auto& trTabS = trialScalarElement.getTabulation(*m_qf);
        const auto& teTabS = testScalarElement.getTabulation(*m_qf);

        static thread_local std::vector<Math::SpatialVector<ScalarType>> GtrS;
        static thread_local std::vector<Math::SpatialVector<ScalarType>> GteS;

        if (GtrS.size() < trialScalarCount) GtrS.resize(trialScalarCount);
        if (GteS.size() < testScalarCount)  GteS.resize(testScalarCount);
        for (size_t a = 0; a < trialScalarCount; ++a) GtrS[a].resize(static_cast<std::uint8_t>(d));
        for (size_t b = 0; b < testScalarCount;  ++b) GteS[b].resize(static_cast<std::uint8_t>(d));

        for (size_t qp = 0; qp < m_ps.size(); ++qp)
        {
          const auto& p = m_ps[qp];
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto Jinv = p.getJacobianInverse();

          // Map scalar reference gradients to physical gradients (once per scalar DOF)
          if (d == 3)
          {
            const ScalarType a00 = Jinv(0,0), a10 = Jinv(1,0), a20 = Jinv(2,0);
            const ScalarType a01 = Jinv(0,1), a11 = Jinv(1,1), a21 = Jinv(2,1);
            const ScalarType a02 = Jinv(0,2), a12 = Jinv(1,2), a22 = Jinv(2,2);

            for (size_t a = 0; a < trialScalarCount; ++a)
            {
              const auto g = trTabS.getGradient(qp, a);
              const ScalarType gx = g[0], gy = g[1], gz = g[2];
              GtrS[a][0] = a00*gx + a10*gy + a20*gz;
              GtrS[a][1] = a01*gx + a11*gy + a21*gz;
              GtrS[a][2] = a02*gx + a12*gy + a22*gz;
            }
            for (size_t b = 0; b < testScalarCount; ++b)
            {
              const auto g = teTabS.getGradient(qp, b);
              const ScalarType gx = g[0], gy = g[1], gz = g[2];
              GteS[b][0] = a00*gx + a10*gy + a20*gz;
              GteS[b][1] = a01*gx + a11*gy + a21*gz;
              GteS[b][2] = a02*gx + a12*gy + a22*gz;
            }
          }
          else if (d == 2)
          {
            const ScalarType a00 = Jinv(0,0), a10 = Jinv(1,0);
            const ScalarType a01 = Jinv(0,1), a11 = Jinv(1,1);

            for (size_t a = 0; a < trialScalarCount; ++a)
            {
              const auto g = trTabS.getGradient(qp, a);
              const ScalarType gx = g[0], gy = g[1];
              GtrS[a][0] = a00*gx + a10*gy;
              GtrS[a][1] = a01*gx + a11*gy;
            }
            for (size_t b = 0; b < testScalarCount; ++b)
            {
              const auto g = teTabS.getGradient(qp, b);
              const ScalarType gx = g[0], gy = g[1];
              GteS[b][0] = a00*gx + a10*gy;
              GteS[b][1] = a01*gx + a11*gy;
            }
          }
          else if (d == 1)
          {
            const ScalarType a00 = Jinv(0,0);
            for (size_t a = 0; a < trialScalarCount; ++a)
            {
              const auto g = trTabS.getGradient(qp, a);
              GtrS[a][0] = a00 * g[0];
            }
            for (size_t b = 0; b < testScalarCount; ++b)
            {
              const auto g = teTabS.getGradient(qp, b);
              GteS[b][0] = a00 * g[0];
            }
          }
          else
          {
            assert(false);
          }

          // Assemble: for each component c, add the same scalar matrix into block (c,c)
          // Local DOF mapping is exactly your convention: dof = scalarIndex*vdim + c.
          if (symmetric)
          {
            // lower triangle in vector-dof space corresponds to scalar lower + component blocks
            for (size_t ib = 0; ib < testScalarCount; ++ib)
            {
              const auto& gb = GteS[ib];

              for (size_t ia = 0; ia <= ib; ++ia)
              {
                const ScalarType kij = wdet * Math::dot(gb, GtrS[ia]);
                if (kij == ScalarType(0))
                  continue;

                for (size_t c = 0; c < vdim; ++c)
                {
                  const size_t row = ib * vdim + c;
                  const size_t col = ia * vdim + c;
                  A[row * ntr + col] += kij;
                }
              }
            }
          }
          else
          {
            for (size_t ib = 0; ib < testScalarCount; ++ib)
            {
              const auto& gb = GteS[ib];

              for (size_t ia = 0; ia < trialScalarCount; ++ia)
              {
                const ScalarType kij = wdet * Math::dot(gb, GtrS[ia]);
                if (kij == ScalarType(0))
                  continue;

                for (size_t c = 0; c < vdim; ++c)
                {
                  const size_t row = ib * vdim + c;
                  const size_t col = ia * vdim + c;
                  A[row * ntr + col] += kij;
                }
              }
            }
          }
        }

        if (symmetric)
        {
          m_mat.template triangularView<Eigen::Upper>() =
            m_mat.transpose().template triangularView<Eigen::Upper>();
        }

        return *this;
      }

      inline ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_mat(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      std::vector<Geometry::Point> m_ps;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_mat;
  };

  // CTAD helper
  template <size_t KTrial, size_t KTest, class LHSDerived, class RHSDerived, class Scalar, class Mesh>
  QuadratureRule(
    const Dot<
      ShapeFunctionBase<
        Jacobian<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
        H1<KTrial, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        Jacobian<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
        H1<KTest, Scalar, Mesh>, TestSpace>>&)
  -> QuadratureRule<
       Dot<
         ShapeFunctionBase<
           Jacobian<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
           H1<KTrial, Scalar, Mesh>, TrialSpace>,
         ShapeFunctionBase<
           Jacobian<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
           H1<KTest, Scalar, Mesh>, TestSpace>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Specialization for @f$\int (\mathbf{J}\,u \cdot f) \cdot v \ dx@f$ with H1 trial and test shape functions.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int (\mathbf{J} \: u \cdot f) \cdot v \ dx \: ,
   * @f]
   * where @f$u@f$ and @f$v@f$ are H1 shape functions and @f$f@f$ is a
   * coefficient function (typically a vector-valued GridFunction).
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int (\mathbf{J} \: u \cdot f) \cdot v \ dx :
   * \texttt{QuadratureRule}}
   * {\vdash u : \texttt{H1}<K_{\mathrm{trial}}>, \ \vdash v : \texttt{H1}<K_{\mathrm{test}}>}
   * @f]
   *
   * Implements the linearized convection term for Navier-Stokes using
   * multi-point quadrature appropriate for the polynomial degree.
   */
  template <
    size_t KTrial, size_t KTest,
    class CoefficientDerived, class LHSDerived, class RHSDerived,
    class Scalar, class Mesh>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Mult<
          ShapeFunctionBase<
            Jacobian<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
            H1<KTrial, Scalar, Mesh>, TrialSpace>,
          FunctionBase<CoefficientDerived>>,
        H1<KTrial, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
        H1<KTest, Scalar, Mesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<
          Dot<
            ShapeFunctionBase<
              Mult<
                ShapeFunctionBase<
                  Jacobian<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
                  H1<KTrial, Scalar, Mesh>, TrialSpace>,
                FunctionBase<CoefficientDerived>>,
              H1<KTrial, Scalar, Mesh>, TrialSpace>,
            ShapeFunctionBase<
              ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
              H1<KTest, Scalar, Mesh>, TestSpace>>>::ScalarType>
  {
    public:
      using TrialFESType = H1<KTrial, Scalar, Mesh>;
      using TestFESType  = H1<KTest, Scalar, Mesh>;

      using TrialSFType =
        ShapeFunctionBase<
          Jacobian<ShapeFunction<LHSDerived, TrialFESType, TrialSpace>>,
          TrialFESType, TrialSpace>;

      using CoefficientType = FunctionBase<CoefficientDerived>;

      using LHSType =
        ShapeFunctionBase<
          Mult<TrialSFType, CoefficientType>,
          TrialFESType, TrialSpace>;

      using RHSType =
        ShapeFunctionBase<
          ShapeFunction<RHSDerived, TestFESType, TestSpace>,
          TestFESType, TestSpace>;

      using IntegrandType = Dot<LHSType, RHSType>;
      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_ps(std::move(other.m_ps)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_mat(std::move(other.m_mat))
      {}

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();

        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes  = rhs.getFiniteElementSpace();

        // The coefficient is the RHS of the Mult node
        const auto& coeff = lhs.getDerived().getRHS();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes .getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe .getOrder();
        // Gradient of trial (k_tr - 1) times value of test (k_te)
        const size_t order = (k_tr == 0) ? k_te : (k_tr + k_te - 1);

        const auto geometry = polytope.getGeometry();
        const bool recompute_qf = (!m_set || m_order != order || m_geometry != geometry);

        if (recompute_qf)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::GenericPolytopeQuadrature::get(order, geometry);

          m_ps.clear();
          m_ps.reserve(m_qf->getSize());
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps.emplace_back(polytope, m_qf->getPoint(qp));
        }
        else
        {
          assert(m_qf);
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps[qp].setPolytope(polytope);
        }

        const size_t ntr = lhs.getDOFs(polytope);
        const size_t nte = rhs.getDOFs(polytope);

        const H1Element<KTrial, ScalarType> trialScalarFE(geometry);
        const H1Element<KTest,  ScalarType> testScalarFE(geometry);
        const size_t trialScalarCount = trialScalarFE.getCount();
        const size_t testScalarCount  = testScalarFE.getCount();

        assert(trialScalarCount > 0 && ntr % trialScalarCount == 0);
        assert(testScalarCount  > 0 && nte % testScalarCount  == 0);
        const size_t vdim = ntr / trialScalarCount;
        assert(vdim == nte / testScalarCount);

        const auto& trTab = trialScalarFE.getTabulation(*m_qf);
        const auto& teTab = testScalarFE .getTabulation(*m_qf);

        m_mat.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_mat.setZero();
        ScalarType* A = m_mat.data();

        static thread_local std::vector<Math::SpatialVector<ScalarType>> GtrS;
        if (GtrS.size() < trialScalarCount) GtrS.resize(trialScalarCount);
        for (size_t a = 0; a < trialScalarCount; ++a)
          GtrS[a].resize(static_cast<std::uint8_t>(d));

        for (size_t qp = 0; qp < m_ps.size(); ++qp)
        {
          const auto& p = m_ps[qp];
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto Jinv = p.getJacobianInverse();

          // Map scalar reference gradients to physical gradients
          if (d == 3)
          {
            const ScalarType a00 = Jinv(0,0), a10 = Jinv(1,0), a20 = Jinv(2,0);
            const ScalarType a01 = Jinv(0,1), a11 = Jinv(1,1), a21 = Jinv(2,1);
            const ScalarType a02 = Jinv(0,2), a12 = Jinv(1,2), a22 = Jinv(2,2);

            for (size_t a = 0; a < trialScalarCount; ++a)
            {
              const auto g = trTab.getGradient(qp, a);
              const ScalarType gx = g[0], gy = g[1], gz = g[2];
              GtrS[a][0] = a00*gx + a10*gy + a20*gz;
              GtrS[a][1] = a01*gx + a11*gy + a21*gz;
              GtrS[a][2] = a02*gx + a12*gy + a22*gz;
            }
          }
          else if (d == 2)
          {
            const ScalarType a00 = Jinv(0,0), a10 = Jinv(1,0);
            const ScalarType a01 = Jinv(0,1), a11 = Jinv(1,1);

            for (size_t a = 0; a < trialScalarCount; ++a)
            {
              const auto g = trTab.getGradient(qp, a);
              const ScalarType gx = g[0], gy = g[1];
              GtrS[a][0] = a00*gx + a10*gy;
              GtrS[a][1] = a01*gx + a11*gy;
            }
          }
          else if (d == 1)
          {
            const ScalarType a00 = Jinv(0,0);
            for (size_t a = 0; a < trialScalarCount; ++a)
            {
              const auto g = trTab.getGradient(qp, a);
              GtrS[a][0] = a00 * g[0];
            }
          }
          else
          {
            assert(false);
          }

          // Evaluate coefficient at this quadrature point
          const auto fval = coeff.getValue(p);

          // Assemble: K(b*vdim+c, a*vdim+c) += wdet * (∇φ_a · f) * ψ_b
          for (size_t ib = 0; ib < testScalarCount; ++ib)
          {
            const ScalarType phi_te = teTab.getBasis(qp, ib);
            for (size_t ia = 0; ia < trialScalarCount; ++ia)
            {
              const ScalarType gradDotF = Math::dot(GtrS[ia], fval);
              const ScalarType kij = wdet * gradDotF * phi_te;
              if (kij == ScalarType(0))
                continue;
              for (size_t c = 0; c < vdim; ++c)
                A[(ib * vdim + c) * ntr + (ia * vdim + c)] += kij;
            }
          }
        }

        return *this;
      }

      inline ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_mat(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      std::vector<Geometry::Point> m_ps;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_mat;
  };

  // CTAD helper
  template <
    size_t KTrial, size_t KTest,
    class CoefficientDerived, class LHSDerived, class RHSDerived,
    class Scalar, class Mesh>
  QuadratureRule(
    const Dot<
      ShapeFunctionBase<
        Mult<
          ShapeFunctionBase<
            Jacobian<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
            H1<KTrial, Scalar, Mesh>, TrialSpace>,
          FunctionBase<CoefficientDerived>>,
        H1<KTrial, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
        H1<KTest, Scalar, Mesh>, TestSpace>>&)
    -> QuadratureRule<
         Dot<
           ShapeFunctionBase<
             Mult<
               ShapeFunctionBase<
                 Jacobian<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
                 H1<KTrial, Scalar, Mesh>, TrialSpace>,
               FunctionBase<CoefficientDerived>>,
             H1<KTrial, Scalar, Mesh>, TrialSpace>,
           ShapeFunctionBase<
             ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
             H1<KTest, Scalar, Mesh>, TestSpace>>>;
}

#endif
