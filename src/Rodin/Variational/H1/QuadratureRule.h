/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_QUADRATURERULE_H
#define RODIN_VARIATIONAL_H1_QUADRATURERULE_H

/**
 * @file
 * @brief H1 quadrature rules with point-local reference-table borrowing.
 *
 * ## Coefficient evaluation and reference data
 *
 * A coefficient callback may query unrelated reference formulas and evict a
 * bounded reference-table slot. Coefficient-bearing rules evaluate owned
 * coefficient values before borrowing the tables consumed by that stage.
 * No table contents are read across a later arbitrary coefficient callback.
 * This protocol preserves the bounded store and requires no new cache,
 * coordinate comparisons or distributed communication.
 */

#include "Rodin/QF/PolytopeQuadratureFormula.h"

#include "Rodin/FormLanguage/Traits.h"
#include "Rodin/Math/Common.h"
#include "Rodin/Variational/IntegrationPoint.h"
#include "Rodin/Variational/ShapeFunction.h"

#include "H1.h"
#include "H1Element.h"
#include "Rodin/Math/Traits.h"

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
    ShapeFunctionBase<ShapeFunction<NestedDerived, H1<K, Scalar, Mesh>, TestSpace>,
      H1<K, Scalar, Mesh>, TestSpace>>
    : public LinearFormIntegratorBase<typename FormLanguage::Traits<
        ShapeFunctionBase<ShapeFunction<NestedDerived, H1<K, Scalar, Mesh>, TestSpace>,
          H1<K, Scalar, Mesh>, TestSpace>>::ScalarType>
  {
    public:
      /// @brief Reports this handler as an optimized specialization.
      static constexpr bool Specialized = true;

      /// @brief Finite element space type.
      using FESType = H1<K, Scalar, Mesh>;
      /// @brief Integrand expression type.
      using IntegrandType =
        ShapeFunctionBase<ShapeFunction<NestedDerived, FESType, TestSpace>, FESType,
          TestSpace>;
      /// @brief Range type of the integrand.
      using IntegrandRangeType = typename FormLanguage::Traits<IntegrandType>::RangeType;
      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;
      /// @brief Parent class type.
      using Parent = LinearFormIntegratorBase<ScalarType>;

      /**
       * @brief Constructs the integrator for the given integrand.
       * @param integrand Expression to integrate.
       */
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr),
          m_quadrature(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr),
          m_quadrature(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_quadrature(std::move(other.m_quadrature)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_vec(std::move(other.m_vec))
      {}

      /**
       * @brief Gets the integrand.
       * @returns The integrand.
       */
      constexpr const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      /**
       * @brief Returns the polytope the integrator is bound to.
       * @returns The polytope the integrator is bound to.
       */
      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      /**
       * @brief Binds the integrator to a polytope and tabulates the quadrature on it.
       * @param polytope Mesh entity used by this operation.
       * @returns Reference to this object after the operation.
       */
      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& fes = integrand.getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, idx);

        const size_t order = this->getOrder(polytope).value_or(
          integrand.getOrder(polytope).value_or(fe.getOrder()));

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set = true;
          m_order = order;
          m_geometry = geometry;

          m_qf = &QF::PolytopeQuadratureFormula::get(order, geometry);
        }

        assert(m_qf);
        m_quadrature = std::make_unique<Geometry::PolytopeQuadrature>(polytope, *m_qf);

        const size_t nte = integrand.getDOFs(polytope);

        const H1Element<K, ScalarType> scalarFE(geometry);
        const auto& tab = scalarFE.getTabulation(*m_qf);

        m_vec.resize(static_cast<Eigen::Index>(nte));
        m_vec.setZero();

        assert(m_quadrature);
        const auto& q = *m_quadrature;
        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          if constexpr (std::is_same_v<IntegrandRangeType, ScalarType>)
          {
            for (size_t local = 0; local < nte; ++local)
              m_vec(local) += wdet * tab.getBasis(qp, local);
          }
          else if constexpr (FormLanguage::IsVectorRange<IntegrandRangeType>::Value)
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
            static_assert(std::is_same_v<IntegrandRangeType, ScalarType> ||
                FormLanguage::IsVectorRange<IntegrandRangeType>::Value,
              "Unsupported H1 Integral(v) range type. Expected scalar or vector-valued "
              "shape function.");
          }
        }

        return *this;
      }

      /**
       * @brief Returns an entry of the element vector.
       * @param local Index in the local numbering.
       * @returns Integral computed by the quadrature rule.
       */
      inline ScalarType integrate(size_t local) final override
      {
        return m_vec(local);
      }

      /**
       * @brief Returns the integration region.
       * @returns The integration region.
       */
      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      std::unique_ptr<Geometry::PolytopeQuadrature> m_quadrature;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      Math::Vector<ScalarType> m_vec;
  };

  /**
   * @brief Deduction guide for @c QuadratureRule.
   * @param integrand Expression to integrate.
   */
  template <size_t K, class NestedDerived, class Scalar, class Mesh>
  QuadratureRule(
    const ShapeFunctionBase<ShapeFunction<NestedDerived, H1<K, Scalar, Mesh>, TestSpace>,
      H1<K, Scalar, Mesh>, TestSpace>& integrand)
    -> QuadratureRule<
      ShapeFunctionBase<ShapeFunction<NestedDerived, H1<K, Scalar, Mesh>, TestSpace>,
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
  class QuadratureRule<ShapeFunctionBase<
    Dot<FunctionBase<LHSDerived>,
      ShapeFunctionBase<ShapeFunction<RHSDerived, H1<K, Scalar, Mesh>, TestSpace>,
        H1<K, Scalar, Mesh>, TestSpace>>,
    H1<K, Scalar, Mesh>, TestSpace>>
    : public LinearFormIntegratorBase<typename FormLanguage::Traits<ShapeFunctionBase<
        Dot<FunctionBase<LHSDerived>,
          ShapeFunctionBase<ShapeFunction<RHSDerived, H1<K, Scalar, Mesh>, TestSpace>,
            H1<K, Scalar, Mesh>, TestSpace>>,
        H1<K, Scalar, Mesh>, TestSpace>>::ScalarType>
  {
    public:
      /// @brief Reports this handler as an optimized specialization.
      static constexpr bool Specialized = true;

      /// @brief Finite element space type.
      using FESType = H1<K, Scalar, Mesh>;

      /// @brief Left-hand side operand type.
      using LHSType = FunctionBase<LHSDerived>;

      /// @brief Right-hand side operand type.
      using RHSType = ShapeFunctionBase<ShapeFunction<RHSDerived, FESType, TestSpace>,
        FESType, TestSpace>;

      /// @brief Range type of the left-hand side operand.
      using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;
      /// @brief Range type of the right-hand side operand.
      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

      /// @brief Integrand expression type.
      using IntegrandType = ShapeFunctionBase<Dot<LHSType, RHSType>, FESType, TestSpace>;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;

      /// @brief Parent class type.
      using Parent = LinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);

      /**
       * @brief Constructs the integrator for the given integrand.
       * @param integrand Expression to integrate.
       */
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr),
          m_quadrature(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr),
          m_quadrature(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_quadrature(std::move(other.m_quadrature)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_vec(std::move(other.m_vec))
      {}

      /**
       * @brief Gets the integrand.
       * @returns The integrand.
       */
      constexpr const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      /**
       * @brief Returns the polytope the integrator is bound to.
       * @returns The polytope the integrator is bound to.
       */
      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      /**
       * @brief Binds the integrator to a polytope and tabulates the quadrature on it.
       * @param polytope Mesh entity used by this operation.
       * @returns Reference to this object after the operation.
       */
      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& f = integrand.getDerived().getLHS();
        const auto& fes = integrand.getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, idx);

        const size_t order = this->getOrder(polytope).value_or(
          integrand.getOrder(polytope).value_or(fe.getOrder()));

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set = true;
          m_order = order;
          m_geometry = geometry;

          m_qf = &QF::PolytopeQuadratureFormula::get(order, geometry);
        }

        assert(m_qf);
        m_quadrature = std::make_unique<Geometry::PolytopeQuadrature>(polytope, *m_qf);

        const size_t nte = integrand.getDOFs(polytope);

        const H1Element<K, ScalarType> scalarFE(geometry);

        m_vec.resize(static_cast<Eigen::Index>(nte));
        m_vec.setZero();

        assert(m_quadrature);
        const auto& q = *m_quadrature;
        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          const IntegrationPoint ip(q, qp);
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const LHSRangeType fval = f(ip);
          const auto& tab = scalarFE.getTabulation(*m_qf);

          if constexpr (std::is_same_v<RHSRangeType, ScalarType>)
          {
            for (size_t local = 0; local < nte; ++local)
              m_vec(local) += wdet * fval * tab.getBasis(qp, local);
          }
          else if constexpr (FormLanguage::IsVectorRange<RHSRangeType>::Value)
          {
            const size_t vdim = fes.getVectorDimension();
            assert(nte == scalarFE.getCount() * vdim);

            for (size_t local = 0; local < nte; ++local)
            {
              const size_t scalarLocal = local / vdim;
              const size_t comp = local % vdim;
              m_vec(local) += wdet * fval(comp) * tab.getBasis(qp, scalarLocal);
            }
          }
          else if constexpr (FormLanguage::IsMatrixRange<RHSRangeType>::Value)
          {
            const size_t components = fes.getVectorDimension();
            const size_t cols = fes.getColumns();
            assert(nte == scalarFE.getCount() * components);
            assert(fval.rows() == fes.getRows() && fval.cols() == cols);
            for (size_t local = 0; local < nte; ++local)
            {
              const size_t comp = local % components;
              m_vec(local) += wdet * fval(comp / cols, comp % cols) *
                tab.getBasis(qp, local / components);
            }
          }
          else
          {
            static_assert(std::is_same_v<RHSRangeType, ScalarType> ||
                FormLanguage::IsVectorRange<RHSRangeType>::Value ||
                FormLanguage::IsMatrixRange<RHSRangeType>::Value,
              "Unsupported H1 Integral(f,v) RHS range type. Expected scalar or "
              "vector-valued shape function.");
          }
        }

        return *this;
      }

      /**
       * @brief Returns an entry of the element vector.
       * @param local Index in the local numbering.
       * @returns Integral computed by the quadrature rule.
       */
      inline ScalarType integrate(size_t local) final override
      {
        return m_vec(local);
      }

      /**
       * @brief Returns the integration region.
       * @returns The integration region.
       */
      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      std::unique_ptr<Geometry::PolytopeQuadrature> m_quadrature;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      Math::Vector<ScalarType> m_vec;
  };

  /**
   * @brief Deduction guide for @c QuadratureRule.
   * @param integrand Expression to integrate.
   */
  template <size_t K, class LHSDerived, class RHSDerived, class Scalar, class Mesh>
  QuadratureRule(const ShapeFunctionBase<
    Dot<FunctionBase<LHSDerived>,
      ShapeFunctionBase<ShapeFunction<RHSDerived, H1<K, Scalar, Mesh>, TestSpace>,
        H1<K, Scalar, Mesh>, TestSpace>>,
    H1<K, Scalar, Mesh>, TestSpace>& integrand)
    -> QuadratureRule<ShapeFunctionBase<
      Dot<FunctionBase<LHSDerived>,
        ShapeFunctionBase<ShapeFunction<RHSDerived, H1<K, Scalar, Mesh>, TestSpace>,
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
   * ## Component addressing
   *
   * For @f$m@f$ components and @f$N_r@f$ scalar trial basis functions,
   * let @f$n_r=mN_r@f$ be the local trial DOF count. A matching-component
   * entry of the row-major operator has scalar address
   * @f[
   * j_c=(bm+c)n_r+(am+c)=bm n_r+am+c(n_r+1).
   * @f]
   * Real storage advances an integer address by @f$n_r+1@f$ between
   * components. Complex storage retains the original row/column expression:
   * equivalent integer formulas do not produce equally efficient kernels for
   * both scalar representations. The component sequence and floating-point
   * accumulation remain unchanged. Mixed degrees retain the trial count as
   * their row stride; no coordinate matching or distributed operation is used.
   *
   * ## Empty component ranges
   *
   * With zero components, the local matrix has shape @f$0\times0@f$.
   * Binding retains the mapped-quadrature lifecycle but returns before
   * tabulating or contracting scalar bases. Equal empty DOF counts do not
   * establish equality of the trial and test scalar polynomial degrees.
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
  template <size_t KTrial, size_t KTest, class LHSDerived, class RHSDerived, class Scalar,
    class Mesh>
  class QuadratureRule<
    Dot<ShapeFunctionBase<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>,
          H1<KTrial, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
        H1<KTest, Scalar, Mesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<typename FormLanguage::Traits<Dot<
        ShapeFunctionBase<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>,
          H1<KTrial, Scalar, Mesh>, TrialSpace>,
        ShapeFunctionBase<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
          H1<KTest, Scalar, Mesh>, TestSpace>>>::ScalarType>
  {
    public:
      /// @brief Reports this handler as an optimized specialization.
      static constexpr bool Specialized = true;

      /// @brief Trial finite element space type.
      using TrialFESType = H1<KTrial, Scalar, Mesh>;
      /// @brief Test finite element space type.
      using TestFESType = H1<KTest, Scalar, Mesh>;

      /// @brief Left-hand side operand type.
      using LHSType =
        ShapeFunctionBase<ShapeFunction<LHSDerived, TrialFESType, TrialSpace>,
          TrialFESType, TrialSpace>;

      /// @brief Right-hand side operand type.
      using RHSType = ShapeFunctionBase<ShapeFunction<RHSDerived, TestFESType, TestSpace>,
        TestFESType, TestSpace>;

      /// @brief Integrand expression type.
      using IntegrandType = Dot<LHSType, RHSType>;
      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;
      /// @brief Parent class type.
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      /**
       * @brief Constructs the integrator for the given integrand.
       * @param integrand Expression to integrate.
       */
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr),
          m_quadrature(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr),
          m_quadrature(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_quadrature(std::move(other.m_quadrature)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_mat(std::move(other.m_mat))
      {}

      /**
       * @brief Gets the integrand.
       * @returns The integrand.
       */
      constexpr const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      /**
       * @brief Returns the polytope the integrator is bound to.
       * @returns The polytope the integrator is bound to.
       */
      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      /**
       * @brief Binds the integrator to a polytope and tabulates the quadrature on it.
       * @param polytope Mesh entity used by this operation.
       * @returns Reference to this object after the operation.
       */
      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();

        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe = testfes.getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe.getOrder();
        const size_t order = this->getOrder(polytope).value_or(
          getIntegrand().getOrder(polytope).value_or(k_tr + k_te));

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set = true;
          m_order = order;
          m_geometry = geometry;

          m_qf = &QF::PolytopeQuadratureFormula::get(order, geometry);
        }

        assert(m_qf);
        m_quadrature = std::make_unique<Geometry::PolytopeQuadrature>(polytope, *m_qf);

        const size_t ntr = lhs.getDOFs(polytope);
        const size_t nte = rhs.getDOFs(polytope);

        const H1Element<KTrial, ScalarType> trialScalarFE(geometry);
        const H1Element<KTest, ScalarType> testScalarFE(geometry);
        const size_t scalarCountTr = trialScalarFE.getCount();
        const size_t scalarCountTe = testScalarFE.getCount();

        assert(scalarCountTr > 0 && ntr % scalarCountTr == 0);
        assert(scalarCountTe > 0 && nte % scalarCountTe == 0);
        const size_t vdim = ntr / scalarCountTr;
        assert(vdim == nte / scalarCountTe);
        if constexpr (FormLanguage::IsMatrixRange<
                        typename FormLanguage::Traits<TrialFESType>::RangeType>::Value)
        {
          assert(trialfes.getRows() == testfes.getRows());
          assert(trialfes.getColumns() == testfes.getColumns());
        }

        const bool symmetric =
          (&trialfes.getMesh() == &testfes.getMesh()) && (ntr == nte);

        m_mat.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_mat.setZero();
        if (vdim == 0)
          return *this;

        const auto& trTab = trialScalarFE.getTabulation(*m_qf);
        const auto& teTab = testScalarFE.getTabulation(*m_qf);
        ScalarType* A = m_mat.data();

        assert(m_quadrature);
        const auto& q = *m_quadrature;
        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          if (symmetric)
          {
            for (size_t ib = 0; ib < scalarCountTe; ++ib)
            {
              const ScalarType phi_te = teTab.getBasis(qp, ib);

              {
                const ScalarType kii = wdet * phi_te * trTab.getBasis(qp, ib);
                size_t entry = (ib * vdim) * ntr + ib * vdim;
                for (size_t c = 0; c < vdim; ++c)
                {
                  const size_t index = std::is_same_v<ScalarType, Real>
                    ? entry
                    : (ib * vdim + c) * ntr + (ib * vdim + c);
                  A[index] += kii;
                  entry += ntr + 1;
                }
              }

              for (size_t ia = 0; ia < ib; ++ia)
              {
                const ScalarType kij = wdet * phi_te * trTab.getBasis(qp, ia);
                size_t entry = (ib * vdim) * ntr + ia * vdim;
                for (size_t c = 0; c < vdim; ++c)
                {
                  const size_t index = std::is_same_v<ScalarType, Real>
                    ? entry
                    : (ib * vdim + c) * ntr + (ia * vdim + c);
                  A[index] += kij;
                  entry += ntr + 1;
                }
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
                size_t entry = (ib * vdim) * ntr + ia * vdim;
                for (size_t c = 0; c < vdim; ++c)
                {
                  const size_t index = std::is_same_v<ScalarType, Real>
                    ? entry
                    : (ib * vdim + c) * ntr + (ia * vdim + c);
                  A[index] += kij;
                  entry += ntr + 1;
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

      /**
       * @brief Returns an entry of the element matrix.
       * @param tr Trial shape-function expression.
       * @param te Test shape-function expression.
       * @returns Integral computed by the quadrature rule.
       */
      inline ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_mat(te, tr);
      }

      /**
       * @brief Returns the integration region.
       * @returns The integration region.
       */
      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      std::unique_ptr<Geometry::PolytopeQuadrature> m_quadrature;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_mat;
  };

  /**
   * @brief Deduction guide for @c QuadratureRule.
   * @param integrand Expression to integrate.
   */
  template <size_t KTrial, size_t KTest, class LHSDerived, class RHSDerived, class Scalar,
    class Mesh>
  QuadratureRule(const Dot<
    ShapeFunctionBase<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>,
      H1<KTrial, Scalar, Mesh>, TrialSpace>,
    ShapeFunctionBase<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
      H1<KTest, Scalar, Mesh>, TestSpace>>& integrand)
    -> QuadratureRule<Dot<
      ShapeFunctionBase<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>,
        H1<KTrial, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
        H1<KTest, Scalar, Mesh>, TestSpace>>>;

  /**
   * @internal
   * @brief Certifies the plain-leaf operands of an H1 weighted mass expression.
   *
   * Only a direct H1 trial leaf has the component support used by the mass
   * rule. Nested products, traces and transformed shapes retain generic
   * integration. This metadata neither rewrites nor owns an expression.
   *
   * @tparam Expression Multiplication expression to classify.
   */
  template <class Expression>
  struct H1MassOperands
  {
    /// @brief Whether the expression has the required component structure.
    static constexpr bool Value = false;
  };

  /**
   * @internal
   * @brief Retains the existing left-coefficient mass applicability.
   * @tparam K Trial degree.
   * @tparam CoefficientDerived Coefficient expression.
   * @tparam Derived Plain trial-leaf reference.
   * @tparam Range Trial field range.
   * @tparam Mesh Mesh type.
   */
  template <size_t K, class CoefficientDerived, class Derived, class Range, class Mesh>
  struct H1MassOperands<Mult<FunctionBase<CoefficientDerived>,
    ShapeFunctionBase<ShapeFunction<Derived, H1<K, Range, Mesh>, TrialSpace>,
      H1<K, Range, Mesh>, TrialSpace>>>
  {
    /// @brief Original coefficient operand type.
    using CoefficientType = FunctionBase<CoefficientDerived>;
    /// @brief Direct H1 trial operand type.
    using MultiplicandType =
      ShapeFunctionBase<ShapeFunction<Derived, H1<K, Range, Mesh>, TrialSpace>,
        H1<K, Range, Mesh>, TrialSpace>;
    /// @brief Whether the original left-coefficient rule accepts this range.
    static constexpr bool Value = !FormLanguage::IsTensorRange<
      typename FormLanguage::Traits<CoefficientType>::RangeType>::Value;
    /// @brief Whether the coefficient precedes the trial operand.
    static constexpr bool CoefficientOnLeft = true;
  };

  /**
   * @internal
   * @brief Certifies scalar right scaling of aggregate H1 carriers.
   *
   * Scalar-valued leaves retain the generic rule, which reuses each scaled
   * trial value across every test entry. Component separation is reserved for
   * vector and matrix carriers; no second scalar cache is introduced.
   *
   * @tparam K Trial degree.
   * @tparam Derived Plain trial-leaf reference.
   * @tparam Range Trial field range.
   * @tparam Mesh Mesh type.
   * @tparam CoefficientDerived Coefficient expression.
   */
  template <size_t K, class Derived, class Range, class Mesh, class CoefficientDerived>
  struct H1MassOperands<
    Mult<ShapeFunctionBase<ShapeFunction<Derived, H1<K, Range, Mesh>, TrialSpace>,
           H1<K, Range, Mesh>, TrialSpace>,
      FunctionBase<CoefficientDerived>>>
  {
    /// @brief Original coefficient operand type.
    using CoefficientType = FunctionBase<CoefficientDerived>;
    /// @brief Direct H1 trial operand type.
    using MultiplicandType =
      ShapeFunctionBase<ShapeFunction<Derived, H1<K, Range, Mesh>, TrialSpace>,
        H1<K, Range, Mesh>, TrialSpace>;
    /// @brief Whether the scalar coefficient acts on an aggregate carrier.
    static constexpr bool Value = (FormLanguage::IsVectorRange<Range>::Value ||
                                    FormLanguage::IsMatrixRange<Range>::Value) &&
      (std::is_same_v<typename FormLanguage::Traits<CoefficientType>::RangeType,
         Real> ||
        std::is_same_v<typename FormLanguage::Traits<CoefficientType>::RangeType,
          Complex>);
    /// @brief Whether the coefficient precedes the trial operand.
    static constexpr bool CoefficientOnLeft = false;
  };

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Weighted mass quadrature for plain H1 trial and test leaves.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int (A u) \cdot v \ dx \: ,
   * @f]
   * where @f$u@f$ and @f$v@f$ are H1 shape functions and @f$A@f$ is a coefficient function.
   *
   * ## Architecture
   *
   * Mapped quadrature samples and the completed local matrix are owned by the
   * rule. Reference tables are borrowed only during binding, after the current
   * coefficient value has been materialized. Public entry access retains test
   * rows and trial columns independently of the private storage order.
   *
   * ### Scalar coefficients
   *
   * A scalar coefficient gives @f$M_{bm+e,am+c}=\delta_{ec}S_{ba}@f$.
   * For left scaling, row-major addresses of successive matching components
   * differ by
   * @f$n_r+1@f$, where @f$n_r@f$ counts trial DOFs, including for rectangular
   * mixed-degree operators. Scalar fields have @f$m=1@f$ by their range type
   * and update @f$M_{b,a}@f$ directly, without runtime component traversal.
   * Right scaling uses trial-major traversal and column-major storage, with
   * matching-component stride @f$n_t+1@f$, where @f$n_t@f$ counts test DOFs.
   * For aggregate ranges, real storage uses the corresponding integer stride;
   * complex storage retains explicit row/column addressing. The public entry
   * convention remains test rows and trial columns, including mixed degrees.
   * Matrix-valued coefficients remain coupled and do not use this component
   * diagonal structure.
   * Scalar right scaling retains the expression @f$u\rho@f$. Its contribution
   * is evaluated in the original order
   * @f$w_qJ_q[(\phi_a\rho_q)\overline{\phi_b}]@f$, not by commuting or
   * regrouping the coefficient. Right-side reference values retain the
   * finite-element scalar type before coefficient promotion. Coefficients
   * are evaluated before either reference table is borrowed.
   * Each scaled right-side trial value is reused across the test bases at the
   * current sample. Fixed one-element views use the existing
   * @ref Rodin::Math::dot "dot product" without constructing component-expanded
   * values. The complex
   * weighted contribution is formed once before component replication. Real
   * multiplication remains in the entry update to preserve the established
   * floating-point association. Quadrature accumulation order is unchanged,
   * and no additional persistent cache or assembly-backend condition is introduced.
   * Right-scaled scalar fields retain generic cached-trial quadrature; this
   * component rule accepts right scaling only for vector and matrix carriers.
   *
   * @tparam KTrial Trial polynomial degree.
   * @tparam KTest Test polynomial degree.
   * @tparam LeftOperand Original left multiplication operand.
   * @tparam RightOperand Original right multiplication operand.
   * @tparam RHSDerived Plain test-leaf reference.
   * @tparam Scalar Finite-element field range.
   * @tparam Mesh Mesh type.
   *
   * ## Judgement
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int (A u) \cdot v \ dx : \texttt{QuadratureRule}}
   * {\vdash u : \texttt{H1}<K_{\mathrm{trial}}>, \ \vdash v : \texttt{H1}<K_{\mathrm{test}}>}
   * @f]
   */
  template <size_t KTrial, size_t KTest, class LeftOperand, class RightOperand,
    class RHSDerived, class Scalar, class Mesh>
    requires(H1MassOperands<Mult<LeftOperand, RightOperand>>::Value &&
      std::is_same_v<typename H1MassOperands<
                       Mult<LeftOperand, RightOperand>>::MultiplicandType::FESType,
        H1<KTrial, Scalar, Mesh>>)
  class QuadratureRule<Dot<ShapeFunctionBase<Mult<LeftOperand, RightOperand>,
                             H1<KTrial, Scalar, Mesh>, TrialSpace>,
    ShapeFunctionBase<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
      H1<KTest, Scalar, Mesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<typename FormLanguage::Traits<
        Dot<ShapeFunctionBase<Mult<LeftOperand, RightOperand>, H1<KTrial, Scalar, Mesh>,
              TrialSpace>,
          ShapeFunctionBase<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
            H1<KTest, Scalar, Mesh>, TestSpace>>>::ScalarType>
  {
    public:
      /// @brief Reports this handler as an optimized specialization.
      static constexpr bool Specialized = true;

      /// @brief Trial finite element space type.
      using TrialFESType = H1<KTrial, Scalar, Mesh>;
      /// @brief Test finite element space type.
      using TestFESType = H1<KTest, Scalar, Mesh>;

      /// @brief Certified operand positions without expression normalization.
      using Operands = H1MassOperands<Mult<LeftOperand, RightOperand>>;

      /// @brief Coefficient type appearing in the integrand.
      using CoefficientType = typename Operands::CoefficientType;
      /// @brief Range type of the coefficient.
      using CoefficientRangeType =
        typename FormLanguage::Traits<CoefficientType>::RangeType;

      /// @brief Type of the multiplicand in the integrand.
      using MultiplicandType = typename Operands::MultiplicandType;

      /// @brief Left-hand side operand type.
      using LHSType =
        ShapeFunctionBase<Mult<LeftOperand, RightOperand>, TrialFESType, TrialSpace>;

      /// @brief Right-hand side operand type.
      using RHSType = ShapeFunctionBase<ShapeFunction<RHSDerived, TestFESType, TestSpace>,
        TestFESType, TestSpace>;

      /// @brief Integrand expression type.
      using IntegrandType = Dot<LHSType, RHSType>;
      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;
      /// @brief Parent class type.
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      /**
       * @brief Constructs the integrator for the given integrand.
       * @param integrand Expression to integrate.
       */
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr),
          m_quadrature(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr),
          m_quadrature(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_quadrature(std::move(other.m_quadrature)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_cmv(std::move(other.m_cmv)),
          m_mat(std::move(other.m_mat))
      {}

      /**
       * @brief Gets the integrand.
       * @returns The integrand.
       */
      constexpr const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      /**
       * @brief Returns the polytope the integrator is bound to.
       * @returns The polytope the integrator is bound to.
       */
      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      /**
       * @brief Binds the integrator to a polytope and tabulates the quadrature on it.
       * @param polytope Mesh entity used by this operation.
       * @returns Reference to this object after the operation.
       */
      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& lhs = integrand.getLHS();
        const auto& coeff = [&]() -> const auto& {
          if constexpr (Operands::CoefficientOnLeft)
            return lhs.getDerived().getLHS();
          else
            return lhs.getDerived().getRHS();
        }();
        const auto& rhs = integrand.getRHS();

        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe = testfes.getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe.getOrder();
        const size_t order = this->getOrder(polytope).value_or(
          getIntegrand().getOrder(polytope).value_or(k_tr + k_te));

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set = true;
          m_order = order;
          m_geometry = geometry;

          m_qf = &QF::PolytopeQuadratureFormula::get(order, geometry);
        }

        assert(m_qf);
        m_quadrature = std::make_unique<Geometry::PolytopeQuadrature>(polytope, *m_qf);

        const size_t ntr = lhs.getDOFs(polytope);
        const size_t nte = rhs.getDOFs(polytope);

        using ReferenceScalar = std::conditional_t<Operands::CoefficientOnLeft,
          ScalarType, typename FormLanguage::Traits<TrialFESType>::ScalarType>;
        const H1Element<KTrial, ReferenceScalar> trialScalarFE(geometry);
        const H1Element<KTest, ReferenceScalar> testScalarFE(geometry);
        const size_t scalarCountTr = trialScalarFE.getCount();
        const size_t scalarCountTe = testScalarFE.getCount();

        assert(scalarCountTr > 0 && ntr % scalarCountTr == 0);
        assert(scalarCountTe > 0 && nte % scalarCountTe == 0);
        const size_t vdim = ntr / scalarCountTr;
        assert(vdim == nte / scalarCountTe);
        if constexpr (FormLanguage::IsMatrixRange<
                        typename FormLanguage::Traits<TrialFESType>::RangeType>::Value)
        {
          assert(trialfes.getRows() == testfes.getRows());
          assert(trialfes.getColumns() == testfes.getColumns());
        }

        m_mat.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_mat.setZero();
        ScalarType* A = m_mat.data();

        assert(m_quadrature);
        const auto& q = *m_quadrature;
        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          const IntegrationPoint ip(q, qp);
          const ScalarType wdet = [&]() {
            if constexpr (Operands::CoefficientOnLeft)
              return static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());
            else
              return static_cast<ScalarType>(m_qf->getWeight(qp)) *
                static_cast<ScalarType>(p.getDistortion());
          }();

          const CoefficientRangeType coefficientValue = coeff.getValue(ip);
          if constexpr (!Operands::CoefficientOnLeft)
          {
            if (vdim == 0)
              continue;
          }
          const auto& trTab = trialScalarFE.getTabulation(*m_qf);
          const auto& teTab = testScalarFE.getTabulation(*m_qf);

          if constexpr (std::is_same_v<CoefficientRangeType, ScalarType> ||
            std::is_same_v<CoefficientRangeType, Real>)
          {
            using WeightScalar = std::conditional_t<Operands::CoefficientOnLeft,
              ScalarType, CoefficientRangeType>;
            const WeightScalar csv = coefficientValue;
            if constexpr (!Operands::CoefficientOnLeft)
            {
              using ScalarView = Eigen::Map<const Eigen::Matrix<ScalarType, 1, 1>>;
              for (size_t ia = 0; ia < scalarCountTr; ++ia)
              {
                const auto value = trTab.getBasis(qp, ia) * csv;
                for (size_t ib = 0; ib < scalarCountTe; ++ib)
                {
                  const auto phi_te = teTab.getBasis(qp, ib);
                  const ScalarType contraction =
                    Math::dot(ScalarView(&value), ScalarView(&phi_te));
                  const ScalarType contribution = [&]() {
                    if constexpr (std::is_same_v<ScalarType, Complex>)
                      return wdet * contraction;
                    else
                      return contraction;
                  }();
                  size_t entry = (ia * vdim) * nte + ib * vdim;
                  for (size_t c = 0; c < vdim; ++c)
                  {
                    const size_t index = std::is_same_v<ScalarType, Real>
                      ? entry
                      : (ia * vdim + c) * nte + (ib * vdim + c);
                    if constexpr (std::is_same_v<ScalarType, Complex>)
                      A[index] += contribution;
                    else
                      A[index] += wdet * contraction;
                    entry += nte + 1;
                  }
                }
              }
            }
            else
            {
              for (size_t ib = 0; ib < scalarCountTe; ++ib)
              {
                const auto phi_te = teTab.getBasis(qp, ib);
                for (size_t ia = 0; ia < scalarCountTr; ++ia)
                {
                  const ScalarType kij = wdet * csv * phi_te * trTab.getBasis(qp, ia);
                  if constexpr (std::is_same_v<Scalar, Real> ||
                    std::is_same_v<Scalar, Complex>)
                    A[ib * ntr + ia] += kij;
                  else
                  {
                    size_t entry = (ib * vdim) * ntr + ia * vdim;
                    for (size_t c = 0; c < vdim; ++c)
                    {
                      const size_t index = std::is_same_v<ScalarType, Real>
                        ? entry
                        : (ib * vdim + c) * ntr + (ia * vdim + c);
                      A[index] += kij;
                      entry += ntr + 1;
                    }
                  }
                }
              }
            }
          }
          else if constexpr (FormLanguage::IsMatrixRange<CoefficientRangeType>::Value)
          {
            m_cmv = coefficientValue;
            for (size_t ib = 0; ib < scalarCountTe; ++ib)
            {
              const ScalarType phi_te = teTab.getBasis(qp, ib);
              for (size_t ia = 0; ia < scalarCountTr; ++ia)
              {
                const ScalarType basisProduct = wdet * phi_te * trTab.getBasis(qp, ia);
                if constexpr (FormLanguage::IsMatrixRange<typename FormLanguage::Traits<
                                TrialFESType>::RangeType>::Value)
                {
                  const size_t rows = trialfes.getRows(), cols = trialfes.getColumns();
                  assert(m_cmv.rows() == rows && m_cmv.cols() == rows);
                  for (size_t column = 0; column < cols; ++column)
                  {
                    for (size_t rowTest = 0; rowTest < rows; ++rowTest)
                    {
                      for (size_t rowTrial = 0; rowTrial < rows; ++rowTrial)
                      {
                        A[(ib * vdim + rowTest * cols + column) * ntr +
                          (ia * vdim + rowTrial * cols + column)] +=
                          basisProduct * m_cmv(rowTest, rowTrial);
                      }
                    }
                  }
                }
                else
                {
                  for (size_t dd = 0; dd < vdim; ++dd)
                  {
                    for (size_t cc = 0; cc < vdim; ++cc)
                    {
                      A[(ib * vdim + dd) * ntr + (ia * vdim + cc)] +=
                        basisProduct * m_cmv(dd, cc);
                    }
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

        return *this;
      }

      /**
       * @brief Returns an entry of the element matrix.
       * @param tr Trial shape-function expression.
       * @param te Test shape-function expression.
       * @returns Integral computed by the quadrature rule.
       */
      inline ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_mat(te, tr);
      }

      /**
       * @brief Returns the integration region.
       * @returns The integration region.
       */
      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      std::unique_ptr<Geometry::PolytopeQuadrature> m_quadrature;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      Math::SpatialMatrix<ScalarType> m_cmv;
      Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic,
        Operands::CoefficientOnLeft ? Eigen::RowMajor : Eigen::ColMajor>
        m_mat;
  };

  /**
   * @brief Deduction guide for @c QuadratureRule.
   * @param integrand Expression to integrate.
   */
  template <size_t KTrial, size_t KTest, class CoefficientDerived, class LHSDerived,
    class RHSDerived, class Scalar, class Mesh>
  QuadratureRule(const Dot<
    ShapeFunctionBase<
      Mult<FunctionBase<CoefficientDerived>,
        ShapeFunctionBase<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>,
          H1<KTrial, Scalar, Mesh>, TrialSpace>>,
      H1<KTrial, Scalar, Mesh>, TrialSpace>,
    ShapeFunctionBase<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
      H1<KTest, Scalar, Mesh>, TestSpace>>& integrand)
    -> QuadratureRule<
      Dot<ShapeFunctionBase<Mult<FunctionBase<CoefficientDerived>,
                              ShapeFunctionBase<ShapeFunction<LHSDerived,
                                                  H1<KTrial, Scalar, Mesh>, TrialSpace>,
                                H1<KTrial, Scalar, Mesh>, TrialSpace>>,
            H1<KTrial, Scalar, Mesh>, TrialSpace>,
        ShapeFunctionBase<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
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
  template <size_t KTrial, size_t KTest, class CoefficientDerived, class LHSDerived,
    class RHSDerived, class Scalar, class Mesh>
    requires(((!FormLanguage::IsMatrixRange<Scalar>::Value)) &&
      (!FormLanguage::IsTensorRange<typename FormLanguage::Traits<
          FunctionBase<CoefficientDerived>>::RangeType>::Value))
  class QuadratureRule<Dot<
    ShapeFunctionBase<Mult<FunctionBase<CoefficientDerived>,
                        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived,
                                            H1<KTrial, Scalar, Mesh>, TrialSpace>>,
                          H1<KTrial, Scalar, Mesh>, TrialSpace>>,
      H1<KTrial, Scalar, Mesh>, TrialSpace>,
    ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
      H1<KTest, Scalar, Mesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<typename FormLanguage::Traits<Dot<
        ShapeFunctionBase<Mult<FunctionBase<CoefficientDerived>,
                            ShapeFunctionBase<Grad<ShapeFunction<LHSDerived,
                                                H1<KTrial, Scalar, Mesh>, TrialSpace>>,
                              H1<KTrial, Scalar, Mesh>, TrialSpace>>,
          H1<KTrial, Scalar, Mesh>, TrialSpace>,
        ShapeFunctionBase<
          Grad<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
          H1<KTest, Scalar, Mesh>, TestSpace>>>::ScalarType>
  {
    public:
      /// @brief Reports this handler as an optimized specialization.
      static constexpr bool Specialized = true;

      /// @brief Trial finite element space type.
      using TrialFESType = H1<KTrial, Scalar, Mesh>;
      /// @brief Test finite element space type.
      using TestFESType = H1<KTest, Scalar, Mesh>;

      /// @brief Coefficient type appearing in the integrand.
      using CoefficientType = FunctionBase<CoefficientDerived>;
      /// @brief Range type of the coefficient.
      using CoefficientRangeType =
        typename FormLanguage::Traits<CoefficientType>::RangeType;

      /// @brief Type of the multiplicand in the integrand.
      using MultiplicandType =
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, TrialFESType, TrialSpace>>,
          TrialFESType, TrialSpace>;

      /// @brief Left-hand side operand type.
      using LHSType = ShapeFunctionBase<Mult<CoefficientType, MultiplicandType>,
        TrialFESType, TrialSpace>;

      /// @brief Right-hand side operand type.
      using RHSType =
        ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, TestFESType, TestSpace>>,
          TestFESType, TestSpace>;

      /// @brief Integrand expression type.
      using IntegrandType = Dot<LHSType, RHSType>;
      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;
      /// @brief Parent class type.
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      /**
       * @brief Constructs the integrator for the given integrand.
       * @param integrand Expression to integrate.
       */
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr),
          m_quadrature(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr),
          m_quadrature(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_quadrature(std::move(other.m_quadrature)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_cmv(std::move(other.m_cmv)),
          m_mat(std::move(other.m_mat))
      {}

      /**
       * @brief Gets the integrand.
       * @returns The integrand.
       */
      constexpr const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      /**
       * @brief Returns the polytope the integrator is bound to.
       * @returns The polytope the integrator is bound to.
       */
      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      /**
       * @brief Binds the integrator to a polytope and tabulates the quadrature on it.
       * @param polytope Mesh entity used by this operation.
       * @returns Reference to this object after the operation.
       */
      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& lhs = integrand.getLHS();
        const auto& coeff = lhs.getDerived().getLHS();
        const auto& rhs = integrand.getRHS();

        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe = testfes.getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe.getOrder();
        const size_t order =
          this->getOrder(polytope).value_or(getIntegrand().getOrder(polytope).value_or(
            ((k_tr == 0 || k_te == 0) ? 0 : (k_tr + k_te - 2))));

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set = true;
          m_order = order;
          m_geometry = geometry;

          m_qf = &QF::PolytopeQuadratureFormula::get(order, geometry);
        }

        assert(m_qf);
        m_quadrature = std::make_unique<Geometry::PolytopeQuadrature>(polytope, *m_qf);

        const size_t ntr = lhs.getDOFs(polytope);
        const size_t nte = rhs.getDOFs(polytope);

        m_mat.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_mat.setZero();
        ScalarType* A_data = m_mat.data();

        static thread_local std::vector<Math::SpatialVector<ScalarType>> Gtr;
        static thread_local std::vector<Math::SpatialVector<ScalarType>> Gte;

        if (Gtr.size() < ntr)
          Gtr.resize(ntr);
        if (Gte.size() < nte)
          Gte.resize(nte);
        for (size_t a = 0; a < ntr; ++a)
          Gtr[a].resize(static_cast<std::uint8_t>(d));
        for (size_t b = 0; b < nte; ++b)
          Gte[b].resize(static_cast<std::uint8_t>(d));

        assert(m_quadrature);
        const auto& q = *m_quadrature;
        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          const IntegrationPoint ip(q, qp);
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto Jinv = p.getJacobianInverse();

          const CoefficientRangeType coefficientValue = coeff.getValue(ip);
          const auto& trTab = trialfe.getTabulation(*m_qf);
          const auto& teTab = testfe.getTabulation(*m_qf);

          // Build physical gradients for trial and test
          if (d == 3)
          {
            const ScalarType a00 = Jinv(0, 0), a10 = Jinv(1, 0), a20 = Jinv(2, 0);
            const ScalarType a01 = Jinv(0, 1), a11 = Jinv(1, 1), a21 = Jinv(2, 1);
            const ScalarType a02 = Jinv(0, 2), a12 = Jinv(1, 2), a22 = Jinv(2, 2);

            for (size_t a = 0; a < ntr; ++a)
            {
              const auto g = trTab.getGradient(qp, a);
              const ScalarType gx = g[0], gy = g[1], gz = g[2];
              Gtr[a][0] = a00 * gx + a10 * gy + a20 * gz;
              Gtr[a][1] = a01 * gx + a11 * gy + a21 * gz;
              Gtr[a][2] = a02 * gx + a12 * gy + a22 * gz;
            }
            for (size_t b = 0; b < nte; ++b)
            {
              const auto g = teTab.getGradient(qp, b);
              const ScalarType gx = g[0], gy = g[1], gz = g[2];
              Gte[b][0] = a00 * gx + a10 * gy + a20 * gz;
              Gte[b][1] = a01 * gx + a11 * gy + a21 * gz;
              Gte[b][2] = a02 * gx + a12 * gy + a22 * gz;
            }
          }
          else if (d == 2)
          {
            const ScalarType a00 = Jinv(0, 0), a10 = Jinv(1, 0);
            const ScalarType a01 = Jinv(0, 1), a11 = Jinv(1, 1);

            for (size_t a = 0; a < ntr; ++a)
            {
              const auto g = trTab.getGradient(qp, a);
              const ScalarType gx = g[0], gy = g[1];
              Gtr[a][0] = a00 * gx + a10 * gy;
              Gtr[a][1] = a01 * gx + a11 * gy;
            }
            for (size_t b = 0; b < nte; ++b)
            {
              const auto g = teTab.getGradient(qp, b);
              const ScalarType gx = g[0], gy = g[1];
              Gte[b][0] = a00 * gx + a10 * gy;
              Gte[b][1] = a01 * gx + a11 * gy;
            }
          }
          else if (d == 1)
          {
            const ScalarType a00 = Jinv(0, 0);
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

          if constexpr (std::is_same_v<CoefficientRangeType, ScalarType> ||
            std::is_same_v<CoefficientRangeType, Real>)
          {
            const ScalarType csv = coefficientValue;
            for (size_t b = 0; b < nte; ++b)
            {
              ScalarType* row = A_data + b * ntr;
              const auto& gb = Gte[b];
              for (size_t a = 0; a < ntr; ++a)
                row[a] += wdet * csv * Math::dot(gb, Gtr[a]);
            }
          }
          else if constexpr (FormLanguage::IsMatrixRange<CoefficientRangeType>::Value)
          {
            m_cmv = coefficientValue;
            for (size_t b = 0; b < nte; ++b)
            {
              ScalarType* row = A_data + b * ntr;
              const auto& gb = Gte[b];
              for (size_t a = 0; a < ntr; ++a)
              {
                Math::SpatialVector<ScalarType> aGtrA = m_cmv * Gtr[a];
                row[a] += wdet * Math::dot(gb, aGtrA);
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

      /**
       * @brief Returns an entry of the element matrix.
       * @param tr Trial shape-function expression.
       * @param te Test shape-function expression.
       * @returns Integral computed by the quadrature rule.
       */
      inline ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_mat(te, tr);
      }

      /**
       * @brief Returns the integration region.
       * @returns The integration region.
       */
      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      std::unique_ptr<Geometry::PolytopeQuadrature> m_quadrature;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      Math::SpatialMatrix<ScalarType> m_cmv;
      Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_mat;
  };

  /**
   * @brief Deduction guide for @c QuadratureRule.
   * @param integrand Expression to integrate.
   */
  template <size_t KTrial, size_t KTest, class CoefficientDerived, class LHSDerived,
    class RHSDerived, class Scalar, class Mesh>
  QuadratureRule(const Dot<
    ShapeFunctionBase<Mult<FunctionBase<CoefficientDerived>,
                        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived,
                                            H1<KTrial, Scalar, Mesh>, TrialSpace>>,
                          H1<KTrial, Scalar, Mesh>, TrialSpace>>,
      H1<KTrial, Scalar, Mesh>, TrialSpace>,
    ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
      H1<KTest, Scalar, Mesh>, TestSpace>>& integrand)
    -> QuadratureRule<
      Dot<ShapeFunctionBase<Mult<FunctionBase<CoefficientDerived>,
                              ShapeFunctionBase<Grad<ShapeFunction<LHSDerived,
                                                  H1<KTrial, Scalar, Mesh>, TrialSpace>>,
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
  template <size_t KTrial, size_t KTest, class CoefficientDerived, class LHSDerived,
    class RHSDerived, class Scalar, class Mesh>
    requires(!FormLanguage::IsTensorRange<typename FormLanguage::Traits<
               FunctionBase<CoefficientDerived>>::RangeType>::Value)
  class QuadratureRule<Mult<FunctionBase<CoefficientDerived>,
    Dot<ShapeFunctionBase<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>,
          H1<KTrial, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
        H1<KTest, Scalar, Mesh>, TestSpace>>>>
    : public LocalBilinearFormIntegratorBase<typename FormLanguage::Traits<Mult<
        FunctionBase<CoefficientDerived>,
        Dot<ShapeFunctionBase<
              ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>,
              H1<KTrial, Scalar, Mesh>, TrialSpace>,
          ShapeFunctionBase<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
            H1<KTest, Scalar, Mesh>, TestSpace>>>>::ScalarType>
  {
    public:
      /// @brief Reports this handler as an optimized specialization.
      static constexpr bool Specialized = true;

      /// @brief Trial finite element space type.
      using TrialFESType = H1<KTrial, Scalar, Mesh>;
      /// @brief Test finite element space type.
      using TestFESType = H1<KTest, Scalar, Mesh>;

      /// @brief Coefficient type appearing in the integrand.
      using CoefficientType = FunctionBase<CoefficientDerived>;

      /// @brief Left-hand side operand type.
      using LHSType =
        ShapeFunctionBase<ShapeFunction<LHSDerived, TrialFESType, TrialSpace>,
          TrialFESType, TrialSpace>;

      /// @brief Right-hand side operand type.
      using RHSType = ShapeFunctionBase<ShapeFunction<RHSDerived, TestFESType, TestSpace>,
        TestFESType, TestSpace>;

      /// @brief Type of the inner integrand.
      using InnerIntegrandType = Dot<LHSType, RHSType>;
      /// @brief Integrand expression type.
      using IntegrandType = Mult<CoefficientType, InnerIntegrandType>;
      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;
      /// @brief Parent class type.
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      /**
       * @brief Constructs the integrator for the given integrand.
       * @param integrand Expression to integrate.
       */
      QuadratureRule(const IntegrandType& integrand)
        : Parent(
            integrand.getRHS().getLHS().getLeaf(), integrand.getRHS().getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr),
          m_quadrature(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr),
          m_quadrature(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_quadrature(std::move(other.m_quadrature)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_mat(std::move(other.m_mat))
      {}

      /**
       * @brief Gets the integrand.
       * @returns The integrand.
       */
      constexpr const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      /**
       * @brief Returns the polytope the integrator is bound to.
       * @returns The polytope the integrator is bound to.
       */
      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      /**
       * @brief Binds the integrator to a polytope and tabulates the quadrature on it.
       * @param polytope Mesh entity used by this operation.
       * @returns Reference to this object after the operation.
       */
      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& coeff = integrand.getLHS();
        const auto& inner = integrand.getRHS();
        const auto& lhs = inner.getLHS();
        const auto& rhs = inner.getRHS();

        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe = testfes.getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe.getOrder();
        const size_t order = this->getOrder(polytope).value_or(
          getIntegrand().getOrder(polytope).value_or(k_tr + k_te));

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set = true;
          m_order = order;
          m_geometry = geometry;

          m_qf = &QF::PolytopeQuadratureFormula::get(order, geometry);
        }

        assert(m_qf);
        m_quadrature = std::make_unique<Geometry::PolytopeQuadrature>(polytope, *m_qf);

        const size_t ntr = lhs.getDOFs(polytope);
        const size_t nte = rhs.getDOFs(polytope);

        const H1Element<KTrial, ScalarType> trialScalarFE(geometry);
        const H1Element<KTest, ScalarType> testScalarFE(geometry);
        const size_t scalarCountTr = trialScalarFE.getCount();
        const size_t scalarCountTe = testScalarFE.getCount();

        assert(scalarCountTr > 0 && ntr % scalarCountTr == 0);
        assert(scalarCountTe > 0 && nte % scalarCountTe == 0);
        const size_t vdim = ntr / scalarCountTr;
        assert(vdim == nte / scalarCountTe);

        const auto& trTab = trialScalarFE.getTabulation(*m_qf);
        const auto& teTab = testScalarFE.getTabulation(*m_qf);

        m_mat.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_mat.setZero();
        ScalarType* A = m_mat.data();

        assert(m_quadrature);
        const auto& q = *m_quadrature;
        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          const IntegrationPoint ip(q, qp);
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const ScalarType csv = coeff(ip);

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

      /**
       * @brief Returns an entry of the element matrix.
       * @param tr Trial shape-function expression.
       * @param te Test shape-function expression.
       * @returns Integral computed by the quadrature rule.
       */
      inline ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_mat(te, tr);
      }

      /**
       * @brief Returns the integration region.
       * @returns The integration region.
       */
      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      std::unique_ptr<Geometry::PolytopeQuadrature> m_quadrature;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_mat;
  };

  /**
   * @brief Deduction guide for @c QuadratureRule.
   * @param integrand Expression to integrate.
   */
  template <size_t KTrial, size_t KTest, class CoefficientDerived, class LHSDerived,
    class RHSDerived, class Scalar, class Mesh>
  QuadratureRule(const Mult<FunctionBase<CoefficientDerived>,
    Dot<ShapeFunctionBase<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>,
          H1<KTrial, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
        H1<KTest, Scalar, Mesh>, TestSpace>>>& integrand)
    -> QuadratureRule<Mult<FunctionBase<CoefficientDerived>,
      Dot<
        ShapeFunctionBase<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>,
          H1<KTrial, Scalar, Mesh>, TrialSpace>,
        ShapeFunctionBase<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
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
  template <size_t KTrial, size_t KTest, class LHSDerived, class RHSDerived,
    class TrialRange, class TestRange, class Mesh>
    requires(!FormLanguage::IsMatrixRange<TrialRange>::Value &&
      !FormLanguage::IsMatrixRange<TestRange>::Value)
  class QuadratureRule<
    Dot<ShapeFunctionBase<
          Div<ShapeFunction<LHSDerived, H1<KTrial, TrialRange, Mesh>, TrialSpace>>,
          H1<KTrial, TrialRange, Mesh>, TrialSpace>,
      ShapeFunctionBase<ShapeFunction<RHSDerived, H1<KTest, TestRange, Mesh>, TestSpace>,
        H1<KTest, TestRange, Mesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<typename FormLanguage::Traits<
        Dot<ShapeFunctionBase<
              Div<ShapeFunction<LHSDerived, H1<KTrial, TrialRange, Mesh>, TrialSpace>>,
              H1<KTrial, TrialRange, Mesh>, TrialSpace>,
          ShapeFunctionBase<
            ShapeFunction<RHSDerived, H1<KTest, TestRange, Mesh>, TestSpace>,
            H1<KTest, TestRange, Mesh>, TestSpace>>>::ScalarType>
  {
    public:
      /// @brief Reports this handler as an optimized specialization.
      static constexpr bool Specialized = true;

      /// @brief Trial finite element space type.
      using TrialFESType = H1<KTrial, TrialRange, Mesh>;
      /// @brief Test finite element space type.
      using TestFESType = H1<KTest, TestRange, Mesh>;

      /// @brief Left-hand side operand type.
      using LHSType =
        ShapeFunctionBase<Div<ShapeFunction<LHSDerived, TrialFESType, TrialSpace>>,
          TrialFESType, TrialSpace>;

      /// @brief Right-hand side operand type.
      using RHSType = ShapeFunctionBase<ShapeFunction<RHSDerived, TestFESType, TestSpace>,
        TestFESType, TestSpace>;

      /// @brief Integrand expression type.
      using IntegrandType = Dot<LHSType, RHSType>;
      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;
      /// @brief Parent class type.
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      /**
       * @brief Constructs the integrator for the given integrand.
       * @param integrand Expression to integrate.
       */
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr),
          m_quadrature(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr),
          m_quadrature(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_quadrature(std::move(other.m_quadrature)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_mat(std::move(other.m_mat))
      {}

      /**
       * @brief Gets the integrand.
       * @returns The integrand.
       */
      constexpr const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      /**
       * @brief Returns the polytope the integrator is bound to.
       * @returns The polytope the integrator is bound to.
       */
      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      /**
       * @brief Binds the integrator to a polytope and tabulates the quadrature on it.
       * @param polytope Mesh entity used by this operation.
       * @returns Reference to this object after the operation.
       */
      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();

        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe = testfes.getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe.getOrder();
        const size_t order =
          this->getOrder(polytope).value_or(getIntegrand().getOrder(polytope).value_or(
            (k_tr == 0) ? 0 : (k_tr + k_te - 1)));

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set = true;
          m_order = order;
          m_geometry = geometry;

          m_qf = &QF::PolytopeQuadratureFormula::get(order, geometry);
        }

        assert(m_qf);
        m_quadrature = std::make_unique<Geometry::PolytopeQuadrature>(polytope, *m_qf);

        const size_t ntr = lhs.getDOFs(polytope);
        const size_t nte = rhs.getDOFs(polytope);

        const H1Element<KTrial, ScalarType> trialScalarFE(geometry);
        const H1Element<KTest, ScalarType> testScalarFE(geometry);
        const size_t scalarCountTr = trialScalarFE.getCount();

        assert(scalarCountTr > 0 && ntr % scalarCountTr == 0);
        const size_t vdim = ntr / scalarCountTr;

        const auto& trTabS = trialScalarFE.getTabulation(*m_qf);
        const auto& teTab = testScalarFE.getTabulation(*m_qf);

        m_mat.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_mat.setZero();
        ScalarType* A = m_mat.data();

        static thread_local std::vector<Math::SpatialVector<ScalarType>> GtrS;
        if (GtrS.size() < scalarCountTr)
          GtrS.resize(scalarCountTr);
        for (size_t a = 0; a < scalarCountTr; ++a)
          GtrS[a].resize(static_cast<std::uint8_t>(d));

        assert(m_quadrature);
        const auto& q = *m_quadrature;
        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto Jinv = p.getJacobianInverse();

          // Physical gradients of scalar trial basis functions
          if (d == 3)
          {
            const ScalarType j00 = Jinv(0, 0), j10 = Jinv(1, 0), j20 = Jinv(2, 0);
            const ScalarType j01 = Jinv(0, 1), j11 = Jinv(1, 1), j21 = Jinv(2, 1);
            const ScalarType j02 = Jinv(0, 2), j12 = Jinv(1, 2), j22 = Jinv(2, 2);

            for (size_t a = 0; a < scalarCountTr; ++a)
            {
              const auto g = trTabS.getGradient(qp, a);
              const ScalarType gx = g[0], gy = g[1], gz = g[2];
              GtrS[a][0] = j00 * gx + j10 * gy + j20 * gz;
              GtrS[a][1] = j01 * gx + j11 * gy + j21 * gz;
              GtrS[a][2] = j02 * gx + j12 * gy + j22 * gz;
            }
          }
          else if (d == 2)
          {
            const ScalarType j00 = Jinv(0, 0), j10 = Jinv(1, 0);
            const ScalarType j01 = Jinv(0, 1), j11 = Jinv(1, 1);

            for (size_t a = 0; a < scalarCountTr; ++a)
            {
              const auto g = trTabS.getGradient(qp, a);
              const ScalarType gx = g[0], gy = g[1];
              GtrS[a][0] = j00 * gx + j10 * gy;
              GtrS[a][1] = j01 * gx + j11 * gy;
            }
          }
          else if (d == 1)
          {
            const ScalarType j00 = Jinv(0, 0);
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

      /**
       * @brief Returns an entry of the element matrix.
       * @param tr Trial shape-function expression.
       * @param te Test shape-function expression.
       * @returns Integral computed by the quadrature rule.
       */
      inline ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_mat(te, tr);
      }

      /**
       * @brief Returns the integration region.
       * @returns The integration region.
       */
      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      std::unique_ptr<Geometry::PolytopeQuadrature> m_quadrature;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_mat;
  };

  /**
   * @brief Deduction guide for @c QuadratureRule.
   * @param integrand Expression to integrate.
   */
  template <size_t KTrial, size_t KTest, class LHSDerived, class RHSDerived,
    class TrialRange, class TestRange, class Mesh>
  QuadratureRule(
    const Dot<ShapeFunctionBase<
                Div<ShapeFunction<LHSDerived, H1<KTrial, TrialRange, Mesh>, TrialSpace>>,
                H1<KTrial, TrialRange, Mesh>, TrialSpace>,
      ShapeFunctionBase<ShapeFunction<RHSDerived, H1<KTest, TestRange, Mesh>, TestSpace>,
        H1<KTest, TestRange, Mesh>, TestSpace>>& integrand)
    -> QuadratureRule<Dot<ShapeFunctionBase<Div<ShapeFunction<LHSDerived,
                                              H1<KTrial, TrialRange, Mesh>, TrialSpace>>,
                            H1<KTrial, TrialRange, Mesh>, TrialSpace>,
      ShapeFunctionBase<ShapeFunction<RHSDerived, H1<KTest, TestRange, Mesh>, TestSpace>,
        H1<KTest, TestRange, Mesh>, TestSpace>>>;

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
  template <size_t KTrial, size_t KTest, class LHSDerived, class RHSDerived,
    class TrialRange, class TestRange, class Mesh>
    requires(!FormLanguage::IsMatrixRange<TrialRange>::Value &&
      !FormLanguage::IsMatrixRange<TestRange>::Value)
  class QuadratureRule<Dot<
    ShapeFunctionBase<ShapeFunction<LHSDerived, H1<KTrial, TrialRange, Mesh>, TrialSpace>,
      H1<KTrial, TrialRange, Mesh>, TrialSpace>,
    ShapeFunctionBase<
      Div<ShapeFunction<RHSDerived, H1<KTest, TestRange, Mesh>, TestSpace>>,
      H1<KTest, TestRange, Mesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<typename FormLanguage::Traits<
        Dot<ShapeFunctionBase<
              ShapeFunction<LHSDerived, H1<KTrial, TrialRange, Mesh>, TrialSpace>,
              H1<KTrial, TrialRange, Mesh>, TrialSpace>,
          ShapeFunctionBase<
            Div<ShapeFunction<RHSDerived, H1<KTest, TestRange, Mesh>, TestSpace>>,
            H1<KTest, TestRange, Mesh>, TestSpace>>>::ScalarType>
  {
    public:
      /// @brief Reports this handler as an optimized specialization.
      static constexpr bool Specialized = true;

      /// @brief Trial finite element space type.
      using TrialFESType = H1<KTrial, TrialRange, Mesh>;
      /// @brief Test finite element space type.
      using TestFESType = H1<KTest, TestRange, Mesh>;

      /// @brief Left-hand side operand type.
      using LHSType =
        ShapeFunctionBase<ShapeFunction<LHSDerived, TrialFESType, TrialSpace>,
          TrialFESType, TrialSpace>;

      /// @brief Right-hand side operand type.
      using RHSType =
        ShapeFunctionBase<Div<ShapeFunction<RHSDerived, TestFESType, TestSpace>>,
          TestFESType, TestSpace>;

      /// @brief Integrand expression type.
      using IntegrandType = Dot<LHSType, RHSType>;
      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;
      /// @brief Parent class type.
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      /**
       * @brief Constructs the integrator for the given integrand.
       * @param integrand Expression to integrate.
       */
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr),
          m_quadrature(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr),
          m_quadrature(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_quadrature(std::move(other.m_quadrature)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_mat(std::move(other.m_mat))
      {}

      /**
       * @brief Gets the integrand.
       * @returns The integrand.
       */
      constexpr const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      /**
       * @brief Returns the polytope the integrator is bound to.
       * @returns The polytope the integrator is bound to.
       */
      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      /**
       * @brief Binds the integrator to a polytope and tabulates the quadrature on it.
       * @param polytope Mesh entity used by this operation.
       * @returns Reference to this object after the operation.
       */
      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();

        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe = testfes.getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe.getOrder();
        const size_t order =
          this->getOrder(polytope).value_or(getIntegrand().getOrder(polytope).value_or(
            (k_te == 0) ? 0 : (k_tr + k_te - 1)));

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set = true;
          m_order = order;
          m_geometry = geometry;

          m_qf = &QF::PolytopeQuadratureFormula::get(order, geometry);
        }

        assert(m_qf);
        m_quadrature = std::make_unique<Geometry::PolytopeQuadrature>(polytope, *m_qf);

        const size_t ntr = lhs.getDOFs(polytope);
        const size_t nte = rhs.getDOFs(polytope);

        const H1Element<KTrial, ScalarType> trialScalarFE(geometry);
        const H1Element<KTest, ScalarType> testScalarFE(geometry);
        const size_t scalarCountTe = testScalarFE.getCount();

        assert(scalarCountTe > 0 && nte % scalarCountTe == 0);
        const size_t vdim = nte / scalarCountTe;

        const auto& trTab = trialScalarFE.getTabulation(*m_qf);
        const auto& teTabS = testScalarFE.getTabulation(*m_qf);

        m_mat.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_mat.setZero();
        ScalarType* A = m_mat.data();

        static thread_local std::vector<Math::SpatialVector<ScalarType>> GteS;
        if (GteS.size() < scalarCountTe)
          GteS.resize(scalarCountTe);
        for (size_t b = 0; b < scalarCountTe; ++b)
          GteS[b].resize(static_cast<std::uint8_t>(d));

        assert(m_quadrature);
        const auto& q = *m_quadrature;
        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto Jinv = p.getJacobianInverse();

          // Physical gradients of scalar test basis functions
          if (d == 3)
          {
            const ScalarType j00 = Jinv(0, 0), j10 = Jinv(1, 0), j20 = Jinv(2, 0);
            const ScalarType j01 = Jinv(0, 1), j11 = Jinv(1, 1), j21 = Jinv(2, 1);
            const ScalarType j02 = Jinv(0, 2), j12 = Jinv(1, 2), j22 = Jinv(2, 2);

            for (size_t b = 0; b < scalarCountTe; ++b)
            {
              const auto g = teTabS.getGradient(qp, b);
              const ScalarType gx = g[0], gy = g[1], gz = g[2];
              GteS[b][0] = j00 * gx + j10 * gy + j20 * gz;
              GteS[b][1] = j01 * gx + j11 * gy + j21 * gz;
              GteS[b][2] = j02 * gx + j12 * gy + j22 * gz;
            }
          }
          else if (d == 2)
          {
            const ScalarType j00 = Jinv(0, 0), j10 = Jinv(1, 0);
            const ScalarType j01 = Jinv(0, 1), j11 = Jinv(1, 1);

            for (size_t b = 0; b < scalarCountTe; ++b)
            {
              const auto g = teTabS.getGradient(qp, b);
              const ScalarType gx = g[0], gy = g[1];
              GteS[b][0] = j00 * gx + j10 * gy;
              GteS[b][1] = j01 * gx + j11 * gy;
            }
          }
          else if (d == 1)
          {
            const ScalarType j00 = Jinv(0, 0);
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

      /**
       * @brief Returns an entry of the element matrix.
       * @param tr Trial shape-function expression.
       * @param te Test shape-function expression.
       * @returns Integral computed by the quadrature rule.
       */
      inline ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_mat(te, tr);
      }

      /**
       * @brief Returns the integration region.
       * @returns The integration region.
       */
      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      std::unique_ptr<Geometry::PolytopeQuadrature> m_quadrature;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_mat;
  };

  /**
   * @brief Deduction guide for @c QuadratureRule.
   * @param integrand Expression to integrate.
   */
  template <size_t KTrial, size_t KTest, class LHSDerived, class RHSDerived,
    class TrialRange, class TestRange, class Mesh>
  QuadratureRule(const Dot<
    ShapeFunctionBase<ShapeFunction<LHSDerived, H1<KTrial, TrialRange, Mesh>, TrialSpace>,
      H1<KTrial, TrialRange, Mesh>, TrialSpace>,
    ShapeFunctionBase<
      Div<ShapeFunction<RHSDerived, H1<KTest, TestRange, Mesh>, TestSpace>>,
      H1<KTest, TestRange, Mesh>, TestSpace>>& integrand)
    -> QuadratureRule<Dot<ShapeFunctionBase<ShapeFunction<LHSDerived,
                                              H1<KTrial, TrialRange, Mesh>, TrialSpace>,
                            H1<KTrial, TrialRange, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        Div<ShapeFunction<RHSDerived, H1<KTest, TestRange, Mesh>, TestSpace>>,
        H1<KTest, TestRange, Mesh>, TestSpace>>>;

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
   * ## Architecture
   *
   * Reference gradients have @f$d@f$ entries and are mapped through the
   * @f$d\times s@f$ geometric inverse into all @f$s@f$ ambient components.
   * A point has an empty reference derivative and a zero local operator.
   * Scalar coefficients permit plain-transpose mirroring for identical trial
   * and test bases, including complex scalars because the basis is real.
   * Their component separation gives
   * @f$A_{bm+e,am+c}=\delta_{ec}S_{ba}@f$. The first component block retains
   * the original quadrature accumulation order; its completed entries are
   * replicated into the other component-diagonal blocks at the final sample,
   * before mirroring. A zero contribution is omitted unless the final sample
   * must copy earlier accumulated entries. At index zero no earlier contributions
   * exist, so a zero first contribution requires no component copy.
   * An empty component range has a zero-sized operator and requires neither
   * coefficient evaluation nor scalar reference preparation.
   * Matrix coefficients act on field components and need not be symmetric;
   * their complete component blocks are retained without mirroring.
   * Equal degrees use the same scalar element type, formula and geometry, so
   * physical-gradient scratch is shared within the current quadrature point.
   * Otherwise separate trial/test storage is used. Gradients are recomputed
   * for every point and polytope; no table pointer or physical derivative is
   * cached across binds.
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
  template <size_t KTrial, size_t KTest, class CoefficientDerived, class LHSDerived,
    class RHSDerived, class Scalar, class Mesh>
    requires(!FormLanguage::IsTensorRange<typename FormLanguage::Traits<
               FunctionBase<CoefficientDerived>>::RangeType>::Value &&
      !FormLanguage::IsMatrixRange<Scalar>::Value)
  class QuadratureRule<
    Dot<ShapeFunctionBase<Mult<FunctionBase<CoefficientDerived>,
                            ShapeFunctionBase<Jacobian<ShapeFunction<LHSDerived,
                                                H1<KTrial, Scalar, Mesh>, TrialSpace>>,
                              H1<KTrial, Scalar, Mesh>, TrialSpace>>,
          H1<KTrial, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        Jacobian<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
        H1<KTest, Scalar, Mesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<typename FormLanguage::Traits<Dot<
        ShapeFunctionBase<Mult<FunctionBase<CoefficientDerived>,
                            ShapeFunctionBase<Jacobian<ShapeFunction<LHSDerived,
                                                H1<KTrial, Scalar, Mesh>, TrialSpace>>,
                              H1<KTrial, Scalar, Mesh>, TrialSpace>>,
          H1<KTrial, Scalar, Mesh>, TrialSpace>,
        ShapeFunctionBase<
          Jacobian<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
          H1<KTest, Scalar, Mesh>, TestSpace>>>::ScalarType>
  {
    public:
      /// @brief Reports this handler as an optimized specialization.
      static constexpr bool Specialized = true;

      /// @brief Trial finite element space type.
      using TrialFESType = H1<KTrial, Scalar, Mesh>;
      /// @brief Test finite element space type.
      using TestFESType = H1<KTest, Scalar, Mesh>;

      /// @brief Coefficient type appearing in the integrand.
      using CoefficientType = FunctionBase<CoefficientDerived>;
      /// @brief Range type of the coefficient.
      using CoefficientRangeType =
        typename FormLanguage::Traits<CoefficientType>::RangeType;

      /// @brief Type of the multiplicand in the integrand.
      using MultiplicandType =
        ShapeFunctionBase<Jacobian<ShapeFunction<LHSDerived, TrialFESType, TrialSpace>>,
          TrialFESType, TrialSpace>;

      /// @brief Left-hand side operand type.
      using LHSType = ShapeFunctionBase<Mult<CoefficientType, MultiplicandType>,
        TrialFESType, TrialSpace>;

      /// @brief Right-hand side operand type.
      using RHSType =
        ShapeFunctionBase<Jacobian<ShapeFunction<RHSDerived, TestFESType, TestSpace>>,
          TestFESType, TestSpace>;

      /// @brief Integrand expression type.
      using IntegrandType = Dot<LHSType, RHSType>;
      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;
      /// @brief Parent class type.
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      /**
       * @brief Constructs the integrator for the given integrand.
       * @param integrand Expression to integrate.
       */
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr),
          m_quadrature(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr),
          m_quadrature(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_quadrature(std::move(other.m_quadrature)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_cmv(std::move(other.m_cmv)),
          m_mat(std::move(other.m_mat))
      {}

      /**
       * @brief Gets the integrand.
       * @returns The integrand.
       */
      constexpr const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      /**
       * @brief Returns the polytope the integrator is bound to.
       * @returns The polytope the integrator is bound to.
       */
      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      /**
       * @brief Binds the integrator to a polytope and tabulates the quadrature on it.
       * @param polytope Mesh entity used by this operation.
       * @returns Reference to this object after the operation.
       */
      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& lhs = integrand.getLHS();
        const auto& coeff = lhs.getDerived().getLHS();
        const auto& rhs = integrand.getRHS();

        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe = testfes.getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe.getOrder();
        const size_t order =
          this->getOrder(polytope).value_or(getIntegrand().getOrder(polytope).value_or(
            ((k_tr == 0 || k_te == 0) ? 0 : (k_tr + k_te - 2))));

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set = true;
          m_order = order;
          m_geometry = geometry;

          m_qf = &QF::PolytopeQuadratureFormula::get(order, geometry);
        }

        assert(m_qf);
        m_quadrature = std::make_unique<Geometry::PolytopeQuadrature>(polytope, *m_qf);

        const size_t ntr = lhs.getDOFs(polytope);
        const size_t nte = rhs.getDOFs(polytope);

        const H1Element<KTrial, ScalarType> trialScalarElement(geometry);
        const H1Element<KTest, ScalarType> testScalarElement(geometry);
        const size_t trialScalarCount = trialScalarElement.getCount();
        const size_t testScalarCount = testScalarElement.getCount();

        assert(trialScalarCount > 0);
        assert(testScalarCount > 0);
        assert(ntr % trialScalarCount == 0);
        assert(nte % testScalarCount == 0);
        const size_t vdim = ntr / trialScalarCount;
        assert(vdim == nte / testScalarCount);

        const bool symmetric =
          !FormLanguage::IsMatrixRange<CoefficientRangeType>::Value && KTrial == KTest &&
          (&trialfes.getMesh() == &testfes.getMesh()) && (ntr == nte);

        m_mat.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_mat.setZero();
        if (d == 0 || ntr == 0 || nte == 0)
          return *this;
        const size_t spaceDimension = polytope.getMesh().getSpaceDimension();
        ScalarType* A_data = m_mat.data();

        constexpr bool sharedGradients = KTrial == KTest;

        static thread_local std::vector<Math::SpatialVector<ScalarType>> GtrS;
        static thread_local std::vector<Math::SpatialVector<ScalarType>> GteStorage;
        auto& GteS = sharedGradients ? GtrS : GteStorage;

        if (GtrS.size() < trialScalarCount)
          GtrS.resize(trialScalarCount);
        if (GteS.size() < testScalarCount)
          GteS.resize(testScalarCount);

        assert(m_quadrature);
        const auto& q = *m_quadrature;
        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          const IntegrationPoint ip(q, qp);
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto Jinv = p.getJacobianInverse();

          const CoefficientRangeType coefficientValue = coeff.getValue(ip);
          // A nested coefficient can change this same-type thread-local scratch.
          // Restore active dimensions after the callback, before any entry access.
          for (size_t a = 0; a < trialScalarCount; ++a)
            GtrS[a].resize(static_cast<std::uint8_t>(spaceDimension));
          for (size_t b = 0; !sharedGradients && b < testScalarCount; ++b)
            GteS[b].resize(static_cast<std::uint8_t>(spaceDimension));
          const auto& trTabS = trialScalarElement.getTabulation(*m_qf);
          const auto& teTabS = testScalarElement.getTabulation(*m_qf);

          // Map scalar reference gradients to physical gradients
          if (spaceDimension != d)
          {
            for (size_t a = 0; a < trialScalarCount; ++a)
            {
              const auto g = trTabS.getGradient(qp, a);
              GtrS[a].setZero();
              for (size_t axis = 0; axis < spaceDimension; ++axis)
              {
                for (size_t ref = 0; ref < d; ++ref)
                  GtrS[a][axis] += Jinv(ref, axis) * g[ref];
              }
            }
            for (size_t b = 0; !sharedGradients && b < testScalarCount; ++b)
            {
              const auto g = teTabS.getGradient(qp, b);
              GteS[b].setZero();
              for (size_t axis = 0; axis < spaceDimension; ++axis)
              {
                for (size_t ref = 0; ref < d; ++ref)
                  GteS[b][axis] += Jinv(ref, axis) * g[ref];
              }
            }
          }
          else if (d == 3)
          {
            const ScalarType a00 = Jinv(0, 0), a10 = Jinv(1, 0), a20 = Jinv(2, 0);
            const ScalarType a01 = Jinv(0, 1), a11 = Jinv(1, 1), a21 = Jinv(2, 1);
            const ScalarType a02 = Jinv(0, 2), a12 = Jinv(1, 2), a22 = Jinv(2, 2);

            for (size_t a = 0; a < trialScalarCount; ++a)
            {
              const auto g = trTabS.getGradient(qp, a);
              const ScalarType gx = g[0], gy = g[1], gz = g[2];
              GtrS[a][0] = a00 * gx + a10 * gy + a20 * gz;
              GtrS[a][1] = a01 * gx + a11 * gy + a21 * gz;
              GtrS[a][2] = a02 * gx + a12 * gy + a22 * gz;
            }
            for (size_t b = 0; !sharedGradients && b < testScalarCount; ++b)
            {
              const auto g = teTabS.getGradient(qp, b);
              const ScalarType gx = g[0], gy = g[1], gz = g[2];
              GteS[b][0] = a00 * gx + a10 * gy + a20 * gz;
              GteS[b][1] = a01 * gx + a11 * gy + a21 * gz;
              GteS[b][2] = a02 * gx + a12 * gy + a22 * gz;
            }
          }
          else if (d == 2)
          {
            const ScalarType a00 = Jinv(0, 0), a10 = Jinv(1, 0);
            const ScalarType a01 = Jinv(0, 1), a11 = Jinv(1, 1);

            for (size_t a = 0; a < trialScalarCount; ++a)
            {
              const auto g = trTabS.getGradient(qp, a);
              const ScalarType gx = g[0], gy = g[1];
              GtrS[a][0] = a00 * gx + a10 * gy;
              GtrS[a][1] = a01 * gx + a11 * gy;
            }
            for (size_t b = 0; !sharedGradients && b < testScalarCount; ++b)
            {
              const auto g = teTabS.getGradient(qp, b);
              const ScalarType gx = g[0], gy = g[1];
              GteS[b][0] = a00 * gx + a10 * gy;
              GteS[b][1] = a01 * gx + a11 * gy;
            }
          }
          else if (d == 1)
          {
            const ScalarType a00 = Jinv(0, 0);
            for (size_t a = 0; a < trialScalarCount; ++a)
            {
              const auto g = trTabS.getGradient(qp, a);
              GtrS[a][0] = a00 * g[0];
            }
            for (size_t b = 0; !sharedGradients && b < testScalarCount; ++b)
            {
              const auto g = teTabS.getGradient(qp, b);
              GteS[b][0] = a00 * g[0];
            }
          }
          else
          {
            assert(false);
          }

          if constexpr (std::is_same_v<CoefficientRangeType, ScalarType> ||
            std::is_same_v<CoefficientRangeType, Real>)
          {
            const ScalarType csv = coefficientValue;
            const bool finalize = qp + 1 == q.getSize();

            // Block-diagonal assembly with scalar coefficient
            if (symmetric)
            {
              for (size_t ib = 0; ib < testScalarCount; ++ib)
              {
                const auto& gb = GteS[ib];
                for (size_t ia = 0; ia <= ib; ++ia)
                {
                  const ScalarType kij = wdet * csv * Math::dot(gb, GtrS[ia]);
                  if (kij == ScalarType(0) && (!finalize || qp == 0))
                    continue;
                  const size_t first = ib * vdim * ntr + ia * vdim;
                  if (kij != ScalarType(0))
                    A_data[first] += kij;
                  if (finalize)
                  {
                    const ScalarType value = A_data[first];
                    size_t entry = first + ntr + 1;
                    for (size_t c = 1; c < vdim; ++c)
                    {
                      A_data[entry] = value;
                      entry += ntr + 1;
                    }
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
                  if (kij == ScalarType(0) && (!finalize || qp == 0))
                    continue;
                  const size_t first = ib * vdim * ntr + ia * vdim;
                  if (kij != ScalarType(0))
                    A_data[first] += kij;
                  if (finalize)
                  {
                    const ScalarType value = A_data[first];
                    size_t entry = first + ntr + 1;
                    for (size_t c = 1; c < vdim; ++c)
                    {
                      A_data[entry] = value;
                      entry += ntr + 1;
                    }
                  }
                }
              }
            }
          }
          else if constexpr (FormLanguage::IsMatrixRange<CoefficientRangeType>::Value)
          {
            m_cmv = coefficientValue;

            // With matrix coefficient, the Frobenius inner product becomes:
            // (A * J_trial) : J_test = sum_{c,j} A(c,c2) * grad_trial_c2[j] * grad_test_c[j]
            // This breaks block-diagonal structure.
            for (size_t ib = 0; ib < testScalarCount; ++ib)
            {
              const auto& gb = GteS[ib];
              for (size_t ia = 0; ia < trialScalarCount; ++ia)
              {
                // grad_tr . grad_te (scalar dot of physical gradients)
                const ScalarType gradDot = Math::dot(gb, GtrS[ia]);

                for (size_t dComp = 0; dComp < vdim; ++dComp)
                {
                  const size_t row = ib * vdim + dComp;
                  for (size_t cComp = 0; cComp < vdim; ++cComp)
                  {
                    const size_t col = ia * vdim + cComp;
                    A_data[row * ntr + col] += wdet * m_cmv(dComp, cComp) * gradDot;
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

      /**
       * @brief Returns an entry of the element matrix.
       * @param tr Trial shape-function expression.
       * @param te Test shape-function expression.
       * @returns Integral computed by the quadrature rule.
       */
      inline ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_mat(te, tr);
      }

      /**
       * @brief Returns the integration region.
       * @returns The integration region.
       */
      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      std::unique_ptr<Geometry::PolytopeQuadrature> m_quadrature;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      Math::SpatialMatrix<ScalarType> m_cmv;
      Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_mat;
  };

  /**
   * @brief Deduction guide for @c QuadratureRule.
   * @param integrand Expression to integrate.
   */
  template <size_t KTrial, size_t KTest, class CoefficientDerived, class LHSDerived,
    class RHSDerived, class Scalar, class Mesh>
  QuadratureRule(const Dot<
    ShapeFunctionBase<Mult<FunctionBase<CoefficientDerived>,
                        ShapeFunctionBase<Jacobian<ShapeFunction<LHSDerived,
                                            H1<KTrial, Scalar, Mesh>, TrialSpace>>,
                          H1<KTrial, Scalar, Mesh>, TrialSpace>>,
      H1<KTrial, Scalar, Mesh>, TrialSpace>,
    ShapeFunctionBase<
      Jacobian<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
      H1<KTest, Scalar, Mesh>, TestSpace>>& integrand)
    -> QuadratureRule<
      Dot<ShapeFunctionBase<Mult<FunctionBase<CoefficientDerived>,
                              ShapeFunctionBase<Jacobian<ShapeFunction<LHSDerived,
                                                  H1<KTrial, Scalar, Mesh>, TrialSpace>>,
                                H1<KTrial, Scalar, Mesh>, TrialSpace>>,
            H1<KTrial, Scalar, Mesh>, TrialSpace>,
        ShapeFunctionBase<
          Jacobian<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
          H1<KTest, Scalar, Mesh>, TestSpace>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Specialization for @f$\int \nabla u \cdot \nabla v \ dx@f$ with H1 trial and test shape functions.
   *
   * ## Physical-gradient preparation
   *
   * At each quadrature point, reference gradients are transformed by the
   * current inverse Jacobian. If both operands borrow the same reference
   * table object, their physical gradients are identical and share one
   * scratch buffer. Distinct tables retain independent buffers. This identity
   * is tested only while both tables are borrowed by the current bind;
   * physical gradients are recomputed on every cell, and no table pointer
   * is stored beyond the bind.
   * Transform expressions and local matrix summation order are unchanged.
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
  template <size_t KTrial, size_t KTest, class LHSDerived, class RHSDerived, class Scalar,
    class Mesh>
  class QuadratureRule<Dot<ShapeFunctionBase<Grad<ShapeFunction<LHSDerived,
                                               H1<KTrial, Scalar, Mesh>, TrialSpace>>,
                             H1<KTrial, Scalar, Mesh>, TrialSpace>,
    ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
      H1<KTest, Scalar, Mesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<typename FormLanguage::Traits<
        Dot<ShapeFunctionBase<
              Grad<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
              H1<KTrial, Scalar, Mesh>, TrialSpace>,
          ShapeFunctionBase<
            Grad<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
            H1<KTest, Scalar, Mesh>, TestSpace>>>::ScalarType>
  {
    public:
      /// @brief Reports this handler as an optimized specialization.
      static constexpr bool Specialized = true;

      /// @brief Trial finite element space type.
      using TrialFESType = H1<KTrial, Scalar, Mesh>;
      /// @brief Test finite element space type.
      using TestFESType = H1<KTest, Scalar, Mesh>;

      /// @brief Left-hand side operand type.
      using LHSType =
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, TrialFESType, TrialSpace>>,
          TrialFESType, TrialSpace>;

      /// @brief Right-hand side operand type.
      using RHSType =
        ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, TestFESType, TestSpace>>,
          TestFESType, TestSpace>;

      /// @brief Integrand expression type.
      using IntegrandType = Dot<LHSType, RHSType>;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;

      /// @brief Parent class type.
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      /**
       * @brief Constructs the integrator for the given integrand.
       * @param integrand Expression to integrate.
       */
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

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
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

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_quadrature(std::move(other.m_quadrature)),
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

      /**
       * @brief Gets the integrand.
       * @returns The integrand.
       */
      constexpr const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      /**
       * @brief Returns the polytope the integrator is bound to.
       * @returns The polytope the integrator is bound to.
       */
      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      /**
       * @brief Binds the integrator to a polytope and tabulates the quadrature on it.
       * @param polytope Mesh entity used by this operation.
       * @returns Reference to this object after the operation.
       */
      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();

        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe = testfes.getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe.getOrder();
        const size_t order =
          this->getOrder(polytope).value_or(integrand.getOrder(polytope).value_or(
            (k_tr == 0 || k_te == 0) ? 0 : (k_tr + k_te - 2)));

        const auto geometry = polytope.getGeometry();
        const bool recomputeQf = (!m_set || m_order != order || m_geometry != geometry);

        if (recomputeQf)
        {
          m_set = true;
          m_order = order;
          m_geometry = geometry;

          m_qf = &QF::PolytopeQuadratureFormula::get(order, geometry);
        }

        assert(m_qf);
        m_quadrature = std::make_unique<Geometry::PolytopeQuadrature>(polytope, *m_qf);

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

        // A point has no reference derivatives.
        if (d == 0)
          return *this;
        const size_t spaceDimension = polytope.getMesh().getSpaceDimension();

        const auto& scalarTrial = [&]() -> const auto& {
          if constexpr (FormLanguage::IsMatrixRange<Scalar>::Value)
            return trialfe.getScalarElement();
          else
            return trialfe;
        }();
        const auto& scalarTest = [&]() -> const auto& {
          if constexpr (FormLanguage::IsMatrixRange<Scalar>::Value)
            return testfe.getScalarElement();
          else
            return testfe;
        }();
        const size_t scalarTrialCount = scalarTrial.getCount();
        const size_t scalarTestCount = scalarTest.getCount();
        const size_t components = ntr / scalarTrialCount;
        assert(components == nte / scalarTestCount);
        if constexpr (FormLanguage::IsMatrixRange<Scalar>::Value)
        {
          assert(trialfes.getRows() == testfes.getRows());
          assert(trialfes.getColumns() == testfes.getColumns());
        }
        const auto& trTab = scalarTrial.getTabulation(*m_qf);
        const auto& teTab = scalarTest.getTabulation(*m_qf);
        const bool sharedGradients =
          static_cast<const void*>(&trTab) == static_cast<const void*>(&teTab);

        // Scratch for physical gradients at this qp.
        static thread_local std::vector<Math::SpatialVector<ScalarType>> Gtr;
        static thread_local std::vector<Math::SpatialVector<ScalarType>> GteStorage;
        auto& Gte = sharedGradients ? Gtr : GteStorage;

        if (Gtr.size() < scalarTrialCount)
          Gtr.resize(scalarTrialCount);
        if (Gte.size() < scalarTestCount)
          Gte.resize(scalarTestCount);
        for (size_t a = 0; a < scalarTrialCount; ++a)
          Gtr[a].resize(static_cast<std::uint8_t>(spaceDimension));
        for (size_t b = 0; b < scalarTestCount; ++b)
          Gte[b].resize(static_cast<std::uint8_t>(spaceDimension));

        assert(m_quadrature);
        const auto& q = *m_quadrature;
        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto Jinv = p.getJacobianInverse();

          // Embedded cells have one physical derivative per ambient axis.
          if (spaceDimension != d)
          {
            for (size_t a = 0; a < scalarTrialCount; ++a)
            {
              const auto g = trTab.getGradient(qp, a);
              Gtr[a].setZero();
              for (size_t k = 0; k < spaceDimension; ++k)
              {
                for (size_t l = 0; l < d; ++l)
                  Gtr[a][k] += Jinv(l, k) * g[l];
              }
            }
            for (size_t b = 0; !sharedGradients && b < scalarTestCount; ++b)
            {
              const auto g = teTab.getGradient(qp, b);
              Gte[b].setZero();
              for (size_t k = 0; k < spaceDimension; ++k)
              {
                for (size_t l = 0; l < d; ++l)
                  Gte[b][k] += Jinv(l, k) * g[l];
              }
            }
          }
          else if (d == 3)
          {
            const ScalarType a00 = Jinv(0, 0), a10 = Jinv(1, 0), a20 = Jinv(2, 0);
            const ScalarType a01 = Jinv(0, 1), a11 = Jinv(1, 1), a21 = Jinv(2, 1);
            const ScalarType a02 = Jinv(0, 2), a12 = Jinv(1, 2), a22 = Jinv(2, 2);

            for (size_t a = 0; a < scalarTrialCount; ++a)
            {
              const auto g = trTab.getGradient(qp, a);
              const ScalarType gx = g[0], gy = g[1], gz = g[2];
              Gtr[a][0] = a00 * gx + a10 * gy + a20 * gz;
              Gtr[a][1] = a01 * gx + a11 * gy + a21 * gz;
              Gtr[a][2] = a02 * gx + a12 * gy + a22 * gz;
            }
            for (size_t b = 0; !sharedGradients && b < scalarTestCount; ++b)
            {
              const auto g = teTab.getGradient(qp, b);
              const ScalarType gx = g[0], gy = g[1], gz = g[2];
              Gte[b][0] = a00 * gx + a10 * gy + a20 * gz;
              Gte[b][1] = a01 * gx + a11 * gy + a21 * gz;
              Gte[b][2] = a02 * gx + a12 * gy + a22 * gz;
            }
          }
          else if (d == 2)
          {
            const ScalarType a00 = Jinv(0, 0), a10 = Jinv(1, 0);
            const ScalarType a01 = Jinv(0, 1), a11 = Jinv(1, 1);

            for (size_t a = 0; a < scalarTrialCount; ++a)
            {
              const auto g = trTab.getGradient(qp, a);
              const ScalarType gx = g[0], gy = g[1];
              Gtr[a][0] = a00 * gx + a10 * gy;
              Gtr[a][1] = a01 * gx + a11 * gy;
            }
            for (size_t b = 0; !sharedGradients && b < scalarTestCount; ++b)
            {
              const auto g = teTab.getGradient(qp, b);
              const ScalarType gx = g[0], gy = g[1];
              Gte[b][0] = a00 * gx + a10 * gy;
              Gte[b][1] = a01 * gx + a11 * gy;
            }
          }
          else if (d == 1)
          {
            const ScalarType a00 = Jinv(0, 0);
            for (size_t a = 0; a < scalarTrialCount; ++a)
            {
              const auto g = trTab.getGradient(qp, a);
              Gtr[a][0] = a00 * g[0];
            }
            for (size_t b = 0; !sharedGradients && b < scalarTestCount; ++b)
            {
              const auto g = teTab.getGradient(qp, b);
              Gte[b][0] = a00 * g[0];
            }
          }
          else
          {
            assert(false);
          }

          if constexpr (FormLanguage::IsMatrixRange<Scalar>::Value)
          {
            for (size_t b = 0; b < scalarTestCount; ++b)
            {
              for (size_t a = 0; a < (symmetric ? b + 1 : scalarTrialCount); ++a)
              {
                const ScalarType entry = wdet * Math::dot(Gtr[a], Gte[b]);
                for (size_t c = 0; c < components; ++c)
                  A[(b * components + c) * ntr + a * components + c] += entry;
              }
            }
          }
          else if (symmetric)
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

      /**
       * @brief Returns an entry of the element matrix.
       * @param tr Trial shape-function expression.
       * @param te Test shape-function expression.
       * @returns Integral computed by the quadrature rule.
       */
      inline ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_mat(te, tr);
      }

      /**
       * @brief Returns the integration region.
       * @returns The integration region.
       */
      virtual Geometry::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      std::unique_ptr<Geometry::PolytopeQuadrature> m_quadrature;

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

  /**
   * @brief Deduction guide for @c QuadratureRule.
   * @param integrand Expression to integrate.
   */
  template <size_t KTrial, size_t KTest, class LHSDerived, class RHSDerived, class Scalar,
    class Mesh>
  QuadratureRule(const Dot<ShapeFunctionBase<Grad<ShapeFunction<LHSDerived,
                                               H1<KTrial, Scalar, Mesh>, TrialSpace>>,
                             H1<KTrial, Scalar, Mesh>, TrialSpace>,
    ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
      H1<KTest, Scalar, Mesh>, TestSpace>>& integrand)
    -> QuadratureRule<
      Dot<ShapeFunctionBase<
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
   * ## Physical-gradient preparation
   *
   * The scalar basis gradients are mapped at the current physical sample by
   * @f$G_{aj}=\sum_{r=0}^{d-1}\partial_r\widehat\phi_a\,(J_K^{\dagger})_{rj}@f$.
   * If both operands borrow the same reference-table object, these images
   * share one scratch buffer. Distinct tables retain separate buffers even
   * when their sizes agree. Identity is checked during the current borrow;
   * no table pointer is retained and physical gradients are never reused
   * across binds. Each cell and quadrature point recomputes its physical
   * gradients. Mapping expressions and local-entry summation order remain
   * unchanged.
   *
   * ## Component-block contraction
   *
   * Let @f$\phi_a@f$ and @f$\psi_b@f$ denote physical trial and test scalar
   * bases, and let @f$m@f$ count scalar components of the value range. The
   * component basis @f$\Phi_{a,c}=\phi_a E_c@f$ has local index @f$am+c@f$.
   * The underlying scalar bases are real-valued in both real and complex H1
   * spaces. Their mapped gradients define the scalar block
   * @f[
   * S_{ba}=\sum_q w_q\,\mathrm{distortion}_K(\widehat x_q)
   *   \nabla\phi_a(x_q)\cdot\nabla\psi_b(x_q),
   * \qquad A_{bm+c,am+e}=\delta_{ce}S_{ba}.
   * @f]
   * Here @f$x_q=F_K(\widehat x_q)@f$ is the mapped reference sample and
   * @f$w_q@f$ is its reference quadrature weight.
   * Only the first component block is accumulated over quadrature points.
   * Its entries are copied to the remaining component-diagonal blocks after
   * the sweep; cross-component entries retain their initialized zeros.
   * This preserves each entry's quadrature summation order without allocating
   * a separate scalar matrix. Matrix-valued ranges use their flattened scalar
   * components. Component-coupling weighted forms have a separate contract.
   *
   * Zero component ranges produce a @f$0\times0@f$ matrix before scalar
   * tabulation or contraction, including for unequal polynomial degrees.
   * On a point, the physical derivative is zero and no scalar gradient
   * contraction is performed, independently of the ambient dimension.
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
  template <size_t KTrial, size_t KTest, class LHSDerived, class RHSDerived, class Scalar,
    class Mesh>
  class QuadratureRule<
    Dot<ShapeFunctionBase<
          Jacobian<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
          H1<KTrial, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        Jacobian<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
        H1<KTest, Scalar, Mesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<typename FormLanguage::Traits<
        Dot<ShapeFunctionBase<
              Jacobian<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
              H1<KTrial, Scalar, Mesh>, TrialSpace>,
          ShapeFunctionBase<
            Jacobian<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
            H1<KTest, Scalar, Mesh>, TestSpace>>>::ScalarType>
  {
    public:
      /// @brief Reports this handler as an optimized specialization.
      static constexpr bool Specialized = true;

      /// @brief Trial finite element space type.
      using TrialFESType = H1<KTrial, Scalar, Mesh>;
      /// @brief Test finite element space type.
      using TestFESType = H1<KTest, Scalar, Mesh>;

      /// @brief Left-hand side operand type.
      using LHSType =
        ShapeFunctionBase<Jacobian<ShapeFunction<LHSDerived, TrialFESType, TrialSpace>>,
          TrialFESType, TrialSpace>;

      /// @brief Right-hand side operand type.
      using RHSType =
        ShapeFunctionBase<Jacobian<ShapeFunction<RHSDerived, TestFESType, TestSpace>>,
          TestFESType, TestSpace>;

      /// @brief Integrand expression type.
      using IntegrandType = Dot<LHSType, RHSType>;
      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;
      /// @brief Parent class type.
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      /**
       * @brief Constructs the integrator for the given integrand.
       * @param integrand Expression to integrate.
       */
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_quadrature(std::move(other.m_quadrature)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_mat(std::move(other.m_mat))
      {}

      /**
       * @brief Returns the polytope the integrator is bound to.
       * @returns The polytope the integrator is bound to.
       */
      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      /**
       * @brief Binds the integrator to a polytope and tabulates the quadrature on it.
       * @param polytope Mesh entity used by this operation.
       * @returns Reference to this object after the operation.
       */
      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();

        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe = testfes.getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe.getOrder();
        const size_t order =
          this->getOrder(polytope).value_or(integrand.getOrder(polytope).value_or(
            ((k_tr == 0 || k_te == 0) ? 0 : (k_tr + k_te - 2))));

        const auto geometry = polytope.getGeometry();
        const bool recomputeQf = (!m_set || m_order != order || m_geometry != geometry);

        if (recomputeQf)
        {
          m_set = true;
          m_order = order;
          m_geometry = geometry;

          m_qf = &QF::PolytopeQuadratureFormula::get(order, geometry);
        }

        assert(m_qf);
        m_quadrature = std::make_unique<Geometry::PolytopeQuadrature>(polytope, *m_qf);

        const size_t ntr = lhs.getDOFs(polytope);
        const size_t nte = rhs.getDOFs(polytope);

        // --- infer vdim from H1Element<K, Scalar> conventions ---
        // Each element's scalar count depends on its polynomial order.
        const H1Element<KTrial, ScalarType> trialScalarElement(polytope.getGeometry());
        const H1Element<KTest, ScalarType> testScalarElement(polytope.getGeometry());
        const size_t trialScalarCount = trialScalarElement.getCount();
        const size_t testScalarCount = testScalarElement.getCount();

        assert(trialScalarCount > 0);
        assert(testScalarCount > 0);
        assert(ntr % trialScalarCount == 0);
        assert(nte % testScalarCount == 0);

        const size_t vdimTr = ntr / trialScalarCount;
        assert(vdimTr == nte / testScalarCount);
        const size_t vdim = vdimTr;

        const bool symmetric =
          (&trialfes.getMesh() == &testfes.getMesh()) && (ntr == nte);

        m_mat.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_mat.setZero();
        ScalarType* A = m_mat.data(); // row-major (rows=test, cols=trial)

        if (d == 0 || vdim == 0)
          return *this;
        const size_t spaceDimension = polytope.getMesh().getSpaceDimension();

        // Use scalar tabulations (fast, cached in H1Element<K, Scalar>::getTabulation)
        const auto& trTabS = trialScalarElement.getTabulation(*m_qf);
        const auto& teTabS = testScalarElement.getTabulation(*m_qf);
        const bool sharedGradients =
          static_cast<const void*>(&trTabS) == static_cast<const void*>(&teTabS);

        static thread_local std::vector<Math::SpatialVector<ScalarType>> GtrS;
        static thread_local std::vector<Math::SpatialVector<ScalarType>> GteStorage;
        auto& GteS = sharedGradients ? GtrS : GteStorage;

        if (GtrS.size() < trialScalarCount)
          GtrS.resize(trialScalarCount);
        if (GteS.size() < testScalarCount)
          GteS.resize(testScalarCount);
        for (size_t a = 0; a < trialScalarCount; ++a)
          GtrS[a].resize(static_cast<std::uint8_t>(spaceDimension));
        for (size_t b = 0; b < testScalarCount; ++b)
          GteS[b].resize(static_cast<std::uint8_t>(spaceDimension));

        assert(m_quadrature);
        const auto& q = *m_quadrature;
        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto Jinv = p.getJacobianInverse();

          // Physical gradients have an entry for every ambient coordinate,
          // including on embedded cells where the inverse Jacobian is rectangular.
          for (size_t axis = 0; axis < spaceDimension; ++axis)
          {
            for (size_t a = 0; a < trialScalarCount; ++a)
            {
              ScalarType value = 0;
              for (size_t r = 0; r < d; ++r)
                value += trTabS.getGradient(qp, a)[r] * Jinv(r, axis);
              GtrS[a][axis] = value;
            }
            for (size_t b = 0; !sharedGradients && b < testScalarCount; ++b)
            {
              ScalarType value = 0;
              for (size_t r = 0; r < d; ++r)
                value += teTabS.getGradient(qp, b)[r] * Jinv(r, axis);
              GteS[b][axis] = value;
            }
          }

          if (symmetric)
          {
            for (size_t ib = 0; ib < testScalarCount; ++ib)
            {
              const auto& gb = GteS[ib];

              for (size_t ia = 0; ia <= ib; ++ia)
              {
                const ScalarType kij = wdet * Math::dot(gb, GtrS[ia]);
                if (kij == ScalarType(0))
                  continue;

                const size_t row = ib * vdim;
                const size_t col = ia * vdim;
                A[row * ntr + col] += kij;
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

                const size_t row = ib * vdim;
                const size_t col = ia * vdim;
                A[row * ntr + col] += kij;
              }
            }
          }
        }

        for (size_t c = 1; c < vdim; ++c)
        {
          for (size_t ib = 0; ib < testScalarCount; ++ib)
          {
            for (size_t ia = 0; ia < trialScalarCount; ++ia)
            {
              A[(ib * vdim + c) * ntr + ia * vdim + c] = A[(ib * vdim) * ntr + ia * vdim];
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

      /**
       * @brief Returns an entry of the element matrix.
       * @param tr Trial shape-function expression.
       * @param te Test shape-function expression.
       * @returns Integral computed by the quadrature rule.
       */
      inline ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_mat(te, tr);
      }

      /**
       * @brief Returns the integration region.
       * @returns The integration region.
       */
      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      std::unique_ptr<Geometry::PolytopeQuadrature> m_quadrature;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_mat;
  };

  /**
   * @brief Deduction guide for @c QuadratureRule.
   * @param integrand Expression to integrate.
   */
  template <size_t KTrial, size_t KTest, class LHSDerived, class RHSDerived, class Scalar,
    class Mesh>
  QuadratureRule(
    const Dot<ShapeFunctionBase<
                Jacobian<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
                H1<KTrial, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        Jacobian<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
        H1<KTest, Scalar, Mesh>, TestSpace>>& integrand)
    -> QuadratureRule<
      Dot<ShapeFunctionBase<
            Jacobian<ShapeFunction<LHSDerived, H1<KTrial, Scalar, Mesh>, TrialSpace>>,
            H1<KTrial, Scalar, Mesh>, TrialSpace>,
        ShapeFunctionBase<
          Jacobian<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>>,
          H1<KTest, Scalar, Mesh>, TestSpace>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Integration of a rank-one form in the H1 Jacobian of ShapeFunction.
   *
   * Assembles
   * @f[
   *   \int_K (A : \mathbf{J}\,u)\,\overline{(B : \mathbf{J}\,v)}\,dx.
   * @f]
   * Each coefficient is evaluated once per quadrature point. The
   * reference gradients of the scalar H1 basis are obtained from the shared
   * element tabulation; unlike the P1 case, they are evaluated at every
   * quadrature point because they vary for polynomial degrees greater than one.
   */
  template <size_t KTrial, size_t KTest, class CoeffDerived, class LHSDerived,
    class RHSDerived, class Scalar, class Mesh>
    requires(!FormLanguage::IsMatrixRange<Scalar>::Value)
  class QuadratureRule<
    Dot<ShapeFunctionBase<Dot<FunctionBase<CoeffDerived>,
                            ShapeFunctionBase<Jacobian<ShapeFunction<LHSDerived,
                                                H1<KTrial, Scalar, Mesh>, TrialSpace>>,
                              H1<KTrial, Scalar, Mesh>, TrialSpace>>,
          H1<KTrial, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<Dot<FunctionBase<CoeffDerived>,
                          ShapeFunctionBase<Jacobian<ShapeFunction<RHSDerived,
                                              H1<KTest, Scalar, Mesh>, TestSpace>>,
                            H1<KTest, Scalar, Mesh>, TestSpace>>,
        H1<KTest, Scalar, Mesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<FunctionBase<CoeffDerived>>::ScalarType>
  {
    public:
      /// @brief Trial finite element space type.
      using TrialFESType = H1<KTrial, Scalar, Mesh>;
      /// @brief Test finite element space type.
      using TestFESType = H1<KTest, Scalar, Mesh>;
      /// @brief Scalar coefficient expression type.
      using CoefficientType = FunctionBase<CoeffDerived>;
      /// @brief Scalar type of the coefficient expression.
      using ScalarType = typename FormLanguage::Traits<CoefficientType>::ScalarType;

      /// @brief Left-hand shape-function expression type.
      using LHSType = ShapeFunctionBase<
        Dot<CoefficientType,
          ShapeFunctionBase<Jacobian<ShapeFunction<LHSDerived, TrialFESType, TrialSpace>>,
            TrialFESType, TrialSpace>>,
        TrialFESType, TrialSpace>;

      /// @brief Right-hand shape-function expression type.
      using RHSType = ShapeFunctionBase<
        Dot<CoefficientType,
          ShapeFunctionBase<Jacobian<ShapeFunction<RHSDerived, TestFESType, TestSpace>>,
            TestFESType, TestSpace>>,
        TestFESType, TestSpace>;

      /// @brief Rank-one Jacobian integrand type.
      using IntegrandType = Dot<LHSType, RHSType>;
      /// @brief Base local bilinear-form integrator type.
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      /// @brief Marks this specialization for compile-time selection tests.
      static constexpr bool IsRankOneJacobianForm = true;

      /**
       * @brief Constructs an integrator for @p integrand.
       * @param integrand Expression to integrate.
       */
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr),
          m_quadrature(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Copy-constructs the integrator from @p other.
       * @param other Object to copy from.
       */
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr),
          m_quadrature(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Move-constructs the integrator from @p other.
       * @param other Object to move from.
       */
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_quadrature(std::move(other.m_quadrature)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_trialAction(std::move(other.m_trialAction)),
          m_testAction(std::move(other.m_testAction)),
          m_matrix(std::move(other.m_matrix))
      {}

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& integrand = *m_integrand;
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialCoefficient = lhs.getDerived().getLHS();
        const auto& testCoefficient = rhs.getDerived().getLHS();
        const bool sharedCoefficient =
          trialCoefficient.getUUID() == testCoefficient.getUUID();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();
        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe = testfes.getFiniteElement(d, idx);

        const size_t kTrial = trialfe.getOrder();
        const size_t kTest = testfe.getOrder();
        const size_t order = this->getOrder(polytope).value_or(
          (kTrial == 0 || kTest == 0) ? 0 : (kTrial + kTest - 2));

        const auto geometry = polytope.getGeometry();
        if (!m_set || m_order != order || m_geometry != geometry)
        {
          m_set = true;
          m_order = order;
          m_geometry = geometry;
          m_qf = &QF::PolytopeQuadratureFormula::get(order, geometry);
        }

        assert(m_qf);
        m_quadrature = std::make_unique<Geometry::PolytopeQuadrature>(polytope, *m_qf);

        const H1Element<KTrial, ScalarType> trialScalarElement(geometry);
        const H1Element<KTest, ScalarType> testScalarElement(geometry);
        const size_t trialScalarCount = trialScalarElement.getCount();
        const size_t testScalarCount = testScalarElement.getCount();
        const size_t trialVdim = trialfes.getVectorDimension();
        const size_t testVdim = testfes.getVectorDimension();

        assert(trialVdim == testVdim);
        assert(trialfe.getCount() == trialScalarCount * trialVdim);
        assert(testfe.getCount() == testScalarCount * testVdim);

        const size_t nTrial = trialfe.getCount();
        const size_t nTest = testfe.getCount();
        m_trialAction.resize(nTrial);
        m_testAction.resize(nTest);
        m_matrix.resize(
          static_cast<Eigen::Index>(nTest), static_cast<Eigen::Index>(nTrial));
        m_matrix.setZero();

        assert(m_quadrature);
        const auto& quadrature = *m_quadrature;

        for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
        {
          const auto& point = quadrature.getPoint(qp);
          const IntegrationPoint ip(quadrature, qp);
          const ScalarType weight =
            static_cast<ScalarType>(m_qf->getWeight(qp) * point.getDistortion());
          const auto Jinv = point.getJacobianInverse();
          const auto A = trialCoefficient.getValue(ip);
          const auto& trialTabulation = trialScalarElement.getTabulation(*m_qf);

          for (size_t node = 0; node < trialScalarCount; ++node)
          {
            const auto gradient = trialTabulation.getGradient(qp, node);
            for (size_t component = 0; component < trialVdim; ++component)
            {
              const size_t local = node * trialVdim + component;
              ScalarType action = 0;
              for (size_t j = 0; j < d; ++j)
              {
                ScalarType physicalDerivative = 0;
                for (size_t r = 0; r < d; ++r)
                  physicalDerivative += gradient[r] * Jinv(r, j);
                action += A(component, j) * physicalDerivative;
              }
              m_trialAction[local] = action;
            }
          }

          if (trialfes == testfes && sharedCoefficient)
          {
            m_testAction = m_trialAction;
          }
          else
          {
            const auto B = testCoefficient.getValue(ip);
            const auto& testTabulation = testScalarElement.getTabulation(*m_qf);
            for (size_t node = 0; node < testScalarCount; ++node)
            {
              const auto gradient = testTabulation.getGradient(qp, node);
              for (size_t component = 0; component < testVdim; ++component)
              {
                const size_t local = node * testVdim + component;
                ScalarType action = 0;
                for (size_t j = 0; j < d; ++j)
                {
                  ScalarType physicalDerivative = 0;
                  for (size_t r = 0; r < d; ++r)
                    physicalDerivative += gradient[r] * Jinv(r, j);
                  action += B(component, j) * physicalDerivative;
                }
                m_testAction[local] = action;
              }
            }
          }

          if (trialfes == testfes && sharedCoefficient)
          {
            for (size_t test = 0; test < nTest; ++test)
            {
              const ScalarType scaledTest = weight * Math::conj(m_trialAction[test]);
              for (size_t trial = 0; trial <= test; ++trial)
                m_matrix(test, trial) += scaledTest * m_trialAction[trial];
            }
          }
          else
          {
            for (size_t test = 0; test < nTest; ++test)
            {
              const ScalarType scaledTest = weight * Math::conj(m_testAction[test]);
              for (size_t trial = 0; trial < nTrial; ++trial)
                m_matrix(test, trial) += scaledTest * m_trialAction[trial];
            }
          }
        }

        if (trialfes == testfes && sharedCoefficient)
        {
          m_matrix.template triangularView<Eigen::Upper>() =
            m_matrix.adjoint().template triangularView<Eigen::Upper>();
        }

        return *this;
      }

      ScalarType integrate(size_t trial, size_t test) final override
      {
        return m_matrix(test, trial);
      }

      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;
      const QF::QuadratureFormulaBase* m_qf;
      std::unique_ptr<Geometry::PolytopeQuadrature> m_quadrature;
      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;
      std::vector<ScalarType> m_trialAction;
      std::vector<ScalarType> m_testAction;
      Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_matrix;
  };

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
   *
   * The scalar nodal bases are real-valued over both real and complex fields.
   * Hence @f$(\nabla\phi_a)^T b=b\cdot\nabla\phi_a@f$ with Rodin's
   * test-second conjugated dot product. The coefficient remains linear;
   * conjugating @f$b@f$ would change the represented expression.
   * Physical gradients have the ambient dimension, including rectangular
   * inverse Jacobians on embedded cells. Point cells and empty component
   * ranges return the zero operator before evaluating the coefficient.
   *
   * ## Working storage
   *
   * Each sample stores @f$d_a=b^T\nabla\phi_a@f$ once per trial scalar basis,
   * rather than retaining a full gradient and contracting it again for each
   * test basis. Matrix traversal and per-entry accumulation remain unchanged.
   * The thread-local scalar buffer only grows. A nested same-type integrator
   * may grow or overwrite it during coefficient evaluation, but the outer
   * invocation then computes every active entry before use. Physical gradients
   * are local bounded values constructed after that callback, so no ambient
   * dimension metadata can be inherited from the nested operation. Scalar
   * reference tables are also acquired after the callback.
   */
  template <size_t KTrial, size_t KTest, class CoefficientDerived, class LHSDerived,
    class RHSDerived, class Scalar, class Mesh>
    requires(!FormLanguage::IsMatrixRange<Scalar>::Value)
  class QuadratureRule<
    Dot<ShapeFunctionBase<Mult<ShapeFunctionBase<Jacobian<ShapeFunction<LHSDerived,
                                                   H1<KTrial, Scalar, Mesh>, TrialSpace>>,
                                 H1<KTrial, Scalar, Mesh>, TrialSpace>,
                            FunctionBase<CoefficientDerived>>,
          H1<KTrial, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
        H1<KTest, Scalar, Mesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<typename FormLanguage::Traits<Dot<
        ShapeFunctionBase<Mult<ShapeFunctionBase<Jacobian<ShapeFunction<LHSDerived,
                                                   H1<KTrial, Scalar, Mesh>, TrialSpace>>,
                                 H1<KTrial, Scalar, Mesh>, TrialSpace>,
                            FunctionBase<CoefficientDerived>>,
          H1<KTrial, Scalar, Mesh>, TrialSpace>,
        ShapeFunctionBase<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
          H1<KTest, Scalar, Mesh>, TestSpace>>>::ScalarType>
  {
    public:
      /// @brief Reports this handler as an optimized specialization.
      static constexpr bool Specialized = true;

      /// @brief Trial finite element space type.
      using TrialFESType = H1<KTrial, Scalar, Mesh>;
      /// @brief Test finite element space type.
      using TestFESType = H1<KTest, Scalar, Mesh>;

      /// @brief Trial shape function type.
      using TrialSFType =
        ShapeFunctionBase<Jacobian<ShapeFunction<LHSDerived, TrialFESType, TrialSpace>>,
          TrialFESType, TrialSpace>;

      /// @brief Coefficient type appearing in the integrand.
      using CoefficientType = FunctionBase<CoefficientDerived>;

      /// @brief Left-hand side operand type.
      using LHSType =
        ShapeFunctionBase<Mult<TrialSFType, CoefficientType>, TrialFESType, TrialSpace>;

      /// @brief Right-hand side operand type.
      using RHSType = ShapeFunctionBase<ShapeFunction<RHSDerived, TestFESType, TestSpace>,
        TestFESType, TestSpace>;

      /// @brief Integrand expression type.
      using IntegrandType = Dot<LHSType, RHSType>;
      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;
      /// @brief Parent class type.
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      /**
       * @brief Constructs the integrator for the given integrand.
       * @param integrand Expression to integrate.
       */
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr),
          m_quadrature(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr),
          m_quadrature(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_quadrature(std::move(other.m_quadrature)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_mat(std::move(other.m_mat))
      {}

      /**
       * @brief Returns the polytope the integrator is bound to.
       * @returns The polytope the integrator is bound to.
       */
      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      /**
       * @brief Binds the integrator to a polytope and tabulates the quadrature on it.
       * @param polytope Mesh entity used by this operation.
       * @returns Reference to this object after the operation.
       */
      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();

        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();

        // The coefficient is the RHS of the Mult node
        const auto& coeff = lhs.getDerived().getRHS();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe = testfes.getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe.getOrder();
        // Gradient of trial (k_tr - 1) times value of test (k_te)
        const size_t order =
          this->getOrder(polytope).value_or(integrand.getOrder(polytope).value_or(
            ((k_tr == 0) ? k_te : (k_tr + k_te - 1))));

        const auto geometry = polytope.getGeometry();
        const bool recomputeQf = (!m_set || m_order != order || m_geometry != geometry);

        if (recomputeQf)
        {
          m_set = true;
          m_order = order;
          m_geometry = geometry;

          m_qf = &QF::PolytopeQuadratureFormula::get(order, geometry);
        }

        assert(m_qf);
        m_quadrature = std::make_unique<Geometry::PolytopeQuadrature>(polytope, *m_qf);

        const size_t ntr = lhs.getDOFs(polytope);
        const size_t nte = rhs.getDOFs(polytope);

        const H1Element<KTrial, ScalarType> trialScalarFE(geometry);
        const H1Element<KTest, ScalarType> testScalarFE(geometry);
        const size_t trialScalarCount = trialScalarFE.getCount();
        const size_t testScalarCount = testScalarFE.getCount();

        assert(trialScalarCount > 0 && ntr % trialScalarCount == 0);
        assert(testScalarCount > 0 && nte % testScalarCount == 0);
        const size_t vdim = ntr / trialScalarCount;
        assert(vdim == nte / testScalarCount);

        m_mat.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_mat.setZero();
        if (d == 0 || vdim == 0)
          return *this;
        ScalarType* A = m_mat.data();

        const size_t spaceDimension = polytope.getMesh().getSpaceDimension();
        static thread_local std::vector<ScalarType> trialDirectional;
        if (trialDirectional.size() < trialScalarCount)
          trialDirectional.resize(trialScalarCount);

        assert(m_quadrature);
        const auto& q = *m_quadrature;
        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          const IntegrationPoint ip(q, qp);
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto Jinv = p.getJacobianInverse();

          const auto fval = coeff.getValue(ip);
          const auto& trTab = trialScalarFE.getTabulation(*m_qf);
          const auto& teTab = testScalarFE.getTabulation(*m_qf);

          // Rectangular inverse Jacobians map into every ambient coordinate.
          if (d != spaceDimension)
          {
            for (size_t a = 0; a < trialScalarCount; ++a)
            {
              const auto g = trTab.getGradient(qp, a);
              Math::SpatialVector<ScalarType> gradient(
                static_cast<std::uint8_t>(spaceDimension));
              for (size_t axis = 0; axis < spaceDimension; ++axis)
              {
                ScalarType value = 0;
                for (size_t r = 0; r < d; ++r)
                  value += Jinv(r, axis) * g[r];
                gradient[axis] = value;
              }
              trialDirectional[a] = Math::dot(fval, gradient);
            }
          }
          else if (d == 3)
          {
            const ScalarType a00 = Jinv(0, 0), a10 = Jinv(1, 0), a20 = Jinv(2, 0);
            const ScalarType a01 = Jinv(0, 1), a11 = Jinv(1, 1), a21 = Jinv(2, 1);
            const ScalarType a02 = Jinv(0, 2), a12 = Jinv(1, 2), a22 = Jinv(2, 2);

            for (size_t a = 0; a < trialScalarCount; ++a)
            {
              const auto g = trTab.getGradient(qp, a);
              const ScalarType gx = g[0], gy = g[1], gz = g[2];
              Math::SpatialVector<ScalarType> gradient(
                static_cast<std::uint8_t>(spaceDimension));
              gradient[0] = a00 * gx + a10 * gy + a20 * gz;
              gradient[1] = a01 * gx + a11 * gy + a21 * gz;
              gradient[2] = a02 * gx + a12 * gy + a22 * gz;
              trialDirectional[a] = Math::dot(fval, gradient);
            }
          }
          else if (d == 2)
          {
            const ScalarType a00 = Jinv(0, 0), a10 = Jinv(1, 0);
            const ScalarType a01 = Jinv(0, 1), a11 = Jinv(1, 1);

            for (size_t a = 0; a < trialScalarCount; ++a)
            {
              const auto g = trTab.getGradient(qp, a);
              const ScalarType gx = g[0], gy = g[1];
              Math::SpatialVector<ScalarType> gradient(
                static_cast<std::uint8_t>(spaceDimension));
              gradient[0] = a00 * gx + a10 * gy;
              gradient[1] = a01 * gx + a11 * gy;
              trialDirectional[a] = Math::dot(fval, gradient);
            }
          }
          else if (d == 1)
          {
            const ScalarType a00 = Jinv(0, 0);
            for (size_t a = 0; a < trialScalarCount; ++a)
            {
              const auto g = trTab.getGradient(qp, a);
              Math::SpatialVector<ScalarType> gradient(
                static_cast<std::uint8_t>(spaceDimension));
              gradient[0] = a00 * g[0];
              trialDirectional[a] = Math::dot(fval, gradient);
            }
          }
          else
          {
            assert(false);
          }

          // Assemble: K(b*vdim+c, a*vdim+c) += wdet * (\nablaφ_a \cdot f) * ψ_b
          for (size_t ib = 0; ib < testScalarCount; ++ib)
          {
            const ScalarType phi_te = teTab.getBasis(qp, ib);
            for (size_t ia = 0; ia < trialScalarCount; ++ia)
            {
              const ScalarType gradDotF = trialDirectional[ia];
              const ScalarType kij = wdet * gradDotF * phi_te;
              if (kij == ScalarType(0))
                continue;
              size_t entry = ib * vdim * ntr + ia * vdim;
              for (size_t c = 0; c < vdim; ++c)
              {
                A[entry] += kij;
                entry += ntr + 1;
              }
            }
          }
        }

        return *this;
      }

      /**
       * @brief Returns an entry of the element matrix.
       * @param tr Trial shape-function expression.
       * @param te Test shape-function expression.
       * @returns Integral computed by the quadrature rule.
       */
      inline ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_mat(te, tr);
      }

      /**
       * @brief Returns the integration region.
       * @returns The integration region.
       */
      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      std::unique_ptr<Geometry::PolytopeQuadrature> m_quadrature;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_mat;
  };

  /**
   * @brief Deduction guide for @c QuadratureRule.
   * @param integrand Expression to integrate.
   */
  template <size_t KTrial, size_t KTest, class CoefficientDerived, class LHSDerived,
    class RHSDerived, class Scalar, class Mesh>
  QuadratureRule(const Dot<
    ShapeFunctionBase<Mult<ShapeFunctionBase<Jacobian<ShapeFunction<LHSDerived,
                                               H1<KTrial, Scalar, Mesh>, TrialSpace>>,
                             H1<KTrial, Scalar, Mesh>, TrialSpace>,
                        FunctionBase<CoefficientDerived>>,
      H1<KTrial, Scalar, Mesh>, TrialSpace>,
    ShapeFunctionBase<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
      H1<KTest, Scalar, Mesh>, TestSpace>>& integrand)
    -> QuadratureRule<Dot<
      ShapeFunctionBase<Mult<ShapeFunctionBase<Jacobian<ShapeFunction<LHSDerived,
                                                 H1<KTrial, Scalar, Mesh>, TrialSpace>>,
                               H1<KTrial, Scalar, Mesh>, TrialSpace>,
                          FunctionBase<CoefficientDerived>>,
        H1<KTrial, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<ShapeFunction<RHSDerived, H1<KTest, Scalar, Mesh>, TestSpace>,
        H1<KTest, Scalar, Mesh>, TestSpace>>>;
}

#endif
