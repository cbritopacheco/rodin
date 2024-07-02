/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_POTENTIAL_H
#define RODIN_VARIATIONAL_POTENTIAL_H

#include "Rodin/FormLanguage/Base.h"
#include "Rodin/FormLanguage/List.h"
#include "Rodin/QF/QuadratureFormula.h"

#include "ForwardDecls.h"
#include "Mult.h"
#include "Dot.h"
#include "Function.h"
#include "Integral.h"
#include "ShapeFunction.h"
#include "QuadratureRule.h"
#include "LinearFormIntegrator.h"

namespace Rodin::FormLanguage
{
  /**
   * @ingroup TraitsSpecializations
   * @brief Traits for Potential.
   */
  template <class LHS, class RHSDerived>
  struct Traits<Variational::Potential<LHS, Variational::FunctionBase<RHSDerived>>>
  {
    using NumberType = Real;

    using LHSType = LHS;

    using RHSType = Variational::FunctionBase<Variational::FunctionBase<RHSDerived>>;

    using KernelType = LHSType;

    using OperandType = RHSType;

    using RHSRangeType =
      typename FormLanguage::Traits<RHSType>::RangeType;

    using LHSRangeType =
      std::conditional_t<
      // If
      std::is_same_v<RHSRangeType, NumberType>,
      // Then
      NumberType,
      // Else
      std::conditional_t<
        // If
        std::is_same_v<RHSRangeType, Math::Vector<NumberType>>,
        // Then
        Math::Matrix<NumberType>,
        // Else
        void>>;

    using RangeType = RHSRangeType;
  };

  /**
   * @ingroup TraitsSpecializations
   * @brief Traits for Potential.
   */
  template <class LHS, class RHSDerived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Potential<
      LHS,
      Variational::ShapeFunctionBase<Variational::ShapeFunction<RHSDerived, FES, Space>>>>
  {
    using NumberType = Real;

    using FESType = FES;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;

    using LHSType = LHS;

    using RHSType = Variational::ShapeFunctionBase<Variational::ShapeFunction<RHSDerived, FES, Space>, FES, Space>;

    using KernelType = LHS;

    using OperandType = RHSType;

    using RHSRangeType =
      typename FormLanguage::Traits<RHSType>::RangeType;

    using LHSRangeType =
      std::conditional_t<
      // If
      std::is_same_v<RHSRangeType, NumberType>,
      // Then
      NumberType,
      // Else
      std::conditional_t<
        // If
        std::is_same_v<RHSRangeType, Math::Vector<NumberType>>,
        // Then
        Math::Matrix<NumberType>,
        // Else
        void>>;

    using RangeType = RHSRangeType;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup PotentialSpecializations Potential Template Specializations
   * @brief Template specializations of the Potential class.
   * @see Potential
   */

  /**
   * @ingroup PotentialSpecializations
   */
  template <class LHS, class RHSDerived>
  class Potential<LHS, FunctionBase<RHSDerived>> final
    : public FunctionBase<Potential<LHS, FunctionBase<RHSDerived>>>
  {
    public:
      using NumberType = Real;

      using LHSType = LHS;

      using KernelType = LHSType;

      using RHSType = FunctionBase<RHSDerived>;

      using OperandType = RHSType;

      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

      using LHSRangeType =
        std::conditional_t<
        // If
        std::is_same_v<RHSRangeType, NumberType>,
        // Then
        NumberType,
        // Else
        std::conditional_t<
          // If
          std::is_same_v<RHSRangeType, Math::Vector<NumberType>>,
          // Then
          Math::Matrix<NumberType>,
          // Else
          void>>;

      using Parent = FunctionBase<Potential<LHSType, RHSType>>;

      Potential(const KernelType& kernel, const OperandType& u)
        : m_kernel(kernel), m_u(u.copy())
      {}

      Potential(const Potential& other)
        : Parent(other),
          m_kernel(other.m_kernel),
          m_u(other.m_u->copy())
      {}

      Potential(Potential&& other)
        : Parent(std::move(other)),
          m_kernel(std::move(other.m_kernel)),
          m_u(std::move(other.m_u))
      {}

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        if constexpr (std::is_same_v<LHSRangeType, NumberType>)
        {
          static_assert(std::is_same_v<RHSRangeType, NumberType>);
          return { 1, 1 };
        }
        else if constexpr (std::is_same_v<LHSRangeType, Math::Matrix<NumberType>>)
        {
          static_assert(std::is_same_v<RHSRangeType, Math::Vector<NumberType>>);
          return getOperand().getRangeShape();
        }
        else
        {
          assert(false);
          return { 0, 0 };
        }
      }

      inline
      const auto& getKernel() const
      {
        return m_kernel.get();
      }

      inline
      const auto& getOperand() const
      {
        assert(m_u);
        return *m_u;
      }

      auto getValue(const Geometry::Point& p) const
      {
        const auto& kernel = getKernel();
        const auto& operand = getOperand();
        const auto& mesh = p.getPolytope().getMesh();
        if (m_qf.has_value())
        {
          if constexpr (std::is_same_v<RHSRangeType, NumberType>)
          {
            NumberType res = 0;
            for (auto it = mesh.getCell(); it; ++it)
            {
              const auto& polytope = *it;
              const auto& qf = m_qf.value()(polytope);
              const auto& trans = polytope.getTransformation();
              for (size_t i = 0; i < qf.getSize(); i++)
              {
                const Geometry::Point y(polytope, trans, std::ref(qf.getPoint(i)));
                res += qf.getWeight(i) * y.getDistortion() * kernel(p, y) * operand(y);
              }
            }
            return res;
          }
          else if constexpr (std::is_same_v<RHSRangeType, Math::Vector<NumberType>>)
          {
            Math::Vector<NumberType> res;
            getValue(res, p);
            return res;
          }
          else
          {
            assert(false);
            return void();
          }
        }
        else
        {
          if constexpr (std::is_same_v<RHSRangeType, NumberType>)
          {
            NumberType res = 0;
            for (auto it = mesh.getCell(); it; ++it)
            {
              const auto& polytope = *it;
              const QF::GenericPolytopeQuadrature qf(polytope.getGeometry());
              const auto& trans = polytope.getTransformation();
              for (size_t i = 0; i < qf.getSize(); i++)
              {
                const Geometry::Point y(polytope, trans, std::cref(qf.getPoint(i)));
                res += qf.getWeight(i) * y.getDistortion() * kernel(p, y) * operand(y);
              }
            }
            return res;
          }
          else if constexpr (std::is_same_v<RHSRangeType, Math::Vector<NumberType>>)
          {
            Math::Vector<NumberType> res;
            getValue(res, p);
            return res;
          }
          else
          {
            assert(false);
            return void();
          }
        }
      }

      void getValue(Math::Vector<NumberType>& res, const Geometry::Point& p) const
      {
        const auto& kernel = getKernel();
        const auto& operand = getOperand();
        const auto& mesh = p.getPolytope().getMesh();
        res.resize(getRangeShape().height());
        res.setZero();
        Math::Matrix<NumberType> kxy;
        if (m_qf.has_value())
        {
          for (auto it = mesh.getCell(); it; ++it)
          {
            const auto& polytope = *it;
            const auto& qf = m_qf.value()(polytope);
            const auto& trans = polytope.getTransformation();
            for (size_t i = 0; i < qf.getSize(); i++)
            {
              const Geometry::Point y(polytope, trans, std::cref(qf.getPoint(i)));
              kernel(kxy, p, y);
              res += qf.getWeight(i) * y.getDistortion() * kxy * operand(y);
            }
          }
        }
        else
        {
          for (auto it = mesh.getCell(); it; ++it)
          {
            const auto& polytope = *it;
            const QF::GenericPolytopeQuadrature qf(polytope.getGeometry());
            const auto& trans = polytope.getTransformation();
            for (size_t i = 0; i < qf.getSize(); i++)
            {
              const Geometry::Point y(polytope, trans, std::cref(qf.getPoint(i)));
              kernel(kxy, p, y);
              res += qf.getWeight(i) * y.getDistortion() * kxy * operand(y);
            }
          }
        }
      }

      Potential& setQuadratureFormula(
          const std::function<const QF::QuadratureFormulaBase&(const Geometry::Polytope&)>& qf)
      {
        m_qf.emplace(qf);
        return *this;
      }

      inline
      const auto& getQuadratureFormula() const
      {
        return m_qf;
      }

    private:
      std::reference_wrapper<const KernelType> m_kernel;
      std::unique_ptr<OperandType> m_u;
      std::optional<
        std::function<const QF::QuadratureFormulaBase&(const Geometry::Polytope&)>> m_qf;
  };

  /**
   * @brief CTAD for Potential.
   */
  template <class LHSType, class RHSDerived>
  Potential(const LHSType&, const FunctionBase<RHSDerived>&)
    -> Potential<LHSType, FunctionBase<RHSDerived>>;

  /**
   * @ingroup PotentialSpecializations
   *
   * Represents the expression:
   * @f[
   *  (K u)(x) = \sum_{\tau \in \mathcal{T}_h} (K_\tau u)(x)
   *  = \sum_{\tau \in \mathcal{T}_h} \sum^{n(\tau)}_{\ell = 1} w_{\tau, \ell}
   *  (K_\tau \phi_{\tau, \ell})(x), \quad w_{\tau, \ell} \in \mathbb{R}
   *  @f]
   *  where:
   *  @f[
   *    (K_\tau \phi_{\tau, \ell})(x) := \ \mathrm{p.v.} \int_{\tau} k(x, y) \phi_{\tau, \ell}(y) \ dy
   *  @f]
   *  are the basis functions.
   */
  template <class LHS, class RHSDerived, class FES, ShapeFunctionSpaceType SpaceType>
  class Potential<
    LHS,
    ShapeFunctionBase<ShapeFunction<RHSDerived, FES, SpaceType>, FES, SpaceType>> final
      : public FormLanguage::Base
  {
    public:
      using Parent = FormLanguage::Base;

      using FESType = FES;
      static constexpr ShapeFunctionSpaceType Space = SpaceType;

      using NumberType = Real;

      using LHSType = LHS;

      using KernelType = LHS;

      using RHSType = ShapeFunctionBase<ShapeFunction<RHSDerived, FES, SpaceType>, FES, Space>;

      using OperandType = RHSType;

      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

      using LHSRangeType =
        std::conditional_t<
        // If
        std::is_same_v<RHSRangeType, NumberType>,
        // Then
        NumberType,
        // Else
        std::conditional_t<
          // If
          std::is_same_v<RHSRangeType, Math::Vector<NumberType>>,
          // Then
          Math::Matrix<NumberType>,
          // Else
          void>>;

      Potential(const KernelType& kernel, const OperandType& u)
        : m_kernel(kernel), m_u(u)
      {}

      Potential(const Potential& other)
        : Parent(other),
          m_kernel(other.m_kernel), m_u(other.m_u)
      {}

      Potential(Potential&& other)
        : Parent(std::move(other)),
          m_kernel(std::move(other.m_kernel)), m_u(std::move(other.m_u))
      {}

      inline
      const KernelType& getKernel() const
      {
        return m_kernel;
      }

      inline
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      inline
      Integrator::Region getRegion() const
      {
        return Integrator::Region::Cells;
      }

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        if constexpr (std::is_same_v<LHSRangeType, NumberType>)
        {
          static_assert(std::is_same_v<RHSRangeType, NumberType>);
          return { 1, 1 };
        }
        else if constexpr (std::is_same_v<LHSRangeType, Math::Matrix<NumberType>>)
        {
          static_assert(std::is_same_v<RHSRangeType, Math::Vector<NumberType>>);
          return getOperand().getRangeShape();
        }
        else
        {
          assert(false);
          return { 0, 0 };
        }
      }

      inline Potential* copy() const noexcept override
      {
        return new Potential(*this);
      }

    private:
      std::reference_wrapper<KernelType> m_kernel;
      std::reference_wrapper<const OperandType> m_u;
  };

  /**
   * @brief CTAD for Potential.
   */
  template <class LHSType, class RHSDerived, class FESType, ShapeFunctionSpaceType SpaceType>
  Potential(const LHSType&, const ShapeFunctionBase<ShapeFunction<RHSDerived, FESType, SpaceType>, FESType, SpaceType>&)
    -> Potential<LHSType, ShapeFunctionBase<ShapeFunction<RHSDerived, FESType, SpaceType>, FESType, SpaceType>>;

  template <class Kernel, class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  class QuadratureRule<
    Dot<
      Potential<Kernel, ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>>,
      ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>
        : public GlobalBilinearFormIntegratorBase<Real>
  {
    public:
      using NumberType = Real;

      using KernelType = Kernel;

      using LHSType = Potential<KernelType, ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>>;

      using RHSType = ShapeFunctionBase<RHSDerived, TestFES, TestSpace>;

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

      QuadratureRule& setPolytope(const Geometry::Polytope& trp, const Geometry::Polytope& tep)
      {
        const auto& trptrans = trp.getTransformation();
        const auto& teptrans = tep.getTransformation();
        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = lhs.getOperand().getFiniteElementSpace();
        const auto& trialfe = trialfes.getFiniteElement(trp.getDimension(), trp.getIndex());
        const auto& testfes = rhs.getFiniteElementSpace();
        const auto& testfe = testfes.getFiniteElement(tep.getDimension(), tep.getIndex());
        const auto& kernel = lhs.getKernel();
        if (tep == trp)
        {
          assert(trp.getGeometry() == tep.getGeometry());
          switch (tep.getGeometry())
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
              const size_t order = std::max(trialfe.getOrder(), testfe.getOrder());
              m_qf.reset(new QF::GenericPolytopeQuadrature(order, Geometry::Polytope::Type::Segment));
              m_rx.resize(6);
              m_rz.resize(6);
              break;
            }
          }
        }
        else
        {
          m_qftr.reset(new QF::GenericPolytopeQuadrature(trialfe.getOrder(), trp.getGeometry()));
          const auto& qftr = *m_qftr;
          for (size_t i = 0; i < qftr.getSize(); i++)
            m_y.emplace_back(trp, trptrans, std::cref(qftr.getPoint(i)));

          m_qfte.reset(new QF::GenericPolytopeQuadrature(testfe.getOrder(), tep.getGeometry()));
          const auto& qfte = *m_qfte;
          for (size_t i = 0; i < qfte.getSize(); i++)
            m_x.emplace_back(tep, teptrans, std::cref(qfte.getPoint(i)));
        }

        return *this;
      }

      NumberType integrate(size_t tr, size_t te) override
      {
        const auto& trp = m_trp.value().get();
        const auto& tep = m_trp.value().get();
        const auto& trptrans = trp.getTransformation();
        const auto& teptrans = tep.getTransformation();
        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = lhs.getOperand().getFiniteElementSpace();
        const auto& trialfe = trialfes.getFiniteElement(trp.getDimension(), trp.getIndex());
        const auto& testfes = rhs.getFiniteElementSpace();
        const auto& testfe = testfes.getFiniteElement(tep.getDimension(), tep.getIndex());
        const auto& kernel = lhs.getKernel();
        NumberType res = 0;
        if (tep == trp)
        {
          assert(trp.getGeometry() == tep.getGeometry());
          switch (tep.getGeometry())
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
              assert(m_qf);
              const auto& qf = *m_qf;
              for (size_t i = 0; i < qf.getSize(); i++)
              {
                const NumberType xi = 1 - qf.getPoint(i).value();
                for (size_t j = 0; j < qf.getSize(); j++)
                {
                  const NumberType eta1 = qf.getPoint(j).value();
                  for (size_t h = 0; h < qf.getSize(); h++)
                  {
                    const NumberType eta2 = qf.getPoint(h).value();
                    m_rx[0] << xi, xi * (1 - eta1 + eta1 * eta2);
                    m_rz[1] << xi, xi * (1 - eta1 + eta1 * eta2);
                    m_rz[2] << xi * (1 - eta1 * eta2), xi * eta1 * (1 - eta2);
                    m_rx[3] << xi * (1 - eta1 * eta2), xi * eta1 * (1 - eta2);
                    m_rz[4] << xi, xi * eta1 * (1 - eta2);
                    m_rx[5] << xi, xi * eta1 * (1 - eta2);
                    const Geometry::Point x0(tep, teptrans, std::cref(m_rx[0]));
                    const Geometry::Point z1(tep, teptrans, std::cref(m_rz[1]));
                    const Geometry::Point z2(tep, teptrans, std::cref(m_rz[2]));
                    const Geometry::Point x3(tep, teptrans, std::cref(m_rx[3]));
                    const Geometry::Point z4(tep, teptrans, std::cref(m_rz[4]));
                    const Geometry::Point x5(tep, teptrans, std::cref(m_rx[5]));
                    const NumberType d = xi * xi * xi * eta1 * eta1 * eta2;
                    for (size_t k = 0; k < qf.getSize(); k++)
                    {
                      const NumberType w = qf.getWeight(i) * qf.getWeight(j) * qf.getWeight(h) * qf.getWeight(k);
                      const NumberType eta3 = qf.getPoint(k).value();
                      m_rz[0] << xi * (1 - eta1 * eta2 * eta3), xi * (1 - eta1);
                      m_rx[1] << xi * (1 - eta1 * eta2 * eta3), xi * (1 - eta1);
                      m_rx[2] << xi, xi * eta1 * (1 - eta2 + eta2 * eta3);
                      m_rz[3] << xi, xi * eta1 * (1 - eta2 + eta2 * eta3);
                      m_rx[4] << xi * (1 - eta1 * eta2 * eta3), xi * eta1 * (1 - eta2 * eta3);
                      m_rz[5] << xi * (1 - eta1 * eta2 * eta3), xi * eta1 * (1 - eta2 * eta3);
                      const Geometry::Point z0(tep, teptrans, std::cref(m_rz[0]));
                      const Geometry::Point x1(tep, teptrans, std::cref(m_rx[1]));
                      const Geometry::Point x2(tep, teptrans, std::cref(m_rx[2]));
                      const Geometry::Point z3(tep, teptrans, std::cref(m_rz[3]));
                      const Geometry::Point x4(tep, teptrans, std::cref(m_rx[4]));
                      const Geometry::Point z5(tep, teptrans, std::cref(m_rz[5]));
                      NumberType s0, s1, s2, s3, s4, s5;
                      if constexpr (std::is_same_v<LHSRangeType, NumberType>)
                      {
                        s0 = kernel(x0, z0) * x0.getDistortion() * z0.getDistortion();
                        s1 = kernel(x1, z1) * x1.getDistortion() * z1.getDistortion();
                        s2 = kernel(x2, z2) * x2.getDistortion() * z2.getDistortion();
                        s3 = kernel(x3, z3) * x3.getDistortion() * z3.getDistortion();
                        s4 = kernel(x4, z4) * x4.getDistortion() * z4.getDistortion();
                        s5 = kernel(x5, z5) * x5.getDistortion() * z5.getDistortion();
                      }
                      else if constexpr (std::is_same_v<LHSRangeType, Math::Vector<NumberType>>)
                      {
                        kernel(m_k0, x0, z0);
                        kernel(m_k1, x1, z1);
                        kernel(m_k2, x2, z2);
                        kernel(m_k3, x3, z3);
                        kernel(m_k4, x4, z4);
                        kernel(m_k5, x5, z5);
                        s0 = x0.getDistortion() * z0.getDistortion();
                        s1 = x1.getDistortion() * z1.getDistortion();
                        s2 = x2.getDistortion() * z2.getDistortion();
                        s3 = x3.getDistortion() * z3.getDistortion();
                        s4 = x4.getDistortion() * z4.getDistortion();
                        s5 = x5.getDistortion() * z5.getDistortion();
                      }
                      else
                      {
                        assert(false);
                        s0 = s1 = s2 = s3 = s4 = s5 = NAN;
                      }

                      for (size_t l = 0; l < testfe.getCount(); l++)
                      {
                        const auto& teb = testfe.getBasis(te);
                        for (size_t m = 0; m < trialfe.getCount(); m++)
                        {
                          const auto& trb = trialfe.getBasis(m);
                          if constexpr (std::is_same_v<LHSRangeType, NumberType>)
                          {
                            res += d * w * s0 * trb(m_rx[0]) * teb(m_rz[0]);
                            assert(std::isfinite(s0));
                            res += d * w * s1 * trb(m_rx[1]) * teb(m_rz[1]);
                            assert(std::isfinite(s1));
                            res += d * w * s2 * trb(m_rx[2]) * teb(m_rz[2]);
                            assert(std::isfinite(s2));
                            res += d * w * s3 * trb(m_rx[3]) * teb(m_rz[3]);
                            assert(std::isfinite(s3));
                            res += d * w * s4 * trb(m_rx[4]) * teb(m_rz[4]);
                            assert(std::isfinite(s4));
                            res += d * w * s5 * trb(m_rx[5]) * teb(m_rz[5]);
                            assert(std::isfinite(s5));
                          }
                          else if constexpr (std::is_same_v<LHSRangeType, Math::Vector<NumberType>>)
                          {
                            assert(std::isfinite(s0));
                            trb(m_trv, m_rx[0]);
                            teb(m_tev, m_rz[0]);
                            res += d * w * s0 * m_tev.dot(m_k0 * m_trv);
                            assert(std::isfinite(s1));
                            trb(m_trv, m_rx[1]);
                            teb(m_tev, m_rz[1]);
                            res += d * w * s1 * m_tev.dot(m_k1 * m_trv);
                            assert(std::isfinite(s2));
                            trb(m_trv, m_rx[2]);
                            teb(m_tev, m_rz[2]);
                            res += d * w * s2 * m_tev.dot(m_k2 * m_trv);
                            assert(std::isfinite(s3));
                            trb(m_trv, m_rx[3]);
                            teb(m_tev, m_rz[3]);
                            res += d * w * s3 * m_tev.dot(m_k3 * m_trv);
                            assert(std::isfinite(s4));
                            trb(m_trv, m_rx[4]);
                            teb(m_tev, m_rz[4]);
                            res += d * w * s4 * m_tev.dot(m_k4 * m_trv);
                            assert(std::isfinite(s5));
                            trb(m_trv, m_rx[5]);
                            teb(m_tev, m_rz[5]);
                            res += d * w * s5 * m_tev.dot(m_k5 * m_trv);
                          }
                          else
                          {
                            assert(false);
                            res = NAN;
                          }
                        }
                      }
                    }
                  }
                }
              }
              break;
            }
          }
        }
        else
        {
          const auto& qftr = *m_qftr;
          const auto& qfte = *m_qfte;
          for (size_t i = 0; i < qfte.getSize(); i++)
          {
            const auto& x = m_x[i];
            for (size_t j = 0; j < qftr.getSize(); j++)
            {
              const auto& y = m_y[i];

              Real d;
              if constexpr (std::is_same_v<LHSRangeType, NumberType>)
              {
                d = kernel(x, y) * x.getDistortion() * y.getDistortion();
              }
              else if constexpr (std::is_same_v<LHSRangeType, Math::Vector<NumberType>>)
              {
                kernel(m_k, x, y);
                d = x.getDistortion() * y.getDistortion();
              }
              else
              {
                d = NAN;
              }

              const Real w = qfte.getWeight(i) * qftr.getWeight(j);
              for (size_t l = 0; l < testfe.getCount(); l++)
              {
                const auto& teb = testfe.getBasis(l);
                NumberType tev;
                if constexpr (std::is_same_v<LHSRangeType, NumberType>)
                {
                  tev = teb(qfte.getPoint(i));
                }
                else if constexpr (std::is_same_v<LHSRangeType, Math::Vector<NumberType>>)
                {
                  teb(m_tev, qfte.getPoint(i));
                }
                for (size_t m = 0; m < trialfe.getCount(); m++)
                {
                  const auto& trb = trialfe.getBasis(m);
                  if constexpr (std::is_same_v<LHSRangeType, NumberType>)
                  {
                    res += w * d * tev * trb(qftr.getPoint(j));
                  }
                  else if constexpr (std::is_same_v<LHSRangeType, Math::Vector<NumberType>>)
                  {
                    trb(m_trv, qftr.getPoint(j));
                    res += w * d * m_tev.dot(m_k * m_trv);
                  }
                }
              }
            }
          }
        }
        return res;
      }

      inline
      Region getTrialRegion() const override
      {
        return getIntegrand().getLHS().getRegion();
      }

      virtual Region getTestRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      std::optional<std::reference_wrapper<Geometry::Polytope>> m_trp;
      std::optional<std::reference_wrapper<Geometry::Polytope>> m_tep;

      std::unique_ptr<QF::QuadratureFormulaBase> m_qf;
      std::unique_ptr<QF::QuadratureFormulaBase> m_qftr;
      std::unique_ptr<QF::QuadratureFormulaBase> m_qfte;

      std::vector<Geometry::Point> m_x;
      std::vector<Geometry::Point> m_y;

      std::vector<Math::SpatialVector<NumberType>> m_rx;
      std::vector<Math::SpatialVector<NumberType>> m_rz;

      Math::Matrix<NumberType> m_k;
      Math::Vector<NumberType> m_trv, m_tev;
      Math::Matrix<NumberType> m_k0, m_k1, m_k2, m_k3, m_k4, m_k5;
  };

  template <class Kernel, class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  class Integral<
    Dot<
      Potential<Kernel, ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>>,
      ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>> final
    : public QuadratureRule<
        Dot<
          Potential<Kernel, ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>>,
          ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>
  {
    public:
      using KernelType = Kernel;

      using LHSType = Potential<KernelType, ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>>;

      using RHSType = ShapeFunctionBase<RHSDerived, TestFES, TestSpace>;

      using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;

      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

      using IntegrandType = Dot<LHSType, RHSType>;

      using Parent =
        QuadratureRule<
          Dot<
            Potential<KernelType, ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>>,
            ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>;

      Integral(const LHSType& lhs, const RHSType& rhs)
        : Integral(Dot(lhs, rhs))
      {}

      Integral(const IntegrandType& integrand)
        : Parent(integrand)
      {}

      Integral(const Integral& other)
        : Parent(other)
      {}

      Integral(Integral&& other)
        : Parent(std::move(other))
      {}

      inline
      Integrator::Region getTestRegion() const override
      {
        return Integrator::Region::Cells;
      }

      inline
      Integral* copy() const noexcept override
      {
        return new Integral(*this);
      }
  };

  template <class KernelType, class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  Integral(const
      Dot<
        Potential<KernelType, ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>>,
        ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>&)
    -> Integral<
          Dot<
            Potential<KernelType, ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>>,
            ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>;

  template <class KernelType, class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  Integral(const Potential<KernelType, ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>>&,
        const ShapeFunctionBase<RHSDerived, TestFES, TestSpace>&)
    -> Integral<
          Dot<
            Potential<KernelType, ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>>,
            ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>;
}

#endif

