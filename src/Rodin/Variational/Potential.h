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
  struct Traits<
    Variational::Potential<
      LHS,
      Variational::FunctionBase<RHSDerived>>>
  {
    using LHSType = LHS;

    using RHSType = Variational::FunctionBase<Variational::FunctionBase<RHSDerived>>;

    using KernelType = LHSType;

    using OperandType = RHSType;

    using RHSRangeType =
      typename FormLanguage::Traits<RHSType>::RangeType;

    using LHSRangeType =
      std::conditional_t<
      // If
      std::is_same_v<RHSRangeType, Scalar>,
      // Then
      Scalar,
      // Else
      std::conditional_t<
        // If
        std::is_same_v<RHSRangeType, Math::Vector<Scalar>>,
        // Then
        Math::Matrix<Scalar>,
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
      std::is_same_v<RHSRangeType, Scalar>,
      // Then
      Scalar,
      // Else
      std::conditional_t<
        // If
        std::is_same_v<RHSRangeType, Math::Vector<Scalar>>,
        // Then
        Math::Matrix<Scalar>,
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
      using LHSType = LHS;

      using KernelType = LHSType;

      using RHSType = FunctionBase<RHSDerived>;

      using OperandType = RHSType;

      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

      using LHSRangeType =
        std::conditional_t<
        // If
        std::is_same_v<RHSRangeType, Scalar>,
        // Then
        Scalar,
        // Else
        std::conditional_t<
          // If
          std::is_same_v<RHSRangeType, Math::Vector<Scalar>>,
          // Then
          Math::Matrix<Scalar>,
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
        if constexpr (std::is_same_v<LHSRangeType, Scalar>)
        {
          static_assert(std::is_same_v<RHSRangeType, Scalar>);
          return { 1, 1 };
        }
        else if constexpr (std::is_same_v<LHSRangeType, Math::Matrix<Scalar>>)
        {
          static_assert(std::is_same_v<RHSRangeType, Math::Vector<Scalar>>);
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
          if constexpr (std::is_same_v<RHSRangeType, Scalar>)
          {
            Scalar res = 0;
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
          else if constexpr (std::is_same_v<RHSRangeType, Math::Vector<Scalar>>)
          {
            Math::Vector<Scalar> res;
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
          if constexpr (std::is_same_v<RHSRangeType, Scalar>)
          {
            Scalar res = 0;
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
          else if constexpr (std::is_same_v<RHSRangeType, Math::Vector<Scalar>>)
          {
            Math::Vector<Scalar> res;
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

      void getValue(Math::Vector<Scalar>& res, const Geometry::Point& p) const
      {
        const auto& kernel = getKernel();
        const auto& operand = getOperand();
        const auto& mesh = p.getPolytope().getMesh();
        res.resize(getRangeShape().height());
        res.setZero();
        Math::Matrix<Scalar> kxy;
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

      using LHSType = LHS;

      using KernelType = LHS;

      using RHSType = ShapeFunctionBase<ShapeFunction<RHSDerived, FES, SpaceType>, FES, Space>;

      using OperandType = RHSType;

      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

      using LHSRangeType =
        std::conditional_t<
        // If
        std::is_same_v<RHSRangeType, Scalar>,
        // Then
        Scalar,
        // Else
        std::conditional_t<
          // If
          std::is_same_v<RHSRangeType, Math::Vector<Scalar>>,
          // Then
          Math::Matrix<Scalar>,
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
        if constexpr (std::is_same_v<LHSRangeType, Scalar>)
        {
          static_assert(std::is_same_v<RHSRangeType, Scalar>);
          return { 1, 1 };
        }
        else if constexpr (std::is_same_v<LHSRangeType, Math::Matrix<Scalar>>)
        {
          static_assert(std::is_same_v<RHSRangeType, Math::Vector<Scalar>>);
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
        : public GlobalBilinearFormIntegratorBase
  {
    public:
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

      void assemble(const Geometry::Polytope& trp, const Geometry::Polytope& tep) final override
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
        auto& res = getMatrix();
        res.resize(testfe.getCount(), trialfe.getCount());
        res.setZero();
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
              const QF::GenericPolytopeQuadrature qf(order, Geometry::Polytope::Type::Segment);
              Math::SpatialVector<Scalar> rx1(2), rz1(2),
                                  rx2(2), rz2(2),
                                  rx3(2), rz3(2),
                                  rx4(2), rz4(2),
                                  rx5(2), rz5(2),
                                  rx6(2), rz6(2);
              for (size_t i = 0; i < qf.getSize(); i++)
              {
                const Scalar xi = 1 - qf.getPoint(i).value();
                for (size_t j = 0; j < qf.getSize(); j++)
                {
                  const Scalar eta1 = qf.getPoint(j).value();
                  for (size_t h = 0; h < qf.getSize(); h++)
                  {
                    const Scalar eta2 = qf.getPoint(h).value();
                    rx1 << xi, xi * (1 - eta1 + eta1 * eta2);
                    rz2 << xi, xi * (1 - eta1 + eta1 * eta2);
                    rz3 << xi * (1 - eta1 * eta2), xi * eta1 * (1 - eta2);
                    rx4 << xi * (1 - eta1 * eta2), xi * eta1 * (1 - eta2);
                    rz5 << xi, xi * eta1 * (1 - eta2);
                    rx6 << xi, xi * eta1 * (1 - eta2);
                    const Geometry::Point x1(tep, teptrans, std::cref(rx1));
                    const Geometry::Point z2(tep, teptrans, std::cref(rz2));
                    const Geometry::Point z3(tep, teptrans, std::cref(rz3));
                    const Geometry::Point x4(tep, teptrans, std::cref(rx4));
                    const Geometry::Point z5(tep, teptrans, std::cref(rz5));
                    const Geometry::Point x6(tep, teptrans, std::cref(rx6));
                    const Scalar d = xi * xi * xi * eta1 * eta1 * eta2;
                    for (size_t k = 0; k < qf.getSize(); k++)
                    {
                      const Scalar w = qf.getWeight(i) * qf.getWeight(j) * qf.getWeight(h) * qf.getWeight(k);
                      const Scalar eta3 = qf.getPoint(k).value();
                      rz1 << xi * (1 - eta1 * eta2 * eta3), xi * (1 - eta1);
                      rx2 << xi * (1 - eta1 * eta2 * eta3), xi * (1 - eta1);
                      rx3 << xi, xi * eta1 * (1 - eta2 + eta2 * eta3);
                      rz4 << xi, xi * eta1 * (1 - eta2 + eta2 * eta3);
                      rx5 << xi * (1 - eta1 * eta2 * eta3), xi * eta1 * (1 - eta2 * eta3);
                      rz6 << xi * (1 - eta1 * eta2 * eta3), xi * eta1 * (1 - eta2 * eta3);
                      const Geometry::Point z1(tep, teptrans, std::cref(rz1));
                      const Geometry::Point x2(tep, teptrans, std::cref(rx2));
                      const Geometry::Point x3(tep, teptrans, std::cref(rx3));
                      const Geometry::Point z4(tep, teptrans, std::cref(rz4));
                      const Geometry::Point x5(tep, teptrans, std::cref(rx5));
                      const Geometry::Point z6(tep, teptrans, std::cref(rz6));

                      Scalar s1, s2, s3, s4, s5, s6;
                      if constexpr (std::is_same_v<LHSRangeType, Scalar>)
                      {
                        s1 = kernel(x1, z1) * x1.getDistortion() * z1.getDistortion();
                        s2 = kernel(x2, z2) * x2.getDistortion() * z2.getDistortion();
                        s3 = kernel(x3, z3) * x3.getDistortion() * z3.getDistortion();
                        s4 = kernel(x4, z4) * x4.getDistortion() * z4.getDistortion();
                        s5 = kernel(x5, z5) * x5.getDistortion() * z5.getDistortion();
                        s6 = kernel(x6, z6) * x6.getDistortion() * z6.getDistortion();
                      }
                      else if constexpr (std::is_same_v<LHSRangeType, Math::Vector<Scalar>>)
                      {
                        kernel(m_k1, x1, z1);
                        kernel(m_k2, x2, z2);
                        kernel(m_k3, x3, z3);
                        kernel(m_k4, x4, z4);
                        kernel(m_k5, x5, z5);
                        kernel(m_k6, x6, z6);
                        s1 = x1.getDistortion() * z1.getDistortion();
                        s2 = x2.getDistortion() * z2.getDistortion();
                        s3 = x3.getDistortion() * z3.getDistortion();
                        s4 = x4.getDistortion() * z4.getDistortion();
                        s5 = x5.getDistortion() * z5.getDistortion();
                        s6 = x6.getDistortion() * z6.getDistortion();
                      }
                      else
                      {
                        assert(false);
                        s1 = s2 = s3 = s4 = s5 = s6 = NAN;
                      }

                      for (size_t l = 0; l < testfe.getCount(); l++)
                      {
                        const auto& teb = testfe.getBasis(l);
                        for (size_t m = 0; m < trialfe.getCount(); m++)
                        {
                          const auto& trb = trialfe.getBasis(m);
                          if constexpr (std::is_same_v<LHSRangeType, Scalar>)
                          {
                            res(l, m) += d * w * s1 * trb(rx1) * teb(rz1);
                            assert(std::isfinite(s1));
                            res(l, m) += d * w * s2 * trb(rx2) * teb(rz2);
                            assert(std::isfinite(s2));
                            res(l, m) += d * w * s3 * trb(rx3) * teb(rz3);
                            assert(std::isfinite(s3));
                            res(l, m) += d * w * s4 * trb(rx4) * teb(rz4);
                            assert(std::isfinite(s4));
                            res(l, m) += d * w * s5 * trb(rx5) * teb(rz5);
                            assert(std::isfinite(s5));
                            res(l, m) += d * w * s6 * trb(rx6) * teb(rz6);
                            assert(std::isfinite(s6));
                          }
                          else if constexpr (std::is_same_v<LHSRangeType, Math::Vector<Scalar>>)
                          {
                            assert(std::isfinite(s1));
                            trb(m_trv, rx1);
                            teb(m_tev, rz1);
                            res(l, m) += d * w * s1 * m_tev.dot(m_k1 * m_trv);
                            assert(std::isfinite(s2));
                            trb(m_trv, rx2);
                            teb(m_tev, rz2);
                            res(l, m) += d * w * s2 * m_tev.dot(m_k2 * m_trv);
                            assert(std::isfinite(s3));
                            trb(m_trv, rx3);
                            teb(m_tev, rz3);
                            res(l, m) += d * w * s3 * m_tev.dot(m_k3 * m_trv);
                            assert(std::isfinite(s4));
                            trb(m_trv, rx4);
                            teb(m_tev, rz4);
                            res(l, m) += d * w * s4 * m_tev.dot(m_k4 * m_trv);
                            assert(std::isfinite(s5));
                            trb(m_trv, rx5);
                            teb(m_tev, rz5);
                            res(l, m) += d * w * s5 * m_tev.dot(m_k5 * m_trv);
                            assert(std::isfinite(s6));
                            trb(m_trv, rx6);
                            teb(m_tev, rz6);
                            res(l, m) += d * w * s6 * m_tev.dot(m_k6 * m_trv);
                          }
                          else
                          {
                            assert(false);
                            res(l, m) = NAN;
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
          const QF::GenericPolytopeQuadrature qftr(trialfe.getOrder(), trp.getGeometry());
          const QF::GenericPolytopeQuadrature qfte(testfe.getOrder(), tep.getGeometry());
          for (size_t i = 0; i < qfte.getSize(); i++)
          {
            const Geometry::Point x(tep, teptrans, std::cref(qfte.getPoint(i)));
            for (size_t j = 0; j < qftr.getSize(); j++)
            {
              const Geometry::Point y(trp, trptrans, std::cref(qftr.getPoint(j)));
              Scalar d;
              if constexpr (std::is_same_v<LHSRangeType, Scalar>)
              {
                d = kernel(x, y) * x.getDistortion() * y.getDistortion();
              }
              else if constexpr (std::is_same_v<LHSRangeType, Math::Vector<Scalar>>)
              {
                kernel(m_k, x, y);
                d = x.getDistortion() * y.getDistortion();
              }
              const Scalar w = qfte.getWeight(i) * qftr.getWeight(j);
              for (size_t l = 0; l < testfe.getCount(); l++)
              {
                const auto& teb = testfe.getBasis(l);
                Scalar tev;
                if constexpr (std::is_same_v<LHSRangeType, Scalar>)
                {
                  tev = teb(qfte.getPoint(i));
                }
                else if constexpr (std::is_same_v<LHSRangeType, Math::Vector<Scalar>>)
                {
                  teb(m_tev, qfte.getPoint(i));
                }
                for (size_t m = 0; m < trialfe.getCount(); m++)
                {
                  const auto& trb = trialfe.getBasis(m);
                  if constexpr (std::is_same_v<LHSRangeType, Scalar>)
                  {
                    res(l, m) += w * d * tev * trb(qftr.getPoint(j));
                  }
                  else if constexpr (std::is_same_v<LHSRangeType, Math::Vector<Scalar>>)
                  {
                    trb(m_trv, qftr.getPoint(j));
                    res(l, m) += w * d * m_tev.dot(m_k * m_trv);
                  }
                }
              }
            }
          }
        }
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

      mutable Math::Matrix<Scalar> m_k;
      mutable Math::Vector<Scalar> m_trv, m_tev;
      mutable Math::Matrix<Scalar> m_k1, m_k2, m_k3, m_k4, m_k5, m_k6;
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

