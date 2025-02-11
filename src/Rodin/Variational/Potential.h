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
    using ScalarType = Real;

    using LHSType = LHS;

    using RHSType = Variational::FunctionBase<Variational::FunctionBase<RHSDerived>>;

    using KernelType = LHSType;

    using OperandType = RHSType;

    using RHSRangeType =
      typename FormLanguage::Traits<RHSType>::RangeType;

    using LHSRangeType =
      std::conditional_t<
      // If
      std::is_same_v<RHSRangeType, ScalarType>,
      // Then
      ScalarType,
      // Else
      std::conditional_t<
        // If
        std::is_same_v<RHSRangeType, Math::Vector<ScalarType>>,
        // Then
        Math::Matrix<ScalarType>,
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
    using ScalarType = Real;

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
      std::is_same_v<RHSRangeType, ScalarType>,
      // Then
      ScalarType,
      // Else
      std::conditional_t<
        // If
        std::is_same_v<RHSRangeType, Math::Vector<ScalarType>>,
        // Then
        Math::Matrix<ScalarType>,
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
      using ScalarType = Real;

      using LHSType = LHS;

      using KernelType = LHSType;

      using RHSType = FunctionBase<RHSDerived>;

      using OperandType = RHSType;

      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

      using LHSRangeType =
        std::conditional_t<
        // If
        std::is_same_v<RHSRangeType, ScalarType>,
        // Then
        ScalarType,
        // Else
        std::conditional_t<
          // If
          std::is_same_v<RHSRangeType, Math::Vector<ScalarType>>,
          // Then
          Math::Matrix<ScalarType>,
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

      constexpr
      RangeShape getRangeShape() const
      {
        if constexpr (std::is_same_v<LHSRangeType, ScalarType>)
        {
          static_assert(std::is_same_v<RHSRangeType, ScalarType>);
          return { 1, 1 };
        }
        else if constexpr (std::is_same_v<LHSRangeType, Math::Matrix<ScalarType>>)
        {
          static_assert(std::is_same_v<RHSRangeType, Math::Vector<ScalarType>>);
          return getOperand().getRangeShape();
        }
        else
        {
          assert(false);
          return { 0, 0 };
        }
      }

      const auto& getKernel() const
      {
        return m_kernel.get();
      }

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
          if constexpr (std::is_same_v<RHSRangeType, ScalarType>)
          {
            ScalarType res = 0;
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
          else if constexpr (std::is_same_v<RHSRangeType, Math::Vector<ScalarType>>)
          {
            Math::Vector<ScalarType> res;
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
          if constexpr (std::is_same_v<RHSRangeType, ScalarType>)
          {
            ScalarType res = 0;
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
          else if constexpr (std::is_same_v<RHSRangeType, Math::Vector<ScalarType>>)
          {
            Math::Vector<ScalarType> res;
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

      void getValue(Math::Vector<ScalarType>& res, const Geometry::Point& p) const
      {
        const auto& kernel = getKernel();
        const auto& operand = getOperand();
        const auto& mesh = p.getPolytope().getMesh();
        res.resize(getRangeShape().height());
        res.setZero();
        Math::Matrix<ScalarType> kxy;
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

      const auto& getQuadratureFormula() const
      {
        return m_qf;
      }

      Potential* copy() const noexcept override
      {
        return new Potential(*this);
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

      using ScalarType = Real;

      using LHSType = LHS;

      using KernelType = LHS;

      using RHSType = ShapeFunctionBase<ShapeFunction<RHSDerived, FES, SpaceType>, FES, Space>;

      using OperandType = RHSType;

      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

      using LHSRangeType =
        std::conditional_t<
        // If
        std::is_same_v<RHSRangeType, ScalarType>,
        // Then
        ScalarType,
        // Else
        std::conditional_t<
          // If
          std::is_same_v<RHSRangeType, Math::Vector<ScalarType>>,
          // Then
          Math::Matrix<ScalarType>,
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

      const KernelType& getKernel() const
      {
        return m_kernel;
      }

      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      Integrator::Region getRegion() const
      {
        return Integrator::Region::Cells;
      }

      constexpr
      RangeShape getRangeShape() const
      {
        if constexpr (std::is_same_v<LHSRangeType, ScalarType>)
        {
          static_assert(std::is_same_v<RHSRangeType, ScalarType>);
          return { 1, 1 };
        }
        else if constexpr (std::is_same_v<LHSRangeType, Math::Matrix<ScalarType>>)
        {
          static_assert(std::is_same_v<RHSRangeType, Math::Vector<ScalarType>>);
          return getOperand().getRangeShape();
        }
        else
        {
          assert(false);
          return { 0, 0 };
        }
      }

      Potential* copy() const noexcept override
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

      Integrator::Region getTestRegion() const override
      {
        return Integrator::Region::Cells;
      }

      Integral* copy() const noexcept override
      {
        return new Integral(*this);
      }
  };

  template <class KernelType, class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  Integral(
      const Dot<Potential<KernelType, ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>>,
      ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>&)
    -> Integral<
          Dot<
            Potential<KernelType, ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>>,
            ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>;

  template <class KernelType, class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  Integral(
      const Potential<KernelType, ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>>&,
      const ShapeFunctionBase<RHSDerived, TestFES, TestSpace>&)
    -> Integral<
          Dot<
            Potential<KernelType, ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>>,
            ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>;
}

#endif

