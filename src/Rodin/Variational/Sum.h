/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Addition operation for functions, shape functions, and form integrators.
 */

#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_SUM_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_SUM_H

#include "Rodin/FormLanguage/List.h"
#include "Rodin/Math/Traits.h"

#include "ForwardDecls.h"
#include "Function.h"
#include "Rodin/Variational/IntegrationPoint.h"
#include "ShapeFunction.h"

#include "RealFunction.h"
#include "ComplexFunction.h"

#include "LinearFormIntegrator.h"
#include "BilinearFormIntegrator.h"

namespace Rodin::FormLanguage
{
  template <class LHSDerived, class RHSDerived>
  struct Traits<
    Variational::Sum<Variational::FunctionBase<LHSDerived>, Variational::FunctionBase<RHSDerived>>>
  {
    using LHSType = Variational::FunctionBase<LHSDerived>;
    using RHSType = Variational::FunctionBase<RHSDerived>;
  };

  template <class LHSDerived, class RHSDerived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Sum<
      Variational::ShapeFunctionBase<LHSDerived, FES, Space>,
      Variational::ShapeFunctionBase<RHSDerived, FES, Space>>>
  {
    using FESType = FES;
    using LHSType = Variational::FunctionBase<LHSDerived>;
    using RHSType = Variational::FunctionBase<RHSDerived>;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;
  };

  template <class LHSNumber, class RHSNumber>
  struct Traits<
    Variational::Sum<
      Variational::LinearFormIntegratorBase<LHSNumber>,
      Variational::LinearFormIntegratorBase<RHSNumber>>>
  {
    using LHSScalarType = LHSNumber;

    using RHSScalarType = RHSNumber;

    using LHSType = Variational::LinearFormIntegratorBase<LHSScalarType>;

    using RHSType = Variational::LinearFormIntegratorBase<RHSScalarType>;

    using ScalarType =
      typename FormLanguage::Sum<LHSScalarType, RHSScalarType>::Type;
  };

}

namespace Rodin::Variational
{
  /**
   * @defgroup SumSpecializations Sum Template Specializations
   * @brief Template specializations of the Sum class.
   * @see Sum
   */

  /**
   * @brief Addition of two functions.
   *
   * Computes pointwise sum of two functions:
   * @f[
   *    (f + g)(x) = f(x) + g(x)
   * @f]
   *
   * Both functions must have the same range type (scalar, vector, or matrix).
   *
   * @tparam LHSDerived Type of left operand function
   * @tparam RHSDerived Type of right operand function
   * @ingroup SumSpecializations
   */
  template <class LHSDerived, class RHSDerived>
  class Sum<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>> final
    : public FunctionBase<Sum<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>
  {
    public:
      using LHSType = FunctionBase<LHSDerived>;
      using RHSType = FunctionBase<RHSDerived>;

      using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;
      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

      using Parent = FunctionBase<Sum<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>;
      static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);

      /**
       * @brief Constructs sum of two functions.
       * @param lhs Left operand function
       * @param rhs Right operand function
       */
      constexpr
      Sum(const LHSType& lhs, const RHSType& rhs)
        : m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {}

      constexpr
      Sum(const Sum& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

      constexpr
      Sum(Sum&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      /**
       * @brief Gets the left operand.
       * @returns Reference to left function
       */
      constexpr
      const auto& getLHS() const
      {
        assert(m_lhs);
        return *m_lhs;
      }

      /**
       * @brief Gets the right operand.
       * @returns Reference to right function
       */
      constexpr
      const auto& getRHS() const
      {
        assert(m_rhs);
        return *m_rhs;
      }

      /**
       * @brief Restricts both functions to a trace.
       * @param args Arguments for trace restriction
       * @returns Reference to this object
       */
      template <class ... Args>
      constexpr
      Sum& traceOf(const Args& ... args)
      {
        m_lhs->traceOf(args...);
        m_rhs->traceOf(args...);
        return *this;
      }

      /**
       * @brief Evaluates the sum at a point.
       * @param p Point at which to evaluate
       * @returns Sum @f$ f(p) + g(p) @f$
       */
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return this->object(this->getLHS().getValue(p)) + this->object(this->getRHS().getValue(p));
      }

      constexpr
      std::optional<size_t> getOrder(const Geometry::Polytope& poly) const noexcept
      {
        const auto lo = getLHS().getOrder(poly);
        const auto ro = getRHS().getOrder(poly);

        if (lo && ro)
          return std::max(*lo, *ro);
        else
          return std::nullopt;
      }

      Sum* copy() const noexcept override
      {
        return new Sum(*this);
      }

    private:
      std::unique_ptr<LHSType> m_lhs;
      std::unique_ptr<RHSType> m_rhs;
  };

  /**
   * @brief Deduction guide for function sum.
   */
  template <class LHSDerived, class RHSDerived>
  Sum(const FunctionBase<LHSDerived>&, const FunctionBase<RHSDerived>&)
    -> Sum<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>;

  /**
   * @brief Addition operator for two functions.
   * @param lhs Left operand function
   * @param rhs Right operand function
   * @returns Sum object
   * @relates Sum
   */
  template <class LHSDerived, class RHSDerived>
  constexpr
  auto
  operator+(const FunctionBase<LHSDerived>& lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return Sum(lhs, rhs);
  }

  /**
   * @brief Addition of function and scalar.
   * @param lhs Function operand
   * @param rhs Real scalar
   * @returns Sum object
   * @relates Sum
   */
  template <class LHSDerived>
  constexpr
  auto
  operator+(const FunctionBase<LHSDerived>& lhs, Real rhs)
  {
    return Sum(lhs, RealFunction(rhs));
  }

  /**
   * @brief Addition of scalar and function.
   * @param lhs Real scalar
   * @param rhs Function operand
   * @returns Sum object
   * @relates Sum
   */
  template <class RHSDerived>
  constexpr
  auto
  operator+(Real lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return Sum(RealFunction(lhs), rhs);
  }

  /**
   * @brief Addition of function and complex scalar.
   * @param lhs Function operand
   * @param rhs Complex scalar
   * @returns Sum object
   * @relates Sum
   */
  template <class LHSDerived>
  constexpr
  auto
  operator+(const FunctionBase<LHSDerived>& lhs, Complex rhs)
  {
    return Sum(lhs, ComplexFunction(rhs));
  }

  /**
   * @brief Addition of complex scalar and function.
   * @param lhs Complex scalar
   * @param rhs Function operand
   * @returns Sum object
   * @relates Sum
   */
  template <class RHSDerived>
  constexpr
  auto
  operator+(Complex lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return Sum(ComplexFunction(lhs), rhs);
  }

  /**
   * @brief Addition of two shape functions.
   *
   * Computes basis-wise sum of two trial or test functions in the same
   * finite element space:
   * @f[
   *    (u + v)_i = u_i + v_i
   * @f]
   *
   * @tparam LHSDerived Type of left shape function
   * @tparam RHSDerived Type of right shape function
   * @tparam FES Finite element space type
   * @tparam Space Trial or test function space
   * @ingroup SumSpecializations
   */
  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  class Sum<ShapeFunctionBase<LHSDerived, FES, Space>, ShapeFunctionBase<RHSDerived, FES, Space>> final
    : public ShapeFunctionBase<Sum<ShapeFunctionBase<LHSDerived, FES, Space>, ShapeFunctionBase<RHSDerived, FES, Space>>>
  {
    public:
      using FESType = FES;
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      using LHSType = ShapeFunctionBase<LHSDerived, FES, Space>;

      using RHSType = ShapeFunctionBase<RHSDerived, FES, Space>;

      using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;

      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

      using Parent = ShapeFunctionBase<Sum<LHSType, RHSType>, FES, Space>;
      static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);

      constexpr
      Sum(const LHSType& lhs, const RHSType& rhs)
        : Parent(lhs.getFiniteElementSpace()),
          m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {
        assert(lhs.getLeaf().getUUID() == rhs.getLeaf().getUUID());
      }

      constexpr
      Sum(const Sum& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

      constexpr
      Sum(Sum&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      constexpr
      const LHSType& getLHS() const
      {
        assert(m_lhs);
        return *m_lhs;
      }

      constexpr
      const RHSType& getRHS() const
      {
        assert(m_rhs);
        return *m_rhs;
      }

      constexpr
      const auto& getLeaf() const
      {
        return getRHS().getLeaf();
      }

      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        assert(getLHS().getDOFs(element) == getRHS().getDOFs(element));
        return getLHS().getDOFs(element);
      }

      Sum& setIntegrationPoint(const IntegrationPoint& ip)
      {
        m_lhs->setIntegrationPoint(ip);
        m_rhs->setIntegrationPoint(ip);
        return *this;
      }

      const IntegrationPoint& getIntegrationPoint() const
      {
        return m_lhs->getIntegrationPoint();
      }

      constexpr
      auto getBasis(size_t local) const
      {
        return this->object(getLHS().getBasis(local)) + this->object(getRHS().getBasis(local));
      }

      constexpr
      const auto& getFiniteElementSpace() const
      {
        return getLHS().getFiniteElementSpace();
      }

      constexpr
      std::optional<size_t> getOrder(const Geometry::Polytope& poly) const noexcept
      {
        const auto lo = getLHS().getOrder(poly);
        const auto ro = getRHS().getOrder(poly);

        if (lo && ro)
          return std::max(*lo, *ro);
        else
          return std::nullopt;
      }

      Sum* copy() const noexcept override
      {
        return new Sum(*this);
      }

    private:
      std::unique_ptr<LHSType> m_lhs;
      std::unique_ptr<RHSType> m_rhs;
  };

  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  Sum(const ShapeFunctionBase<LHSDerived, FES, Space>&, const ShapeFunctionBase<RHSDerived, FES, Space>&)
    -> Sum<ShapeFunctionBase<LHSDerived, FES, Space>, ShapeFunctionBase<RHSDerived, FES, Space>>;

  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  constexpr
  auto
  operator+(const ShapeFunctionBase<LHSDerived, FES, Space>& lhs,
            const ShapeFunctionBase<RHSDerived, FES, Space>& rhs)
  {
    return Sum(lhs, rhs);
  }

  template <class LHSNumber, class RHSNumber>
  class Sum<LinearFormIntegratorBase<LHSNumber>, LinearFormIntegratorBase<RHSNumber>>
    : public FormLanguage::List<
        LinearFormIntegratorBase<typename FormLanguage::Sum<LHSNumber, RHSNumber>::Type>>
  {
    public:
      using LHSScalarType = LHSNumber;

      using RHSScalarType = RHSNumber;

      using LHSType = LinearFormIntegratorBase<LHSScalarType>;

      using RHSType = LinearFormIntegratorBase<RHSScalarType>;

      using ScalarType = typename FormLanguage::Sum<LHSNumber, RHSNumber>::Type;

      using Parent = FormLanguage::List<LinearFormIntegratorBase<ScalarType>>;

      Sum(const LHSType& lhs, const RHSType& rhs)
      {
        this->add(lhs);
        this->add(rhs);
      }

      Sum(const Sum& other)
        : Parent(other)
      {}

      Sum(Sum&& other)
        : Parent(std::move(other))
      {}
  };

  template <class LHSNumber, class RHSNumber>
  Sum(const LinearFormIntegratorBase<LHSNumber>&, const LinearFormIntegratorBase<RHSNumber>&)
    -> Sum<LinearFormIntegratorBase<LHSNumber>, LinearFormIntegratorBase<RHSNumber>>;

  template <class LHSNumber, class RHSNumber>
  constexpr
  auto
  operator+(
      const LinearFormIntegratorBase<LHSNumber>& lhs,
      const LinearFormIntegratorBase<RHSNumber>& rhs)
  {
    return Sum(lhs, rhs);
  }

  template <class LHSNumber, class RHSNumber>
  class Sum<LinearFormIntegratorBase<LHSNumber>, FormLanguage::List<LinearFormIntegratorBase<RHSNumber>>>
    : public FormLanguage::List<
              LinearFormIntegratorBase<typename FormLanguage::Sum<LHSNumber, RHSNumber>::Type>>
  {
    public:
      using LHSScalarType = LHSNumber;

      using RHSScalarType = RHSNumber;

      using LHSType = LinearFormIntegratorBase<LHSScalarType>;

      using RHSType = FormLanguage::List<LinearFormIntegratorBase<RHSNumber>>;

      using ScalarType = typename FormLanguage::Sum<LHSNumber, RHSNumber>::Type;

      using Parent = FormLanguage::List<LinearFormIntegratorBase<ScalarType>>;

      Sum(const LHSType& lhs, const RHSType& rhs)
      {
        this->add(lhs);
        this->add(rhs);
      }

      Sum(const Sum& other)
        : Parent(other)
      {}

      Sum(Sum&& other)
        : Parent(std::move(other))
      {}
  };

  template <class LHSNumber, class RHSNumber>
  Sum(const LinearFormIntegratorBase<LHSNumber>&, const FormLanguage::List<LinearFormIntegratorBase<RHSNumber>&>)
    -> Sum<LinearFormIntegratorBase<LHSNumber>, FormLanguage::List<LinearFormIntegratorBase<RHSNumber>>>;

  template <class LHSNumber, class RHSNumber>
  constexpr
  auto
  operator+(
      const LinearFormIntegratorBase<LHSNumber>& lhs,
      const FormLanguage::List<LinearFormIntegratorBase<RHSNumber>>& rhs)
  {
    return Sum(lhs, rhs);
  }

  template <class LHSNumber, class RHSNumber>
  class Sum<
    FormLanguage::List<LinearFormIntegratorBase<LHSNumber>>, LinearFormIntegratorBase<RHSNumber>>
    : public FormLanguage::List<
              LinearFormIntegratorBase<typename FormLanguage::Sum<LHSNumber, RHSNumber>::Type>>
  {
    public:
      using LHSScalarType = LHSNumber;

      using RHSScalarType = RHSNumber;

      using LHSType = FormLanguage::List<LinearFormIntegratorBase<LHSNumber>>;

      using RHSType = LinearFormIntegratorBase<RHSScalarType>;

      using ScalarType = typename FormLanguage::Sum<LHSNumber, RHSNumber>::Type;

      using Parent = FormLanguage::List<LinearFormIntegratorBase<ScalarType>>;

      Sum(const LHSType& lhs, const RHSType& rhs)
      {
        this->add(lhs);
        this->add(rhs);
      }

      Sum(const Sum& other)
        : Parent(other)
      {}

      Sum(Sum&& other)
        : Parent(std::move(other))
      {}
  };

  template <class LHSNumber, class RHSNumber>
  Sum(const FormLanguage::List<LinearFormIntegratorBase<RHSNumber>>& lhs,
      const LinearFormIntegratorBase<LHSNumber>& rhs)
    -> Sum<FormLanguage::List<LinearFormIntegratorBase<LHSNumber>>, LinearFormIntegratorBase<RHSNumber>>;

  template <class LHSNumber, class RHSNumber>
  constexpr
  auto
  operator+(
      const FormLanguage::List<LinearFormIntegratorBase<RHSNumber>>& lhs,
      const LinearFormIntegratorBase<LHSNumber>& rhs)
  {
    return Sum(lhs, rhs);
  }

  template <class LHSNumber, class RHSNumber>
  class Sum<
    FormLanguage::List<LinearFormIntegratorBase<LHSNumber>>,
    FormLanguage::List<LinearFormIntegratorBase<RHSNumber>>>
      : public FormLanguage::List<
          LinearFormIntegratorBase<typename FormLanguage::Sum<LHSNumber, RHSNumber>::Type>>
  {
    public:
      using LHSScalarType = LHSNumber;

      using RHSScalarType = RHSNumber;

      using LHSType = FormLanguage::List<LinearFormIntegratorBase<RHSNumber>>;

      using RHSType = FormLanguage::List<LinearFormIntegratorBase<RHSNumber>>;

      using ScalarType = typename FormLanguage::Sum<LHSNumber, RHSNumber>::Type;

      using Parent = FormLanguage::List<LinearFormIntegratorBase<ScalarType>>;

      Sum(const LHSType& lhs, const RHSType& rhs)
      {
        this->add(lhs);
        this->add(rhs);
      }

      Sum(const Sum& other)
        : Parent(other)
      {}

      Sum(Sum&& other)
        : Parent(std::move(other))
      {}
  };

  template <class LHSNumber, class RHSNumber>
  Sum(const FormLanguage::List<LinearFormIntegratorBase<LHSNumber>>&,
      const FormLanguage::List<LinearFormIntegratorBase<RHSNumber>>&)
    -> Sum<
        FormLanguage::List<LinearFormIntegratorBase<LHSNumber>>,
        FormLanguage::List<LinearFormIntegratorBase<RHSNumber>>>;

  template <class LHSNumber, class RHSNumber>
  constexpr
  auto
  operator+(
      const FormLanguage::List<LinearFormIntegratorBase<LHSNumber>>& lhs,
      const FormLanguage::List<LinearFormIntegratorBase<RHSNumber>>& rhs)
  {
    return Sum(lhs, rhs);
  }

  template <class LHSNumber, class RHSNumber>
  class Sum<LocalBilinearFormIntegratorBase<LHSNumber>, LocalBilinearFormIntegratorBase<RHSNumber>>
    : public FormLanguage::List<
              LocalBilinearFormIntegratorBase<typename FormLanguage::Sum<LHSNumber, RHSNumber>::Type>>
  {
    public:
      using LHSScalarType = LHSNumber;

      using RHSScalarType = RHSNumber;

      using LHSType = LocalBilinearFormIntegratorBase<LHSScalarType>;

      using RHSType = LocalBilinearFormIntegratorBase<RHSScalarType>;

      using ScalarType = typename FormLanguage::Sum<LHSNumber, RHSNumber>::Type;

      using Parent = FormLanguage::List<LocalBilinearFormIntegratorBase<ScalarType>>;

      Sum(const LHSType& lhs, const RHSType& rhs)
      {
        this->add(lhs);
        this->add(rhs);
      }

      Sum(const Sum& other)
        : Parent(other)
      {}

      Sum(Sum&& other)
        : Parent(std::move(other))
      {}
  };

  template <class LHSNumber, class RHSNumber>
  Sum(const LocalBilinearFormIntegratorBase<LHSNumber>& lhs,
      const LocalBilinearFormIntegratorBase<RHSNumber>& rhs)
    -> Sum<LocalBilinearFormIntegratorBase<LHSNumber>, LocalBilinearFormIntegratorBase<RHSNumber>>;

  template <class LHSNumber, class RHSNumber>
  constexpr
  auto
  operator+(
    const LocalBilinearFormIntegratorBase<LHSNumber>& lhs,
    const LocalBilinearFormIntegratorBase<RHSNumber>& rhs)
  {
    return Sum(lhs, rhs);
  }

  template <class LHSNumber, class RHSNumber>
  class Sum<LocalBilinearFormIntegratorBase<LHSNumber>, FormLanguage::List<LocalBilinearFormIntegratorBase<RHSNumber>>>
    : public FormLanguage::List<
              LocalBilinearFormIntegratorBase<typename FormLanguage::Sum<LHSNumber, RHSNumber>::Type>>
  {
    public:
      using LHSScalarType = LHSNumber;

      using RHSScalarType = RHSNumber;

      using LHSType = LocalBilinearFormIntegratorBase<LHSScalarType>;

      using RHSType = FormLanguage::List<LocalBilinearFormIntegratorBase<RHSNumber>>;

      using ScalarType = typename FormLanguage::Sum<LHSNumber, RHSNumber>::Type;

      using Parent = FormLanguage::List<LocalBilinearFormIntegratorBase<ScalarType>>;

      Sum(const LHSType& lhs, const RHSType& rhs)
      {
        this->add(lhs);
        this->add(rhs);
      }

      Sum(const Sum& other)
        : Parent(other)
      {}

      Sum(Sum&& other)
        : Parent(std::move(other))
      {}
  };

  template <class LHSNumber, class RHSNumber>
  Sum(const LocalBilinearFormIntegratorBase<LHSNumber>&,
      const FormLanguage::List<LocalBilinearFormIntegratorBase<RHSNumber>>&)
    -> Sum<
        LocalBilinearFormIntegratorBase<LHSNumber>,
        FormLanguage::List<LocalBilinearFormIntegratorBase<RHSNumber>>>;

  template <class LHSNumber, class RHSNumber>
  constexpr
  auto
  operator+(
      const LocalBilinearFormIntegratorBase<LHSNumber>& lhs,
      const FormLanguage::List<LocalBilinearFormIntegratorBase<RHSNumber>>& rhs)
  {
    return Sum(lhs, rhs);
  }

  template <class LHSNumber, class RHSNumber>
  class Sum<
    FormLanguage::List<LocalBilinearFormIntegratorBase<LHSNumber>>, LocalBilinearFormIntegratorBase<RHSNumber>>
    : public FormLanguage::List<
              LocalBilinearFormIntegratorBase<typename FormLanguage::Sum<LHSNumber, RHSNumber>::Type>>
  {
    public:
      using LHSScalarType = LHSNumber;

      using RHSScalarType = RHSNumber;

      using LHSType = FormLanguage::List<LocalBilinearFormIntegratorBase<LHSNumber>>;

      using RHSType = LocalBilinearFormIntegratorBase<RHSScalarType>;

      using ScalarType = typename FormLanguage::Sum<LHSNumber, RHSNumber>::Type;

      using Parent = FormLanguage::List<LocalBilinearFormIntegratorBase<ScalarType>>;

      Sum(const LHSType& lhs, const RHSType& rhs)
      {
        this->add(lhs);
        this->add(rhs);
      }

      Sum(const Sum& other)
        : Parent(other)
      {}

      Sum(Sum&& other)
        : Parent(std::move(other))
      {}
  };

  template <class LHSNumber, class RHSNumber>
  Sum(const FormLanguage::List<LocalBilinearFormIntegratorBase<LHSNumber>>&,
      const LocalBilinearFormIntegratorBase<RHSNumber>&)
    -> Sum<
        FormLanguage::List<LocalBilinearFormIntegratorBase<LHSNumber>>,
        LocalBilinearFormIntegratorBase<RHSNumber>>;

  template <class LHSNumber, class RHSNumber>
  constexpr
  auto
  operator+(
      const FormLanguage::List<LocalBilinearFormIntegratorBase<LHSNumber>>& lhs,
      const LocalBilinearFormIntegratorBase<RHSNumber>& rhs)
  {
    return Sum(lhs, rhs);
  }

  template <class LHSNumber, class RHSNumber>
  class Sum<
    FormLanguage::List<LocalBilinearFormIntegratorBase<LHSNumber>>,
    FormLanguage::List<LocalBilinearFormIntegratorBase<RHSNumber>>>
      : public FormLanguage::List<
          LocalBilinearFormIntegratorBase<typename FormLanguage::Sum<LHSNumber, RHSNumber>::Type>>
  {
    public:
      using LHSScalarType = LHSNumber;

      using RHSScalarType = RHSNumber;

      using LHSType = FormLanguage::List<LocalBilinearFormIntegratorBase<RHSNumber>>;

      using RHSType = FormLanguage::List<LocalBilinearFormIntegratorBase<RHSNumber>>;

      using ScalarType = typename FormLanguage::Sum<LHSNumber, RHSNumber>::Type;

      using Parent = FormLanguage::List<LocalBilinearFormIntegratorBase<ScalarType>>;

      Sum(const LHSType& lhs, const RHSType& rhs)
      {
        this->add(lhs);
        this->add(rhs);
      }

      Sum(const Sum& other)
        : Parent(other)
      {}

      Sum(Sum&& other)
        : Parent(std::move(other))
      {}
  };

  template <class LHSNumber, class RHSNumber>
  Sum(const FormLanguage::List<LocalBilinearFormIntegratorBase<LHSNumber>>&,
      const FormLanguage::List<LocalBilinearFormIntegratorBase<RHSNumber>>&)
    -> Sum<
        FormLanguage::List<LocalBilinearFormIntegratorBase<LHSNumber>>,
        FormLanguage::List<LocalBilinearFormIntegratorBase<RHSNumber>>>;

  template <class LHSNumber, class RHSNumber>
  constexpr
  auto
  operator+(
      const FormLanguage::List<LocalBilinearFormIntegratorBase<LHSNumber>>& lhs,
      const FormLanguage::List<LocalBilinearFormIntegratorBase<RHSNumber>>& rhs)
  {
    return Sum(lhs, rhs);
  }

  template <class Operator>
  class Sum<BilinearFormBase<Operator>, BilinearFormBase<Operator>>
    : public FormLanguage::List<BilinearFormBase<Operator>>
  {
    public:
      using LHSType = BilinearFormBase<Operator>;

      using RHSType = BilinearFormBase<Operator>;

      using Parent = FormLanguage::List<BilinearFormBase<Operator>>;

      Sum(const LHSType& lhs, const RHSType& rhs)
      {
        this->add(lhs);
        this->add(rhs);
      }

      Sum(const Sum& other)
        : Parent(other)
      {}

      Sum(Sum&& other)
        : Parent(std::move(other))
      {}
  };

  template <class Operator>
  Sum(const BilinearFormBase<Operator>& lhs, const BilinearFormBase<Operator>& rhs)
    -> Sum<BilinearFormBase<Operator>, BilinearFormBase<Operator>>;

  template <class Operator>
  constexpr
  auto
  operator+(const BilinearFormBase<Operator>& lhs, const BilinearFormBase<Operator>& rhs)
  {
    return Sum(lhs, rhs);
  }
}

#endif
