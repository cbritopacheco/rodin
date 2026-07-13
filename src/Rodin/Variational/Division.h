/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DIVISION_H
#define RODIN_VARIATIONAL_DIVISION_H

/**
 * @file
 * @brief Division operation for functions.
 */

#include "ForwardDecls.h"
#include "Function.h"

namespace Rodin::FormLanguage
{
  template <class LHSDerived, class RHSDerived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Division<
      Variational::ShapeFunctionBase<LHSDerived, FES, Space>,
      Variational::FunctionBase<RHSDerived>>>
  {
    /// @brief Finite element space type.
      using FESType = FES;
      static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;

    /// @brief Left-hand side operand type.
      using LHSType = Variational::ShapeFunctionBase<LHSDerived, FES, Space>;

    /// @brief Right-hand side operand type.
      using RHSType = Variational::FunctionBase<RHSDerived>;

    /// @brief Range (evaluation value) type.
      using RangeType = typename FormLanguage::Traits<LHSType>::RangeType;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup DivisionSpecializations Division Template Specializations
   * @brief Template specializations of the Division class.
   * @see Division
   */

  /**
   * @ingroup DivisionSpecializations
   * @brief Division of a function by another function.
   *
   * Represents the mathematical expression:
   * @f[
   *    \left(\frac{f}{g}\right)(x) = \frac{f(x)}{g(x)}
   * @f]
   * where @f$ f @f$ and @f$ g @f$ are functions with compatible ranges.
   *
   * @tparam LHSDerived Type of the numerator function
   * @tparam RHSDerived Type of the denominator function
   *
   * @note The denominator function @f$ g @f$ must not evaluate to zero on the
   * domain of interest.
   */
  template <class LHSDerived, class RHSDerived>
  class Division<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>
    : public FunctionBase<Division<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>
  {
    public:
      /// @brief Left-hand side operand type.
      using LHSType = FunctionBase<LHSDerived>;

      /// @brief Right-hand side operand type.
      using RHSType = FunctionBase<RHSDerived>;

      /// @brief Range type of the left-hand side operand.
      using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;

      /// @brief Range type of the right-hand side operand.
      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

      /// @brief Parent class type.
      using Parent =
        FunctionBase<Division<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>;

      /**
       * @brief Constructs a division from numerator and denominator functions.
       * @param lhs Numerator function @f$ f @f$
       * @param rhs Denominator function @f$ g @f$
       */
      Division(const FunctionBase<LHSDerived>& lhs, const FunctionBase<RHSDerived>& rhs)
        : m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {}

      /**
       * @brief Copy constructor.
       * @param other Division object to copy
       */
      Division(const Division& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

      /**
       * @brief Move constructor.
       * @param other Division object to move from
       */
      Division(Division&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      /**
       * @brief Restricts both functions to the trace of a mesh object.
       * @param args Arguments specifying the trace
       * @returns Reference to this object
       */
      template <class ... Args>
      Division& traceOf(const Args& ... args)
      {
        m_lhs->traceOf(args...);
        m_rhs->traceOf(args...);
        return *this;
      }

      /**
       * @brief Gets the numerator function.
       * @returns Reference to numerator function @f$ f @f$
       */
      constexpr
      const LHSType& getLHS() const
      {
        assert(m_lhs);
        return *m_lhs;
      }

      /**
       * @brief Gets the denominator function.
       * @returns Reference to denominator function @f$ g @f$
       */
      constexpr
      const RHSType& getRHS() const
      {
        assert(m_rhs);
        return *m_rhs;
      }

      /**
       * @brief Evaluates the division at a point.
       * @param p Point at which to evaluate
       * @returns Value @f$ \frac{f(x)}{g(x)} @f$ at point @f$ x @f$
       */
      template <class Point>
      constexpr
      auto getValue(const Point& p) const
      {
        const auto lhs = getLHS().getValue(p);
        const auto rhs = getRHS().getValue(p);
        return lhs / rhs;
      }

      Optional<size_t> getOrder(const Geometry::Polytope& polytope) const noexcept
      {
        const auto lo = getLHS().getOrder(polytope);
        const auto ro = getRHS().getOrder(polytope);
        // Only polynomial if denominator is polynomial of order 0 (constant)
        if (!lo || !ro || *ro != 0)
          return std::nullopt;
        return lo;
      }

      /**
       * @brief Creates a copy of the division operation.
       * @returns Pointer to copied object
       */
      Division* copy() const noexcept final override
      {
        return new Division(*this);
      }

    private:
      std::unique_ptr<FunctionBase<LHSDerived>> m_lhs;
      std::unique_ptr<FunctionBase<RHSDerived>> m_rhs;
  };
  /**
   * @brief Deduction guide for Division.
   */
  template <class LHSDerived, class RHSDerived>
  Division(const FunctionBase<LHSDerived>&, const FunctionBase<RHSDerived>&)
    -> Division<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>;

  /**
   * @brief Division operator for two functions.
   * @param lhs Numerator function
   * @param rhs Denominator function
   * @returns Division object representing @f$ \frac{f}{g} @f$
   */
  template <class LHSDerived, class RHSDerived>
  auto
  operator/(const FunctionBase<LHSDerived>& lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return Division(lhs, rhs);
  }

  /**
   * @brief Division of a function by a number.
   * @param lhs Function to divide
   * @param rhs Number to divide by
   * @returns Division object representing @f$ \frac{f}{c} @f$
   */
  template <class LHSDerived, class Number,
    typename = std::enable_if_t<
      std::is_arithmetic_v<Number>, Division<LHSDerived, RealFunction<Number>>>>
  auto
  operator/(const FunctionBase<LHSDerived>& lhs, Number rhs)
  {
    return Division(lhs, RealFunction(rhs));
  }

  /**
   * @brief Division of a number by a function.
   * @param lhs Number in numerator
   * @param rhs Function in denominator
   * @returns Division object representing @f$ \frac{c}{g} @f$
   */
  template <class Number, class RHSDerived,
    typename = std::enable_if_t<
      std::is_arithmetic_v<Number>, Division<RHSDerived, RealFunction<Number>>>>
  auto
  operator/(Number lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return Division(RealFunction(lhs), rhs);
  }

  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  class Division<ShapeFunctionBase<LHSDerived, FES, Space>, FunctionBase<RHSDerived>> final
    : public ShapeFunctionBase<
        Division<ShapeFunctionBase<LHSDerived, FES, Space>, FunctionBase<RHSDerived>>,
        FES,
        Space>
  {
    public:
      /// @brief Finite element space type.
      using FESType = FES;
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      /// @brief Left-hand side operand type.
      using LHSType = ShapeFunctionBase<LHSDerived, FES, Space>;
      /// @brief Right-hand side operand type.
      using RHSType = FunctionBase<RHSDerived>;

      /// @brief Parent class type.
      using Parent =
        ShapeFunctionBase<Division<LHSType, RHSType>, FES, Space>;

      Division(const LHSType& lhs, const RHSType& rhs)
        : Parent(lhs.getFiniteElementSpace()),
          m_lhs(lhs.copy()),
          m_rhs(rhs.copy())
      {}

      Division(const Division& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()),
          m_rhs(other.m_rhs->copy())
      {}

      Division(Division&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)),
          m_rhs(std::move(other.m_rhs))
      {}

      const auto& getLeaf() const
      {
        return getLHS().getLeaf();
      }

      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getLHS().getDOFs(element);
      }

      const auto& getFiniteElementSpace() const
      {
        return getLHS().getFiniteElementSpace();
      }

      const LHSType& getLHS() const
      {
        assert(m_lhs);
        return *m_lhs;
      }

      const RHSType& getRHS() const
      {
        assert(m_rhs);
        return *m_rhs;
      }

      const IntegrationPoint& getIntegrationPoint() const
      {
        return m_lhs->getIntegrationPoint();
      }

      Division& setIntegrationPoint(const IntegrationPoint& ip)
      {
        m_lhs->setIntegrationPoint(ip);
        return *this;
      }

      auto getBasis(size_t local) const
      {
        const auto& ip = getIntegrationPoint();
        const auto lhs = getLHS().getBasis(local);
        const auto rhs = getRHS().getValue(ip);
        return lhs / rhs;
      }

      Optional<size_t> getOrder(const Geometry::Polytope& polytope) const
      {
        const auto lo = getLHS().getOrder(polytope);
        const auto ro = getRHS().getOrder(polytope);

        if (!lo || !ro || *ro != 0)
          return std::nullopt;

        return lo;
      }

      Division* copy() const noexcept override
      {
        return new Division(*this);
      }

    private:
      std::unique_ptr<LHSType> m_lhs;
      std::unique_ptr<RHSType> m_rhs;
  };

  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  Division(
      const ShapeFunctionBase<LHSDerived, FES, Space>&,
      const FunctionBase<RHSDerived>&)
    -> Division<ShapeFunctionBase<LHSDerived, FES, Space>, FunctionBase<RHSDerived>>;

  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  auto operator/(
      const ShapeFunctionBase<LHSDerived, FES, Space>& lhs,
      const FunctionBase<RHSDerived>& rhs)
  {
    return Division(lhs, rhs);
  }
}
#endif
