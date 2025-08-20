/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FUNCTION_H
#define RODIN_VARIATIONAL_FUNCTION_H

#include "Rodin/Cast.h"

#include "Rodin/Geometry/Point.h"

#include "Rodin/Variational/Traits.h"

#include "Rodin/FormLanguage/Base.h"
#include "Rodin/FormLanguage/Traits.h"

#include "ForwardDecls.h"

namespace Rodin::FormLanguage
{
  template <class Derived>
  struct Traits<Variational::FunctionBase<Derived>>
  {
    using ResultType = typename ResultOf<Variational::FunctionBase<Derived>>::Type;

    using RangeType = typename RangeOf<Variational::FunctionBase<Derived>>::Type;

    using ScalarType = typename FormLanguage::Traits<RangeType>::ScalarType;
  };
}

namespace Rodin::Variational
{
  /**
   * @brief Base class for functions defined on a mesh.
   */
  template <class Derived>
  class FunctionBase : public FormLanguage::Base
  {
    public:
      using Parent = FormLanguage::Base;

      using TraceDomain = FlatSet<Geometry::Attribute>;

      FunctionBase() = default;

      FunctionBase(const FunctionBase& other)
        : Parent(other),
          m_traceDomain(other.m_traceDomain)
      {}

      FunctionBase(FunctionBase&& other)
        : Parent(std::move(other)),
          m_traceDomain(std::move(other.m_traceDomain))
      {}

      virtual ~FunctionBase() = default;

      FunctionBase& operator=(FunctionBase&& other)
      {
        m_traceDomain = std::move(other.m_traceDomain);
        return *this;
      }

      /**
       * @brief Evaluates the function on a Point belonging to the mesh.
       *
       * This calls the function get getValue(const Geometry::Point&).
       */
      constexpr
      decltype(auto) operator()(const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).getValue(p);
      }

      auto operator()(size_t i) const
      {
        return Component(*this, i);
      }

      auto operator()(size_t i, size_t j) const
      {
        return Component(*this, i, j);
      }

      constexpr
      auto x() const
      {
        return Component(*this, 0);
      }

      constexpr
      auto y() const
      {
        return Component(*this, 1);
      }

      constexpr
      auto z() const
      {
        return Component(*this, 2);
      }

      constexpr
      auto T() const
      {
        return Transpose(*this);
      }

      /**
       * @brief Sets an attribute which will be interpreted as the domain to
       * trace.
       *
       * Convenience function to call traceOf(FlatSet<int>) with only one
       * attribute.
       *
       * @returns Reference to self (for method chaining)
       */
      constexpr
      Derived& traceOf(const Geometry::Attribute& attr)
      {
        return static_cast<Derived&>(*this).traceOf(FlatSet<Geometry::Attribute>{ attr });
      }

      template <class A1, class A2, class ... As>
      constexpr
      Derived& traceOf(const A1& a1, const A2& a2, const As& ... as)
      {
        return static_cast<Derived&>(*this).traceOf(FlatSet<Geometry::Attribute>{ a1, a2, as... });
      }

      constexpr
      Derived& traceOf(const FlatSet<Geometry::Attribute>& attr)
      {
        m_traceDomain = attr;
        return static_cast<Derived&>(*this);
      }

      /**
       * @brief Gets the set of attributes which will be interpreted as the
       * domains to "trace".
       *
       * The domains to trace are interpreted as the domains where there
       * shall be a continuous extension from values to the interior
       * boundaries. If the trace domain is empty, then this has the
       * semantic value that it has not been specified yet.
       */
      constexpr
      const TraceDomain& getTraceDomain() const
      {
        return m_traceDomain;
      }

      /**
       * @brief Evaluates the function on a Point belonging to the mesh.
       * @note CRTP function to be overriden in Derived class.
       */
      constexpr
      decltype(auto) getValue(const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).getValue(p);
      }

      const Derived& getDerived() const noexcept
      {
        return static_cast<const Derived&>(*this);
      }

      virtual FunctionBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      FlatSet<Geometry::Attribute> m_traceDomain;
  };
}

#endif
