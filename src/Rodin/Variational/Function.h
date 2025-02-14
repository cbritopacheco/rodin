/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FUNCTION_H
#define RODIN_VARIATIONAL_FUNCTION_H

#include <set>
#include <variant>
#include <type_traits>

#include "Rodin/Cast.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/Matrix.h"

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Polytope.h"

#include "Rodin/FormLanguage/Base.h"
#include "Rodin/FormLanguage/Traits.h"
#include "Rodin/Utility/Overloaded.h"

#include "ForwardDecls.h"

#include "RangeType.h"
#include "RangeShape.h"

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
  namespace Internal
  {
    template <typename T, class ... Args>
    struct HasGetValueMethod
    {
      template<typename U, typename = decltype(std::declval<U>().getValue(std::declval<Args>()...))>
      static std::true_type Test(int);

      template<typename U>
      static std::false_type Test(...);

      using Type = decltype(Test<T>(0));
      static constexpr bool Value = Type::value;
    };

    template <typename T, typename... Args>
    struct HasGetValueMethod<T, Args&...>
    {
      template <typename U, typename = decltype(std::declval<U>().getValue(std::declval<Args&>()...))>
      static std::true_type Test(int);

      template <typename U>
      static std::false_type Test(...);

      using Type = decltype(Test<T>(0));
      static constexpr bool Value = Type::value;
    };

    template <typename T, class... Args>
    struct HasGetValueMethodR
    {
        template<typename U, typename = decltype(std::declval<U>().getValue(std::declval<Args>()...))>
        static auto Test(int) ->
          decltype(std::is_same<typename std::invoke_result<decltype(&U::getValue)(U, Args...)>::type, T>::value, std::true_type{});

        template<typename U>
        static std::false_type Test(...);

        using Type = decltype(Test<T>(0));
        static constexpr bool Value = Type::value;
    };
  }

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

      constexpr
      Transpose<FunctionBase> T() const
      {
        return Transpose<FunctionBase>(*this);
      }

      constexpr
      RangeShape getRangeShape() const
      {
        return static_cast<const Derived&>(*this).getRangeShape();
      }

      /**
       * @brief Evaluates the function on a Point belonging to the mesh.
       * @note CRTP function to be overriden in Derived class.
       */
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).getValue(p);
      }

      template <class T>
      constexpr
      void getValue(T& res, const Geometry::Point& p) const
      {
        if constexpr (Internal::HasGetValueMethod<Derived, T&, const Geometry::Point&>::Value)
        {
          static_cast<const Derived&>(*this).getValue(res, p);
        }
        else
        {
          res = getValue(p);
        }
      }

      /**
       * @brief Evaluates the function on a Point belonging to the mesh.
       *
       * This calls the function get getValue(const Geometry::Point&).
       */
      constexpr
      auto operator()(const Geometry::Point& p) const
      {
        return getValue(p);
      }

      template <class T>
      constexpr
      void operator()(T& res, const Geometry::Point& p) const
      {
        getValue(res, p);
      }

      auto coeff(size_t i, size_t j) const
      {
        return Component(*this, i, j);
      }

      auto coeff(size_t i) const
      {
        return Component(*this, i);
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
      Derived& traceOf(Geometry::Attribute attr)
      {
        return traceOf(FlatSet<Geometry::Attribute>{ attr });
      }

      template <class A1, class A2, class ... As>
      constexpr
      Derived& traceOf(A1 a1, A2 a2, As ... as)
      {
        return traceOf(FlatSet<Geometry::Attribute>{ a1, a2, as... });
      }

      constexpr
      Derived& traceOf(const FlatSet<Geometry::Attribute>& attr)
      {
        m_traceDomain = attr;
        return static_cast<Derived&>(*this);
      }

      Derived& getDerived()
      {
        return static_cast<Derived&>(*this);
      }

      const Derived& getDerived() const
      {
        return static_cast<const Derived&>(*this);
      }

      template <class ToRange>
      constexpr
      auto cast() const
      {
        return Cast<FunctionBase, ToRange>(*this);
      }

      virtual FunctionBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      FlatSet<Geometry::Attribute> m_traceDomain;
  };
}

namespace Rodin
{
  template <class FromDerived, class ToRange>
  class Cast<Variational::FunctionBase<FromDerived>, ToRange> final
    : public Variational::FunctionBase<Cast<Variational::FunctionBase<FromDerived>, ToRange>>
  {
    public:
      using FromType = Variational::FunctionBase<FromDerived>;

      using FromRangeType = typename FormLanguage::Traits<FromType>::RangeType;

      using ToRangeType = ToRange;

      using TraceDomain = FlatSet<Geometry::Attribute>;

      using Parent = Variational::FunctionBase<Cast<Variational::FunctionBase<FromDerived>, ToRange>>;


      Cast(const FromType& from)
        : m_from(from.copy())
      {}

      Cast(const Cast& other)
        : Parent(other),
          m_from(other.m_from->copy())
      {}

      Cast(Cast&& other)
        : Parent(std::move(other)),
          m_from(std::move(other.m_from))
      {}

      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return static_cast<ToRangeType>(m_from->getValue(p));
      }

      Cast* copy() const noexcept override
      {
        return new Cast(*this);
      }

    private:
      std::unique_ptr<FromType> m_from;
  };
}

#include "Function.hpp"

#endif
