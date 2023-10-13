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
      template <typename U, typename = decltype(std::declval<U>().myMethodRef(std::declval<Args&>()...))>
      static std::true_type Test(int);

      template <typename U>
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

      inline
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
      inline
      constexpr
      const TraceDomain& getTraceDomain() const
      {
        return m_traceDomain;
      }

      inline
      constexpr
      Transpose<FunctionBase> T() const
      {
        return Transpose<FunctionBase>(*this);
      }

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        return static_cast<const Derived&>(*this).getRangeShape();
      }

      inline
      constexpr
      RangeType getRangeType() const
      {
        using R = typename FormLanguage::Traits<FunctionBase<Derived>>::RangeType;
        if constexpr (std::is_same_v<R, Boolean>)
        {
          return RangeType::Boolean;
        }
        else if constexpr (std::is_same_v<R, Integer>)
        {
          return RangeType::Integer;
        }
        else if constexpr (std::is_same_v<R, Scalar>)
        {
          return RangeType::Scalar;
        }
        else if constexpr (std::is_same_v<R, Math::Vector>)
        {
          return RangeType::Vector;
        }
        else if constexpr (std::is_same_v<R, Math::Matrix>)
        {
          return RangeType::Matrix;
        }
        else
        {
          assert(false);
          static_assert(Utility::DependentFalse<R>::Value);
        }
      }

      /**
       * @brief Evaluates the function on a Point belonging to the mesh.
       * @note CRTP function to be overriden in Derived class.
       */
      inline
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).getValue(p);
      }

      inline
      constexpr
      void getValue(Math::Vector& res, const Geometry::Point& p) const
      {
        if constexpr (Internal::HasGetValueMethod<Derived, Math::Vector&>::Value)
        {
          return static_cast<const Derived&>(*this).getValue(res, p);
        }
        else
        {
          res = getValue(p);
        }
      }

      inline
      constexpr
      void getValue(Math::Matrix& res, const Geometry::Point& p) const
      {
        if constexpr (Internal::HasGetValueMethod<Derived, Math::Matrix&>::Value)
        {
          return static_cast<const Derived&>(*this).getValue(res, p);
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
      inline
      constexpr
      auto operator()(const Geometry::Point& p) const
      {
        return getValue(p);
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
      inline
      constexpr
      FunctionBase& traceOf(Geometry::Attribute attr)
      {
        m_traceDomain = FlatSet<Geometry::Attribute>{attr};
        return *this;
      }

      inline
      constexpr
      FunctionBase& traceOf(const FlatSet<Geometry::Attribute>& attr)
      {
        m_traceDomain = attr;
        return *this;
      }

      virtual FunctionBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      FlatSet<Geometry::Attribute> m_traceDomain;
  };
}

#include "Function.hpp"

#endif
