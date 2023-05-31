/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_VECTORFUNCTION_H
#define RODIN_VARIATIONAL_VECTORFUNCTION_H

#include <memory>
#include <optional>
#include <type_traits>

#include <mfem.hpp>

#include "ForwardDecls.h"

#include "Rodin/Alert.h"
#include "Rodin/Utility/ForConstexpr.h"

#include "Function.h"
#include "ScalarFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup VectorFunctionSpecializations VectorFunction Template Specializations
   * @brief Template specializations of the VectorFunction class.
   * @see VectorFunction
   */

  template <class Derived>
  class VectorFunctionBase : public FunctionBase<VectorFunctionBase<Derived>>
  {
    public:
      using Parent = FunctionBase<VectorFunctionBase<Derived>>;

      VectorFunctionBase() = default;

      VectorFunctionBase(const VectorFunctionBase& other)
        : Parent(other)
      {}

      VectorFunctionBase(VectorFunctionBase&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Convenience function to access the 1st component of the
       * vector.
       */
      inline
      constexpr
      auto x() const
      {
        assert(getDimension() >= 1);
        return operator()(0);
      }

      /**
       * @brief Convenience function to access the 2nd component of the
       * vector.
       */
      inline
      constexpr
      auto y() const
      {
        assert(getDimension() >= 2);
        return operator()(1);
      }

      inline
      constexpr
      auto z() const
      {
        assert(getDimension() >= 3);
        return operator()(2);
      }

      inline
      constexpr
      auto operator()(const Geometry::Point& p) const
      {
        return getValue(p);
      }

      /**
       * @brief Access the ith component of the vector function.
       * @returns Object of type Component<VectorFunctionBase> representing
       * the ith component of the VectorFunction.
       */
      inline
      constexpr
      auto operator()(size_t i) const
      {
        assert(0 <= i);
        assert(i < getDimension());
        return Component(*this, i);
      }

      virtual ~VectorFunctionBase() = default;

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        return { getDimension(), 1 };
      }

      inline
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).getValue(p);
      }

      /**
       * @brief Gets the dimension of the vector object.
       * @returns Dimension of vector.
       */
      inline
      constexpr
      size_t getDimension() const
      {
        return static_cast<const Derived&>(*this).getDimension();
      }

      virtual VectorFunctionBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }
  };

  /**
   * @ingroup VectorFunctionSpecializations
   * @tparam V Type of first value
   * @tparam Values Parameter pack of remaining values
   * @brief Represents a vector function which may be constructed from values
   * which can be converted to objects of type ScalarFunction.
   *
   * In general one may construct any VectorFunction by specifying its values
   * in a uniform initialization manner. For example, to construct a
   * VectorFunction with constant entries (1, 2, 3) :
   * @code{.cpp}
   * auto v = VectorFunction{1, 2, 3};
   * @endcode
   * Alternatively, we may construct instances of VectorFunction from any type
   * which is convertible to specializations of ScalarFunction:
   * @code{.cpp}
   * auto s = ScalarFunction(3.1416);
   * auto v = VectorFunction{Dx(s), 42, s};
   * @endcode
   */
  template <class V, class ... Values>
  class VectorFunction<V, Values...> final
    : public VectorFunctionBase<VectorFunction<V, Values...>>
  {
    public:
      using Parent = VectorFunctionBase<VectorFunction<V, Values...>>;
      /**
       * @brief Constructs a vector with the given values.
       * @param[in] values Parameter pack of values
       *
       * Each value passed must be convertible to any specialization of
       * ScalarFunction.
       */
      VectorFunction(const V v, const Values... values)
        : m_fs(ScalarFunction(v), ScalarFunction(values)...)
      {}

      VectorFunction(const VectorFunction& other)
        : Parent(other),
          m_fs(other.m_fs)
      {}

      VectorFunction(VectorFunction&& other)
        : Parent(std::move(other)),
          m_fs(std::move(other.m_fs))
      {}

      inline
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        Math::FixedSizeVector<1 + sizeof...(Values)> value;
        Utility::ForIndex<1 + sizeof...(Values)>(
            [&](auto i){ value(static_cast<Eigen::Index>(i)) = std::get<i>(m_fs).getValue(p); });
        return value;
      }

      inline
      constexpr
      size_t getDimension() const
      {
        return 1 + sizeof...(Values);
      }

      inline
      constexpr
      VectorFunction& traceOf(Geometry::Attribute attrs)
      {
        std::apply([&](auto& s) { s.traceOf(attrs); }, m_fs);
        return *this;
      }

      inline VectorFunction* copy() const noexcept override
      {
        return new VectorFunction(*this);
      }

    private:
      std::tuple<ScalarFunction<V>, ScalarFunction<Values>...> m_fs;
  };

  template <class V, class ... Values>
  VectorFunction(V, Values...) -> VectorFunction<V, Values...>;
}

#endif
