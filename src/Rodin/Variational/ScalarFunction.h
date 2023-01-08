/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_SCALARFUNCTION_H
#define RODIN_VARIATIONAL_SCALARFUNCTION_H

#include <map>
#include <set>
#include <memory>
#include <optional>
#include <type_traits>

#include <mfem.hpp>

#include "Rodin/Geometry/Element.h"
#include "ForwardDecls.h"
#include "Function.h"
#include "RangeShape.h"
#include "Exceptions.h"

namespace Rodin::Variational
{
   /**
    * @defgroup ScalarFunctionSpecializations ScalarFunction Template Specializations
    * @brief Template specializations of the ScalarFunction class.
    * @see ScalarFunction
    */

   class ScalarFunctionBase : public FunctionBase
   {
      public:
         ScalarFunctionBase() = default;

         ScalarFunctionBase(const ScalarFunctionBase& other)
            : FunctionBase(other)
         {}

         ScalarFunctionBase(ScalarFunctionBase&& other)
            : FunctionBase(std::move(other))
         {}

         virtual ~ScalarFunctionBase() = default;

         FunctionValue::Scalar operator()(const Geometry::Point& p) const
         {
            return getValue(p).scalar();
         }

         RangeShape getRangeShape() const override
         {
            return {1, 1};
         }

         RangeType getRangeType() const override
         {
            return RangeType::Scalar;
         }

         virtual ScalarFunctionBase* copy() const noexcept override = 0;
   };

   /**
    * @ingroup ScalarFunctionSpecializations
    */
   template <>
   class ScalarFunction<FunctionBase> : public ScalarFunctionBase
   {
      public:
         ScalarFunction(const FunctionBase& nested)
            : m_nested(nested.copy())
         {
            if (nested.getRangeType() != RangeType::Scalar)
               UnexpectedRangeTypeException(RangeType::Scalar, nested.getRangeType()).raise();
         }

         ScalarFunction(const ScalarFunction& other)
            :  ScalarFunctionBase(other),
               m_nested(other.m_nested->copy())
         {}

         ScalarFunction(ScalarFunction&& other)
            : ScalarFunctionBase(std::move(other)),
              m_nested(std::move(other.m_nested))
         {}

         ScalarFunction& traceOf(Geometry::Attribute attrs) override
         {
            ScalarFunctionBase::traceOf(attrs);
            m_nested->traceOf(attrs);
            return *this;
         }

         FunctionValue getValue(const Geometry::Point& v) const override
         {
            return m_nested->getValue(v);
         }

         ScalarFunction* copy() const noexcept override
         {
            return new ScalarFunction(*this);
         }

      private:
         std::unique_ptr<FunctionBase> m_nested;
   };
   ScalarFunction(const FunctionBase&) -> ScalarFunction<FunctionBase>;

   /**
    * @ingroup ScalarFunctionSpecializations
    * @brief Represents a ScalarFunction of arithmetic type `Number`.
    *
    * @tparam Number Arithmetic type
    * @see [std::is_arithmetic](https://en.cppreference.com/w/cpp/types/is_arithmetic)
    */
   template <class Number>
   class ScalarFunction<Number, std::enable_if_t<std::is_arithmetic_v<Number>>>
      : public ScalarFunctionBase
   {
      public:
         /**
          * @brief Constructs a ScalarFunction from an arithmetic value.
          * @param[in] x Arithmetic value
          */
         constexpr
         ScalarFunction(Number x)
            : m_x(x)
         {}

         constexpr
         ScalarFunction(const ScalarFunction& other)
            : ScalarFunctionBase(other),
              m_x(other.m_x)
         {}

         constexpr
         ScalarFunction(ScalarFunction&& other)
            : ScalarFunctionBase(std::move(other)),
              m_x(other.m_x)
         {}

         FunctionValue getValue(const Geometry::Point& p) const override
         {
            return FunctionValue(static_cast<FunctionValue::Scalar>(m_x));
         }

         ScalarFunction* copy() const noexcept override
         {
            return new ScalarFunction(*this);
         }

         Internal::MFEMFunction build(const Geometry::MeshBase&) const override
         {
            return Internal::MFEMFunction(new mfem::ConstantCoefficient(m_x));
         }

      private:
         const Number m_x;
   };
   template <class Number>
   ScalarFunction(Number)
      -> ScalarFunction<Number, std::enable_if_t<std::is_arithmetic_v<Number>>>;

   /**
    * @ingroup ScalarFunctionSpecializations
    * @brief Represents a scalar function given by an arbitrary scalar function.
    */
   template <>
   class ScalarFunction<std::function<double(const Geometry::Point&)>>
      : public ScalarFunctionBase
   {
      public:
         template <class T>
         ScalarFunction(T&& f)
            : ScalarFunction(
                  std::function<double(const Geometry::Point&)>(std::forward<T>(f)))
         {}

         /**
          * @brief Constructs a ScalarFunction from an std::function.
          */
         ScalarFunction(std::function<double(const Geometry::Point&)> f)
            : m_f(f)
         {}

         ScalarFunction(const ScalarFunction& other)
            : ScalarFunctionBase(other),
              m_f(other.m_f)
         {}

         ScalarFunction(ScalarFunction&& other)
            : ScalarFunctionBase(std::move(other)),
              m_f(std::move(other.m_f))
         {}

         FunctionValue getValue(const Geometry::Point& v) const override
         {
            return m_f(v);
         }

         ScalarFunction* copy() const noexcept override
         {
            return new ScalarFunction(*this);
         }

      private:
         const std::function<double(const Geometry::Point&)> m_f;
   };

   ScalarFunction(std::function<double(const Geometry::Point&)>)
      -> ScalarFunction<std::function<double(const Geometry::Point&)>>;

   template <class T>
   ScalarFunction(T)
      -> ScalarFunction<
         std::enable_if_t<std::is_invocable_r_v<double, T, const Geometry::Point&>,
         std::function<double(const Geometry::Point&)>>>;
}

#endif
