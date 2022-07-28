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

#include "ForwardDecls.h"
#include "Function.h"
#include "RangeShape.h"

namespace Rodin::Variational
{
   /**
    * @brief Abstract base class for objects representing scalar functions.
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

         void getValue(
               mfem::DenseMatrix& value,
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip) const override
         {
            value.SetSize(1, 1);
            value(0, 0) = getValue(trans, ip);
         }

         RangeShape getRangeShape() const override
         {
            return {1, 1};
         }

         RangeType getRangeType() const override
         {
            return RangeType::Scalar;
         }

         /**
          * @brief Computes the value at the given transformation and
          * integration point.
          * @returns Value at given transformation and integration point.
          */
         virtual double getValue(
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const = 0;

         virtual ScalarFunctionBase* copy() const noexcept override = 0;
   };

   template <>
   class ScalarFunction<FunctionBase> : public ScalarFunctionBase
   {
      public:
         ScalarFunction(const FunctionBase& nested);

         ScalarFunction(const ScalarFunction& other);

         ScalarFunction(ScalarFunction&& other);

         ScalarFunction& traceOf(const std::set<int>& attrs) override
         {
            ScalarFunctionBase::traceOf(attrs);
            m_nested->traceOf(attrs);
            return *this;
         }

         double getValue(
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip) const override;

         ScalarFunction* copy() const noexcept override
         {
            return new ScalarFunction(*this);
         }

      private:
         std::unique_ptr<FunctionBase> m_nested;
   };
   ScalarFunction(const FunctionBase&) -> ScalarFunction<FunctionBase>;

   /**
    * @brief Represents a ScalarFunction of arithmetic type `T`.
    *
    * @tparam T Arithmetic type
    * @see [std::is_arithmetic](https://en.cppreference.com/w/cpp/types/is_arithmetic)
    */
   template <class Number>
   class ScalarFunction : public ScalarFunctionBase
   {
      public:
         static_assert(std::is_arithmetic_v<Number>, "T must be an arithmetic type");

         /**
          * @brief Constructs a ScalarFunction from an arithmetic value.
          * @param[in] x Arithmetic value
          */
         template <typename U = Number>
         constexpr
         ScalarFunction(typename std::enable_if_t<std::is_arithmetic_v<U>, U> x)
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

         double getValue(
               mfem::ElementTransformation&,
               const mfem::IntegrationPoint&) const override
         {
            return m_x;
         }

         ScalarFunction* copy() const noexcept override
         {
            return new ScalarFunction(*this);
         }

      private:
         const Number m_x;
   };
   template <class T>
   ScalarFunction(const T&)
      -> ScalarFunction<std::enable_if_t<std::is_arithmetic_v<T>, T>>;

   /**
    * @brief Represents a scalar function given by an arbitrary function.
    *
    * This class represents a ScalarFunction whose values are given by any
    * function taking `double*` data array and `int` dimension parameter.
    */
   template <>
   class ScalarFunction<std::function<double(const double*, int)>>
      : public ScalarFunctionBase
   {
      public:
         template <class T>
         ScalarFunction(T&& f)
            : ScalarFunction(
                  std::function<double(const double*, int)>(std::forward<T>(f)))
         {}

         /**
          * @brief Constructs a ScalarFunction from an std::function.
          */
         ScalarFunction(std::function<double(const double*, int)> f)
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

         double getValue(mfem::ElementTransformation& trans, const
               mfem::IntegrationPoint& ip) const override
         {
            double x[3];
            mfem::Vector transip(x, 3);
            trans.Transform(ip, transip);
            return m_f(transip.GetData(), transip.Size());
         }

         ScalarFunction* copy() const noexcept override
         {
            return new ScalarFunction(*this);
         }

      private:
         const std::function<double(const double*, int)> m_f;
   };
   ScalarFunction(std::function<double(const double*, int)>)
      -> ScalarFunction<std::function<double(const double*, int)>>;
   template <class T>
   ScalarFunction(T)
      -> ScalarFunction<
      std::enable_if_t<std::is_invocable_r_v<double, T, const double*, int>,
      std::function<double(const double*, int)>>>;
}

#endif
