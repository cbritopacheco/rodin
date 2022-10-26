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

#include "Rodin/Mesh/Element.h"
#include "ForwardDecls.h"
#include "Function.h"
#include "RangeShape.h"

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

   /**
    * @ingroup ScalarFunctionSpecializations
    */
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
    * @ingroup ScalarFunctionSpecializations
    * @brief Represents a ScalarFunction of arithmetic type `Number`.
    *
    * @tparam Number Arithmetic type
    * @see [std::is_arithmetic](https://en.cppreference.com/w/cpp/types/is_arithmetic)
    */
   template <class Number>
   class ScalarFunction<Number> : public ScalarFunctionBase
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

         Internal::MFEMFunction build() const override
         {
            return Internal::MFEMFunction(new mfem::ConstantCoefficient(m_x));
         }

      private:
         const Number m_x;
   };
   template <class T>
   ScalarFunction(const T&)
      -> ScalarFunction<std::enable_if_t<std::is_arithmetic_v<T>, T>>;

   /**
    * @ingroup ScalarFunctionSpecializations
    * @brief Represents a scalar function given by an arbitrary scalar function.
    */
   template <>
   class ScalarFunction<std::function<double(const Vertex&)>>
      : public ScalarFunctionBase
   {
      public:
         template <class T>
         ScalarFunction(T&& f)
            : ScalarFunction(
                  std::function<double(const Vertex&)>(std::forward<T>(f)))
         {}

         /**
          * @brief Constructs a ScalarFunction from an std::function.
          */
         ScalarFunction(std::function<double(const Vertex&)> f)
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
            mfem::Vector transip;
            trans.Transform(ip, transip);
            Vertex v(std::move(transip));
            v.setElementTransformation(&trans).setIntegrationPoint(&ip);
            return m_f(v);
         }

         double operator()(const Vertex& v) const
         {
            return getValue(*v.getElementTransformation(), *v.getIntegrationPoint());
         }

         ScalarFunction* copy() const noexcept override
         {
            return new ScalarFunction(*this);
         }

      private:
         const std::function<double(const Vertex&)> m_f;
   };

   ScalarFunction(std::function<double(const Vertex&)>)
      -> ScalarFunction<std::function<double(const Vertex&)>>;

   template <class T>
   ScalarFunction(T)
      -> ScalarFunction<
         std::enable_if_t<std::is_invocable_r_v<double, T, const Vertex&>,
         std::function<double(const Vertex&)>>>;
}

#endif
