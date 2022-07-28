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

#include "Utility.h"
#include "Function.h"
#include "ScalarFunction.h"

namespace Rodin::Variational
{
   /**
    * @brief Abstract base class for objects representing vector coefficients.
    *
    * @note Vectors are zero indexed. This means that the 0-index corresponds
    * to the 1st entry of the vector.
    */
   class VectorFunctionBase : public FunctionBase
   {
      public:
         VectorFunctionBase() = default;

         VectorFunctionBase(const VectorFunctionBase& other)
            : FunctionBase(other)
         {}

         VectorFunctionBase(VectorFunctionBase&& other)
            : FunctionBase(std::move(other))
         {}

         /**
          * @brief Convenience function to access the 1st component of the
          * vector.
          */
         Component<FunctionBase> x() const;

         /**
          * @brief Convenience function to access the 2nd component of the
          * vector.
          */
         Component<FunctionBase> y() const;

         /**
          * @brief Convenience function to access the 3rd component of the
          * vector.
          */
         Component<FunctionBase> z() const;

         virtual ~VectorFunctionBase() = default;

         void getValue(
               mfem::DenseMatrix& value,
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip) const override
         {
            mfem::Vector v;
            getValue(v, trans, ip);
            value.UseExternalData(v.StealData(), getDimension(), 1);
            value.GetMemory().SetHostPtrOwner(true);
         }

         RangeShape getRangeShape() const override
         {
            return {getDimension(), 1};
         }

         RangeType getRangeType() const override
         {
            return RangeType::Vector;
         }

         /**
          * @brief Access the ith component of the vector function.
          * @returns Object of type Component<VectorFunctionBase> representing
          * the ith component of the VectorFunction.
          */
         virtual Component<FunctionBase> operator()(int i) const;

         /**
          * @brief Computes the value at the given transformation and
          * integration point.
          * @returns Value at given transformation and integration point.
          */
         virtual void getValue(
               mfem::Vector& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const = 0;

         /**
          * @brief Gets the dimension of the vector object.
          * @returns Dimension of vector.
          */
         virtual int getDimension() const = 0;

         virtual VectorFunctionBase* copy() const noexcept override = 0;
   };

   /**
    * @brief Represents a vector which may be constructed from values which can
    * be converted to objects of type ScalarFunction.
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
   template <class ... Values>
   class VectorFunction : public VectorFunctionBase
   {
      public:
         /**
          * @brief Constructs a vector with the given values.
          * @param[in] values Parameter pack of values
          *
          * Each value passed must be convertible to any specialization of
          * ScalarFunction.
          */
         constexpr
         VectorFunction(Values... values)
         {
            m_fs.reserve(sizeof...(Values));
            makeFsFromTuple(std::forward_as_tuple(values...));
         }

         constexpr
         VectorFunction(const VectorFunction& other)
            : VectorFunctionBase(other)
         {
            m_fs.reserve(sizeof...(Values));
            for (const auto& v : other.m_fs)
               m_fs.emplace_back(v->copy());
         }

         constexpr
         VectorFunction(VectorFunction&& other)
            :  VectorFunctionBase(std::move(other)),
               m_fs(std::move(other.m_fs))
         {}

         void getValue(
               mfem::Vector& value,
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip) const override
         {
            value.SetSize(static_cast<int>(sizeof...(Values)));
            for (size_t i = 0; i < sizeof...(Values); i++)
               value(i) = m_fs[i]->getValue(trans, ip);
         }

         int getDimension() const override
         {
            return sizeof...(Values);
         }

         VectorFunction& traceOf(const std::set<int>& attrs) override
         {
            VectorFunctionBase::traceOf(attrs);
            for (auto& f : m_fs)
               f->traceOf(attrs);
            return *this;
         }

         VectorFunction* copy() const noexcept override
         {
            return new VectorFunction(*this);
         }

      private:
         template<std::size_t I = 0, class ... Tp>
         typename std::enable_if_t<I == sizeof...(Tp)>
         makeFsFromTuple(const std::tuple<Tp...>&)
         {}

         template<std::size_t I = 0, class ... Tp>
         typename std::enable_if_t<I < sizeof...(Tp)>
         makeFsFromTuple(const std::tuple<Tp...>& t)
         {
            m_fs.emplace_back(new ScalarFunction(std::get<I>(t)));
            makeFsFromTuple<I + 1, Tp...>(t);
         }

         std::vector<std::unique_ptr<ScalarFunctionBase>> m_fs;
   };
   template <class ... Values>
   VectorFunction(Values&&...) -> VectorFunction<Values...>;
}

#endif
