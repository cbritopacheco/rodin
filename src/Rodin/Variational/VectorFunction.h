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
#include "FormLanguage/Base.h"

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
   class VectorFunctionBase
      :  public FunctionBase,
         public FormLanguage::Buildable<mfem::VectorCoefficient>
   {
      public:
         VectorFunctionBase() = default;

         VectorFunctionBase(const VectorFunctionBase&) = default;

         /**
          * @brief Sets an attribute which will be interpreted as the domain to
          * trace.
          *
          * Convenience function to call traceOf(std::set<int>) with only one
          * attribute.
          *
          * @returns Reference to self (for method chaining)
          */
         VectorFunctionBase& traceOf(int attr)
         {
            return traceOf(std::set<int>{attr});
         }

         /**
          * @brief Sets which attributes will be interpreted as the domain to
          * trace.
          * @returns Reference to self (for method chaining)
          *
          * When integrating along interior boundaries sometimes it is
          * necessary to specify which attributes should be interpreted as the
          * respective "interior" domain, since it is not clear which domain
          * attribute can be used to extend the value continuously up to the
          * boundary. To resolve this ambiguity the trace domain is interpreted
          * as the domain which shall be used to make this continuous
          * extension.
          *
          * @note Setting the trace domain of a VectorFunctionBase instance
          * does not guarantee that it will taken into consideration when
          * computing its value. That said, it is up to the subclass to decide
          * how it will use this information which can be obtained via the
          * getTraceDomain() method.
          *
          * @see @ref VectorFunctionBase::getTraceDomain() "getTraceDomain()"
          *
          */
         VectorFunctionBase& traceOf(std::set<int> attrs)
         {
            m_traceDomain = attrs;
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
         const std::set<int>& getTraceDomain() const
         {
            return m_traceDomain;
         }

         /**
          * @brief Convenience function to access the 1st component of the
          * vector.
          */
         Component<VectorFunctionBase> x() const;

         /**
          * @brief Convenience function to access the 2nd component of the
          * vector.
          */
         Component<VectorFunctionBase> y() const;

         /**
          * @brief Convenience function to access the 3rd component of the
          * vector.
          */
         Component<VectorFunctionBase> z() const;

         virtual ~VectorFunctionBase() = default;

         /**
          * @brief Access the ith component of the vector function.
          * @returns Object of type Component<VectorFunctionBase> representing
          * the ith component of the VectorFunction.
          */
         virtual Component<VectorFunctionBase> operator()(int i) const;

         /**
          * @brief Computes the value at the given transformation and
          * integration point.
          * @returns Value at given transformation and integration point.
          */
         virtual void getValue(
               mfem::Vector& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const = 0;

         void getValue(
               mfem::DenseMatrix& value,
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip) const override
         {
            mfem::Vector v;
            getValue(v, trans, ip);
            value.SetSize(getDimension(), 1);
            value.GetFromVector(0, v);
         }

         std::tuple<int, int> getRangeShape() const override
         {
            return {getDimension(), 1};
         }

         RangeType getRangeType() const override
         {
            return RangeType::Vector;
         }

         /**
          * @brief Gets the dimension of the vector object.
          * @returns Dimension of vector.
          */
         virtual int getDimension() const = 0;

         virtual VectorFunctionBase* copy() const noexcept override = 0;

         std::unique_ptr<mfem::VectorCoefficient> build() const override;

      private:
         std::set<int> m_traceDomain;
   };

   /**
    * @brief Represents from a GridFunction with vector
    * values.
    *
    * Represents a vector which may be constructed from a GridFunction which
    * takes on vector values at the mesh vertices.
    */
   template <>
   class VectorFunction<GridFunctionBase>
      : public VectorFunctionBase
   {
      public:
         /**
          * @brief Constructs a VectorFunction from a vector valued GridFunction.
          */
         VectorFunction(const GridFunctionBase& u);

         VectorFunction(const VectorFunction& other)
            :  VectorFunctionBase(other),
               m_dimension(other.m_dimension),
               m_u(other.m_u)
         {}

         int getDimension() const override
         {
            return m_dimension;
         }

         void getValue(
               mfem::Vector& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const override;

         VectorFunction* copy() const noexcept override
         {
            return new VectorFunction(*this);
         }

      private:
         const size_t m_dimension;
         const GridFunctionBase& m_u;
   };
   VectorFunction(GridFunctionBase&) -> VectorFunction<GridFunctionBase>;
   VectorFunction(const GridFunctionBase&) -> VectorFunction<GridFunctionBase>;

   template <class FEC, class Trait>
   VectorFunction(GridFunction<FEC, Trait>&) -> VectorFunction<GridFunctionBase>;
   template <class FEC, class Trait>
   VectorFunction(const GridFunction<FEC, Trait>&) -> VectorFunction<GridFunctionBase>;

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
            m_coeffs.reserve(sizeof...(Values));
            makeCoefficientsFromTuple(std::forward_as_tuple(values...));
         }

         constexpr
         VectorFunction(const VectorFunction& other)
            : VectorFunctionBase(other)
         {
            m_coeffs.reserve(sizeof...(Values));
            for (const auto& v : other.m_coeffs)
               m_coeffs.emplace_back(v->copy());
         }

         constexpr
         VectorFunction(VectorFunction&& other)
            : VectorFunctionBase(std::move(other)),
               m_coeffs(std::move(other.m_coeffs))
         {}

         void getValue(
               mfem::Vector& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const override
         {
            value.SetSize(static_cast<int>(sizeof...(Values)));
            for (size_t i = 0; i < sizeof...(Values); i++)
               value(i) = m_coeffs[i]->getValue(trans, ip);
         }

         int getDimension() const override
         {
            return sizeof...(Values);
         }

         VectorFunction* copy() const noexcept override
         {
            return new VectorFunction(*this);
         }

      private:
         template<std::size_t I = 0, class ... Tp>
         typename std::enable_if_t<I == sizeof...(Tp)>
         makeCoefficientsFromTuple(const std::tuple<Tp...>&)
         {}

         template<std::size_t I = 0, class ... Tp>
         typename std::enable_if_t<I < sizeof...(Tp)>
         makeCoefficientsFromTuple(const std::tuple<Tp...>& t)
         {
            m_coeffs.emplace_back(new ScalarFunction(std::get<I>(t)));
            makeCoefficientsFromTuple<I + 1, Tp...>(t);
         }

         std::vector<std::unique_ptr<ScalarFunctionBase>> m_coeffs;
   };
   template <class ... Values>
   VectorFunction(Values&&...) -> VectorFunction<Values...>;
}

namespace Rodin::Variational::Internal
{
   class ProxyVectorFunction : public mfem::VectorCoefficient
   {
      public:
         ProxyVectorFunction(const VectorFunctionBase& v)
            :  mfem::VectorCoefficient(v.getDimension()),
               m_v(v)
         {}

         ProxyVectorFunction(const ProxyVectorFunction& other)
            :  mfem::VectorCoefficient(other),
               m_v(other.m_v)
         {}

         ProxyVectorFunction(ProxyVectorFunction&& other)
            :  mfem::VectorCoefficient(std::move(other)),
               m_v(other.m_v)
         {}

         void Eval(mfem::Vector& value, mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip) override
         {
            m_v.getValue(value, trans, ip);
         }

      private:
         const VectorFunctionBase& m_v;
   };
}

#endif
