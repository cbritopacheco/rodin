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
      : public FormLanguage::Buildable<mfem::VectorCoefficient>
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
          * respective "interior" domain. For example, coefficients which
          * involve the derivatives of a GridFunction need to know the element
          * to "trace".
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
          * @brief Convenience function to access the 3nd component of the
          * vector.
          */
         Component<VectorFunctionBase> z() const;

         virtual ~VectorFunctionBase() = default;

         virtual Component<VectorFunctionBase> operator()(int i) const;

         virtual void getValue(
               mfem::Vector& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const = 0;

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
    * @brief Variadic vector of values
    * Represents a vector which may be constructed from values which can be
    * converted to objects of type ScalarFunction.
    */
   template <class ... Values>
   class VectorFunction : public VectorFunctionBase
   {
      public:
         /**
          * @brief Constructs a vector with the given values.
          * @param[in] values Parameter pack of values
          */
         constexpr
         VectorFunction(Values... values)
         {
            m_coeffs.reserve(sizeof...(Values));
            makeCoefficientsFromTuple(std::forward_as_tuple(values...));
         }

         constexpr
         VectorFunction(const VectorFunction& other)
         {
            m_coeffs.reserve(sizeof...(Values));
            for (const auto& v : other.m_coeffs)
               m_coeffs.emplace_back(v->copy());
         }

         constexpr
         VectorFunction(VectorFunction&&) = default;

         void getValue(
               mfem::Vector& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const override
         {
            value.SetSize(static_cast<int>(sizeof...(Values)));
            for (int i = 0; i < sizeof...(Values); i++)
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

   /**
    * @brief Vector which can be constructed from a GridFunction with vector
    * values.
    *
    * Represents a vector which may be constructed from a GridFunction which
    * takes on vector values at the mesh vertices.
    */
   template <class FEC, class Trait>
   class VectorFunction<GridFunction<FEC, Trait>>
      : public VectorFunctionBase
   {
      public:
         /**
          * @brief Constructs a VectorFunction from a vector valued GridFunction.
          */
         constexpr
         VectorFunction(GridFunction<FEC, Trait>& u)
            :  m_dimension(u.getFiniteElementSpace().getVectorDimension()),
               m_u(u),
               m_mfemVectorFunction(&u.getHandle())
         {}

         constexpr
         VectorFunction(const VectorFunction& other)
            :  m_dimension(other.m_dimension),
               m_u(other.m_u)
         {}

         int getDimension() const override
         {
            return m_dimension;
         }

         void getValue(
               mfem::Vector& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const override
         {
            mfem::Mesh* gfMesh = m_u.getHandle().FESpace()->GetMesh();
            if (trans.mesh == gfMesh)
            {
               m_u.getHandle().GetVectorValue(trans, ip, value);
            }
            else
            {
               mfem::IntegrationPoint coarseIp;
               mfem::ElementTransformation* coarseTrans = refinedToCoarse(*gfMesh, trans, ip, coarseIp);
               m_u.getHandle().GetVectorValue(*coarseTrans, coarseIp, value);
            }
         }

         VectorFunction* copy() const noexcept override
         {
            return new VectorFunction(*this);
         }

      private:
         const size_t m_dimension;
         GridFunction<FEC, Trait>& m_u;

         mfem::VectorGridFunctionCoefficient m_mfemVectorFunction;
   };
   template <class FEC, class Trait>
   VectorFunction(GridFunction<FEC, Trait>&)
      -> VectorFunction<GridFunction<FEC, Trait>>;
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
