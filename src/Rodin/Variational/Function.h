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

#include <mfem.hpp>

#include "Rodin/FormLanguage/Base.h"
#include "Rodin/Utility/Overloaded.h"

#include "ForwardDecls.h"
#include "Rodin/Geometry/Element.h"

namespace Rodin::Variational::Internal
{
   template <RangeType R>
   class ProxyFunction;

   using ScalarProxyFunction = ProxyFunction<RangeType::Scalar>;
   using VectorProxyFunction = ProxyFunction<RangeType::Vector>;
   using MatrixProxyFunction = ProxyFunction<RangeType::Matrix>;

   class MFEMFunction
   {
      public:
         constexpr
         MFEMFunction(nullptr_t)
         {}

         MFEMFunction(mfem::Coefficient* p)
            :  m_range(RangeType::Scalar),
               m_v(std::unique_ptr<mfem::Coefficient>(p))
         {}

         MFEMFunction(mfem::VectorCoefficient* p)
            :  m_range(RangeType::Vector),
               m_v(std::unique_ptr<mfem::VectorCoefficient>(p))
         {}

         MFEMFunction(mfem::MatrixCoefficient* p)
            :  m_range(RangeType::Matrix),
               m_v(std::unique_ptr<mfem::MatrixCoefficient>(p))
         {}

         MFEMFunction(MFEMFunction&& other)
            :  m_range(std::move(other.m_range)),
               m_v(std::move(other.m_v))
         {}

         constexpr
         RangeType getRangeType() const
         {
            assert(m_range);
            return *m_range;
         }

         template <RangeType R>
         constexpr
         std::enable_if_t<R == RangeType::Scalar, mfem::Coefficient&>
         get()
         {
            assert(m_range);
            assert(*m_range == R);
            assert(std::get<std::unique_ptr<mfem::Coefficient>>(m_v));
            return *std::get<std::unique_ptr<mfem::Coefficient>>(m_v);
         }

         template <RangeType R>
         constexpr
         std::enable_if_t<R == RangeType::Scalar, const mfem::Coefficient&>
         get() const
         {
            assert(m_range);
            assert(*m_range == R);
            assert(std::get<std::unique_ptr<mfem::Coefficient>>(m_v));
            return *std::get<std::unique_ptr<mfem::Coefficient>>(m_v);
         }

         template <RangeType R>
         constexpr
         std::enable_if_t<R == RangeType::Vector, mfem::VectorCoefficient&>
         get()
         {
            assert(m_range);
            assert(*m_range == R);
            assert(std::get<std::unique_ptr<mfem::VectorCoefficient>>(m_v));
            return *std::get<std::unique_ptr<mfem::VectorCoefficient>>(m_v);
         }

         template <RangeType R>
         constexpr
         std::enable_if_t<R == RangeType::Vector, const mfem::VectorCoefficient&>
         get() const
         {
            assert(m_range);
            assert(*m_range == R);
            assert(std::get<std::unique_ptr<mfem::VectorCoefficient>>(m_v));
            return *std::get<std::unique_ptr<mfem::VectorCoefficient>>(m_v);
         }

         template <RangeType R>
         constexpr
         std::enable_if_t<R == RangeType::Matrix, mfem::MatrixCoefficient&>
         get()
         {
            assert(m_range);
            assert(*m_range == R);
            assert(std::get<std::unique_ptr<mfem::MatrixCoefficient>>(m_v));
            return *std::get<std::unique_ptr<mfem::MatrixCoefficient>>(m_v);
         }

         template <RangeType R>
         constexpr
         std::enable_if_t<R == RangeType::Matrix, const mfem::MatrixCoefficient&>
         get() const
         {
            assert(m_range);
            assert(*m_range == R);
            assert(std::get<std::unique_ptr<mfem::MatrixCoefficient>>(m_v));
            return *std::get<std::unique_ptr<mfem::MatrixCoefficient>>(m_v);
         }

      private:
         std::optional<RangeType> m_range;
         std::variant<
            std::unique_ptr<mfem::Coefficient>,
            std::unique_ptr<mfem::VectorCoefficient>,
            std::unique_ptr<mfem::MatrixCoefficient>> m_v;
   };
}

namespace Rodin::Variational
{
   /**
    * @brief Represents the value of a FunctionBase object when evaluated
    * on a Mesh Vertex.
    *
    * In Rodin, valuations of instances of FunctionBase on mesh vertices
    * can take the form three possible value types:
    * - Scalar
    * - Vector
    * - Matrix
    *
    * The type resolution of the actual type is performed at runtime.
    */
   class FunctionValue
   {
      public:
         using Scalar   = double; ///< Scalar value type
         using Boolean  = bool; ///< Scalar value type
         using Vector   = mfem::Vector; ///< Vector value type
         using Matrix   = mfem::DenseMatrix; ///< Matrix value type

         constexpr
         FunctionValue() = delete;

         constexpr
         FunctionValue(Scalar v)
            : m_v(v)
         {}

         constexpr
         FunctionValue(Boolean v)
            : m_v(v)
         {}

         constexpr
         FunctionValue(Vector&& v)
         {
            m_v.emplace<Vector>();
            std::get<Vector>(m_v).Swap(v);
         }

         constexpr
         FunctionValue(Matrix&& v)
         {
            m_v.emplace<mfem::DenseMatrix>();
            std::get<mfem::DenseMatrix>(m_v).Swap(v);
         }

         FunctionValue(const FunctionValue&) = default;

         FunctionValue(FunctionValue&&) = default;

         /**
          * @brief Queries which type the FunctionValue %holds.
          */
         template <class T>
         inline
         constexpr
         bool holds() const
         {
            return std::holds_alternative<T>(m_v);
         }

         inline
         constexpr
         operator Boolean() const
         {
            assert(holds<Boolean>());
            return std::get<Boolean>(m_v);
         }

         inline
         constexpr
         operator Scalar() const
         {
            assert(holds<Scalar>());
            return std::get<Scalar>(m_v);
         }

         inline
         constexpr
         operator Vector&() &
         {
            assert(holds<Vector>());
            return std::get<Vector>(m_v);
         }

         inline
         constexpr
         operator Vector&&() &&
         {
            assert(holds<Vector>());
            return std::move(std::get<Vector>(m_v));
         }

         inline
         constexpr
         operator Matrix&() &
         {
            assert(holds<Matrix>());
            return std::get<Matrix>(m_v);
         }

         inline
         constexpr
         operator Matrix&&() &&
         {
            assert(holds<Matrix>());
            return std::move(std::get<Matrix>(m_v));
         }

         inline
         constexpr
         Boolean boolean() const
         {
            assert(holds<Boolean>());
            return std::get<Boolean>(m_v);
         }

         inline
         constexpr
         Scalar scalar() const
         {
            assert(holds<Scalar>());
            return std::get<Scalar>(m_v);
         }

         inline
         constexpr
         Vector& vector() &
         {
            assert(holds<Vector>());
            return std::get<Vector>(m_v);
         }

         inline
         constexpr
         const Vector& vector() const &
         {
            assert(holds<Vector>());
            return std::get<Vector>(m_v);
         }

         inline
         constexpr
         Vector&& vector() &&
         {
            assert(holds<Vector>());
            return std::move(std::get<Vector>(m_v));
         }

         inline
         constexpr
         Matrix& matrix() &
         {
            assert(holds<Matrix>());
            return std::get<Matrix>(m_v);
         }

         inline
         constexpr
         const Matrix& matrix() const &
         {
            assert(holds<Matrix>());
            return std::get<Matrix>(m_v);
         }

         inline
         constexpr
         Matrix&& matrix() &&
         {
            assert(holds<Matrix>());
            return std::move(std::get<Matrix>(m_v));
         }

         inline
         constexpr
         FunctionValue& operator-()
         {
            *this *= -1.0;
            return *this;
         }

         inline
         constexpr
         FunctionValue& operator+=(const FunctionValue& s)
         {
            std::visit(Utility::Overloaded{
                  [&](Scalar& v) { v += s.scalar(); },
                  [&](Vector& v) { v += s.vector(); },
                  [&](Matrix& v) { v += s.matrix(); },
                  [&](Boolean& v) { v += s.boolean(); }} , m_v);
            return *this;
         }

         inline
         constexpr
         FunctionValue& operator*=(Scalar s)
         {
            std::visit(Utility::Overloaded{
                  [&](Scalar& v) { v *= s; },
                  [&](Vector& v) { v *= s; },
                  [&](Matrix& v) { v *= s; },
                  [&](Boolean& v) { v *= s; }} , m_v);
            return *this;
         }

         inline
         constexpr
         FunctionValue& operator/=(Scalar s)
         {
            std::visit(Utility::Overloaded{
                  [&](Scalar& v) { v *= (1.0 / s); },
                  [&](Vector& v) { v *= (1.0 / s); },
                  [&](Matrix& v) { v *= (1.0 / s); },
                  [&](Boolean& v) { v *= (1.0 / s); }} , m_v);
            return *this;
         }

      private:
         std::variant<Scalar, Boolean, Vector, Matrix> m_v;
   };

   class FunctionBase : public FormLanguage::Base
   {
      public:
         FunctionBase() = default;

         FunctionBase(const FunctionBase& other)
            : FormLanguage::Base(other),
              m_traceDomain(other.m_traceDomain)
         {}

         FunctionBase(FunctionBase&& other)
            : FormLanguage::Base(std::move(other)),
              m_traceDomain(std::move(other.m_traceDomain))
         {}

         FunctionBase& operator=(FunctionBase&& other)
         {
            m_traceDomain = std::move(other.m_traceDomain);
            return *this;
         }

         /**
          * @brief Sets an attribute which will be interpreted as the domain to
          * trace.
          *
          * Convenience function to call traceOf(std::set<int>) with only one
          * attribute.
          *
          * @returns Reference to self (for method chaining)
          */
         virtual FunctionBase& traceOf(int attr)
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
          * @note Setting the trace domain of a FunctionBase instance
          * does not guarantee that it will taken into consideration when
          * computing its value. That said, it is up to the subclass to decide
          * how it will use this information which can be obtained via the
          * getTraceDomain() method.
          *
          * @see @ref getTraceDomain() "getTraceDomain()"
          */
         virtual FunctionBase& traceOf(const std::set<int>& attrs)
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

         virtual Transpose<FunctionBase> T() const;

         virtual RangeShape getRangeShape() const = 0;

         virtual RangeType getRangeType() const;

         /**
          * @note It is not necessary to set the size beforehand.
          */
         virtual FunctionValue getValue(const Geometry::Point& p) const = 0;

         /**
          * @brief Evaluates the function on a vertex of the mesh.
          * @param[in] v Vertex belonging to the mesh
          */
         FunctionValue operator()(const Geometry::Point& v) const
         {
            return getValue(v);
         }

         virtual FunctionBase* copy() const noexcept override = 0;

         virtual Internal::MFEMFunction build(const Geometry::MeshBase& mesh) const;

      private:
         std::set<int> m_traceDomain;
   };
}

namespace Rodin::Variational::Internal
{
   template <>
   class ProxyFunction<RangeType::Scalar> : public mfem::Coefficient
   {
      public:
         ProxyFunction(const Geometry::MeshBase& mesh, const FunctionBase& s);

         ProxyFunction(const ProxyFunction& other);

         ProxyFunction(ProxyFunction&& other);

         double Eval(
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip) override;

      private:
         const Geometry::MeshBase& m_mesh;
         const FunctionBase& m_s;
   };

   template <>
   class ProxyFunction<RangeType::Vector> : public mfem::VectorCoefficient
   {
      public:
         ProxyFunction(const Geometry::MeshBase& mesh, const FunctionBase& s);

         ProxyFunction(const ProxyFunction& other);

         ProxyFunction(ProxyFunction&& other);

         void Eval(
               mfem::Vector& value,
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip) override;

      private:
         const Geometry::MeshBase& m_mesh;
         const FunctionBase& m_s;
   };

   template <>
   class ProxyFunction<RangeType::Matrix> : public mfem::MatrixCoefficient
   {
      public:
         ProxyFunction(const Geometry::MeshBase& mesh, const FunctionBase& s);

         ProxyFunction(const ProxyFunction& other);

         ProxyFunction(ProxyFunction&& other);

         void Eval(
               mfem::DenseMatrix& value,
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip) override;
      private:
         const Geometry::MeshBase& m_mesh;
         const FunctionBase& m_s;
   };
}

#endif
