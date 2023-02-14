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

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/DenseMatrix.h"

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Simplex.h"

#include "Rodin/FormLanguage/Base.h"
#include "Rodin/FormLanguage/Traits.h"
#include "Rodin/Utility/Overloaded.h"

#include "ForwardDecls.h"

#include "RangeShape.h"

namespace Rodin::Variational::Internal
{
  template <class Derived>
  class ScalarProxyFunction;

  template <class Derived>
  class VectorProxyFunction;

  template <class Derived>
  class MatrixProxyFunction;

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
    using Scalar  = double; ///< Scalar value type
    using Boolean = bool; ///< Scalar value type
    using Vector  = Math::Vector; ///< Vector value type
    using Matrix  = Math::Matrix; ///< Matrix value type

    public:
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

      template <int Size>
      constexpr
      FunctionValue(Eigen::Vector<Scalar, Size>&& v)
        : m_v(Math::Vector(std::move(v)))
      {}

      template <int Rows, int Cols>
      constexpr
      FunctionValue(Eigen::Matrix<Scalar, Rows, Cols>&& v)
        : m_v(Math::Matrix(std::move(v)))
      {}

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
            [&](Boolean& v) { v = v && static_cast<Boolean>(s); }} , m_v);
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
            [&](Boolean& v) { v = v && static_cast<Boolean>(1.0 / s); }} , m_v);
        return *this;
      }

    private:
      std::variant<Scalar, Boolean, Vector, Matrix> m_v;
  };

  template <class Derived>
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

      virtual ~FunctionBase() = default;

      FunctionBase& operator=(FunctionBase&& other)
      {
        m_traceDomain = other.m_traceDomain;
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
      Geometry::Attribute getTraceDomain() const
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
       * @note It is not necessary to set the size beforehand.
       */
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).getValue(p);
      }

      /**
       * @brief Evaluates the function on a vertex of the mesh.
       * @param[in] v Vertex belonging to the mesh
       */
      constexpr
      auto operator()(const Geometry::Point& v) const
      {
        return getValue(v);
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
      virtual FunctionBase& traceOf(Geometry::Attribute attr)
      {
        m_traceDomain = attr;
        return *this;
      }

      virtual FunctionBase* copy() const noexcept override = 0;

      virtual Internal::MFEMFunction build(const Geometry::MeshBase& mesh) const;

    private:
      Geometry::Attribute m_traceDomain;
  };
}

namespace Rodin::FormLanguage
{
  template <class Derived>
  struct Traits<Variational::FunctionBase<Derived>>
  {
    using ResultType = std::result_of_t<
      decltype(&Variational::FunctionBase<Derived>::getValue())(const Geometry::Point&)>;
  };
}

namespace Rodin::Variational::Internal
{
  template <class Derived>
  class ScalarProxyFunction : public mfem::Coefficient
  {
    public:
      ScalarProxyFunction(const Geometry::MeshBase& mesh, const FunctionBase<Derived>& s)
         : m_mesh(mesh), m_s(s)
      {}

      ScalarProxyFunction(const ScalarProxyFunction& other)
        : mfem::Coefficient(other),
          m_mesh(other.m_mesh),
          m_s(other.m_s)
      {}

      ScalarProxyFunction(ScalarProxyFunction&& other)
        : mfem::Coefficient(std::move(other)),
          m_mesh(other.m_mesh),
          m_s(other.m_s)
      {}

      inline
      double Eval(mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip)
      final override
      {
        Scalar res = 0;
        switch (trans.ElementType)
        {
          case mfem::ElementTransformation::ELEMENT:
          {
            res = m_s.getValue(Geometry::Point(*m_mesh.get().getElement(trans.ElementNo), ip));
            break;
          }
          case mfem::ElementTransformation::BDR_ELEMENT:
          {
            int faceIdx = m_mesh.get().getHandle().GetBdrFace(trans.ElementNo);
            res = m_s.getValue(Geometry::Point(*m_mesh.get().getFace(faceIdx), ip));
            break;
          }
          case mfem::ElementTransformation::FACE:
          {
            res = m_s.getValue(Geometry::Point(*m_mesh.get().getFace(trans.ElementNo), ip));
            break;
          }
        }
        return res;
      }

    private:
      std::reference_wrapper<const Geometry::MeshBase> m_mesh;
      std::reference_wrapper<const FunctionBase<Derived>> m_s;
  };

  template <class Derived>
  class VectorProxyFunction : public mfem::VectorCoefficient
  {
    public:
      VectorProxyFunction(const Geometry::MeshBase& mesh, const FunctionBase<Derived>& s)
        : mfem::VectorCoefficient(
            s.getRangeShape().height() == 1 ?
              s.getRangeShape().width() : s.getRangeShape().height()),
          m_mesh(mesh),
          m_s(s)
      {}

      VectorProxyFunction(const VectorProxyFunction& other)
        : mfem::VectorCoefficient(other),
          m_mesh(other.m_mesh),
          m_s(other.m_s)
      {}

      VectorProxyFunction(VectorProxyFunction&& other)
        : mfem::VectorCoefficient(std::move(other)),
          m_mesh(other.mesh),
          m_s(other.m_s)
      {}

      inline
      void Eval(
          mfem::Vector& value, mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip)
      final override
      {
        Math::Vector vec;
        switch (trans.ElementType)
        {
          case mfem::ElementTransformation::ELEMENT:
          {
            vec = m_s.getValue(Geometry::Point(*m_mesh.get().getElement(trans.ElementNo), ip));
            break;
          }
          case mfem::ElementTransformation::BDR_ELEMENT:
          {
            int faceIdx = m_mesh.get().getHandle().GetBdrFace(trans.ElementNo);
            vec = m_s.getValue(Geometry::Point(*m_mesh.get().getFace(faceIdx), ip));
            break;
          }
          case mfem::ElementTransformation::FACE:
          {
            vec = m_s.getValue(Geometry::Point(*m_mesh.get().getFace(trans.ElementNo), ip));
            break;
          }
        }
        value.SetSize(vec.size());
        std::copy(vec.begin(), vec.end(), value.begin());
      }

    private:
      std::reference_wrapper<const Geometry::MeshBase> m_mesh;
      std::reference_wrapper<const FunctionBase<Derived>> m_s;
  };

  template <class Derived>
  class MatrixProxyFunction : public mfem::MatrixCoefficient
  {
    public:
      MatrixProxyFunction(const Geometry::MeshBase& mesh, const FunctionBase<Derived>& s)
        : mfem::MatrixCoefficient(s.getRangeShape().height(), s.getRangeShape().width()),
          m_mesh(mesh),
          m_s(s)
      {}

      MatrixProxyFunction(const MatrixProxyFunction& other)
        : mfem::MatrixCoefficient(other),
          m_mesh(other.m_mesh),
          m_s(other.m_s)
      {}

      MatrixProxyFunction(MatrixProxyFunction&& other)
        : mfem::MatrixCoefficient(std::move(other)),
          m_mesh(other.m_mesh),
          m_s(other.m_s)
      {}

      inline
      void Eval(
          mfem::DenseMatrix& value, mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip)
      final override
      {
        Math::Matrix mat;
        switch (trans.ElementType)
        {
          case mfem::ElementTransformation::ELEMENT:
          {
            mat = m_s.getValue(Geometry::Point(*m_mesh.get().getElement(trans.ElementNo), ip));
            break;
          }
          case mfem::ElementTransformation::BDR_ELEMENT:
          {
            int faceIdx = m_mesh.get().getHandle().GetBdrFace(trans.ElementNo);
            mat = m_s.getValue(Geometry::Point(*m_mesh.get().getFace(faceIdx), ip));
            break;
          }
          case mfem::ElementTransformation::FACE:
          {
            mat = m_s.getValue(Geometry::Point(*m_mesh.get().getFace(trans.ElementNo), ip));
            break;
          }
        }
        value.SetSize(mat.rows(), mat.cols());
        value = mat.data();
      }
    private:
      std::reference_wrapper<const Geometry::MeshBase> m_mesh;
      std::reference_wrapper<const FunctionBase<Derived>> m_s;
  };
}

#endif
