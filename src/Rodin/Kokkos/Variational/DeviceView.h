/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_DEVICEVIEW_H
#define RODIN_KOKKOS_VARIATIONAL_DEVICEVIEW_H

/**
 * @file
 * @brief Lowering of form language expressions to device-callable views.
 */

#include <string>
#include <tuple>
#include <type_traits>

#include <Kokkos_Core.hpp>

#include "Rodin/Configure.h"

#include "Rodin/Variational/Function.h"
#include "Rodin/Variational/RealFunction.h"
#include "Rodin/Variational/VectorFunction.h"
#include "Rodin/Variational/GridFunction.h"
#include "Rodin/Variational/ShapeFunction.h"
#include "Rodin/Variational/Sum.h"
#include "Rodin/Variational/Mult.h"
#include "Rodin/Variational/Dot.h"
#include "Rodin/Variational/UnaryMinus.h"
#include "Rodin/Variational/Transpose.h"
#include "Rodin/Variational/P1/Div.h"
#include "Rodin/Variational/P1/Grad.h"
#include "Rodin/Variational/P1/Jacobian.h"
#include "Rodin/Variational/P1/ShapeFunction.h"

#include "Rodin/Kokkos/Geometry/MeshView.h"

#include "P1.h"

namespace Rodin::FormLanguage
{
  /**
   * @brief Claims callables taking a device point for the Kokkos backend.
   *
   * Callables invocable on a host Geometry::Point are left to the host,
   * which is tested first so that generic host lambdas are never
   * instantiated with a device point.
   *
   * @tparam F Callable type.
   */
  template <class F>
  struct IsDevicePointCallable<F,
    std::enable_if_t<
      std::conjunction_v<std::negation<std::is_invocable<F, const Geometry::Point&>>,
        std::is_invocable_r<Real, F, const ::Rodin::Kokkos::Geometry::PointView&>>>>
  {
    /// @brief True, the Kokkos backend evaluating the callable.
    static constexpr bool Value = true;
  };
}

namespace Rodin::Kokkos::Variational
{
  /// @brief P1 space of the given range on a resident mesh with the given backend.
  template <class Range, class Backend>
  using DeviceP1Space = ::Rodin::Variational::P1<Range,
    ::Rodin::Geometry::Mesh<::Rodin::Context::Kokkos<::Rodin::Context::Local, Backend>>>;

  /// @brief Scalar P1 space on a resident mesh with the given backend.
  template <class Backend>
  using DeviceP1 = DeviceP1Space<Real, Backend>;

  /// @brief Vector P1 space on a resident mesh with the given backend.
  template <class Backend>
  using DeviceVectorP1 = DeviceP1Space<Math::SpatialVector<Real>, Backend>;

  /**
   * @brief Multiplies two scalar values.
   * @param a First factor.
   * @param b Second factor.
   * @return Product @f$ a b @f$.
   */
  KOKKOS_INLINE_FUNCTION
  Real mult(Real a, Real b)
  {
    return a * b;
  }

  /**
   * @brief Scales a device vector.
   * @param s Scaling factor.
   * @param v Device vector.
   * @return Vector @f$ s v @f$.
   */
  KOKKOS_INLINE_FUNCTION
  Geometry::Vector mult(Real s, const Geometry::Vector& v)
  {
    Geometry::Vector result(v.size());
    for (Index i = 0; i < v.size(); ++i)
      result(i) = s * v(i);
    return result;
  }

  /**
   * @brief Scales a device vector.
   * @param v Device vector.
   * @param s Scaling factor.
   * @return Vector @f$ v s @f$.
   */
  KOKKOS_INLINE_FUNCTION
  Geometry::Vector mult(const Geometry::Vector& v, Real s)
  {
    return mult(s, v);
  }

  /**
   * @brief Scales a device matrix.
   * @param s Scaling factor.
   * @param m Device matrix.
   * @return Matrix @f$ s m @f$.
   */
  KOKKOS_INLINE_FUNCTION
  Geometry::Matrix mult(Real s, const Geometry::Matrix& m)
  {
    Geometry::Matrix result(m.rows(), m.cols());
    for (Index i = 0; i < m.rows(); ++i)
    {
      for (Index j = 0; j < m.cols(); ++j)
        result(i, j) = s * m(i, j);
    }
    return result;
  }

  /**
   * @brief Scales a device matrix.
   * @param m Device matrix.
   * @param s Scaling factor.
   * @return Matrix @f$ m s @f$.
   */
  KOKKOS_INLINE_FUNCTION
  Geometry::Matrix mult(const Geometry::Matrix& m, Real s)
  {
    return mult(s, m);
  }

  /**
   * @brief Multiplies a device vector by a device matrix.
   * @param m Device matrix.
   * @param v Device vector, with as many entries as @p m has columns.
   * @return Vector @f$ m v @f$.
   */
  KOKKOS_INLINE_FUNCTION
  Geometry::Vector mult(const Geometry::Matrix& m, const Geometry::Vector& v)
  {
    KOKKOS_ASSERT(m.cols() == v.size());
    Geometry::Vector result(m.rows());
    for (Index i = 0; i < m.rows(); ++i)
    {
      for (Index j = 0; j < m.cols(); ++j)
        result(i) += m(i, j) * v(j);
    }
    return result;
  }

  /**
   * @brief Multiplies two device matrices.
   * @param a First matrix.
   * @param b Second matrix, with as many rows as @p a has columns.
   * @return Matrix @f$ a b @f$.
   */
  KOKKOS_INLINE_FUNCTION
  Geometry::Matrix mult(const Geometry::Matrix& a, const Geometry::Matrix& b)
  {
    KOKKOS_ASSERT(a.cols() == b.rows());
    Geometry::Matrix result(a.rows(), b.cols());
    for (Index i = 0; i < a.rows(); ++i)
    {
      for (Index j = 0; j < b.cols(); ++j)
      {
        for (Index k = 0; k < a.cols(); ++k)
          result(i, j) += a(i, k) * b(k, j);
      }
    }
    return result;
  }

  /**
   * @brief Adds two scalar values.
   * @param a First value.
   * @param b Second value.
   * @return Sum @f$ a + b @f$.
   */
  KOKKOS_INLINE_FUNCTION
  Real add(Real a, Real b)
  {
    return a + b;
  }

  /**
   * @brief Adds two device vectors.
   * @param a First vector.
   * @param b Second vector, with the dimension of @p a.
   * @return Sum @f$ a + b @f$.
   */
  KOKKOS_INLINE_FUNCTION
  Geometry::Vector add(const Geometry::Vector& a, const Geometry::Vector& b)
  {
    KOKKOS_ASSERT(a.size() == b.size());
    Geometry::Vector result(a.size());
    for (Index i = 0; i < a.size(); ++i)
      result(i) = a(i) + b(i);
    return result;
  }

  /**
   * @brief Adds two device matrices.
   * @param a First matrix.
   * @param b Second matrix, with the shape of @p a.
   * @return Sum @f$ a + b @f$.
   */
  KOKKOS_INLINE_FUNCTION
  Geometry::Matrix add(const Geometry::Matrix& a, const Geometry::Matrix& b)
  {
    KOKKOS_ASSERT(a.rows() == b.rows() && a.cols() == b.cols());
    Geometry::Matrix result(a.rows(), a.cols());
    for (Index i = 0; i < a.rows(); ++i)
    {
      for (Index j = 0; j < a.cols(); ++j)
        result(i, j) = a(i, j) + b(i, j);
    }
    return result;
  }

  /**
   * @brief Transposes a device matrix.
   * @param m Device matrix.
   * @return Transpose @f$ m^T @f$.
   */
  KOKKOS_INLINE_FUNCTION
  Geometry::Matrix transpose(const Geometry::Matrix& m)
  {
    Geometry::Matrix result(m.cols(), m.rows());
    for (Index i = 0; i < m.rows(); ++i)
    {
      for (Index j = 0; j < m.cols(); ++j)
        result(j, i) = m(i, j);
    }
    return result;
  }

  /**
   * @brief Multiplies two scalar values.
   * @param a First value.
   * @param b Second value.
   * @return Product @f$ a b @f$.
   */
  KOKKOS_INLINE_FUNCTION
  Real dot(Real a, Real b)
  {
    return a * b;
  }

  /**
   * @brief Computes the Euclidean inner product of two device vectors.
   * @param a First vector.
   * @param b Second vector, with the dimension of @p a.
   * @return Inner product @f$ a \cdot b @f$.
   */
  KOKKOS_INLINE_FUNCTION
  Real dot(const Geometry::Vector& a, const Geometry::Vector& b)
  {
    KOKKOS_ASSERT(a.size() == b.size());
    Real result = 0;
    for (Index i = 0; i < a.size(); ++i)
      result += a(i) * b(i);
    return result;
  }

  /**
   * @brief Computes the Frobenius inner product of two device matrices.
   * @param a First matrix.
   * @param b Second matrix, with the shape of @p a.
   * @return Inner product @f$ a : b @f$.
   */
  KOKKOS_INLINE_FUNCTION
  Real dot(const Geometry::Matrix& a, const Geometry::Matrix& b)
  {
    KOKKOS_ASSERT(a.rows() == b.rows() && a.cols() == b.cols());
    Real result = 0;
    for (Index i = 0; i < a.rows(); ++i)
    {
      for (Index j = 0; j < a.cols(); ++j)
        result += a(i, j) * b(i, j);
    }
    return result;
  }

  /**
   * @brief P1 basis values and physical gradients at one quadrature point.
   *
   * The point is evaluated once per quadrature point and shared by every
   * shape function of the integrand. Gradients are mapped to physical space
   * by @f$ \nabla \phi = J^{-T} \hat{\nabla} \hat{\phi} @f$, as on the host.
   */
  class IntegrationPointView
  {
    public:
      /// @brief Maximum number of local basis functions.
      static constexpr Index Capacity = RODIN_MAXIMUM_POLYTOPE_VERTICES;

      /**
       * @brief Evaluates the P1 basis of the point's polytope.
       * @param point Point on a mesh entity.
       */
      KOKKOS_INLINE_FUNCTION
      explicit IntegrationPointView(const Geometry::PointView& point)
        : m_point(point),
          m_count(point.getPolytope().getVertices().size())
      {
        KOKKOS_ASSERT(m_count <= Capacity);
        const Geometry::Matrix inverse = point.getJacobianInverse();
        const Index rdim = inverse.rows();
        const Index sdim = inverse.cols();
        for (Index a = 0; a < m_count; ++a)
        {
          Geometry::Vector reference(rdim);
          point.evaluateBasis(a, m_basis[a], reference);
          Geometry::Vector gradient(sdim);
          for (Index c = 0; c < sdim; ++c)
          {
            for (Index k = 0; k < rdim; ++k)
              gradient(c) += inverse(k, c) * reference(k);
          }
          m_gradients[a] = gradient;
        }
      }

      /// @brief Returns the evaluation point. @return Device point.
      KOKKOS_INLINE_FUNCTION
      const Geometry::PointView& getPoint() const
      {
        return m_point;
      }

      /// @brief Returns the number of local basis functions. @return Vertex count.
      KOKKOS_INLINE_FUNCTION
      Index getCount() const
      {
        return m_count;
      }

      /**
       * @brief Returns a basis value.
       * @param local Local vertex index.
       * @return Basis value @f$ \phi_a @f$.
       */
      KOKKOS_INLINE_FUNCTION
      Real getBasis(Index local) const
      {
        KOKKOS_ASSERT(local < m_count);
        return m_basis[local];
      }

      /**
       * @brief Returns a physical basis gradient.
       * @param local Local vertex index.
       * @return Gradient @f$ \nabla \phi_a @f$.
       */
      KOKKOS_INLINE_FUNCTION
      const Geometry::Vector& getGradient(Index local) const
      {
        KOKKOS_ASSERT(local < m_count);
        return m_gradients[local];
      }

    private:
      Geometry::PointView m_point;
      Index m_count;
      ::Kokkos::Array<Real, Capacity> m_basis;
      ::Kokkos::Array<Geometry::Vector, Capacity> m_gradients;
  };

  /**
   * @brief Device-callable view of a form language expression.
   *
   * A view is built on the host from the expression and copied by value into
   * kernels. It holds no host pointer and no virtual function. Function views
   * provide `getValue(const Geometry::PointView&)`, while shape function views
   * provide `getBasis(const IntegrationPointView&, Index)`.
   *
   * The primary template is left undefined so that an expression without a
   * device lowering is rejected at compile time.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | `DeviceView<FunctionBase<Derived>>` | Forwards to the derived function. |
   * | `DeviceView<ScalarFunctionBase<Scalar, Derived>>` | Forwards to the derived function. |
   * | `DeviceView<RealFunctionBase<Derived>>` | Forwards to the derived function. |
   * | `DeviceView<RealFunction<F>>` | Callable taking a device point. |
   * | `DeviceView<RealFunction<Real>>` | Real constant. |
   * | `DeviceView<Sum<FunctionBase<L>, FunctionBase<R>>>` | Sum of functions. |
   * | `DeviceView<Mult<FunctionBase<L>, FunctionBase<R>>>` | Product of functions. |
   * | `DeviceView<UnaryMinus<FunctionBase<N>>>` | Negated function. |
   * | `DeviceView<RealFunction<Integer>>` | Integer constant. |
   * | `DeviceView<VectorFunctionBase<Scalar, Derived>>` | Forwards to the derived function. |
   * | `DeviceView<MatrixFunctionBase<Scalar, Derived>>` | Forwards to the derived function. |
   * | `DeviceView<GradBase<Operand, Derived>>` | Forwards to the derived function. |
   * | `DeviceView<JacobianBase<Operand, Derived>>` | Forwards to the derived function. |
   * | `DeviceView<DivBase<Operand, Derived>>` | Forwards to the derived function. |
   * | `DeviceView<VectorFunction<V, Values...>>` | Vector function given component-wise. |
   * | `DeviceView<GridFunctionBaseReference<GridFunction<DeviceP1Space<R, B>, Data>>>` | Referenced P1 grid function. |
   * | `DeviceView<GridFunction<DeviceP1Space<R, B>, Data>>` | P1 grid function. |
   * | `DeviceView<Grad<GridFunction<DeviceP1<B>, Data>>>` | Physical gradient of a P1 grid function. |
   * | `DeviceView<Jacobian<GridFunction<DeviceVectorP1<B>, Data>>>` | Physical Jacobian of a vector P1 grid function. |
   * | `DeviceView<Div<GridFunction<DeviceVectorP1<B>, Data>>>` | Divergence of a vector P1 grid function. |
   * | `DeviceView<ShapeFunctionBase<Derived, FES, Space>>` | Forwards to the derived shape function. |
   * | `DeviceView<ShapeFunction<Derived, DeviceP1<B>, Space>>` | P1 basis value. |
   * | `DeviceView<Grad<ShapeFunction<Derived, DeviceP1<B>, Space>>>` | Physical P1 gradient. |
   * | `DeviceView<ShapeFunction<Derived, DeviceVectorP1<B>, Space>>` | Vector P1 basis value. |
   * | `DeviceView<Jacobian<ShapeFunction<Derived, DeviceVectorP1<B>, Space>>>` | Physical vector P1 Jacobian. |
   * | `DeviceView<Div<ShapeFunction<Derived, DeviceVectorP1<B>, Space>>>` | Vector P1 divergence. |
   * | `DeviceView<Sum<ShapeFunctionBase<L, FES, Space>, ShapeFunctionBase<R, FES, Space>>>` | Sum of shape functions. |
   * | `DeviceView<UnaryMinus<ShapeFunctionBase<N, FES, Space>>>` | Negated shape function. |
   * | `DeviceView<Transpose<ShapeFunctionBase<N, FES, Space>>>` | Transposed shape function. |
   * | `DeviceView<Mult<FunctionBase<L>, ShapeFunctionBase<R, FES, Space>>>` | Function times shape function. |
   * | `DeviceView<Mult<ShapeFunctionBase<L, FES, Space>, FunctionBase<R>>>` | Shape function times function. |
   * | `DeviceView<Dot<FunctionBase<L>, ShapeFunctionBase<R, FES, Space>>>` | Dot product of a function and a shape function. |
   *
   * @tparam Node Form language expression type.
   */
  template <class Node>
  struct DeviceView;

  /// @brief Forwards a function to its derived type.
  template <class Derived>
  struct DeviceView<::Rodin::Variational::FunctionBase<Derived>> : DeviceView<Derived>
  {
    /// @param node Function expression.
    DeviceView(const ::Rodin::Variational::FunctionBase<Derived>& node)
      : DeviceView<Derived>(static_cast<const Derived&>(node))
    {}
  };

  /// @brief Forwards a scalar function to its derived type.
  template <class Scalar, class Derived>
  struct DeviceView<::Rodin::Variational::ScalarFunctionBase<Scalar, Derived>>
    : DeviceView<Derived>
  {
    /// @param node Function expression.
    DeviceView(const ::Rodin::Variational::ScalarFunctionBase<Scalar, Derived>& node)
      : DeviceView<Derived>(static_cast<const Derived&>(node))
    {}
  };

  /// @brief Forwards a real function to its derived type.
  template <class Derived>
  struct DeviceView<::Rodin::Variational::RealFunctionBase<Derived>> : DeviceView<Derived>
  {
    /// @param node Function expression.
    DeviceView(const ::Rodin::Variational::RealFunctionBase<Derived>& node)
      : DeviceView<Derived>(static_cast<const Derived&>(node))
    {}
  };

  /// @brief Real function given by a callable taking a device point.
  template <class F>
  struct DeviceView<::Rodin::Variational::RealFunction<F>>
  {
    static_assert(std::is_invocable_r_v<Real, const F&, const Geometry::PointView&>,
      "A function evaluated on a Kokkos mesh must be a KOKKOS_LAMBDA taking a "
      "const Rodin::Kokkos::Geometry::PointView&.");

    /// @param node Function expression.
    DeviceView(const ::Rodin::Variational::RealFunction<F>& node)
      : f(node.getFunction())
    {}

    /// @param p Device point. @return Function value.
    KOKKOS_INLINE_FUNCTION
    Real getValue(const Geometry::PointView& p) const
    {
      return f(p);
    }

    /// @brief Device copy of the callable.
    F f;
  };

  /// @brief Real constant.
  template <>
  struct DeviceView<::Rodin::Variational::RealFunction<Real>>
  {
    /// @param node Function expression.
    DeviceView(const ::Rodin::Variational::RealFunction<Real>& node)
      : value(node.getValue())
    {}

    /// @param p Device point, unused. @return Constant value.
    KOKKOS_INLINE_FUNCTION
    Real getValue([[maybe_unused]] const Geometry::PointView& p) const
    {
      return value;
    }

    /// @brief Constant value.
    Real value;
  };

  /// @brief Sum of two functions.
  template <class LHSDerived, class RHSDerived>
  struct DeviceView<
    ::Rodin::Variational::Sum<::Rodin::Variational::FunctionBase<LHSDerived>,
      ::Rodin::Variational::FunctionBase<RHSDerived>>>
  {
    /// @param node Function expression.
    DeviceView(
      const ::Rodin::Variational::Sum<::Rodin::Variational::FunctionBase<LHSDerived>,
        ::Rodin::Variational::FunctionBase<RHSDerived>>& node)
      : lhs(node.getLHS()),
        rhs(node.getRHS())
    {}

    /// @param p Device point. @return Function value.
    KOKKOS_INLINE_FUNCTION
    auto getValue(const Geometry::PointView& p) const
    {
      return add(lhs.getValue(p), rhs.getValue(p));
    }

    /// @brief Left operand.
    DeviceView<::Rodin::Variational::FunctionBase<LHSDerived>> lhs;

    /// @brief Right operand.
    DeviceView<::Rodin::Variational::FunctionBase<RHSDerived>> rhs;
  };

  /// @brief Product of two functions.
  template <class LHSDerived, class RHSDerived>
  struct DeviceView<
    ::Rodin::Variational::Mult<::Rodin::Variational::FunctionBase<LHSDerived>,
      ::Rodin::Variational::FunctionBase<RHSDerived>>>
  {
    /// @param node Function expression.
    DeviceView(
      const ::Rodin::Variational::Mult<::Rodin::Variational::FunctionBase<LHSDerived>,
        ::Rodin::Variational::FunctionBase<RHSDerived>>& node)
      : lhs(node.getLHS()),
        rhs(node.getRHS())
    {}

    /// @param p Device point. @return Function value.
    KOKKOS_INLINE_FUNCTION
    auto getValue(const Geometry::PointView& p) const
    {
      return mult(lhs.getValue(p), rhs.getValue(p));
    }

    /// @brief Left operand.
    DeviceView<::Rodin::Variational::FunctionBase<LHSDerived>> lhs;

    /// @brief Right operand.
    DeviceView<::Rodin::Variational::FunctionBase<RHSDerived>> rhs;
  };

  /// @brief Negated function.
  template <class NestedDerived>
  struct DeviceView<
    ::Rodin::Variational::UnaryMinus<::Rodin::Variational::FunctionBase<NestedDerived>>>
  {
    /// @param node Function expression.
    DeviceView(const ::Rodin::Variational::UnaryMinus<
      ::Rodin::Variational::FunctionBase<NestedDerived>>& node)
      : operand(node.getOperand())
    {}

    /// @param p Device point. @return Function value.
    KOKKOS_INLINE_FUNCTION
    auto getValue(const Geometry::PointView& p) const
    {
      return mult(Real(-1), operand.getValue(p));
    }

    /// @brief Operand.
    DeviceView<::Rodin::Variational::FunctionBase<NestedDerived>> operand;
  };

  /// @brief Integer constant.
  template <>
  struct DeviceView<::Rodin::Variational::RealFunction<Integer>>
  {
    /// @param node Function expression.
    DeviceView(const ::Rodin::Variational::RealFunction<Integer>& node)
      : value(node.getValue())
    {}

    /// @param p Device point, unused. @return Constant value.
    KOKKOS_INLINE_FUNCTION
    Real getValue([[maybe_unused]] const Geometry::PointView& p) const
    {
      return value;
    }

    /// @brief Constant value.
    Real value;
  };

  /// @brief Forwards a vector function to its derived type.
  template <class Scalar, class Derived>
  struct DeviceView<::Rodin::Variational::VectorFunctionBase<Scalar, Derived>>
    : DeviceView<Derived>
  {
    /// @param node Function expression.
    DeviceView(const ::Rodin::Variational::VectorFunctionBase<Scalar, Derived>& node)
      : DeviceView<Derived>(static_cast<const Derived&>(node))
    {}
  };

  /// @brief Forwards a matrix function to its derived type.
  template <class Scalar, class Derived>
  struct DeviceView<::Rodin::Variational::MatrixFunctionBase<Scalar, Derived>>
    : DeviceView<Derived>
  {
    /// @param node Function expression.
    DeviceView(const ::Rodin::Variational::MatrixFunctionBase<Scalar, Derived>& node)
      : DeviceView<Derived>(static_cast<const Derived&>(node))
    {}
  };

  /// @brief Forwards a gradient of a function to its derived type.
  template <class Operand, class Derived>
  struct DeviceView<::Rodin::Variational::GradBase<Operand, Derived>>
    : DeviceView<Derived>
  {
    /// @param node Function expression.
    DeviceView(const ::Rodin::Variational::GradBase<Operand, Derived>& node)
      : DeviceView<Derived>(static_cast<const Derived&>(node))
    {}
  };

  /// @brief Forwards a Jacobian of a function to its derived type.
  template <class Operand, class Derived>
  struct DeviceView<::Rodin::Variational::JacobianBase<Operand, Derived>>
    : DeviceView<Derived>
  {
    /// @param node Function expression.
    DeviceView(const ::Rodin::Variational::JacobianBase<Operand, Derived>& node)
      : DeviceView<Derived>(static_cast<const Derived&>(node))
    {}
  };

  /// @brief Forwards a divergence of a function to its derived type.
  template <class Operand, class Derived>
  struct DeviceView<::Rodin::Variational::DivBase<Operand, Derived>> : DeviceView<Derived>
  {
    /// @param node Function expression.
    DeviceView(const ::Rodin::Variational::DivBase<Operand, Derived>& node)
      : DeviceView<Derived>(static_cast<const Derived&>(node))
    {}
  };
}

namespace Rodin::Kokkos::Variational::Internal
{
  /**
   * @brief Device views of the components of a vector function.
   *
   * Recursive storage, so that the views of components of distinct types are
   * held by value without std::tuple in device code.
   *
   * @tparam Nodes Component expression types.
   */
  template <class... Nodes>
  struct DeviceViewPack
  {
    /// @brief Evaluates no component.
    KOKKOS_INLINE_FUNCTION
    void evaluate([[maybe_unused]] Geometry::Vector& out,
      [[maybe_unused]] const Geometry::PointView& p, [[maybe_unused]] Index i) const
    {}
  };

  /// @brief Device views of a first component and of the remaining ones.
  template <class Head, class... Tail>
  struct DeviceViewPack<Head, Tail...>
  {
    /// @param head First component. @param tail Remaining components.
    DeviceViewPack(const Head& head, const Tail&... tail)
      : first(head),
        rest(tail...)
    {}

    /**
     * @brief Evaluates the components from a given position.
     * @param out Vector receiving the values.
     * @param p Device point.
     * @param i Position of the first component.
     */
    KOKKOS_INLINE_FUNCTION
    void evaluate(Geometry::Vector& out, const Geometry::PointView& p, Index i) const
    {
      out(i) = first.getValue(p);
      rest.evaluate(out, p, i + 1);
    }

    /// @brief First component.
    DeviceView<Head> first;

    /// @brief Remaining components.
    DeviceViewPack<Tail...> rest;
  };
}

namespace Rodin::Kokkos::Variational
{
  /**
   * @brief Vector function given component-wise.
   *
   * Each component is a real constant or a callable taking a device point.
   */
  template <class V, class... Values>
  struct DeviceView<::Rodin::Variational::VectorFunction<V, Values...>>
  {
    static_assert(1 + sizeof...(Values) <= 3,
      "Vector functions on Kokkos meshes have at most three components.");

    /// @brief Component views.
    using PackType = Internal::DeviceViewPack<::Rodin::Variational::RealFunction<V>,
      ::Rodin::Variational::RealFunction<Values>...>;

    /// @param node Function expression.
    DeviceView(const ::Rodin::Variational::VectorFunction<V, Values...>& node)
      : components(std::make_from_tuple<PackType>(node.getFunctions()))
    {}

    /// @param p Device point. @return Function value.
    KOKKOS_INLINE_FUNCTION
    Geometry::Vector getValue(const Geometry::PointView& p) const
    {
      Geometry::Vector result(1 + sizeof...(Values));
      components.evaluate(result, p, 0);
      return result;
    }

    /// @brief Component views.
    PackType components;
  };

  /**
   * @brief Coefficients of a P1 grid function in device memory.
   *
   * The coefficients are copied to the mesh memory space when the view is
   * built, so the view reflects the grid function at lowering time. Values and
   * derivatives are interpolated from the vertices of the polytope holding the
   * point, as on the host.
   *
   * @tparam Range Range type of the space, Real or a spatial vector.
   * @tparam Backend Kokkos backend of the mesh context.
   */
  template <class Range, class Backend>
  class P1GridFunctionView
  {
    public:
      /// @brief Finite element space.
      using FESType = DeviceP1Space<Range, Backend>;

      /// @brief Memory space of the coefficients.
      using MemorySpace = typename FESType::MeshType::MemorySpace;

      /**
       * @brief Copies the coefficients of a grid function to the device.
       * @param fes Space of the grid function.
       * @param data Coefficients, numbered as by the space.
       */
      template <class Data>
      P1GridFunctionView(const FESType& fes, const Data& data)
        : m_vdim(fes.getVectorDimension()),
          m_vertexCount(fes.getMesh().getVertexCount())
      {
        const auto& execution = fes.getMesh().getContext().getExecutionSpace();
        const Index size = data.size();
        m_data = ::Kokkos::View<Real*, MemorySpace>(
          ::Kokkos::view_alloc(execution, std::string("Rodin grid function")), size);
        auto host = ::Kokkos::create_mirror_view(m_data);
        for (Index i = 0; i < size; ++i)
          host(i) = data.coeff(i);
        ::Kokkos::deep_copy(execution, m_data, host);
        execution.fence("Rodin grid function transfer");
      }

      /// @brief Returns the vector dimension. @return Component count.
      KOKKOS_INLINE_FUNCTION
      Index getVectorDimension() const
      {
        return m_vdim;
      }

      /**
       * @brief Interpolates one component.
       * @param p Device point.
       * @param r Component index.
       * @return Component @f$ u_r(p) @f$.
       */
      KOKKOS_INLINE_FUNCTION
      Real interpolate(const Geometry::PointView& p, Index r) const
      {
        const auto vertices = p.getPolytope().getVertices();
        Geometry::Vector reference(p.getReferenceCoordinates().size());
        Real result = 0;
        for (Index a = 0; a < vertices.size(); ++a)
        {
          Real phi;
          p.evaluateBasis(a, phi, reference);
          result += phi * m_data(vertices[a] + r * m_vertexCount);
        }
        return result;
      }

      /**
       * @brief Interpolates the physical Jacobian.
       * @param p Device point.
       * @return Matrix with one row per component, row @f$ r @f$ holding
       * @f$ \nabla u_r(p) @f$.
       */
      KOKKOS_INLINE_FUNCTION
      Geometry::Matrix interpolateJacobian(const Geometry::PointView& p) const
      {
        const auto vertices = p.getPolytope().getVertices();
        const Geometry::Matrix inverse = p.getJacobianInverse();
        const Index rdim = inverse.rows();
        const Index sdim = inverse.cols();
        Geometry::Matrix result(m_vdim, sdim);
        for (Index a = 0; a < vertices.size(); ++a)
        {
          Real phi;
          Geometry::Vector reference(rdim);
          p.evaluateBasis(a, phi, reference);
          for (Index c = 0; c < sdim; ++c)
          {
            Real gradient = 0;
            for (Index k = 0; k < rdim; ++k)
              gradient += inverse(k, c) * reference(k);
            for (Index r = 0; r < m_vdim; ++r)
              result(r, c) += gradient * m_data(vertices[a] + r * m_vertexCount);
          }
        }
        return result;
      }

    private:
      Index m_vdim;
      Index m_vertexCount;
      ::Kokkos::View<Real*, MemorySpace> m_data;
  };

  /// @brief P1 grid function.
  template <class Range, class Backend, class Data>
  struct DeviceView<
    ::Rodin::Variational::GridFunction<DeviceP1Space<Range, Backend>, Data>>
  {
    /// @param node Grid function.
    DeviceView(
      const ::Rodin::Variational::GridFunction<DeviceP1Space<Range, Backend>, Data>&
        node)
      : u(node.getFiniteElementSpace(), node.getData())
    {}

    /// @param fes Space of the grid function. @param data Coefficients.
    DeviceView(const DeviceP1Space<Range, Backend>& fes, const Data& data)
      : u(fes, data)
    {}

    /// @param p Device point. @return Function value.
    KOKKOS_INLINE_FUNCTION
    auto getValue(const Geometry::PointView& p) const
    {
      if constexpr (std::is_same_v<Range, Real>)
      {
        return u.interpolate(p, 0);
      }
      else
      {
        Geometry::Vector result(u.getVectorDimension());
        for (Index r = 0; r < u.getVectorDimension(); ++r)
          result(r) = u.interpolate(p, r);
        return result;
      }
    }

    /// @brief Device coefficients.
    P1GridFunctionView<Range, Backend> u;
  };

  /// @brief Grid function referenced by an expression.
  template <class Range, class Backend, class Data>
  struct DeviceView<::Rodin::Variational::GridFunctionBaseReference<
    ::Rodin::Variational::GridFunction<DeviceP1Space<Range, Backend>, Data>>>
    : DeviceView<::Rodin::Variational::GridFunction<DeviceP1Space<Range, Backend>, Data>>
  {
    /// @param node Reference to a grid function.
    DeviceView(const ::Rodin::Variational::GridFunctionBaseReference<
      ::Rodin::Variational::GridFunction<DeviceP1Space<Range, Backend>, Data>>& node)
      : DeviceView<
          ::Rodin::Variational::GridFunction<DeviceP1Space<Range, Backend>, Data>>(
          node.getFiniteElementSpace(), node.getData())
    {}
  };

  /// @brief Physical gradient of a scalar P1 grid function.
  template <class Backend, class Data>
  struct DeviceView<::Rodin::Variational::Grad<
    ::Rodin::Variational::GridFunction<DeviceP1<Backend>, Data>>>
  {
    /// @param node Gradient expression.
    DeviceView(const ::Rodin::Variational::Grad<
      ::Rodin::Variational::GridFunction<DeviceP1<Backend>, Data>>& node)
      : u(node.getOperand().getFiniteElementSpace(), node.getOperand().getData())
    {}

    /// @param p Device point. @return Gradient @f$ \nabla u(p) @f$.
    KOKKOS_INLINE_FUNCTION
    Geometry::Vector getValue(const Geometry::PointView& p) const
    {
      const Geometry::Matrix jacobian = u.interpolateJacobian(p);
      Geometry::Vector result(jacobian.cols());
      for (Index c = 0; c < jacobian.cols(); ++c)
        result(c) = jacobian(0, c);
      return result;
    }

    /// @brief Device coefficients.
    P1GridFunctionView<Real, Backend> u;
  };

  /// @brief Physical Jacobian of a vector P1 grid function.
  template <class Backend, class Data>
  struct DeviceView<::Rodin::Variational::Jacobian<
    ::Rodin::Variational::GridFunction<DeviceVectorP1<Backend>, Data>>>
  {
    /// @param node Jacobian expression.
    DeviceView(const ::Rodin::Variational::Jacobian<
      ::Rodin::Variational::GridFunction<DeviceVectorP1<Backend>, Data>>& node)
      : u(node.getOperand().getFiniteElementSpace(), node.getOperand().getData())
    {}

    /// @param p Device point. @return Jacobian @f$ \nabla u(p) @f$.
    KOKKOS_INLINE_FUNCTION
    Geometry::Matrix getValue(const Geometry::PointView& p) const
    {
      return u.interpolateJacobian(p);
    }

    /// @brief Device coefficients.
    P1GridFunctionView<Math::SpatialVector<Real>, Backend> u;
  };

  /// @brief Divergence of a vector P1 grid function.
  template <class Backend, class Data>
  struct DeviceView<::Rodin::Variational::Div<
    ::Rodin::Variational::GridFunction<DeviceVectorP1<Backend>, Data>>>
  {
    /// @param node Divergence expression.
    DeviceView(const ::Rodin::Variational::Div<
      ::Rodin::Variational::GridFunction<DeviceVectorP1<Backend>, Data>>& node)
      : u(node.getOperand().getFiniteElementSpace(), node.getOperand().getData())
    {}

    /// @param p Device point. @return Divergence @f$ \nabla \cdot u(p) @f$.
    KOKKOS_INLINE_FUNCTION
    Real getValue(const Geometry::PointView& p) const
    {
      const Geometry::Matrix jacobian = u.interpolateJacobian(p);
      KOKKOS_ASSERT(jacobian.rows() == jacobian.cols());
      Real result = 0;
      for (Index r = 0; r < jacobian.rows(); ++r)
        result += jacobian(r, r);
      return result;
    }

    /// @brief Device coefficients.
    P1GridFunctionView<Math::SpatialVector<Real>, Backend> u;
  };

  /// @brief Forwards a shape function to its derived type.
  template <class Derived, class FES, ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct DeviceView<::Rodin::Variational::ShapeFunctionBase<Derived, FES, Space>>
    : DeviceView<Derived>
  {
    /// @param node Shape function expression.
    DeviceView(const ::Rodin::Variational::ShapeFunctionBase<Derived, FES, Space>& node)
      : DeviceView<Derived>(static_cast<const Derived&>(node))
    {}
  };

  /// @brief Scalar P1 basis.
  template <class Derived, class Backend,
    ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct DeviceView<
    ::Rodin::Variational::ShapeFunction<Derived, DeviceP1<Backend>, Space>>
  {
    /// @param node Shape function expression.
    DeviceView([[maybe_unused]] const ::Rodin::Variational::ShapeFunction<Derived,
      DeviceP1<Backend>, Space>& node)
    {}

    /// @param ip Quadrature point. @param local Local index. @return Basis value.
    KOKKOS_INLINE_FUNCTION
    Real getBasis(const IntegrationPointView& ip, Index local) const
    {
      return ip.getBasis(local);
    }
  };

  /// @brief Physical gradient of the scalar P1 basis.
  template <class Derived, class Backend,
    ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct DeviceView<::Rodin::Variational::Grad<
    ::Rodin::Variational::ShapeFunction<Derived, DeviceP1<Backend>, Space>>>
  {
    /// @param node Shape function expression.
    DeviceView([[maybe_unused]] const ::Rodin::Variational::Grad<
      ::Rodin::Variational::ShapeFunction<Derived, DeviceP1<Backend>, Space>>& node)
    {}

    /// @param ip Quadrature point. @param local Local index. @return Gradient.
    KOKKOS_INLINE_FUNCTION
    const Geometry::Vector& getBasis(const IntegrationPointView& ip, Index local) const
    {
      return ip.getGradient(local);
    }
  };

  /// @brief Function multiplying a shape function on the left.
  template <class LHSDerived, class RHSDerived, class FES,
    ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct DeviceView<
    ::Rodin::Variational::Mult<::Rodin::Variational::FunctionBase<LHSDerived>,
      ::Rodin::Variational::ShapeFunctionBase<RHSDerived, FES, Space>>>
  {
    /// @param node Shape function expression.
    DeviceView(
      const ::Rodin::Variational::Mult<::Rodin::Variational::FunctionBase<LHSDerived>,
        ::Rodin::Variational::ShapeFunctionBase<RHSDerived, FES, Space>>& node)
      : lhs(node.getLHS()),
        rhs(node.getRHS())
    {}

    /// @param ip Quadrature point. @param local Local index. @return Scaled basis value.
    KOKKOS_INLINE_FUNCTION
    auto getBasis(const IntegrationPointView& ip, Index local) const
    {
      return mult(lhs.getValue(ip.getPoint()), rhs.getBasis(ip, local));
    }

    /// @brief Function factor.
    DeviceView<::Rodin::Variational::FunctionBase<LHSDerived>> lhs;

    /// @brief Shape function factor.
    DeviceView<::Rodin::Variational::ShapeFunctionBase<RHSDerived, FES, Space>> rhs;
  };

  /// @brief Function multiplying a shape function on the right.
  template <class LHSDerived, class RHSDerived, class FES,
    ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct DeviceView<::Rodin::Variational::Mult<
    ::Rodin::Variational::ShapeFunctionBase<LHSDerived, FES, Space>,
    ::Rodin::Variational::FunctionBase<RHSDerived>>>
  {
    /// @param node Shape function expression.
    DeviceView(const ::Rodin::Variational::Mult<
      ::Rodin::Variational::ShapeFunctionBase<LHSDerived, FES, Space>,
      ::Rodin::Variational::FunctionBase<RHSDerived>>& node)
      : lhs(node.getLHS()),
        rhs(node.getRHS())
    {}

    /// @param ip Quadrature point. @param local Local index. @return Scaled basis value.
    KOKKOS_INLINE_FUNCTION
    auto getBasis(const IntegrationPointView& ip, Index local) const
    {
      return mult(lhs.getBasis(ip, local), rhs.getValue(ip.getPoint()));
    }

    /// @brief Shape function factor.
    DeviceView<::Rodin::Variational::ShapeFunctionBase<LHSDerived, FES, Space>> lhs;

    /// @brief Function factor.
    DeviceView<::Rodin::Variational::FunctionBase<RHSDerived>> rhs;
  };

  /// @brief Dot product of a function and a shape function.
  template <class LHSDerived, class RHSDerived, class FES,
    ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct DeviceView<
    ::Rodin::Variational::Dot<::Rodin::Variational::FunctionBase<LHSDerived>,
      ::Rodin::Variational::ShapeFunctionBase<RHSDerived, FES, Space>>>
  {
    /// @param node Shape function expression.
    DeviceView(
      const ::Rodin::Variational::Dot<::Rodin::Variational::FunctionBase<LHSDerived>,
        ::Rodin::Variational::ShapeFunctionBase<RHSDerived, FES, Space>>& node)
      : lhs(node.getLHS()),
        rhs(node.getRHS())
    {}

    /// @param ip Quadrature point. @param local Local index. @return Dot product value.
    KOKKOS_INLINE_FUNCTION
    Real getBasis(const IntegrationPointView& ip, Index local) const
    {
      return dot(lhs.getValue(ip.getPoint()), rhs.getBasis(ip, local));
    }

    /// @brief Function factor.
    DeviceView<::Rodin::Variational::FunctionBase<LHSDerived>> lhs;

    /// @brief Shape function factor.
    DeviceView<::Rodin::Variational::ShapeFunctionBase<RHSDerived, FES, Space>> rhs;
  };

  /**
   * @brief Vector P1 basis.
   *
   * Local index @f$ a\,d + r @f$ is the basis function
   * @f$ \phi_a e_r @f$, as on the host.
   */
  template <class Derived, class Backend,
    ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct DeviceView<
    ::Rodin::Variational::ShapeFunction<Derived, DeviceVectorP1<Backend>, Space>>
  {
    /// @param node Shape function expression.
    DeviceView(const ::Rodin::Variational::ShapeFunction<Derived,
      DeviceVectorP1<Backend>, Space>& node)
      : vdim(node.getFiniteElementSpace().getVectorDimension())
    {}

    /// @param ip Quadrature point. @param local Local index. @return Basis value.
    KOKKOS_INLINE_FUNCTION
    Geometry::Vector getBasis(const IntegrationPointView& ip, Index local) const
    {
      Geometry::Vector result(vdim);
      result(local % vdim) = ip.getBasis(local / vdim);
      return result;
    }

    /// @brief Vector dimension.
    Index vdim;
  };

  /**
   * @brief Physical Jacobian of the vector P1 basis.
   *
   * Row @f$ r @f$ of the Jacobian of @f$ \phi_a e_r @f$ is
   * @f$ \nabla \phi_a @f$, the other rows vanish.
   */
  template <class Derived, class Backend,
    ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct DeviceView<::Rodin::Variational::Jacobian<
    ::Rodin::Variational::ShapeFunction<Derived, DeviceVectorP1<Backend>, Space>>>
  {
    /// @param node Shape function expression.
    DeviceView(
      const ::Rodin::Variational::Jacobian<::Rodin::Variational::ShapeFunction<Derived,
        DeviceVectorP1<Backend>, Space>>& node)
      : vdim(node.getFiniteElementSpace().getVectorDimension())
    {}

    /// @param ip Quadrature point. @param local Local index. @return Jacobian.
    KOKKOS_INLINE_FUNCTION
    Geometry::Matrix getBasis(const IntegrationPointView& ip, Index local) const
    {
      const auto& gradient = ip.getGradient(local / vdim);
      const Index r = local % vdim;
      Geometry::Matrix result(vdim, gradient.size());
      for (Index c = 0; c < gradient.size(); ++c)
        result(r, c) = gradient(c);
      return result;
    }

    /// @brief Vector dimension.
    Index vdim;
  };

  /**
   * @brief Divergence of the vector P1 basis.
   *
   * The divergence of @f$ \phi_a e_r @f$ is @f$ \partial_r \phi_a @f$.
   */
  template <class Derived, class Backend,
    ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct DeviceView<::Rodin::Variational::Div<
    ::Rodin::Variational::ShapeFunction<Derived, DeviceVectorP1<Backend>, Space>>>
  {
    /// @param node Shape function expression.
    DeviceView(
      const ::Rodin::Variational::Div<::Rodin::Variational::ShapeFunction<Derived,
        DeviceVectorP1<Backend>, Space>>& node)
      : vdim(node.getFiniteElementSpace().getVectorDimension())
    {}

    /// @param ip Quadrature point. @param local Local index. @return Divergence.
    KOKKOS_INLINE_FUNCTION
    Real getBasis(const IntegrationPointView& ip, Index local) const
    {
      return ip.getGradient(local / vdim)(local % vdim);
    }

    /// @brief Vector dimension.
    Index vdim;
  };

  /// @brief Sum of two shape functions.
  template <class LHSDerived, class RHSDerived, class FES,
    ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct DeviceView<::Rodin::Variational::Sum<
    ::Rodin::Variational::ShapeFunctionBase<LHSDerived, FES, Space>,
    ::Rodin::Variational::ShapeFunctionBase<RHSDerived, FES, Space>>>
  {
    /// @param node Shape function expression.
    DeviceView(const ::Rodin::Variational::Sum<
      ::Rodin::Variational::ShapeFunctionBase<LHSDerived, FES, Space>,
      ::Rodin::Variational::ShapeFunctionBase<RHSDerived, FES, Space>>& node)
      : lhs(node.getLHS()),
        rhs(node.getRHS())
    {}

    /// @param ip Quadrature point. @param local Local index. @return Basis value.
    KOKKOS_INLINE_FUNCTION
    auto getBasis(const IntegrationPointView& ip, Index local) const
    {
      return add(lhs.getBasis(ip, local), rhs.getBasis(ip, local));
    }

    /// @brief Left operand.
    DeviceView<::Rodin::Variational::ShapeFunctionBase<LHSDerived, FES, Space>> lhs;

    /// @brief Right operand.
    DeviceView<::Rodin::Variational::ShapeFunctionBase<RHSDerived, FES, Space>> rhs;
  };

  /// @brief Negated shape function.
  template <class NestedDerived, class FES,
    ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct DeviceView<::Rodin::Variational::UnaryMinus<
    ::Rodin::Variational::ShapeFunctionBase<NestedDerived, FES, Space>>>
  {
    /// @param node Shape function expression.
    DeviceView(const ::Rodin::Variational::UnaryMinus<
      ::Rodin::Variational::ShapeFunctionBase<NestedDerived, FES, Space>>& node)
      : operand(node.getOperand())
    {}

    /// @param ip Quadrature point. @param local Local index. @return Basis value.
    KOKKOS_INLINE_FUNCTION
    auto getBasis(const IntegrationPointView& ip, Index local) const
    {
      return mult(Real(-1), operand.getBasis(ip, local));
    }

    /// @brief Operand.
    DeviceView<::Rodin::Variational::ShapeFunctionBase<NestedDerived, FES, Space>>
      operand;
  };

  /// @brief Transposed matrix valued shape function.
  template <class NestedDerived, class FES,
    ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct DeviceView<::Rodin::Variational::Transpose<
    ::Rodin::Variational::ShapeFunctionBase<NestedDerived, FES, Space>>>
  {
    /// @param node Shape function expression.
    DeviceView(const ::Rodin::Variational::Transpose<
      ::Rodin::Variational::ShapeFunctionBase<NestedDerived, FES, Space>>& node)
      : operand(node.getOperand())
    {}

    /// @param ip Quadrature point. @param local Local index. @return Basis value.
    KOKKOS_INLINE_FUNCTION
    Geometry::Matrix getBasis(const IntegrationPointView& ip, Index local) const
    {
      return transpose(operand.getBasis(ip, local));
    }

    /// @brief Operand.
    DeviceView<::Rodin::Variational::ShapeFunctionBase<NestedDerived, FES, Space>>
      operand;
  };
}

#endif
