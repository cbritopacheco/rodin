/*
 *          Copyright Carlos BRITO PACHECO 2026.
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_VIEW_COMMON_H
#define RODIN_KOKKOS_VARIATIONAL_VIEW_COMMON_H

/**
 * @file
 * @brief Expression-view declarations, value algebra and bound-basis composition.
 */

#include <concepts>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include <Kokkos_Core.hpp>

#include "Rodin/Configure.h"

#include "Rodin/Variational/Function.h"
#include "Rodin/Variational/RealFunction.h"
#include "Rodin/Variational/VectorFunction.h"
#include "Rodin/Variational/MatrixFunction.h"
#include "Rodin/Variational/GridFunction.h"
#include "Rodin/Variational/ShapeFunction.h"
#include "Rodin/Variational/Sum.h"
#include "Rodin/Variational/Mult.h"
#include "Rodin/Variational/Dot.h"
#include "Rodin/Variational/UnaryMinus.h"
#include "Rodin/Variational/Transpose.h"
#include "Rodin/Variational/BoundaryNormal.h"
#include "Rodin/Variational/Zero.h"
#include "Rodin/Variational/P1/Div.h"
#include "Rodin/Variational/P1/Grad.h"
#include "Rodin/Variational/P1/Jacobian.h"
#include "Rodin/Variational/P1/ShapeFunction.h"

#include "Rodin/Kokkos/Geometry/View.h"

#include "Rodin/Kokkos/Variational/P1.h"
#include "Rodin/Kokkos/Variational/GridFunctionData.h"

#include "IntegrationPoint.h"

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
        std::is_invocable_r<Real, F,
          const ::Rodin::Kokkos::View<::Rodin::Geometry::Point>&>>>>
  {
    /// @brief True, the Kokkos backend evaluating the callable.
    static constexpr bool Value = true;
  };
}

namespace Rodin::Kokkos::Variational
{
  /// @brief P1 space of the given range on a resident mesh with the given backend.
  template <class Range, class Backend>
  using KokkosP1Space = ::Rodin::Variational::P1<Range,
    ::Rodin::Geometry::Mesh<::Rodin::Context::Kokkos<::Rodin::Context::Local, Backend>>>;

  /// @brief Scalar P1 space on a resident mesh with the given backend.
  template <class Backend>
  using KokkosP1 = KokkosP1Space<Real, Backend>;

  /// @brief Vector P1 space on a resident mesh with the given backend.
  template <class Backend>
  using KokkosVectorP1 = KokkosP1Space<::Rodin::Math::SpatialVector<Real>, Backend>;

  /// @brief Tells whether a space is a real valued device space.
  template <class FES>
  concept ScalarSpace = Space<FES> && std::is_same_v<typename FES::RangeType, Real>;

  /// @brief Tells whether a space is a vector valued device space.
  template <class FES>
  concept VectorSpace = Space<FES> &&
    std::is_same_v<typename FES::RangeType, ::Rodin::Math::SpatialVector<Real>>;

  /**
   * @brief Multiplies two scalar values.
   * @param a First factor.
   * @param b Second factor.
   * @return Product @f$ a b @f$.
   */
  KOKKOS_INLINE_FUNCTION Real mult(Real a, Real b)
  {
    return a * b;
  }

  /**
   * @brief Scales a device vector.
   * @param s Scaling factor.
   * @param v Device vector.
   * @return Vector @f$ s v @f$.
   */
  KOKKOS_INLINE_FUNCTION Geometry::Vector mult(Real s, const Geometry::Vector& v)
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
  KOKKOS_INLINE_FUNCTION Geometry::Vector mult(const Geometry::Vector& v, Real s)
  {
    return mult(s, v);
  }

  /**
   * @brief Scales a device matrix.
   * @param s Scaling factor.
   * @param m Device matrix.
   * @return Matrix @f$ s m @f$.
   */
  KOKKOS_INLINE_FUNCTION Geometry::Matrix mult(Real s, const Geometry::Matrix& m)
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
  KOKKOS_INLINE_FUNCTION Geometry::Matrix mult(const Geometry::Matrix& m, Real s)
  {
    return mult(s, m);
  }

  /**
   * @brief Multiplies a device vector by a device matrix.
   * @param m Device matrix.
   * @param v Device vector, with as many entries as @p m has columns.
   * @return Vector @f$ m v @f$.
   */
  KOKKOS_INLINE_FUNCTION Geometry::Vector mult(
    const Geometry::Matrix& m, const Geometry::Vector& v)
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
  KOKKOS_INLINE_FUNCTION Geometry::Matrix mult(
    const Geometry::Matrix& a, const Geometry::Matrix& b)
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
  KOKKOS_INLINE_FUNCTION Real add(Real a, Real b)
  {
    return a + b;
  }

  /**
   * @brief Adds two device vectors.
   * @param a First vector.
   * @param b Second vector, with the dimension of @p a.
   * @return Sum @f$ a + b @f$.
   */
  KOKKOS_INLINE_FUNCTION Geometry::Vector add(
    const Geometry::Vector& a, const Geometry::Vector& b)
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
  KOKKOS_INLINE_FUNCTION Geometry::Matrix add(
    const Geometry::Matrix& a, const Geometry::Matrix& b)
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
  KOKKOS_INLINE_FUNCTION Geometry::Matrix transpose(const Geometry::Matrix& m)
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
  KOKKOS_INLINE_FUNCTION Real dot(Real a, Real b)
  {
    return a * b;
  }

  /**
   * @brief Computes the Euclidean inner product of two device vectors.
   * @param a First vector.
   * @param b Second vector, with the dimension of @p a.
   * @return Inner product @f$ a \cdot b @f$.
   */
  KOKKOS_INLINE_FUNCTION Real dot(const Geometry::Vector& a, const Geometry::Vector& b)
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
  KOKKOS_INLINE_FUNCTION Real dot(const Geometry::Matrix& a, const Geometry::Matrix& b)
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
   * @brief Device vector with a single nonzero component.
   *
   * Value of a vector Lagrange basis @f$ \phi_a e_r @f$. The operations without
   * a dedicated overload use its conversion to a dense device vector.
   */
  struct UnitVector
  {
    /// @return Dense device vector.
    KOKKOS_INLINE_FUNCTION operator Geometry::Vector() const
    {
      Geometry::Vector result(size);
      result(component) = value;
      return result;
    }

    /// @brief Vector dimension.
    Index size;

    /// @brief Index of the nonzero component.
    Index component;

    /// @brief Value of the nonzero component.
    Real value;
  };

  /**
   * @brief Device matrix with a single nonzero row.
   *
   * Jacobian of a vector Lagrange basis @f$ \phi_a e_r @f$, whose row
   * @f$ r @f$ is @f$ \nabla \phi_a @f$. The operations without a dedicated
   * overload use its conversion to a dense device matrix.
   */
  struct RowMatrix
  {
    /// @return Dense device matrix.
    KOKKOS_INLINE_FUNCTION operator Geometry::Matrix() const
    {
      Geometry::Matrix result(rows, values.size());
      for (Index c = 0; c < values.size(); ++c)
        result(row, c) = values(c);
      return result;
    }

    /// @brief Number of rows.
    Index rows;

    /// @brief Index of the nonzero row.
    Index row;

    /// @brief Entries of the nonzero row.
    Geometry::Vector values;
  };

  /// @param s Scaling factor. @param v Unit vector. @return Vector @f$ s v @f$.
  KOKKOS_INLINE_FUNCTION UnitVector mult(Real s, const UnitVector& v)
  {
    return {v.size, v.component, s * v.value};
  }

  /// @param v Unit vector. @param s Scaling factor. @return Vector @f$ v s @f$.
  KOKKOS_INLINE_FUNCTION UnitVector mult(const UnitVector& v, Real s)
  {
    return mult(s, v);
  }

  /// @param s Scaling factor. @param m Row matrix. @return Matrix @f$ s m @f$.
  KOKKOS_INLINE_FUNCTION RowMatrix mult(Real s, const RowMatrix& m)
  {
    return {m.rows, m.row, mult(s, m.values)};
  }

  /// @param m Row matrix. @param s Scaling factor. @return Matrix @f$ m s @f$.
  KOKKOS_INLINE_FUNCTION RowMatrix mult(const RowMatrix& m, Real s)
  {
    return mult(s, m);
  }

  /**
   * @brief Multiplies a device vector by a row matrix.
   * @param m Row matrix.
   * @param v Device vector, with as many entries as @p m has columns.
   * @return Vector @f$ m v @f$.
   */
  KOKKOS_INLINE_FUNCTION UnitVector mult(const RowMatrix& m, const Geometry::Vector& v)
  {
    return {m.rows, m.row, dot(m.values, v)};
  }

  /// @param a First unit vector. @param b Second unit vector. @return @f$ a \cdot b @f$.
  KOKKOS_INLINE_FUNCTION Real dot(const UnitVector& a, const UnitVector& b)
  {
    KOKKOS_ASSERT(a.size == b.size);
    return a.component == b.component ? a.value * b.value : Real(0);
  }

  /// @param a Device vector. @param b Unit vector. @return @f$ a \cdot b @f$.
  KOKKOS_INLINE_FUNCTION Real dot(const Geometry::Vector& a, const UnitVector& b)
  {
    KOKKOS_ASSERT(a.size() == b.size);
    return a(b.component) * b.value;
  }

  /// @param a Unit vector. @param b Device vector. @return @f$ a \cdot b @f$.
  KOKKOS_INLINE_FUNCTION Real dot(const UnitVector& a, const Geometry::Vector& b)
  {
    return dot(b, a);
  }

  /// @param a First row matrix. @param b Second row matrix. @return @f$ a : b @f$.
  KOKKOS_INLINE_FUNCTION Real dot(const RowMatrix& a, const RowMatrix& b)
  {
    KOKKOS_ASSERT(a.rows == b.rows);
    return a.row == b.row ? dot(a.values, b.values) : Real(0);
  }

  /// @param a Device matrix. @param b Row matrix. @return @f$ a : b @f$.
  KOKKOS_INLINE_FUNCTION Real dot(const Geometry::Matrix& a, const RowMatrix& b)
  {
    KOKKOS_ASSERT(a.rows() == b.rows && a.cols() == b.values.size());
    Real result = 0;
    for (Index c = 0; c < b.values.size(); ++c)
      result += a(b.row, c) * b.values(c);
    return result;
  }

  /// @param a Row matrix. @param b Device matrix. @return @f$ a : b @f$.
  KOKKOS_INLINE_FUNCTION Real dot(const RowMatrix& a, const Geometry::Matrix& b)
  {
    return dot(b, a);
  }

  /**
   * @brief Tells whether the inner product of two values may be nonzero.
   *
   * Kernels skip the pairs that do not overlap. Values without sparsity
   * always overlap.
   *
   * @param lhs Left operand, ignored.
   * @param rhs Right operand, ignored.
   * @return True.
   */
  template <class LHS, class RHS>
  KOKKOS_INLINE_FUNCTION constexpr bool overlaps(
    [[maybe_unused]] const LHS& lhs, [[maybe_unused]] const RHS& rhs)
  {
    return true;
  }

  /// @param a First unit vector. @param b Second unit vector. @return Whether the components match.
  KOKKOS_INLINE_FUNCTION bool overlaps(const UnitVector& a, const UnitVector& b)
  {
    return a.component == b.component;
  }

  /// @param a First row matrix. @param b Second row matrix. @return Whether the rows match.
  KOKKOS_INLINE_FUNCTION bool overlaps(const RowMatrix& a, const RowMatrix& b)
  {
    return a.row == b.row;
  }

  /**
   * @brief Tells whether a basis value type is sparse by component.
   *
   * Values of these types come from the vector basis @f$ \phi_a e_r @f$,
   * scaled by scalars or multiplied on the right by vectors: the value of
   * the scalar degree of freedom a and component r has its nonzero
   * component or row at r, and its entries do not depend on r. The inner
   * product of two such values vanishes across components and, within one
   * component, equals the product of the values of component 0, so the
   * bilinear kernels integrate one scalar matrix and replicate it on the
   * diagonal of the component blocks. An operation producing these types
   * keeps this property.
   *
   * @tparam T Basis value type.
   */
  template <class T>
  inline constexpr bool IsComponentSparse =
    std::is_same_v<T, UnitVector> || std::is_same_v<T, RowMatrix>;

  /// @param v Unit vector. @return Value of its nonzero component.
  KOKKOS_INLINE_FUNCTION Real getPayload(const UnitVector& v)
  {
    return v.value;
  }

  /// @param m Row matrix. @return Entries of its nonzero row.
  KOKKOS_INLINE_FUNCTION Geometry::Vector getPayload(const RowMatrix& m)
  {
    return m.values;
  }

  /**
   * @brief Gets the value stored by a bilinear kernel.
   *
   * A component diagonal kernel (IsComponentSparse values on both sides)
   * stores the payload of the values of component 0 only: their inner
   * product is the inner product of the payloads, and the payload is
   * smaller, so more quadrature points fit the team scratch memory.
   *
   * @tparam Reduce Whether the kernel is component diagonal.
   * @param v Basis value.
   * @return Payload of @p v if Reduce, @p v otherwise.
   */
  template <bool Reduce, class T>
  KOKKOS_INLINE_FUNCTION auto getStored(const T& v)
  {
    if constexpr (Reduce)
      return getPayload(v);
    else
      return v;
  }

  /**
   * @brief Device-callable view of a form language expression.
   *
   * A view is built on the host from the expression and copied by value into
   * kernels. It holds no host pointer and no virtual function. Function views
   * provide `getValue(const ::Rodin::Kokkos::View<::Rodin::Geometry::Point>&)`, while shape function views
   * provide `getBasis(const IntegrationPointView&, Index, Index)`, which
   * returns the value of component r of the basis function of the local scalar
   * degree of freedom a. A scalar view has the single component 0.
   *
   * The primary template is left undefined so that an expression without a
   * device lowering is rejected at compile time.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | `View<FunctionBase<Derived>>` | Forwards to the derived function. |
   * | `View<ScalarFunctionBase<Scalar, Derived>>` | Forwards to the derived function. |
   * | `View<RealFunctionBase<Derived>>` | Forwards to the derived function. |
   * | `View<RealFunction<F>>` | Callable taking a device point. |
   * | `View<RealFunction<Real>>` | Real constant. |
   * | `View<Sum<FunctionBase<L>, FunctionBase<R>>>` | Sum of functions. |
   * | `View<Mult<FunctionBase<L>, FunctionBase<R>>>` | Product of functions. |
   * | `View<Dot<FunctionBase<L>, FunctionBase<R>>>` | Dot product of functions. |
   * | `View<UnaryMinus<FunctionBase<N>>>` | Negated function. |
   * | `View<RealFunction<Integer>>` | Integer constant. |
   * | `View<VectorFunctionBase<Scalar, Derived>>` | Forwards to the derived function. |
   * | `View<MatrixFunctionBase<Scalar, Derived>>` | Forwards to the derived function. |
   * | `View<MatrixFunction<Matrix>>` | Constant matrix. |
   * | `View<GradBase<Operand, Derived>>` | Forwards to the derived function. |
   * | `View<JacobianBase<Operand, Derived>>` | Forwards to the derived function. |
   * | `View<DivBase<Operand, Derived>>` | Forwards to the derived function. |
   * | `View<VectorFunction<V, Values...>>` | Vector function given component-wise. |
   * | `View<Zero<Real>>` | Scalar zero. |
   * | `View<Zero<SpatialVector<Real>>>` | Vector zero. |
   * | `View<BoundaryNormal>` | Outward unit normal of the boundary faces. |
   * | `View<GridFunctionBaseReference<GridFunction<FES, Data>>>` | Referenced grid function of a device space. |
   * | `View<GridFunction<FES, Data>>` | Grid function of a device space. |
   * | `View<Grad<GridFunction<FES, Data>>>` | Physical gradient of a scalar grid function. |
   * | `View<Jacobian<GridFunction<FES, Data>>>` | Physical Jacobian of a vector grid function. |
   * | `View<Div<GridFunction<FES, Data>>>` | Divergence of a vector grid function. |
   * | `View<ShapeFunctionBase<Derived, FES, Space>>` | Forwards to the derived shape function. |
   * | `View<ShapeFunction<Derived, FES, Space>>` | Scalar or vector basis value of a device space. |
   * | `View<Grad<ShapeFunction<Derived, FES, Space>>>` | Physical gradient of a scalar basis. |
   * | `View<Jacobian<ShapeFunction<Derived, FES, Space>>>` | Physical Jacobian of a vector basis. |
   * | `View<Div<ShapeFunction<Derived, FES, Space>>>` | Divergence of a vector basis. |
   * | `View<Sum<ShapeFunctionBase<L, FES, Space>, ShapeFunctionBase<R, FES, Space>>>` | Sum of shape functions. |
   * | `View<UnaryMinus<ShapeFunctionBase<N, FES, Space>>>` | Negated shape function. |
   * | `View<Transpose<ShapeFunctionBase<N, FES, Space>>>` | Transposed shape function. |
   * | `View<Mult<FunctionBase<L>, ShapeFunctionBase<R, FES, Space>>>` | Function times shape function. |
   * | `View<Mult<ShapeFunctionBase<L, FES, Space>, FunctionBase<R>>>` | Shape function times function. |
   * | `View<Dot<FunctionBase<L>, ShapeFunctionBase<R, FES, Space>>>` | Dot product of a function and a shape function. |
   *
   * @tparam Node Form language expression type.
   */
  template <class Node>
  struct View;

  /**
   * @brief Evaluates a function view at a point.
   *
   * Views that use the basis already evaluated at an integration point
   * provide `getValue(const IntegrationPointView&)`. The others are
   * evaluated at the underlying device point.
   *
   * @param view Function view.
   * @param p Device point or integration point.
   * @return Function value.
   */
  template <class View, class Point>
  KOKKOS_INLINE_FUNCTION auto evaluate(const View& view, const Point& p)
  {
    if constexpr (requires { view.getValue(p); })
      return view.getValue(p);
    else
      return view.getValue(p.getPoint());
  }

  namespace Internal
  {
    /// @brief Shape function view evaluated at a fixed integration point.
    template <class View, class Point = IntegrationPointView>
    struct BoundBasis
    {
      /**
       * @param dof Local scalar degree of freedom.
       * @param component Vector component. @return Basis value.
       */
      KOKKOS_INLINE_FUNCTION decltype(auto) getBasis(Index dof, Index component) const
      {
        return view.getBasis(ip, dof, component);
      }

      /// @brief Shape function view.
      const View& view;

      /// @brief Integration point.
      const Point& ip;
    };

    /// @brief Basis multiplied on the left by a bound function value.
    template <class Value, class Basis>
    struct LeftScaledBasis
    {
      /**
       * @param dof Local scalar degree of freedom.
       * @param component Vector component. @return Scaled basis value.
       */
      KOKKOS_INLINE_FUNCTION auto getBasis(Index dof, Index component) const
      {
        return mult(value, basis.getBasis(dof, component));
      }

      /// @brief Function value.
      Value value;

      /// @brief Bound basis.
      Basis basis;
    };

    /// @brief Basis multiplied on the right by a bound function value.
    template <class Basis, class Value>
    struct RightScaledBasis
    {
      /**
       * @param dof Local scalar degree of freedom.
       * @param component Vector component. @return Scaled basis value.
       */
      KOKKOS_INLINE_FUNCTION auto getBasis(Index dof, Index component) const
      {
        return mult(basis.getBasis(dof, component), value);
      }

      /// @brief Bound basis.
      Basis basis;

      /// @brief Function value.
      Value value;
    };

    /// @brief Dot product of a bound function value and a basis.
    template <class Value, class Basis>
    struct DotBasis
    {
      /**
       * @param dof Local scalar degree of freedom.
       * @param component Vector component. @return Dot product value.
       */
      KOKKOS_INLINE_FUNCTION Real getBasis(Index dof, Index component) const
      {
        return dot(value, basis.getBasis(dof, component));
      }

      /// @brief Function value.
      Value value;

      /// @brief Bound basis.
      Basis basis;
    };

    /// @brief Sum of two bound bases.
    template <class LHS, class RHS>
    struct SumBasis
    {
      /**
       * @param dof Local scalar degree of freedom.
       * @param component Vector component. @return Sum of the basis values.
       */
      KOKKOS_INLINE_FUNCTION auto getBasis(Index dof, Index component) const
      {
        return add(lhs.getBasis(dof, component), rhs.getBasis(dof, component));
      }

      /// @brief Left bound basis.
      LHS lhs;

      /// @brief Right bound basis.
      RHS rhs;
    };

    /// @brief Opposite of a bound basis.
    template <class Basis>
    struct NegatedBasis
    {
      /**
       * @param dof Local scalar degree of freedom.
       * @param component Vector component. @return Opposite basis value.
       */
      KOKKOS_INLINE_FUNCTION auto getBasis(Index dof, Index component) const
      {
        return mult(Real(-1), basis.getBasis(dof, component));
      }

      /// @brief Bound basis.
      Basis basis;
    };

    /// @brief Transpose of a bound basis.
    template <class Basis>
    struct TransposedBasis
    {
      /**
       * @param dof Local scalar degree of freedom.
       * @param component Vector component. @return Transposed basis value.
       */
      KOKKOS_INLINE_FUNCTION Geometry::Matrix getBasis(Index dof, Index component) const
      {
        return transpose(basis.getBasis(dof, component));
      }

      /// @brief Bound basis.
      Basis basis;
    };
  }

  /**
   * @brief Binds a shape function view to an integration point.
   *
   * The functions multiplying the basis are evaluated once at the point,
   * instead of once per local basis function. Views that hold no function
   * are evaluated through their `getBasis(const IntegrationPointView&,
   * Index, Index)`.
   *
   * @param view Shape function view.
   * @param ip Integration point, which must outlive the result. A
   * CellIntegrationPointView is passed on to the functions.
   * @return Object providing `getBasis(Index, Index)`.
   */
  template <class View, class Point>
    requires std::derived_from<Point, IntegrationPointView>
  KOKKOS_INLINE_FUNCTION auto bindBasis(const View& view, const Point& ip)
  {
    if constexpr (requires { view.bind(ip); })
      return view.bind(ip);
    else
      return Internal::BoundBasis<View, Point>{view, ip};
  }

}

#endif
