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

#include <type_traits>

#include <Kokkos_Core.hpp>

#include "Rodin/Configure.h"

#include "Rodin/Variational/Function.h"
#include "Rodin/Variational/RealFunction.h"
#include "Rodin/Variational/ShapeFunction.h"
#include "Rodin/Variational/Sum.h"
#include "Rodin/Variational/Mult.h"
#include "Rodin/Variational/UnaryMinus.h"
#include "Rodin/Variational/P1/Grad.h"
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
  /// @brief Scalar P1 space on a resident mesh with the given backend.
  template <class Backend>
  using DeviceP1 = ::Rodin::Variational::P1<Real,
    ::Rodin::Geometry::Mesh<::Rodin::Context::Kokkos<::Rodin::Context::Local, Backend>>>;

  /**
   * @brief Scales a scalar value.
   * @param s Scaling factor.
   * @param v Scalar value.
   * @return Product @f$ s v @f$.
   */
  KOKKOS_INLINE_FUNCTION
  Real scale(Real s, Real v)
  {
    return s * v;
  }

  /**
   * @brief Scales a device vector.
   * @param s Scaling factor.
   * @param v Device vector.
   * @return Vector @f$ s v @f$.
   */
  KOKKOS_INLINE_FUNCTION
  Geometry::Vector scale(Real s, const Geometry::Vector& v)
  {
    Geometry::Vector result(v.size());
    for (Index i = 0; i < v.size(); ++i)
      result(i) = s * v(i);
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
   * | `DeviceView<ShapeFunctionBase<Derived, FES, Space>>` | Forwards to the derived shape function. |
   * | `DeviceView<ShapeFunction<Derived, DeviceP1<B>, Space>>` | P1 basis value. |
   * | `DeviceView<Grad<ShapeFunction<Derived, DeviceP1<B>, Space>>>` | Physical P1 gradient. |
   * | `DeviceView<Mult<FunctionBase<L>, ShapeFunctionBase<R, FES, Space>>>` | Function times shape function. |
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
    Real getValue(const Geometry::PointView& p) const
    {
      return lhs.getValue(p) + rhs.getValue(p);
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
    Real getValue(const Geometry::PointView& p) const
    {
      return lhs.getValue(p) * rhs.getValue(p);
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
    Real getValue(const Geometry::PointView& p) const
    {
      return -operand.getValue(p);
    }

    /// @brief Operand.
    DeviceView<::Rodin::Variational::FunctionBase<NestedDerived>> operand;
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
      return scale(lhs.getValue(ip.getPoint()), rhs.getBasis(ip, local));
    }

    /// @brief Function factor.
    DeviceView<::Rodin::Variational::FunctionBase<LHSDerived>> lhs;

    /// @brief Shape function factor.
    DeviceView<::Rodin::Variational::ShapeFunctionBase<RHSDerived, FES, Space>> rhs;
  };
}

#endif
