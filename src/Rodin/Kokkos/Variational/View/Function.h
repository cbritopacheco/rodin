/*
 *          Copyright Carlos BRITO PACHECO 2026.
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_VIEW_FUNCTION_H
#define RODIN_KOKKOS_VARIATIONAL_VIEW_FUNCTION_H

/**
 * @file
 * @brief Device evaluation of known function expressions.
 */

#include "Common.h"

namespace Rodin::Kokkos::Variational
{
  /// @brief Forwards a function to its derived type.
  template <class Derived>
  struct View<::Rodin::Variational::FunctionBase<Derived>> : View<Derived>
  {
    /// @param node Function expression.
    View(const ::Rodin::Variational::FunctionBase<Derived>& node)
      : View<Derived>(static_cast<const Derived&>(node))
    {}
  };

  /// @brief Forwards a scalar function to its derived type.
  template <class Scalar, class Derived>
  struct View<::Rodin::Variational::ScalarFunctionBase<Scalar, Derived>> : View<Derived>
  {
    /// @param node Function expression.
    View(const ::Rodin::Variational::ScalarFunctionBase<Scalar, Derived>& node)
      : View<Derived>(static_cast<const Derived&>(node))
    {}
  };

  /// @brief Forwards a real function to its derived type.
  template <class Derived>
  struct View<::Rodin::Variational::RealFunctionBase<Derived>> : View<Derived>
  {
    /// @param node Function expression.
    View(const ::Rodin::Variational::RealFunctionBase<Derived>& node)
      : View<Derived>(static_cast<const Derived&>(node))
    {}
  };

  /// @brief Real function given by a callable taking a device point.
  template <class F>
  struct View<::Rodin::Variational::RealFunction<F>>
  {
    static_assert(std::is_invocable_r_v<Real, const F&,
                    const ::Rodin::Kokkos::View<::Rodin::Geometry::Point>&>,
      "A function evaluated on a Kokkos mesh must be a KOKKOS_LAMBDA taking a "
      "const ::Rodin::Kokkos::View<::Rodin::Geometry::Point>&.");

    /// @param node Function expression.
    View(const ::Rodin::Variational::RealFunction<F>& node)
      : f(node.getFunction())
    {}

    /// @param p Device point. @return Function value.
    KOKKOS_INLINE_FUNCTION Real getValue(
      const ::Rodin::Kokkos::View<::Rodin::Geometry::Point>& p) const
    {
      return f(p);
    }

    /// @brief Device copy of the callable.
    F f;
  };

  /// @brief Real function wrapping a real function, as a vector component.
  template <class F>
  struct View<::Rodin::Variational::RealFunction<::Rodin::Variational::RealFunction<F>>>
    : View<::Rodin::Variational::RealFunction<F>>
  {
    /// @param node Function expression.
    View(
      const ::Rodin::Variational::RealFunction<::Rodin::Variational::RealFunction<F>>&
        node)
      : View<::Rodin::Variational::RealFunction<F>>(node.getFunction())
    {}
  };

  /// @brief Real constant.
  template <>
  struct View<::Rodin::Variational::RealFunction<Real>>
  {
    /// @param node Function expression.
    View(const ::Rodin::Variational::RealFunction<Real>& node)
      : value(node.getValue())
    {}

    /// @param p Device point, unused. @return Constant value.
    KOKKOS_INLINE_FUNCTION Real getValue(
      [[maybe_unused]] const ::Rodin::Kokkos::View<::Rodin::Geometry::Point>& p) const
    {
      return value;
    }

    /// @brief Constant value.
    Real value;
  };

  /// @brief Sum of two functions.
  template <class LHSDerived, class RHSDerived>
  struct View<::Rodin::Variational::Sum<::Rodin::Variational::FunctionBase<LHSDerived>,
    ::Rodin::Variational::FunctionBase<RHSDerived>>>
  {
    /// @param node Function expression.
    View(const ::Rodin::Variational::Sum<::Rodin::Variational::FunctionBase<LHSDerived>,
      ::Rodin::Variational::FunctionBase<RHSDerived>>& node)
      : lhs(node.getLHS()),
        rhs(node.getRHS())
    {}

    /// @param p Device or integration point. @return Function value.
    template <class Point>
    KOKKOS_INLINE_FUNCTION auto getValue(const Point& p) const
    {
      return add(evaluate(lhs, p), evaluate(rhs, p));
    }

    /// @brief Left operand.
    View<::Rodin::Variational::FunctionBase<LHSDerived>> lhs;

    /// @brief Right operand.
    View<::Rodin::Variational::FunctionBase<RHSDerived>> rhs;
  };

  /// @brief Product of two functions.
  template <class LHSDerived, class RHSDerived>
  struct View<::Rodin::Variational::Mult<::Rodin::Variational::FunctionBase<LHSDerived>,
    ::Rodin::Variational::FunctionBase<RHSDerived>>>
  {
    /// @param node Function expression.
    View(
      const ::Rodin::Variational::Mult<::Rodin::Variational::FunctionBase<LHSDerived>,
        ::Rodin::Variational::FunctionBase<RHSDerived>>& node)
      : lhs(node.getLHS()),
        rhs(node.getRHS())
    {}

    /// @param p Device or integration point. @return Function value.
    template <class Point>
    KOKKOS_INLINE_FUNCTION auto getValue(const Point& p) const
    {
      return mult(evaluate(lhs, p), evaluate(rhs, p));
    }

    /// @brief Left operand.
    View<::Rodin::Variational::FunctionBase<LHSDerived>> lhs;

    /// @brief Right operand.
    View<::Rodin::Variational::FunctionBase<RHSDerived>> rhs;
  };

  /// @brief Dot product of two functions.
  template <class LHSDerived, class RHSDerived>
  struct View<::Rodin::Variational::Dot<::Rodin::Variational::FunctionBase<LHSDerived>,
    ::Rodin::Variational::FunctionBase<RHSDerived>>>
  {
    /// @param node Function expression.
    View(const ::Rodin::Variational::Dot<::Rodin::Variational::FunctionBase<LHSDerived>,
      ::Rodin::Variational::FunctionBase<RHSDerived>>& node)
      : lhs(node.getLHS()),
        rhs(node.getRHS())
    {}

    /// @param p Device or integration point. @return Function value.
    template <class Point>
    KOKKOS_INLINE_FUNCTION Real getValue(const Point& p) const
    {
      return dot(evaluate(lhs, p), evaluate(rhs, p));
    }

    /// @brief Left operand.
    View<::Rodin::Variational::FunctionBase<LHSDerived>> lhs;

    /// @brief Right operand.
    View<::Rodin::Variational::FunctionBase<RHSDerived>> rhs;
  };

  /// @brief Negated function.
  template <class NestedDerived>
  struct View<
    ::Rodin::Variational::UnaryMinus<::Rodin::Variational::FunctionBase<NestedDerived>>>
  {
    /// @param node Function expression.
    View(const ::Rodin::Variational::UnaryMinus<
      ::Rodin::Variational::FunctionBase<NestedDerived>>& node)
      : operand(node.getOperand())
    {}

    /// @param p Device or integration point. @return Function value.
    template <class Point>
    KOKKOS_INLINE_FUNCTION auto getValue(const Point& p) const
    {
      return mult(Real(-1), evaluate(operand, p));
    }

    /// @brief Operand.
    View<::Rodin::Variational::FunctionBase<NestedDerived>> operand;
  };

  /// @brief Integer constant.
  template <>
  struct View<::Rodin::Variational::RealFunction<Integer>>
  {
    /// @param node Function expression.
    View(const ::Rodin::Variational::RealFunction<Integer>& node)
      : value(node.getValue())
    {}

    /// @param p Device point, unused. @return Constant value.
    KOKKOS_INLINE_FUNCTION Real getValue(
      [[maybe_unused]] const ::Rodin::Kokkos::View<::Rodin::Geometry::Point>& p) const
    {
      return value;
    }

    /// @brief Constant value.
    Real value;
  };

  /// @brief Forwards a vector function to its derived type.
  template <class Scalar, class Derived>
  struct View<::Rodin::Variational::VectorFunctionBase<Scalar, Derived>> : View<Derived>
  {
    /// @param node Function expression.
    View(const ::Rodin::Variational::VectorFunctionBase<Scalar, Derived>& node)
      : View<Derived>(static_cast<const Derived&>(node))
    {}
  };

  /// @brief Forwards a matrix function to its derived type.
  template <class Scalar, class Derived>
  struct View<::Rodin::Variational::MatrixFunctionBase<Scalar, Derived>> : View<Derived>
  {
    /// @param node Function expression.
    View(const ::Rodin::Variational::MatrixFunctionBase<Scalar, Derived>& node)
      : View<Derived>(static_cast<const Derived&>(node))
    {}
  };

  /// @brief Constant matrix of at most three rows and three columns.
  template <class Matrix>
    requires std::is_same_v<Matrix, ::Rodin::Math::Matrix<Real>> ||
    std::is_same_v<Matrix, ::Rodin::Math::SpatialMatrix<Real>>
  struct View<::Rodin::Variational::MatrixFunction<Matrix>>
  {
    /// @param node Function expression.
    View(const ::Rodin::Variational::MatrixFunction<Matrix>& node)
    {
      const auto& matrix = node.getValue();
      if (matrix.rows() > 3 || matrix.cols() > 3)
      {
        Alert::Exception()
          << "Matrix functions on Kokkos meshes have at most three rows and three "
             "columns."
          << Alert::Raise;
      }
      value = Geometry::Matrix(matrix.rows(), matrix.cols());
      for (Index r = 0; r < value.rows(); ++r)
      {
        for (Index c = 0; c < value.cols(); ++c)
          value(r, c) = matrix(r, c);
      }
    }

    /// @param p Device or integration point, unused. @return Constant value.
    template <class Point>
    KOKKOS_INLINE_FUNCTION Geometry::Matrix getValue(
      [[maybe_unused]] const Point& p) const
    {
      return value;
    }

    /// @brief Constant value.
    Geometry::Matrix value;
  };

  /// @brief Forwards a gradient of a function to its derived type.
  template <class Operand, class Derived>
  struct View<::Rodin::Variational::GradBase<Operand, Derived>> : View<Derived>
  {
    /// @param node Function expression.
    View(const ::Rodin::Variational::GradBase<Operand, Derived>& node)
      : View<Derived>(static_cast<const Derived&>(node))
    {}
  };

  /// @brief Forwards a Jacobian of a function to its derived type.
  template <class Operand, class Derived>
  struct View<::Rodin::Variational::JacobianBase<Operand, Derived>> : View<Derived>
  {
    /// @param node Function expression.
    View(const ::Rodin::Variational::JacobianBase<Operand, Derived>& node)
      : View<Derived>(static_cast<const Derived&>(node))
    {}
  };

  /// @brief Forwards a divergence of a function to its derived type.
  template <class Operand, class Derived>
  struct View<::Rodin::Variational::DivBase<Operand, Derived>> : View<Derived>
  {
    /// @param node Function expression.
    View(const ::Rodin::Variational::DivBase<Operand, Derived>& node)
      : View<Derived>(static_cast<const Derived&>(node))
    {}
  };

  /// @brief Scalar zero.
  template <>
  struct View<::Rodin::Variational::Zero<Real>>
  {
    /// @param node Function expression, unused.
    View([[maybe_unused]] const ::Rodin::Variational::Zero<Real>& node) {}

    /// @param p Device point, unused. @return Zero.
    KOKKOS_INLINE_FUNCTION Real getValue(
      [[maybe_unused]] const ::Rodin::Kokkos::View<::Rodin::Geometry::Point>& p) const
    {
      return 0;
    }
  };

  /// @brief Vector zero.
  template <>
  struct View<::Rodin::Variational::Zero<::Rodin::Math::SpatialVector<Real>>>
  {
    /// @param node Function expression.
    View(const ::Rodin::Variational::Zero<::Rodin::Math::SpatialVector<Real>>& node)
      : dimension(node.getDimension())
    {}

    /// @param p Device point, unused. @return Zero vector.
    KOKKOS_INLINE_FUNCTION Geometry::Vector getValue(
      [[maybe_unused]] const ::Rodin::Kokkos::View<::Rodin::Geometry::Point>& p) const
    {
      return Geometry::Vector(dimension);
    }

    /// @brief Number of components.
    size_t dimension;
  };

  /**
   * @brief Outward unit normal of the boundary faces.
   *
   * The normal of a face of a mesh of dimension two or three, embedded in a
   * space of the same dimension, is computed from the face Jacobian as on
   * the host: @f$ (J_{10}, -J_{00}) @f$ for an edge and the cross product of
   * the Jacobian columns for a surface. It is oriented away from the cell
   * incident to the face, whose vertex mean stands for its interior, then
   * normalized. Evaluating it needs the face to cell incidence. Surface
   * meshes embedded in three dimensions and codimension two entities are not
   * supported.
   */
  template <>
  struct View<::Rodin::Variational::BoundaryNormal>
  {
    /// @param node Function expression.
    View([[maybe_unused]] const ::Rodin::Variational::BoundaryNormal& node) {}

    /// @param p Point on a boundary face. @return Outward unit normal.
    KOKKOS_INLINE_FUNCTION Geometry::Vector getValue(
      const ::Rodin::Kokkos::View<::Rodin::Geometry::Point>& p) const
    {
      const auto& face = p.getPolytope();
      const ::Rodin::Kokkos::View<::Rodin::Geometry::MeshBase>& mesh = face.mesh;
      const Index d = face.getDimension();
      const Index sdim = mesh.getSpaceDimension();
      KOKKOS_ASSERT(d + 1 == mesh.getDimension() && sdim == d + 1);
      const Geometry::Matrix jacobian = p.getJacobian();
      Geometry::Vector normal(sdim);
      if (sdim == 2)
      {
        normal(0) = jacobian(1, 0);
        normal(1) = -jacobian(0, 0);
      }
      else
      {
        normal(0) = jacobian(1, 0) * jacobian(2, 1) - jacobian(2, 0) * jacobian(1, 1);
        normal(1) = jacobian(2, 0) * jacobian(0, 1) - jacobian(0, 0) * jacobian(2, 1);
        normal(2) = jacobian(0, 0) * jacobian(1, 1) - jacobian(1, 0) * jacobian(0, 1);
      }
      const auto cells = mesh.getConnectivity().getIncidence(d, d + 1, face.getIndex());
      KOKKOS_ASSERT(cells.size() > 0);
      const auto vertices = mesh.getPolytope(d + 1, cells[0]).getVertices();
      const Geometry::Vector& x = p.getPhysicalCoordinates();
      Geometry::Vector centroid(sdim);
      for (Index a = 0; a < vertices.size(); ++a)
        centroid = add(centroid, mesh.getVertexCoordinates(vertices[a]));
      Real side = 0;
      for (Index c = 0; c < sdim; ++c)
        side += normal(c) * (centroid(c) / vertices.size() - x(c));
      const Real norm = ::Kokkos::sqrt(dot(normal, normal));
      const Real scale = (side > 0 ? Real(-1) : Real(1)) / norm;
      for (Index c = 0; c < sdim; ++c)
        normal(c) *= scale;
      return normal;
    }
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
  struct ViewPack
  {
    /**
     * @brief Evaluates no component.
     * @param out Vector receiving the values, left unchanged.
     * @param p Device or integration point.
     * @param i Position past the last component.
     */
    template <class Point>
    KOKKOS_INLINE_FUNCTION void evaluate([[maybe_unused]] Geometry::Vector& out,
      [[maybe_unused]] const Point& p, [[maybe_unused]] Index i) const
    {}
  };

  /// @brief Device views of a first component and of the remaining ones.
  template <class Head, class... Tail>
  struct ViewPack<Head, Tail...>
  {
    /// @param head First component. @param tail Remaining components.
    ViewPack(const Head& head, const Tail&... tail)
      : first(head),
        rest(tail...)
    {}

    /**
     * @brief Evaluates the components from a given position.
     * @param out Vector receiving the values.
     * @param p Device or integration point.
     * @param i Position of the first component.
     */
    template <class Point>
    KOKKOS_INLINE_FUNCTION void evaluate(
      Geometry::Vector& out, const Point& p, Index i) const
    {
      out(i) = Variational::evaluate(first, p);
      rest.evaluate(out, p, i + 1);
    }

    /// @brief First component.
    View<Head> first;

    /// @brief Remaining components.
    ViewPack<Tail...> rest;
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
  struct View<::Rodin::Variational::VectorFunction<V, Values...>>
  {
    static_assert(1 + sizeof...(Values) <= 3,
      "Vector functions on Kokkos meshes have at most three components.");

    /// @brief Component views.
    using PackType = Internal::ViewPack<::Rodin::Variational::RealFunction<V>,
      ::Rodin::Variational::RealFunction<Values>...>;

    /// @param node Function expression.
    View(const ::Rodin::Variational::VectorFunction<V, Values...>& node)
      : components(std::make_from_tuple<PackType>(node.getFunctions()))
    {}

    /// @param p Device or integration point. @return Function value.
    template <class Point>
    KOKKOS_INLINE_FUNCTION Geometry::Vector getValue(const Point& p) const
    {
      Geometry::Vector result(1 + sizeof...(Values));
      components.evaluate(result, p, 0);
      return result;
    }

    /// @brief Component views.
    PackType components;
  };

}

#endif
