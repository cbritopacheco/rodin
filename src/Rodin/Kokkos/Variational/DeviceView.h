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

  /// @brief Tells whether a space is a real valued device space.
  template <class FES>
  concept ScalarDeviceSpace =
    DeviceSpace<FES> && std::is_same_v<typename FES::RangeType, Real>;

  /// @brief Tells whether a space is a vector valued device space.
  template <class FES>
  concept VectorDeviceSpace = DeviceSpace<FES> &&
    std::is_same_v<typename FES::RangeType, Math::SpatialVector<Real>>;

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
   * @brief Device vector with a single nonzero component.
   *
   * Value of a vector Lagrange basis @f$ \phi_a e_r @f$. The operations without
   * a dedicated overload use its conversion to a dense device vector.
   */
  struct UnitVector
  {
    /// @return Dense device vector.
    KOKKOS_INLINE_FUNCTION
    operator Geometry::Vector() const
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
    KOKKOS_INLINE_FUNCTION
    operator Geometry::Matrix() const
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
  KOKKOS_INLINE_FUNCTION
  UnitVector mult(Real s, const UnitVector& v)
  {
    return {v.size, v.component, s * v.value};
  }

  /// @param v Unit vector. @param s Scaling factor. @return Vector @f$ v s @f$.
  KOKKOS_INLINE_FUNCTION
  UnitVector mult(const UnitVector& v, Real s)
  {
    return mult(s, v);
  }

  /// @param s Scaling factor. @param m Row matrix. @return Matrix @f$ s m @f$.
  KOKKOS_INLINE_FUNCTION
  RowMatrix mult(Real s, const RowMatrix& m)
  {
    return {m.rows, m.row, mult(s, m.values)};
  }

  /// @param m Row matrix. @param s Scaling factor. @return Matrix @f$ m s @f$.
  KOKKOS_INLINE_FUNCTION
  RowMatrix mult(const RowMatrix& m, Real s)
  {
    return mult(s, m);
  }

  /**
   * @brief Multiplies a device vector by a row matrix.
   * @param m Row matrix.
   * @param v Device vector, with as many entries as @p m has columns.
   * @return Vector @f$ m v @f$.
   */
  KOKKOS_INLINE_FUNCTION
  UnitVector mult(const RowMatrix& m, const Geometry::Vector& v)
  {
    return {m.rows, m.row, dot(m.values, v)};
  }

  /// @param a First unit vector. @param b Second unit vector. @return @f$ a \cdot b @f$.
  KOKKOS_INLINE_FUNCTION
  Real dot(const UnitVector& a, const UnitVector& b)
  {
    KOKKOS_ASSERT(a.size == b.size);
    return a.component == b.component ? a.value * b.value : Real(0);
  }

  /// @param a Device vector. @param b Unit vector. @return @f$ a \cdot b @f$.
  KOKKOS_INLINE_FUNCTION
  Real dot(const Geometry::Vector& a, const UnitVector& b)
  {
    KOKKOS_ASSERT(a.size() == b.size);
    return a(b.component) * b.value;
  }

  /// @param a Unit vector. @param b Device vector. @return @f$ a \cdot b @f$.
  KOKKOS_INLINE_FUNCTION
  Real dot(const UnitVector& a, const Geometry::Vector& b)
  {
    return dot(b, a);
  }

  /// @param a First row matrix. @param b Second row matrix. @return @f$ a : b @f$.
  KOKKOS_INLINE_FUNCTION
  Real dot(const RowMatrix& a, const RowMatrix& b)
  {
    KOKKOS_ASSERT(a.rows == b.rows);
    return a.row == b.row ? dot(a.values, b.values) : Real(0);
  }

  /// @param a Device matrix. @param b Row matrix. @return @f$ a : b @f$.
  KOKKOS_INLINE_FUNCTION
  Real dot(const Geometry::Matrix& a, const RowMatrix& b)
  {
    KOKKOS_ASSERT(a.rows() == b.rows && a.cols() == b.values.size());
    Real result = 0;
    for (Index c = 0; c < b.values.size(); ++c)
      result += a(b.row, c) * b.values(c);
    return result;
  }

  /// @param a Row matrix. @param b Device matrix. @return @f$ a : b @f$.
  KOKKOS_INLINE_FUNCTION
  Real dot(const RowMatrix& a, const Geometry::Matrix& b)
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
  KOKKOS_INLINE_FUNCTION
  bool overlaps(const UnitVector& a, const UnitVector& b)
  {
    return a.component == b.component;
  }

  /// @param a First row matrix. @param b Second row matrix. @return Whether the rows match.
  KOKKOS_INLINE_FUNCTION
  bool overlaps(const RowMatrix& a, const RowMatrix& b)
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
  KOKKOS_INLINE_FUNCTION
  Real getPayload(const UnitVector& v)
  {
    return v.value;
  }

  /// @param m Row matrix. @return Entries of its nonzero row.
  KOKKOS_INLINE_FUNCTION
  Geometry::Vector getPayload(const RowMatrix& m)
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

  /// @brief Tag of polytopes whose geometry is known at run time only.
  struct DynamicGeometry
  {};

  /**
   * @brief Tag of polytopes of geometry G whose dimension is the spatial dimension.
   * @tparam G Reference geometry, as supported by Geometry::StaticBasis.
   */
  template <::Rodin::Geometry::Polytope::Type G>
  struct StaticGeometry
  {
    /// @brief Reference geometry.
    static constexpr ::Rodin::Geometry::Polytope::Type Value = G;
  };

  /**
   * @brief P1 basis values and physical gradients at one quadrature point.
   *
   * The point is evaluated once per quadrature point and shared by every
   * shape function of the integrand. Gradients are mapped to physical space
   * by @f$ \nabla \phi = J^{-T} \hat{\nabla} \hat{\phi} @f$, as on the host.
   *
   * A StaticGeometry tag evaluates the point with static bounds, getting the
   * physical point and the Jacobian in one pass over the vertices.
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
        evaluate(point);
      }

      /**
       * @brief Evaluates the P1 basis at a reference point of a polytope.
       * @param polytope Device polytope.
       * @param reference Reference coordinates.
       * @param geometry Tag of a geometry known at run time only.
       */
      KOKKOS_INLINE_FUNCTION
      IntegrationPointView(const Geometry::PolytopeView& polytope,
        const Geometry::Vector& reference, [[maybe_unused]] DynamicGeometry geometry)
        : IntegrationPointView(polytope.getPoint(reference))
      {}

      /**
       * @brief Evaluates the P1 basis at a reference point of a polytope of geometry G.
       * @param polytope Device polytope of geometry G, in a mesh of spatial
       * dimension the dimension of G.
       * @param reference Reference coordinates.
       * @param geometry Tag of the geometry G.
       */
      template <::Rodin::Geometry::Polytope::Type G>
      KOKKOS_INLINE_FUNCTION IntegrationPointView(const Geometry::PolytopeView& polytope,
        const Geometry::Vector& reference, [[maybe_unused]] StaticGeometry<G> geometry)
        : m_point(polytope, reference, reference),
          m_count(Geometry::StaticBasis<G>::Count)
      {
        KOKKOS_ASSERT(polytope.getGeometry() == G);
        KOKKOS_ASSERT(
          polytope.mesh.spatialDimension == Geometry::StaticBasis<G>::Dimension);
        evaluate<G, true>(polytope, reference);
      }

      /**
       * @brief Moves the integration point within its polytope.
       *
       * The P1 map of a simplex is affine, so its Jacobian, distortion and
       * basis gradients are kept and only the basis values are evaluated
       * again. Other polytopes are evaluated again in full.
       *
       * @param point Point on the polytope of the current point.
       */
      KOKKOS_INLINE_FUNCTION
      void setPoint(const Geometry::PointView& point)
      {
        using Type = ::Rodin::Geometry::Polytope::Type;
        const Type type = point.getPolytope().getGeometry();
        m_point = point;
        if (type != Type::Point && type != Type::Segment && type != Type::Triangle &&
          type != Type::Tetrahedron)
        {
          evaluate(point);
          return;
        }
        Geometry::Vector reference(point.getReferenceCoordinates().size());
        for (Index a = 0; a < m_count; ++a)
          point.evaluateBasis(a, m_basis[a], reference);
      }

      /**
       * @brief Moves the integration point within its polytope.
       * @param reference Reference coordinates in the current polytope.
       * @param geometry Tag of a geometry known at run time only.
       */
      KOKKOS_INLINE_FUNCTION
      void setPoint(
        const Geometry::Vector& reference, [[maybe_unused]] DynamicGeometry geometry)
      {
        setPoint(m_point.getPolytope().getPoint(reference));
      }

      /**
       * @brief Moves the integration point within its polytope of geometry G.
       *
       * As for a dynamic geometry, a simplex keeps its Jacobian, distortion
       * and basis gradients.
       *
       * @param reference Reference coordinates in the current polytope.
       * @param geometry Tag of the geometry G.
       */
      template <::Rodin::Geometry::Polytope::Type G>
      KOKKOS_INLINE_FUNCTION void setPoint(
        const Geometry::Vector& reference, [[maybe_unused]] StaticGeometry<G> geometry)
      {
        const Geometry::PolytopeView polytope = m_point.getPolytope();
        evaluate<G, !Geometry::StaticBasis<G>::Affine>(polytope, reference);
      }

      /// @brief Returns the evaluation point. @return Device point.
      KOKKOS_INLINE_FUNCTION
      const Geometry::PointView& getPoint() const
      {
        return m_point;
      }

      /// @brief Returns the integration measure. @return Distortion of the point.
      KOKKOS_INLINE_FUNCTION
      Real getDistortion() const
      {
        return m_distortion;
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
      /**
       * @brief Evaluates the basis, the gradients and the distortion at a point of the polytope.
       * @param point Point on a mesh entity.
       */
      KOKKOS_INLINE_FUNCTION
      void evaluate(const Geometry::PointView& point)
      {
        // The inverse and the measure of a square map share one Jacobian.
        const Geometry::Matrix jacobian = point.getJacobian();
        Geometry::Matrix inverse;
        if (jacobian.cols() > 0 && jacobian.rows() == jacobian.cols())
        {
          inverse = jacobian.inverse();
          m_distortion = ::Kokkos::abs(jacobian.determinant());
        }
        else
        {
          inverse = point.getJacobianInverse();
          m_distortion = point.getDistortion();
        }
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

      /**
       * @brief Evaluates the point on a polytope of geometry G.
       * @tparam G Reference geometry.
       * @tparam Full Whether to evaluate the Jacobian, the distortion and the
       * gradients, or only the basis values and the physical point.
       * @param polytope Polytope of geometry G filling its space.
       * @param reference Reference coordinates.
       */
      template <::Rodin::Geometry::Polytope::Type G, bool Full>
      KOKKOS_INLINE_FUNCTION void evaluate(
        const Geometry::PolytopeView& polytope, const Geometry::Vector& reference)
      {
        using Basis = Geometry::StaticBasis<G>;
        constexpr Index D = Basis::Dimension;
        constexpr Index N = Basis::Count;
        ::Kokkos::Array<::Kokkos::Array<Real, D>, N> gradients;
        if constexpr (Full)
          Basis::evaluate(reference, m_basis, gradients);
        else
          Basis::evaluate(reference, m_basis);
        const auto& mesh = polytope.mesh;
        const Index entity = mesh.entityOffsets[polytope.dimension] + polytope.index;
        const Geometry::GeometryKind kind = mesh.geometryKinds[entity];
        Geometry::Vector physical(D);
        Geometry::Matrix jacobian(D, D);
        if (kind == Geometry::GeometryKind::Identity)
        {
          physical = reference;
          for (Index d = 0; d < D; ++d)
            jacobian(d, d) = 1;
        }
        else
        {
          const bool vertexMap = kind == Geometry::GeometryKind::VertexP1;
          const Index begin =
            vertexMap ? mesh.vertexOffsets[entity] : mesh.controlOffsets[entity];
          for (Index a = 0; a < N; ++a)
          {
            const Real* control = vertexMap
              ? mesh.coordinates + mesh.vertexIndices[begin + a] * D
              : mesh.controls + (begin + a) * D;
            for (Index c = 0; c < D; ++c)
            {
              physical(c) += control[c] * m_basis[a];
              if constexpr (Full)
              {
                for (Index d = 0; d < D; ++d)
                  jacobian(c, d) += control[c] * gradients[a][d];
              }
            }
          }
        }
        m_point = Geometry::PointView(polytope, reference, physical);
        if constexpr (Full)
        {
          const Geometry::Matrix inverse = jacobian.inverse();
          m_distortion = ::Kokkos::abs(jacobian.determinant());
          for (Index a = 0; a < N; ++a)
          {
            Geometry::Vector gradient(D);
            for (Index c = 0; c < D; ++c)
            {
              for (Index k = 0; k < D; ++k)
                gradient(c) += inverse(k, c) * gradients[a][k];
            }
            m_gradients[a] = gradient;
          }
        }
      }

      Geometry::PointView m_point;
      Real m_distortion;
      Index m_count;
      ::Kokkos::Array<Real, Capacity> m_basis;
      ::Kokkos::Array<Geometry::Vector, Capacity> m_gradients;
  };

  /**
   * @brief Integration point known at compile time to lie in a cell of the
   * mesh.
   *
   * The kernels of a static geometry integrate polytopes filling the space,
   * hence cells, and build their points with this type. The derivatives of
   * grid functions then leave out the trace on interfaces, whose registers
   * would lower the occupancy of the kernel.
   */
  class CellIntegrationPointView : public IntegrationPointView
  {
    public:
      using IntegrationPointView::IntegrationPointView;
  };

  /**
   * @brief Device-callable view of a form language expression.
   *
   * A view is built on the host from the expression and copied by value into
   * kernels. It holds no host pointer and no virtual function. Function views
   * provide `getValue(const Geometry::PointView&)`, while shape function views
   * provide `getBasis(const IntegrationPointView&, Index, Index)`, which
   * returns the value of component r of the basis function of the local scalar
   * degree of freedom a. A scalar view has the single component 0.
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
   * | `DeviceView<Dot<FunctionBase<L>, FunctionBase<R>>>` | Dot product of functions. |
   * | `DeviceView<UnaryMinus<FunctionBase<N>>>` | Negated function. |
   * | `DeviceView<RealFunction<Integer>>` | Integer constant. |
   * | `DeviceView<VectorFunctionBase<Scalar, Derived>>` | Forwards to the derived function. |
   * | `DeviceView<MatrixFunctionBase<Scalar, Derived>>` | Forwards to the derived function. |
   * | `DeviceView<MatrixFunction<Matrix>>` | Constant matrix. |
   * | `DeviceView<GradBase<Operand, Derived>>` | Forwards to the derived function. |
   * | `DeviceView<JacobianBase<Operand, Derived>>` | Forwards to the derived function. |
   * | `DeviceView<DivBase<Operand, Derived>>` | Forwards to the derived function. |
   * | `DeviceView<VectorFunction<V, Values...>>` | Vector function given component-wise. |
   * | `DeviceView<Zero<Real>>` | Scalar zero. |
   * | `DeviceView<Zero<SpatialVector<Real>>>` | Vector zero. |
   * | `DeviceView<BoundaryNormal>` | Outward unit normal of the boundary faces. |
   * | `DeviceView<GridFunctionBaseReference<GridFunction<FES, Data>>>` | Referenced grid function of a device space. |
   * | `DeviceView<GridFunction<FES, Data>>` | Grid function of a device space. |
   * | `DeviceView<Grad<GridFunction<FES, Data>>>` | Physical gradient of a scalar grid function. |
   * | `DeviceView<Jacobian<GridFunction<FES, Data>>>` | Physical Jacobian of a vector grid function. |
   * | `DeviceView<Div<GridFunction<FES, Data>>>` | Divergence of a vector grid function. |
   * | `DeviceView<ShapeFunctionBase<Derived, FES, Space>>` | Forwards to the derived shape function. |
   * | `DeviceView<ShapeFunction<Derived, FES, Space>>` | Scalar or vector basis value of a device space. |
   * | `DeviceView<Grad<ShapeFunction<Derived, FES, Space>>>` | Physical gradient of a scalar basis. |
   * | `DeviceView<Jacobian<ShapeFunction<Derived, FES, Space>>>` | Physical Jacobian of a vector basis. |
   * | `DeviceView<Div<ShapeFunction<Derived, FES, Space>>>` | Divergence of a vector basis. |
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
      KOKKOS_INLINE_FUNCTION
      decltype(auto) getBasis(Index dof, Index component) const
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
      KOKKOS_INLINE_FUNCTION
      auto getBasis(Index dof, Index component) const
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
      KOKKOS_INLINE_FUNCTION
      auto getBasis(Index dof, Index component) const
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
      KOKKOS_INLINE_FUNCTION
      Real getBasis(Index dof, Index component) const
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
      KOKKOS_INLINE_FUNCTION
      auto getBasis(Index dof, Index component) const
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
      KOKKOS_INLINE_FUNCTION
      auto getBasis(Index dof, Index component) const
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
      KOKKOS_INLINE_FUNCTION
      Geometry::Matrix getBasis(Index dof, Index component) const
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

  /// @brief Real function wrapping a real function, as a vector component.
  template <class F>
  struct DeviceView<
    ::Rodin::Variational::RealFunction<::Rodin::Variational::RealFunction<F>>>
    : DeviceView<::Rodin::Variational::RealFunction<F>>
  {
    /// @param node Function expression.
    DeviceView(
      const ::Rodin::Variational::RealFunction<::Rodin::Variational::RealFunction<F>>&
        node)
      : DeviceView<::Rodin::Variational::RealFunction<F>>(node.getFunction())
    {}
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

    /// @param p Device or integration point. @return Function value.
    template <class Point>
    KOKKOS_INLINE_FUNCTION auto getValue(const Point& p) const
    {
      return add(evaluate(lhs, p), evaluate(rhs, p));
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

    /// @param p Device or integration point. @return Function value.
    template <class Point>
    KOKKOS_INLINE_FUNCTION auto getValue(const Point& p) const
    {
      return mult(evaluate(lhs, p), evaluate(rhs, p));
    }

    /// @brief Left operand.
    DeviceView<::Rodin::Variational::FunctionBase<LHSDerived>> lhs;

    /// @brief Right operand.
    DeviceView<::Rodin::Variational::FunctionBase<RHSDerived>> rhs;
  };

  /// @brief Dot product of two functions.
  template <class LHSDerived, class RHSDerived>
  struct DeviceView<
    ::Rodin::Variational::Dot<::Rodin::Variational::FunctionBase<LHSDerived>,
      ::Rodin::Variational::FunctionBase<RHSDerived>>>
  {
    /// @param node Function expression.
    DeviceView(
      const ::Rodin::Variational::Dot<::Rodin::Variational::FunctionBase<LHSDerived>,
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

    /// @param p Device or integration point. @return Function value.
    template <class Point>
    KOKKOS_INLINE_FUNCTION auto getValue(const Point& p) const
    {
      return mult(Real(-1), evaluate(operand, p));
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

  /// @brief Constant matrix of at most three rows and three columns.
  template <class Matrix>
    requires std::is_same_v<Matrix, Math::Matrix<Real>> ||
    std::is_same_v<Matrix, Math::SpatialMatrix<Real>>
  struct DeviceView<::Rodin::Variational::MatrixFunction<Matrix>>
  {
    /// @param node Function expression.
    DeviceView(const ::Rodin::Variational::MatrixFunction<Matrix>& node)
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

  /// @brief Scalar zero.
  template <>
  struct DeviceView<::Rodin::Variational::Zero<Real>>
  {
    /// @param node Function expression, unused.
    DeviceView([[maybe_unused]] const ::Rodin::Variational::Zero<Real>& node) {}

    /// @param p Device point, unused. @return Zero.
    KOKKOS_INLINE_FUNCTION
    Real getValue([[maybe_unused]] const Geometry::PointView& p) const
    {
      return 0;
    }
  };

  /// @brief Vector zero.
  template <>
  struct DeviceView<::Rodin::Variational::Zero<Math::SpatialVector<Real>>>
  {
    /// @param node Function expression.
    DeviceView(const ::Rodin::Variational::Zero<Math::SpatialVector<Real>>& node)
      : dimension(node.getDimension())
    {}

    /// @param p Device point, unused. @return Zero vector.
    KOKKOS_INLINE_FUNCTION
    Geometry::Vector getValue([[maybe_unused]] const Geometry::PointView& p) const
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
  struct DeviceView<::Rodin::Variational::BoundaryNormal>
  {
    /// @param node Function expression.
    DeviceView([[maybe_unused]] const ::Rodin::Variational::BoundaryNormal& node) {}

    /// @param p Point on a boundary face. @return Outward unit normal.
    KOKKOS_INLINE_FUNCTION
    Geometry::Vector getValue(const Geometry::PointView& p) const
    {
      const auto& face = p.getPolytope();
      const Geometry::MeshView& mesh = face.mesh;
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
  struct DeviceViewPack
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

  /**
   * @brief Trace domain of a function, in device memory.
   *
   * A derivative of a grid function evaluated on a face is that of an
   * incident cell, as on the host. A boundary face has a single incident
   * cell. On an interface, the cell is the one whose attribute is in the
   * trace domain, and a kernel evaluating it with no such cell aborts.
   */
  class TraceDomainView
  {
    public:
      /// @brief Maximum number of attributes.
      static constexpr Index Capacity = 16;

      /// @brief Constructs the empty trace domain.
      TraceDomainView() = default;

      /**
       * @brief Copies a trace domain.
       * @param domain Attributes of the cells giving the trace.
       */
      explicit TraceDomainView(const FlatSet<::Rodin::Geometry::Attribute>& domain)
        : m_count(domain.size())
      {
        if (m_count > Capacity)
        {
          Alert::Exception() << "Device trace domains hold at most " << Capacity
                             << " attributes." << Alert::Raise;
        }
        Index i = 0;
        for (const auto& attribute : domain)
          m_attributes[i++] = attribute;
      }

      /**
       * @brief Tests whether a point lies on a face.
       * @param p Device point. @return Whether the point is bound to a face.
       */
      KOKKOS_INLINE_FUNCTION
      static bool isFace(const Geometry::PointView& p)
      {
        const auto& polytope = p.getPolytope();
        return polytope.getDimension() + 1 == polytope.mesh.getDimension();
      }

      /**
       * @brief Carries a point of a face to the cell giving the trace.
       * @param p Point on a face; face-to-cell incidence must be computed.
       * @return Same physical point, bound to the cell.
       */
      KOKKOS_INLINE_FUNCTION
      Geometry::PointView getCellPoint(const Geometry::PointView& p) const
      {
        const auto& face = p.getPolytope();
        const Index d = face.mesh.getDimension();
        const auto connectivity = face.mesh.getConnectivity();
        KOKKOS_ASSERT(connectivity.hasIncidence(d - 1, d));
        const auto cells = connectivity.getIncidence(d - 1, d, face.getIndex());
        if (cells.size() == 1)
          return p.getCellPoint(cells[0]);
        for (Index c = 0; c < cells.size(); ++c)
        {
          const auto cell = face.mesh.getCell(cells[c]);
          if (!cell.hasAttribute())
            continue;
          for (Index i = 0; i < m_count; ++i)
          {
            if (m_attributes[i] == cell.getAttribute())
              return p.getCellPoint(cells[c]);
          }
        }
        ::Kokkos::abort(
          "Rodin: derivative at an interface with no cell in the trace domain.");
        return p;
      }

    private:
      Index m_count = 0;
      ::Kokkos::Array<::Rodin::Geometry::Attribute, Capacity> m_attributes = {};
  };

  /**
   * @brief Copies the coefficients of a grid function to the memory space
   * of its mesh, in the local numbering of the device degrees of freedom.
   *
   * This primary template reads containers indexed by coeff(), such as
   * Math::Vector. Modules whose grid functions keep their coefficients in
   * another container specialize it, as Rodin/Kokkos/Assembly/PETSc.h does
   * for PETSc vectors.
   *
   * @tparam Data Coefficient container of the grid function.
   */
  template <class Data>
  struct GridFunctionCoefficients
  {
    /**
     * @param fes Device space of the grid function.
     * @param data Coefficients, numbered as by the space.
     * @returns Device coefficients, after transfer completion.
     */
    template <class FES>
    static ::Kokkos::View<Real*, typename FES::MeshType::MemorySpace> upload(
      const FES& fes, const Data& data)
    {
      const auto& mesh = fes.getMesh();
      const auto& execution = mesh.getContext().getExecutionSpace();
      const Index size = data.size();
      if constexpr (requires { mesh.getTransferCache(); })
      {
        // Lowered at every assembly, so that the same coefficients stay
        // resident instead of being reallocated and copied each time.
        std::vector<Real> values(size);
        for (Index i = 0; i < size; ++i)
          values[i] = data.coeff(i);
        return mesh.getTransferCache().upload(
          execution, std::move(values), "Rodin grid function");
      }
      else
      {
        ::Kokkos::View<Real*, typename FES::MeshType::MemorySpace> res(
          ::Kokkos::view_alloc(execution, std::string("Rodin grid function")), size);
        auto host = ::Kokkos::create_mirror_view(res);
        for (Index i = 0; i < size; ++i)
          host(i) = data.coeff(i);
        ::Kokkos::deep_copy(execution, res, host);
        execution.fence("Rodin grid function transfer");
        return res;
      }
    }
  };

  /**
   * @brief Coefficients of a grid function of a device space in device
   * memory.
   *
   * The coefficients are copied to the mesh memory space when the view is
   * built, so the view reflects the grid function at lowering time. The
   * copy comes from the transfer cache of the mesh, which keeps it
   * resident while the coefficients do not change. Values and
   * derivatives are interpolated from the degrees of freedom of the polytope
   * holding the point, as on the host.
   *
   * @tparam FES Device space of the grid function.
   */
  template <DeviceSpace FES>
  class GridFunctionView
  {
    public:
      /// @brief Finite element space.
      using FESType = FES;

      /// @brief Device degrees of freedom.
      using DOFsType = typename FESType::DOFsType;

      /// @brief Memory space of the coefficients.
      using MemorySpace = typename FESType::MeshType::MemorySpace;

      /**
       * @brief Copies the coefficients of a grid function to the device.
       * @param fes Space of the grid function.
       * @param data Coefficients, numbered as by the space.
       */
      template <class Data>
      GridFunctionView(const FESType& fes, const Data& data)
        : m_vdim(fes.getVectorDimension()),
          m_dofs(fes.getDeviceDOFs()),
          m_data(GridFunctionCoefficients<std::decay_t<Data>>::upload(fes, data))
      {}

      /// @brief Gets the device degrees of freedom. @return Degrees of freedom.
      KOKKOS_INLINE_FUNCTION
      const DOFsType& getDOFs() const
      {
        return m_dofs;
      }

      /// @brief Gets the coefficients. @return Device coefficients.
      KOKKOS_INLINE_FUNCTION
      const ::Kokkos::View<Real*, MemorySpace>& getData() const
      {
        return m_data;
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
       *
       * The P1 basis is evaluated directly at the point, the other bases
       * through an IntegrationPointView of the point.
       */
      KOKKOS_INLINE_FUNCTION
      Real interpolate(const Geometry::PointView& p, Index r) const
      {
        if constexpr (std::is_same_v<DOFsType, VertexDOFs>)
        {
          const auto vertices = p.getPolytope().getVertices();
          Geometry::Vector reference(p.getReferenceCoordinates().size());
          Real result = 0;
          for (Index a = 0; a < vertices.size(); ++a)
          {
            Real phi;
            p.evaluateBasis(a, phi, reference);
            result += phi * m_data(vertices[a] + r * m_dofs.getStride());
          }
          return result;
        }
        else
        {
          return interpolate(IntegrationPointView(p), r);
        }
      }

      /**
       * @brief Interpolates the physical Jacobian.
       * @param p Device point.
       * @param trace Trace domain, used when the point lies on an interface.
       * @return Matrix with one row per component, row @f$ r @f$ holding
       * @f$ \nabla u_r(p) @f$.
       *
       * On a face, the Jacobian is that of an incident cell.
       */
      KOKKOS_INLINE_FUNCTION
      Geometry::Matrix interpolateJacobian(
        const Geometry::PointView& p, const TraceDomainView& trace = {}) const
      {
        return TraceDomainView::isFace(p) ? getCellJacobian(trace.getCellPoint(p))
                                          : getCellJacobian(p);
      }

      /**
       * @brief Interpolates the divergence.
       * @param p Device point.
       * @param trace Trace domain, used when the point lies on an interface.
       * @return Trace of interpolateJacobian(p), summed in the same order.
       */
      KOKKOS_INLINE_FUNCTION
      Real interpolateDivergence(
        const Geometry::PointView& p, const TraceDomainView& trace = {}) const
      {
        return TraceDomainView::isFace(p) ? getCellDivergence(trace.getCellPoint(p))
                                          : getCellDivergence(p);
      }

      /**
       * @brief Interpolates one component with the basis of an integration
       * point.
       * @param ip Integration point on a polytope of the space mesh.
       * @param r Component index.
       * @return Component @f$ u_r(p) @f$.
       */
      KOKKOS_INLINE_FUNCTION
      Real interpolate(const IntegrationPointView& ip, Index r) const
      {
        const auto dofs = m_dofs.bind(ip.getPoint().getPolytope());
        Real result = 0;
        forEachDOF(dofs.size(), [&](Index q) {
          result += m_dofs.getBasis(ip, q) * m_data(dofs[q] + r * m_dofs.getStride());
        });
        return result;
      }

      /**
       * @brief Interpolates every component with the basis of an integration
       * point.
       * @param ip Integration point on a polytope of the space mesh.
       * @return Vector of the components @f$ u_r(p) @f$.
       */
      KOKKOS_INLINE_FUNCTION
      Geometry::Vector interpolateVector(const IntegrationPointView& ip) const
      {
        const auto dofs = m_dofs.bind(ip.getPoint().getPolytope());
        Geometry::Vector result(m_vdim);
        forEachDOF(dofs.size(), [&](Index q) {
          const Real basis = m_dofs.getBasis(ip, q);
          forEachComponent(m_vdim, [&](Index r) {
            result(r) += basis * m_data(dofs[q] + r * m_dofs.getStride());
          });
        });
        return result;
      }

      /**
       * @brief Interpolates every component.
       * @param p Device point.
       * @return Vector of the components @f$ u_r(p) @f$.
       */
      KOKKOS_INLINE_FUNCTION
      Geometry::Vector interpolateVector(const Geometry::PointView& p) const
      {
        Geometry::Vector result(m_vdim);
        for (Index r = 0; r < m_vdim; ++r)
          result(r) = interpolate(p, r);
        return result;
      }

      /**
       * @brief Interpolates the physical Jacobian with the basis gradients of
       * an integration point.
       * @param ip Integration point on a polytope of the space mesh.
       * @param trace Trace domain, used when the point lies on an interface.
       * @return Matrix with one row per component, row @f$ r @f$ holding
       * @f$ \nabla u_r(p) @f$.
       */
      KOKKOS_INLINE_FUNCTION
      Geometry::Matrix interpolateJacobian(
        const IntegrationPointView& ip, const TraceDomainView& trace = {}) const
      {
        return TraceDomainView::isFace(ip.getPoint())
          ? getCellJacobian(trace.getCellPoint(ip.getPoint()))
          : getCellJacobian(ip);
      }

      /**
       * @brief Interpolates the divergence with the basis gradients of an
       * integration point.
       * @param ip Integration point on a polytope of the space mesh.
       * @param trace Trace domain, used when the point lies on an interface.
       * @return Trace of interpolateJacobian(ip), summed in the same order.
       */
      KOKKOS_INLINE_FUNCTION
      Real interpolateDivergence(
        const IntegrationPointView& ip, const TraceDomainView& trace = {}) const
      {
        return TraceDomainView::isFace(ip.getPoint())
          ? getCellDivergence(trace.getCellPoint(ip.getPoint()))
          : getCellDivergence(ip);
      }

      /**
       * @brief Interpolates the physical Jacobian at an integration point of
       * a cell.
       * @param ip Integration point on a cell of the space mesh.
       * @param trace Trace domain, unused since a cell is not an interface.
       * @return Matrix with one row per component, row @f$ r @f$ holding
       * @f$ \nabla u_r(p) @f$.
       */
      KOKKOS_INLINE_FUNCTION
      Geometry::Matrix interpolateJacobian(const CellIntegrationPointView& ip,
        [[maybe_unused]] const TraceDomainView& trace = {}) const
      {
        return getCellJacobian(ip);
      }

      /**
       * @brief Interpolates the divergence at an integration point of a cell.
       * @param ip Integration point on a cell of the space mesh.
       * @param trace Trace domain, unused since a cell is not an interface.
       * @return Trace of interpolateJacobian(ip), summed in the same order.
       */
      KOKKOS_INLINE_FUNCTION
      Real interpolateDivergence(const CellIntegrationPointView& ip,
        [[maybe_unused]] const TraceDomainView& trace = {}) const
      {
        return getCellDivergence(ip);
      }

    private:
      /**
       * @brief Calls a function on each local degree of freedom.
       * @param count Local degree of freedom count.
       * @param f Function of the local index.
       *
       * The P1 loop has a static bound, so that the basis cached by the
       * integration point is indexed statically and stays in registers.
       */
      template <class F>
      KOKKOS_INLINE_FUNCTION static void forEachDOF(Index count, const F& f)
      {
        if constexpr (std::is_same_v<DOFsType, VertexDOFs>)
        {
          RODIN_KOKKOS_UNROLL
          for (Index q = 0; q < IntegrationPointView::Capacity; ++q)
          {
            if (q < count)
              f(q);
          }
        }
        else
        {
          for (Index q = 0; q < count; ++q)
            f(q);
        }
      }

      /**
       * @brief Calls a function on each component, or each coordinate.
       * @param count Component count, at most three.
       * @param f Function of the index.
       *
       * The static bound keeps the device vectors and matrices in registers.
       */
      template <class F>
      KOKKOS_INLINE_FUNCTION static void forEachComponent(Index count, const F& f)
      {
        KOKKOS_ASSERT(count <= 3);
        RODIN_KOKKOS_UNROLL
        for (Index r = 0; r < 3; ++r)
        {
          if (r < count)
            f(r);
        }
      }

      /**
       * @brief Interpolates the physical Jacobian at a point of a cell.
       *
       * The P1 basis is evaluated directly at the point, the other bases
       * through an IntegrationPointView of the point.
       *
       * @param p Point on a cell of the space mesh.
       * @return Matrix with one row per component, row @f$ r @f$ holding
       * @f$ \nabla u_r(p) @f$.
       */
      KOKKOS_INLINE_FUNCTION
      Geometry::Matrix getCellJacobian(const Geometry::PointView& p) const
      {
        if constexpr (std::is_same_v<DOFsType, VertexDOFs>)
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
                result(r, c) += gradient * m_data(vertices[a] + r * m_dofs.getStride());
            }
          }
          return result;
        }
        else
        {
          return getCellJacobian(IntegrationPointView(p));
        }
      }

      /**
       * @brief Interpolates the divergence at a point of a cell.
       * @param p Point on a cell of the space mesh.
       * @return Trace of getCellJacobian(p).
       */
      KOKKOS_INLINE_FUNCTION
      Real getCellDivergence(const Geometry::PointView& p) const
      {
        if constexpr (std::is_same_v<DOFsType, VertexDOFs>)
        {
          const auto vertices = p.getPolytope().getVertices();
          const Geometry::Matrix inverse = p.getJacobianInverse();
          const Index rdim = inverse.rows();
          KOKKOS_ASSERT(m_vdim == inverse.cols());
          Real result = 0;
          for (Index r = 0; r < m_vdim; ++r)
          {
            Real diagonal = 0;
            for (Index a = 0; a < vertices.size(); ++a)
            {
              Real phi;
              Geometry::Vector reference(rdim);
              p.evaluateBasis(a, phi, reference);
              Real gradient = 0;
              for (Index k = 0; k < rdim; ++k)
                gradient += inverse(k, r) * reference(k);
              diagonal += gradient * m_data(vertices[a] + r * m_dofs.getStride());
            }
            result += diagonal;
          }
          return result;
        }
        else
        {
          return getCellDivergence(IntegrationPointView(p));
        }
      }

      /**
       * @brief Interpolates the physical Jacobian at an integration point of a cell.
       * @param ip Integration point on a cell of the space mesh.
       * @return Matrix with one row per component, row @f$ r @f$ holding
       * @f$ \nabla u_r(p) @f$.
       */
      KOKKOS_INLINE_FUNCTION
      Geometry::Matrix getCellJacobian(const IntegrationPointView& ip) const
      {
        const auto dofs = m_dofs.bind(ip.getPoint().getPolytope());
        const Index sdim = ip.getGradient(0).size();
        Geometry::Matrix result(m_vdim, sdim);
        forEachDOF(dofs.size(), [&](Index q) {
          const auto& gradient = m_dofs.getGradient(ip, q);
          forEachComponent(m_vdim, [&](Index r) {
            const Real value = m_data(dofs[q] + r * m_dofs.getStride());
            forEachComponent(sdim, [&](Index c) { result(r, c) += gradient(c) * value; });
          });
        });
        return result;
      }

      /**
       * @brief Interpolates the divergence at an integration point of a cell.
       * @param ip Integration point on a cell of the space mesh.
       * @return Trace of getCellJacobian(ip).
       */
      KOKKOS_INLINE_FUNCTION
      Real getCellDivergence(const IntegrationPointView& ip) const
      {
        const auto dofs = m_dofs.bind(ip.getPoint().getPolytope());
        KOKKOS_ASSERT(m_vdim == ip.getGradient(0).size());
        Real result = 0;
        forEachComponent(m_vdim, [&](Index r) {
          Real diagonal = 0;
          forEachDOF(dofs.size(), [&](Index q) {
            diagonal +=
              m_dofs.getGradient(ip, q)(r) * m_data(dofs[q] + r * m_dofs.getStride());
          });
          result += diagonal;
        });
        return result;
      }
      Index m_vdim;
      DOFsType m_dofs;
      ::Kokkos::View<Real*, MemorySpace> m_data;
  };

  /// @brief Grid function of a device space.
  template <DeviceSpace FES, class Data>
  struct DeviceView<::Rodin::Variational::GridFunction<FES, Data>>
  {
    /// @param node Grid function.
    DeviceView(const ::Rodin::Variational::GridFunction<FES, Data>& node)
      : u(node.getFiniteElementSpace(), node.getData())
    {}

    /// @param fes Space of the grid function. @param data Coefficients.
    DeviceView(const FES& fes, const Data& data)
      : u(fes, data)
    {}

    /// @param p Device or integration point. @return Function value.
    template <class Point>
    KOKKOS_INLINE_FUNCTION auto getValue(const Point& p) const
    {
      if constexpr (std::is_same_v<typename FES::RangeType, Real>)
      {
        return u.interpolate(p, 0);
      }
      else
      {
        return u.interpolateVector(p);
      }
    }

    /// @brief Device coefficients.
    GridFunctionView<FES> u;
  };

  /// @brief Grid function referenced by an expression.
  template <DeviceSpace FES, class Data>
  struct DeviceView<::Rodin::Variational::GridFunctionBaseReference<
    ::Rodin::Variational::GridFunction<FES, Data>>>
    : DeviceView<::Rodin::Variational::GridFunction<FES, Data>>
  {
    /// @param node Reference to a grid function.
    DeviceView(const ::Rodin::Variational::GridFunctionBaseReference<
      ::Rodin::Variational::GridFunction<FES, Data>>& node)
      : DeviceView<::Rodin::Variational::GridFunction<FES, Data>>(
          node.getFiniteElementSpace(), node.getData())
    {}
  };

  /// @brief Physical gradient of a scalar grid function.
  template <ScalarDeviceSpace FES, class Data>
  struct DeviceView<
    ::Rodin::Variational::Grad<::Rodin::Variational::GridFunction<FES, Data>>>
  {
    /// @param node Gradient expression.
    DeviceView(
      const ::Rodin::Variational::Grad<::Rodin::Variational::GridFunction<FES, Data>>&
        node)
      : u(node.getOperand().getFiniteElementSpace(), node.getOperand().getData()),
        trace(node.getTraceDomain())
    {}

    /// @param p Device or integration point. @return Gradient @f$ \nabla u(p) @f$.
    template <class Point>
    KOKKOS_INLINE_FUNCTION Geometry::Vector getValue(const Point& p) const
    {
      const Geometry::Matrix jacobian = u.interpolateJacobian(p, trace);
      Geometry::Vector result(jacobian.cols());
      for (Index c = 0; c < jacobian.cols(); ++c)
        result(c) = jacobian(0, c);
      return result;
    }

    /// @brief Device coefficients.
    GridFunctionView<FES> u;

    /// @brief Trace domain of the derivative on interfaces.
    TraceDomainView trace;
  };

  /// @brief Physical Jacobian of a vector grid function.
  template <VectorDeviceSpace FES, class Data>
  struct DeviceView<
    ::Rodin::Variational::Jacobian<::Rodin::Variational::GridFunction<FES, Data>>>
  {
    /// @param node Jacobian expression.
    DeviceView(const ::Rodin::Variational::Jacobian<
      ::Rodin::Variational::GridFunction<FES, Data>>& node)
      : u(node.getOperand().getFiniteElementSpace(), node.getOperand().getData()),
        trace(node.getTraceDomain())
    {}

    /// @param p Device or integration point. @return Jacobian @f$ \nabla u(p) @f$.
    template <class Point>
    KOKKOS_INLINE_FUNCTION Geometry::Matrix getValue(const Point& p) const
    {
      return u.interpolateJacobian(p, trace);
    }

    /// @brief Device coefficients.
    GridFunctionView<FES> u;

    /// @brief Trace domain of the derivative on interfaces.
    TraceDomainView trace;
  };

  /// @brief Divergence of a vector grid function.
  template <VectorDeviceSpace FES, class Data>
  struct DeviceView<
    ::Rodin::Variational::Div<::Rodin::Variational::GridFunction<FES, Data>>>
  {
    /// @param node Divergence expression.
    DeviceView(
      const ::Rodin::Variational::Div<::Rodin::Variational::GridFunction<FES, Data>>&
        node)
      : u(node.getOperand().getFiniteElementSpace(), node.getOperand().getData()),
        trace(node.getTraceDomain())
    {}

    /// @param p Device or integration point. @return Divergence @f$ \nabla \cdot u(p) @f$.
    template <class Point>
    KOKKOS_INLINE_FUNCTION Real getValue(const Point& p) const
    {
      return u.interpolateDivergence(p, trace);
    }

    /// @brief Device coefficients.
    GridFunctionView<FES> u;

    /// @brief Trace domain of the derivative on interfaces.
    TraceDomainView trace;
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

  /// @brief Scalar basis of a device space.
  template <class Derived, ScalarDeviceSpace FES,
    ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct DeviceView<::Rodin::Variational::ShapeFunction<Derived, FES, Space>>
  {
    /// @param node Shape function expression.
    DeviceView(const ::Rodin::Variational::ShapeFunction<Derived, FES, Space>& node)
      : dofs(node.getFiniteElementSpace().getDeviceDOFs())
    {}

    /**
     * @param ip Quadrature point. @param dof Local scalar degree of freedom.
     * @param component Vector component. @return Basis value.
     */
    KOKKOS_INLINE_FUNCTION
    Real getBasis(const IntegrationPointView& ip, Index dof, Index component) const
    {
      KOKKOS_ASSERT(component == 0);
      return dofs.getBasis(ip, dof);
    }

    /// @brief Device degrees of freedom.
    typename FES::DOFsType dofs;
  };

  /// @brief Physical gradient of the scalar basis of a device space.
  template <class Derived, ScalarDeviceSpace FES,
    ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct DeviceView<
    ::Rodin::Variational::Grad<::Rodin::Variational::ShapeFunction<Derived, FES, Space>>>
  {
    /// @param node Shape function expression.
    DeviceView(const ::Rodin::Variational::Grad<
      ::Rodin::Variational::ShapeFunction<Derived, FES, Space>>& node)
      : dofs(node.getOperand().getFiniteElementSpace().getDeviceDOFs())
    {}

    /**
     * @param ip Quadrature point. @param dof Local scalar degree of freedom.
     * @param component Vector component. @return Gradient.
     */
    KOKKOS_INLINE_FUNCTION
    decltype(auto) getBasis(
      const IntegrationPointView& ip, Index dof, Index component) const
    {
      KOKKOS_ASSERT(component == 0);
      return dofs.getGradient(ip, dof);
    }

    /// @brief Device degrees of freedom.
    typename FES::DOFsType dofs;
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

    /**
     * @param ip Quadrature point. @param dof Local scalar degree of freedom.
     * @param component Vector component. @return Scaled basis value.
     */
    KOKKOS_INLINE_FUNCTION
    auto getBasis(const IntegrationPointView& ip, Index dof, Index component) const
    {
      return mult(evaluate(lhs, ip), rhs.getBasis(ip, dof, component));
    }

    /// @param ip Integration point. @return View bound to the point.
    template <class Point>
    KOKKOS_INLINE_FUNCTION auto bind(const Point& ip) const
    {
      using ValueType = decltype(evaluate(lhs, ip));
      using BasisType = decltype(bindBasis(rhs, ip));
      return Internal::LeftScaledBasis<ValueType, BasisType>{
        evaluate(lhs, ip), bindBasis(rhs, ip)};
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

    /**
     * @param ip Quadrature point. @param dof Local scalar degree of freedom.
     * @param component Vector component. @return Scaled basis value.
     */
    KOKKOS_INLINE_FUNCTION
    auto getBasis(const IntegrationPointView& ip, Index dof, Index component) const
    {
      return mult(lhs.getBasis(ip, dof, component), evaluate(rhs, ip));
    }

    /// @param ip Integration point. @return View bound to the point.
    template <class Point>
    KOKKOS_INLINE_FUNCTION auto bind(const Point& ip) const
    {
      using ValueType = decltype(evaluate(rhs, ip));
      using BasisType = decltype(bindBasis(lhs, ip));
      return Internal::RightScaledBasis<BasisType, ValueType>{
        bindBasis(lhs, ip), evaluate(rhs, ip)};
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

    /**
     * @param ip Quadrature point. @param dof Local scalar degree of freedom.
     * @param component Vector component. @return Dot product value.
     */
    KOKKOS_INLINE_FUNCTION
    Real getBasis(const IntegrationPointView& ip, Index dof, Index component) const
    {
      return dot(evaluate(lhs, ip), rhs.getBasis(ip, dof, component));
    }

    /// @param ip Integration point. @return View bound to the point.
    template <class Point>
    KOKKOS_INLINE_FUNCTION auto bind(const Point& ip) const
    {
      using ValueType = decltype(evaluate(lhs, ip));
      using BasisType = decltype(bindBasis(rhs, ip));
      return Internal::DotBasis<ValueType, BasisType>{
        evaluate(lhs, ip), bindBasis(rhs, ip)};
    }

    /// @brief Function factor.
    DeviceView<::Rodin::Variational::FunctionBase<LHSDerived>> lhs;

    /// @brief Shape function factor.
    DeviceView<::Rodin::Variational::ShapeFunctionBase<RHSDerived, FES, Space>> rhs;
  };

  /**
   * @brief Vector basis of a device space.
   *
   * Local index @f$ a\,d + r @f$ is the basis function
   * @f$ \phi_a e_r @f$, as on the host.
   */
  template <class Derived, VectorDeviceSpace FES,
    ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct DeviceView<::Rodin::Variational::ShapeFunction<Derived, FES, Space>>
  {
    /// @param node Shape function expression.
    DeviceView(const ::Rodin::Variational::ShapeFunction<Derived, FES, Space>& node)
      : vdim(node.getFiniteElementSpace().getVectorDimension()),
        dofs(node.getFiniteElementSpace().getDeviceDOFs())
    {}

    /**
     * @param ip Quadrature point. @param dof Local scalar degree of freedom.
     * @param component Vector component. @return Basis value.
     */
    KOKKOS_INLINE_FUNCTION
    UnitVector getBasis(
      const IntegrationPointView& ip, Index dof, Index component) const
    {
      KOKKOS_ASSERT(component < vdim);
      return {vdim, component, dofs.getBasis(ip, dof)};
    }

    /// @brief Vector dimension.
    Index vdim;

    /// @brief Device degrees of freedom.
    typename FES::DOFsType dofs;
  };

  /**
   * @brief Physical Jacobian of the vector basis of a device space.
   *
   * Row @f$ r @f$ of the Jacobian of @f$ \phi_a e_r @f$ is
   * @f$ \nabla \phi_a @f$, the other rows vanish.
   */
  template <class Derived, VectorDeviceSpace FES,
    ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct DeviceView<::Rodin::Variational::Jacobian<
    ::Rodin::Variational::ShapeFunction<Derived, FES, Space>>>
  {
    /// @param node Shape function expression.
    DeviceView(const ::Rodin::Variational::Jacobian<
      ::Rodin::Variational::ShapeFunction<Derived, FES, Space>>& node)
      : vdim(node.getFiniteElementSpace().getVectorDimension()),
        dofs(node.getFiniteElementSpace().getDeviceDOFs())
    {}

    /**
     * @param ip Quadrature point. @param dof Local scalar degree of freedom.
     * @param component Vector component. @return Jacobian.
     */
    KOKKOS_INLINE_FUNCTION
    RowMatrix getBasis(const IntegrationPointView& ip, Index dof, Index component) const
    {
      KOKKOS_ASSERT(component < vdim);
      return {vdim, component, dofs.getGradient(ip, dof)};
    }

    /// @brief Vector dimension.
    Index vdim;

    /// @brief Device degrees of freedom.
    typename FES::DOFsType dofs;
  };

  /**
   * @brief Divergence of the vector basis of a device space.
   *
   * The divergence of @f$ \phi_a e_r @f$ is @f$ \partial_r \phi_a @f$.
   */
  template <class Derived, VectorDeviceSpace FES,
    ::Rodin::Variational::ShapeFunctionSpaceType Space>
  struct DeviceView<
    ::Rodin::Variational::Div<::Rodin::Variational::ShapeFunction<Derived, FES, Space>>>
  {
    /// @param node Shape function expression.
    DeviceView(const ::Rodin::Variational::Div<
      ::Rodin::Variational::ShapeFunction<Derived, FES, Space>>& node)
      : vdim(node.getFiniteElementSpace().getVectorDimension()),
        dofs(node.getFiniteElementSpace().getDeviceDOFs())
    {}

    /**
     * @param ip Quadrature point. @param dof Local scalar degree of freedom.
     * @param component Vector component. @return Divergence.
     */
    KOKKOS_INLINE_FUNCTION
    Real getBasis(const IntegrationPointView& ip, Index dof, Index component) const
    {
      KOKKOS_ASSERT(component < vdim);
      return dofs.getGradient(ip, dof)(component);
    }

    /// @brief Vector dimension.
    Index vdim;

    /// @brief Device degrees of freedom.
    typename FES::DOFsType dofs;
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

    /**
     * @param ip Quadrature point. @param dof Local scalar degree of freedom.
     * @param component Vector component. @return Basis value.
     */
    KOKKOS_INLINE_FUNCTION
    auto getBasis(const IntegrationPointView& ip, Index dof, Index component) const
    {
      return add(lhs.getBasis(ip, dof, component), rhs.getBasis(ip, dof, component));
    }

    /// @param ip Integration point. @return View bound to the point.
    template <class Point>
    KOKKOS_INLINE_FUNCTION auto bind(const Point& ip) const
    {
      using LHSType = decltype(bindBasis(lhs, ip));
      using RHSType = decltype(bindBasis(rhs, ip));
      return Internal::SumBasis<LHSType, RHSType>{
        bindBasis(lhs, ip), bindBasis(rhs, ip)};
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

    /**
     * @param ip Quadrature point. @param dof Local scalar degree of freedom.
     * @param component Vector component. @return Basis value.
     */
    KOKKOS_INLINE_FUNCTION
    auto getBasis(const IntegrationPointView& ip, Index dof, Index component) const
    {
      return mult(Real(-1), operand.getBasis(ip, dof, component));
    }

    /// @param ip Integration point. @return View bound to the point.
    template <class Point>
    KOKKOS_INLINE_FUNCTION auto bind(const Point& ip) const
    {
      using BasisType = decltype(bindBasis(operand, ip));
      return Internal::NegatedBasis<BasisType>{bindBasis(operand, ip)};
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

    /**
     * @param ip Quadrature point. @param dof Local scalar degree of freedom.
     * @param component Vector component. @return Basis value.
     */
    KOKKOS_INLINE_FUNCTION
    Geometry::Matrix getBasis(
      const IntegrationPointView& ip, Index dof, Index component) const
    {
      return transpose(operand.getBasis(ip, dof, component));
    }

    /// @param ip Integration point. @return View bound to the point.
    template <class Point>
    KOKKOS_INLINE_FUNCTION auto bind(const Point& ip) const
    {
      using BasisType = decltype(bindBasis(operand, ip));
      return Internal::TransposedBasis<BasisType>{bindBasis(operand, ip)};
    }

    /// @brief Operand.
    DeviceView<::Rodin::Variational::ShapeFunctionBase<NestedDerived, FES, Space>>
      operand;
  };
}

#endif
