/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_DOFS_H
#define RODIN_KOKKOS_VARIATIONAL_DOFS_H

/**
 * @file
 * @brief Degrees of freedom of the device spaces.
 *
 * A device space exposes its degrees of freedom by `getDeviceDOFs()`, an
 * object copied by value into kernels that provides:
 *
 * | Member | Description |
 * |--------|-------------|
 * | `static constexpr Index getCapacity(Index vertices)` | Maximum local count on a polytope of that many vertices. |
 * | `Index getStride() const` | Offset between the components of a vector space. |
 * | `bind(const Geometry::PolytopeView&) const` | Scalar global indices of a polytope, in local order. |
 * | `getBasis(const IntegrationPointView&, Index) const` | Basis value at an integration point. |
 * | `getGradient(const IntegrationPointView&, Index) const` | Physical basis gradient at an integration point. |
 * | `Geometry::Vector getNode(const Geometry::PolytopeView&, Index) const` | Reference coordinates of a Lagrange node, when the space interpolates at nodes. |
 *
 * Component @f$ r @f$ of the scalar degree of freedom @f$ q @f$ of a vector
 * space has the global index @f$ q + r s @f$, where @f$ s @f$ is the stride,
 * and the local index @f$ a\,d + r @f$ when @f$ q @f$ is the local index
 * @f$ a @f$ and @f$ d @f$ the vector dimension.
 */

#include "Rodin/Configure.h"
#include "Rodin/Geometry/Polytope.h"

#include "Rodin/Kokkos/Geometry/MeshView.h"

namespace Rodin::Kokkos::Variational
{
  /**
   * @brief Tells whether a space is resident on a Kokkos mesh.
   *
   * Device spaces provide the degrees of freedom read by the kernels.
   */
  template <class FES>
  concept DeviceSpace = requires(const FES& fes) { fes.getDeviceDOFs(); };

  /**
   * @brief Reference coordinates of a vertex of a polytope.
   * @param polytope Device polytope.
   * @param a Local vertex index.
   * @return Reference coordinates of the P1 node of local vertex a.
   */
  KOKKOS_INLINE_FUNCTION
  Geometry::Vector getReferenceVertex(const Geometry::PolytopeView& polytope, Index a)
  {
    const Index rdim = polytope.getDimension();
    Geometry::Vector result(rdim);
    for (Index d = 0; d < rdim; ++d)
      result(d) = Geometry::getReferenceVertexCoordinate(polytope.getGeometry(), a, d);
    return result;
  }

  /**
   * @brief Degrees of freedom of the P1 spaces, one per vertex.
   *
   * The scalar degree of freedom of a vertex is its index, and the basis is
   * the P1 basis cached by the integration point.
   */
  struct VertexDOFs
  {
    /// @param vertices Vertex count. @return Vertex count.
    KOKKOS_INLINE_FUNCTION
    static constexpr Index getCapacity(Index vertices)
    {
      return vertices;
    }

    /// @brief Gets the offset between components. @return Vertex count of the mesh.
    KOKKOS_INLINE_FUNCTION
    Index getStride() const
    {
      return stride;
    }

    /// @param polytope Device polytope. @return Its vertices.
    KOKKOS_INLINE_FUNCTION
    Geometry::IndexView bind(const Geometry::PolytopeView& polytope) const
    {
      return polytope.getVertices();
    }

    /// @param ip Integration point. @param q Local vertex index. @return Basis value.
    template <class Point>
    KOKKOS_INLINE_FUNCTION Real getBasis(const Point& ip, Index q) const
    {
      return ip.getBasis(q);
    }

    /// @param ip Integration point. @param q Local vertex index. @return Physical gradient.
    template <class Point>
    KOKKOS_INLINE_FUNCTION const Geometry::Vector& getGradient(
      const Point& ip, Index q) const
    {
      return ip.getGradient(q);
    }

    /// @param polytope Device polytope. @param q Local vertex index. @return Reference vertex.
    KOKKOS_INLINE_FUNCTION
    Geometry::Vector getNode(const Geometry::PolytopeView& polytope, Index q) const
    {
      return getReferenceVertex(polytope, q);
    }

    /// @brief Number of vertices of the mesh.
    Index stride;
  };

  /// @brief Local to global map of a single global degree of freedom.
  struct ConstantIndex
  {
    /// @return One.
    KOKKOS_INLINE_FUNCTION
    Index size() const
    {
      return 1;
    }

    /// @param q Local degree of freedom, necessarily zero. @return Zero, the global index.
    KOKKOS_INLINE_FUNCTION
    Index operator[]([[maybe_unused]] Index q) const
    {
      return 0;
    }
  };

  /**
   * @brief Degrees of freedom of the global constant space.
   *
   * The single basis function is one on the whole mesh. It is also the
   * implicit trial space of a linear form.
   */
  struct ConstantDOFs
  {
    /// @param vertices Vertex count, ignored. @return One.
    KOKKOS_INLINE_FUNCTION
    static constexpr Index getCapacity([[maybe_unused]] Index vertices)
    {
      return 1;
    }

    /// @return One.
    KOKKOS_INLINE_FUNCTION
    Index getStride() const
    {
      return 1;
    }

    /// @param polytope Device polytope, ignored. @return The single degree of freedom.
    KOKKOS_INLINE_FUNCTION
    ConstantIndex bind([[maybe_unused]] const Geometry::PolytopeView& polytope) const
    {
      return {};
    }

    /// @param ip Integration point, ignored. @param q Local degree of freedom, necessarily zero. @return One.
    template <class Point>
    KOKKOS_INLINE_FUNCTION Real getBasis(
      [[maybe_unused]] const Point& ip, [[maybe_unused]] Index q) const
    {
      return 1;
    }

    /// @param ip Integration point. @param q Local degree of freedom, necessarily zero. @return Zero gradient.
    template <class Point>
    KOKKOS_INLINE_FUNCTION Geometry::Vector getGradient(
      const Point& ip, [[maybe_unused]] Index q) const
    {
      return Geometry::Vector(ip.getGradient(0).size());
    }
  };
}

#endif
