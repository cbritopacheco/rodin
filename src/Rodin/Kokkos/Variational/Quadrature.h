/*
 *          Copyright Carlos BRITO PACHECO 2026.
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_QUADRATURE_H
#define RODIN_KOKKOS_VARIATIONAL_QUADRATURE_H

/**
 * @file
 * @brief Prepared reference quadrature tables for finite-element integration.
 */

#include <algorithm>
#include <array>
#include <vector>
#include "Rodin/QF/QuadratureFormula.h"
#include "Rodin/Kokkos/Utility.h"
#include "View/IntegrationPoint.h"

namespace Rodin::Kokkos::Variational
{
  /**
   * @brief Quadrature tables of the polytopes of one dimension of a resident
   * mesh, one formula per reference geometry.
   *
   * ## Architecture
   *
   * The quadrature order is selected on the host, once per reference
   * geometry, using the first polytope of each geometry. It must therefore
   * depend only on the geometry, which holds for meshes whose polytopes of a
   * given geometry share one transformation degree. The reference
   * coordinates are flattened and padded to three components.
   *
   * @tparam Mesh Resident Kokkos mesh type.
   */
  template <class Mesh>
  class Quadrature
  {
    public:
      /// @brief Kokkos execution space.
      using ExecutionSpace = typename Mesh::ExecutionSpace;

      /// @brief Kokkos memory space.
      using MemorySpace = typename Mesh::MemorySpace;

      /// @brief Number of reference geometries.
      static constexpr size_t GeometryCount = ::Rodin::Geometry::Polytope::Types.size();

      /// @brief Stride of the flattened reference coordinates.
      static constexpr Index ReferenceStride = 3;

      /// @brief Prepares an empty rule for a region with no candidate entities.
      Quadrature() = default;

      /**
       * @brief Builds and uploads the tables of every geometry present.
       * @param mesh Resident mesh.
       * @param dimension Dimension of the integrated polytopes.
       * @param order Callable giving the quadrature order on a host polytope.
       */
      template <class Order>
      Quadrature(const Mesh& mesh, size_t dimension, const Order& order)
      {
        const auto& execution = mesh.getContext().getExecutionSpace();
        m_representatives = getRepresentatives(mesh, dimension);
        const auto& representatives = m_representatives;
        std::vector<Index> offsets(GeometryCount + 1, 0);
        std::vector<Real> weights;
        std::vector<Real> points;
        for (size_t g = 0; g < GeometryCount; ++g)
        {
          offsets[g] = weights.size();
          if (!representatives[g])
            continue;
          const auto it = mesh.getPolytope(dimension, *representatives[g]);
          const ::Rodin::Geometry::Polytope& polytope = *it;
          m_vertices = std::max<Index>(m_vertices, polytope.getVertices().size());
          m_minimum = m_minimum
            ? std::min<Index>(m_minimum, polytope.getVertices().size())
            : polytope.getVertices().size();
          m_orders[g] = order(polytope);
          const auto& qf =
            QF::QuadratureFormula::get(m_orders[g], polytope.getGeometry());
          for (size_t qp = 0; qp < qf.getSize(); ++qp)
          {
            weights.push_back(qf.getWeight(qp));
            const auto& point = qf.getPoint(qp);
            for (size_t k = 0; k < ReferenceStride; ++k)
            {
              points.push_back(
                k < static_cast<size_t>(point.size()) ? point(k) : Real(0));
            }
          }
        }
        offsets[GeometryCount] = weights.size();
        m_offsets =
          Rodin::Kokkos::Internal::upload(execution, offsets, "Rodin quadrature offsets");
        m_weights =
          Rodin::Kokkos::Internal::upload(execution, weights, "Rodin quadrature weights");
        m_points =
          Rodin::Kokkos::Internal::upload(execution, points, "Rodin quadrature points");
      }

      /**
       * @brief Tests whether the tables were built with the given orders.
       * @param mesh Resident mesh the tables were built on.
       * @param dimension Dimension of the integrated polytopes.
       * @param order Callable giving the quadrature order on a host polytope.
       * @return Whether every geometry present has the same order.
       */
      template <class Order>
      bool hasOrders(const Mesh& mesh, size_t dimension, const Order& order) const
      {
        for (size_t g = 0; g < GeometryCount; ++g)
        {
          if (!m_representatives[g])
            continue;
          const auto it = mesh.getPolytope(dimension, *m_representatives[g]);
          if (order(*it) != m_orders[g])
            return false;
        }
        return true;
      }

      /// @brief Gets the largest vertex count of the polytopes. @return Vertex count.
      Index getMaximumVertexCount() const
      {
        return m_vertices;
      }

      /**
       * @brief Returns the smallest vertex count of the geometries present.
       * @return Minimum vertex count, zero on an empty mesh.
       */
      Index getMinimumVertexCount() const
      {
        return m_minimum;
      }

      /// @param g Reference geometry. @return First quadrature point of the geometry.
      KOKKOS_INLINE_FUNCTION Index begin(Index g) const
      {
        return m_offsets(g);
      }

      /// @param g Reference geometry. @return Past the last quadrature point of the geometry.
      KOKKOS_INLINE_FUNCTION Index end(Index g) const
      {
        return m_offsets(g + 1);
      }

      /// @param q Quadrature point. @return Reference weight.
      KOKKOS_INLINE_FUNCTION Real getWeight(Index q) const
      {
        return m_weights(q);
      }

      /**
       * @brief Reads the reference coordinates of a quadrature point.
       * @param cell Device polytope.
       * @param q Quadrature point of the polytope geometry.
       * @return Reference coordinates, of the polytope dimension.
       */
      KOKKOS_INLINE_FUNCTION Rodin::Kokkos::Geometry::Vector getReference(
        const ::Rodin::Kokkos::View<::Rodin::Geometry::Polytope>& cell, Index q) const
      {
        const Index rdim = cell.getDimension();
        Rodin::Kokkos::Geometry::Vector reference(rdim);
        for (Index k = 0; k < rdim; ++k)
          reference(k) = m_points(q * ReferenceStride + k);
        return reference;
      }

      /**
       * @brief Maps a quadrature point onto a polytope.
       * @param cell Device polytope.
       * @param q Quadrature point of the polytope geometry.
       * @return Point on the polytope.
       */
      KOKKOS_INLINE_FUNCTION ::Rodin::Kokkos::View<::Rodin::Geometry::Point> getPoint(
        const ::Rodin::Kokkos::View<::Rodin::Geometry::Polytope>& cell, Index q) const
      {
        return cell.getPoint(getReference(cell, q));
      }

      /**
       * @brief Calls a functor with the static geometry tag of the polytopes.
       *
       * The tag applies when every polytope has one geometry with a static
       * basis, of dimension the spatial dimension.
       *
       * @param spaceDimension Spatial dimension of the mesh.
       * @param f Callable taking a StaticGeometry tag and returning whether it
       * integrated the polytopes.
       * @return Result of @p f, false when no tag applies.
       */
      template <class F>
      bool visitStaticGeometry(size_t spaceDimension, F&& f) const
      {
        using Type = ::Rodin::Geometry::Polytope::Type;
        Optional<Type> geometry;
        for (size_t g = 0; g < GeometryCount; ++g)
        {
          if (!m_representatives[g])
            continue;
          if (geometry)
            return false;
          geometry = static_cast<Type>(g);
        }
        if (!geometry ||
          ::Rodin::Geometry::Polytope::Traits(*geometry).getDimension() != spaceDimension)
          return false;
        switch (*geometry)
        {
          case Type::Triangle:
            return f(StaticGeometry<Type::Triangle>{});
          case Type::Quadrilateral:
            return f(StaticGeometry<Type::Quadrilateral>{});
          case Type::Tetrahedron:
            return f(StaticGeometry<Type::Tetrahedron>{});
          case Type::Hexahedron:
            return f(StaticGeometry<Type::Hexahedron>{});
          default:
            return false;
        }
      }

      /**
       * @brief Builds a range policy on the mesh execution space.
       * @param execution Execution space instance.
       * @param count Iteration count.
       * @return Range policy over @f$ [0, count) @f$.
       */
      static auto policy(const ExecutionSpace& execution, Index count)
      {
        // The kernels capture the mesh view, the quadrature tables and the
        // integrand, more than the 512 bytes CUDA passes as arguments by
        // default; the hint keeps them in the kernel arguments instead of
        // the constant memory copy, which waits for an event at each launch.
        return ::Kokkos::Experimental::require(
          ::Kokkos::RangePolicy<ExecutionSpace, ::Kokkos::IndexType<Index>>(
            execution, 0, count),
          ::Kokkos::Experimental::WorkItemProperty::HintLightWeight);
      }

      /**
       * @brief Finds the first polytope of every reference geometry.
       * @param mesh Resident mesh.
       * @param dimension Dimension of the polytopes.
       * @return Index of the first polytope per geometry, empty when absent.
       */
      static std::array<Optional<Index>, GeometryCount> getRepresentatives(
        const Mesh& mesh, size_t dimension)
      {
        const auto& execution = mesh.getContext().getExecutionSpace();
        constexpr Index None = std::numeric_limits<Index>::max();
        ::Kokkos::View<Index[GeometryCount], MemorySpace> first(
          ::Kokkos::view_alloc(execution, std::string("Rodin geometry representatives")));
        ::Kokkos::deep_copy(execution, first, None);
        const auto view = mesh.getView();
        findRepresentatives(
          execution, view, first, dimension, mesh.getPolytopeCount(dimension));
        auto host = ::Kokkos::create_mirror_view(first);
        ::Kokkos::deep_copy(execution, host, first);
        execution.fence("Rodin geometry representatives");
        std::array<Optional<Index>, GeometryCount> result;
        for (size_t g = 0; g < GeometryCount; ++g)
        {
          if (host(g) != None)
            result[g] = host(g);
        }
        return result;
      }

      /**
       * @brief Records the minimum polytope index of every geometry.
       * @param execution Execution space instance.
       * @param view Device mesh.
       * @param first Per-geometry minimum, initialized to the maximum index.
       * @param dimension Dimension of the polytopes.
       * @param count Polytope count.
       */
      static void findRepresentatives(const ExecutionSpace& execution,
        const ::Rodin::Kokkos::View<::Rodin::Geometry::MeshBase>& view,
        const ::Kokkos::View<Index[GeometryCount], MemorySpace>& first, Index dimension,
        Index count)
      {
        ::Kokkos::parallel_for("Rodin geometry representatives", policy(execution, count),
          Representatives{view, first, dimension});
      }

      /// @brief Records the minimum polytope index of every geometry.
      struct Representatives
      {
        /// @brief Device mesh.
        ::Rodin::Kokkos::View<::Rodin::Geometry::MeshBase> view;

        /// @brief Minimum polytope index of every geometry.
        ::Kokkos::View<Index[GeometryCount], MemorySpace> first;

        /// @brief Dimension of the polytopes.
        Index dimension;

        /**
         * @brief Lowers the minimum of the geometry of a polytope to its index.
         * @param p Polytope index.
         */
        KOKKOS_INLINE_FUNCTION void operator()(Index p) const
        {
          const Index g =
            static_cast<Index>(view.getPolytope(dimension, p).getGeometry());
          ::Kokkos::atomic_min(&first(g), p);
        }
      };

    private:
      ::Kokkos::View<Index*, MemorySpace> m_offsets;
      ::Kokkos::View<Real*, MemorySpace> m_weights;
      ::Kokkos::View<Real*, MemorySpace> m_points;
      std::array<Optional<Index>, GeometryCount> m_representatives;
      std::array<size_t, GeometryCount> m_orders{};
      Index m_vertices = 0;
      Index m_minimum = 0;
  };
}

#endif
