/*
 *          Copyright Carlos BRITO PACHECO 2026.
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_VIEW_INTEGRATIONPOINT_H
#define RODIN_KOKKOS_VARIATIONAL_VIEW_INTEGRATIONPOINT_H

/**
 * @file
 * @brief Prepared geometry and basis values at an integration point.
 */

#include <type_traits>
#include "Rodin/Configure.h"
#include "Rodin/Kokkos/Geometry/View.h"

namespace Rodin::Kokkos::Variational
{
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
      KOKKOS_INLINE_FUNCTION explicit IntegrationPointView(
        const ::Rodin::Kokkos::View<::Rodin::Geometry::Point>& point)
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
      KOKKOS_INLINE_FUNCTION IntegrationPointView(
        const ::Rodin::Kokkos::View<::Rodin::Geometry::Polytope>& polytope,
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
      KOKKOS_INLINE_FUNCTION IntegrationPointView(
        const ::Rodin::Kokkos::View<::Rodin::Geometry::Polytope>& polytope,
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
      KOKKOS_INLINE_FUNCTION void setPoint(
        const ::Rodin::Kokkos::View<::Rodin::Geometry::Point>& point)
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
          point.evaluate(a, m_basis[a], reference);
      }

      /**
       * @brief Moves the integration point within its polytope.
       * @param reference Reference coordinates in the current polytope.
       * @param geometry Tag of a geometry known at run time only.
       */
      KOKKOS_INLINE_FUNCTION void setPoint(
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
        const ::Rodin::Kokkos::View<::Rodin::Geometry::Polytope> polytope =
          m_point.getPolytope();
        evaluate<G, !Geometry::StaticBasis<G>::Affine>(polytope, reference);
      }

      /// @brief Returns the evaluation point. @return Device point.
      KOKKOS_INLINE_FUNCTION const ::Rodin::Kokkos::View<::Rodin::Geometry::Point>&
      getPoint() const
      {
        return m_point;
      }

      /// @brief Returns the integration measure. @return Distortion of the point.
      KOKKOS_INLINE_FUNCTION Real getDistortion() const
      {
        return m_distortion;
      }

      /// @brief Returns the number of local basis functions. @return Vertex count.
      KOKKOS_INLINE_FUNCTION Index getCount() const
      {
        return m_count;
      }

      /**
       * @brief Returns a basis value.
       * @param local Local vertex index.
       * @return Basis value @f$ \phi_a @f$.
       */
      KOKKOS_INLINE_FUNCTION Real getBasis(Index local) const
      {
        KOKKOS_ASSERT(local < m_count);
        return m_basis[local];
      }

      /**
       * @brief Returns a physical basis gradient.
       * @param local Local vertex index.
       * @return Gradient @f$ \nabla \phi_a @f$.
       */
      KOKKOS_INLINE_FUNCTION const Geometry::Vector& getGradient(Index local) const
      {
        KOKKOS_ASSERT(local < m_count);
        return m_gradients[local];
      }

    private:
      /**
       * @brief Evaluates the basis, the gradients and the distortion at a point of the polytope.
       * @param point Point on a mesh entity.
       */
      KOKKOS_INLINE_FUNCTION void evaluate(
        const ::Rodin::Kokkos::View<::Rodin::Geometry::Point>& point)
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
          point.evaluate(a, m_basis[a], reference);
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
        const ::Rodin::Kokkos::View<::Rodin::Geometry::Polytope>& polytope,
        const Geometry::Vector& reference)
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
        m_point =
          ::Rodin::Kokkos::View<::Rodin::Geometry::Point>(polytope, reference, physical);
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

      ::Rodin::Kokkos::View<::Rodin::Geometry::Point> m_point;
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
   * @brief Integration point type of the kernels of a geometry tag.
   *
   * A static geometry fills the space, so its polytopes are cells.
   */
  template <class Shape>
  using KernelPointView = std::conditional_t<std::is_same_v<Shape, DynamicGeometry>,
    IntegrationPointView, CellIntegrationPointView>;
}

#endif
