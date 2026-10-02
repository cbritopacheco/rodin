/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef RODIN_TESTS_CONVERGENCE_CURVED_GEOMETRY_H
#define RODIN_TESTS_CONVERGENCE_CURVED_GEOMETRY_H

#include <functional>
#include <vector>

#include "Convergence.h"

namespace Rodin::Tests::Convergence
{
  /**
   * @brief Installs the quadratic map
   * @f$\Phi(\xi)=\xi+0.1\xi_0^2e_{d-1}@f$ on a unit-box grid.
   * @par Architecture
   * Original vertex positions are retained before any coordinate mutation.
   * Geometry control points are evaluated from the original P1 chart, then
   * mapped once by @f$\Phi@f$. Transformations are installed on every
   * positive-dimensional entity, including traces and MPI halo entities.
   * The same deterministic physical map is used on each shard; topology,
   * ownership and reconciliation metadata are not changed. Entity iterators
   * enumerate local entities (including halos); global MPI counts are never
   * used as local-index bounds. The mesh must
   * remain alive and its topology must not change while this object is used.
   * In dimensions two and three, @f$\det D\Phi=1@f$; in dimension one,
   * @f$\Phi'=1+0.2\xi_0>0@f$. P2 geometry represents the map exactly.
   */
  template <class MeshType>
  class CurvedGeometry
  {
    public:
      explicit CurvedGeometry(MeshType& mesh)
        : m_mesh(mesh)
      {
        for (auto vertex = mesh.getVertex(); vertex; ++vertex)
        {
          assert(vertex->getIndex() == m_vertices.size());
          m_vertices.push_back(mesh.getVertexCoordinates(vertex->getIndex()));
        }
      }

      Math::SpatialPoint mapToPhysical(Math::SpatialPoint point) const
      {
        assert(point.size() >= 1);
        point(point.size() - 1) += Real(0.1) * point(0) * point(0);
        return point;
      }

      Math::SpatialPoint referencePosition(
        const Geometry::Polytope& polytope, const Math::SpatialPoint& rc) const
      {
        Variational::RealP1Element affine(polytope.getGeometry());
        Math::SpatialPoint point(m_mesh.get().getSpaceDimension());
        point.setZero();
        const auto vertices = polytope.getVertices();
        for (size_t local = 0; local < affine.getCount(); ++local)
          point += m_vertices.at(vertices[local]) * affine.getBasis(local)(rc);
        return point;
      }

      template <size_t GeometryOrder>
      void install()
      {
        static_assert(GeometryOrder >= 1);
        auto& mesh = m_mesh.get();
        const size_t dim = mesh.getSpaceDimension();
        for (Index vertex = 0; vertex < m_vertices.size(); ++vertex)
          mesh.setVertexCoordinates(vertex, mapToPhysical(m_vertices.at(vertex)));
        if constexpr (GeometryOrder > 1)
        {
          for (size_t d = 1; d <= dim; ++d)
          {
            for (auto polytope = mesh.getPolytope(d); polytope; ++polytope)
            {
              Variational::RealH1Element<GeometryOrder> element(polytope->getGeometry());
              Geometry::PointCloud nodes(dim, element.getCount());
              for (size_t local = 0; local < element.getCount(); ++local)
              {
                const auto point =
                  mapToPhysical(referencePosition(*polytope, element.getNode(local)));
                for (size_t coordinate = 0; coordinate < dim; ++coordinate)
                  nodes(coordinate, local) = point(coordinate);
              }
              mesh.setPolytopeTransformation({d, polytope->getIndex()},
                new Geometry::ParametricTransformation<
                  Variational::RealH1Element<GeometryOrder>>(
                  std::move(nodes), std::move(element)));
            }
          }
        }
      }

    private:
      std::reference_wrapper<MeshType> m_mesh;
      std::vector<Math::SpatialPoint> m_vertices;
  };
}

#endif
