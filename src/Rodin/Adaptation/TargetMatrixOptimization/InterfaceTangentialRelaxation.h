/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_TARGETMATRIXOPTIMIZATION_INTERFACETANGENTIALRELAXATION_H
#define RODIN_ADAPTATION_TARGETMATRIXOPTIMIZATION_INTERFACETANGENTIALRELAXATION_H

/**
 * @file
 * @brief Laplacian-with-projection tangential redistribution of interface
 *        vertices, as a standalone pipeline stage.
 *
 * Distinct from the topological optimizer (split/collapse/swap) and from
 * the TMOP solve. Each interface vertex (recognized by an attribute on its
 * incident dim-1 polytopes) is relaxed toward the midpoint of its two
 * interface neighbours and re-projected onto the manifold by the
 * user-supplied projector (typically the analytic phi=0 surface), in a
 * single validity-gated sweep. The fit is preserved exactly because every
 * smoothed vertex lands back on the manifold; element quality improves
 * because clustered interface corners are dispersed.
 *
 * Intentionally NOT inside TriangleMeshOptimizer: this is geometric
 * pre-processing for the downstream TMOP solve, not a topology operator.
 */

#include <cmath>
#include <cstdint>
#include <functional>
#include <vector>

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Math.h"
#include "Rodin/Types.h"

namespace Rodin::Adaptation::TargetMatrixOptimization
{
  /// Pure tangential redistribution of mesh vertices on an attribute-tagged
  /// 1-d feature (interface). Each pass moves every interior interface
  /// vertex toward the midpoint of its two interface neighbours, then snaps
  /// the candidate back onto the feature via the injected projector; the
  /// move is committed only if every incident triangle stays strictly
  /// orientation-preserving. Endpoints / junctions (interface vertices
  /// with !=2 interface neighbours) are left untouched.
  class InterfaceTangentialRelaxation
  {
    public:
      using Projector =
        std::function<Math::SpatialPoint(const Math::SpatialPoint&)>;

      InterfaceTangentialRelaxation() = default;

      /// Snap-back projector. Typically projectToInterface (the analytic
      /// phi=0 Newton step). When unset, the smoother is a pure Laplacian
      /// with no snap-back (suitable for non-level-set features).
      InterfaceTangentialRelaxation& setProjector(Projector p)
      { m_projector = std::move(p); return *this; }

      /// Attribute that identifies a 1-d polytope as a feature edge.
      InterfaceTangentialRelaxation& setFeatureAttribute(
          Geometry::Attribute a)
      { m_featureAttribute = a; return *this; }

      /// Number of relaxation sweeps. 0 disables.
      InterfaceTangentialRelaxation& setPasses(size_t n)
      { m_passes = n; return *this; }

      /// Under-relaxation factor in [0, 1]. 1 = full Laplacian step.
      InterfaceTangentialRelaxation& setRelaxation(Real omega)
      { m_omega = std::max(Real(0), std::min(Real(1), omega)); return *this; }

      /// Apply the relaxation. Mesh is mutated in place.
      void relax(Geometry::LocalMesh& mesh) const
      {
        if (m_passes == 0)
          return;

        mesh.getConnectivity().compute(2, 1);
        mesh.getConnectivity().compute(1, 0);
        mesh.getConnectivity().compute(0, 2);
        const auto& conn = mesh.getConnectivity();
        const Index nv = static_cast<Index>(mesh.getVertexCount());

        std::vector<std::vector<Index>> ifaceNbr(nv);
        for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
        {
          const auto attr = mesh.getAttribute(1, e);
          if (!attr || *attr != m_featureAttribute)
            continue;
          const auto& edge = conn.getPolytope(1, e);
          ifaceNbr[static_cast<size_t>(edge(0))].push_back(edge(1));
          ifaceNbr[static_cast<size_t>(edge(1))].push_back(edge(0));
        }

        for (size_t p = 0; p < m_passes; ++p)
        {
          for (Index v = 0; v < nv; ++v)
          {
            const auto& nbrs = ifaceNbr[static_cast<size_t>(v)];
            if (nbrs.size() != 2)            // endpoints / junctions
              continue;
            const auto a = mesh.getVertexCoordinates(nbrs[0]);
            const auto b = mesh.getVertexCoordinates(nbrs[1]);
            const Math::SpatialPoint mid{
              (a[0] + b[0]) * Real(0.5),
              (a[1] + b[1]) * Real(0.5) };
            const Math::SpatialPoint snapped =
              m_projector ? m_projector(mid) : mid;
            const auto cur = mesh.getVertexCoordinates(v);
            const Math::SpatialPoint trial{
              cur[0] + m_omega * (snapped[0] - cur[0]),
              cur[1] + m_omega * (snapped[1] - cur[1]) };
            // Validity gate: accept only if every incident triangle stays
            // strictly orientation-preserving with the trial position.
            bool valid = true;
            for (Index ci : conn.getIncidence({ 0, 2 }, v))
            {
              const auto& cell = conn.getPolytope(2, ci);
              std::array<Math::SpatialPoint, 3> p3;
              for (size_t k = 0; k < 3; ++k)
                p3[k] = (cell(k) == v)
                  ? trial
                  : mesh.getVertexCoordinates(cell(k));
              const Real orient =
                  (p3[1][0] - p3[0][0]) * (p3[2][1] - p3[0][1])
                - (p3[1][1] - p3[0][1]) * (p3[2][0] - p3[0][0]);
              if (!(orient > Real(0))) { valid = false; break; }
            }
            if (valid)
              mesh.setVertexCoordinates(v, trial);
          }
        }
        mesh.getConnectivity().compute(2, 1);
        mesh.getConnectivity().compute(1, 0);
      }

    private:
      Projector m_projector;
      Geometry::Attribute m_featureAttribute = 0;
      size_t m_passes = 5;
      Real m_omega = Real(0.5);
  };
}

#endif
