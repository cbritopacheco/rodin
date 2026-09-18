/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef KELVIN_BALL_ROTATED_BOUNDARY_POLICY_H
#define KELVIN_BALL_ROTATED_BOUNDARY_POLICY_H

#include <array>
#include <functional>
#include <stdexcept>

#include <Rodin/Advection/Lagrangian.h>

#include "Common.h"

namespace KelvinBall
{
  using namespace Rodin;
  using namespace Rodin::Geometry;

  /** Transports a characteristic across a rotational chamber cut. */
  template <class Mesh, class Locator>
  class RotatedBoundaryPolicy
  {
    public:
      RotatedBoundaryPolicy(Real dt, const Mesh& mesh, const Locator& locator,
        const std::array<RotationPair, 2>& pairs)
        : m_mesh(mesh), m_locator(locator), m_pairs(pairs), m_stop(dt, mesh)
      {}

      bool operator()(const BoundaryHit& hit) const
      {
        const auto& mesh = m_mesh.get();
        const size_t dimension = mesh.getDimension();
        const auto& faces =
          mesh.getConnectivity().getIncidence({dimension, dimension - 1}, hit.cell);
        if (hit.face >= faces.size())
          return m_stop(hit);
        const auto face = mesh.getPolytope(dimension - 1, faces[hit.face]);
        if (!face || !face->getAttribute())
          return m_stop(hit);

        Attribute target = 0;
        Math::SpatialMatrix<Real> transform(3, 3);
        bool periodic = false;
        for (const auto& pair : m_pairs.get())
        {
          if (*face->getAttribute() == pair.slave)
          {
            target = pair.master;
            transform = pair.rotation;
            periodic = true;
            break;
          }
          if (*face->getAttribute() == pair.master)
          {
            target = pair.slave;
            transform = pair.rotation.transpose();
            periodic = true;
            break;
          }
        }
        if (!periodic)
          return m_stop(hit);

        Math::SpatialPoint physical;
        mesh.getPolytopeTransformation(dimension, hit.cell).transform(physical, hit.rref);
        const auto mapped = m_locator.get().locate(target, transform * physical);
        if (!mapped)
          throw std::runtime_error(
            "A periodic characteristic could not cross a chamber cut.");
        const auto& incidence = mesh.getConnectivity().getIncidence(
          {dimension - 1, dimension}, mapped->getPolytope().getIndex());
        if (incidence.size() != 1)
          throw std::runtime_error("A periodic target face has no unique incident cell.");
        const Index cell = incidence[0];
        Math::SpatialPoint reference;
        mesh.getPolytopeTransformation(dimension, cell)
          .inverse(reference, mapped->getPhysicalCoordinates());
        const auto geometry = mesh.getGeometry(dimension, cell);
        const auto centroid = Polytope::Traits(geometry).getCentroid();
        reference = (Real(1) - Real(1e-10)) * reference + Real(1e-10) * centroid;
        hit.cell = cell;
        hit.rref = reference;
        return true;
      }

    private:
      std::reference_wrapper<const Mesh> m_mesh;
      std::reference_wrapper<const Locator> m_locator;
      std::reference_wrapper<const std::array<RotationPair, 2>> m_pairs;
      Advection::StopInsideBoundaryPolicy m_stop;
  };
}

#endif
