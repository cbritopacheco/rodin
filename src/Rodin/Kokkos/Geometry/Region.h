/*
 *          Copyright Carlos BRITO PACHECO 2026.
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_GEOMETRY_REGION_H
#define RODIN_KOKKOS_GEOMETRY_REGION_H

/**
 * @file
 * @brief Ownership, region and attribute selection on backend meshes.
 */

#include <vector>
#include "Rodin/Geometry/Region.h"
#include "Rodin/Kokkos/Geometry/View.h"
#include "Rodin/Kokkos/Utility.h"

namespace Rodin::Kokkos::Geometry
{
  /**
   * @brief Tests whether a cell belongs to an attribute selection.
   * @param cell Device cell.
   * @param attributes Selected attributes; empty selects every cell.
   * @return Whether the cell is integrated.
   */
  template <class AttributeView>
  KOKKOS_INLINE_FUNCTION bool isSelected(
    const ::Rodin::Kokkos::View<::Rodin::Geometry::Polytope>& cell,
    const AttributeView& attributes)
  {
    if (attributes.extent(0) == 0)
      return true;
    if (!cell.hasAttribute())
      return false;
    const auto attribute = cell.getAttribute();
    for (Index i = 0; i < attributes.extent(0); ++i)
    {
      if (attributes(i) == attribute)
        return true;
    }
    return false;
  }

  /**
   * @brief Tests whether a polytope is integrated over a region.
   * @param view Device mesh.
   * @param region Integration region.
   * @param polytope Device polytope of the dimension of the region.
   * @param attributes Selected attributes; empty selects every polytope of
   * the region.
   * @return Whether the polytope is integrated.
   *
   * Boundary faces are incident to one cell and interface faces to more
   * than one, as on the host. These regions need the face to cell
   * incidence. Polytopes of a distributed mesh are integrated by their
   * owner only.
   */
  template <class AttributeView>
  KOKKOS_INLINE_FUNCTION bool isIntegrated(
    const ::Rodin::Kokkos::View<::Rodin::Geometry::MeshBase>& view,
    ::Rodin::Geometry::Region region,
    const ::Rodin::Kokkos::View<::Rodin::Geometry::Polytope>& polytope,
    const AttributeView& attributes)
  {
    if (!view.isOwned(polytope.getDimension(), polytope.getIndex()))
      return false;
    if (region == ::Rodin::Geometry::Region::Boundary ||
      region == ::Rodin::Geometry::Region::Interface)
    {
      const Index d = view.getDimension();
      const Index n =
        view.getConnectivity().getIncidence(d - 1, d, polytope.getIndex()).size();
      if (region == ::Rodin::Geometry::Region::Boundary ? n != 1 : n <= 1)
        return false;
    }
    return isSelected(polytope, attributes);
  }

  /**
   * @brief Uploads the attribute selection of an integrator.
   * @param execution Execution space instance.
   * @param attributes Selected attributes; empty selects every cell.
   * @return Device copy of the selection.
   */
  template <class Execution>
  auto uploadAttributes(
    const Execution& execution, const FlatSet<::Rodin::Geometry::Attribute>& attributes)
  {
    const std::vector<::Rodin::Geometry::Attribute> list(
      attributes.begin(), attributes.end());
    return Rodin::Kokkos::Internal::upload(
      execution, list, "Rodin integration attributes");
  }

  /**
   * @brief Dimension of the polytopes of an integration region.
   * @param dimension Mesh dimension.
   * @param region Integration region.
   * @return The mesh dimension for cells, the face dimension otherwise.
   * A zero-dimensional mesh has no face candidates.
   */
  inline size_t getRegionDimension(size_t dimension, ::Rodin::Geometry::Region region)
  {
    return region == ::Rodin::Geometry::Region::Cells || dimension == 0 ? dimension
                                                                        : dimension - 1;
  }

  /**
   * @brief Gets the topological dimension defining the integration regions.
   * @param mesh Local mesh or distributed shard.
   * @return The dimension of the whole mesh, independent of the rank.
   */
  template <class Mesh>
  size_t getDomainDimension(const Mesh& mesh)
  {
    if constexpr (requires { mesh.getGlobalDimension(); })
      return mesh.getGlobalDimension();
    else
      return mesh.getDimension();
  }

  /**
   * @brief Counts candidates of an integration region on this rank.
   * @param mesh Local mesh or distributed shard.
   * @param region Integration region.
   * @return Candidate count before ownership and attribute filtering.
   */
  template <class Mesh>
  size_t getRegionPolytopeCount(const Mesh& mesh, ::Rodin::Geometry::Region region)
  {
    if (region == ::Rodin::Geometry::Region::Cells)
      return mesh.getCellCount();
    const size_t dimension = getDomainDimension(mesh);
    if (dimension == 0 || mesh.getDimension() != dimension)
      return 0;
    return mesh.getPolytopeCount(getRegionDimension(dimension, region));
  }

}

#endif
