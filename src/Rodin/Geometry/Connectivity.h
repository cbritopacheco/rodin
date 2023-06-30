/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_CONNECTIVITY_H
#define RODIN_GEOMETRY_CONNECTIVITY_H

#include <set>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <boost/bimap.hpp>
#include <boost/bimap/unordered_set_of.hpp>

#include "Rodin/Array.h"

#include "ForwardDecls.h"

#include "Types.h"
#include "Simplex.h"

namespace Rodin::Geometry
{
  /**
   * @brief Represents the set of incidence relations.
   *
   * Stores the set of incidence relations:
   * @f[
   *  d \rightarrow d'
   * @f]
   * for a fixed pair of topological dimensions @f$ (d, d') @f$.
   */
  class MeshConnectivity
  {
    public:
      using PolytopeIndex =
        boost::bimap<
          boost::bimaps::unordered_set_of<IndexArray, IndexArraySymmetricHash, IndexArraySymmetricEquality>,
          boost::bimaps::unordered_set_of<Index>
        >;

      struct SubPolytope
      {
        Polytope::Geometry geometry;
        Array<Index> vertices;
      };

      MeshConnectivity() = default;

      MeshConnectivity(const MeshConnectivity&) = default;

      MeshConnectivity(MeshConnectivity&&) = default;

      MeshConnectivity& operator=(const MeshConnectivity&) = default;

      MeshConnectivity& operator=(MeshConnectivity&&) = default;

      MeshConnectivity& initialize(size_t maximalDimension);

      MeshConnectivity& nodes(size_t count);

      MeshConnectivity& reserve(size_t d, size_t count);

      MeshConnectivity& polytope(
          Geometry::Polytope::Geometry t, std::initializer_list<Index> p)
      {
        Array<Index> arr(p.size());
        std::copy(p.begin(), p.end(), arr.begin());
        return polytope(t, std::move(arr));
      }

      MeshConnectivity& polytope(
          Geometry::Polytope::Geometry t, const Array<Index>& polytope);

      MeshConnectivity& polytope(
          Geometry::Polytope::Geometry t, Array<Index>&& polytope);

      MeshConnectivity& compute(size_t d, size_t dp);

      size_t getCount(size_t dim) const;

      size_t getCount(Polytope::Geometry g) const;

      size_t getMeshDimension() const;

      const PolytopeIndex& getIndexMap(size_t dim) const;

      const std::optional<Index> getIndex(size_t dim, const IndexArray& key) const;

      Polytope::Geometry getGeometry(size_t d, Index idx) const;

      const Array<Index>& getPolytope(size_t d, Index idx) const;

      const Incidence& getIncidence(size_t d, size_t dp) const;

      const IndexSet& getIncidence(const std::pair<size_t, size_t> p, Index idx) const;

    protected:
      void local(std::vector<SubPolytope>& out, size_t dim, Index i);

      /**
       * D -> d and D -> 0 from D -> 0 and D -> D
       */
      MeshConnectivity& build(size_t d);

      MeshConnectivity& transpose(size_t d, size_t dp);

      MeshConnectivity& intersection(size_t d, size_t dp, size_t dpp);

    private:
      size_t m_maximalDimension;

      std::vector<std::vector<bool>> m_dirty;

      // std::vector<std::vector<Array<Index>>> m_polytopes;
      std::vector<std::vector<Polytope::Geometry>> m_geometry;

      std::vector<size_t> m_count;
      FlatMap<Polytope::Geometry, size_t> m_gcount;
      std::vector<PolytopeIndex> m_index;
      std::vector<std::vector<Incidence>> m_connectivity;
  };
}

#endif
