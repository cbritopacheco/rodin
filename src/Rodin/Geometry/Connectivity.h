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
#include "Rodin/Context/Sequential.h"

#include "ForwardDecls.h"

#include "Types.h"
#include "Polytope.h"
#include "GeometryIndexed.h"

namespace Rodin::Geometry
{
  class ConnectivityBase
  {
    public:
      virtual Polytope::Type getGeometry(size_t d, Index idx) const = 0;
      virtual const Array<Index>& getPolytope(size_t d, Index idx) const = 0;
      virtual size_t getCount(size_t dim) const = 0;
      virtual size_t getCount(Polytope::Type g) const = 0;
      virtual size_t getMeshDimension() const = 0;
      virtual const Incidence& getIncidence(size_t d, size_t dp) const = 0;
      virtual const IndexSet& getIncidence(const std::pair<size_t, size_t> p, Index idx) const = 0;
  };

  using SequentialConnectivity = Connectivity<Context::Local>;

  /**
   * @brief Represents the set of incidence relations of a Mesh.
   *
   * This class stores the set of incidence relations:
   * @f[
   *  d \longrightarrow d', \quad 0 \leq d, d' \leq D,
   * @f]
   * where @f$ D @f$ represents the topological mesh dimension. This class is
   * based on @cite logg2009efficient.
   *
   */
  template <>
  class Connectivity<Context::Local> final : public ConnectivityBase
  {
    public:
      using PolytopeIndex =
        boost::bimap<
          boost::bimaps::unordered_set_of<IndexArray, IndexArraySymmetricHash, IndexArraySymmetricEquality>,
          boost::bimaps::unordered_set_of<Index>
        >;

      struct SubPolytope
      {
        Polytope::Type geometry;
        Array<Index> vertices;
      };

      Connectivity();

      Connectivity(const Connectivity&) = default;

      Connectivity(Connectivity&&) = default;

      Connectivity& operator=(const Connectivity&) = default;

      Connectivity& operator=(Connectivity&&) = default;

      Connectivity& initialize(size_t maximalDimension);

      /**
       * @brief Sets the number of nodes (vertices) in the mesh.
       */
      Connectivity& nodes(size_t count);

      Connectivity& clear(size_t d, size_t dp);

      /**
       * @brief Reserves space for the polytopes of the given dimension.
       */
      Connectivity& reserve(size_t d, size_t count);

      Connectivity& polytope(
          Geometry::Polytope::Type t, std::initializer_list<Index> p)
      {
        Array<Index> arr(p.size());
        std::copy(p.begin(), p.end(), arr.begin());
        return polytope(t, std::move(arr));
      }

      Connectivity& polytope(
          Geometry::Polytope::Type t, const Array<Index>& polytope);

      Connectivity& polytope(
          Geometry::Polytope::Type t, Array<Index>&& polytope);

      /**
       * @brief Computes the entities of dimension @f$ d @f$ of each cell and
       * for each such entity the vertices of that entity.
       *
       * Computes the connectivities:
       * @f[
       *  D \longrightarrow d \quad \text{and} \quad D \longrightarrow 0, \quad 0 < d < D,
       * @f]
       * from @f$ D \longrightarrow 0 @f$ and @f$ D \longrightarrow D @f$.
       */
      Connectivity& build(size_t d);

      Connectivity& local(size_t i, size_t d);

      Connectivity& compute(size_t d, size_t dp);

      Connectivity& transpose(size_t d, size_t dp);

      Connectivity& intersection(size_t d, size_t dp, size_t dpp);

      void getSubPolytopes(std::vector<SubPolytope>& out, Index i, size_t d);

      const PolytopeIndex& getIndexMap(size_t dim) const;

      const std::optional<Index> getIndex(size_t dim, const IndexArray& key) const;

      Connectivity& setIncidence(const std::pair<size_t, size_t>& p, Incidence&& inc);

      size_t getCount(size_t dim) const override;

      size_t getCount(Polytope::Type g) const override;

      size_t getMeshDimension() const override;

      Polytope::Type getGeometry(size_t d, Index idx) const override;

      const Array<Index>& getPolytope(size_t d, Index idx) const override;

      const Incidence& getIncidence(size_t d, size_t dp) const override;

      const IndexSet& getIncidence(const std::pair<size_t, size_t> p, Index idx) const override;

    private:
      size_t m_maximalDimension;
      std::vector<size_t> m_count;
      GeometryIndexed<size_t> m_gcount;
      std::vector<PolytopeIndex> m_index;
      std::vector<std::vector<bool>> m_dirty;
      std::vector<std::vector<Polytope::Type>> m_geometry;
      std::vector<std::vector<Incidence>> m_connectivity;

  };
}

#endif
