/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_CONNECTIVITY_H
#define RODIN_GEOMETRY_CONNECTIVITY_H

/**
 * @file
 * @brief Mesh connectivity for representing topological relationships.
 */

#include <vector>
#include <boost/bimap.hpp>
#include <boost/bimap/vector_of.hpp>
#include <boost/bimap/unordered_set_of.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/split_member.hpp>
#include "Rodin/Serialization/UnorderedMap.h"

#include "Rodin/Array.h"
#include "Rodin/Context/Local.h"

#include "ForwardDecls.h"

#include "Types.h"
#include "Polytope.h"
#include "GeometryIndexed.h"

namespace Rodin::Geometry
{
  /**
   * @brief Abstract base class for mesh connectivity.
   *
   * Defines the interface for accessing mesh topology and incidence relations
   * between polytopes of different dimensions.
   */
  class ConnectivityBase
  {
    public:
      /**
       * @brief Gets the geometry type of a polytope.
       * @param[in] d Dimension of the polytope
       * @param[in] idx Index of the polytope
       * @returns Geometry type (Triangle, Tetrahedron, etc.)
       */
      virtual Polytope::Type getGeometry(size_t d, Index idx) const = 0;

      /**
       * @brief Gets the vertex indices defining a polytope.
       * @param[in] d Dimension of the polytope
       * @param[in] idx Index of the polytope
       * @returns Array of vertex indices
       */
      virtual const Array<Index>& getPolytope(size_t d, Index idx) const = 0;

      /**
       * @brief Gets the count of polytopes in a dimension.
       * @param[in] dim Dimension
       * @returns Number of polytopes of dimension @p dim
       */
      virtual size_t getCount(size_t dim) const = 0;

      /**
       * @brief Gets the count of polytopes of a specific geometry.
       * @param[in] g Geometry type
       * @returns Number of polytopes of type @p g
       */
      virtual size_t getCount(Polytope::Type g) const = 0;

      /**
       * @brief Gets the topological dimension of the mesh.
       * @returns Maximal dimension @f$ D @f$ of polytopes in the mesh
       */
      virtual size_t getDimension() const = 0;

      /**
       * @brief Gets the full incidence relation @f$ d \rightarrow d' @f$.
       * @param[in] d Source dimension
       * @param[in] dp Target dimension
       * @returns Incidence relation data structure
       */
      virtual const Incidence& getIncidence(size_t d, size_t dp) const = 0;

      /**
       * @brief Gets incident polytopes for a specific polytope.
       * @param[in] p Pair (d, idx) specifying the polytope
       * @param[in] idx Index of the polytope
       * @returns Vector of incident polytope indices
       */
      virtual const IndexVector& getIncidence(const std::pair<size_t, size_t> p, Index idx) const = 0;
  };

  /**
   * @brief Type alias for sequential (non-distributed) connectivity.
   */
  using SequentialConnectivity = Connectivity<Context::Local>;

  /**
   * @brief Mesh connectivity for sequential (local) meshes.
   *
   * This class stores and manages the topological relationships between
   * polytopes in a mesh. It maintains:
   * - Polytope-to-vertex mappings
   * - Incidence relations @f$ d \rightarrow d' @f$ between polytopes
   * - Geometry types for all polytopes
   *
   * # Mathematical Foundation
   *
   * For a mesh with topological dimension @f$ D @f$, the connectivity stores
   * all incidence relations:
   * @f[
   *  d \longrightarrow d', \quad 0 \leq d, d' \leq D
   * @f]
   * where the relation @f$ d \rightarrow d' @f$ maps each @f$ d @f$-dimensional
   * polytope to its incident @f$ d' @f$-dimensional polytopes.
   *
   * # Usage
   *
   * Connectivity must be computed before use:
   * @code{.cpp}
   * mesh.getConnectivity().compute(1, 2); // Compute face-to-cell incidence
   * @endcode
   *
   * # Storage Efficiency
   *
   * The implementation is based on the efficient connectivity storage scheme
   * described in @cite logg2009efficient, which minimizes memory usage while
   * enabling fast queries.
   *
   * # Thread Safety
   *
   * Connectivity objects are not thread-safe during construction. Once
   * finalized, read-only access is thread-safe.
   *
   * @see ConnectivityBase, Mesh
   */
  template <>
  class Connectivity<Context::Local> final : public ConnectivityBase
  {
    friend class boost::serialization::access;

    public:
      /**
       * @brief Bidirectional index mapping for polytope identification.
       *
       * Maintains a bidirectional mapping between vertex arrays (defining
       * polytopes) and their indices in the mesh.
       */
      struct PolytopeIndex
      {
        friend class boost::serialization::access;

        public:
          std::vector<const IndexArray*> left;  ///< Index to vertex array mapping
          UnorderedMap<IndexArray, Index, IndexArraySymmetricHash, IndexArraySymmetricEquality> right;  ///< Vertex array to index mapping

          /**
           * @brief Serialization save method.
           * @param[in,out] ar Archive object
           */
          template <class Archive>
          void save(Archive& ar, const unsigned int /*version*/) const
          {
            ar & right;
            std::vector<IndexArray> left_keys;
            left_keys.reserve(left.size());
            for (const IndexArray* p : left)
              left_keys.push_back(*p);
            ar & left_keys;
          }

          /**
           * @brief Serialization load method.
           * @param[in,out] ar Archive object
           */
          template <class Archive>
          void load(Archive& ar, const unsigned int /*version*/)
          {
            ar & right;
            std::vector<IndexArray> left_keys;
            ar & left_keys;
            left.clear();
            left.reserve(left_keys.size());
            for (const auto& k : left_keys)
            {
              auto it = right.find(k);
              assert(it != right.end());
              left.push_back(&it->first);
            }
          }

          BOOST_SERIALIZATION_SPLIT_MEMBER()
      };

      /**
       * @brief Represents a sub-polytope (lower-dimensional face).
       */
      struct SubPolytope
      {
        Polytope::Type geometry;  ///< Geometry type of the sub-polytope
        Array<Index> vertices;    ///< Vertex indices defining the sub-polytope
      };

      /**
       * @brief Default constructor.
       */
      Connectivity();

      /**
       * @brief Copy constructor.
       */
      Connectivity(const Connectivity&) = default;

      /**
       * @brief Move constructor.
       */
      Connectivity(Connectivity&&) = default;

      /**
       * @brief Copy assignment operator.
       */
      Connectivity& operator=(const Connectivity&) = default;

      /**
       * @brief Move assignment operator.
       */
      Connectivity& operator=(Connectivity&&) = default;

      /**
       * @brief Initializes connectivity for a mesh of given dimension.
       * @param[in] maximalDimension Topological dimension @f$ D @f$ of the mesh
       * @returns Reference to this object for method chaining
       *
       * Must be called before adding polytopes or computing incidence relations.
       */
      Connectivity& initialize(size_t maximalDimension);

      /**
       * @brief Sets the number of vertices in the mesh.
       * @param[in] count Number of vertices (0-dimensional polytopes)
       * @returns Reference to this object for method chaining
       */
      Connectivity& nodes(size_t count);

      /**
       * @brief Clears an incidence relation.
       * @param[in] d Source dimension
       * @param[in] dp Target dimension
       * @returns Reference to this object for method chaining
       *
       * Clears the stored incidence relation @f$ d \rightarrow d' @f$.
       */
      Connectivity& clear(size_t d, size_t dp);

      /**
       * @brief Reserves storage space for polytopes.
       * @param[in] d Dimension of polytopes
       * @param[in] count Number of polytopes to reserve space for
       * @returns Reference to this object for method chaining
       *
       * Pre-allocates memory to avoid repeated reallocations during construction.
       */
      Connectivity& reserve(size_t d, size_t count);

      /**
       * @brief Adds a polytope from an initializer list.
       * @param[in] t Geometry type
       * @param[in] p Initializer list of vertex indices
       * @returns Reference to this object for method chaining
       */
      Connectivity& polytope(
          Geometry::Polytope::Type t, std::initializer_list<Index> p)
      {
        Array<Index> arr(p.size());
        std::copy(p.begin(), p.end(), arr.begin());
        return polytope(t, std::move(arr));
      }

      /**
       * @brief Adds a polytope to the mesh.
       * @param[in] t Geometry type
       * @param[in] polytope Array of vertex indices defining the polytope
       * @returns Reference to this object for method chaining
       */
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

      /**
       * @brief Extracts local connectivity for a specific polytope.
       * @param[in] i Polytope index
       * @param[in] d Dimension of the polytope
       * @returns Reference to this connectivity object
       */
      Connectivity& local(size_t i, size_t d);

      /**
       * @brief Computes connectivity between dimensions.
       * @param[in] d Source dimension
       * @param[in] dp Target dimension
       * @returns Reference to this connectivity object
       *
       * Computes the incidence relation from polytopes of dimension @p d
       * to polytopes of dimension @p dp.
       */
      Connectivity& compute(size_t d, size_t dp);

      /**
       * @brief Transposes a connectivity relation.
       * @param[in] d First dimension
       * @param[in] dp Second dimension
       * @returns Reference to this connectivity object
       *
       * Given connectivity from dimension @p d to @p dp, computes the
       * transpose relation from @p dp to @p d.
       */
      Connectivity& transpose(size_t d, size_t dp);

      /**
       * @brief Computes connectivity through intersection of two relations.
       * @param[in] d Source dimension
       * @param[in] dp Intermediate dimension
       * @param[in] dpp Target dimension
       * @returns Reference to this connectivity object
       *
       * Computes connectivity from @p d to @p dpp via @p dp, i.e.,
       * the composition of @f$ d \to dp @f$ and @f$ dp \to dpp @f$.
       */
      Connectivity& intersection(size_t d, size_t dp, size_t dpp);

      /**
       * @brief Gets all sub-polytopes of a given polytope.
       * @param[out] out Vector to store sub-polytopes
       * @param[in] i Polytope index
       * @param[in] d Polytope dimension
       */
      void getSubPolytopes(std::vector<SubPolytope>& out, Index i, size_t d) const;

      /**
       * @brief Gets the index mapping for a given dimension.
       * @param[in] dim Dimension
       * @returns Reference to the polytope index map for the dimension
       */
      const PolytopeIndex& getIndexMap(size_t dim) const;

      /**
       * @brief Gets the index for a polytope defined by a set of vertices.
       * @param[in] dim Dimension of the polytope
       * @param[in] key Array of vertex indices defining the polytope
       * @returns Optional index if found, empty otherwise
       */
      const Optional<Index> getIndex(size_t dim, const IndexArray& key) const;

      /**
       * @brief Sets the incidence relation for a dimension pair.
       * @param[in] p Dimension pair (from, to)
       * @param[in] inc Incidence relation to set
       * @returns Reference to this connectivity object
       */
      Connectivity& setIncidence(const std::pair<size_t, size_t>& p, Incidence&& inc);

      size_t getCount(size_t dim) const override;

      size_t getCount(Polytope::Type g) const override;

      size_t getDimension() const override;

      Polytope::Type getGeometry(size_t d, Index idx) const override;

      const Array<Index>& getPolytope(size_t d, Index idx) const override;

      const Incidence& getIncidence(size_t d, size_t dp) const override;

      const IndexVector& getIncidence(const std::pair<size_t, size_t> p, Index idx) const override;

      template<class Archive>
      void serialize(Archive& ar, const unsigned int version)
      {
        ar & m_maximalDimension;
        ar & m_count;
        ar & m_gcount;
        ar & m_index;
        ar & m_dirty;
        ar & m_geometry;
        ar & m_connectivity;
      }

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
