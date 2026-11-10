/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_ORDERING_H
#define RODIN_GEOMETRY_ORDERING_H

#include "SpaceFillingCurveOrdering.h"

namespace Rodin::Geometry
{
  /// @brief Available local mesh vertex ordering algorithms.
  enum class OrderingType
  {
    Morton, ///< Morton (Z-curve) ordering.
    Hilbert, ///< Hilbert space-filling curve ordering.
    RCM ///< Reverse Cuthill-McKee graph ordering.
  };

  /**
   * @brief Vertex ordering policy selected by its algorithm.
   * @tparam Type Ordering algorithm; Morton is the default.
   *
   * Each specialization is callable on a local mesh and returns a new-to-old
   * vertex permutation, accepted by Mesh::renumber().
   *
   * | Specialization | Configuration |
   * | --- | --- |
   * | OrderingType::Morton | Bounding-box quantization bit depth |
   * | OrderingType::Hilbert | Bounding-box quantization bit depth |
   * | OrderingType::RCM | Deterministic roots and degree/index tie breaking |
   *
   * @code{.cpp}
   * Ordering<OrderingType::Hilbert> ordering;
   * ordering.setBitsPerCoordinate(21);
   * auto reordered = mesh.renumber(ordering);
   * @endcode
   */
  template <OrderingType Type = OrderingType::Morton>
  class Ordering;

  /**
   * @brief Morton (Z-curve) ordering of local mesh vertices.
   *
   * Integer coordinate bits are interleaved, most significant first.
   * Bounding-box normalization, configurable bit depth and stable tie handling
   * are provided by Internal::SpaceFillingCurveOrdering. Construction uses
   * automatic precision. Apply before constructing finite element spaces or
   * transferring a mesh to a device.
   */
  template <>
  class Ordering<OrderingType::Morton>
    : public Internal::SpaceFillingCurveOrdering<Ordering<OrderingType::Morton>>
  {
    friend class Internal::SpaceFillingCurveOrdering<Ordering<OrderingType::Morton>>;

    private:
      static void transform(std::vector<std::uint64_t>&, size_t) {}
  };

  /**
   * @brief Hilbert ordering of local mesh vertices.
   *
   * Integer coordinates are transformed using Skilling's axes-to-transpose
   * algorithm before their bits are interleaved (doi:10.1063/1.1751381).
   * Bounding-box normalization, configurable bit depth and stable tie handling
   * are provided by Internal::SpaceFillingCurveOrdering. Construction uses
   * automatic precision. Apply before constructing finite element spaces or
   * transferring a mesh to a device.
   */
  template <>
  class Ordering<OrderingType::Hilbert>
    : public Internal::SpaceFillingCurveOrdering<Ordering<OrderingType::Hilbert>>
  {
    friend class Internal::SpaceFillingCurveOrdering<Ordering<OrderingType::Hilbert>>;

    private:
      static void transform(std::vector<std::uint64_t>& coordinates, size_t bits);
  };

  /**
   * @brief Deterministic reverse Cuthill-McKee ordering of local mesh vertices.
   *
   * Vertices sharing a maximal-dimensional cell form a clique in the graph.
   * Each connected component starts from a pseudo-peripheral vertex. Initial
   * seeds, last-level candidates and BFS neighbors are selected by increasing
   * degree, with original vertex index breaking ties. The BFS ordering is
   * reversed separately in each component; components are visited by their
   * minimum (degree, original index) seed. Isolated vertices are included.
   *
   * This is a bandwidth-reduction heuristic for the vertex graph, not an
   * ordering of high-order or mixed-space degrees of freedom. No derived
   * incidence is required or computed, and the input mesh is unchanged.
   */
  template <>
  class Ordering<OrderingType::RCM>
  {
    public:
      /**
       * @brief Computes the reverse Cuthill-McKee vertex ordering.
       * @param[in] mesh Local mesh supplying cell-to-vertex lists.
       * @returns Permutation with newToOld[newIndex] equal to the old vertex index.
       */
      IndexVector operator()(const Mesh<Context::Local>& mesh) const;
  };

  /// @brief Default construction selects the Morton ordering.
  Ordering() -> Ordering<OrderingType::Morton>;

}

#endif
