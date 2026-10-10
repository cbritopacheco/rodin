/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_MINSTCUT_H
#define RODIN_GEOMETRY_MINSTCUT_H

#include <functional>
#include <vector>

#include "Rodin/Types.h"
#include "ForwardDecls.h"

namespace Rodin::Geometry
{
  /**
   * @brief Binary Potts classifier selected by mesh type.
   * @tparam MeshType Mesh type, including its execution context.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref MinSTCut<Mesh<Context::Local>> | Serial cut on a local mesh. |
   *
   * @section minstcut_usage Usage
   * @code{.cpp}
   * const size_t d = mesh.getDimension();
   * mesh.getConnectivity().compute(d - 1, d);
   * MinSTCut classifier(mesh);
   * decltype(classifier)::Parameters parameters;
   * parameters.fidelity = 1;
   * parameters.smoothing = [](const Geometry::Polytope& facet)
   * {
   *   return facet.getAttribute().value_or(0) == 2 ? Real(0.1) : Real(1);
   * };
   * classifier.setParameters(parameters);
   * // Return the signed field average over each cell, not its integral.
   * const auto result = classifier.classify([&](const Geometry::Polytope& cell)
   * {
   *   return averages[cell.getIndex()];
   * });
   * // Applying labels to mesh attributes is left to the caller.
   * @endcode
   *
   * @section minstcut_model Model
   * For cell volumes @f$V_i@f$, signed cell averages @f$m_i@f$ and interior
   * facets @f$\sigma=K_i\cap K_j@f$, classification minimizes
   * @f[
   * E(\ell)=\eta\sum_i U_i(\ell_i)
   *   +\sum_{\sigma=K_i\cap K_j}s(\sigma)|\sigma|
   *     \boldsymbol1_{\ell_i\ne\ell_j},
   * \qquad U_i(-1)=V_i(m_i)_+,\quad U_i(+1)=V_i(-m_i)_+.
   * @f]
   * Here @f$\eta@f$ is fidelity and @f$s@f$ is facet-dependent smoothing.
   * Negative averages favor Inside; positive averages favor Outside.
   * A typical average is @f$m_i=V_i^{-1}\int_{K_i}\tanh(\phi/\varepsilon)\,dx@f$.
   * The supplied callable computes this cell statistic; the classifier does
   * not choose a field, phase map or integration rule on the caller's behalf.
   * Boundary facets contribute no pairwise term and no labels are pinned.
   * In one dimension the facet measure is the counting measure, namely one.
   *
   * @section minstcut_architecture Architecture
   * The classifier borrows a @ref Mesh and reads its face-to-cell
   * @ref Connectivity. Cell and facet measures are obtained from
   * @ref Polytope::getMeasure; no geometry or adjacency is reconstructed.
   * Interior facets become graph connections between their two incident
   * cells. A private Dinic residual network computes the minimum cut.
   * Results use mesh cell and facet indices; neither mesh attributes nor
   * topology are modified. The network and geometric weights are rebuilt
   * for each classification, so mesh edits are not hidden by a stale cache.
   *
   * @note The borrowed mesh must outlive the classifier. Face-to-cell
   * incidence must be explicitly computed before classification. Only local
   * meshes with at most two cells incident to a facet are supported;
   * distributed minimum cuts are not implemented.
   */
  template <class MeshType>
  class MinSTCut;

  /**
   * @brief Serial mesh-aware Potts classifier using local connectivity.
   * @see MinSTCut
   */
  template <>
  class MinSTCut<Mesh<Context::Local>> final
  {
    public:

      /// @brief Mesh type supported by this serial classifier.
      using MeshType = Mesh<Context::Local>;

      /// @brief Cell classification and the separating mesh facets.
      struct Result
      {
          /// Mesh indices of inside cells.
          IndexVector inside;
          /// Mesh indices of outside cells.
          IndexVector outside;
          /// Mesh indices of interior facets separating different labels.
          IndexVector cut;
          /// Value of the minimized Potts objective.
          Real energy = 0;
      };

      /// @brief Fidelity and geometric perimeter weighting.
      struct Parameters
      {
          /// Finite nonnegative multiplier applied to every unary cost.
          Real fidelity = 1;

          /**
           * @brief Finite nonnegative multiplier of an interior facet's measure.
           * Evaluated once per interior facet; defaults to one. Capturing
           * lambdas are supported. Reference captures must remain valid for
           * every classification using these parameters.
           */
          std::function<Real(const Polytope&)> smoothing =
            [](const Polytope&) -> Real { return 1; };
      };

      /**
       * @brief Constructs a classifier borrowing a mesh with default parameters.
       * @param mesh Local mesh, which must outlive the classifier.
       */
      explicit MinSTCut(const MeshType& mesh);

      /**
       * @brief Gets the borrowed mesh.
       * @returns Mesh used to construct adjacency and geometric costs.
       */
      const MeshType& getMesh() const
      {
        return m_mesh.get();
      }

      /**
       * @brief Sets the owned classification parameters.
       * @param parameters Parameters copied into this object.
       * @returns This classifier for chaining.
       */
      MinSTCut& setParameters(const Parameters& parameters)
      {
        m_parameters = parameters;
        return *this;
      }

      /**
       * @brief Gets the owned classification parameters.
       * @returns Parameters used without an explicit per-call override.
       */
      const Parameters& getParameters() const
      {
        return m_parameters;
      }

      /**
       * @brief Classifies mesh cells using the owned parameters.
       * @param average Callable returning a finite signed field average for
       * each cell polytope. It is evaluated once per cell and is not retained.
       * @returns Inside/outside cell indices, cut facets and minimized objective.
       */
      Result classify(const std::function<Real(const Polytope&)>& average) const;

    private:
      /**
       * @brief Dinic maximum-flow algorithm for the private cell graph.
       *
       * Every directed connection has a reverse residual arc. Breadth-first
       * traversal builds source-distance levels; augmentation follows only
       * arcs advancing one level. Blocking flows are repeated until the sink
       * is unreachable. The final levels identify the source side of the
       * minimum cut, without a second graph traversal.
       */
      class Dinic
      {
        public:
          /// @brief Whether a connection carries one or two directed capacities.
          enum class Type
          {
            Directed, ///< Capacity from the first vertex to the second only.
            Undirected ///< Equal directed capacities in both directions.
          };

          /**
           * @brief Allocates an empty residual network.
           * @param count Number of vertices, including source and sink.
           */
          explicit Dinic(Index count);

          /**
           * @brief Adds a connection and its reverse residual arcs.
           * @param first First graph vertex.
           * @param second Second graph vertex.
           * @param capacity Finite nonnegative connection capacity.
           * @param type Directed or undirected connection.
           */
          void add(Index first, Index second, Real capacity, Type type);

          /**
           * @brief Computes maximum flow in the current residual network.
           * @param source Source vertex.
           * @param sink Sink vertex, distinct from the source.
           * @returns Flow added by this call.
           */
          Real getMaximumFlow(Index source, Index sink);

          /**
           * @brief Checks source reachability after maximum-flow completion.
           * @param vertex Graph vertex to inspect.
           * @returns Whether the vertex belongs to the source-side cut.
           */
          Boolean isReachable(Index vertex) const;

        private:
          /// @brief Residual arc paired with its reverse arc.
          struct Arc
          {
              Index to; ///< Destination vertex.
              Index reverse; ///< Index of the reverse arc at the destination.
              Real capacity; ///< Remaining residual capacity.
          };

          /**
           * @brief Builds source-distance levels along positive residual arcs.
           * @param source Source vertex.
           * @param sink Sink vertex.
           * @returns Whether the sink is reachable.
           */
          Boolean build(Index source, Index sink);

          /**
           * @brief Augments flow along the current level graph.
           * @param current Current vertex on the augmenting path.
           * @param sink Sink vertex.
           * @param amount Maximum additional flow allowed on the path.
           * @returns Flow pushed along the path, or zero if blocked.
           */
          Real augment(Index current, Index sink, Real amount);

          std::vector<std::vector<Arc>> m_graph;
          std::vector<Integer> m_level;
          IndexVector m_next;
      };

      /// @brief Interior mesh facet represented in the private cell graph.
      struct Edge
      {
          Index first; ///< First incident cell index.
          Index second; ///< Second incident cell index.
          Real capacity; ///< Smoothing times the facet measure.
          Index face; ///< Mesh facet index.
      };

      /**
       * @brief Solves the assembled cell graph and translates its cut to mesh indices.
       * @param insideCosts Costs of assigning cells to Inside.
       * @param outsideCosts Costs of assigning cells to Outside.
       * @param edges Interior facets with their weighted capacities.
       * @returns Inside/outside sets, cut facets and energy.
       */
      Result solve(const std::vector<Real>& insideCosts,
        const std::vector<Real>& outsideCosts, const std::vector<Edge>& edges) const;

      std::reference_wrapper<const MeshType> m_mesh;
      Parameters m_parameters;
  };

  /**
   * @brief Deduces the classifier specialization from the mesh.
   * @tparam MeshType Mesh type, including its execution context.
   * @param mesh Mesh borrowed by the classifier.
   */
  template <class MeshType>
  MinSTCut(const MeshType& mesh) -> MinSTCut<MeshType>;

}

#endif
