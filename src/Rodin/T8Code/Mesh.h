/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Mesh.h
 * @brief Nonconforming AMR mesh type backed by the t8code library.
 *
 * This file defines the @ref Rodin::T8Code::Mesh class which extends
 * @ref Rodin::Geometry::Mesh with adaptive mesh refinement (AMR) capabilities
 * provided by t8code. The mesh supports:
 * - Local cell refinement and coarsening
 * - 2:1 balance enforcement
 * - Hanging node tracking and constraint queries
 * - Refinement level queries per cell
 *
 * t8code manages the underlying forest-of-trees representation while this
 * class exposes refined meshes as standard Rodin mesh objects for use in
 * finite element computations.
 */
#ifndef RODIN_T8CODE_MESH_H
#define RODIN_T8CODE_MESH_H

#include <vector>
#include <functional>
#include <unordered_map>

#include <t8.h>
#include <t8_cmesh.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>

#include "Rodin/Types.h"
#include "Rodin/Geometry.h"

#include "ForwardDecls.h"

namespace Rodin::T8Code
{
  /**
   * @brief Nonconforming AMR mesh backed by t8code.
   *
   * This class extends @ref Rodin::Geometry::Mesh<Context::Local> with
   * adaptive mesh refinement capabilities provided by t8code's
   * forest-of-trees data structure. It supports:
   *
   * - **Local refinement**: Refine individual cells based on user predicates
   * - **Coarsening**: Merge families of sibling cells back into parents
   * - **2:1 balance**: Ensure neighboring cells differ by at most one
   *   refinement level
   * - **Hanging node tracking**: Query which vertices are hanging nodes
   *   and retrieve their constraining vertices
   * - **Refinement hierarchy**: Query the refinement level of each cell
   *
   * # Usage
   *
   * @code{.cpp}
   * #include <Rodin/Geometry.h>
   * #include <Rodin/T8Code/Mesh.h>
   *
   * using namespace Rodin;
   * using namespace Rodin::Geometry;
   *
   * // Create a base mesh
   * Mesh base = Mesh<Context::Local>::UniformGrid(
   *   Polytope::Type::Quadrilateral, {5, 5});
   *
   * // Create an AMR mesh from the base mesh
   * T8Code::Mesh amr(std::move(base));
   *
   * // Refine cells near a point
   * amr.refine([](const Cell& cell) {
   *   auto v = cell.getVertex();
   *   auto coords = v->getCoordinates();
   *   return coords.norm() < 0.3;
   * });
   *
   * // Enforce 2:1 balance
   * amr.balance();
   * @endcode
   *
   * @see <a href="https://github.com/DLR-AMR/t8code">t8code</a>
   */
  class Mesh : public Geometry::Mesh<Context::Local>
  {
    public:
      /// Parent class type.
      using Parent = Geometry::Mesh<Rodin::Context::Local>;

      /// Mesh context type.
      using Context = typename Parent::Context;

      /**
       * @brief Constructs an empty T8Code mesh.
       */
      Mesh();

      /**
       * @brief Move-constructs a T8Code mesh from a base Rodin mesh.
       *
       * Converts the given mesh into a t8code coarse mesh (cmesh) and
       * initializes a uniform forest at refinement level 0.
       *
       * @param[in] mesh Base mesh to convert.
       */
      explicit Mesh(Parent&& mesh);

      /**
       * @brief Copy constructor.
       */
      Mesh(const Mesh& other);

      /**
       * @brief Move constructor.
       */
      Mesh(Mesh&& other);

      /**
       * @brief Copy assignment.
       */
      Mesh& operator=(const Mesh& other);

      /**
       * @brief Move assignment.
       */
      Mesh& operator=(Mesh&& other);

      /**
       * @brief Move-assigns from a base local mesh.
       *
       * Converts the mesh into a t8code coarse mesh and reinitializes the
       * forest at refinement level 0.
       *
       * @param[in] other Base mesh to convert.
       * @returns Reference to this mesh.
       */
      Mesh& operator=(Parent&& other);

      /**
       * @brief Destructor.
       *
       * Releases the underlying t8code forest and coarse mesh resources.
       */
      ~Mesh() override;

      /**
       * @brief Refines cells matching a predicate.
       *
       * Cells for which the predicate returns @c true are refined using
       * t8code's element-specific refinement rules (e.g., 1:4 for triangles
       * and quadrilaterals, 1:8 for tetrahedra and hexahedra).
       *
       * After refinement, the underlying Rodin mesh is rebuilt from the
       * updated t8code forest. Calling code must re-create finite element
       * spaces and grid functions on the new mesh.
       *
       * @param[in] predicate Function returning @c true for cells to refine.
       * @returns Reference to this mesh.
       */
      Mesh& refine(std::function<bool(const Geometry::Cell&)> predicate);

      /**
       * @brief Uniformly refines all cells by one level.
       * @returns Reference to this mesh.
       */
      Mesh& refine();

      /**
       * @brief Coarsens cells matching a predicate.
       *
       * Families of sibling cells for which the predicate returns @c true on
       * every member are merged back into their parent cell.
       *
       * @param[in] predicate Function returning @c true for cells to coarsen.
       * @returns Reference to this mesh.
       */
      Mesh& coarsen(std::function<bool(const Geometry::Cell&)> predicate);

      /**
       * @brief Enforces 2:1 balance on the forest.
       *
       * Ensures that neighboring cells differ by at most one refinement
       * level. Additional cells may be refined to satisfy this constraint.
       *
       * @returns Reference to this mesh.
       */
      Mesh& balance();

      /**
       * @brief Returns the refinement level of a cell.
       * @param[in] cellIdx Cell index.
       * @returns Refinement level (0 for coarse-mesh cells).
       */
      size_t getRefinementLevel(Index cellIdx) const;

      /**
       * @brief Returns the maximum refinement level across all cells.
       * @returns Maximum refinement level in the forest.
       */
      size_t getMaxRefinementLevel() const;

      /**
       * @brief Tests whether a vertex is a hanging node.
       *
       * A hanging node is a vertex that lies on an edge or face of a
       * neighboring element at a coarser refinement level.
       *
       * @param[in] vertexIdx Vertex index.
       * @returns @c true if the vertex is a hanging node.
       */
      bool isHangingNode(Index vertexIdx) const;

      /**
       * @brief Returns all hanging node vertex indices.
       * @returns Set of vertex indices that are hanging nodes.
       */
      FlatSet<Index> getHangingNodes() const;

      /**
       * @brief Returns the pair of constraining vertices for a hanging node.
       *
       * A hanging node on an edge is constrained by the two endpoints of
       * that edge on the coarser neighbor. The hanging node's value in a
       * conforming finite element space must be the average of these two
       * vertices' values.
       *
       * @param[in] vertexIdx Hanging node vertex index.
       * @returns Pair of vertex indices that constrain the hanging node.
       * @pre @p vertexIdx must be a hanging node.
       */
      std::pair<Index, Index> getConstrainingVertices(Index vertexIdx) const;

      /**
       * @brief Returns the t8code forest handle.
       * @returns Pointer to the underlying t8code forest.
       */
      t8_forest_t getForest() const;

      /**
       * @brief Returns the t8code coarse mesh handle.
       * @returns Pointer to the underlying t8code coarse mesh.
       */
      t8_cmesh_t getCoarseMesh() const;

    private:
      /**
       * @brief Converts a Rodin mesh into a t8code coarse mesh.
       */
      void initializeFromMesh(const Parent& mesh);

      /**
       * @brief Rebuilds the Rodin mesh from the current t8code forest.
       *
       * Extracts vertex coordinates and element connectivity from the
       * t8code forest and populates the parent Mesh data structures.
       * Also recomputes hanging node information.
       */
      void rebuild();

      /**
       * @brief Recomputes the hanging node set and constraint map.
       */
      void computeHangingNodes();

      /**
       * @brief Releases t8code resources.
       */
      void cleanup();

      t8_forest_t  m_forest;  ///< t8code forest handle.
      t8_cmesh_t   m_cmesh;   ///< t8code coarse mesh handle.

      std::vector<size_t> m_refinementLevels; ///< Per-cell refinement level.
      size_t m_maxRefinementLevel;            ///< Maximum refinement level.

      FlatSet<Index> m_hangingNodes; ///< Set of hanging node vertex indices.

      /// Map from hanging node vertex index to constraining vertex pair.
      std::unordered_map<Index, std::pair<Index, Index>> m_constraints;
  };
}

#endif
