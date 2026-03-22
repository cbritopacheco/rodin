/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cassert>
#include <cstring>

#include <t8.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_element_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_vec.h>

#include "Rodin/Alert.h"

#include "Mesh.h"

namespace Rodin::T8Code
{
  // --------------------------------------------------------------------------
  // Helper: hash 3D coordinates for vertex deduplication
  // --------------------------------------------------------------------------
  static
  size_t hashCoordinates(const double coords[3])
  {
    size_t hash = 0;
    for (size_t d = 0; d < 3; ++d)
    {
      auto bits = static_cast<int64_t>(coords[d] * 1e10);
      hash ^= std::hash<int64_t>{}(bits) + 0x9e3779b9
               + (hash << 6) + (hash >> 2);
    }
    return hash;
  }

  // --------------------------------------------------------------------------
  // Helper: map Rodin polytope type -> t8code element class
  // --------------------------------------------------------------------------
  static
  t8_eclass_t polytopeToEclass(Geometry::Polytope::Type g)
  {
    switch (g)
    {
      case Geometry::Polytope::Type::Point:
        return T8_ECLASS_VERTEX;
      case Geometry::Polytope::Type::Segment:
        return T8_ECLASS_LINE;
      case Geometry::Polytope::Type::Triangle:
        return T8_ECLASS_TRIANGLE;
      case Geometry::Polytope::Type::Quadrilateral:
        return T8_ECLASS_QUAD;
      case Geometry::Polytope::Type::Tetrahedron:
        return T8_ECLASS_TET;
      case Geometry::Polytope::Type::Hexahedron:
        return T8_ECLASS_HEX;
      case Geometry::Polytope::Type::Wedge:
        return T8_ECLASS_PRISM;
      default:
        Alert::Exception()
          << "Unsupported polytope type for t8code conversion."
          << Alert::Raise;
        return T8_ECLASS_COUNT; // unreachable
    }
  }

  // --------------------------------------------------------------------------
  // Helper: map t8code element class -> Rodin polytope type
  // --------------------------------------------------------------------------
  static
  Geometry::Polytope::Type eclassToPolytope(t8_eclass_t eclass)
  {
    switch (eclass)
    {
      case T8_ECLASS_VERTEX:
        return Geometry::Polytope::Type::Point;
      case T8_ECLASS_LINE:
        return Geometry::Polytope::Type::Segment;
      case T8_ECLASS_TRIANGLE:
        return Geometry::Polytope::Type::Triangle;
      case T8_ECLASS_QUAD:
        return Geometry::Polytope::Type::Quadrilateral;
      case T8_ECLASS_TET:
        return Geometry::Polytope::Type::Tetrahedron;
      case T8_ECLASS_HEX:
        return Geometry::Polytope::Type::Hexahedron;
      case T8_ECLASS_PRISM:
        return Geometry::Polytope::Type::Wedge;
      default:
        Alert::Exception()
          << "Unsupported t8code element class."
          << Alert::Raise;
        return Geometry::Polytope::Type::Point; // unreachable
    }
  }

  // --------------------------------------------------------------------------
  // Construction / destruction
  // --------------------------------------------------------------------------

  Mesh::Mesh()
    : m_forest(nullptr),
      m_cmesh(nullptr),
      m_maxRefinementLevel(0)
  {}

  Mesh::Mesh(Parent&& mesh)
    : Parent(std::move(mesh)),
      m_forest(nullptr),
      m_cmesh(nullptr),
      m_maxRefinementLevel(0)
  {
    initializeFromMesh(*this);
    rebuild();
  }

  Mesh::Mesh(const Mesh& other)
    : Parent(other),
      m_forest(nullptr),
      m_cmesh(nullptr),
      m_refinementLevels(other.m_refinementLevels),
      m_maxRefinementLevel(other.m_maxRefinementLevel),
      m_hangingNodes(other.m_hangingNodes),
      m_constraints(other.m_constraints)
  {
    if (other.m_cmesh != nullptr)
    {
      t8_cmesh_ref(other.m_cmesh);
      m_cmesh = other.m_cmesh;
    }
    if (other.m_forest != nullptr)
    {
      t8_forest_ref(other.m_forest);
      m_forest = other.m_forest;
    }
  }

  Mesh::Mesh(Mesh&& other)
    : Parent(std::move(other)),
      m_forest(other.m_forest),
      m_cmesh(other.m_cmesh),
      m_refinementLevels(std::move(other.m_refinementLevels)),
      m_maxRefinementLevel(other.m_maxRefinementLevel),
      m_hangingNodes(std::move(other.m_hangingNodes)),
      m_constraints(std::move(other.m_constraints))
  {
    other.m_forest = nullptr;
    other.m_cmesh = nullptr;
    other.m_maxRefinementLevel = 0;
  }

  Mesh& Mesh::operator=(const Mesh& other)
  {
    if (this != &other)
    {
      cleanup();
      Parent::operator=(other);
      m_refinementLevels = other.m_refinementLevels;
      m_maxRefinementLevel = other.m_maxRefinementLevel;
      m_hangingNodes = other.m_hangingNodes;
      m_constraints = other.m_constraints;
      if (other.m_cmesh != nullptr)
      {
        t8_cmesh_ref(other.m_cmesh);
        m_cmesh = other.m_cmesh;
      }
      else
      {
        m_cmesh = nullptr;
      }
      if (other.m_forest != nullptr)
      {
        t8_forest_ref(other.m_forest);
        m_forest = other.m_forest;
      }
      else
      {
        m_forest = nullptr;
      }
    }
    return *this;
  }

  Mesh& Mesh::operator=(Mesh&& other)
  {
    if (this != &other)
    {
      cleanup();
      Parent::operator=(std::move(other));
      m_forest = other.m_forest;
      m_cmesh = other.m_cmesh;
      m_refinementLevels = std::move(other.m_refinementLevels);
      m_maxRefinementLevel = other.m_maxRefinementLevel;
      m_hangingNodes = std::move(other.m_hangingNodes);
      m_constraints = std::move(other.m_constraints);
      other.m_forest = nullptr;
      other.m_cmesh = nullptr;
      other.m_maxRefinementLevel = 0;
    }
    return *this;
  }

  Mesh& Mesh::operator=(Parent&& other)
  {
    cleanup();
    Parent::operator=(std::move(other));
    m_refinementLevels.clear();
    m_maxRefinementLevel = 0;
    m_hangingNodes.clear();
    m_constraints.clear();
    initializeFromMesh(*this);
    rebuild();
    return *this;
  }

  Mesh::~Mesh()
  {
    cleanup();
  }

  void Mesh::cleanup()
  {
    if (m_forest != nullptr)
    {
      t8_forest_unref(&m_forest);
      m_forest = nullptr;
    }
    if (m_cmesh != nullptr)
    {
      t8_cmesh_unref(&m_cmesh);
      m_cmesh = nullptr;
    }
  }

  // --------------------------------------------------------------------------
  // Initialization from a Rodin mesh
  // --------------------------------------------------------------------------

  void Mesh::initializeFromMesh(const Parent& mesh)
  {
    cleanup();

    const size_t numCells = mesh.getCellCount();
    if (numCells == 0)
      return;

    const size_t sdim = mesh.getSpaceDimension();

    // Initialize coarse mesh
    t8_cmesh_init(&m_cmesh);

    // Set tree classes from mesh cells
    for (Index i = 0; i < static_cast<Index>(numCells); ++i)
    {
      auto it = mesh.getCell(i);
      const auto geom = it->getGeometry();
      t8_cmesh_set_tree_class(m_cmesh, i, polytopeToEclass(geom));
    }

    // Set vertex coordinates for each tree
    for (Index i = 0; i < static_cast<Index>(numCells); ++i)
    {
      auto it = mesh.getCell(i);
      const auto& verts = it->getVertices();
      const size_t nv = verts.size();

      // t8code expects 3D coordinates; pad with zeros if needed
      std::vector<double> coords(nv * 3, 0.0);
      for (size_t v = 0; v < nv; ++v)
      {
        auto vcoords = mesh.getVertexCoordinates(verts[v]);
        for (size_t d = 0; d < sdim && d < 3; ++d)
          coords[v * 3 + d] = vcoords[d];
      }
      t8_cmesh_set_tree_vertices(m_cmesh, i, coords.data(), nv);
    }

    // Commit the coarse mesh
    // Note: sc_MPI_COMM_WORLD resolves to a serial communicator when
    // t8code is built without MPI support.
    t8_cmesh_commit(m_cmesh, sc_MPI_COMM_WORLD);

    // Create the forest from the coarse mesh
    t8_forest_init(&m_forest);
    t8_forest_set_cmesh(m_forest, m_cmesh);
    t8_forest_set_scheme(m_forest, t8_scheme_new_default_cxx());
    t8_forest_set_level(m_forest, 0);
    t8_forest_commit(m_forest);

    // Take a reference to the cmesh so it outlives the forest
    t8_cmesh_ref(m_cmesh);

    m_maxRefinementLevel = 0;
  }

  // --------------------------------------------------------------------------
  // Rebuild the Rodin mesh from the t8code forest
  // --------------------------------------------------------------------------

  void Mesh::rebuild()
  {
    if (m_forest == nullptr)
      return;

    const t8_locidx_t numTrees = t8_forest_get_num_local_trees(m_forest);
    const size_t sdim = getSpaceDimension() > 0 ? getSpaceDimension() : 3;

    // Collect all element vertices and build the Rodin mesh
    std::unordered_map<size_t, Index> vertexMap;
    std::vector<Math::SpatialPoint> vertices;
    std::vector<std::pair<Geometry::Polytope::Type, std::vector<Index>>> elements;

    m_refinementLevels.clear();
    m_maxRefinementLevel = 0;

    for (t8_locidx_t itree = 0; itree < numTrees; ++itree)
    {
      const t8_eclass_t eclass = t8_forest_get_tree_class(m_forest, itree);
      const t8_eclass_scheme_c* scheme =
        t8_forest_get_eclass_scheme(m_forest, eclass);
      const t8_locidx_t numElements =
        t8_forest_get_tree_num_elements(m_forest, itree);

      for (t8_locidx_t ielem = 0; ielem < numElements; ++ielem)
      {
        const t8_element_t* element =
          t8_forest_get_element_in_tree(m_forest, itree, ielem);
        const int level = scheme->t8_element_level(element);
        const int numCorners = scheme->t8_element_num_corners(element);

        if (static_cast<size_t>(level) > m_maxRefinementLevel)
          m_maxRefinementLevel = static_cast<size_t>(level);

        m_refinementLevels.push_back(static_cast<size_t>(level));

        // Get the vertex coordinates of this element
        std::vector<Index> elemVertices(numCorners);
        for (int ic = 0; ic < numCorners; ++ic)
        {
          double coords[3] = {0.0, 0.0, 0.0};
          t8_forest_element_coordinate(
            m_forest, itree, element, ic, coords);

          const size_t hash = hashCoordinates(coords);

          auto ins = vertexMap.insert({hash, static_cast<Index>(vertices.size())});
          if (ins.second)
          {
            Math::SpatialPoint pt(sdim);
            for (size_t d = 0; d < sdim && d < 3; ++d)
              pt[d] = coords[d];
            vertices.push_back(std::move(pt));
          }
          elemVertices[ic] = ins.first->second;
        }

        auto polyType = eclassToPolytope(eclass);
        elements.emplace_back(polyType, std::move(elemVertices));
      }
    }

    // Rebuild the parent mesh using the Builder
    auto builder = Parent::Builder();
    builder.initialize(sdim);
    builder.nodes(vertices.size());
    for (auto& v : vertices)
      builder.vertex(v);
    for (auto& [geom, vs] : elements)
      builder.polytope(geom, vs);

    Parent::operator=(builder.finalize());

    // Recompute hanging node information
    computeHangingNodes();
  }

  // --------------------------------------------------------------------------
  // Hanging node detection
  // --------------------------------------------------------------------------

  void Mesh::computeHangingNodes()
  {
    m_hangingNodes.clear();
    m_constraints.clear();

    if (m_forest == nullptr)
      return;

    const t8_locidx_t numTrees = t8_forest_get_num_local_trees(m_forest);
    const size_t sdim = getSpaceDimension() > 0 ? getSpaceDimension() : 3;

    // Build a coordinate-to-vertex-index map for the rebuilt Rodin mesh
    std::unordered_map<size_t, Index> coordToVertex;
    for (Index vi = 0; vi < static_cast<Index>(getVertexCount()); ++vi)
    {
      auto vcoords = getVertexCoordinates(vi);
      double c[3] = {0.0, 0.0, 0.0};
      for (size_t d = 0; d < sdim && d < 3; ++d)
        c[d] = vcoords[d];
      coordToVertex[hashCoordinates(c)] = vi;
    }

    // For each tree and each element, check face neighbors for level
    // differences indicating hanging nodes
    for (t8_locidx_t itree = 0; itree < numTrees; ++itree)
    {
      const t8_eclass_t eclass = t8_forest_get_tree_class(m_forest, itree);
      const t8_eclass_scheme_c* scheme =
        t8_forest_get_eclass_scheme(m_forest, eclass);
      const t8_locidx_t numElements =
        t8_forest_get_tree_num_elements(m_forest, itree);

      for (t8_locidx_t ielem = 0; ielem < numElements; ++ielem)
      {
        const t8_element_t* element =
          t8_forest_get_element_in_tree(m_forest, itree, ielem);
        const int level = scheme->t8_element_level(element);
        const int numFaces = scheme->t8_element_num_faces(element);

        for (int iface = 0; iface < numFaces; ++iface)
        {
          // Query face neighbors
          int numNeighbors = 0;
          t8_element_t** neighbors = nullptr;
          t8_eclass_scheme_c* neighScheme = nullptr;
          int* faceIndices = nullptr;
          t8_locidx_t neighTree;
          t8_eclass_t neighClass;

          t8_forest_leaf_face_neighbors(
            m_forest, itree, element, &neighbors, iface,
            &faceIndices, &numNeighbors, &neighTree,
            &neighScheme, 1);

          if (numNeighbors > 0 && neighScheme != nullptr)
          {
            const int neighLevel = neighScheme->t8_element_level(neighbors[0]);

            // If neighbor is coarser, the midpoint of the face on the
            // fine side is a hanging node
            if (neighLevel < level)
            {
              const int numFaceCorners =
                scheme->t8_element_num_face_corners(element, iface);

              // The midpoint of an edge face (2 corners) is a hanging node
              if (numFaceCorners == 2)
              {
                double c0[3] = {}, c1[3] = {};
                const int fc0 = scheme->t8_element_get_face_corner(
                  element, iface, 0);
                const int fc1 = scheme->t8_element_get_face_corner(
                  element, iface, 1);

                t8_forest_element_coordinate(
                  m_forest, itree, element, fc0, c0);
                t8_forest_element_coordinate(
                  m_forest, itree, element, fc1, c1);

                // Compute midpoint
                double mid[3];
                for (int d = 0; d < 3; ++d)
                  mid[d] = 0.5 * (c0[d] + c1[d]);

                // Find the midpoint vertex in the Rodin mesh
                const size_t midHash = hashCoordinates(mid);

                auto midIt = coordToVertex.find(midHash);
                if (midIt != coordToVertex.end())
                {
                  // Find the coarse edge endpoints on the neighbor
                  double nc0[3] = {}, nc1[3] = {};
                  const int neighFaceCorners =
                    neighScheme->t8_element_num_face_corners(
                      neighbors[0], faceIndices[0]);

                  if (neighFaceCorners >= 2)
                  {
                    const int nfc0 = neighScheme->t8_element_get_face_corner(
                      neighbors[0], faceIndices[0], 0);
                    const int nfc1 = neighScheme->t8_element_get_face_corner(
                      neighbors[0], faceIndices[0], 1);

                    t8_forest_element_coordinate(
                      m_forest, neighTree, neighbors[0], nfc0, nc0);
                    t8_forest_element_coordinate(
                      m_forest, neighTree, neighbors[0], nfc1, nc1);

                    const size_t h0 = hashCoordinates(nc0);
                    const size_t h1 = hashCoordinates(nc1);

                    auto it0 = coordToVertex.find(h0);
                    auto it1 = coordToVertex.find(h1);
                    if (it0 != coordToVertex.end() &&
                        it1 != coordToVertex.end())
                    {
                      m_hangingNodes.insert(midIt->second);
                      m_constraints[midIt->second] = {
                        it0->second, it1->second
                      };
                    }
                  }
                }
              }
            }
          }

          // Clean up neighbor allocations
          if (neighbors != nullptr)
          {
            for (int n = 0; n < numNeighbors; ++n)
              neighScheme->t8_element_destroy(1, &neighbors[n]);
            T8_FREE(neighbors);
          }
          if (faceIndices != nullptr)
            T8_FREE(faceIndices);
        }
      }
    }
  }

  // --------------------------------------------------------------------------
  // Refinement
  // --------------------------------------------------------------------------

  // Adaptation callback context
  struct AdaptContext
  {
    const Mesh* mesh;
    std::function<bool(const Geometry::Cell&)> predicate;
    bool uniformRefine;
    bool coarsen;
  };

  // t8code adaptation callback
  static int
  adaptCallback(
    t8_forest_t forest,
    t8_forest_t forestFrom,
    t8_locidx_t whichTree,
    t8_locidx_t lelement_id,
    t8_eclass_scheme_c* ts,
    const int isFamily,
    const int numElements,
    t8_element_t* elements[])
  {
    (void)forest;
    (void)ts;
    (void)numElements;

    const auto* ctx = static_cast<const AdaptContext*>(
      t8_forest_get_user_data(forestFrom));

    if (ctx->uniformRefine)
      return 1; // Refine

    if (ctx->coarsen && isFamily)
    {
      // For coarsening, the predicate must return true for the first element
      // of the family (representative of the parent)
      const t8_locidx_t globalIdx =
        t8_forest_get_tree_element_offset(forestFrom, whichTree) + lelement_id;
      if (static_cast<size_t>(globalIdx) < ctx->mesh->getCellCount())
      {
        auto cellIt = ctx->mesh->getCell(static_cast<Index>(globalIdx));
        if (ctx->predicate(*cellIt))
          return -1; // Coarsen
      }
      return 0;
    }

    if (!ctx->coarsen)
    {
      const t8_locidx_t globalIdx =
        t8_forest_get_tree_element_offset(forestFrom, whichTree) + lelement_id;
      if (static_cast<size_t>(globalIdx) < ctx->mesh->getCellCount())
      {
        auto cellIt = ctx->mesh->getCell(static_cast<Index>(globalIdx));
        if (ctx->predicate(*cellIt))
          return 1; // Refine
      }
    }

    return 0; // Do nothing
  }

  Mesh& Mesh::refine(std::function<bool(const Geometry::Cell&)> predicate)
  {
    if (m_forest == nullptr)
    {
      Alert::Exception()
        << "Cannot refine: no forest initialized."
        << Alert::Raise;
    }

    AdaptContext ctx;
    ctx.mesh = this;
    ctx.predicate = std::move(predicate);
    ctx.uniformRefine = false;
    ctx.coarsen = false;

    t8_forest_t newForest;
    t8_forest_init(&newForest);
    t8_forest_set_user_data(m_forest, &ctx);
    t8_forest_set_adapt(newForest, m_forest, adaptCallback, 0);
    t8_forest_commit(newForest);

    m_forest = newForest;
    rebuild();
    return *this;
  }

  Mesh& Mesh::refine()
  {
    if (m_forest == nullptr)
    {
      Alert::Exception()
        << "Cannot refine: no forest initialized."
        << Alert::Raise;
    }

    AdaptContext ctx;
    ctx.mesh = this;
    ctx.uniformRefine = true;
    ctx.coarsen = false;

    t8_forest_t newForest;
    t8_forest_init(&newForest);
    t8_forest_set_user_data(m_forest, &ctx);
    t8_forest_set_adapt(newForest, m_forest, adaptCallback, 0);
    t8_forest_commit(newForest);

    m_forest = newForest;
    rebuild();
    return *this;
  }

  // --------------------------------------------------------------------------
  // Coarsening
  // --------------------------------------------------------------------------

  Mesh& Mesh::coarsen(std::function<bool(const Geometry::Cell&)> predicate)
  {
    if (m_forest == nullptr)
    {
      Alert::Exception()
        << "Cannot coarsen: no forest initialized."
        << Alert::Raise;
    }

    AdaptContext ctx;
    ctx.mesh = this;
    ctx.predicate = std::move(predicate);
    ctx.uniformRefine = false;
    ctx.coarsen = true;

    t8_forest_t newForest;
    t8_forest_init(&newForest);
    t8_forest_set_user_data(m_forest, &ctx);
    t8_forest_set_adapt(newForest, m_forest, adaptCallback, 0);
    t8_forest_commit(newForest);

    m_forest = newForest;
    rebuild();
    return *this;
  }

  // --------------------------------------------------------------------------
  // Balance
  // --------------------------------------------------------------------------

  Mesh& Mesh::balance()
  {
    if (m_forest == nullptr)
    {
      Alert::Exception()
        << "Cannot balance: no forest initialized."
        << Alert::Raise;
    }

    t8_forest_t newForest;
    t8_forest_init(&newForest);
    t8_forest_set_balance(newForest, m_forest, 0);
    t8_forest_commit(newForest);

    m_forest = newForest;
    rebuild();
    return *this;
  }

  // --------------------------------------------------------------------------
  // Queries
  // --------------------------------------------------------------------------

  size_t Mesh::getRefinementLevel(Index cellIdx) const
  {
    assert(static_cast<size_t>(cellIdx) < m_refinementLevels.size());
    return m_refinementLevels[static_cast<size_t>(cellIdx)];
  }

  size_t Mesh::getMaxRefinementLevel() const
  {
    return m_maxRefinementLevel;
  }

  bool Mesh::isHangingNode(Index vertexIdx) const
  {
    return m_hangingNodes.count(vertexIdx) > 0;
  }

  FlatSet<Index> Mesh::getHangingNodes() const
  {
    return m_hangingNodes;
  }

  std::pair<Index, Index>
  Mesh::getConstrainingVertices(Index vertexIdx) const
  {
    auto it = m_constraints.find(vertexIdx);
    assert(it != m_constraints.end());
    return it->second;
  }

  t8_forest_t Mesh::getForest() const
  {
    return m_forest;
  }

  t8_cmesh_t Mesh::getCoarseMesh() const
  {
    return m_cmesh;
  }
}
