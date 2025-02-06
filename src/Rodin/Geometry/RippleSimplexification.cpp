#include "GeometryIndexed.h"
#include "RippleSimplexification.h"

namespace Rodin::Geometry
{
  GeometryIndexed<std::vector<std::vector<std::vector<IndexArray>>>>
  RippleSimplexification<Mesh<Context::Local>>::Simplexification::s_simplices =
  {
    // =====================================================
    // Polytope: Triangle
    // =====================================================
    {
      Polytope::Type::Triangle,
      {
        // Level 0: 0d simplices (vertices)
        {
          { IndexArray{{0}},
            IndexArray{{1}},
            IndexArray{{2}}
          }
        },
        // Level 1: 1d simplices (edges)
        {
          { IndexArray{{0, 1}},
            IndexArray{{1, 2}},
            IndexArray{{2, 0}}
          }
        },
        // Level 2: 2d simplex (face)
        {
          { IndexArray{{0, 1, 2}} }
        },
      },
    },
    // =====================================================
    // Polytope: Quadrilateral
    // =====================================================
    {
      Polytope::Type::Quadrilateral,
      {
        // -----------------------------------------------------
        // Triangulation Variant 1:
        // Diagonal chosen from Vertex 0 to Vertex 3.
        // -----------------------------------------------------
        {
          // Level 0: 0d simplices (vertices)
          {
            // The quadrilateral always has four vertices:
            { IndexArray{{0}},
              IndexArray{{1}},
              IndexArray{{2}},
              IndexArray{{3}}
            }
          },
          // Level 1: 1d simplices (edges)
          {
            IndexArray{{0, 1}},
            IndexArray{{1, 3}},
            IndexArray{{3, 2}},
            IndexArray{{2, 0}},
            IndexArray{{0, 3}},  // Interior (diagonal) edge
          },
          // Level 2: 2d simplices (faces)
          {
            IndexArray{{0, 1, 3}},
            IndexArray{{0, 3, 2}},
          },
        },
        // -----------------------------------------------------
        // Triangulation Variant 2:
        // Diagonal chosen from Vertex 1 to Vertex 2.
        // -----------------------------------------------------
        {
          // Level 0: 0d simplices (vertices)
          {
            { IndexArray{{0}},
              IndexArray{{1}},
              IndexArray{{2}},
              IndexArray{{3}}
            }
          },
          // Level 1: 1d simplices (edges)
          {
            IndexArray{{0, 1}},
            IndexArray{{1, 3}},
            IndexArray{{3, 2}},
            IndexArray{{2, 0}},
            IndexArray{{1, 2}}, // Interior (diagonal) edge
          },
          // Level 2: 2d simplices (faces)
          {
            IndexArray{{0, 1, 2}},
            IndexArray{{1, 3, 2}},
          },
        },
      }
    },
    // =====================================================
    // Polytope: Triangular Prism
    // =====================================================
    {
      Polytope::Type::TriangularPrism,
      {
        // =============================================================
        // Triangulation 1:
        // Tetrahedra: {0,1,2,4}, {0,2,4,5}, {0,3,4,5}
        // =============================================================
        {
          // 0d vertices
          {
            { IndexArray{{0}} },
            { IndexArray{{1}} },
            { IndexArray{{2}} },
            { IndexArray{{3}} },
            { IndexArray{{4}} },
            { IndexArray{{5}} },
          },
          // 1d edges (first the fixed boundary edges, then interior)
          {
            // Boundary edges:
            IndexArray{{0,1}}, IndexArray{{1,2}}, IndexArray{{2,0}},
            IndexArray{{3,4}}, IndexArray{{4,5}}, IndexArray{{5,3}},
            IndexArray{{0,3}}, IndexArray{{1,4}}, IndexArray{{2,5}},
            // Interior edges:
            IndexArray{{0,4}}, IndexArray{{2,4}}, IndexArray{{0,5}},
          },
          // 2d faces
          {
            // Faces from T0 = {0,1,2,4}:
            IndexArray{{0,1,2}}, IndexArray{{0,1,4}}, IndexArray{{1,2,4}}, IndexArray{{0,2,4}},
            // Faces from T1 = {0,2,4,5}:
            IndexArray{{0,2,5}}, IndexArray{{2,4,5}}, IndexArray{{4,0,5}},
            // Faces from T2 = {0,3,4,5}:
            IndexArray{{0,3,4}}, IndexArray{{0,3,5}}, IndexArray{{3,4,5}},
          },
          // 3d tetrahedra
          {
            IndexArray{{0,1,2,4}},
            IndexArray{{0,2,4,5}},
            IndexArray{{0,3,4,5}},
          },
        },
        // =============================================================
        // Triangulation 2:
        // Tetrahedra: {0,1,2,5}, {0,1,4,5}, {0,3,4,5}
        // =============================================================
        {
          // 0d vertices (same as above)
          {
            { IndexArray{{0}} },
            { IndexArray{{1}} },
            { IndexArray{{2}} },
            { IndexArray{{3}} },
            { IndexArray{{4}} },
            { IndexArray{{5}} },
          },
          // 1d edges
          {
            IndexArray{{0,1}}, IndexArray{{1,2}}, IndexArray{{2,0}},
            IndexArray{{3,4}}, IndexArray{{4,5}}, IndexArray{{5,3}},
            IndexArray{{0,3}}, IndexArray{{1,4}}, IndexArray{{2,5}},
            IndexArray{{0,5}}, IndexArray{{0,4}}, IndexArray{{1,5}},
          },
          // 2d faces
          {
            // Faces from T0 = {0,1,2,5}:
            IndexArray{{0,1,2}}, IndexArray{{0,1,5}}, IndexArray{{1,2,5}}, IndexArray{{2,0,5}},
            // Faces from T1 = {0,1,4,5}:
            IndexArray{{0,1,4}}, IndexArray{{0,1,5}}, IndexArray{{1,4,5}}, IndexArray{{4,0,5}},
            // Faces from T2 = {0,3,4,5}:
            IndexArray{{0,3,4}}, IndexArray{{0,3,5}}, IndexArray{{3,4,5}},
          },
          // 3d tetrahedra
          {
            IndexArray{{0,1,2,5}},
            IndexArray{{0,1,4,5}},
            IndexArray{{0,3,4,5}},
          },
        },
        // =============================================================
        // Triangulation 3:
        // Tetrahedra: {0,1,2,4}, {1,2,3,4}, {1,3,4,5}
        // =============================================================
        {
          {
            { IndexArray{{0}} },
            { IndexArray{{1}} },
            { IndexArray{{2}} },
            { IndexArray{{3}} },
            { IndexArray{{4}} },
            { IndexArray{{5}} },
          },
          {
            // Boundary edges:
            IndexArray{{0,1}}, IndexArray{{1,2}}, IndexArray{{2,0}},
            IndexArray{{3,4}}, IndexArray{{4,5}}, IndexArray{{5,3}},
            IndexArray{{0,3}}, IndexArray{{1,4}}, IndexArray{{2,5}},
            // Interior edges:
            IndexArray{{0,4}}, IndexArray{{2,4}},
            IndexArray{{1,3}}, IndexArray{{2,3}},
            IndexArray{{1,5}}, IndexArray{{3,5}},
          },
          {
            // T0 = {0,1,2,4}:
            IndexArray{{0,1,2}}, IndexArray{{0,1,4}}, IndexArray{{1,2,4}}, IndexArray{{0,2,4}},
            // T1 = {1,2,3,4}:
            IndexArray{{1,2,3}}, IndexArray{{1,3,4}}, IndexArray{{2,3,4}}, IndexArray{{1,2,4}},
            // T2 = {1,3,4,5}:
            IndexArray{{1,3,4}}, IndexArray{{1,4,5}}, IndexArray{{3,4,5}}, IndexArray{{1,3,5}},
          },
          {
            IndexArray{{0,1,2,4}},
            IndexArray{{1,2,3,4}},
            IndexArray{{1,3,4,5}},
          },
        },
        // =============================================================
        // Triangulation 4:
        // Tetrahedra: {0,1,2,5}, {1,2,3,5}, {1,3,4,5}
        // =============================================================
        {
          {
            { IndexArray{{0}} },
            { IndexArray{{1}} },
            { IndexArray{{2}} },
            { IndexArray{{3}} },
            { IndexArray{{4}} },
            { IndexArray{{5}} },
          },
          {
            IndexArray{{0,1}}, IndexArray{{1,2}}, IndexArray{{2,0}},
            IndexArray{{3,4}}, IndexArray{{4,5}}, IndexArray{{5,3}},
            IndexArray{{0,3}}, IndexArray{{1,4}}, IndexArray{{2,5}},
            // Interior edges:
            IndexArray{{0,5}}, IndexArray{{1,5}},
            IndexArray{{1,3}}, IndexArray{{2,3}},
            IndexArray{{3,5}},
          },
          {
            // T0 = {0,1,2,5}:
            IndexArray{{0,1,2}}, IndexArray{{0,1,5}}, IndexArray{{1,2,5}}, IndexArray{{0,2,5}},
            // T1 = {1,2,3,5}:
            IndexArray{{1,2,3}}, IndexArray{{1,3,5}}, IndexArray{{2,3,5}}, IndexArray{{1,2,5}},
            // T2 = {1,3,4,5}:
            IndexArray{{1,3,4}}, IndexArray{{1,4,5}}, IndexArray{{3,4,5}}, IndexArray{{1,3,5}},
          },
          {
            IndexArray{{0,1,2,5}},
            IndexArray{{1,2,3,5}},
            IndexArray{{1,3,4,5}},
          },
        },
        // =============================================================
        // Triangulation 5:
        // Tetrahedra: {0,1,2,3}, {0,2,3,4}, {0,3,4,5}
        // =============================================================
        {
          {
            { IndexArray{{0}} },
            { IndexArray{{1}} },
            { IndexArray{{2}} },
            { IndexArray{{3}} },
            { IndexArray{{4}} },
            { IndexArray{{5}} },
          },
          {
            IndexArray{{0,1}}, IndexArray{{1,2}}, IndexArray{{2,0}},
            IndexArray{{3,4}}, IndexArray{{4,5}}, IndexArray{{5,3}},
            IndexArray{{0,3}}, IndexArray{{1,4}}, IndexArray{{2,5}},
            // Interior edges:
            IndexArray{{1,3}}, IndexArray{{2,3}},
            IndexArray{{0,4}}, IndexArray{{2,4}},
            IndexArray{{0,5}}, IndexArray{{3,5}},
          },
          {
            // T0 = {0,1,2,3}:
            IndexArray{{0,1,2}}, IndexArray{{0,1,3}}, IndexArray{{1,2,3}}, IndexArray{{0,2,3}},
            // T1 = {0,2,3,4}:
            IndexArray{{0,2,3}}, IndexArray{{0,2,4}}, IndexArray{{2,3,4}}, IndexArray{{0,3,4}},
            // T2 = {0,3,4,5}:
            IndexArray{{0,3,4}}, IndexArray{{0,3,5}}, IndexArray{{3,4,5}}, IndexArray{{0,4,5}},
          },
          {
            IndexArray{{0,1,2,3}},
            IndexArray{{0,2,3,4}},
            IndexArray{{0,3,4,5}},
          },
        },
        // =============================================================
        // Triangulation 6:
        // Tetrahedra: {0,1,2,3}, {1,2,3,5}, {1,3,4,5}
        // =============================================================
        {
          {
            { IndexArray{{0}} },
            { IndexArray{{1}} },
            { IndexArray{{2}} },
            { IndexArray{{3}} },
            { IndexArray{{4}} },
            { IndexArray{{5}} },
          },
          {
            IndexArray{{0,1}}, IndexArray{{1,2}}, IndexArray{{2,0}},
            IndexArray{{3,4}}, IndexArray{{4,5}}, IndexArray{{5,3}},
            IndexArray{{0,3}}, IndexArray{{1,4}}, IndexArray{{2,5}},
            // Interior edges:
            IndexArray{{1,3}}, IndexArray{{2,3}},
            IndexArray{{1,5}}, IndexArray{{3,5}},
          },
          {
            // T0 = {0,1,2,3}:
            IndexArray{{0,1,2}}, IndexArray{{0,1,3}}, IndexArray{{1,2,3}}, IndexArray{{0,2,3}},
            // T1 = {1,2,3,5}:
            IndexArray{{1,2,3}}, IndexArray{{1,3,5}}, IndexArray{{2,3,5}}, IndexArray{{1,2,5}},
            // T2 = {1,3,4,5}:
            IndexArray{{1,3,4}}, IndexArray{{1,4,5}}, IndexArray{{3,4,5}}, IndexArray{{1,3,5}},
          },
          {
            IndexArray{{0,1,2,3}},
            IndexArray{{1,2,3,5}},
            IndexArray{{1,3,4,5}},
          },
        },
      }
    }
  };

  constexpr
  size_t RippleSimplexification<Mesh<Context::Local>>::Simplexification::getCount() const
  {
    switch (getGeometry())
    {
      case Polytope::Type::Point:
      case Polytope::Type::Segment:
      case Polytope::Type::Triangle:
      case Polytope::Type::Tetrahedron:
        return 1;
      case Polytope::Type::Quadrilateral:
        return 2;
      case Polytope::Type::TriangularPrism:
        return 6;
    }
    assert(false);
    return 0;
  }

  constexpr
  size_t RippleSimplexification<Mesh<Context::Local>>::Simplexification::getCount(size_t, size_t d) const
  {
    switch (getGeometry())
    {
      case Polytope::Type::Point:
      {
        assert(d == 0);
        return 1;
      }
      case Polytope::Type::Segment:
      {
        if (d == 0)
        {
          return 2;
        }
        else if (d == 1)
        {
          return 1;
        }
        else
        {
          assert(false);
          return 0;
        }
      }
      case Polytope::Type::Triangle:
      {
        if (d == 0)
        {
          return 3;
        }
        else if (d == 1)
        {
          return 3;
        }
        else if (d == 2)
        {
          return 1;
        }
        else
        {
          assert(false);
          return 0;
        }
      }
      case Polytope::Type::Tetrahedron:
      {
        if (d == 0)
        {
          return 4;
        }
        else if (d == 1)
        {
          return 6;
        }
        else if (d == 2)
        {
          return 4;
        }
        else if (d == 3)
        {
          return 1;
        }
        else
        {
          assert(false);
          return 0;
        }
      }
      case Polytope::Type::Quadrilateral:
      {
        if (d == 0)
        {
          return 4;
        }
        else if (d == 1)
        {
          return 5;
        }
        else if (d == 2)
        {
          return 2;
        }
        else
        {
          assert(false);
          return 0;
        }
      }
      case Polytope::Type::TriangularPrism:
      {
        if (d == 0)
        {
          return 6;
        }
        else if (d == 1)
        {
          return 12;
        }
        else if (d == 2)
        {
          return 8;
        }
        else
        {
          assert(false);
          return 0;
        }
      }
    }
    assert(false);
    return 0;
  }

  constexpr
  const IndexArray& RippleSimplexification<Mesh<Context::Local>>::Simplexification::getSimplex(size_t k, size_t d, size_t i) const
  {
    return s_simplices[getGeometry()][k][d][i];
  }

  constexpr
  Polytope::Type RippleSimplexification<Mesh<Context::Local>>::Simplexification::getGeometry() const
  {
    return m_geometry;
  }

  RippleSimplexification<Mesh<Context::Local>>
  ::RippleSimplexification(const Mesh<Context::Local>& mesh)
    : Parent(mesh)
  {}

  Mesh<Context::Local> RippleSimplexification<Mesh<Context::Local>>::simplexify()
  {
    // const auto& mesh = getMesh();
    // const size_t D = mesh.getDimension();
    // RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, D, D - 1);
    // RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, D - 1, D);
    // std::deque<Index> searchQueue;
    // std::vector<size_t> visitedCount(mesh.getFaceCount(), 0);
    // m_triangulation.resize(mesh.getFaceCount(), 0);
    // for (auto it = mesh.getFace(); it; ++it)
    //   searchQueue.push_back(it->getIndex());
    // while (searchQueue.size() > 0)
    // {
    //   const Index i = searchQueue.back();
    //   searchQueue.pop_back();
    //   visitedCount[i] += 1;
    //   auto it = mesh.getFace(i);

    //   // Check if the face triangulation is compatible with the adjacent cells
    //   bool compatible;

    //   if (it->getGeometry() == Polytope::Type::Triangle)
    //     continue;

    //   const auto& cells = mesh.getConnectivity().getIncidence({D - 1, D}, i);

    //   for (const auto& j : cells)
    //   {
    //     auto cell = mesh.getCell(j);
    //     IndexMap<size_t> count;
    //     count.reserve(6);
    //     switch (cell->getGeometry())
    //     {
    //       case Polytope::Type::TriangularPrism:
    //       {
    //         compatible = false;
    //         const auto& faces = mesh.getConnectivity().getIncidence({D, D - 1}, j);
    //         for (const Index k : faces)
    //         {
    //           const auto& triangulation = m_triangulation[k];
    //           const auto& vs = mesh.getConnectivity().getPolytope(D - 1, k);
    //           const auto& simplices = s_simplices[mesh.getGeometry(D - 1, k)].at(triangulation);
    //           for (const auto& simplex : simplices)
    //           {
    //             for (const Index sv : simplex)
    //               count[vs[sv]] += 1;
    //           }
    //         }

    //         for (const auto& [k, v] : count)
    //         {
    //           // One vertex must be shared by at least 5 triangles
    //           if (v >= 5)
    //           {
    //             compatible = true;
    //             break;
    //           }
    //         }
    //         break;
    //       }
    //     }

    //     if (!compatible)
    //       break;
    //   }

    //   // If not compatible, then change the triangulation and add the faces in
    //   // the incident cells to the search queue
    //   if (!compatible)
    //   {
    //     m_triangulation[i] = visitedCount[i] % s_simplices[it->getGeometry()].size();
    //     for (const auto& j : cells)
    //     {
    //       auto cell = mesh.getCell(j);
    //       const auto& faces = mesh.getConnectivity().getIncidence({D, D - 1}, j);
    //       for (const Index k : faces)
    //       {
    //         if (k == i)
    //           continue;
    //         searchQueue.push_back(k);
    //       }
    //     }
    //   }
    // }

    // // Build the new mesh
    // for (auto it = mesh.getCell(); it; ++it)
    // {
    //   const auto& cell = *it;
    //   const auto& faces = mesh.getConnectivity().getIncidence({D, D - 1}, cell.getIndex());
    //   for (const Index i : faces)
    //   {
    //     const auto& triangulation = m_triangulation[i];
    //     const auto& simplices = s_simplices[mesh.getGeometry(D - 1, i)].at(triangulation);
    //     for (const auto& simplex : simplices)
    //     {
    //     }
    //   }
    // }

    return m_build.finalize();
  }
}
