#ifndef RODIN_VARIATIONAL_H1_H1_HPP
#define RODIN_VARIATIONAL_H1_H1_HPP

#include <cstddef>

#include "Rodin/Alert/Raise.h"
#include "Rodin/Alert/Exception.h"

#include "H1.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Utility/ForConstexpr.h"

namespace Rodin::Variational
{
  /**
   * @brief Local cochain map for nodal H¹ DOFs on a given geometry.
   *
   * This class encodes, at compile time, how the nodal degrees of freedom
   * living on the boundary entities of a reference polytope @p G are
   * injected into the nodal DOFs of @p G itself.
   *
   * Concretely:
   *  - @p G is the codomain polytope (Segment, Triangle, Quadrilateral,
   *    Tetrahedron, or Wedge),
   *  - @p Domain is a container of DOF indices on the boundary entity,
   *  - @p Codomain is a container of DOF indices on the whole polytope @p G.
   *
   * The method Cochain::map writes appropriate entries of @p codomain from
   * the entries of @p domain, respecting:
   *  - the canonical reference-geometry vertex ordering in Polytope::Traits,
   *  - the edge / face orientation conventions in Connectivity::getSubPolytopes,
   *  - the nodal ordering used in:
   *      - GLL<K> / GLL01<K> on segments,
   *      - FeketeTriangle<K> on triangles,
   *      - FeketeTetrahedron<K> on tetrahedra,
   *      - tensor-product GLL nodes on quadrilaterals,
   *      - tensor-product (triangle × GLL in z) on wedges.
   *
   * From the algebraic viewpoint, this realizes the local cochain map
   * associated with the trace operator on the boundary of @p G, written
   * in the nodal basis.
   *
   * @tparam K     Polynomial degree of the H¹ element.
   * @tparam Scalar Scalar type of the H¹ field (Real, Complex, ...).
   * @tparam G     Codomain reference polytope type (Segment, Triangle, ...).
   */
  template <size_t K, class Scalar>
  template <Geometry::Polytope::Type G>
  class H1<K, Scalar, Geometry::Mesh<Context::Local>>::Cochain
  {
    public:

      /**
       * @name Local DOF counts per reference geometry
       * These constants describe how many nodal DOFs are attached to each
       * reference polytope for the degree @p K H¹ element.
       * @{
       */

      /// Number of nodal DOFs on a reference point (single vertex).
      static constexpr size_t PointCount = 1;

      /// Number of nodal DOFs on a reference segment (K+1 GLL nodes).
      static constexpr size_t SegmentCount = GLL<K>::Count;

      /// Number of nodal DOFs on a reference triangle (FeketeTriangle<K>).
      static constexpr size_t TriangleCount = FeketeTriangle<K>::Count;

      /// Number of nodal DOFs on a reference tetrahedron (FeketeTetrahedron<K>).
      static constexpr size_t TetrahedronCount = FeketeTetrahedron<K>::Count;

      /// Number of nodal DOFs on a reference quadrilateral ((K+1)×(K+1) GLL grid).
      static constexpr size_t QuadrilateralCount = (K + 1) * (K + 1);

      static constexpr size_t WedgeCount = (K + 1) * FeketeTriangle<K>::Count;

      /** @} */

      static constexpr size_t Count =
        []()
        {
          if constexpr (G == Geometry::Polytope::Type::Point)
            return PointCount;
          else if constexpr (G == Geometry::Polytope::Type::Segment)
            return SegmentCount;
          else if constexpr (G == Geometry::Polytope::Type::Triangle)
            return TriangleCount;
          else if constexpr (G == Geometry::Polytope::Type::Quadrilateral)
            return QuadrilateralCount;
          else if constexpr (G == Geometry::Polytope::Type::Tetrahedron)
            return TetrahedronCount;
          else if constexpr (G == Geometry::Polytope::Type::Wedge)
            return WedgeCount;
          else
            return 0;
        }();


      /**
       * @brief Injects boundary DOFs into the DOFs of the polytope @p G.
       *
       * This method maps the degrees of freedom associated with a boundary
       * entity of @p G (vertex, edge, or face) into the degrees of freedom
       * associated with the whole polytope @p G, for the Local-th boundary
       * entity.
       *
       * The specific pattern depends on:
       *  - the codomain geometry @p G,
       *  - the Local boundary index (edge/face index),
       *  - the polynomial degree @p K and the chosen nodal sets
       *    (GLL, Fekete, tensor-product).
       *
       * The mapping is consistent with:
       *  - the reference vertex coordinates in Polytope::Traits::getVertex,
       *  - the edge/face orientations used in Connectivity::getSubPolytopes,
       *  - the reference nodal ordering on each geometry (getNodes).
       *
       * @tparam Domain   Container type for the boundary DOFs (e.g. IndexArray).
       * @tparam Codomain Container type for the element DOFs (e.g. IndexArray).
       *
       *
       * @param[in]  local     Index of the boundary entity of @p G
       *                       (vertex, edge, or face, depending on the dimension)
       * @param[out] codomain  DOF container on @p G to be updated
       *                       along the Local-th boundary entity.
       * @param[in]  domain    DOF container on the Local-th boundary entity.
       */
      template <class Domain, class Codomain>
      static constexpr void map(size_t local, Codomain& codomain, const Domain& domain)
      {
        using Type = Geometry::Polytope::Type;

        // -------------------------------------------------------------------
        // G = Segment: from Point (vertex) -> Segment
        // Local = 0 or 1 = which endpoint
        // -------------------------------------------------------------------
        if constexpr (G == Type::Segment)
        {
          assert(local < 2 && "Segment has 2 vertices (Local = 0, 1).");
          if (local == 0)
          {
            // left endpoint
            codomain[0] = domain[0];
          }
          else if (local == 1)
          {
            // right endpoint
            codomain[SegmentCount - 1] = domain[0];
          }
        }
        // -------------------------------------------------------------------
        // G = Triangle: from Segment edge -> Triangle
        //
        // Triangle vertices: (0,0),(1,0),(0,1)
        // Edge locals:
        //   Local 0: (0->1)  bottom edge
        //   Local 1: (1->2)  "hypotenuse" (1,0)->(0,1)
        //   Local 2: (2->0)  left edge
        //
        // Triangle DOF ordering (FeketeTriangle<K>):
        //   for j2 = 0..K:
        //     rowStart = j2*(K+1) - j2*(j2-1)/2
        //     for i2 = 0..K-j2:
        //       idx = rowStart + i2
        // -------------------------------------------------------------------
        else if constexpr (G == Type::Triangle)
        {
          assert(local < 3 && "Triangle has 3 edges (Local = 0, 1, 2).");

          if (local == 0)
          {
            // bottom edge (0->1)
            Utility::ForIndex<SegmentCount>([&](auto ii)
            {
              constexpr size_t r   = ii.value; // 0..K
              constexpr size_t idx = r;        // bottom row contiguous
              codomain[idx] = domain[r];
            });
          }
          else if (local == 1)
          {
            // hypotenuse (1->2)
            Utility::ForIndex<SegmentCount>([&](auto ii)
            {
              constexpr size_t r = ii.value;   // 0..K
              constexpr size_t j = r;          // parameter along hypotenuse
              constexpr size_t rowStart_j =
                  j * (K + 1) - (j * (j - 1)) / 2;
              constexpr size_t idx = rowStart_j + (K - j);
              codomain[idx] = domain[r];
            });
          }
          else if (local == 2)
          {
            // left edge (2->0)
            Utility::ForIndex<SegmentCount>([&](auto ii)
            {
              constexpr size_t r_on_edge = ii.value;  // 0..K
              constexpr size_t j         = K - r_on_edge;
              constexpr size_t rowStart_j =
                  j * (K + 1) - (j * (j - 1)) / 2;
              constexpr size_t idx = rowStart_j;      // i = 0 on row j
              codomain[idx] = domain[r_on_edge];
            });
          }
        }
        // -------------------------------------------------------------------
        // G = Quadrilateral: from Segment edge -> Quadrilateral
        //
        // Quad reference vertices: (0,0),(1,0),(1,1),(0,1)
        // Edges (local):
        //   Local 0: (0->1) bottom
        //   Local 1: (1->2) right
        //   Local 2: (2->3) top
        //   Local 3: (3->0) left
        //
        // Quad DOF ordering:
        //   idx(i,j) = j*(K+1) + i, i,j in 0..K
        // -------------------------------------------------------------------
        else if constexpr (G == Type::Quadrilateral)
        {
          assert(local < 4 && "Quadrilateral has 4 edges (Local = 0, 1, 2, 3).");

          if (local == 0)
          {
            // bottom edge: j = 0, i = 0..K
            Utility::ForIndex<SegmentCount>([&](auto ii)
            {
              constexpr size_t r = ii.value;  // 0..K
              constexpr size_t i = r;
              constexpr size_t j = 0;
              constexpr size_t idx = j * (K + 1) + i;
              codomain[idx] = domain[r];
            });
          }
          else if (local == 1)
          {
            // right edge: i = K, j = 0..K
            Utility::ForIndex<SegmentCount>([&](auto ii)
            {
              constexpr size_t r = ii.value;  // 0..K
              constexpr size_t i = K;
              constexpr size_t j = r;
              constexpr size_t idx = j * (K + 1) + i;
              codomain[idx] = domain[r];
            });
          }
          else if (local == 2)
          {
            // top edge: j = K, i = K..0 (reverse to keep global orientation)
            Utility::ForIndex<SegmentCount>([&](auto ii)
            {
              constexpr size_t r = ii.value;  // 0..K
              constexpr size_t i = K - r;
              constexpr size_t j = K;
              constexpr size_t idx = j * (K + 1) + i;
              codomain[idx] = domain[r];
            });
          }
          else if (local == 3)
          {
            // left edge: i = 0, j = K..0 (reverse)
            Utility::ForIndex<SegmentCount>([&](auto ii)
            {
              constexpr size_t r = ii.value;  // 0..K
              constexpr size_t i = 0;
              constexpr size_t j = K - r;
              constexpr size_t idx = j * (K + 1) + i;
              codomain[idx] = domain[r];
            });
          }
        }
        // -------------------------------------------------------------------
        // G = Tetrahedron: from Triangle face -> Tetrahedron
        //
        // Tetra vertices: 0:(0,0,0), 1:(1,0,0), 2:(0,1,0), 3:(0,0,1)
        //
        // Faces (Connectivity::getSubPolytopes, dim=2):
        //   Local 0: (1,2,3)  // +[1,2,3]
        //   Local 1: (0,3,2)  // -[0,2,3]
        //   Local 2: (0,1,3)  // +[0,1,3]
        //   Local 3: (0,2,1)  // -[0,1,2]
        //
        // Triangle DOF ordering:
        //   triIdx(j2,i2) with
        //   rowStart(j2)=j2*(K+1) - j2*(j2-1)/2, i2=0..K-j2
        //
        // Tetra DOF ordering: equispaced i,j,k with i+j+k<=K; idx computed inline.
        // -------------------------------------------------------------------
        else if constexpr (G == Type::Tetrahedron)
        {
          assert(local < 4 && "Tetrahedron has 4 faces (Local = 0, 1, 2, 3).");

          if (local == 0)
          {
            // Face (1,2,3) opposite vertex 0: tri(0,1,2) -> tet(1,2,3)
            // λ1 = λ_tri0, λ2 = λ_tri1, λ3 = λ_tri2, λ0 = 0
            // (i,j,k) = (K - i2 - j2, i2, j2)
            Utility::ForIndex<K + 1>([&](auto jj)
            {
              constexpr size_t j2 = jj.value;
              Utility::ForIndex<K + 1>([&](auto ii)
              {
                constexpr size_t i2 = ii.value;
                if constexpr (i2 + j2 <= K)
                {
                  constexpr size_t triRowStart =
                      j2 * (K + 1) - (j2 * (j2 - 1)) / 2;
                  constexpr size_t triIdx = triRowStart + i2;
                  constexpr size_t i = K - i2 - j2;
                  constexpr size_t j = i2;
                  constexpr size_t k = j2;
                  constexpr size_t tetraTotal =
                      (K + 1) * (K + 2) * (K + 3) / 6;
                  constexpr size_t m_tail = K - k;
                  constexpr size_t tetraTail =
                      (m_tail + 1) * (m_tail + 2) * (m_tail + 3) / 6;
                  constexpr size_t offset_k = tetraTotal - tetraTail;
                  constexpr size_t offset_j =
                      j * (K - k + 1) - (j * (j - 1)) / 2;
                  constexpr size_t tetIdx = offset_k + offset_j + i;
                  codomain[tetIdx] = domain[triIdx];
                }
              });
            });
          }
          else if (local == 1)
          {
            // Face (0,3,2) opposite vertex 1: tri(0,1,2) -> tet(0,3,2)
            // λ0 = λ_tri0, λ3 = λ_tri1, λ2 = λ_tri2, λ1 = 0
            // (i,j,k) = (0, j2, i2)
            Utility::ForIndex<K + 1>([&](auto jj)
            {
              constexpr size_t j2 = jj.value;
              Utility::ForIndex<K + 1>([&](auto ii)
              {
                constexpr size_t i2 = ii.value;
                if constexpr (i2 + j2 <= K)
                {
                  constexpr size_t triRowStart =
                      j2 * (K + 1) - (j2 * (j2 - 1)) / 2;
                  constexpr size_t triIdx = triRowStart + i2;
                  constexpr size_t i = 0;
                  constexpr size_t j = j2;
                  constexpr size_t k = i2;
                  constexpr size_t tetraTotal =
                      (K + 1) * (K + 2) * (K + 3) / 6;
                  constexpr size_t m_tail = K - k;
                  constexpr size_t tetraTail =
                      (m_tail + 1) * (m_tail + 2) * (m_tail + 3) / 6;
                  constexpr size_t offset_k = tetraTotal - tetraTail;
                  constexpr size_t offset_j =
                      j * (K - k + 1) - (j * (j - 1)) / 2;
                  constexpr size_t tetIdx = offset_k + offset_j + i;
                  codomain[tetIdx] = domain[triIdx];
                }
              });
            });
          }
          else if (local == 2)
          {
            // Face (0,1,3) opposite vertex 2: tri(0,1,2) -> tet(0,1,3)
            // λ0 = λ_tri0, λ1 = λ_tri1, λ3 = λ_tri2, λ2 = 0
            // (i,j,k) = (i2, 0, j2)
            Utility::ForIndex<K + 1>([&](auto jj)
            {
              constexpr size_t j2 = jj.value;
              Utility::ForIndex<K + 1>([&](auto ii)
              {
                constexpr size_t i2 = ii.value;
                if constexpr (i2 + j2 <= K)
                {
                  constexpr size_t triRowStart =
                      j2 * (K + 1) - (j2 * (j2 - 1)) / 2;
                  constexpr size_t triIdx = triRowStart + i2;

                  constexpr size_t i = i2;
                  constexpr size_t j = 0;
                  constexpr size_t k = j2;

                  constexpr size_t tetraTotal =
                      (K + 1) * (K + 2) * (K + 3) / 6;
                  constexpr size_t m_tail = K - k;
                  constexpr size_t tetraTail =
                      (m_tail + 1) * (m_tail + 2) * (m_tail + 3) / 6;
                  constexpr size_t offset_k = tetraTotal - tetraTail;

                  constexpr size_t offset_j =
                      j * (K - k + 1) - (j * (j - 1)) / 2;

                  constexpr size_t tetIdx = offset_k + offset_j + i;

                  codomain[tetIdx] = domain[triIdx];
                }
              });
            });
          }
          else if (local == 3)
          {
            // Face (0,2,1) opposite vertex 3: tri(0,1,2) -> tet(0,2,1)
            // λ0 = λ_tri0, λ2 = λ_tri1, λ1 = λ_tri2, λ3 = 0
            // (i,j,k) = (j2, i2, 0)
            Utility::ForIndex<K + 1>([&](auto jj)
            {
              constexpr size_t j2 = jj.value;
              Utility::ForIndex<K + 1>([&](auto ii)
              {
                constexpr size_t i2 = ii.value;
                if constexpr (i2 + j2 <= K)
                {
                  constexpr size_t triRowStart =
                      j2 * (K + 1) - (j2 * (j2 - 1)) / 2;
                  constexpr size_t triIdx = triRowStart + i2;
                  constexpr size_t i = j2;
                  constexpr size_t j = i2;
                  constexpr size_t k = 0;
                  constexpr size_t tetraTotal =
                      (K + 1) * (K + 2) * (K + 3) / 6;
                  constexpr size_t m_tail = K - k;
                  constexpr size_t tetraTail =
                      (m_tail + 1) * (m_tail + 2) * (m_tail + 3) / 6;
                  constexpr size_t offset_k = tetraTotal - tetraTail;
                  constexpr size_t offset_j =
                      j * (K - k + 1) - (j * (j - 1)) / 2;
                  constexpr size_t tetIdx = offset_k + offset_j + i;
                  codomain[tetIdx] = domain[triIdx];
                }
              });
            });
          }
        }
        // -------------------------------------------------------------------
        // G = Wedge: from Triangle or Quad face -> Wedge
        //
        // Wedge vertices: (0,1,2) bottom tri, (3,4,5) top tri
        //
        // Faces (dim=2, Connectivity::getSubPolytopes):
        //   Local 0 : Triangle (0,1,2)      bottom
        //   Local 1 : Quad     (0,1,4,3)
        //   Local 2 : Quad     (1,2,5,4)
        //   Local 3 : Quad     (2,0,3,5)
        //   Local 4 : Triangle (3,5,4)      top
        //
        // Wedge DOF ordering:
        //   for k = 0..K:
        //     for triIdx = 0..TriangleCount-1:
        //       wedgeIdx = k*TriangleCount + triIdx
        //
        // - For Local 0 and 4: From = Triangle
        // - For Local 1,2,3   : From = Quadrilateral
        // -------------------------------------------------------------------
        else if constexpr (G == Type::Wedge)
        {
          // Triangular faces
          if (local == 0 || local == 4)
          {
            // Triangle -> Wedge
            assert(local == 0 || local == 4 && "Triangle -> Wedge only for Local = 0 (bottom) or 4 (top).");

            if (local == 0)
            {
              // bottom triangular face, k = 0
              Utility::ForIndex<TriangleCount>([&](auto ii)
              {
                constexpr size_t triIdx   = ii.value;
                constexpr size_t wedgeIdx = triIdx; // k=0
                codomain[wedgeIdx] = domain[triIdx];
              });
            }
            else if (local == 4)
            {
              // top triangular face, k = K
              Utility::ForIndex<TriangleCount>([&](auto ii)
              {
                constexpr size_t triIdx   = ii.value;
                constexpr size_t wedgeIdx = K * TriangleCount + triIdx;
                codomain[wedgeIdx] = domain[triIdx];
              });
            }
          }
          // Quadrilateral faces
          else if (local == 1 || local == 2 || local == 3)
          {
            // Quad -> Wedge
            assert(local >= 1 && local <= 3 && "Quadrilateral -> Wedge implemented only for quad faces Local = 1, 2, 3.");

            if (local == 1)
            {
              // Face (0,1,4,3): extrude edge 0->1 in z
              // base i along 0->1, j vertical 0..K
              Utility::ForIndex<K + 1>([&](auto jj)
              {
                constexpr size_t j = jj.value; // layer 0..K
                Utility::ForIndex<K + 1>([&](auto ii)
                {
                  constexpr size_t i = ii.value; // 0..K along edge 0->1
                  constexpr size_t quadIdx = j * (K + 1) + i;
                  // triangle edge 0->1: same as Segment->Triangle Local=0
                  constexpr size_t triEdgeIdx = i;
                  constexpr size_t wedgeIdx = j * TriangleCount + triEdgeIdx;
                  codomain[wedgeIdx] = domain[quadIdx];
                });
              });
            }
            else if (local == 2)
            {
              // Face (1,2,5,4): extrude edge 1->2 in z
              Utility::ForIndex<K + 1>([&](auto jj)
              {
                constexpr size_t j = jj.value;
                Utility::ForIndex<K + 1>([&](auto ii)
                {
                  constexpr size_t i = ii.value; // 0..K along edge 1->2
                  constexpr size_t quadIdx = j * (K + 1) + i;
                  // triangle edge 1->2 (Segment->Triangle Local=1)
                  constexpr size_t r = i;
                  constexpr size_t rowStart =
                      r * (K + 1) - (r * (r - 1)) / 2;
                  constexpr size_t triEdgeIdx = rowStart + (K - r);
                  constexpr size_t wedgeIdx = j * TriangleCount + triEdgeIdx;
                  codomain[wedgeIdx] = domain[quadIdx];
                });
              });
            }
            else if (local == 3)
            {
              // Face (2,0,3,5): extrude edge 2->0 in z
              Utility::ForIndex<K + 1>([&](auto jj)
              {
                constexpr size_t j = jj.value;
                Utility::ForIndex<K + 1>([&](auto ii)
                {
                  constexpr size_t i = ii.value; // 0..K along edge 2->0
                  constexpr size_t quadIdx = j * (K + 1) + i;
                  // triangle edge 2->0 (Segment->Triangle Local=2)
                  constexpr size_t r      = i;
                  constexpr size_t j_edge = K - r;
                  constexpr size_t rowStart =
                      j_edge * (K + 1) - (j_edge * (j_edge - 1)) / 2;
                  constexpr size_t triEdgeIdx = rowStart; // i=0 on that row
                  constexpr size_t wedgeIdx = j * TriangleCount + triEdgeIdx;
                  codomain[wedgeIdx] = domain[quadIdx];
                });
              });
            }
          }
        }
        else
        {
          Alert::Exception()
            << "Cochain for geometry " << G
            << " and local = " << local
            << " not implemented."
            << Alert::Raise;
        }
      }
  };

  template <size_t K, class Scalar>
  void H1<K, Scalar, Geometry::Mesh<Context::Local>>::getClosure(size_t d, Index idx)
  {
    if (m_visited[d][idx])
      return;

    const auto& mesh = m_mesh.get();
    const auto& conn = mesh.getConnectivity();
    const auto g = mesh.getGeometry(d, idx);

    m_visited[d][idx] = 1;

    auto& local = m_closure[d][idx];

    switch (g)
    {
      case Geometry::Polytope::Type::Point:
      {
        // exactly 1 local DOF
        local[0] = m_size++;
        break;
      }

      case Geometry::Polytope::Type::Segment:
      {
        const auto& inc = conn.getIncidence({ d, d - 1 }, idx); // 2 vertices

        // Left endpoint
        this->getClosure(d - 1, inc[0]);
        Cochain<Geometry::Polytope::Type::Segment>::map(0, local, m_closure[d - 1][inc[0]]);

        // Interior DOFs
        for (size_t k = 1; k + 1 < Cochain<Geometry::Polytope::Type::Segment>::Count; ++k)
          local[k] = m_size++;

        // Right endpoint
        this->getClosure(d - 1, inc[1]);
        Cochain<Geometry::Polytope::Type::Segment>::map(1, local, m_closure[d - 1][inc[1]]);

        break;
      }

      case Geometry::Polytope::Type::Triangle:
      {
        // Triangle edges: (0->1), (1->2), (2->0)
        const auto& inc = conn.getIncidence({ d, d - 1 }, idx);
        assert(inc.size() == 3);

        // Mark which triangle nodes belong to edges (boundary DOFs)
        std::array<uint8_t,
                   Cochain<Geometry::Polytope::Type::Triangle>::Count> used{};
        used.fill(0);

        // Convenience aliases
        using TriCochain = Cochain<Geometry::Polytope::Type::Triangle>;
        using SegCochain = Cochain<Geometry::Polytope::Type::Segment>;
        constexpr size_t Ns = SegCochain::Count;

        // Edge 0: (0->1), "bottom edge"
        {
          const Index e = inc[0];
          this->getClosure(d - 1, e);
          const auto& edge = m_closure[d - 1][e];

          Utility::ForIndex<Ns>([&](auto ii)
          {
            constexpr size_t r   = ii.value; // 0..K
            constexpr size_t tId = r;        // bottom row contiguous
            local[tId] = edge[r];
            used[tId]  = 1;
          });
        }

        // Edge 1: (1->2), "hypotenuse"
        {
          const Index e = inc[1];
          this->getClosure(d - 1, e);
          const auto& edge = m_closure[d - 1][e];

          Utility::ForIndex<Ns>([&](auto ii)
          {
            constexpr size_t r = ii.value;   // 0..K
            constexpr size_t j = r;          // param along hypotenuse

            constexpr size_t rowStart_j =
                j * (K + 1) - (j * (j - 1)) / 2;
            constexpr size_t tId = rowStart_j + (K - j);

            local[tId] = edge[r];
            used[tId]  = 1;
          });
        }

        // Edge 2: (2->0), "left edge"
        {
          const Index e = inc[2];
          this->getClosure(d - 1, e);
          const auto& edge = m_closure[d - 1][e];

          Utility::ForIndex<Ns>([&](auto ii)
          {
            constexpr size_t r      = ii.value;  // 0..K along edge
            constexpr size_t j      = K - r;
            constexpr size_t rowStart_j =
                j * (K + 1) - (j * (j - 1)) / 2;
            constexpr size_t tId = rowStart_j;   // i = 0 on row j

            local[tId] = edge[r];
            used[tId]  = 1;
          });
        }

        // Interior triangle DOFs (not on any edge)
        for (size_t tId = 0; tId < TriCochain::Count; ++tId)
        {
          if (!used[tId])
            local[tId] = m_size++;
        }

        break;
      }

      case Geometry::Polytope::Type::Quadrilateral:
      {
        // Quad edges in CCW loop: (0->1) bottom, (1->2) right,
        //                         (2->3) top,    (3->0) left
        const auto& inc = conn.getIncidence({ d, d - 1 }, idx);
        assert(inc.size() == 4);

        using QuadCochain = Cochain<Geometry::Polytope::Type::Quadrilateral>;
        constexpr size_t N1  = K + 1;

        // Map each edge closure into the quad boundary DOFs
        // local-edge indices: 0 bottom, 1 right, 2 top, 3 left
        for (size_t le = 0; le < 4; ++le)
        {
          const Index e = inc[le];
          this->getClosure(d - 1, e);
          const auto& edge = m_closure[d - 1][e];

          QuadCochain::map(le, local, edge);
        }

        // Interior quad DOFs:
        // indices (i,j) with 0 < i < K and 0 < j < K
        for (size_t j = 1; j < K; ++j)
        {
          for (size_t i = 1; i < K; ++i)
          {
            const size_t qId = j * N1 + i;
            local[qId] = m_size++;
          }
        }

        break;
      }

      case Geometry::Polytope::Type::Tetrahedron:
      {
        // Faces (dim=2, Connectivity::getSubPolytopes):
        //   Local 0: (1,2,3)  // +[1,2,3]
        //   Local 1: (0,3,2)  // -[0,2,3]
        //   Local 2: (0,1,3)  // +[0,1,3]
        //   Local 3: (0,2,1)  // -[0,1,2]
        const auto& inc = conn.getIncidence({ d, d - 1 }, idx);
        assert(inc.size() == 4);

        using TetCochain = Cochain<Geometry::Polytope::Type::Tetrahedron>;

        std::array<uint8_t, TetCochain::Count> used{};
        used.fill(0);

        constexpr size_t tetraTotal =
            (K + 1) * (K + 2) * (K + 3) / 6;

        // Helper lambdas to apply each face mapping and mark "used"
        auto face0 = [&](const Index fIdx)
        {
          const auto& face = m_closure[d - 1][fIdx]; // triangle DOFs

          Utility::ForIndex<K + 1>([&](auto jj)
          {
            constexpr size_t j2 = jj.value;
            Utility::ForIndex<K + 1>([&](auto ii)
            {
              constexpr size_t i2 = ii.value;
              if constexpr (i2 + j2 <= K)
              {
                constexpr size_t triRowStart =
                    j2 * (K + 1) - (j2 * (j2 - 1)) / 2;
                constexpr size_t triIdx = triRowStart + i2;

                // (i,j,k) = (K - i2 - j2, i2, j2)
                constexpr size_t i = K - i2 - j2;
                constexpr size_t j = i2;
                constexpr size_t k = j2;

                constexpr size_t m_tail = K - k;
                constexpr size_t tetraTail =
                    (m_tail + 1) * (m_tail + 2) * (m_tail + 3) / 6;
                constexpr size_t offset_k = tetraTotal - tetraTail;

                constexpr size_t offset_j =
                    j * (K - k + 1) - (j * (j - 1)) / 2;

                constexpr size_t tetIdx = offset_k + offset_j + i;

                local[tetIdx] = face[triIdx];
                used[tetIdx]  = 1;
              }
            });
          });
        };

        auto face1 = [&](const Index fIdx)
        {
          const auto& face = m_closure[d - 1][fIdx];

          Utility::ForIndex<K + 1>([&](auto jj)
          {
            constexpr size_t j2 = jj.value;
            Utility::ForIndex<K + 1>([&](auto ii)
            {
              constexpr size_t i2 = ii.value;
              if constexpr (i2 + j2 <= K)
              {
                constexpr size_t triRowStart =
                    j2 * (K + 1) - (j2 * (j2 - 1)) / 2;
                constexpr size_t triIdx = triRowStart + i2;

                // (i,j,k) = (0, j2, i2)
                constexpr size_t i = 0;
                constexpr size_t j = j2;
                constexpr size_t k = i2;

                constexpr size_t m_tail = K - k;
                constexpr size_t tetraTail =
                    (m_tail + 1) * (m_tail + 2) * (m_tail + 3) / 6;
                constexpr size_t offset_k = tetraTotal - tetraTail;

                constexpr size_t offset_j =
                    j * (K - k + 1) - (j * (j - 1)) / 2;

                constexpr size_t tetIdx = offset_k + offset_j + i;

                local[tetIdx] = face[triIdx];
                used[tetIdx]  = 1;
              }
            });
          });
        };

        auto face2 = [&](const Index fIdx)
        {
          const auto& face = m_closure[d - 1][fIdx];

          Utility::ForIndex<K + 1>([&](auto jj)
          {
            constexpr size_t j2 = jj.value;
            Utility::ForIndex<K + 1>([&](auto ii)
            {
              constexpr size_t i2 = ii.value;
              if constexpr (i2 + j2 <= K)
              {
                constexpr size_t triRowStart =
                    j2 * (K + 1) - (j2 * (j2 - 1)) / 2;
                constexpr size_t triIdx = triRowStart + i2;

                // (i,j,k) = (i2, 0, j2)
                constexpr size_t i = i2;
                constexpr size_t j = 0;
                constexpr size_t k = j2;

                constexpr size_t m_tail = K - k;
                constexpr size_t tetraTail =
                    (m_tail + 1) * (m_tail + 2) * (m_tail + 3) / 6;
                constexpr size_t offset_k = tetraTotal - tetraTail;

                constexpr size_t offset_j =
                    j * (K - k + 1) - (j * (j - 1)) / 2;

                constexpr size_t tetIdx = offset_k + offset_j + i;

                local[tetIdx] = face[triIdx];
                used[tetIdx]  = 1;
              }
            });
          });
        };

        auto face3 = [&](const Index fIdx)
        {
          const auto& face = m_closure[d - 1][fIdx];

          Utility::ForIndex<K + 1>([&](auto jj)
          {
            constexpr size_t j2 = jj.value;
            Utility::ForIndex<K + 1>([&](auto ii)
            {
              constexpr size_t i2 = ii.value;
              if constexpr (i2 + j2 <= K)
              {
                constexpr size_t triRowStart =
                    j2 * (K + 1) - (j2 * (j2 - 1)) / 2;
                constexpr size_t triIdx = triRowStart + i2;

                // (i,j,k) = (j2, i2, 0)
                constexpr size_t i = j2;
                constexpr size_t j = i2;
                constexpr size_t k = 0;

                constexpr size_t m_tail = K - k;
                constexpr size_t tetraTail =
                    (m_tail + 1) * (m_tail + 2) * (m_tail + 3) / 6;
                constexpr size_t offset_k = tetraTotal - tetraTail;

                constexpr size_t offset_j =
                    j * (K - k + 1) - (j * (j - 1)) / 2;

                constexpr size_t tetIdx = offset_k + offset_j + i;

                local[tetIdx] = face[triIdx];
                used[tetIdx]  = 1;
              }
            });
          });
        };

        // Recurse on faces and apply mappings
        // Face 0: (1,2,3)
        this->getClosure(d - 1, inc[0]);
        face0(inc[0]);

        // Face 1: (0,3,2)
        this->getClosure(d - 1, inc[1]);
        face1(inc[1]);

        // Face 2: (0,1,3)
        this->getClosure(d - 1, inc[2]);
        face2(inc[2]);

        // Face 3: (0,2,1)
        this->getClosure(d - 1, inc[3]);
        face3(inc[3]);

        // Interior tetra DOFs (not on any face)
        for (size_t tId = 0; tId < TetCochain::Count; ++tId)
        {
          if (!used[tId])
            local[tId] = m_size++;
        }

        break;
      }

      case Geometry::Polytope::Type::Wedge:
      {
        // Faces (dim=2, getSubPolytopes):
        //   Local 0 : Triangle (0,1,2)      bottom
        //   Local 1 : Quad     (0,1,4,3)
        //   Local 2 : Quad     (1,2,5,4)
        //   Local 3 : Quad     (2,0,3,5)
        //   Local 4 : Triangle (3,5,4)      top
        const auto& inc = conn.getIncidence({ d, d - 1 }, idx);
        assert(inc.size() == 5);

        using WedgeCochain = Cochain<Geometry::Polytope::Type::Wedge>;
        using TriCochain   = Cochain<Geometry::Polytope::Type::Triangle>;
        constexpr size_t TriCount = TriCochain::Count;

        std::array<uint8_t, WedgeCochain::Count> used{};
        used.fill(0);

        // Triangular faces ----------------------------------------------------
        // Local 0: bottom triangle (k = 0)
        {
          const Index f = inc[0];
          this->getClosure(d - 1, f);
          const auto& tri = m_closure[d - 1][f];

          Utility::ForIndex<TriCount>([&](auto ii)
          {
            constexpr size_t triIdx   = ii.value;
            constexpr size_t wedgeIdx = triIdx; // k = 0
            local[wedgeIdx] = tri[triIdx];
            used[wedgeIdx]  = 1;
          });
        }

        // Local 4: top triangle (k = K)
        {
          const Index f = inc[4];
          this->getClosure(d - 1, f);
          const auto& tri = m_closure[d - 1][f];

          Utility::ForIndex<TriCount>([&](auto ii)
          {
            constexpr size_t triIdx   = ii.value;
            constexpr size_t wedgeIdx = K * TriCount + triIdx;
            local[wedgeIdx] = tri[triIdx];
            used[wedgeIdx]  = 1;
          });
        }

        // Quadrilateral faces -------------------------------------------------
        // We follow the same formulas as Cochain<Wedge>::map, but we also
        // mark which wedge DOFs are used.

        // Local 1: Quad (0,1,4,3): extrude edge 0->1
        {
          const Index f = inc[1];
          this->getClosure(d - 1, f);
          const auto& quad = m_closure[d - 1][f];

          Utility::ForIndex<K + 1>([&](auto jj)
          {
            constexpr size_t j = jj.value; // layer 0..K
            Utility::ForIndex<K + 1>([&](auto ii)
            {
              constexpr size_t i = ii.value; // along edge 0->1

              constexpr size_t quadIdx = j * (K + 1) + i;
              // triangle edge 0->1 corresponds to triIdx = i
              constexpr size_t triEdgeIdx = i;
              constexpr size_t wedgeIdx   = j * TriCount + triEdgeIdx;

              local[wedgeIdx] = quad[quadIdx];
              used[wedgeIdx]  = 1;
            });
          });
        }

        // Local 2: Quad (1,2,5,4): extrude edge 1->2
        {
          const Index f = inc[2];
          this->getClosure(d - 1, f);
          const auto& quad = m_closure[d - 1][f];

          Utility::ForIndex<K + 1>([&](auto jj)
          {
            constexpr size_t j = jj.value;
            Utility::ForIndex<K + 1>([&](auto ii)
            {
              constexpr size_t i = ii.value; // along edge 1->2

              constexpr size_t quadIdx = j * (K + 1) + i;

              // triangle edge 1->2 (Segment->Triangle Local=1):
              //  r = i
              //  rowStart = r*(K+1) - r*(r-1)/2
              //  triEdgeIdx = rowStart + (K - r)
              constexpr size_t r        = i;
              constexpr size_t rowStart =
                  r * (K + 1) - (r * (r - 1)) / 2;
              constexpr size_t triEdgeIdx = rowStart + (K - r);
              constexpr size_t wedgeIdx   = j * TriCount + triEdgeIdx;

              local[wedgeIdx] = quad[quadIdx];
              used[wedgeIdx]  = 1;
            });
          });
        }

        // Local 3: Quad (2,0,3,5): extrude edge 2->0
        {
          const Index f = inc[3];
          this->getClosure(d - 1, f);
          const auto& quad = m_closure[d - 1][f];

          Utility::ForIndex<K + 1>([&](auto jj)
          {
            constexpr size_t j = jj.value;
            Utility::ForIndex<K + 1>([&](auto ii)
            {
              constexpr size_t i = ii.value; // along edge 2->0

              constexpr size_t quadIdx = j * (K + 1) + i;

              // triangle edge 2->0 (Segment->Triangle Local=2):
              //  r      = i
              //  j_edge = K - r
              //  rowStart = j_edge*(K+1) - j_edge*(j_edge-1)/2
              //  triEdgeIdx = rowStart (i=0)
              constexpr size_t r      = i;
              constexpr size_t j_edge = K - r;
              constexpr size_t rowStart =
                  j_edge * (K + 1) - (j_edge * (j_edge - 1)) / 2;
              constexpr size_t triEdgeIdx = rowStart;
              constexpr size_t wedgeIdx   = j * TriCount + triEdgeIdx;

              local[wedgeIdx] = quad[quadIdx];
              used[wedgeIdx]  = 1;
            });
          });
        }

        // Interior wedge DOFs: not touched by any face mapping
        for (size_t wId = 0; wId < WedgeCochain::Count; ++wId)
        {
          if (!used[wId])
            local[wId] = m_size++;
        }

        break;
      }
    }
  }


  template <size_t K, class Scalar>
  H1<K, Scalar, Geometry::Mesh<Context::Local>>::H1(
      std::integral_constant<size_t, K>, const MeshType& mesh)
    : m_mesh(mesh),
      m_size(0)
  {
    const size_t D = mesh.getDimension();

    m_visited.resize(D + 1);
    m_closure.resize(D + 1);

    // Pre-size closure arrays by geometry
    for (size_t d = 0; d <= D; ++d)
    {
      const size_t count = mesh.getPolytopeCount(d);
      m_closure[d].resize(count);

      for (Index i = 0; i < static_cast<Index>(count); ++i)
      {
        const auto g = mesh.getGeometry(d, i);
        switch (g)
        {
          case Geometry::Polytope::Type::Point:
          {
            m_closure[d][i].resize(
              Cochain<Geometry::Polytope::Type::Point>::Count);
            break;
          }
          case Geometry::Polytope::Type::Segment:
          {
            m_closure[d][i].resize(
              Cochain<Geometry::Polytope::Type::Segment>::Count);
            break;
          }
          case Geometry::Polytope::Type::Triangle:
          {
            m_closure[d][i].resize(
              Cochain<Geometry::Polytope::Type::Triangle>::Count);
            break;
          }
          case Geometry::Polytope::Type::Quadrilateral:
          {
            m_closure[d][i].resize(
              Cochain<Geometry::Polytope::Type::Quadrilateral>::Count);
            break;
          }
          case Geometry::Polytope::Type::Tetrahedron:
          {
            m_closure[d][i].resize(
              Cochain<Geometry::Polytope::Type::Tetrahedron>::Count);
            break;
          }
          case Geometry::Polytope::Type::Wedge:
          {
            m_closure[d][i].resize(
              Cochain<Geometry::Polytope::Type::Wedge>::Count);
            break;
          }
        }
      }

      // Now initialize visited flags for this dimension
      m_visited[d].assign(count, 0);
    }

    // Ensure we have incidence d -> d-1 for d >= 1
    for (size_t d = 1; d <= D; ++d)
      RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, d, d - 1);

    // Build closure starting from top-dimensional cells
    const size_t nCells = mesh.getPolytopeCount(D);
    for (Index c = 0; c < static_cast<Index>(nCells); ++c)
      this->getClosure(D, c);
  }

  template <size_t K, class Scalar>
  H1<K, Math::Vector<Scalar>, Geometry::Mesh<Context::Local>>::
  H1(std::integral_constant<size_t, K>, const Geometry::Mesh<ContextType>& mesh,
     size_t vdim)
    : m_mesh(mesh),
      m_vdim(vdim),
      m_size(0)
  {
    using MeshType    = Geometry::Mesh<ContextType>;
    using ScalarSpace = H1<K, Scalar, MeshType>;

    const size_t D = mesh.getDimension();

    // 1. Build the scalar H1 space on the same mesh
    ScalarSpace scalar(std::integral_constant<size_t, K>{}, mesh);

    const size_t scalarSize = scalar.getSize(); // total scalar DOFs
    m_size                  = scalarSize * vdim; // total vector DOFs

    // 2. Lift scalar closure to vector closure
    m_closure.resize(D + 1);

    for (size_t d = 0; d <= D; ++d)
    {
      const size_t count = mesh.getPolytopeCount(d);
      m_closure[d].resize(count);

      for (Index i = 0; i < static_cast<Index>(count); ++i)
      {
        const IndexArray& scalarLocal = scalar.getDOFs(d, i);
        const size_t nLocalScalar     = scalarLocal.size();

        IndexArray& vecLocal = m_closure[d][i];
        vecLocal.resize(nLocalScalar * vdim);

        // Local layout: (node q, component c) -> q*vdim + c
        // Global layout: block by component
        //   component c lives in [c*scalarSize, (c+1)*scalarSize)
        for (size_t q = 0; q < nLocalScalar; ++q)
        {
          const Index sIdx = scalarLocal(q); // scalar global index

          for (size_t c = 0; c < vdim; ++c)
          {
            const Index vIdx = sIdx + static_cast<Index>(c * scalarSize);
            vecLocal(q * vdim + c) = vIdx;
          }
        }
      }
    }
  }
}

#endif // RODIN_VARIATIONAL_H1_H1_HPP
