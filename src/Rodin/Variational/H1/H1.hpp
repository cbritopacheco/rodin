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
        using SegCochain = Cochain<Geometry::Polytope::Type::Segment>;
        constexpr size_t Ns = SegCochain::Count; // = K+1 GLL nodes

        // Segment vertices as stored in the connectivity (dimension 0)
        assert(d == 1);
        assert(d - 1 == 0);
        const auto& segVertsIA = conn.getIncidence({ d, d - 1 }, idx);
        assert(segVertsIA.size() == 2);

        const Index gv0 = segVertsIA[0]; // endpoint "0" in this cell
        const Index gv1 = segVertsIA[1]; // endpoint "1" in this cell

        // Ensure closure for the vertices is built
        this->getClosure(0, gv0);
        this->getClosure(0, gv1);

        // We adopt the canonical reference orientation for the segment:
        //   vertex 0 -> vertex 1
        //
        // GLL<K>::getNodes() is assumed ordered monotonically along
        // this 1D parametric direction, so:
        //   - local[0]      is the DOF at gv0,
        //   - local[Ns - 1] is the DOF at gv1,
        //   - intermediate indices are interior edge DOFs.

        // Left endpoint DOF (canonical vertex 0)
        local[0] = m_closure[0][gv0][0];

        // Interior segment DOFs
        for (size_t k = 1; k + 1 < Ns; ++k)
          local[k] = m_size++;

        // Right endpoint DOF (canonical vertex 1)
        local[Ns - 1] = m_closure[0][gv1][0];

        break;
      }

      case Geometry::Polytope::Type::Triangle:
      {
        // Incident edges (dimension 1)
        const auto& inc = conn.getIncidence({ d, d - 1 }, idx);
        assert(inc.size() == 3);

        using TriCochain = Cochain<Geometry::Polytope::Type::Triangle>;
        using SegCochain = Cochain<Geometry::Polytope::Type::Segment>;
        constexpr size_t Ns = SegCochain::Count; // K+1

        // Mark which triangle nodes belong to edges
        std::array<uint8_t, TriCochain::Count> used{};
        used.fill(0);

        // Triangle vertices in this cell's local order (0,1,2)
        const auto& triVertsIA = conn.getPolytope(d, idx);
        assert(triVertsIA.size() == 3);
        std::array<Index,3> v = {
          triVertsIA(0),
          triVertsIA(1),
          triVertsIA(2)
        };

        // Canonical local edge vertices in global indices:
        //  edge 0: (0->1)
        //  edge 1: (1->2)
        //  edge 2: (2->0)
        auto canonicalEdgeVerts = [&](size_t le) -> std::array<Index,2>
        {
          switch (le)
          {
            case 0: return { v[0], v[1] }; // (0->1)
            case 1: return { v[1], v[2] }; // (1->2)
            case 2: return { v[2], v[0] }; // (2->0)
            default:
              assert(false && "Invalid triangle local edge index.");
              return { 0, 0 };
          }
        };

        // For a given local edge 'le', find which segment entity in 'inc'
        // corresponds to it, and whether it is oriented forward or backward
        // compared to the canonical orientation.
        auto findEdgeEntity = [&](size_t le, bool& forward) -> Index
        {
          const auto wanted  = canonicalEdgeVerts(le); // [a,b]
          const Index a = wanted[0];
          const Index b = wanted[1];

          for (Index e : inc)
          {
            const auto& eVertsIA = conn.getPolytope(d - 1, e);
            assert(eVertsIA.size() == 2);
            const Index ev0 = eVertsIA(0);
            const Index ev1 = eVertsIA(1);

            if (ev0 == a && ev1 == b)
            {
              forward = true;  // stored with canonical orientation
              return e;
            }
            if (ev0 == b && ev1 == a)
            {
              forward = false; // stored reversed
              return e;
            }
          }

          assert(false && "Could not match triangle local edge to incident segment entity.");
          return -1;
        };

        // --- Edge 0: canonical (0->1), "bottom edge" -----------------------
        {
          bool forward = true;
          const Index e = findEdgeEntity(0, forward);
          this->getClosure(d - 1, e);
          const auto& edge = m_closure[d - 1][e];

          Utility::ForIndex<Ns>([&](auto ii)
          {
            constexpr size_t r = ii.value;       // 0..K
            const size_t rr    = forward ? r : (Ns - 1 - r);

            // bottom row contiguous in Fekete ordering: idx = r
            constexpr size_t tId = r;

            local[tId] = edge[rr];
            used[tId]  = 1;
          });
        }

        // --- Edge 1: canonical (1->2), "hypotenuse" ------------------------
        {
          bool forward = true;
          const Index e = findEdgeEntity(1, forward);
          this->getClosure(d - 1, e);
          const auto& edge = m_closure[d - 1][e];

          Utility::ForIndex<Ns>([&](auto ii)
          {
            constexpr size_t r = ii.value;       // 0..K
            const size_t rr    = forward ? r : (Ns - 1 - r);

            constexpr size_t j = r;              // parameter along hypotenuse
            constexpr size_t rowStart_j =
                j * (K + 1) - (j * (j - 1)) / 2;
            constexpr size_t tId = rowStart_j + (K - j);

            local[tId] = edge[rr];
            used[tId]  = 1;
          });
        }

        // --- Edge 2: canonical (2->0), "left edge" -------------------------
        {
          bool forward = true;
          const Index e = findEdgeEntity(2, forward);
          this->getClosure(d - 1, e);
          const auto& edge = m_closure[d - 1][e];

          Utility::ForIndex<Ns>([&](auto ii)
          {
            constexpr size_t r = ii.value;       // 0..K
            const size_t rr    = forward ? r : (Ns - 1 - r);

            constexpr size_t j = K - r;
            constexpr size_t rowStart_j =
                j * (K + 1) - (j * (j - 1)) / 2;
            constexpr size_t tId = rowStart_j;   // i = 0 on row j

            local[tId] = edge[rr];
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
        // Incidence: cell (d) -> edges (d-1)
        const auto& inc = conn.getIncidence({ d, d - 1 }, idx);
        assert(inc.size() == 4);

        using QuadCochain = Cochain<Geometry::Polytope::Type::Quadrilateral>;
        using SegCochain  = Cochain<Geometry::Polytope::Type::Segment>;
        constexpr size_t Ns = SegCochain::Count;   // K+1
        constexpr size_t N1 = K + 1;

        // Cell vertices in its local order 0..3
        const auto& quadVertsIA = conn.getPolytope(d, idx);
        assert(quadVertsIA.size() == 4);
        std::array<Index,4> v = {
          quadVertsIA(0),
          quadVertsIA(1),
          quadVertsIA(2),
          quadVertsIA(3)
        };

        // Canonical local edges of the reference quad:
        //   0: (0->1) bottom
        //   1: (1->2) right
        //   2: (2->3) top
        //   3: (3->0) left
        auto canonicalEdgeVerts = [&](size_t le) -> std::array<Index,2>
        {
          switch (le)
          {
            case 0: return { v[0], v[1] }; // bottom
            case 1: return { v[1], v[2] }; // right
            case 2: return { v[2], v[3] }; // top
            case 3: return { v[3], v[0] }; // left
            default:
              assert(false && "Invalid local edge index for quadrilateral.");
              return { 0, 0 };
          }
        };

        // For a given local edge le, find which segment in 'inc' it corresponds
        // to, and whether that segment is oriented forward or backward w.r.t.
        // the canonical (0->1), (1->2), (2->3), (3->0).
        auto findEdgeEntity = [&](size_t le, bool& forward) -> Index
        {
          const auto canon = canonicalEdgeVerts(le);
          const Index a = canon[0];
          const Index b = canon[1];

          for (Index e : inc)
          {
            const auto& eVertsIA = conn.getPolytope(d - 1, e);
            assert(eVertsIA.size() == 2);
            const Index ev0 = eVertsIA(0);
            const Index ev1 = eVertsIA(1);

            if (ev0 == a && ev1 == b)
            {
              forward = true;  // stored in canonical orientation
              return e;
            }
            if (ev0 == b && ev1 == a)
            {
              forward = false; // stored in opposite orientation
              return e;
            }
          }

          assert(false && "Could not match quad local edge to incident segment.");
          return -1;
        };

        // Map all four edges, correcting orientation via a small local buffer.
        for (size_t le = 0; le < 4; ++le)
        {
          bool forward = true;
          const Index e = findEdgeEntity(le, forward);

          // Ensure the segment closure is built
          this->getClosure(d - 1, e);
          const auto& edge = m_closure[d - 1][e]; // size = Ns

          // Re-orient edge DOFs into canonical 0..K ordering along the local edge
          std::array<Index, Ns> orientedEdge{};
          for (size_t r = 0; r < Ns; ++r)
          {
            const size_t rr = forward ? r : (Ns - 1 - r);
            orientedEdge[r] = edge[rr];
          }

          // Inject into quad using the canonical edge index 'le'
          QuadCochain::map(le, local, orientedEdge);
        }

        // Interior quad DOFs: (i,j) with 0 < i < K, 0 < j < K
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
        const auto& mesh = m_mesh.get();
        const auto& conn = mesh.getConnectivity();

        // Faces incident to this tetrahedron (dimension 2)
        const auto& inc = conn.getIncidence({ d, d - 1 }, idx);
        assert(inc.size() == 4);

        using TetCochain = Cochain<Geometry::Polytope::Type::Tetrahedron>;
        using TriCochain = Cochain<Geometry::Polytope::Type::Triangle>;

        std::array<uint8_t, TetCochain::Count> used{};
        used.fill(0);

        constexpr size_t tetraTotal =
            (K + 1) * (K + 2) * (K + 3) / 6;

        // Cell vertices in local order 0,1,2,3
        const auto& cellVertsIA = conn.getPolytope(d, idx);
        assert(cellVertsIA.size() == 4);
        std::array<Index, 4> v = {
          cellVertsIA(0),
          cellVertsIA(1),
          cellVertsIA(2),
          cellVertsIA(3)
        };

        // Canonical faces in terms of the cell vertices:
        // match Connectivity::getSubPolytopes convention:
        //   0: (1,2,3)  // +[1,2,3]
        //   1: (0,3,2)  // -[0,2,3]
        //   2: (0,1,3)  // +[0,1,3]
        //   3: (0,2,1)  // -[0,1,2]
        auto canonicalFaceVerts = [&](size_t lf) -> std::array<Index, 3>
        {
          switch (lf)
          {
            case 0: return { v[1], v[2], v[3] };
            case 1: return { v[0], v[3], v[2] };
            case 2: return { v[0], v[1], v[3] };
            case 3: return { v[0], v[2], v[1] };
            default:
              assert(false && "Invalid local face index for tetrahedron.");
              return { 0, 0, 0 };
          }
        };

        // All 6 permutations of a triangle's 3 local vertices
        static constexpr int perms[6][3] =
        {
          {0,1,2}, {1,2,0}, {2,0,1},
          {0,2,1}, {2,1,0}, {1,0,2}
        };

        // For a given local tetra face lf (0..3), find the triangle entity
        // in 'inc' which matches it, and the permutation that maps
        // canonical face vertices -> triangle-local vertices.
        auto getFaceEntityAndPerm =
          [&](size_t lf, std::array<int, 3>& canonToTri) -> Index
        {
          const auto wanted = canonicalFaceVerts(lf); // global vertex ids (a,b,c)

          for (Index f : inc)
          {
            const auto& fVertsIA = conn.getPolytope(d - 1, f);
            assert(fVertsIA.size() == 3);
            std::array<Index, 3> tv = {
              fVertsIA(0),
              fVertsIA(1),
              fVertsIA(2)
            };

            // Try all permutations of the triangle's local vertices
            for (int pi = 0; pi < 6; ++pi)
            {
              const int a = perms[pi][0];
              const int b = perms[pi][1];
              const int c = perms[pi][2];

              if (tv[a] == wanted[0] &&
                  tv[b] == wanted[1] &&
                  tv[c] == wanted[2])
              {
                // canonical vertex 0 -> triangle local a
                // canonical vertex 1 -> triangle local b
                // canonical vertex 2 -> triangle local c
                canonToTri[0] = a;
                canonToTri[1] = b;
                canonToTri[2] = c;
                return f;
              }
            }
          }

          assert(false && "Could not match tetra face to incident triangle entity.");
          return -1;
        };

        // Build a face DOF array in *canonical* FeketeTriangle<K> ordering
        // for this tetra face: 'faceCanon' corresponds to the face whose
        // vertices are the canonical ones used above, with the standard
        // FeketeTriangle enumeration.
        auto buildCanonicalFace =
          [&](Index fIdx,
              const std::array<int, 3>& canonToTri,
              IndexArray& faceCanon)
        {
          const auto& faceLocal = m_closure[d - 1][fIdx]; // triangle DOFs, as built for that triangle
          faceCanon.resize(TriCochain::Count);

          // Inverse permutation: triangle local vertex q -> canonical vertex p
          std::array<int, 3> triToCanon;
          for (int p = 0; p < 3; ++p)
          {
            const int q = canonToTri[p];
            triToCanon[q] = p;
          }

          size_t canonIdx = 0;
          for (size_t j2 = 0; j2 <= K; ++j2)
          {
            for (size_t i2 = 0; i2 <= K - j2; ++i2, ++canonIdx)
            {
              // Canonical FeketeTriangle integer barycentric (a,b,c) w.r.t.
              // vertices (0,1,2) of the reference triangle:
              //   a = K - i2 - j2,  b = i2,  c = j2
              const size_t a = K - i2 - j2;
              const size_t b = i2;
              const size_t c = j2;
              const size_t abc[3] = { a, b, c };

              // Triangle-local barycentric integers (a_t,b_t,c_t)
              // for vertices 0,1,2 of that triangle:
              const size_t a_t = abc[ triToCanon[0] ];
              const size_t b_t = abc[ triToCanon[1] ];
              const size_t c_t = abc[ triToCanon[2] ];

              // Back to the triangle's Fekete enumeration:
              //  i_loc = b_t, j_loc = c_t
              const size_t i_loc = b_t;
              const size_t j_loc = c_t;

              const size_t rowStartLoc =
                  j_loc * (K + 1) - (j_loc * (j_loc - 1)) / 2;
              const size_t locIdx = rowStartLoc + i_loc;

              faceCanon[canonIdx] = faceLocal[locIdx];
            }
          }
        };

        // Now apply the *canonical* tetra face maps using faceCanon.
        auto applyFace0 = [&](const IndexArray& faceCanon)
        {
          // Face 0: (1,2,3), opposite vertex 0
          // integer embedding: (i,j,k) = (K - i2 - j2, i2, j2)
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

                constexpr size_t m_tail   = K - k;
                constexpr size_t tetraTail =
                    (m_tail + 1) * (m_tail + 2) * (m_tail + 3) / 6;
                constexpr size_t offset_k = tetraTotal - tetraTail;

                constexpr size_t offset_j =
                    j * (K - k + 1) - (j * (j - 1)) / 2;

                constexpr size_t tetIdx = offset_k + offset_j + i;

                local[tetIdx] = faceCanon[triIdx];
                used[tetIdx]  = 1;
              }
            });
          });
        };

        auto applyFace1 = [&](const IndexArray& faceCanon)
        {
          // Face 1: (0,3,2), opposite vertex 1
          // embedding: (i,j,k) = (0, j2, i2)
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

                constexpr size_t m_tail   = K - k;
                constexpr size_t tetraTail =
                    (m_tail + 1) * (m_tail + 2) * (m_tail + 3) / 6;
                constexpr size_t offset_k = tetraTotal - tetraTail;

                constexpr size_t offset_j =
                    j * (K - k + 1) - (j * (j - 1)) / 2;

                constexpr size_t tetIdx = offset_k + offset_j + i;

                local[tetIdx] = faceCanon[triIdx];
                used[tetIdx]  = 1;
              }
            });
          });
        };

        auto applyFace2 = [&](const IndexArray& faceCanon)
        {
          // Face 2: (0,1,3), opposite vertex 2
          // embedding: (i,j,k) = (i2, 0, j2)
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

                constexpr size_t m_tail   = K - k;
                constexpr size_t tetraTail =
                    (m_tail + 1) * (m_tail + 2) * (m_tail + 3) / 6;
                constexpr size_t offset_k = tetraTotal - tetraTail;

                constexpr size_t offset_j =
                    j * (K - k + 1) - (j * (j - 1)) / 2;

                constexpr size_t tetIdx = offset_k + offset_j + i;

                local[tetIdx] = faceCanon[triIdx];
                used[tetIdx]  = 1;
              }
            });
          });
        };

        auto applyFace3 = [&](const IndexArray& faceCanon)
        {
          // Face 3: (0,2,1), opposite vertex 3
          // embedding: (i,j,k) = (j2, i2, 0)
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

                constexpr size_t m_tail   = K - k;
                constexpr size_t tetraTail =
                    (m_tail + 1) * (m_tail + 2) * (m_tail + 3) / 6;
                constexpr size_t offset_k = tetraTotal - tetraTail;

                constexpr size_t offset_j =
                    j * (K - k + 1) - (j * (j - 1)) / 2;

                constexpr size_t tetIdx = offset_k + offset_j + i;

                local[tetIdx] = faceCanon[triIdx];
                used[tetIdx]  = 1;
              }
            });
          });
        };

        // ---- Build and apply each face in canonical orientation -----------

        // Face 0: (1,2,3)
        {
          std::array<int, 3> canonToTri{};
          const Index f = getFaceEntityAndPerm(0, canonToTri);
          this->getClosure(d - 1, f);

          IndexArray faceCanon;
          buildCanonicalFace(f, canonToTri, faceCanon);
          applyFace0(faceCanon);
        }

        // Face 1: (0,3,2)
        {
          std::array<int, 3> canonToTri{};
          const Index f = getFaceEntityAndPerm(1, canonToTri);
          this->getClosure(d - 1, f);

          IndexArray faceCanon;
          buildCanonicalFace(f, canonToTri, faceCanon);
          applyFace1(faceCanon);
        }

        // Face 2: (0,1,3)
        {
          std::array<int, 3> canonToTri{};
          const Index f = getFaceEntityAndPerm(2, canonToTri);
          this->getClosure(d - 1, f);

          IndexArray faceCanon;
          buildCanonicalFace(f, canonToTri, faceCanon);
          applyFace2(faceCanon);
        }

        // Face 3: (0,2,1)
        {
          std::array<int, 3> canonToTri{};
          const Index f = getFaceEntityAndPerm(3, canonToTri);
          this->getClosure(d - 1, f);

          IndexArray faceCanon;
          buildCanonicalFace(f, canonToTri, faceCanon);
          applyFace3(faceCanon);
        }

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
        const auto& mesh = m_mesh.get();
        const auto& conn = mesh.getConnectivity();

        // Faces (dim=2, Connectivity::getSubPolytopes):
        //   Local 0 : Triangle (0,1,2)      bottom
        //   Local 1 : Quad     (0,1,4,3)
        //   Local 2 : Quad     (1,2,5,4)
        //   Local 3 : Quad     (2,0,3,5)
        //   Local 4 : Triangle (3,4,5)  (canonical wedge top)
        //
        // The actual Triangle / Quad entities in 'inc' may have these
        // vertex sets but in a permuted order. We match them combinatorially.
        const auto& inc = conn.getIncidence({ d, d - 1 }, idx);
        assert(inc.size() == 5);

        using WedgeCochain = Cochain<Geometry::Polytope::Type::Wedge>;
        using TriCochain   = Cochain<Geometry::Polytope::Type::Triangle>;
        using SegCochain   = Cochain<Geometry::Polytope::Type::Segment>;

        constexpr size_t TriCount = TriCochain::Count;
        constexpr size_t Ns       = SegCochain::Count; // = K+1
        constexpr size_t N1       = K + 1;

        std::array<uint8_t, WedgeCochain::Count> used{};
        used.fill(0);

        // Cell vertices in local order 0..5
        const auto& cellVertsIA = conn.getPolytope(d, idx);
        assert(cellVertsIA.size() == 6);
        std::array<Index,6> v = {
          cellVertsIA(0),
          cellVertsIA(1),
          cellVertsIA(2),
          cellVertsIA(3),
          cellVertsIA(4),
          cellVertsIA(5)
        };

        auto sort3 = [](std::array<Index,3> a)
        {
          std::sort(a.begin(), a.end());
          return a;
        };

        auto sort4 = [](std::array<Index,4> a)
        {
          std::sort(a.begin(), a.end());
          return a;
        };

        // Canonical triangular faces in terms of the cell vertices:
        //   lf = 0 -> bottom: (0,1,2)
        //   lf = 4 -> top   : (3,4,5)
        auto canonicalTriFaceVerts = [&](size_t lf) -> std::array<Index,3>
        {
          switch (lf)
          {
            case 0: // bottom
              return { v[0], v[1], v[2] };
            case 4: // top (canonical orientation (3,4,5))
              return { v[3], v[4], v[5] };
            default:
              assert(false && "Invalid local triangular face index for wedge.");
              return { 0, 0, 0 };
          }
        };

        // Canonical quadrilateral faces:
        //   lf = 1 -> (0,1,4,3)
        //   lf = 2 -> (1,2,5,4)
        //   lf = 3 -> (2,0,3,5)
        auto canonicalQuadFaceVerts = [&](size_t lf) -> std::array<Index,4>
        {
          switch (lf)
          {
            case 1:
              return { v[0], v[1], v[4], v[3] };
            case 2:
              return { v[1], v[2], v[5], v[4] };
            case 3:
              return { v[2], v[0], v[3], v[5] };
            default:
              assert(false && "Invalid local quadrilateral face index for wedge.");
              return { 0, 0, 0, 0 };
          }
        };

        // All 6 permutations of 3 vertices
        static constexpr int perms3[6][3] =
        {
          {0,1,2}, {1,2,0}, {2,0,1},
          {0,2,1}, {2,1,0}, {1,0,2}
        };

        // All 24 permutations of 4 vertices
        static constexpr int perms4[24][4] =
        {
          {0,1,2,3}, {0,1,3,2}, {0,2,1,3}, {0,2,3,1},
          {0,3,1,2}, {0,3,2,1}, {1,0,2,3}, {1,0,3,2},
          {1,2,0,3}, {1,2,3,0}, {1,3,0,2}, {1,3,2,0},
          {2,0,1,3}, {2,0,3,1}, {2,1,0,3}, {2,1,3,0},
          {2,3,0,1}, {2,3,1,0}, {3,0,1,2}, {3,0,2,1},
          {3,1,0,2}, {3,1,2,0}, {3,2,0,1}, {3,2,1,0}
        };

        // For lf in {0,4}: find the triangle entity in 'inc' that matches
        // the canonical face, and the permutation mapping canonical
        // face vertices -> triangle local vertices.
        auto getTriFaceEntityAndPerm =
          [&](size_t lf, std::array<int,3>& canonToTri) -> Index
        {
          const auto wanted  = canonicalTriFaceVerts(lf);
          const auto wantedS = sort3(wanted);

          for (Index f : inc)
          {
            if (conn.getGeometry(d - 1, f) != Geometry::Polytope::Type::Triangle)
              continue;

            const auto& fVertsIA = conn.getPolytope(d - 1, f);
            assert(fVertsIA.size() == 3);
            std::array<Index,3> tv = {
              fVertsIA(0),
              fVertsIA(1),
              fVertsIA(2)
            };

            // Match up to permutation
            for (int pi = 0; pi < 6; ++pi)
            {
              const int a = perms3[pi][0];
              const int b = perms3[pi][1];
              const int c = perms3[pi][2];

              if (tv[a] == wanted[0] &&
                  tv[b] == wanted[1] &&
                  tv[c] == wanted[2])
              {
                // canonical vertex 0 -> triangle local a
                // canonical vertex 1 -> triangle local b
                // canonical vertex 2 -> triangle local c
                canonToTri[0] = a;
                canonToTri[1] = b;
                canonToTri[2] = c;
                return f;
              }
            }
          }

          assert(false && "Could not match wedge triangular face to incident triangle entity.");
          return -1;
        };

        // Build a face DOF array in *canonical* FeketeTriangle<K> ordering
        // for this triangular face, given the triangle entity index fIdx
        // and the permutation canonToTri (canonical -> tri-local).
        auto buildCanonicalFace =
          [&](Index fIdx,
              const std::array<int,3>& canonToTri,
              IndexArray& faceCanon)
        {
          const auto& faceLocal = m_closure[d - 1][fIdx]; // triangle DOFs as built
          faceCanon.resize(TriCochain::Count);

          // Inverse permutation: triangle local vertex q -> canonical vertex p
          std::array<int,3> triToCanon;
          for (int p = 0; p < 3; ++p)
          {
            const int q = canonToTri[p];
            triToCanon[q] = p;
          }

          size_t canonIdx = 0;
          for (size_t j2 = 0; j2 <= K; ++j2)
          {
            for (size_t i2 = 0; i2 <= K - j2; ++i2, ++canonIdx)
            {
              // Canonical integer barycentric (a,b,c) for FeketeTriangle
              // w.r.t. vertices (0,1,2) of the reference triangle:
              //   a = K - i2 - j2,  b = i2,  c = j2
              const size_t a = K - i2 - j2;
              const size_t b = i2;
              const size_t c = j2;
              const size_t abc[3] = { a, b, c };

              // Triangle-local barycentric integers for its vertices 0,1,2:
              const size_t a_t = abc[ triToCanon[0] ];
              const size_t b_t = abc[ triToCanon[1] ];
              const size_t c_t = abc[ triToCanon[2] ];

              // Back to triangle's Fekete enumeration:
              //   i_loc = b_t,  j_loc = c_t
              const size_t i_loc = b_t;
              const size_t j_loc = c_t;

              const size_t rowStartLoc =
                  j_loc * (K + 1) - (j_loc * (j_loc - 1)) / 2;
              const size_t locIdx = rowStartLoc + i_loc;

              faceCanon[canonIdx] = faceLocal[locIdx];
            }
          }
        };

        // -------------------------------------------------------------------
        // 1) Triangular faces: bottom (lf=0, k=0) and top (lf=4, k=K)
        // -------------------------------------------------------------------

        // Bottom triangular face (k = 0)
        {
          std::array<int,3> canonToTri{};
          const Index f = getTriFaceEntityAndPerm(0, canonToTri);
          this->getClosure(d - 1, f);

          IndexArray faceCanon;
          buildCanonicalFace(f, canonToTri, faceCanon);

          Utility::ForIndex<TriCount>([&](auto ii)
          {
            constexpr size_t triIdx   = ii.value;
            constexpr size_t wedgeIdx = triIdx; // k = 0
            local[wedgeIdx] = faceCanon[triIdx];
            used[wedgeIdx]  = 1;
          });
        }

        // Top triangular face (k = K)
        {
          std::array<int,3> canonToTri{};
          const Index f = getTriFaceEntityAndPerm(4, canonToTri);
          this->getClosure(d - 1, f);

          IndexArray faceCanon;
          buildCanonicalFace(f, canonToTri, faceCanon);

          Utility::ForIndex<TriCount>([&](auto ii)
          {
            constexpr size_t triIdx   = ii.value;
            constexpr size_t wedgeIdx = K * TriCount + triIdx;
            local[wedgeIdx] = faceCanon[triIdx];
            used[wedgeIdx]  = 1;
          });
        }

        // -------------------------------------------------------------------
        // 2) Quadrilateral faces: fill side faces without overriding
        //    DOFs already set by the triangular faces (especially vertices).
        // -------------------------------------------------------------------

        // Local 1: Quad (0,1,4,3): extrude edge 0->1 in z
        {
          const auto wanted = canonicalQuadFaceVerts(1);
          const auto wantedS = sort4(wanted);

          Index f = -1;
          for (Index cand : inc)
          {
            if (conn.getGeometry(d - 1, cand) != Geometry::Polytope::Type::Quadrilateral)
              continue;
            const auto& qVertsIA = conn.getPolytope(d - 1, cand);
            assert(qVertsIA.size() == 4);
            std::array<Index,4> qv = {
              qVertsIA(0),
              qVertsIA(1),
              qVertsIA(2),
              qVertsIA(3)
            };
            if (sort4(qv) == wantedS)
            {
              f = cand;
              break;
            }
          }

          this->getClosure(d - 1, f);
          const auto& quad = m_closure[d - 1][f];

          Utility::ForIndex<K + 1>([&](auto jj)
          {
            constexpr size_t j = jj.value; // layer 0..K
            Utility::ForIndex<K + 1>([&](auto ii)
            {
              constexpr size_t i = ii.value; // 0..K along edge 0->1

              constexpr size_t quadIdx = j * (K + 1) + i;
              // triangle edge 0->1 corresponds to triEdgeIdx = i
              constexpr size_t triEdgeIdx = i;
              constexpr size_t wedgeIdx   = j * TriCount + triEdgeIdx;

              if (!used[wedgeIdx])
              {
                local[wedgeIdx] = quad[quadIdx];
                used[wedgeIdx]  = 1;
              }
            });
          });
        }

        // Local 2: Quad (1,2,5,4): extrude edge 1->2
        {
          const auto wanted = canonicalQuadFaceVerts(2);
          const auto wantedS = sort4(wanted);

          Index f = -1;
          for (Index cand : inc)
          {
            if (conn.getGeometry(d - 1, cand) != Geometry::Polytope::Type::Quadrilateral)
              continue;
            const auto& qVertsIA = conn.getPolytope(d - 1, cand);
            assert(qVertsIA.size() == 4);
            std::array<Index,4> qv = {
              qVertsIA(0),
              qVertsIA(1),
              qVertsIA(2),
              qVertsIA(3)
            };
            if (sort4(qv) == wantedS)
            {
              f = cand;
              break;
            }
          }
          this->getClosure(d - 1, f);
          const auto& quad = m_closure[d - 1][f];

          Utility::ForIndex<K + 1>([&](auto jj)
          {
            constexpr size_t j = jj.value;
            Utility::ForIndex<K + 1>([&](auto ii)
            {
              constexpr size_t i = ii.value; // 0..K along edge 1->2

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

              if (!used[wedgeIdx])
              {
                local[wedgeIdx] = quad[quadIdx];
                used[wedgeIdx]  = 1;
              }
            });
          });
        }

        // Local 3: Quad (2,0,3,5): extrude edge 2->0
        {
          const auto wanted = canonicalQuadFaceVerts(3);
          const auto wantedS = sort4(wanted);

          Index f = -1;
          for (Index cand : inc)
          {
            if (conn.getGeometry(d - 1, cand) != Geometry::Polytope::Type::Quadrilateral)
              continue;
            const auto& qVertsIA = conn.getPolytope(d - 1, cand);
            assert(qVertsIA.size() == 4);
            std::array<Index,4> qv = {
              qVertsIA(0),
              qVertsIA(1),
              qVertsIA(2),
              qVertsIA(3)
            };
            if (sort4(qv) == wantedS)
            {
              f = cand;
              break;
            }
          }

          this->getClosure(d - 1, f);
          const auto& quad = m_closure[d - 1][f];

          Utility::ForIndex<K + 1>([&](auto jj)
          {
            constexpr size_t j = jj.value;
            Utility::ForIndex<K + 1>([&](auto ii)
            {
              constexpr size_t i = ii.value; // 0..K along edge 2->0

              constexpr size_t quadIdx = j * (K + 1) + i;

              // triangle edge 2->0 (Segment->Triangle Local=2):
              //  r      = i
              //  j_edge = K - r
              //  rowStart = j_edge*(K+1) - j_edge*(j_edge-1)/2
              //  triEdgeIdx = rowStart (i=0 on that row)
              constexpr size_t r      = i;
              constexpr size_t j_edge = K - r;
              constexpr size_t rowStart =
                  j_edge * (K + 1) - (j_edge * (j_edge - 1)) / 2;
              constexpr size_t triEdgeIdx = rowStart;
              constexpr size_t wedgeIdx   = j * TriCount + triEdgeIdx;

              if (!used[wedgeIdx])
              {
                local[wedgeIdx] = quad[quadIdx];
                used[wedgeIdx]  = 1;
              }
            });
          });
        }

        // -------------------------------------------------------------------
        // 3) Interior wedge DOFs (not on any face)
        // -------------------------------------------------------------------
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
