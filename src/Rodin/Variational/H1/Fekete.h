/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_FEKETE_H
#define RODIN_VARIATIONAL_H1_FEKETE_H

#include <array>
#include <cstddef>

#include "Rodin/Types.h"
#include "Rodin/Math/Vector.h"

#include "WarpBlend.h"

#define RODIN_VARIATIONAL_H1_FEKETE_TOLERANCE 1e-12

namespace Rodin::Variational
{
  /**
   * @brief Cached Fekete-type nodes on the reference triangle for degree K.
   *
   * The nodes are constructed as:
   *  - start from equispaced barycentric nodes on the reference triangle
   *    with vertices (0,0), (1,0), (0,1),
   *  - apply the warp–blend algorithm to move them towards Fekete positions.
   *
   * For a given polynomial degree K, the total number of nodes is
   * \f$ (K+1)(K+2)/2 \f$. The nodes are computed once per template
   * instantiation and cached in a static array.
   *
   * @tparam K Polynomial degree on the triangle.
   */
  template <size_t K>
  class FeketeTriangle
  {
  public:
    static constexpr size_t Count = (K + 1) * (K + 2) / 2;

    static const std::array<Math::SpatialPoint, Count>& getNodes()
    {
      static const std::array<Math::SpatialPoint, Count> s_nodes = compute();
      return s_nodes;
    }

  private:
    static std::array<Math::SpatialPoint, Count> compute()
    {
      std::array<Math::SpatialPoint, Count> nodes{};
      size_t idx = 0;

      if constexpr (K == 0)
      {
        nodes[0] = Math::SpatialPoint{{0.0, 0.0}};
        return nodes;
      }

      // 1. Equispaced nodes in (s,t) on ref triangle
      for (size_t j = 0; j <= K; ++j)
      {
        for (size_t i = 0; i <= K - j; ++i, ++idx)
        {
          const Real s = static_cast<Real>(i) / static_cast<Real>(K);
          const Real t = static_cast<Real>(j) / static_cast<Real>(K);
          nodes[idx] = Math::SpatialPoint{{s, t}};
        }
      }

      // 2. Warp–blend → Fekete-type positions (does not change indexing)
      WarpBlendTriangle<K>::template apply<Count>(nodes);

      // 3. Reorder into H1 canonical ordering
      constexpr Real tol = static_cast<Real>(RODIN_VARIATIONAL_H1_FEKETE_TOLERANCE);

      std::array<Math::SpatialPoint, Count> reordered{};
      std::array<bool, Count> used{};
      used.fill(false);

      // --- 3.1 vertices: indices v0 (L1≈1), v1 (L2≈1), v2 (L3≈1) ---
      size_t v0 = Count, v1 = Count, v2 = Count;
      for (size_t i = 0; i < Count; ++i)
      {
        const Real x = nodes[i].x();
        const Real y = nodes[i].y();
        const Real L1 = static_cast<Real>(1.0) - x - y;
        const Real L2 = x;
        const Real L3 = y;

        if (v0 == Count && Math::abs(L1 - static_cast<Real>(1.0)) < tol)
          v0 = i;
        else if (v1 == Count && Math::abs(L2 - static_cast<Real>(1.0)) < tol)
          v1 = i;
        else if (v2 == Count && Math::abs(L3 - static_cast<Real>(1.0)) < tol)
          v2 = i;
      }

      // Sanity
      assert(v0 < Count && v1 < Count && v2 < Count);

      reordered[0] = nodes[v0]; used[v0] = true; // vertex 0 = (0,0)
      reordered[1] = nodes[v1]; used[v1] = true; // vertex 1 = (1,0)
      reordered[2] = nodes[v2]; used[v2] = true; // vertex 2 = (0,1)

      // --- 3.2 classify remaining nodes into edges and interior ---
      std::array<size_t, Count> edge01{};
      std::array<size_t, Count> edge12{};
      std::array<size_t, Count> edge20{};
      std::array<size_t, Count> interior{};

      size_t n01 = 0, n12 = 0, n20 = 0, ni = 0;

      for (size_t i = 0; i < Count; ++i)
      {
        if (used[i]) continue;

        const Real x = nodes[i].x();
        const Real y = nodes[i].y();
        const Real L1 = static_cast<Real>(1.0) - x - y;
        const Real L2 = x;
        const Real L3 = y;

        // edge 0–1 : y = 0 → L3 ≈ 0
        if (Math::abs(L3) < tol)
          edge01[n01++] = i;
        // edge 1–2 : x + y = 1 → L1 ≈ 0
        else if (Math::abs(L1) < tol)
          edge12[n12++] = i;
        // edge 2–0 : x = 0 → L2 ≈ 0
        else if (Math::abs(L2) < tol)
          edge20[n20++] = i;
        else
          interior[ni++] = i;
      }

      // --- 3.3 sort each edge by parameter t consistent with Connectivity ---
      // edge (0->1): t = L2 increasing from 0 to 1
      auto sort_edge01 = [&](size_t* idxs, size_t n)
      {
        std::sort(idxs, idxs + n, [&](size_t a, size_t b)
        {
          return nodes[a].x() < nodes[b].x();
        });
      };

      // edge (1->2): t = L3 increasing from 0 to 1 along 1→2
      auto sort_edge12 = [&](size_t* idxs, size_t n)
      {
        std::sort(idxs, idxs + n, [&](size_t a, size_t b)
        {
          return nodes[a].y() < nodes[b].y();
        });
      };

      // edge (2->0): t = L1 increasing from 0 to 1 along 2→0
      auto sort_edge20 = [&](size_t* idxs, size_t n)
      {
        std::sort(idxs, idxs + n, [&](size_t a, size_t b)
        {
          const Real xa = nodes[a].x();
          const Real ya = nodes[a].y();
          const Real xb = nodes[b].x();
          const Real yb = nodes[b].y();
          const Real L1a = static_cast<Real>(1.0) - xa - ya;
          const Real L1b = static_cast<Real>(1.0) - xb - yb;
          return L1a < L1b;
        });
      };

      sort_edge01(edge01.data(), n01); // corresponds to (0->1)
      sort_edge12(edge12.data(), n12); // (1->2)
      sort_edge20(edge20.data(), n20); // (2->0)

      // --- 3.4 fill reordered ---
      size_t out = 3;
      for (size_t k = 0; k < n01; ++k) reordered[out++] = nodes[edge01[k]];
      for (size_t k = 0; k < n12; ++k) reordered[out++] = nodes[edge12[k]];
      for (size_t k = 0; k < n20; ++k) reordered[out++] = nodes[edge20[k]];
      for (size_t k = 0; k < ni;  ++k) reordered[out++] = nodes[interior[k]];

      assert(out == Count);
      return reordered;
    }
  };

  /**
   * @brief Cached Fekete-type nodes on the reference tetrahedron for degree K.
   *
   * Reference tetrahedron:
   *   (0,0,0), (1,0,0), (0,1,0), (0,0,1).
   *
   * Construction:
   *  - start from equispaced barycentric nodes (integer lattice with
   *    i + j + k <= K, normalized by K),
   *  - apply warp–blend algorithm (WarpBlendTetrahedron<K>) to move them
   *    towards Fekete-type positions.
   *
   * For a given polynomial degree K, the number of nodes is
   * @f$ (K + 1)(K + 2)(K + 3) / 6 @f$.
   *
   * Nodes are computed once per template instantiation and cached.
   *
   * @tparam K Polynomial degree on the tetrahedron.
   */
  template <size_t K>
  class FeketeTetrahedron
  {
    public:
      /// Total number of nodes: (K+1)(K+2)(K+3)/6.
      static constexpr size_t Count = (K + 1) * (K + 2) * (K + 3) / 6;

      /// Return cached nodes as a std::array.
      static const std::array<Math::SpatialPoint, Count>& getNodes()
      {
        static const std::array<Math::SpatialPoint, Count> s_nodes = compute();
        return s_nodes;
      }

    private:
      /// Build equispaced nodes and apply warp–blend once.
      static std::array<Math::SpatialPoint, Count> compute()
      {
        std::array<Math::SpatialPoint, Count> nodes;
        size_t idx = 0;

        if constexpr (K == 0)
        {
          // Single vertex node; choose (0,0,0) by convention.
          nodes[0] = Math::SpatialPoint{{0.0, 0.0, 0.0}};
          return nodes;
        }

        // Equispaced nodes on the reference tetrahedron:
        // vertices (0,0,0), (1,0,0), (0,1,0), (0,0,1)
        for (size_t k = 0; k <= K; ++k)
        {
          for (size_t j = 0; j <= K - k; ++j)
          {
            for (size_t i = 0; i <= K - j - k; ++i, ++idx)
            {
              const Real r = static_cast<Real>(i) / static_cast<Real>(K);
              const Real s = static_cast<Real>(j) / static_cast<Real>(K);
              const Real t = static_cast<Real>(k) / static_cast<Real>(K);
              nodes[idx] = Math::SpatialPoint{{r, s, t}};
            }
          }
        }

        // 2) WarpBlendTetrahedron<K>::apply(nodes);
        WarpBlendTetrahedron<K>::template apply<Count>(nodes);

        // 3) Reorder into H1 convention
        constexpr Real tol = static_cast<Real>(RODIN_VARIATIONAL_H1_FEKETE_TOLERANCE);

        std::array<Math::SpatialPoint, Count> reordered;
        std::array<bool, Count> used{};
        used.fill(false);

        // 3.1 vertices: indices v0..v3
        size_t v[4] = {Count, Count, Count, Count};
        for (size_t i = 0; i < Count; ++i)
        {
          const Real x = nodes[i].x();
          const Real y = nodes[i].y();
          const Real z = nodes[i].z();

          const Real L1 = Real(1) - x - y - z;
          const Real L2 = x;
          const Real L3 = y;
          const Real L4 = z;

          if (v[0]==Count && Math::abs(L1 - Real(1)) < tol) v[0] = i;
          else if (v[1]==Count && Math::abs(L2 - Real(1)) < tol) v[1] = i;
          else if (v[2]==Count && Math::abs(L3 - Real(1)) < tol) v[2] = i;
          else if (v[3]==Count && Math::abs(L4 - Real(1)) < tol) v[3] = i;
        }
        assert(v[0]<Count && v[1]<Count && v[2]<Count && v[3]<Count);

        size_t out = 0;
        for (int a = 0; a < 4; ++a)
        {
          reordered[out++] = nodes[v[a]];
          used[v[a]] = true;
        }

        // 3.2 edge, face, interior buckets
        std::vector<size_t> edge[6];
        std::vector<size_t> face[4];
        std::vector<size_t> vol;

        for (size_t i = 0; i < Count; ++i)
        {
          if (used[i]) continue;

          const Real x = nodes[i].x();
          const Real y = nodes[i].y();
          const Real z = nodes[i].z();

          const Real L1 = Real(1) - x - y - z;
          const Real L2 = x;
          const Real L3 = y;
          const Real L4 = z;

          const bool z1 = (Math::abs(L1) < tol);
          const bool z2 = (Math::abs(L2) < tol);
          const bool z3 = (Math::abs(L3) < tol);
          const bool z4 = (Math::abs(L4) < tol);

          const int nZero = (z1?1:0)+(z2?1:0)+(z3?1:0)+(z4?1:0);

          if (nZero == 2)
          {
            // edge: 2 Li ~ 0, 2 > 0
            // map non-zero pair {ia,ib} → edge index 0..5
            int nz[2], cnt=0;
            if (!z1) nz[cnt++] = 0;
            if (!z2) nz[cnt++] = 1;
            if (!z3) nz[cnt++] = 2;
            if (!z4) nz[cnt++] = 3;
            assert(cnt==2);

            int a = nz[0], b = nz[1];

            int e = -1;
            if ((a==0 && b==1) || (a==1 && b==0)) e = 0; // 0-1
            else if ((a==0 && b==2) || (a==2 && b==0)) e = 1; // 0-2
            else if ((a==0 && b==3) || (a==3 && b==0)) e = 2; // 0-3
            else if ((a==1 && b==2) || (a==2 && b==1)) e = 3; // 1-2
            else if ((a==1 && b==3) || (a==3 && b==1)) e = 4; // 1-3
            else if ((a==2 && b==3) || (a==3 && b==2)) e = 5; // 2-3

            assert(e>=0);
            edge[e].push_back(i);
          }
          else if (nZero == 1)
          {
            // face interior node: exactly one Li ≈ 0
            int f = -1;
            if (z1) f = 0; // face opposite vertex 0: (1,2,3)
            else if (z2) f = 1; // opposite 1
            else if (z3) f = 2; // opposite 2
            else if (z4) f = 3; // opposite 3
            assert(f>=0);
            face[f].push_back(i);
          }
          else
          {
            // interior (volume) node
            vol.push_back(i);
          }
        }

        // 3.3 sort edge nodes along each edge consistent with orientation
        // Connectivity edges: (0,1),(0,2),(0,3),(1,2),(1,3),(2,3)

        auto sort_by_t = [&](std::vector<size_t>& eNodes, auto get_t)
        {
          std::sort(eNodes.begin(), eNodes.end(),
                    [&](size_t a, size_t b)
                    {
                      return get_t(nodes[a]) < get_t(nodes[b]);
                    });
        };

        // parameter t in [0,1] along each edge:
        // 0->1: t = L2 (= x)
        sort_by_t(edge[0], [](const Math::SpatialPoint& p){ return p.x(); });
        // 0->2: t = L3 (= y)
        sort_by_t(edge[1], [](const Math::SpatialPoint& p){ return p.y(); });
        // 0->3: t = L4 (= z)
        sort_by_t(edge[2], [](const Math::SpatialPoint& p){ return p.z(); });
        // 1->2: z=0, L1=0, t = L3 (= y)
        sort_by_t(edge[3], [](const Math::SpatialPoint& p){ return p.y(); });
        // 1->3: y=0, L1=0, t = L4 (= z)
        sort_by_t(edge[4], [](const Math::SpatialPoint& p){ return p.z(); });
        // 2->3: x=0, L1=0, t = L4 (= z)
        sort_by_t(edge[5], [](const Math::SpatialPoint& p){ return p.z(); });

        // 3.4 fill reordered: vertices, edges, faces, volume
        for (int e = 0; e < 6; ++e)
          for (auto idx : edge[e])
            reordered[out++] = nodes[idx];

        for (int f = 0; f < 4; ++f)
          for (auto idx : face[f])
            reordered[out++] = nodes[idx];

        for (auto idx : vol)
          reordered[out++] = nodes[idx];

        assert(out == Count);
        return reordered;
      }
  };
}

#endif
