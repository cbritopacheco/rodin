/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_WARPBLEND_H
#define RODIN_VARIATIONAL_H1_WARPBLEND_H

#include <cstddef>
#include <array>
#include <algorithm>
#include <cmath>

#include "Rodin/Types.h"
#include "Rodin/Math/Common.h"
#include "Rodin/Math/SpatialVector.h"

#include "LagrangeBasis.h" // LagrangeBasis1D<K>
#include "GLL.h"           // GLL<K>

#define RODIN_VARIATIONAL_H1_WARPBLEND_TOLERANCE 1e-14

namespace Rodin::Variational
{
  /**
   * @brief Computes the 1D warp factor for moving equispaced nodes to GLL positions.
   *
   * Given an evaluation point @f$ r \in [-1,1] @f$, this class computes the
   * warp displacement that interpolates the difference between GLL nodes and
   * equispaced nodes using Lagrange interpolation on the equispaced nodes.
   *
   * The warp factor is normalized by @f$ (1 - r^2) @f$ to ensure smooth
   * blending near the endpoints.
   *
   * @tparam K Polynomial degree.
   *
   * @see GLL for the target GLL node positions.
   */
  template <size_t K>
  class WarpFactor1D
  {
    public:
      /**
       * @brief Computes the warp factor at point r.
       *
       * @param r Evaluation point in [-1,1].
       * @return The warp displacement factor.
       */
      static Real get(Real r)
      {
        if constexpr (K <= 1)
          return static_cast<Real>(0.0);

        // Equispaced nodes req_j on [-1,1]
        static const std::array<Real, K + 1> s_req = []()
        {
          std::array<Real, K + 1> req{};
          for (size_t j = 0; j <= K; ++j)
          {
            req[j] = static_cast<Real>(-1.0)
                   + static_cast<Real>(2.0) * static_cast<Real>(j)
                   / static_cast<Real>(K);
          }
          return req;
        }();

        // Legendre–GLL nodes on [-1,1]
        const auto& gll = GLL<K>::getNodes();

        // Interpolate (gll - req) at r using Lagrange basis on req
        Real warp = static_cast<Real>(0.0);
        for (size_t i = 0; i <= K; ++i)
        {
          const Real Li = LagrangeBasis1D<K>::getBasis(i, r, s_req);
          warp += Li * (gll[i] - s_req[i]);
        }

        if (Math::abs(r) < 1.0 - RODIN_VARIATIONAL_H1_WARPBLEND_TOLERANCE)
        {
          const Real sf = 1.0 - r * r;
          warp /= sf;
        }
        else
        {
          warp = static_cast<Real>(0.0);
        }

        return warp;
      }
  };

  /**
   * @brief Optimized blending parameter α for triangular warp-blend.
   *
   * Provides the blending parameter α for different polynomial degrees,
   * empirically optimized to minimize the Lebesgue constant.
   * Values follow Hesthaven & Warburton, "Nodal DG Methods" (2008).
   *
   * @tparam K Polynomial degree.
   */
  template <size_t K>
  class TriangleBlend
  {
    public:
      static constexpr Real getAlpha()
      {
        if constexpr (K <= 2)
          return 0.0;
        else if constexpr (K == 3)
          return 1.4152;
        else if constexpr (K == 4)
          return 0.1001;
        else if constexpr (K == 5)
          return 0.2751;
        else if constexpr (K == 6)
          return 0.9800;
        else if constexpr (K == 7)
          return 1.0999;
        else if constexpr (K == 8)
          return 1.2832;
        else if constexpr (K == 9)
          return 1.3648;
        else if constexpr (K == 10)
          return 1.4773;
        else if constexpr (K == 11)
          return 1.4959;
        else if constexpr (K == 12)
          return 1.5743;
        else if constexpr (K == 13)
          return 1.5770;
        else if constexpr (K == 14)
          return 1.6223;
        else if constexpr (K == 15)
          return 1.6258;
        else
          return 5.0 / 3.0; // HW choice for higher orders
      }
  };

  /**
   * @brief Optimized blending parameter α for tetrahedral warp-blend.
   *
   * Provides the blending parameter α for different polynomial degrees,
   * empirically optimized for tetrahedra.
   *
   * @tparam K Polynomial degree.
   */
  template <size_t K>
  class TetrahedronBlend
  {
    public:
      static constexpr Real getAlpha()
      {
        if constexpr (K <= 3)
          return 0.0;
        else if constexpr (K == 4)
          return 0.1002;
        else if constexpr (K == 5)
          return 1.1332;
        else if constexpr (K == 6)
          return 1.5608;
        else if constexpr (K == 7)
          return 1.3413;
        else if constexpr (K == 8)
          return 1.2577;
        else if constexpr (K == 9)
          return 1.1603;
        else if constexpr (K == 10)
          return 1.10153;
        else if constexpr (K == 11)
          return 0.6080;
        else if constexpr (K == 12)
          return 0.4523;
        else if constexpr (K == 13)
          return 0.8856;
        else if constexpr (K == 14)
          return 0.8717;
        else if constexpr (K == 15)
          return 0.9655;
        else
          return 1.0; // HW choice for higher orders
      }
  };

  /**
   * @brief Computes 2D warp-blend shift for a triangular face.
   *
   * Given barycentric coordinates @f$ (L_1, L_2, L_3) @f$ on an equilateral
   * triangle, computes the displacement @f$ (\Delta x, \Delta y) @f$ that
   * moves the point towards its Fekete position.
   *
   * @tparam K Polynomial degree.
   */
  template <size_t K>
  class WarpShiftFace2D
  {
    public:
      static std::array<Real, 2> apply(Real L1, Real L2, Real L3, Real alpha)
      {
        if constexpr (K <= 1)
          return { Real(0), Real(0) };

        // 2) blending per edge
        const Real blend1 = L2 * L3; // edge opposite L1
        const Real blend2 = L1 * L3; // edge opposite L2
        const Real blend3 = L1 * L2; // edge opposite L3

        // 3) edge coordinates in [-1,1]
        const Real r1 = L3 - L2;
        const Real r2 = L1 - L3;
        const Real r3 = L2 - L1;

        const Real warpf1 = static_cast<Real>(4.0) * WarpFactor1D<K>::get(r1);
        const Real warpf2 = static_cast<Real>(4.0) * WarpFactor1D<K>::get(r2);
        const Real warpf3 = static_cast<Real>(4.0) * WarpFactor1D<K>::get(r3);

        const Real aL1 = alpha * L1;
        const Real aL2 = alpha * L2;
        const Real aL3 = alpha * L3;

        const Real warp1 = blend1 * warpf1 * (static_cast<Real>(1.0) + aL1 * aL1);
        const Real warp2 = blend2 * warpf2 * (static_cast<Real>(1.0) + aL2 * aL2);
        const Real warp3 = blend3 * warpf3 * (static_cast<Real>(1.0) + aL3 * aL3);

        // 5) shift in equilateral triangle
        const Real cos2pi3 = static_cast<Real>(-0.5);
        const Real sin2pi3 = static_cast<Real>(0.5) * std::sqrt(static_cast<Real>(3.0));
        const Real cos4pi3 = cos2pi3;
        const Real sin4pi3 = -sin2pi3;

        const Real dx = warp1
           + cos2pi3 * warp2
           + cos4pi3 * warp3;

        const Real dy = static_cast<Real>(0.0)
           + sin2pi3 * warp2
           + sin4pi3 * warp3;

        return { dx, dy };
      }
  };

  /**
   * @brief Computes 3D warp-blend shift for a tetrahedral face.
   *
   * Extends the 2D face warp to 3D by computing displacement in the
   * tangent plane of a tetrahedral face.
   *
   * @tparam K Polynomial degree.
   */
  template <size_t K>
  class WarpShiftFace3D
  {
    public:
      static std::array<Real, 2> apply(Real La, Real Lb, Real Lc, Real Ld, Real alpha)
      {
        (void) La;
        return WarpShiftFace2D<K>::apply(Lb, Lc, Ld, alpha);
      }
  };

  /**
   * @brief Applies warp-blend algorithm to move triangle nodes toward Fekete positions.
   *
   * Given an array of equispaced nodes on the reference triangle, this class
   * applies the warp-blend algorithm from Hesthaven & Warburton to move them
   * toward optimal Fekete-type positions that minimize the Lebesgue constant.
   *
   * The algorithm:
   * 1. Converts reference coordinates to equilateral triangle coordinates
   * 2. Computes edge-based warp contributions using 1D GLL warping
   * 3. Blends contributions using barycentric weights with parameter α
   * 4. Converts back to reference triangle coordinates
   *
   * @tparam K Polynomial degree.
   *
   * @see WarpFactor1D, TriangleBlend, FeketeTriangle
   */
  template <size_t K>
  class WarpBlendTriangle
  {
    public:
      template <size_t N>
      static void apply(std::array<Math::SpatialPoint, N>& nodes)
      {
        static_assert(
          N == (K + 1) * (K + 2) / 2,
          "WarpBlendTriangle<K>::apply: N must be (K+1)(K+2)/2."
        );

        for (auto& p : nodes)
          p.resize(2);

        if constexpr (K <= 1)
          return;

        // -------------------------------------------------------------------
        // Helpers
        // -------------------------------------------------------------------

        // Row offset for the (i,j) enumeration used in FeketeTriangle:
        // for j in [0..K], i in [0..K-j]:
        //   idx(j,i) = sum_{m=0}^{j-1} (K+1-m) + i
        static constexpr auto rowOffset = [](size_t j) -> size_t
        {
          size_t off = 0;
          for (size_t m = 0; m < j; ++m)
            off += (K + 1 - m);
          return off;
        };

        constexpr Real SQRT3     = static_cast<Real>(1.7320508075688772);
        constexpr Real INV_SQRT3 = static_cast<Real>(1.0) / SQRT3;

        const Real alpha = TriangleBlend<K>::getAlpha();

        // -------------------------------------------------------------------
        // 1. Warp only strictly interior nodes with the 2D warp–blend
        // -------------------------------------------------------------------
        // Reference triangle vertices: v1=(0,0), v2=(1,0), v3=(0,1).
        //
        // Equally spaced enumeration from FeketeTriangle:
        //   for j=0..K:
        //     for i=0..K-j:
        //       nodes[idx] = (i/K, j/K)
        //
        // Vertices:
        //   v1: (i=0,   j=0)   -> idx = 0
        //   v2: (i=K,   j=0)   -> idx = K
        //   v3: (i=0,   j=K)   -> idx = rowOffset(K)
        //
        // Edges:
        //   e12: j=0, i=0..K
        //   e31: i=0, j=0..K
        //   e23: i+j=K

        size_t idx = 0;
        for (size_t j = 0; j <= K; ++j)
        {
          const size_t rowOff = rowOffset(j);
          for (size_t i = 0; i <= K - j; ++i, ++idx)
          {
            // Skip vertices
            const bool is_v1 = (i == 0   && j == 0);
            const bool is_v2 = (i == K   && j == 0);
            const bool is_v3 = (i == 0   && j == K);
            if (is_v1 || is_v2 || is_v3)
              continue;

            // Skip edge nodes (we will overwrite them with exact GLL nodes)
            const bool on_e12 = (j == 0);         // v1–v2
            const bool on_e31 = (i == 0);         // v3–v1
            const bool on_e23 = (i + j == K);     // v2–v3
            if (on_e12 || on_e31 || on_e23)
              continue;

            // Strictly interior node: apply full 2D warp–blend.

            // Start from equispaced coordinates (i/K, j/K)
            const Real xref = static_cast<Real>(i) / static_cast<Real>(K);
            const Real yref = static_cast<Real>(j) / static_cast<Real>(K);

            // Barycentric on reference triangle
            Real L1 = static_cast<Real>(1.0) - xref - yref; // at (0,0)
            Real L2 = xref;                                 // at (1,0)
            Real L3 = yref;                                 // at (0,1)

            // To equilateral coordinates
            Real x = -L2 + L3;
            Real y = (-L2 - L3 + static_cast<Real>(2.0) * L1) * INV_SQRT3;

            // 2D warp–blend shift
            const auto d = WarpShiftFace2D<K>::apply(L1, L2, L3, alpha);
            const Real dx = d[0];
            const Real dy = d[1];

            x += dx;
            y += dy;

            // Back: equilateral → barycentric (inverse mapping)
            L1 =  y * INV_SQRT3 + static_cast<Real>(1.0) / static_cast<Real>(3.0);
            L2 = -static_cast<Real>(0.5) * x
               - static_cast<Real>(0.5) * y * INV_SQRT3
               + static_cast<Real>(1.0) / static_cast<Real>(3.0);
            L3 =  static_cast<Real>(0.5) * x
               - static_cast<Real>(0.5) * y * INV_SQRT3
               + static_cast<Real>(1.0) / static_cast<Real>(3.0);

            // Normalize barycentric (should be very close to 1)
            const Real sumL = L1 + L2 + L3;
            const Real invSum = static_cast<Real>(1.0) / sumL;
            L1 *= invSum;
            L2 *= invSum;
            L3 *= invSum;

            // Back to reference triangle: (x,y) = (L2,L3)
            assert(rowOff + i < nodes.size());
            nodes[rowOff + i].resize(2);
            assert(nodes[rowOff + i].size() == 2);
            nodes[rowOff + i][0] = L2;
            nodes[rowOff + i][1] = L3;
          }
        }

        // -------------------------------------------------------------------
        // 2. Overwrite edge nodes with exact GLL01<K> positions
        // -------------------------------------------------------------------

        const auto& gll01 = GLL01<K>::getNodes();

        // Edge e12: from v1=(0,0) to v2=(1,0), j=0, i=0..K
        // Param t ∈ [0,1], t = x = L2; L1 = 1 - t, L3 = 0.
        for (size_t i_edge = 0; i_edge <= K; ++i_edge)
        {
          const size_t idx_edge = i_edge; // rowOffset(0) + i_edge = i_edge
          const Real t = gll01[i_edge];
          const Real L2 = t;
          const Real L3 = static_cast<Real>(0.0);

          assert(idx_edge < nodes.size());
          assert(nodes[idx_edge].size() == 2);
          nodes[idx_edge][0] = L2;
          nodes[idx_edge][1] = L3;
        }

        // Edge e31: from v1=(0,0) to v3=(0,1), i=0, j=0..K
        // Param t ∈ [0,1], t = y = L3; L1 = 1 - t, L2 = 0.
        for (size_t j_edge = 0; j_edge <= K; ++j_edge)
        {
          const size_t idx_edge = rowOffset(j_edge); // i=0
          const Real t = gll01[j_edge];

          const Real L2 = static_cast<Real>(0.0);
          const Real L3 = t;

          assert(idx_edge < nodes.size());
          assert(nodes[idx_edge].size() == 2);
          nodes[idx_edge][0] = L2;
          nodes[idx_edge][1] = L3;
        }

        // Edge e23: from v2=(1,0) to v3=(0,1), nodes with i+j=K.
        // Param t ∈ [0,1] along v2→v3; we take t = L3.
        // On this edge: L1 = 0, L2 = 1 - t, L3 = t.
        //
        // Equispaced enumeration: for each j in [0..K], i = K - j.
        for (size_t j_edge = 0; j_edge <= K; ++j_edge)
        {
          const size_t i_edge = K - j_edge;
          const size_t idx_edge = rowOffset(j_edge) + i_edge;
          const Real t = gll01[j_edge];

          const Real L3 = t;
          const Real L2 = static_cast<Real>(1.0) - t;

          assert(idx_edge < nodes.size());
          assert(nodes[idx_edge].size() == 2);
          nodes[idx_edge][0] = L2;
          nodes[idx_edge][1] = L3;
        }
      }
  };

  /**
   * @brief Applies warp-blend algorithm to move tetrahedron nodes toward Fekete positions.
   *
   * Extends the triangular warp-blend to 3D tetrahedra. The algorithm:
   * 1. Converts to equilateral tetrahedron coordinates
   * 2. Computes face-based warp contributions from all 4 faces
   * 3. Special handling for edge nodes (2 zero barycentric coords) using 1D GLL warp
   * 4. Blends contributions using 4D barycentric weights
   * 5. Converts back to reference tetrahedron coordinates
   *
   * Edge nodes are handled separately to ensure they match GLL01 positions
   * exactly, which is essential for H1 conformity between adjacent elements.
   *
   * @tparam K Polynomial degree.
   *
   * @see WarpBlendTriangle, TetrahedronBlend, FeketeTetrahedron
   */
  template <size_t K>
  class WarpBlendTetrahedron
  {
    private:
      // -------------------------------------------------------------------
      // Helper: index mapping (i,j,k) -> flat idx (no lambdas; avoids ASan
      // stack-use-after-scope reports from short-lifetime closure objects).
      // -------------------------------------------------------------------
      static constexpr size_t layerOffset(size_t kk)
      {
        size_t off = 0;
        for (size_t m = 0; m < kk; ++m)
        {
          const size_t n = K - m;
          off += (n + 1) * (n + 2) / 2;
        }
        return off;
      }

      static constexpr size_t rowOffsetWithinLayer(size_t kk, size_t jj)
      {
        size_t off = 0;
        const size_t n = K - kk;
        for (size_t r = 0; r < jj; ++r)
          off += (n + 1 - r);
        return off;
      }

      static constexpr size_t idxOf(size_t i, size_t j, size_t k)
      {
        return layerOffset(k) + rowOffsetWithinLayer(k, j) + i;
      }

      static inline void set_from_bary(
        Math::SpatialPoint& p, Real L2, Real L3, Real L4)
      {
        assert(p.size() == 3);
        p[0] = L2;
        p[1] = L3;
        p[2] = L4;
      }

    public:
      template <size_t N>
      static void apply(std::array<Math::SpatialPoint, N>& nodes)
      {
        static_assert(
          N == (K + 1) * (K + 2) * (K + 3) / 6,
          "WarpBlendTetrahedron<K>::apply: N must be (K + 1)(K + 2)(K + 3) / 6."
        );

        for (auto& p : nodes)
          p.resize(3);

        if constexpr (K <= 1)
          return;

        using Real = Rodin::Real;
        constexpr Real TOL = static_cast<Real>(RODIN_VARIATIONAL_H1_WARPBLEND_TOLERANCE);

        const Real alpha  = TetrahedronBlend<K>::getAlpha();
        const Real alphaT = TriangleBlend<K>::getAlpha();
        const Real invK   = static_cast<Real>(1.0) / static_cast<Real>(K);

        // -------------------------------------------------------------------
        // Equilateral tetrahedron geometry (as before)
        // -------------------------------------------------------------------
        const Real invSqrt3 = static_cast<Real>(1.0) / std::sqrt(static_cast<Real>(3.0));
        const Real invSqrt6 = static_cast<Real>(1.0) / std::sqrt(static_cast<Real>(6.0));

        const Real v1x = static_cast<Real>(-1.0);
        const Real v1y = -invSqrt3;
        const Real v1z = -invSqrt6;

        const Real v2x = static_cast<Real>(1.0);
        const Real v2y = -invSqrt3;
        const Real v2z = -invSqrt6;

        const Real v3x = static_cast<Real>(0.0);
        const Real v3y = static_cast<Real>(2.0) * invSqrt3;
        const Real v3z = -invSqrt6;

        const Real v4x = static_cast<Real>(0.0);
        const Real v4y = static_cast<Real>(0.0);
        const Real v4z = static_cast<Real>(3.0) * invSqrt6;

        // Tangent directions for each face
        Real t1[4][3];
        Real t2[4][3];

        // Face 0: opposite v4, vertices (v1,v2,v3)
        t1[0][0] = v2x - v1x;
        t1[0][1] = v2y - v1y;
        t1[0][2] = v2z - v1z;
        t2[0][0] = v3x - static_cast<Real>(0.5) * (v1x + v2x);
        t2[0][1] = v3y - static_cast<Real>(0.5) * (v1y + v2y);
        t2[0][2] = v3z - static_cast<Real>(0.5) * (v1z + v2z);

        // Face 1: opposite v3, vertices (v1,v2,v4)
        t1[1][0] = v2x - v1x;
        t1[1][1] = v2y - v1y;
        t1[1][2] = v2z - v1z;
        t2[1][0] = v4x - static_cast<Real>(0.5) * (v1x + v2x);
        t2[1][1] = v4y - static_cast<Real>(0.5) * (v1y + v2y);
        t2[1][2] = v4z - static_cast<Real>(0.5) * (v1z + v2z);

        // Face 2: opposite v2, vertices (v1,v3,v4)
        t1[2][0] = v3x - v2x;
        t1[2][1] = v3y - v2y;
        t1[2][2] = v3z - v2z;
        t2[2][0] = v4x - static_cast<Real>(0.5) * (v2x + v3x);
        t2[2][1] = v4y - static_cast<Real>(0.5) * (v2y + v3y);
        t2[2][2] = v4z - static_cast<Real>(0.5) * (v2z + v3z);

        // Face 3: opposite v1, vertices (v2,v3,v4)
        t1[3][0] = v3x - v1x;
        t1[3][1] = v3y - v1y;
        t1[3][2] = v3z - v1z;
        t2[3][0] = v4x - static_cast<Real>(0.5) * (v1x + v3x);
        t2[3][1] = v4y - static_cast<Real>(0.5) * (v1y + v3y);
        t2[3][2] = v4z - static_cast<Real>(0.5) * (v1z + v3z);

        // Normalize t1,t2
        for (int f = 0; f < 4; ++f)
        {
          const Real n1 = std::sqrt(t1[f][0]*t1[f][0] + t1[f][1]*t1[f][1] + t1[f][2]*t1[f][2]);
          const Real n2 = std::sqrt(t2[f][0]*t2[f][0] + t2[f][1]*t2[f][1] + t2[f][2]*t2[f][2]);

          t1[f][0] /= n1; t1[f][1] /= n1; t1[f][2] /= n1;
          t2[f][0] /= n2; t2[f][1] /= n2; t2[f][2] /= n2;
        }

        // Inverse of A = [v1-v4, v2-v4, v3-v4] for barycentric recovery
        const Real a11 = static_cast<Real>(-0.5);
        const Real a12 = -static_cast<Real>(0.5) * invSqrt3;
        const Real a13 = -static_cast<Real>(0.5) * invSqrt6;

        const Real a21 = static_cast<Real>(0.5);
        const Real a22 = -static_cast<Real>(0.5) * invSqrt3;
        const Real a23 = -static_cast<Real>(0.5) * invSqrt6;

        const Real a31 = static_cast<Real>(0.0);
        const Real a32 = invSqrt3;
        const Real a33 = -static_cast<Real>(0.5) * invSqrt6;

        // -------------------------------------------------------------------
        // 1. Warp strictly interior volume nodes (4 barycentric ints > 0)
        // -------------------------------------------------------------------
        for (size_t k = 0; k <= K; ++k)
        {
          const size_t lo  = layerOffset(k);
          const size_t nJK = K - k;

          for (size_t j = 0; j <= nJK; ++j)
          {
            const size_t ro = rowOffsetWithinLayer(k, j);

            for (size_t i = 0; i <= nJK - j; ++i)
            {
              const size_t idx = lo + ro + i;

              const size_t l2 = i;
              const size_t l3 = j;
              const size_t l4 = k;
              const size_t l1 = K - i - j - k;

              // Hoist these to the full loop-iteration scope:
              // avoids ASan stack-use-after-scope reports caused by
              // compiler lifetime shortening and later stack spills.
              Real l1n = static_cast<Real>(0.0);
              Real l2n = static_cast<Real>(0.0);
              Real l3n = static_cast<Real>(0.0);
              Real l4n = static_cast<Real>(0.0);

              const int nz =
                (l1 > 0 ? 1 : 0) +
                (l2 > 0 ? 1 : 0) +
                (l3 > 0 ? 1 : 0) +
                (l4 > 0 ? 1 : 0);

              // interior volume node: all four > 0
              if (nz != 4)
                continue;

              Real L1 = static_cast<Real>(l1) * invK;
              Real L2 = static_cast<Real>(l2) * invK;
              Real L3 = static_cast<Real>(l3) * invK;
              Real L4 = static_cast<Real>(l4) * invK;

              // Equilateral coordinates of undeformed point
              const Real rx0 = L1 * v1x + L2 * v2x + L3 * v3x + L4 * v4x;
              const Real ry0 = L1 * v1y + L2 * v2y + L3 * v3y + L4 * v4y;
              const Real rz0 = L1 * v1z + L2 * v2z + L3 * v3z + L4 * v4z;

              Real shiftx = static_cast<Real>(0.0);
              Real shifty = static_cast<Real>(0.0);
              Real shiftz = static_cast<Real>(0.0);

              // Face contributions (volume-blended)
              for (int face = 0; face < 4; ++face)
              {
                Real La, Lb, Lc, Ld;

                if (face == 0) { La = L1; Lb = L2; Lc = L3; Ld = L4; }
                else if (face == 1) { La = L2; Lb = L1; Lc = L3; Ld = L4; }
                else if (face == 2) { La = L3; Lb = L1; Lc = L4; Ld = L2; }
                else { /* face == 3 */ La = L4; Lb = L1; Lc = L3; Ld = L2; }

                // Avoid structured bindings temporaries (keeps lifetime simple under ASan).
                const auto w = WarpShiftFace3D<K>::apply(La, Lb, Lc, Ld, alpha);
                const Real warp1 = w[0];
                const Real warp2 = w[1];

                Real blend = Lb * Lc * Ld;
                const Real denom = (Lb + static_cast<Real>(0.5) * La)
                                 * (Lc + static_cast<Real>(0.5) * La)
                                 * (Ld + static_cast<Real>(0.5) * La);

                if (denom > TOL)
                {
                  const Real alphaLa = alpha * La;
                  blend = (static_cast<Real>(1.0) + alphaLa * alphaLa) * blend / denom;
                }
                else
                {
                  blend = static_cast<Real>(0.0);
                }

                const Real bw1 = blend * warp1;
                const Real bw2 = blend * warp2;

                shiftx += bw1 * t1[face][0] + bw2 * t2[face][0];
                shifty += bw1 * t1[face][1] + bw2 * t2[face][1];
                shiftz += bw1 * t1[face][2] + bw2 * t2[face][2];
              }

              const Real rx = rx0 + shiftx;
              const Real ry = ry0 + shifty;
              const Real rz = rz0 + shiftz;

              // Recover barycentric w.r.t v1..v4
              const Real dxv = rx - v4x;
              const Real dyv = ry - v4y;
              const Real dzv = rz - v4z;

              l1n = a11 * dxv + a12 * dyv + a13 * dzv;
              l2n = a21 * dxv + a22 * dyv + a23 * dzv;
              l3n = a31 * dxv + a32 * dyv + a33 * dzv;
              l4n = static_cast<Real>(1.0) - l1n - l2n - l3n;

              l1n = std::max(static_cast<Real>(0.0), l1n);
              l2n = std::max(static_cast<Real>(0.0), l2n);
              l3n = std::max(static_cast<Real>(0.0), l3n);
              l4n = std::max(static_cast<Real>(0.0), l4n);

              Real sumL = l1n + l2n + l3n + l4n;
              if (sumL > TOL)
              {
                const Real invSum = static_cast<Real>(1.0) / sumL;
                l1n *= invSum;
                l2n *= invSum;
                l3n *= invSum;
                l4n *= invSum;
              }

              // Back to reference tetra: (x,y,z) = (L2,L3,L4)
              assert(idx < nodes.size());
              assert(nodes[idx].size() == 3);
              nodes[idx][0] = l2n;
              nodes[idx][1] = l3n;
              nodes[idx][2] = l4n;
            }
          }
        }

        // -------------------------------------------------------------------
        // 2. Per-face 2D warp for strictly face-interior nodes
        //    (exactly one barycentric integer is 0, others > 0)
        // -------------------------------------------------------------------
        constexpr Real SQRT3     = static_cast<Real>(1.7320508075688772);
        constexpr Real INV_SQRT3 = static_cast<Real>(1.0) / SQRT3;

        for (size_t k = 0; k <= K; ++k)
        {
          const size_t lo  = layerOffset(k);
          const size_t nJK = K - k;

          for (size_t j = 0; j <= nJK; ++j)
          {
            const size_t ro = rowOffsetWithinLayer(k, j);

            for (size_t i = 0; i <= nJK - j; ++i)
            {
              const size_t idx = lo + ro + i;

              const size_t l2 = i;
              const size_t l3 = j;
              const size_t l4 = k;
              const size_t l1 = K - i - j - k;

              const int nzero =
                (l1 == 0 ? 1 : 0) +
                (l2 == 0 ? 1 : 0) +
                (l3 == 0 ? 1 : 0) +
                (l4 == 0 ? 1 : 0);

              const int nz =
                (l1 > 0 ? 1 : 0) +
                (l2 > 0 ? 1 : 0) +
                (l3 > 0 ? 1 : 0) +
                (l4 > 0 ? 1 : 0);

              // Face-interior node: exactly one zero, three strictly positive
              if (nzero != 1 || nz != 3)
                continue;

              Real L[4];
              L[0] = static_cast<Real>(l1) * invK;
              L[1] = static_cast<Real>(l2) * invK;
              L[2] = static_cast<Real>(l3) * invK;
              L[3] = static_cast<Real>(l4) * invK;

              int i_zero = -1;
              if (l1 == 0) i_zero = 0;
              else if (l2 == 0) i_zero = 1;
              else if (l3 == 0) i_zero = 2;
              else i_zero = 3;

              int ia, ib, ic;
              if (i_zero == 0) { ia = 1; ib = 2; ic = 3; }       // face v2-v3-v4 (L1=0)
              else if (i_zero == 1) { ia = 0; ib = 2; ic = 3; }  // face v1-v3-v4 (L2=0)
              else if (i_zero == 2) { ia = 0; ib = 1; ic = 3; }  // face v1-v2-v4 (L3=0)
              else { ia = 0; ib = 1; ic = 2; }                   // face v1-v2-v3 (L4=0)

              Real La = L[ia];
              Real Lb = L[ib];
              Real Lc = L[ic];

              // (La,Lb,Lc) already sum to 1, but normalize for safety
              const Real sum = La + Lb + Lc;
              const Real invSum = static_cast<Real>(1.0) / sum;
              La *= invSum;
              Lb *= invSum;
              Lc *= invSum;

              // Map to equilateral triangle coordinates
              Real x = -Lb + Lc;
              Real y = (-Lb - Lc + static_cast<Real>(2.0) * La) * INV_SQRT3;

              // 2D warp–blend on face
              const auto d = WarpShiftFace2D<K>::apply(La, Lb, Lc, alphaT);
              const Real dx = d[0];
              const Real dy = d[1];

              x += dx;
              y += dy;

              // Back: equilateral → barycentric (La_new,Lb_new,Lc_new)
              Real La_new = y * INV_SQRT3 + static_cast<Real>(1.0) / static_cast<Real>(3.0);
              Real Lb_new = -static_cast<Real>(0.5) * x
                          - static_cast<Real>(0.5) * y * INV_SQRT3
                          + static_cast<Real>(1.0) / static_cast<Real>(3.0);
              Real Lc_new =  static_cast<Real>(0.5) * x
                          - static_cast<Real>(0.5) * y * INV_SQRT3
                          + static_cast<Real>(1.0) / static_cast<Real>(3.0);

              La_new = std::max(static_cast<Real>(0.0), La_new);
              Lb_new = std::max(static_cast<Real>(0.0), Lb_new);
              Lc_new = std::max(static_cast<Real>(0.0), Lc_new);

              Real sumL = La_new + Lb_new + Lc_new;
              if (sumL > TOL)
              {
                const Real invSumL = static_cast<Real>(1.0) / sumL;
                La_new *= invSumL;
                Lb_new *= invSumL;
                Lc_new *= invSumL;
              }

              // Update barycentrics: zero stays 0, others get (La_new,Lb_new,Lc_new)
              Real Lnew[4] = {static_cast<Real>(0.0),
                              static_cast<Real>(0.0),
                              static_cast<Real>(0.0),
                              static_cast<Real>(0.0)};

              Lnew[ia] = La_new;
              Lnew[ib] = Lb_new;
              Lnew[ic] = Lc_new;

              // Back to reference tetra: (x,y,z) = (L2,L3,L4)
              assert(idx < nodes.size());
              assert(nodes[idx].size() == 3);
              nodes[idx][0] = Lnew[1];
              nodes[idx][1] = Lnew[2];
              nodes[idx][2] = Lnew[3];
            }
          }
        }

        // -------------------------------------------------------------------
        // 3. Overwrite edge nodes with exact GLL01<K> positions (combinatorial)
        // -------------------------------------------------------------------
        const auto& gll01 = GLL01<K>::getNodes();

        // Edge v1-v2: (0,0,0)–(1,0,0): j=0, k=0, i=0..K
        for (size_t i = 0; i <= K; ++i)
        {
          const size_t id = idxOf(i, 0, 0);
          const Real t  = gll01[i];
          const Real L2 = t;
          const Real L3 = static_cast<Real>(0.0);
          const Real L4 = static_cast<Real>(0.0);

          assert(id < nodes.size());
          set_from_bary(nodes[id], L2, L3, L4);
        }

        // Edge v1-v3: (0,0,0)–(0,1,0): i=0, k=0, j=0..K
        for (size_t j = 0; j <= K; ++j)
        {
          const size_t id = idxOf(0, j, 0);
          const Real t  = gll01[j];
          const Real L2 = static_cast<Real>(0.0);
          const Real L3 = t;
          const Real L4 = static_cast<Real>(0.0);

          assert(id < nodes.size());
          set_from_bary(nodes[id], L2, L3, L4);
        }

        // Edge v1-v4: (0,0,0)–(0,0,1): i=0, j=0, k=0..K
        for (size_t k = 0; k <= K; ++k)
        {
          const size_t id = idxOf(0, 0, k);
          const Real t  = gll01[k];
          const Real L2 = static_cast<Real>(0.0);
          const Real L3 = static_cast<Real>(0.0);
          const Real L4 = t;

          assert(id < nodes.size());
          set_from_bary(nodes[id], L2, L3, L4);
        }

        // Edge v2-v3: (1,0,0)–(0,1,0): k=0, i+j=K
        //   L1 = 0, L2 = 1-t, L3 = t, L4 = 0
        for (size_t j = 0; j <= K; ++j)
        {
          const size_t i  = K - j;
          const size_t id = idxOf(i, j, 0);

          const Real t  = gll01[j]; // parameter from v2→v3: t = L3
          const Real L2 = static_cast<Real>(1.0) - t;
          const Real L3 = t;
          const Real L4 = static_cast<Real>(0.0);

          assert(id < nodes.size());
          set_from_bary(nodes[id], L2, L3, L4);
        }

        // Edge v2-v4: (1,0,0)–(0,0,1): j=0, i+k=K
        //   L1 = 0, L2 = 1-t, L4 = t, L3 = 0
        for (size_t k = 0; k <= K; ++k)
        {
          const size_t i  = K - k;
          const size_t id = idxOf(i, 0, k);

          const Real t  = gll01[k]; // parameter from v2→v4: t = L4
          const Real L2 = static_cast<Real>(1.0) - t;
          const Real L3 = static_cast<Real>(0.0);
          const Real L4 = t;

          assert(id < nodes.size());
          set_from_bary(nodes[id], L2, L3, L4);
        }

        // Edge v3-v4: (0,1,0)–(0,0,1): i=0, j+k=K
        //   L1 = 0, L3 = 1-t, L4 = t, L2 = 0
        for (size_t j = 0; j <= K; ++j)
        {
          const size_t k  = K - j;
          const size_t id = idxOf(0, j, k);

          const Real t  = gll01[k]; // parameter from v3→v4: t = L4
          const Real L2 = static_cast<Real>(0.0);
          const Real L3 = static_cast<Real>(1.0) - t;
          const Real L4 = t;

          assert(id < nodes.size());
          set_from_bary(nodes[id], L2, L3, L4);
        }
      }
  };
}

#endif // RODIN_VARIATIONAL_H1_WARPBLEND_H
