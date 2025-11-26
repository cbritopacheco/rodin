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

#include "Rodin/Types.h"
#include "Rodin/Math/Common.h"
#include "Rodin/Math/Vector.h"

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

        const Real one = static_cast<Real>(1.0);
        const Real tol = static_cast<Real>(1.0e-10);

        if (Math::abs(r) < one - tol)
        {
          const Real sf = one - r * r;
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
      static constexpr void apply(Real& dx, Real& dy,
          Real L1, Real L2, Real L3, Real alpha)
      {
        if constexpr (K <= 1)
        {
          dx = dy = static_cast<Real>(0.0);
          return;
        }

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

        dx = warp1
           + cos2pi3 * warp2
           + cos4pi3 * warp3;

        dy = static_cast<Real>(0.0)
           + sin2pi3 * warp2
           + sin4pi3 * warp3;
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
      static constexpr void apply(
          Real& warpx, Real& warpy,
          Real La, Real Lb, Real Lc, Real Ld, Real alpha)
      {
        (void) Ld;
        WarpShiftFace2D<K>::apply(warpx, warpy, Lb, Lc, Ld, alpha);
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
        if constexpr (K <= 1)
          return;

        constexpr Real TOL       = static_cast<Real>(RODIN_VARIATIONAL_H1_WARPBLEND_TOLERANCE);
        constexpr Real SQRT3     = static_cast<Real>(1.7320508075688772);
        constexpr Real INV_SQRT3 = static_cast<Real>(1.0) / SQRT3;

        const Real alpha = TriangleBlend<K>::getAlpha();

        for (auto& p : nodes)
        {
          Real xref = p.x();
          Real yref = p.y();

          // Barycentric on reference triangle
          Real L1 = static_cast<Real>(1.0) - xref - yref;
          Real L2 = xref;
          Real L3 = yref;

          // Vertices unchanged
          if ( (Math::abs(L1 - 1.0) < TOL) ||
               (Math::abs(L2 - 1.0) < TOL) ||
               (Math::abs(L3 - 1.0) < TOL) )
          {
            continue;
          }

          // To equilateral coordinates
          Real x = -L2 + L3;
          Real y = (-L2 - L3 + static_cast<Real>(2.0) * L1) * INV_SQRT3;

          // 2D warp–blend shift
          Real dx, dy;
          WarpShiftFace2D<K>::apply(dx, dy, L1, L2, L3, alpha);

          x += dx;
          y += dy;

          // Back: equilateral → barycentric (inverse of above)
          L1 = y * INV_SQRT3 + static_cast<Real>(1.0) / static_cast<Real>(3.0);
          L2 = -static_cast<Real>(0.5) * x
             - static_cast<Real>(0.5) * y * INV_SQRT3
             + static_cast<Real>(1.0) / static_cast<Real>(3.0);
          L3 =  static_cast<Real>(0.5) * x
             - static_cast<Real>(0.5) * y * INV_SQRT3
             + static_cast<Real>(1.0) / static_cast<Real>(3.0);

          // Clamp and renormalize (robustness)
          L1 = std::max(static_cast<Real>(0.0), L1);
          L2 = std::max(static_cast<Real>(0.0), L2);
          L3 = std::max(static_cast<Real>(0.0), L3);

          Real sumL = L1 + L2 + L3;
          if (sumL > TOL)
          {
            L1 /= sumL;
            L2 /= sumL;
            L3 /= sumL;
          }

          // Back to reference triangle: (x,y) = (L2,L3)
          p = Math::SpatialPoint{{L2, L3}};
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
    public:
      template <size_t N>
      static void apply(std::array<Math::SpatialPoint, N>& nodes)
      {
        if constexpr (K <= 1)
          return;

        constexpr Real TOL = static_cast<Real>(RODIN_VARIATIONAL_H1_WARPBLEND_TOLERANCE);

        const Real alpha = TetrahedronBlend<K>::getAlpha();

        // Equilateral tetra vertices (same as HW)
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

        // Tangent directions t1(face,:), t2(face,:)
        Real t1[4][3];
        Real t2[4][3];

        // Face 1
        t1[0][0] = v2x - v1x;
        t1[0][1] = v2y - v1y;
        t1[0][2] = v2z - v1z;
        t2[0][0] = v3x - static_cast<Real>(0.5) * (v1x + v2x);
        t2[0][1] = v3y - static_cast<Real>(0.5) * (v1y + v2y);
        t2[0][2] = v3z - static_cast<Real>(0.5) * (v1z + v2z);

        // Face 2
        t1[1][0] = v2x - v1x;
        t1[1][1] = v2y - v1y;
        t1[1][2] = v2z - v1z;
        t2[1][0] = v4x - static_cast<Real>(0.5) * (v1x + v2x);
        t2[1][1] = v4y - static_cast<Real>(0.5) * (v1y + v2y);
        t2[1][2] = v4z - static_cast<Real>(0.5) * (v1z + v2z);

        // Face 3
        t1[2][0] = v3x - v2x;
        t1[2][1] = v3y - v2y;
        t1[2][2] = v3z - v2z;
        t2[2][0] = v4x - static_cast<Real>(0.5) * (v2x + v3x);
        t2[2][1] = v4y - static_cast<Real>(0.5) * (v2y + v3y);
        t2[2][2] = v4z - static_cast<Real>(0.5) * (v2z + v3z);

        // Face 4
        t1[3][0] = v3x - v1x;
        t1[3][1] = v3y - v1y;
        t1[3][2] = v3z - v1z;
        t2[3][0] = v4x - static_cast<Real>(0.5) * (v1x + v3x);
        t2[3][1] = v4y - static_cast<Real>(0.5) * (v1y + v3y);
        t2[3][2] = v4z - static_cast<Real>(0.5) * (v1z + v3z);

        // Normalize tangents
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

        for (auto& p : nodes)
        {
          Real xr = p.x();
          Real yr = p.y();
          Real zr = p.z();

          // Barycentric on reference tetra
          Real L1 = static_cast<Real>(1.0) - xr - yr - zr; // vertex (0,0,0)
          Real L2 = xr;                                   // vertex (1,0,0)
          Real L3 = yr;                                   // vertex (0,1,0)
          Real L4 = zr;                                   // vertex (0,0,1)

          // Keep vertices fixed
          const bool near_v1 = (L1 > static_cast<Real>(1.0) - TOL);
          const bool near_v2 = (L2 > static_cast<Real>(1.0) - TOL);
          const bool near_v3 = (L3 > static_cast<Real>(1.0) - TOL);
          const bool near_v4 = (L4 > static_cast<Real>(1.0) - TOL);
          if (near_v1 || near_v2 || near_v3 || near_v4)
            continue;

          // Count how many barycentric coords are essentially zero
          const int nzero = (L1 < TOL ? 1 : 0) + (L2 < TOL ? 1 : 0)
                          + (L3 < TOL ? 1 : 0) + (L4 < TOL ? 1 : 0);

          // Edge nodes: exactly 2 barycentric coords are zero
          // For edge nodes, apply 1D warp along the edge direction only
          if (nzero == 2)
          {
            // Determine which edge this node is on
            // Edges connect pairs of vertices where exactly 2 L's are nonzero
            // Apply 1D GLL warp along that edge

            // Find the two non-zero barycentric coordinates
            Real La = 0, Lb = 0;
            int ia = -1, ib = -1;
            Real Lvals[4] = {L1, L2, L3, L4};
            for (int i = 0; i < 4; ++i)
            {
              if (Lvals[i] >= TOL)
              {
                if (ia < 0) { ia = i; La = Lvals[i]; }
                else { ib = i; Lb = Lvals[i]; }
              }
            }

            // Edge parameter in [0,1]: t = Lb / (La + Lb)
            Real t = Lb / (La + Lb);

            // Map to [-1,1] for warp calculation
            Real r = static_cast<Real>(2.0) * t - static_cast<Real>(1.0);

            // Get the RAW interpolated displacement (gll[i] - equispaced[i])
            // We need to compute this without the (1-r²) division that WarpFactor1D does
            // Use the same approach as WarpFactor1D but without the division
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

            const auto& gll = GLL<K>::getNodes();

            // Interpolate (gll - req) at r using Lagrange basis on req
            Real warp = static_cast<Real>(0.0);
            for (size_t i = 0; i <= K; ++i)
            {
              const Real Li = LagrangeBasis1D<K>::getBasis(i, r, s_req);
              warp += Li * (gll[i] - s_req[i]);
            }

            // The warp is directly in [-1,1] space, so apply it
            Real r_new = r + warp;

            // Clamp to valid range
            r_new = std::max(static_cast<Real>(-1.0), std::min(static_cast<Real>(1.0), r_new));

            // Map back to [0,1]
            Real t_new = (r_new + static_cast<Real>(1.0)) * static_cast<Real>(0.5);

            // Update barycentric coordinates
            Real La_new = static_cast<Real>(1.0) - t_new;
            Real Lb_new = t_new;

            // Reconstruct position (l1n not directly used but needed for completeness)
            [[maybe_unused]] Real l1n = 0;
            Real l2n = 0, l3n = 0, l4n = 0;
            if (ia == 0) l1n = La_new; else if (ia == 1) l2n = La_new;
            else if (ia == 2) l3n = La_new; else l4n = La_new;
            if (ib == 0) l1n = Lb_new; else if (ib == 1) l2n = Lb_new;
            else if (ib == 2) l3n = Lb_new; else l4n = Lb_new;

            p.x() = l2n;
            p.y() = l3n;
            p.z() = l4n;
            continue;
          }

          // Face nodes: exactly 1 barycentric coord is zero
          // Apply 2D warp-blend within the face plane using triangle warp-blend
          // This ensures face nodes match FeketeTriangle nodes for conformity
          if (nzero == 1)
          {
            // Determine which face this node is on (which L is zero)
            // Face opposite to vertex i has Li = 0
            Real Lvals[4] = {L1, L2, L3, L4};
            int i_zero = -1;
            for (int i = 0; i < 4; ++i)
            {
              if (Lvals[i] < TOL)
              {
                i_zero = i;
                break;
              }
            }

            // Get the three non-zero barycentric coordinates for this face
            // These form the barycentric coordinates within the triangular face
            Real La, Lb, Lc;
            int ia, ib, ic;

            // Map face vertices to reference triangle vertices
            // Face opposite L1 (i_zero=0): vertices 2,3,4 with L2,L3,L4
            // Face opposite L2 (i_zero=1): vertices 1,3,4 with L1,L3,L4
            // Face opposite L3 (i_zero=2): vertices 1,2,4 with L1,L2,L4
            // Face opposite L4 (i_zero=3): vertices 1,2,3 with L1,L2,L3 (base face z=0)
            if (i_zero == 0) { ia = 1; ib = 2; ic = 3; }
            else if (i_zero == 1) { ia = 0; ib = 2; ic = 3; }
            else if (i_zero == 2) { ia = 0; ib = 1; ic = 3; }
            else { ia = 0; ib = 1; ic = 2; } // i_zero == 3: base face

            La = Lvals[ia];
            Lb = Lvals[ib];
            Lc = Lvals[ic];

            // Normalize to get face barycentric coordinates
            Real sum = La + Lb + Lc;
            La /= sum;
            Lb /= sum;
            Lc /= sum;

            // Use the same warp-blend as triangle (FeketeTriangle uses WarpBlendTriangle)
            // But we need to apply it in face coordinates and then map back

            // For the reference triangle with vertices at (0,0), (1,0), (0,1):
            // - La corresponds to vertex (0,0) i.e., 1 - x - y
            // - Lb corresponds to vertex (1,0) i.e., x
            // - Lc corresponds to vertex (0,1) i.e., y
            // So on the face: xface = Lb, yface = Lc

            // Apply 2D warp-blend (same as WarpBlendTriangle)
            const Real alphaT = TriangleBlend<K>::getAlpha();

            // Convert to equilateral triangle coordinates
            constexpr Real SQRT3     = static_cast<Real>(1.7320508075688772);
            constexpr Real INV_SQRT3 = static_cast<Real>(1.0) / SQRT3;

            Real x = -Lb + Lc;
            Real y = (-Lb - Lc + static_cast<Real>(2.0) * La) * INV_SQRT3;

            // 2D warp–blend shift
            Real dx, dy;
            WarpShiftFace2D<K>::apply(dx, dy, La, Lb, Lc, alphaT);

            x += dx;
            y += dy;

            // Back: equilateral → barycentric (inverse of above)
            Real La_new = y * INV_SQRT3 + static_cast<Real>(1.0) / static_cast<Real>(3.0);
            Real Lb_new = -static_cast<Real>(0.5) * x
                        - static_cast<Real>(0.5) * y * INV_SQRT3
                        + static_cast<Real>(1.0) / static_cast<Real>(3.0);
            Real Lc_new =  static_cast<Real>(0.5) * x
                        - static_cast<Real>(0.5) * y * INV_SQRT3
                        + static_cast<Real>(1.0) / static_cast<Real>(3.0);

            // Clamp and renormalize
            La_new = std::max(static_cast<Real>(0.0), La_new);
            Lb_new = std::max(static_cast<Real>(0.0), Lb_new);
            Lc_new = std::max(static_cast<Real>(0.0), Lc_new);

            Real sumL = La_new + Lb_new + Lc_new;
            if (sumL > TOL)
            {
              La_new /= sumL;
              Lb_new /= sumL;
              Lc_new /= sumL;
            }

            // Map back to 3D tetrahedron barycentric coordinates
            // The zero coordinate stays zero (node stays on face)
            // l1n is assigned for code symmetry but its value is not used since
            // reference tetrahedron coordinates are (x,y,z) = (L2,L3,L4)
            Real l2n = 0, l3n = 0, l4n = 0;

            // Assign new face barycentric coords to the correct tet barycentric coords
            // Note: assignments to index 0 (L1) are no-ops since l1n isn't used
            if (ia == 1) l2n = La_new;
            else if (ia == 2) l3n = La_new;
            else if (ia == 3) l4n = La_new;

            if (ib == 1) l2n = Lb_new;
            else if (ib == 2) l3n = Lb_new;
            else if (ib == 3) l4n = Lb_new;

            if (ic == 1) l2n = Lc_new;
            else if (ic == 2) l3n = Lc_new;
            else if (ic == 3) l4n = Lc_new;

            // The i_zero coordinate stays 0 (already initialized to 0)

            // Back to reference tetra: (x,y,z) = (L2,L3,L4)
            p.x() = l2n;
            p.y() = l3n;
            p.z() = l4n;
            continue;
          }

          // Equilateral coordinates of undeformed point
          const Real rx0 = L1 * v1x + L2 * v2x + L3 * v3x + L4 * v4x;
          const Real ry0 = L1 * v1y + L2 * v2y + L3 * v3y + L4 * v4y;
          const Real rz0 = L1 * v1z + L2 * v2z + L3 * v3z + L4 * v4z;

          Real shiftx = 0.0;
          Real shifty = 0.0;
          Real shiftz = 0.0;

          // Loop over faces 1..4 (indices 0..3)
          for (int face = 0; face < 4; ++face)
          {
            Real La, Lb, Lc, Ld;

            if (face == 0) { La = L1; Lb = L2; Lc = L3; Ld = L4; }
            else if (face == 1) { La = L2; Lb = L1; Lc = L3; Ld = L4; }
            else if (face == 2) { La = L3; Lb = L1; Lc = L4; Ld = L2; }
            else { /* face == 3 */ La = L4; Lb = L1; Lc = L3; Ld = L2; }

            Real warp1, warp2;
            WarpShiftFace3D<K>::apply(warp1, warp2, La, Lb, Lc, Ld, alpha);

            // Volume blending
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

            // Boundary fix: pure face warp if La≈0 and not all (Lb,Lc,Ld) vertices
            const bool boundary_face =
              (La < TOL) &&
              ( (Lb > TOL ? 1 : 0)
              + (Lc > TOL ? 1 : 0)
              + (Ld > TOL ? 1 : 0) < 3);

            Real sx_face, sy_face, sz_face;
            if (boundary_face)
            {
              // Use unblended tangential warp
              sx_face = warp1 * t1[face][0] + warp2 * t2[face][0];
              sy_face = warp1 * t1[face][1] + warp2 * t2[face][1];
              sz_face = warp1 * t1[face][2] + warp2 * t2[face][2];
            }
            else
            {
              const Real bw1 = blend * warp1;
              const Real bw2 = blend * warp2;
              sx_face = bw1 * t1[face][0] + bw2 * t2[face][0];
              sy_face = bw1 * t1[face][1] + bw2 * t2[face][1];
              sz_face = bw1 * t1[face][2] + bw2 * t2[face][2];
            }

            shiftx += sx_face;
            shifty += sy_face;
            shiftz += sz_face;
          }

          // New equilateral coords
          const Real rx = rx0 + shiftx;
          const Real ry = ry0 + shifty;
          const Real rz = rz0 + shiftz;

          // Recover barycentric w.r.t v1..v4 via A^{-1}
          const Real dxv = rx - v4x;
          const Real dyv = ry - v4y;
          const Real dzv = rz - v4z;

          Real l1n = a11 * dxv + a12 * dyv + a13 * dzv;
          Real l2n = a21 * dxv + a22 * dyv + a23 * dzv;
          Real l3n = a31 * dxv + a32 * dyv + a33 * dzv;
          Real l4n = static_cast<Real>(1.0) - l1n - l2n - l3n;

          // Clamp and renormalize barycentric
          l1n = std::max(static_cast<Real>(0.0), l1n);
          l2n = std::max(static_cast<Real>(0.0), l2n);
          l3n = std::max(static_cast<Real>(0.0), l3n);
          l4n = std::max(static_cast<Real>(0.0), l4n);

          Real sumL = l1n + l2n + l3n + l4n;
          if (sumL > TOL)
          {
            l1n /= sumL;
            l2n /= sumL;
            l3n /= sumL;
            l4n /= sumL;
          }

          // Back to reference tetra: (x,y,z) = (λ2,λ3,λ4)
          p.x() = l2n;
          p.y() = l3n;
          p.z() = l4n;
        }
      }
  };
}

#endif // RODIN_VARIATIONAL_H1_WARPBLEND_H
