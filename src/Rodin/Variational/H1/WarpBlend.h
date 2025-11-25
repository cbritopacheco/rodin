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
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/Common.h"

#include "GLL.h"

#define RODIN_VARIATIONAL_H1_WARPBLEND_TOLERANCE 1e-14

namespace Rodin::Variational
{
  //======================================================================
  // 1D warp factor: equispaced -> GLL (→ GLL01 on [0,1])
  //
  // Concept:
  //   - req_j: equispaced nodes in [-1,1]
  //   - rGL_j: GLL nodes in [-1,1]
  //   - difference d_j = rGL_j - req_j
  //   - warp(r) = ∑_j L_j(r) d_j, where L_j are Lagrange polynomials on req_j
  //
  // Using GLL01 on [0,1] is consistent because GLL01<K> is defined
  // from GLL<K> by an affine map.
  //======================================================================

  template <size_t K>
  class WarpFactor1D
  {
    public:
      static Real get(Real r)
      {
        if constexpr (K <= 1)
          return 0.0;

        // Equispaced nodes on [-1, 1]
        static constexpr auto buildReq = []()
        {
          std::array<Real, K + 1> req{};
          for (size_t j = 0; j <= K; ++j)
            req[j] = static_cast<Real>(-1.0)
                   + static_cast<Real>(2.0) * static_cast<Real>(j) / static_cast<Real>(K);
          return req;
        };

        static const std::array<Real, K + 1> req = buildReq();

        // GLL nodes on [-1,1]
        const auto& rGL = GLL<K>::getNodes();

        // Interpolate the displacement (rGL - req) at r
        Real warp = 0.0;
        for (size_t i = 0; i <= K; ++i)
        {
          Real Li = 1.0;
          const Real xi = req[i];
          for (size_t m = 0; m <= K; ++m)
          {
            if (m == i) continue;
            const Real xm = req[m];
            Li *= (r - xm) / (xi - xm);
          }
          warp += Li * (rGL[i] - req[i]);
        }

        // This warp already vanishes at the endpoints (because rGL and req coincide).
        // No extra scaling is strictly necessary; if you want, you can
        // introduce the usual (1 - r^2) “compression” here.
        return warp;
      }
  };

  //======================================================================
  // Triangle warp–blend
  //
  // Input nodes are in reference triangle with vertices:
  //   V1 = (0,0), V2 = (1,0), V3 = (0,1).
  //
  // Steps:
  //  1. Compute barycentric L1,L2,L3.
  //  2. For edge nodes (one Li ≈ 0), apply pure 1D GLL01 edge warp tangentially.
  //  3. Convert barycentric -> equilateral coordinates (x,y).
  //  4. For all nodes, compute 1D warp values in edge coordinates:
  //       r1 = L3 - L2, r2 = L1 - L3, r3 = L2 - L1,
  //     and get w1,w2,w3 from WarpFactor1D.
  //  5. Blend these edge warps into the interior with a symmetric
  //     blending function B = 4 L1 L2 L3.
  //  6. Update equilateral coordinates (x,y), then map back to
  //     barycentric and then to reference (s,t).
  //
  // Notes:
  //   - Edge warps ensure that edge nodes are on the 1D GLL01 grid.
  //   - The interior warp provides “Fekete-type” clustering.
  //======================================================================

  template <size_t K>
  class WarpBlendTriangle
  {
    public:
      template <size_t N>
      static void apply(std::array<Math::SpatialPoint, N>& nodes, Real alpha = 0)
      {
        constexpr Real TOL        = static_cast<Real>(RODIN_VARIATIONAL_H1_WARPBLEND_TOLERANCE);
        constexpr Real SQRT3      = static_cast<Real>(1.7320508075688772); // ≈ sqrt(3)
        constexpr Real INV_SQRT3  = static_cast<Real>(1.0) / SQRT3;

        if constexpr (K <= 1)
        {
          // K=0: single vertex, K=1: standard P1, no warp.
          return;
        }

        // -------------------------------------------------------------------
        // Step A: 1D warp along each boundary edge to put edge nodes on GLL01
        // -------------------------------------------------------------------
        for (auto& p : nodes)
        {
          Real s = p.x();
          Real t = p.y();
          Real L1 = static_cast<Real>(1.0) - s - t;
          Real L2 = s;
          Real L3 = t;

          // Classify edges (exclude vertices)
          const bool on_e23 = (Math::abs(L1) < TOL) && (L2 > TOL) && (L3 > TOL); // edge V2–V3
          const bool on_e13 = (Math::abs(L2) < TOL) && (L1 > TOL) && (L3 > TOL); // edge V1–V3
          const bool on_e12 = (Math::abs(L3) < TOL) && (L1 > TOL) && (L2 > TOL); // edge V1–V2

          if (on_e23)
          {
            // Edge between V2 (L2=1,L3=0) and V3 (L2=0,L3=1)
            // Param in [-1,1]: r = L3 - L2
            Real r = L3 - L2;
            Real warp = WarpFactor1D<K>::get(r);
            Real r_new = r + warp;

            // Recover L2,L3 on that edge (L1=0, L2+L3=1, L3-L2 = r_new)
            L2 = (static_cast<Real>(1.0) - r_new) * static_cast<Real>(0.5);
            L3 = (static_cast<Real>(1.0) + r_new) * static_cast<Real>(0.5);
            L1 = static_cast<Real>(0.0);
          }
          else if (on_e13)
          {
            // Edge between V1 (L1=1,L3=0) and V3 (L1=0,L3=1)
            // Param in [-1,1]: r = L1 - L3
            Real r = L1 - L3;
            Real warp = WarpFactor1D<K>::get(r);
            Real r_new = r + warp;

            // L1+L3=1, L1-L3=r_new
            L1 = (static_cast<Real>(1.0) + r_new) * static_cast<Real>(0.5);
            L3 = (static_cast<Real>(1.0) - r_new) * static_cast<Real>(0.5);
            L2 = static_cast<Real>(0.0);
          }
          else if (on_e12)
          {
            // Edge between V1 (L1=1,L2=0) and V2 (L1=0,L2=1)
            // Param in [-1,1]: r = L2 - L1
            Real r = L2 - L1;
            Real warp = WarpFactor1D<K>::get(r);
            Real r_new = r + warp;

            // L1+L2=1, L2-L1=r_new
            L2 = (static_cast<Real>(1.0) + r_new) * static_cast<Real>(0.5);
            L1 = (static_cast<Real>(1.0) - r_new) * static_cast<Real>(0.5);
            L3 = static_cast<Real>(0.0);
          }

          // Project back to reference coordinates
          p = Math::SpatialPoint{{L2, L3}};
        }

        // -------------------------------------------------------------------
        // Step B: interior warp–blend (Fekete-type clustering + α factor)
        // -------------------------------------------------------------------
        for (auto& p : nodes)
        {
          Real s = p.x();
          Real t = p.y();
          Real L1 = static_cast<Real>(1.0) - s - t;
          Real L2 = s;
          Real L3 = t;

          // Vertices: do nothing
          if ( (Math::abs(L1 - 1.0) < TOL) ||
               (Math::abs(L2 - 1.0) < TOL) ||
               (Math::abs(L3 - 1.0) < TOL) )
          {
            continue;
          }

          // Convert to equilateral coordinates
          // v1 = (0,  2/√3), v2 = (-1,-1/√3), v3 = (1,-1/√3)
          Real x = -L2 + L3;
          Real y = (-L2 - L3 + static_cast<Real>(2.0)*L1) * INV_SQRT3;

          // 1D edge coordinates in [-1,1]
          Real r1 = L3 - L2; // edge 2–3 (opposite vertex 1)
          Real r2 = L1 - L3; // edge 1–3 (opposite vertex 2)
          Real r3 = L2 - L1; // edge 1–2 (opposite vertex 3)

          Real w1 = WarpFactor1D<K>::get(r1);
          Real w2 = WarpFactor1D<K>::get(r2);
          Real w3 = WarpFactor1D<K>::get(r3);

          // Base blending factor: vanishes on edges, max at interior
          Real blend = static_cast<Real>(4.0) * L1 * L2 * L3;

          // α-dependent quadratic factors using opposite vertex barycentric
          Real a1 = static_cast<Real>(1.0) + (alpha * L1) * (alpha * L1); // edge 2–3, opp. L1
          Real a2 = static_cast<Real>(1.0) + (alpha * L2) * (alpha * L2); // edge 1–3, opp. L2
          Real a3 = static_cast<Real>(1.0) + (alpha * L3) * (alpha * L3); // edge 1–2, opp. L3

          // Apply blending + α
          Real dw1 = blend * a1 * w1;
          Real dw2 = blend * a2 * w2;
          Real dw3 = blend * a3 * w3;

          // Symmetric combination (equilateral)
          Real dx =  dw2 - dw3;
          Real dy = (static_cast<Real>(2.0) * dw1 - dw2 - dw3) * INV_SQRT3;

          x += dx;
          y += dy;

          // Map back: equilateral -> barycentric
          L1 = y * INV_SQRT3 + static_cast<Real>(1.0) / static_cast<Real>(3.0);
          L2 = -static_cast<Real>(0.5) * x
             - static_cast<Real>(0.5) * y * INV_SQRT3
             + static_cast<Real>(1.0) / static_cast<Real>(3.0);
          L3 =  static_cast<Real>(0.5) * x
             - static_cast<Real>(0.5) * y * INV_SQRT3
             + static_cast<Real>(1.0) / static_cast<Real>(3.0);

          // Clamp and renormalize
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

          p = Math::SpatialPoint{{L2, L3}};
        }
      }
  };

  //----------------------------------------------------------------------
  // Helper: degree-dependent α for tetrahedron (Table 10.1 in H&W)
  //----------------------------------------------------------------------

  template <size_t K>
  static constexpr Real getTetraAlpha()
  {
    if constexpr (K <= 3)
    {
      return static_cast<Real>(0.0);
    }
    else
    {
      // Optimized α for N = 1..15 (N=degree)
      //  N :   1      2      3      4       5       6       7       8
      // αN : 0.0,   0.0,   0.0, 0.1002, 1.1332, 1.5608, 1.3413, 1.2577,
      //   9      10      11      12      13      14      15
      // 1.1603, 1.0153, 0.6080, 0.4523, 0.8856, 0.8717, 0.9655
      constexpr Real table[16] =
      {
        static_cast<Real>(0.0), // dummy for index 0
        static_cast<Real>(0.0), // 1
        static_cast<Real>(0.0), // 2
        static_cast<Real>(0.0), // 3
        static_cast<Real>(0.1002), // 4
        static_cast<Real>(1.1332), // 5
        static_cast<Real>(1.5608), // 6
        static_cast<Real>(1.3413), // 7
        static_cast<Real>(1.2577), // 8
        static_cast<Real>(1.1603), // 9
        static_cast<Real>(1.0153), // 10
        static_cast<Real>(0.6080), // 11
        static_cast<Real>(0.4523), // 12
        static_cast<Real>(0.8856), // 13
        static_cast<Real>(0.8717), // 14
        static_cast<Real>(0.9655)  // 15
      };

      if constexpr (K < 16)
        return table[K];
      else
        // For K > 15, just freeze α at the last optimized value.
        return table[15];
    }
  }

  //----------------------------------------------------------------------
  // Helper: 2D Warp–Blend shift on a triangle face (evalshift.m)
  //
  // Input: barycentric (L1,L2,L3) on the face, and α (triangle parameter).
  // Output: (dx,dy) in equilateral triangle coordinates.
  //----------------------------------------------------------------------

  template <size_t K>
  inline void evalShiftTriangle(Real L1, Real L2, Real L3,
                                Real alpha, Real& dx, Real& dy)
  {
    if constexpr (K <= 1)
    {
      dx = dy = static_cast<Real>(0.0);
      return;
    }

    // Blending per edge
    const Real blend1 = L2 * L3;
    const Real blend2 = L1 * L3;
    const Real blend3 = L1 * L2;

    // Edge coordinates in [-1,1]
    const Real r1 = L3 - L2;
    const Real r2 = L1 - L3;
    const Real r3 = L2 - L1;

    // 1D warp (equispaced -> GLL) along each edge, with factor 4
    const Real warpf1 = static_cast<Real>(4.0) * WarpFactor1D<K>::get(r1);
    const Real warpf2 = static_cast<Real>(4.0) * WarpFactor1D<K>::get(r2);
    const Real warpf3 = static_cast<Real>(4.0) * WarpFactor1D<K>::get(r3);

    const Real aL1 = alpha * L1;
    const Real aL2 = alpha * L2;
    const Real aL3 = alpha * L3;

    const Real warp1 = blend1 * warpf1 * (static_cast<Real>(1.0) + aL1 * aL1);
    const Real warp2 = blend2 * warpf2 * (static_cast<Real>(1.0) + aL2 * aL2);
    const Real warp3 = blend3 * warpf3 * (static_cast<Real>(1.0) + aL3 * aL3);

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

  //======================================================================
  // Tetrahedral Warp–Blend (Hesthaven–Warburton style)
  //
  // - Uses 2D warp–blend on each face (via evalShiftTriangle)
  // - Builds 3D face warps using orthonormal face tangents t_{f,1}, t_{f,2}
  // - Blends the four face warps with b_f(λ) and β (α_3D) in barycentric space
  // - Works in an origin-centered equilateral tetrahedron and maps back to
  //   the standard reference tetrahedron with vertices (0,0,0),(1,0,0),(0,1,0),(0,0,1).
  //======================================================================

  template <size_t K>
  class WarpBlendTetrahedron
  {
    public:
      template <size_t N>
      static void apply(std::array<Math::SpatialPoint, N>& nodes)
      {
        constexpr Real TOL   = static_cast<Real>(RODIN_VARIATIONAL_H1_WARPBLEND_TOLERANCE);

        if constexpr (K <= 1)
          return;

        // α for triangle faces (usual choice 5/3 for p>=4)
        const Real alphaTri =
          (K >= 4 ? static_cast<Real>(5.0) / static_cast<Real>(3.0)
                  : static_cast<Real>(0.0));

        // β for 3D blending (degree-dependent)
        const Real beta = getTetraAlpha<K>();

        // Geometry of the equilateral tetrahedron (Hesthaven–Warburton / Göbel)
        const Real invSqrt2 = static_cast<Real>(1.0) / std::sqrt(static_cast<Real>(2.0));
        const Real invSqrt3 = static_cast<Real>(1.0) / std::sqrt(static_cast<Real>(3.0));
        const Real invSqrt6 = static_cast<Real>(1.0) / std::sqrt(static_cast<Real>(6.0));

        // Vertices v1..v4 in R^3 (origin-centered equilateral tetrahedron)
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

        // Inverse of A = [v1-v4, v2-v4, v3-v4], for barycentric recovery
        const Real a11 = static_cast<Real>(-0.5);
        const Real a12 = -static_cast<Real>(0.5) * invSqrt3;
        const Real a13 = -static_cast<Real>(0.5) * invSqrt6;

        const Real a21 = static_cast<Real>(0.5);
        const Real a22 = -static_cast<Real>(0.5) * invSqrt3;
        const Real a23 = -static_cast<Real>(0.5) * invSqrt6;

        const Real a31 = static_cast<Real>(0.0);
        const Real a32 = invSqrt3;
        const Real a33 = -static_cast<Real>(0.5) * invSqrt6;

        // Orthonormal face tangents t_{f,1}, t_{f,2} (Göbel 2024, Sec. 4.2)
        const Real t11x = static_cast<Real>(1.0);
        const Real t11y = static_cast<Real>(0.0);
        const Real t11z = static_cast<Real>(0.0);
        const Real t12x = static_cast<Real>(0.0);
        const Real t12y = static_cast<Real>(1.0);
        const Real t12z = static_cast<Real>(0.0);

        const Real t21x = static_cast<Real>(1.0);
        const Real t21y = static_cast<Real>(0.0);
        const Real t21z = static_cast<Real>(0.0);
        const Real t22x = static_cast<Real>(0.0);
        const Real t22y = static_cast<Real>(1.0) / static_cast<Real>(3.0);
        const Real t22z = static_cast<Real>(4.0)
                        / (static_cast<Real>(3.0) * std::sqrt(static_cast<Real>(2.0)));

        const Real t31x = -static_cast<Real>(0.5);
        const Real t31y = static_cast<Real>(3.0)
                        / (static_cast<Real>(2.0) * std::sqrt(static_cast<Real>(3.0)));
        const Real t31z = static_cast<Real>(0.0);
        const Real t32x = -static_cast<Real>(0.5) * invSqrt3;
        const Real t32y = -static_cast<Real>(1.0) / static_cast<Real>(6.0);
        const Real t32z = static_cast<Real>(4.0)
                        / (static_cast<Real>(3.0) * std::sqrt(static_cast<Real>(2.0)));

        const Real t41x =  static_cast<Real>(0.5);
        const Real t41y =  static_cast<Real>(3.0)
                         / (static_cast<Real>(2.0) * std::sqrt(static_cast<Real>(3.0)));
        const Real t41z =  static_cast<Real>(0.0);
        const Real t42x =  static_cast<Real>(0.5) * invSqrt3;
        const Real t42y = -static_cast<Real>(1.0) / static_cast<Real>(6.0);
        const Real t42z =  static_cast<Real>(4.0)
                         / (static_cast<Real>(3.0) * std::sqrt(static_cast<Real>(2.0)));

        for (auto& p : nodes)
        {
          // Standard reference tetra coordinates
          Real xr = p.x();
          Real yr = p.y();
          Real zr = p.z();

          // Barycentric coordinates w.r.t. standard simplex
          Real l1 = static_cast<Real>(1.0) - xr - yr - zr; // vertex (0,0,0)
          Real l2 = xr;                                   // vertex (1,0,0)
          Real l3 = yr;                                   // vertex (0,1,0)
          Real l4 = zr;                                   // vertex (0,0,1)

          // Skip vertices (keep them fixed)
          const bool near_v1 = (l1 > static_cast<Real>(1.0) - TOL);
          const bool near_v2 = (l2 > static_cast<Real>(1.0) - TOL);
          const bool near_v3 = (l3 > static_cast<Real>(1.0) - TOL);
          const bool near_v4 = (l4 > static_cast<Real>(1.0) - TOL);
          if (near_v1 || near_v2 || near_v3 || near_v4)
            continue;

          // Equilateral coordinates of the current point
          const Real rx0 = l1 * v1x + l2 * v2x + l3 * v3x + l4 * v4x;
          const Real ry0 = l1 * v1y + l2 * v2y + l3 * v3y + l4 * v4y;
          const Real rz0 = l1 * v1z + l2 * v2z + l3 * v3z + l4 * v4z;

          //------------------------------------------------------------------
          // Face warps w1..w4 using 2D Warp–Blend on appropriate triples
          //------------------------------------------------------------------
          Real dxF1 = 0, dyF1 = 0;
          Real dxF2 = 0, dyF2 = 0;
          Real dxF3 = 0, dyF3 = 0;
          Real dxF4 = 0, dyF4 = 0;

          // Face F1 (opposite λ1), triple (λ2,λ3,λ4)
          evalShiftTriangle<K>(l2, l3, l4, alphaTri, dxF1, dyF1);

          // Face F2 (opposite λ2), triple (λ1,λ3,λ4)
          evalShiftTriangle<K>(l1, l3, l4, alphaTri, dxF2, dyF2);

          // Face F3 (opposite λ3), triple (λ1,λ4,λ2)
          evalShiftTriangle<K>(l1, l4, l2, alphaTri, dxF3, dyF3);

          // Face F4 (opposite λ4), triple (λ1,λ3,λ2)
          evalShiftTriangle<K>(l1, l3, l2, alphaTri, dxF4, dyF4);

          // Build 3D warps w_f
          const Real w1x = dxF1 * t11x + dyF1 * t12x;
          const Real w1y = dxF1 * t11y + dyF1 * t12y;
          const Real w1z = dxF1 * t11z + dyF1 * t12z;

          const Real w2x = dxF2 * t21x + dyF2 * t22x;
          const Real w2y = dxF2 * t21y + dyF2 * t22y;
          const Real w2z = dxF2 * t21z + dyF2 * t22z;

          const Real w3x = dxF3 * t31x + dyF3 * t32x;
          const Real w3y = dxF3 * t31y + dyF3 * t32y;
          const Real w3z = dxF3 * t31z + dyF3 * t32z;

          const Real w4x = dxF4 * t41x + dyF4 * t42x;
          const Real w4y = dxF4 * t41y + dyF4 * t42y;
          const Real w4z = dxF4 * t41z + dyF4 * t42z;

          //------------------------------------------------------------------
          // Face blending b_f(λ) and β-dependent factors (Göbel (4.7))
          //------------------------------------------------------------------
          auto safe_ratio = [&](Real num, Real den) -> Real
          {
            if (Math::abs(den) < TOL)
              return static_cast<Real>(0.0);
            return num / den;
          };

          const Real b1 =
            safe_ratio(static_cast<Real>(2.0) * l2, static_cast<Real>(2.0) * l2 + l1) *
            safe_ratio(static_cast<Real>(2.0) * l3, static_cast<Real>(2.0) * l3 + l1) *
            safe_ratio(static_cast<Real>(2.0) * l4, static_cast<Real>(2.0) * l4 + l1);

          const Real b2 =
            safe_ratio(static_cast<Real>(2.0) * l1, static_cast<Real>(2.0) * l1 + l2) *
            safe_ratio(static_cast<Real>(2.0) * l3, static_cast<Real>(2.0) * l3 + l2) *
            safe_ratio(static_cast<Real>(2.0) * l4, static_cast<Real>(2.0) * l4 + l2);

          const Real b3 =
            safe_ratio(static_cast<Real>(2.0) * l1, static_cast<Real>(2.0) * l1 + l3) *
            safe_ratio(static_cast<Real>(2.0) * l2, static_cast<Real>(2.0) * l2 + l3) *
            safe_ratio(static_cast<Real>(2.0) * l4, static_cast<Real>(2.0) * l4 + l3);

          const Real b4 =
            safe_ratio(static_cast<Real>(2.0) * l1, static_cast<Real>(2.0) * l1 + l4) *
            safe_ratio(static_cast<Real>(2.0) * l2, static_cast<Real>(2.0) * l2 + l4) *
            safe_ratio(static_cast<Real>(2.0) * l3, static_cast<Real>(2.0) * l3 + l4);

          const Real betaL1 = beta * l1;
          const Real betaL2 = beta * l2;
          const Real betaL3 = beta * l3;
          const Real betaL4 = beta * l4;

          const Real f1 = b1 * (static_cast<Real>(1.0) + betaL1 * betaL1);
          const Real f2 = b2 * (static_cast<Real>(1.0) + betaL2 * betaL2);
          const Real f3 = b3 * (static_cast<Real>(1.0) + betaL3 * betaL3);
          const Real f4 = b4 * (static_cast<Real>(1.0) + betaL4 * betaL4);

          const Real gx = f1 * w1x + f2 * w2x + f3 * w3x + f4 * w4x;
          const Real gy = f1 * w1y + f2 * w2y + f3 * w3y + f4 * w4y;
          const Real gz = f1 * w1z + f2 * w2z + f3 * w3z + f4 * w4z;

          // New equilateral coordinates
          const Real rx = rx0 + gx;
          const Real ry = ry0 + gy;
          const Real rz = rz0 + gz;

          //------------------------------------------------------------------
          // Recover barycentric λ' from equilateral coordinates
          //------------------------------------------------------------------
          const Real dxv = rx - v4x;
          const Real dyv = ry - v4y;
          const Real dzv = rz - v4z;

          Real l1n = a11 * dxv + a12 * dyv + a13 * dzv;
          Real l2n = a21 * dxv + a22 * dyv + a23 * dzv;
          Real l3n = a31 * dxv + a32 * dyv + a33 * dzv;
          Real l4n = static_cast<Real>(1.0) - l1n - l2n - l3n;

          // Clamp and renormalize barycentric coordinates
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

          // Back to standard reference tetra: (x,y,z) = (λ2,λ3,λ4)
          xr = l2n;
          yr = l3n;
          zr = l4n;

          p = Math::SpatialPoint{{xr, yr, zr}};
        }
      }
  };
}

#endif
