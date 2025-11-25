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

  //======================================================================
  // Tetrahedral warp–blend (face-based heuristic)
  //
  // Here I keep essentially the structure you already had, but use
  // WarpFactor1D<K> for face-based warps so that edge traces are
  // consistent with the 1D GLL01 distribution.
  //
  // This is a reasonable “Hesthaven–Warburton style” warp–blend,
  // but not guaranteed to match Nodes3D.m exactly.
  //======================================================================

  template <size_t K>
  class WarpBlendTetrahedron
  {
    public:
      template <size_t N>
      static void apply(std::array<Math::SpatialPoint, N>& nodes)
      {
        constexpr Real TOL      = static_cast<Real>(RODIN_VARIATIONAL_H1_WARPBLEND_TOLERANCE);
        constexpr Real ONE      = static_cast<Real>(1.0);
        constexpr Real HALF     = static_cast<Real>(0.5);

        if constexpr (K <= 1)
          return;

        for (auto& p : nodes)
        {
          Real x = p.x();
          Real y = p.y();
          Real z = p.z();
          Real w = ONE - x - y - z; // 4th barycentric

          // Skip vertices
          const bool is_vertex =
              ((Math::abs(x)     < TOL) && (Math::abs(y) < TOL) && (Math::abs(z) < TOL)) ||
              ((Math::abs(x-ONE) < TOL) && (Math::abs(y) < TOL) && (Math::abs(z) < TOL)) ||
              ((Math::abs(y-ONE) < TOL) && (Math::abs(x) < TOL) && (Math::abs(z) < TOL)) ||
              ((Math::abs(z-ONE) < TOL) && (Math::abs(x) < TOL) && (Math::abs(y) < TOL));
          if (is_vertex)
            continue;

          // Simple face-based warp: treat each face as a triangle and
          // apply a 1D-like warp in the appropriate direction.
          //
          // Face opposite (1,0,0): (y,z,w)
          Real warp_x = 0.0;
          if (y + z + w > TOL)
          {
            // Parameter in [-1,1] on that face: r_face ≈ (z - y)/(y+z+w)
            Real r_face = (z - y) / (y + z + w);
            Real wf = WarpFactor1D<K>::get(r_face);
            // scale: small, to avoid large distortions
            warp_x += HALF * wf * y * z * w;
          }

          // Face opposite (0,1,0): (x,z,w)
          Real warp_y = 0.0;
          if (x + z + w > TOL)
          {
            Real r_face = (z - x) / (x + z + w);
            Real wf = WarpFactor1D<K>::get(r_face);
            warp_y += HALF * wf * x * z * w;
          }

          // Face opposite (0,0,1): (x,y,w)
          Real warp_z = 0.0;
          if (x + y + w > TOL)
          {
            Real r_face = (y - x) / (x + y + w);
            Real wf = WarpFactor1D<K>::get(r_face);
            warp_z += HALF * wf * x * y * w;
          }

          // Global blending factor to keep nodes near faces mostly on faces.
          Real blend = x * y * z * w;
          Real scale = ONE / (ONE + static_cast<Real>(0.1) * static_cast<Real>(K * K));

          x += scale * blend * warp_x;
          y += scale * blend * warp_y;
          z += scale * blend * warp_z;

          // Clamp to tetra
          x = std::max(static_cast<Real>(0.0), std::min(ONE, x));
          y = std::max(static_cast<Real>(0.0), std::min(ONE, y));
          z = std::max(static_cast<Real>(0.0), std::min(ONE, z));

          Real sum = x + y + z;
          if (sum > ONE)
          {
            Real excess = sum - ONE;
            x -= excess / static_cast<Real>(3.0);
            y -= excess / static_cast<Real>(3.0);
            z -= excess / static_cast<Real>(3.0);
          }

          p = Math::SpatialPoint{{x, y, z}};
        }
      }
  };
}

#endif
