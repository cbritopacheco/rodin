/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_WARPBLEND_H
#define RODIN_VARIATIONAL_H1_WARPBLEND_H

#include <cstddef>

#include "Rodin/Math/Vector.h"
#include "Rodin/Types.h"

#include "GLL.h"

#include "LagrangeBasis.h"

#define RODIN_VARIATIONAL_H1_WARPBLEND_TOLERANCE 1e-14

namespace Rodin::Variational
{
  template <size_t K, size_t MaxItGLL = 25>
  class WarpBlend
  {
    public:
      /**
       * @brief Computes the warp factor for a given edge in the triangle
       * This is used to move equispaced nodes toward optimal Fekete positions.
       */
      static Real getFactor(Real r, Real s)
      {
        constexpr Real TOL = RODIN_VARIATIONAL_H1_WARPBLEND_TOLERANCE;

        if constexpr (K == 0)
        {
          return 0.0;
        }
        else
        {
          // Get GLL nodes on [-1, 1]
          const auto& gll_nodes = GLL<K, MaxItGLL>::getNodes();

          // Evaluate warp based on distance from edge
          Real warp = 0.0;

          // Compute barycentric coordinate along edge
          const Real lambda = r / (r + s + TOL);  // Add small epsilon to avoid division by zero

          // Map to [-1, 1]
          const Real xi = 2.0 * lambda - 1.0;

          // Compute warp as weighted sum of Lagrange polynomials
          for (size_t i = 1; i < K; ++i)  // Skip endpoints
          {
            const Real L = LagrangeBasis1D<K>(GLL<K, MaxItGLL>::getNodes()).getBasis(i, xi);
            const Real target = gll_nodes[i];
            const Real equi = 2.0 * static_cast<Real>(i) / static_cast<Real>(K) - 1.0;
            warp += L * (target - equi);
          }

          return warp;
        }
      }
  };

  template <size_t K, size_t MaxItGLL = 25>
  class WarpBlendTriangle
  {
    public:
      template <size_t N>
      static void apply(std::array<Math::SpatialPoint, N>& nodes)
      {
        constexpr Real TOL = RODIN_VARIATIONAL_H1_WARPBLEND_TOLERANCE;

        if constexpr (K <= 1)
        {
          return;  // No warping needed for linear elements
        }
        else
        {
          // Convert to equilateral triangle coordinates
          std::array<std::pair<Real, Real>, N> equi_coords;
          for (size_t idx = 0; idx < N; ++idx)
          {
            const Real x = nodes[idx].x();
            const Real y = nodes[idx].y();

            // Convert from reference triangle (x,y) to equilateral coordinates (r,s)
            // Reference triangle: (0,0), (1,0), (0,1)
            // Equilateral triangle: centered coordinate system
            const Real r = 2.0 * x - 1.0 + y;
            const Real s = Math::sqrt(3.0) * y - 1.0;

            equi_coords[idx] = { r, s };
          }

          // Apply warp for each node
          for (size_t idx = 0; idx < N; ++idx)
          {
            // Note: equi_coords computed above but not used in simplified warp implementation
            // (void)equi_coords;  // Suppress unused warning

            const Real x = nodes[idx].x();
            const Real y = nodes[idx].y();
            const Real z = 1.0 - x - y;  // Third barycentric coordinate

            // Skip vertices (they should remain fixed)
            if (x < TOL && y < TOL)  // Vertex (0,0)
              continue;
            if (x > 1.0 - TOL && y < TOL)  // Vertex (1,0)
              continue;
            if (x < TOL && y > 1.0 - TOL)  // Vertex (0,1)
              continue;

            // Compute warp contributions from each edge
            Real warp1 = 0.0, warp2 = 0.0, warp3 = 0.0;

            // Edge 1: from (1,0) to (0,1), perpendicular direction
            if (x + y > TOL)
              warp1 = WarpBlend<K, MaxItGLL>::getFactor(x, y);

            // Edge 2: from (0,0) to (0,1), perpendicular direction
            if (y + z > TOL)
              warp2 = WarpBlend<K, MaxItGLL>::getFactor(y, z);

            // Edge 3: from (0,0) to (1,0), perpendicular direction
            if (z + x > TOL)
              warp3 = WarpBlend<K, MaxItGLL>::getFactor(z, x);

            // Blend the warp contributions using barycentric coordinates as weights
            // This ensures smooth transition and maintains symmetry
            Real blend1 = y * z;
            Real blend2 = z * x;
            Real blend3 = x * y;
            Real blend_sum = blend1 + blend2 + blend3 + TOL;

            // Apply scaled warp in each direction
            Real scale = 1.0;  // Scaling factor for warp magnitude
            Real dx = scale * (blend1 * warp1 + blend2 * warp2 + blend3 * warp3) / blend_sum;

            // Update node position (warp is applied in x-direction as approximation)
            nodes[idx] = Math::SpatialPoint{{x + dx * 0.5, y - dx * 0.5 * std::sqrt(3.0)}};

            // Clamp to valid triangle domain
            Real new_x = nodes[idx].x();
            Real new_y = nodes[idx].y();
            new_x = std::max(0.0, std::min(1.0, new_x));
            new_y = std::max(0.0, std::min(1.0, new_y));
            if (new_x + new_y > 1.0)
            {
              Real excess = new_x + new_y - 1.0;
              new_x -= excess * 0.5;
              new_y -= excess * 0.5;
            }
            nodes[idx] = Math::SpatialPoint{{ new_x, new_y }};
          }
        }
      }
  };

  template <size_t K, size_t MaxItGLL = 25>
  class WarpBlendTetrahedron
  {
    public:
      template <size_t N>
      static void apply(std::array<Math::SpatialPoint, N>& nodes)
      {
        constexpr Real TOL = RODIN_VARIATIONAL_H1_WARPBLEND_TOLERANCE;

        if constexpr (K <= 1)
        {
          return;  // No warping needed for linear elements
        }
        else
        {
          // Apply warp for each node
          for (size_t idx = 0; idx < N; ++idx)
          {
            Real x = nodes[idx].x();
            Real y = nodes[idx].y();
            Real z = nodes[idx].z();
            Real w = 1.0 - x - y - z;  // Fourth barycentric coordinate

            // Skip vertices (they should remain fixed)
            bool is_vertex =
              (x < TOL && y < TOL && z < TOL) ||  // Vertex (0,0,0)
              (x > 1.0 - TOL && y < TOL && z < TOL) ||  // Vertex (1,0,0)
              (x < TOL && y > 1.0-TOL && z < TOL) ||  // Vertex (0,1,0)
              (x < TOL && y < TOL && z > 1.0 - TOL);    // Vertex (0,0,1)

            if (is_vertex)
              continue;

            // Compute warp contributions from each face/edge
            // Simplified approach: warp based on distance from faces
            Real warp_x = 0.0, warp_y = 0.0, warp_z = 0.0;

            // Face 1: opposite to vertex (1,0,0)
            if (y + z + w > TOL)
            {
              warp_x += warpFactor<K>(y, z + w) * 0.3;
            }

            // Face 2: opposite to vertex (0,1,0)
            if (x + z + w > TOL)
            {
              warp_y += warpFactor<K>(x, z + w) * 0.3;
            }

            // Face 3: opposite to vertex (0,0,1)
            if (x + y + w > TOL)
            {
              warp_z += warpFactor<K>(x, y + w) * 0.3;
            }

            // Blend using barycentric coordinates
            Real blend_factor = x * y * z * w;
            Real scale = 1.0 / (1.0 + K * K * 0.1);  // Reduce warp for higher orders

            // Apply warp
            Real new_x = x + warp_x * scale * blend_factor;
            Real new_y = y + warp_y * scale * blend_factor;
            Real new_z = z + warp_z * scale * blend_factor;

            // Clamp to valid tetrahedron domain
            new_x = std::max(0.0, std::min(1.0, new_x));
            new_y = std::max(0.0, std::min(1.0, new_y));
            new_z = std::max(0.0, std::min(1.0, new_z));

            if (new_x + new_y + new_z > 1.0)
            {
              Real excess = new_x + new_y + new_z - 1.0;
              new_x -= excess / 3.0;
              new_y -= excess / 3.0;
              new_z -= excess / 3.0;
            }

            nodes[idx] = Math::SpatialPoint{{new_x, new_y, new_z}};
          }
        }
      }
  };
}

#endif
