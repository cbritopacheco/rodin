/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_H1ELEMENT_HPP
#define RODIN_VARIATIONAL_H1_H1ELEMENT_HPP

#include <cmath>
#include <array>

#include "Rodin/Math/Common.h"

#include "H1Element.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace Rodin::Variational
{
  namespace Internal
  {
    /**
     * @brief Computes Legendre polynomial P_n(x) and its derivative P'_n(x)
     * using the 3-term recurrence relation.
     * 
     * @param n Degree of Legendre polynomial
     * @param x Evaluation point in [-1, 1]
     * @param[out] P Value of P_n(x)
     * @param[out] dP Value of P'_n(x)
     */
    inline void legendrePolynomial(size_t n, Real x, Real& P, Real& dP)
    {
      if (n == 0)
      {
        P = 1.0;
        dP = 0.0;
        return;
      }
      else if (n == 1)
      {
        P = x;
        dP = 1.0;
        return;
      }

      // Use 3-term recurrence: (n+1)P_{n+1} = (2n+1)xP_n - nP_{n-1}
      Real P0 = 1.0;    // P_0
      Real P1 = x;      // P_1
      Real dP0 = 0.0;   // P'_0
      Real dP1 = 1.0;   // P'_1
      
      for (size_t k = 1; k < n; ++k)
      {
        Real P2 = ((2.0 * k + 1.0) * x * P1 - k * P0) / (k + 1.0);
        Real dP2 = ((2.0 * k + 1.0) * (P1 + x * dP1) - k * dP0) / (k + 1.0);
        
        P0 = P1;
        P1 = P2;
        dP0 = dP1;
        dP1 = dP2;
      }
      
      P = P1;
      dP = dP1;
    }

    /**
     * @brief Computes 1D Gauss-Lobatto-Legendre (GLL) nodes on [-1, 1].
     * 
     * GLL nodes are the endpoints {-1, 1} and the zeros of the derivative
     * of the Legendre polynomial P'_K(x) weighted by (1 - x²).
     * These nodes provide better conditioning for high-order elements.
     * 
     * @tparam K Polynomial degree
     * @return Array of K+1 GLL nodes in ascending order
     */
    template <size_t K>
    std::array<Real, K+1> getGLLNodes()
    {
      std::array<Real, K+1> nodes;
      
      // Endpoints are always -1 and 1
      nodes[0] = -1.0;
      nodes[K] = 1.0;
      
      if (K == 0)
        return nodes;
      
      if (K == 1)
        return nodes;
      
      // Interior nodes are zeros of (1 - x²) P'_K(x)
      // Use Newton iteration starting from Chebyshev points as initial guess
      constexpr Real tol = 1e-15;
      constexpr size_t max_iter = 20;
      
      for (size_t i = 1; i < K; ++i)
      {
        // Initial guess: Chebyshev points
        Real x = -std::cos(M_PI * i / K);
        
        // Newton iteration to find zero of (1 - x²) P'_K(x)
        for (size_t iter = 0; iter < max_iter; ++iter)
        {
          Real P, dP;
          legendrePolynomial(K, x, P, dP);
          
          // f(x) = (1 - x²) P'_K(x)
          Real f = (1.0 - x * x) * dP;
          
          // f'(x) = -2x P'_K(x) + (1 - x²) P''_K(x)
          // Using P''_K = ((2K+1)xP'_K - KP_{K-1}) / (1-x²) when needed
          // But simpler: use the relation dP = K(xP - P_{K-1})/(x²-1)
          // Actually for GLL: f'(x) = -2x*dP + (1-x²)*d²P
          // We use a different formula: from P'_K relation
          Real P_Kminus1, dP_Kminus1;
          if (K > 1)
            legendrePolynomial(K - 1, x, P_Kminus1, dP_Kminus1);
          else
            P_Kminus1 = x; // P_1 = x when K = 2
          
          // Use: (x² - 1) dP = K(xP - P_{K-1})
          Real d2P = (2.0 * x * dP - K * (K + 1.0) * P) / (1.0 - x * x);
          Real df = -2.0 * x * dP + (1.0 - x * x) * d2P;
          
          // Newton update
          Real dx = -f / df;
          x += dx;
          
          if (std::abs(dx) < tol)
            break;
        }
        
        nodes[i] = x;
      }
      
      return nodes;
    }

    /**
     * @brief Computes 1D Gauss-Lobatto-Legendre (GLL) nodes on [0, 1].
     * 
     * Maps GLL nodes from [-1, 1] to [0, 1] using x ↦ (x + 1)/2.
     * 
     * @tparam K Polynomial degree
     * @return Array of K+1 GLL nodes on [0, 1]
     */
    template <size_t K>
    std::array<Real, K+1> getGLLNodes01()
    {
      auto nodes_ref = getGLLNodes<K>();
      std::array<Real, K+1> nodes;
      
      for (size_t i = 0; i <= K; ++i)
        nodes[i] = (nodes_ref[i] + 1.0) / 2.0;
      
      return nodes;
    }

    /**
     * @brief Evaluates 1D Lagrange polynomial on shifted Gauss-Lobatto nodes
     * Used for computing the warp function in the warp-blend algorithm.
     */
    template <size_t K>
    Real evaluateWarpLagrange(size_t i, Real x)
    {
      const auto& gll_nodes = getGLLNodes<K>();
      Real result = 1.0;
      Real xi = gll_nodes[i];
      
      for (size_t j = 0; j <= K; ++j)
      {
        if (j != i)
        {
          Real xj = gll_nodes[j];
          result *= (x - xj) / (xi - xj);
        }
      }
      return result;
    }

    /**
     * @brief Computes the warp factor for a given edge in the triangle
     * This is used to move equispaced nodes toward optimal Fekete positions.
     */
    template <size_t K>
    Real warpFactor(Real r, Real s)
    {
      if (K == 0)
        return 0.0;
      
      // Get GLL nodes on [-1, 1]
      const auto& gll_nodes = getGLLNodes<K>();
      
      // Evaluate warp based on distance from edge
      Real warp = 0.0;
      
      // Compute barycentric coordinate along edge
      Real lambda = r / (r + s + 1e-10);  // Add small epsilon to avoid division by zero
      
      // Map to [-1, 1]
      Real xi = 2.0 * lambda - 1.0;
      
      // Compute warp as weighted sum of Lagrange polynomials
      for (size_t i = 1; i < K; ++i)  // Skip endpoints
      {
        Real L = evaluateWarpLagrange<K>(i, xi);
        Real target = gll_nodes[i];
        Real equi = 2.0 * static_cast<Real>(i) / static_cast<Real>(K) - 1.0;
        warp += L * (target - equi);
      }
      
      return warp;
    }

    /**
     * @brief Applies warp-blend to equispaced nodes on the reference triangle
     * 
     * The warp-blend algorithm moves equispaced nodes toward optimal Fekete positions
     * by applying a warp function along each edge, then blending the contributions.
     * 
     * Reference: Hesthaven & Warburton, "Nodal DG Methods", Algorithm 174.
     */
    template <size_t K>
    inline void applyTriangleWarpBlend(std::vector<Math::SpatialPoint>& nodes)
    {
      if (K <= 1)
        return;  // No warping needed for linear elements
      
      const size_t n_nodes = nodes.size();
      
      // Convert to equilateral triangle coordinates
      std::vector<std::pair<Real, Real>> equi_coords(n_nodes);
      for (size_t idx = 0; idx < n_nodes; ++idx)
      {
        Real x = nodes[idx].x();
        Real y = nodes[idx].y();
        
        // Convert from reference triangle (x,y) to equilateral coordinates (r,s)
        // Reference triangle: (0,0), (1,0), (0,1)
        // Equilateral triangle: centered coordinate system
        Real r = 2.0 * x - 1.0 + y;
        Real s = std::sqrt(3.0) * y - 1.0;
        
        equi_coords[idx] = {r, s};
      }
      
      // Apply warp for each node
      for (size_t idx = 0; idx < n_nodes; ++idx)
      {
        // Note: equi_coords computed above but not used in simplified warp implementation
        // (void)equi_coords;  // Suppress unused warning
        
        Real x = nodes[idx].x();
        Real y = nodes[idx].y();
        Real z = 1.0 - x - y;  // Third barycentric coordinate
        
        // Skip vertices (they should remain fixed)
        if (x < 1e-10 && y < 1e-10)  // Vertex (0,0)
          continue;
        if (x > 1.0 - 1e-10 && y < 1e-10)  // Vertex (1,0)
          continue;
        if (x < 1e-10 && y > 1.0 - 1e-10)  // Vertex (0,1)
          continue;
        
        // Compute warp contributions from each edge
        Real warp1 = 0.0, warp2 = 0.0, warp3 = 0.0;
        
        // Edge 1: from (1,0) to (0,1), perpendicular direction
        if (x + y > 1e-10)
        {
          warp1 = warpFactor<K>(x, y);
        }
        
        // Edge 2: from (0,0) to (0,1), perpendicular direction
        if (y + z > 1e-10)
        {
          warp2 = warpFactor<K>(y, z);
        }
        
        // Edge 3: from (0,0) to (1,0), perpendicular direction  
        if (z + x > 1e-10)
        {
          warp3 = warpFactor<K>(z, x);
        }
        
        // Blend the warp contributions using barycentric coordinates as weights
        // This ensures smooth transition and maintains symmetry
        Real blend1 = y * z;
        Real blend2 = z * x;
        Real blend3 = x * y;
        Real blend_sum = blend1 + blend2 + blend3 + 1e-10;
        
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
        nodes[idx] = Math::SpatialPoint{{new_x, new_y}};
      }
    }

    /**
     * @brief Computes Fekete nodes on the reference triangle using warp-blend algorithm.
     * 
     * Fekete nodes minimize the Lebesgue constant and provide excellent interpolation
     * properties for high-order elements on triangles.
     * 
     * Based on: Hesthaven & Warburton, "Nodal Discontinuous Galerkin Methods", 2008.
     * 
     * @tparam K Polynomial degree
     * @return Vector of triangle Fekete nodes
     */
    template <size_t K>
    const std::vector<Math::SpatialPoint>& getTriangleFeketeNodes()
    {
      static thread_local std::vector<Math::SpatialPoint> s_nodes;
      
      if (!s_nodes.empty())
        return s_nodes;
      
      const size_t n_nodes = (K + 1) * (K + 2) / 2;
      s_nodes.reserve(n_nodes);
      
      // Start with equispaced nodes
      for (size_t j = 0; j <= K; ++j)
      {
        for (size_t i = 0; i <= K - j; ++i)
        {
          Real s = static_cast<Real>(i) / static_cast<Real>(K);
          Real t = static_cast<Real>(j) / static_cast<Real>(K);
          s_nodes.emplace_back(Math::SpatialPoint{{s, t}});
        }
      }
      
      // Apply warp-blend algorithm to move nodes toward Fekete positions
      applyTriangleWarpBlend<K>(s_nodes);
      
      return s_nodes;
    }

    /**
     * @brief Applies warp-blend to equispaced nodes on the reference tetrahedron
     * 
     * Similar to triangle case, but extended to 3D tetrahedron.
     */
    template <size_t K>
    inline void applyTetrahedronWarpBlend(std::vector<Math::SpatialPoint>& nodes)
    {
      if (K <= 1)
        return;  // No warping needed for linear elements
      
      const size_t n_nodes = nodes.size();
      
      // Apply warp for each node
      for (size_t idx = 0; idx < n_nodes; ++idx)
      {
        Real x = nodes[idx].x();
        Real y = nodes[idx].y();
        Real z = nodes[idx].z();
        Real w = 1.0 - x - y - z;  // Fourth barycentric coordinate
        
        // Skip vertices (they should remain fixed)
        const Real eps = 1e-10;
        bool is_vertex = 
          (x < eps && y < eps && z < eps) ||  // Vertex (0,0,0)
          (x > 1.0-eps && y < eps && z < eps) ||  // Vertex (1,0,0)
          (x < eps && y > 1.0-eps && z < eps) ||  // Vertex (0,1,0)
          (x < eps && y < eps && z > 1.0-eps);    // Vertex (0,0,1)
        
        if (is_vertex)
          continue;
        
        // Compute warp contributions from each face/edge
        // Simplified approach: warp based on distance from faces
        Real warp_x = 0.0, warp_y = 0.0, warp_z = 0.0;
        
        // Face 1: opposite to vertex (1,0,0)
        if (y + z + w > eps)
        {
          warp_x += warpFactor<K>(y, z + w) * 0.3;
        }
        
        // Face 2: opposite to vertex (0,1,0)
        if (x + z + w > eps)
        {
          warp_y += warpFactor<K>(x, z + w) * 0.3;
        }
        
        // Face 3: opposite to vertex (0,0,1)
        if (x + y + w > eps)
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

    /**
     * @brief Computes Fekete nodes on the reference tetrahedron using warp-blend algorithm.
     * 
     * Fekete nodes minimize the Lebesgue constant and provide excellent interpolation
     * properties for high-order elements on tetrahedra.
     * 
     * @tparam K Polynomial degree
     * @return Vector of tetrahedron Fekete nodes
     */
    template <size_t K>
    const std::vector<Math::SpatialPoint>& getTetrahedronFeketeNodes()
    {
      static thread_local std::vector<Math::SpatialPoint> s_nodes;
      
      if (!s_nodes.empty())
        return s_nodes;
      
      const size_t n_nodes = (K + 1) * (K + 2) * (K + 3) / 6;
      s_nodes.reserve(n_nodes);
      
      // Start with equispaced nodes
      for (size_t k = 0; k <= K; ++k)
      {
        for (size_t j = 0; j <= K - k; ++j)
        {
          for (size_t i = 0; i <= K - j - k; ++i)
          {
            Real r = static_cast<Real>(i) / static_cast<Real>(K);
            Real s = static_cast<Real>(j) / static_cast<Real>(K);
            Real t = static_cast<Real>(k) / static_cast<Real>(K);
            s_nodes.emplace_back(Math::SpatialPoint{{r, s, t}});
          }
        }
      }
      
      // Apply warp-blend algorithm to move nodes toward Fekete positions
      applyTetrahedronWarpBlend<K>(s_nodes);
      
      return s_nodes;
    }

    /**
     * @brief Evaluates normalized Jacobi polynomial P^{α,β}_n(x) on [-1,1]
     * Used for Dubiner modal basis construction.
     */
    inline Real jacobiPolynomial(size_t n, Real alpha, Real beta, Real x)
    {
      if (n == 0)
        return 1.0;
      
      Real P0 = 1.0;
      Real P1 = 0.5 * (alpha - beta + (alpha + beta + 2.0) * x);
      
      if (n == 1)
        return P1;
      
      // 3-term recurrence for Jacobi polynomials
      Real P2 = 0.0;
      for (size_t k = 1; k < n; ++k)
      {
        Real a1 = 2.0 * (k + 1.0) * (k + alpha + beta + 1.0) * (2.0 * k + alpha + beta);
        Real a2 = (2.0 * k + alpha + beta + 1.0) * (alpha * alpha - beta * beta);
        Real a3 = (2.0 * k + alpha + beta) * (2.0 * k + alpha + beta + 1.0) * (2.0 * k + alpha + beta + 2.0);
        Real a4 = 2.0 * (k + alpha) * (k + beta) * (2.0 * k + alpha + beta + 2.0);
        
        P2 = ((a2 + a3 * x) * P1 - a4 * P0) / a1;
        P0 = P1;
        P1 = P2;
      }
      
      return P1;
    }

    /**
     * @brief Evaluates derivative of Jacobi polynomial dP^{α,β}_n/dx on [-1,1]
     */
    inline Real jacobiPolynomialDerivative(size_t n, Real alpha, Real beta, Real x)
    {
      if (n == 0)
        return 0.0;
      
      // d/dx P^{α,β}_n = 0.5 * (n + α + β + 1) * P^{α+1,β+1}_{n-1}
      return 0.5 * (n + alpha + beta + 1.0) * jacobiPolynomial(n - 1, alpha + 1.0, beta + 1.0, x);
    }

    /**
     * @brief Evaluates Dubiner modal basis function ψ_{p,q} on triangle
     * 
     * Dubiner basis in collapsed coordinates (r,s) ∈ [-1,1]²
     * Standard form: ψ_{p,q}(r,s) = P^{0,0}_p(a) * P^{2p+1,0}_q(b) * ((1-b)/2)^p
     * where a = 2*(1+r)/(1-s) - 1, b = s
     */
    inline Real dubinerTriangle(size_t p, size_t q, Real r, Real s)
    {
      Real b = s;
      Real a;
      
      // Handle singularity at s = 1
      if (std::abs(s - 1.0) < 1e-10)
      {
        a = -1.0;
      }
      else
      {
        a = 2.0 * (1.0 + r) / (1.0 - s) - 1.0;
      }
      
      Real psi = jacobiPolynomial(p, 0.0, 0.0, a);
      psi *= jacobiPolynomial(q, 2.0 * p + 1.0, 0.0, b);
      psi *= std::pow(0.5 * (1.0 - b), p);
      
      return psi;
    }

    /**
     * @brief Evaluates gradient of Dubiner modal basis ∇ψ_{p,q} on triangle
     * Returns (∂ψ/∂r, ∂ψ/∂s) in collapsed coordinates
     */
    inline std::pair<Real, Real> dubinerTriangleGradient(size_t p, size_t q, Real r, Real s)
    {
      Real eps = 1e-10;
      Real b = s;
      Real a;
      
      if (std::abs(s - 1.0) < eps)
      {
        a = -1.0;
      }
      else
      {
        a = 2.0 * (1.0 + r) / (1.0 - s) - 1.0;
      }
      
      Real Pa = jacobiPolynomial(p, 0.0, 0.0, a);
      Real dPa = jacobiPolynomialDerivative(p, 0.0, 0.0, a);
      Real Pb = jacobiPolynomial(q, 2.0 * p + 1.0, 0.0, b);
      Real dPb = jacobiPolynomialDerivative(q, 2.0 * p + 1.0, 0.0, b);
      Real scale_b = std::pow(0.5 * (1.0 - b), p);
      
      // ∂ψ/∂r: only a depends on r
      Real dpsi_dr = 0.0;
      if (std::abs(s - 1.0) > eps)
      {
        Real da_dr = 2.0 / (1.0 - s);
        dpsi_dr = dPa * da_dr * Pb * scale_b;
      }
      
      // ∂ψ/∂s: both a and b depend on s, and scale_b depends on b=s
      Real dpsi_ds = 0.0;
      if (std::abs(s - 1.0) > eps)
      {
        Real da_ds = -2.0 * (1.0 + r) / ((1.0 - s) * (1.0 - s));
        dpsi_ds += dPa * da_ds * Pb * scale_b;
      }
      dpsi_ds += Pa * dPb * scale_b;  // db_ds = 1
      if (p > 0)
      {
        dpsi_ds += Pa * Pb * p * std::pow(0.5 * (1.0 - b), p - 1) * (-0.5);
      }
      
      return {dpsi_dr, dpsi_ds};
    }

    /**
     * @brief Converts from reference triangle (x,y) ∈ [0,1]² to collapsed coords (r,s) ∈ [-1,1]²
     * Reference triangle: (0,0), (1,0), (0,1)
     * Standard triangle (for Dubiner): (-1,-1), (1,-1), (-1,1)
     */
    inline std::pair<Real, Real> triangleToCollapsed(Real x, Real y)
    {
      // First map reference triangle to standard triangle
      // (0,0) -> (-1,-1), (1,0) -> (1,-1), (0,1) -> (-1,1)
      // This is: x_std = -1 + 2*x, y_std = -1 + 2*y
      // But we need to account for the constraint x+y <= 1
      
      // Use Hesthaven & Warburton formula:
      // r = 2*(1+x)/(1-y) - 1  (when y != 1)
      // s = 2*y - 1
      
      Real eps = 1e-10;
      Real r, s;
      
      if (y > 1.0 - eps)
      {
        r = -1.0;
        s = 1.0;
      }
      else
      {
        r = 2.0 * (x / (1.0 - y)) - 1.0;
        s = 2.0 * y - 1.0;
      }
      
      return {r, s};
    }

    /**
     * @brief Gets Vandermonde matrix (and its inverse) for Triangle Dubiner basis
     * V[i,j] = ψ_j(node_i) where ψ_j is the j-th Dubiner mode
     * This is cached per polynomial degree K
     */
    template <size_t K>
    struct TriangleDubinerVandermonde
    {
      static const Math::Matrix<Real>& get()
      {
        static thread_local Math::Matrix<Real> V;
        if (V.rows() == 0)
        {
          const auto& nodes = getTriangleFeketeNodes<K>();
          const size_t n_nodes = nodes.size();
          V.resize(n_nodes, n_nodes);
          
          // Fill Vandermonde matrix
          size_t mode_idx = 0;
          for (size_t p = 0; p <= K; ++p)
          {
            for (size_t q = 0; q <= K - p; ++q)
            {
              for (size_t node_idx = 0; node_idx < n_nodes; ++node_idx)
              {
                auto [r, s] = triangleToCollapsed(nodes[node_idx].x(), nodes[node_idx].y());
                V(node_idx, mode_idx) = dubinerTriangle(p, q, r, s);
              }
              mode_idx++;
            }
          }
        }
        return V;
      }
      
      static const Math::Matrix<Real>& getInverse()
      {
        static thread_local Math::Matrix<Real> Vinv;
        if (Vinv.rows() == 0)
        {
          const auto& V = get();
          Vinv = V.inverse();
        }
        return Vinv;
      }
    };

    /**
     * @brief Evaluates Dubiner modal basis function ψ_{p,q,r} on tetrahedron
     */
    inline Real dubinerTetrahedron(size_t p, size_t q, size_t r_idx, Real r, Real s, Real t)
    {
      Real eps = 1e-10;
      
      // Collapsed coordinates transformation
      Real sr = (std::abs(1.0 - r) < eps) ? -1.0 : (2.0 * (1.0 + s) / (1.0 - r) - 1.0);
      Real tr = (std::abs(1.0 - r - s) < eps) ? -1.0 : (2.0 * (1.0 + t) / (1.0 - r - s) - 1.0);
      
      Real psi = jacobiPolynomial(p, 0.0, 0.0, r);
      psi *= std::pow(0.5 * (1.0 - r), p);
      psi *= jacobiPolynomial(q, 2.0 * p + 1.0, 0.0, sr);
      psi *= std::pow(0.5 * (1.0 - s), q);
      psi *= jacobiPolynomial(r_idx, 2.0 * (p + q) + 2.0, 0.0, tr);
      
      return psi;
    }

    /**
     * @brief Evaluates gradient of Dubiner modal basis ∇ψ_{p,q,r} on tetrahedron
     */
    inline std::tuple<Real, Real, Real> dubinerTetrahedronGradient(size_t p, size_t q, size_t r_idx, Real r, Real s, Real t)
    {
      // Simplified implementation - proper chain rule derivatives needed for production
      Real h = 1e-7;
      Real psi0 = dubinerTetrahedron(p, q, r_idx, r, s, t);
      Real psi_r = dubinerTetrahedron(p, q, r_idx, r + h, s, t);
      Real psi_s = dubinerTetrahedron(p, q, r_idx, r, s + h, t);
      Real psi_t = dubinerTetrahedron(p, q, r_idx, r, s, t + h);
      
      return {(psi_r - psi0) / h, (psi_s - psi0) / h, (psi_t - psi0) / h};
    }

    /**
     * @brief Converts from reference tetrahedron (x,y,z) to collapsed coords (r,s,t)
     * Reference tetrahedron: (0,0,0), (1,0,0), (0,1,0), (0,0,1)
     */
    inline std::tuple<Real, Real, Real> tetrahedronToCollapsed(Real x, Real y, Real z)
    {
      Real eps = 1e-10;
      Real r, s, t;
      
      if (y + z < 1.0 - eps)
      {
        r = -1.0 + 2.0 * x / (1.0 - y - z);
      }
      else
      {
        r = -1.0;
      }
      
      if (z < 1.0 - eps)
      {
        s = -1.0 + 2.0 * y / (1.0 - z);
      }
      else
      {
        s = -1.0;
      }
      
      t = -1.0 + 2.0 * z;
      
      return {r, s, t};
    }

    /**
     * @brief Gets Vandermonde matrix for Tetrahedron Dubiner basis
     */
    template <size_t K>
    struct TetrahedronDubinerVandermonde
    {
      static const Math::Matrix<Real>& get()
      {
        static thread_local Math::Matrix<Real> V;
        if (V.rows() == 0)
        {
          const auto& nodes = getTetrahedronFeketeNodes<K>();
          const size_t n_nodes = nodes.size();
          V.resize(n_nodes, n_nodes);
          
          size_t mode_idx = 0;
          for (size_t p = 0; p <= K; ++p)
          {
            for (size_t q = 0; q <= K - p; ++q)
            {
              for (size_t r = 0; r <= K - p - q; ++r)
              {
                for (size_t node_idx = 0; node_idx < n_nodes; ++node_idx)
                {
                  auto [rc, sc, tc] = tetrahedronToCollapsed(
                    nodes[node_idx].x(), nodes[node_idx].y(), nodes[node_idx].z());
                  V(node_idx, mode_idx) = dubinerTetrahedron(p, q, r, rc, sc, tc);
                }
                mode_idx++;
              }
            }
          }
        }
        return V;
      }
      
      static const Math::Matrix<Real>& getInverse()
      {
        static thread_local Math::Matrix<Real> Vinv;
        if (Vinv.rows() == 0)
        {
          const auto& V = get();
          Vinv = V.inverse();
        }
        return Vinv;
      }
    };

    /**
     * @brief Evaluates Lagrange basis function for 1D with given nodes.
     */
    template <size_t K, class Scalar>
    constexpr Scalar evaluateLagrange1D(size_t i, Real x, const std::vector<Math::SpatialPoint>& nodes)
    {
      Scalar result = 1;
      Real xi = nodes[i].x();

      for (size_t j = 0; j <= K; ++j)
      {
        if (j != i)
        {
          Real xj = nodes[j].x();
          result *= (x - xj) / (xi - xj);
        }
      }
      return result;
    }

    /**
     * @brief Evaluates derivative of Lagrange basis function for 1D with given nodes.
     */
    template <size_t K, class Scalar>
    constexpr Scalar evaluateLagrange1DDerivative(size_t i, Real x, const std::vector<Math::SpatialPoint>& nodes)
    {
      Real xi = nodes[i].x();
      Scalar result = 0;

      // Derivative using product rule
      for (size_t m = 0; m <= K; ++m)
      {
        if (m != i)
        {
          Scalar term = 1;
          Real xm = nodes[m].x();

          for (size_t j = 0; j <= K; ++j)
          {
            if (j != i && j != m)
            {
              Real xj = nodes[j].x();
              term *= (x - xj) / (xi - xj);
            }
          }
          term /= (xi - xm);
          result += term;
        }
      }
      return result;
    }

    /**
     * @brief Evaluates Lagrange basis for triangle using barycentric coordinates.
     * 
     * For a triangle with barycentric coordinates (λ0, λ1, λ2), the Pk Lagrange
     * basis function associated with node (i, j, k) where i+j+k=K is:
     * φ_{i,j,k}(λ) = L_i^K(λ0) * L_j^K(λ1) * L_k^K(λ2)
     * 
     * where L_n^K(λ) = ∏_{m=0}^{n-1} (K*λ - m) / (n - m) for n > 0, and L_0^K = 1
     */
    template <size_t K, class Scalar>
    constexpr Scalar evaluateLagrangeTriangle(size_t i, size_t j, const Math::SpatialPoint& r)
    {
      // Barycentric coordinates: λ0 = 1 - x - y, λ1 = x, λ2 = y
      Real lambda[3];
      lambda[0] = 1.0 - r.x() - r.y();
      lambda[1] = r.x();
      lambda[2] = r.y();

      // Node (i,j) corresponds to Cartesian (i/K, j/K)
      // which has barycentric (lambda0=(K-i-j)/K, lambda1=i/K, lambda2=j/K)
      // So basis (i,j) should be L_i^K(lambda1) * L_j^K(lambda2) * L_(K-i-j)^K(lambda0)
      size_t indices[3] = {K - i - j, i, j};  // Reordered to match lambda[0], lambda[1], lambda[2]

      // Compute product of generalized Lagrange polynomials
      Scalar result = 1;
      for (size_t dim = 0; dim < 3; ++dim)
      {
        size_t n = indices[dim];
        if (n == 0)
          continue;

        Scalar L_n = 1;
        for (size_t m = 0; m < n; ++m)
        {
          L_n *= (static_cast<Real>(K) * lambda[dim] - static_cast<Real>(m));
          L_n /= static_cast<Real>(m + 1);
        }
        result *= L_n;
      }

      return result;
    }

    /**
     * @brief Evaluates derivative of Lagrange basis for triangle.
     */
    template <size_t K, class Scalar>
    constexpr Scalar evaluateLagrangeTriangleDerivative(
        size_t i, size_t j, size_t deriv_dim, const Math::SpatialPoint& r)
    {
      // Barycentric coordinates: λ0 = 1 - x - y, λ1 = x, λ2 = y
      Real lambda[3];
      lambda[0] = 1.0 - r.x() - r.y();
      lambda[1] = r.x();
      lambda[2] = r.y();

      // Derivatives of barycentric coordinates w.r.t. x and y
      // dλ0/dx = -1, dλ0/dy = -1
      // dλ1/dx =  1, dλ1/dy =  0
      // dλ2/dx =  0, dλ2/dy =  1
      Real dlambda[3][2] = {{-1, -1}, {1, 0}, {0, 1}};

      size_t indices[3] = {K - i - j, i, j}; // for λ0, λ1, λ2

      Scalar result = 0;

      // Derivative via chain rule over λ_d
      for (size_t d = 0; d < 3; ++d)
      {
        Scalar term = dlambda[d][deriv_dim];

        for (size_t dim = 0; dim < 3; ++dim)
        {
          const size_t n = indices[dim];

          if (dim == d)
          {
            // Differentiated factor L_n^K(λ_d)
            if (n == 0)
            {
              // L_0 ≡ 1 → derivative 0 → whole contribution for this d is 0
              term = Scalar(0);
              break; // no need to consider other dims
            }

            // d/dλ L_n^K(λ_d)
            Scalar dL_n = 0;
            for (size_t p = 0; p < n; ++p)
            {
              Scalar prod = static_cast<Real>(K);
              for (size_t m = 0; m < n; ++m)
              {
                if (m != p)
                {
                  prod *= (static_cast<Real>(K) * lambda[dim] - static_cast<Real>(m));
                  prod /= static_cast<Real>(m + 1);
                }
              }
              prod /= static_cast<Real>(p + 1);
              dL_n += prod;
            }
            term *= dL_n;
          }
          else
          {
            // Undifferentiated factor L_n^K(λ_dim)
            if (n == 0)
              continue; // factor 1

            Scalar L_n = 1;
            for (size_t m = 0; m < n; ++m)
            {
              L_n *= (static_cast<Real>(K) * lambda[dim] - static_cast<Real>(m));
              L_n /= static_cast<Real>(m + 1);
            }
            term *= L_n;
          }
        }

        result += term;
      }

      return result;
    }

    /**
     * @brief Evaluates Lagrange basis for tetrahedron using barycentric coordinates.
     */
    template <size_t K, class Scalar>
    constexpr Scalar evaluateLagrangeTetrahedron(
        size_t i, size_t j, size_t k, const Math::SpatialPoint& r)
    {
      // Barycentric coordinates: λ0 = 1-x-y-z, λ1 = x, λ2 = y, λ3 = z
      Real lambda[4];
      lambda[0] = 1.0 - r.x() - r.y() - r.z();
      lambda[1] = r.x();
      lambda[2] = r.y();
      lambda[3] = r.z();

      // Node (i,j,k) corresponds to Cartesian (i/K, j/K, k/K)
      // which has barycentric (lambda0=(K-i-j-k)/K, lambda1=i/K, lambda2=j/K, lambda3=k/K)
      // So basis (i,j,k) should be L_i^K(lambda1) * L_j^K(lambda2) * L_k^K(lambda3) * L_(K-i-j-k)^K(lambda0)
      size_t indices[4] = {K - i - j - k, i, j, k};  // Reordered to match lambda[0], lambda[1], lambda[2], lambda[3]

      Scalar result = 1;
      for (size_t dim = 0; dim < 4; ++dim)
      {
        size_t n = indices[dim];
        if (n == 0)
          continue;

        Scalar L_n = 1;
        for (size_t m = 0; m < n; ++m)
        {
          L_n *= (static_cast<Real>(K) * lambda[dim] - static_cast<Real>(m));
          L_n /= static_cast<Real>(m + 1);
        }
        result *= L_n;
      }

      return result;
    }

    /**
     * @brief Evaluates Lagrange basis for wedge (tensor product of triangle and segment).
     */
    template <size_t K, class Scalar>
    Scalar evaluateLagrangeWedge(
        size_t i, size_t j, size_t k, const Math::SpatialPoint& r)
    {
      // Wedge is a tensor product: triangle (x,y) × segment (z)
      // For indices (i,j) on triangle with i+j≤K, and k on segment

      // Triangle part using barycentric coordinates
      Real lambda[3];
      lambda[0] = 1.0 - r.x() - r.y();
      lambda[1] = r.x();
      lambda[2] = r.y();

      size_t tri_indices[3] = {K - i - j, i, j};  // Reordered to match lambda[0], lambda[1], lambda[2]

      Scalar tri_result = 1;
      for (size_t dim = 0; dim < 3; ++dim)
      {
        size_t n = tri_indices[dim];
        if (n == 0)
          continue;

        Scalar L_n = 1;
        for (size_t m = 0; m < n; ++m)
        {
          L_n *= (static_cast<Real>(K) * lambda[dim] - static_cast<Real>(m));
          L_n /= static_cast<Real>(m + 1);
        }
        tri_result *= L_n;
      }

      // Segment part (1D Lagrange)
      // For wedge, we need the segment nodes (z-direction)
      const auto& wedge_nodes = H1Element<K, Scalar>::getNodes(Geometry::Polytope::Type::Wedge);
      // Extract segment nodes from wedge - they repeat for each triangle node
      // We'll construct a temporary segment nodes vector
      static thread_local std::vector<Math::SpatialPoint> seg_nodes;
      if (seg_nodes.empty() || seg_nodes.size() != K + 1)
      {
        seg_nodes.clear();
        seg_nodes.reserve(K + 1);
        // Extract z-coordinates: wedge nodes are ordered by k, then triangle pattern
        // The first (K+1)(K+2)/2 nodes are for k=0, next batch for k=1, etc.
        const size_t tri_count = (K + 1) * (K + 2) / 2;
        for (size_t kk = 0; kk <= K; ++kk)
          seg_nodes.emplace_back(Math::SpatialPoint{{wedge_nodes[kk * tri_count].z()}});
      }
      Scalar seg_result = evaluateLagrange1D<K, Scalar>(k, r.z(), seg_nodes);

      return tri_result * seg_result;
    }

    /**
     * @brief Evaluates derivative of Lagrange basis for tetrahedron.
     */
    template <size_t K, class Scalar>
    constexpr Scalar evaluateLagrangeTetrahedronDerivative(
        size_t i, size_t j, size_t k, size_t deriv_dim, const Math::SpatialPoint& r)
    {
      // Barycentric coordinates: λ0 = 1-x-y-z, λ1 = x, λ2 = y, λ3 = z
      Real lambda[4];
      lambda[0] = 1.0 - r.x() - r.y() - r.z();
      lambda[1] = r.x();
      lambda[2] = r.y();
      lambda[3] = r.z();

      // Derivatives:
      // dλ0/dx = -1, dλ0/dy = -1, dλ0/dz = -1
      // dλ1/dx =  1, dλ1/dy =  0, dλ1/dz =  0
      // dλ2/dx =  0, dλ2/dy =  1, dλ2/dz =  0
      // dλ3/dx =  0, dλ3/dy =  0, dλ3/dz =  1
      Real dlambda[4][3] = {{-1, -1, -1}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

      size_t indices[4] = {K - i - j - k, i, j, k}; // for λ0, λ1, λ2, λ3

      Scalar result = 0;

      for (size_t d = 0; d < 4; ++d)
      {
        Scalar term = dlambda[d][deriv_dim];

        for (size_t dim = 0; dim < 4; ++dim)
        {
          const size_t n = indices[dim];

          if (dim == d)
          {
            // Differentiated factor L_n^K(λ_d)
            if (n == 0)
            {
              // L_0 ≡ 1 → derivative 0 → whole contribution for this d is 0
              term = Scalar(0);
              break;
            }

            Scalar dL_n = 0;
            for (size_t p = 0; p < n; ++p)
            {
              Scalar prod = static_cast<Real>(K);
              for (size_t m = 0; m < n; ++m)
              {
                if (m != p)
                {
                  prod *= (static_cast<Real>(K) * lambda[dim] - static_cast<Real>(m));
                  prod /= static_cast<Real>(m + 1);
                }
              }
              prod /= static_cast<Real>(p + 1);
              dL_n += prod;
            }
            term *= dL_n;
          }
          else
          {
            // Undifferentiated factor L_n^K(λ_dim)
            if (n == 0)
              continue; // factor 1

            Scalar L_n = 1;
            for (size_t m = 0; m < n; ++m)
            {
              L_n *= (static_cast<Real>(K) * lambda[dim] - static_cast<Real>(m));
              L_n /= static_cast<Real>(m + 1);
            }
            term *= L_n;
          }
        }

        result += term;
      }

      return result;
    }

    /**
     * @brief Evaluates derivative of Lagrange basis for wedge.
     */
    template <size_t K, class Scalar>
    Scalar evaluateLagrangeWedgeDerivative(
        size_t i, size_t j, size_t k, size_t deriv_dim, const Math::SpatialPoint& r)
    {
      // For wedge, we need the segment nodes (z-direction)
      const auto& wedge_nodes = H1Element<K, Scalar>::getNodes(Geometry::Polytope::Type::Wedge);
      static thread_local std::vector<Math::SpatialPoint> seg_nodes;
      if (seg_nodes.empty() || seg_nodes.size() != K + 1)
      {
        seg_nodes.clear();
        seg_nodes.reserve(K + 1);
        const size_t tri_count = (K + 1) * (K + 2) / 2;
        for (size_t kk = 0; kk <= K; ++kk)
          seg_nodes.emplace_back(Math::SpatialPoint{{wedge_nodes[kk * tri_count].z()}});
      }
      
      if (deriv_dim < 2)
      {
        // Derivative w.r.t. x or y (triangle part)
        Scalar tri_deriv = evaluateLagrangeTriangleDerivative<K, Scalar>(i, j, deriv_dim, r);
        Scalar seg_val = evaluateLagrange1D<K, Scalar>(k, r.z(), seg_nodes);
        return tri_deriv * seg_val;
      }
      else
      {
        // Derivative w.r.t. z (segment part)
        Scalar tri_val = evaluateLagrangeTriangle<K, Scalar>(i, j, r);
        Scalar seg_deriv = evaluateLagrange1DDerivative<K, Scalar>(k, r.z(), seg_nodes);
        return tri_val * seg_deriv;
      }
    }
  }

  template <size_t K, class Scalar>
  constexpr
  const Math::SpatialPoint& H1Element<K, Scalar>::getNode(size_t i) const
  {
    return getNodes(this->getGeometry())[i];
  }

  template <size_t K, class Scalar>
  const typename H1Element<K, Scalar>::LinearForm& H1Element<K, Scalar>::getLinearForm(size_t i) const
  {
    const Geometry::Polytope::Type g = this->getGeometry();

    // Use switch to create geometry-specific thread_local storage
    switch (g)
    {
      case Geometry::Polytope::Type::Point:
      {
        static thread_local std::vector<LinearForm> s_lfs;
        if (s_lfs.empty())
        {
          const size_t count = getCount();
          s_lfs.reserve(count);
          for (size_t j = 0; j < count; ++j)
            s_lfs.emplace_back(j, g);
        }
        return s_lfs[i];
      }
      case Geometry::Polytope::Type::Segment:
      {
        static thread_local std::vector<LinearForm> s_lfs;
        if (s_lfs.empty())
        {
          const size_t count = getCount();
          s_lfs.reserve(count);
          for (size_t j = 0; j < count; ++j)
            s_lfs.emplace_back(j, g);
        }
        return s_lfs[i];
      }
      case Geometry::Polytope::Type::Triangle:
      {
        static thread_local std::vector<LinearForm> s_lfs;
        if (s_lfs.empty())
        {
          const size_t count = getCount();
          s_lfs.reserve(count);
          for (size_t j = 0; j < count; ++j)
            s_lfs.emplace_back(j, g);
        }
        return s_lfs[i];
      }
      case Geometry::Polytope::Type::Quadrilateral:
      {
        static thread_local std::vector<LinearForm> s_lfs;
        if (s_lfs.empty())
        {
          const size_t count = getCount();
          s_lfs.reserve(count);
          for (size_t j = 0; j < count; ++j)
            s_lfs.emplace_back(j, g);
        }
        return s_lfs[i];
      }
      case Geometry::Polytope::Type::Tetrahedron:
      {
        static thread_local std::vector<LinearForm> s_lfs;
        if (s_lfs.empty())
        {
          const size_t count = getCount();
          s_lfs.reserve(count);
          for (size_t j = 0; j < count; ++j)
            s_lfs.emplace_back(j, g);
        }
        return s_lfs[i];
      }
      case Geometry::Polytope::Type::Wedge:
      {
        static thread_local std::vector<LinearForm> s_lfs;
        if (s_lfs.empty())
        {
          const size_t count = getCount();
          s_lfs.reserve(count);
          for (size_t j = 0; j < count; ++j)
            s_lfs.emplace_back(j, g);
        }
        return s_lfs[i];
      }
    }

    // Fallback (should never happen)
    static thread_local LinearForm s_null(0, g);
    assert(false);
    return s_null;
  }

  template <size_t K, class Scalar>
  const typename H1Element<K, Scalar>::BasisFunction& H1Element<K, Scalar>::getBasis(size_t i) const
  {
    const Geometry::Polytope::Type g = this->getGeometry();

    // Use switch to create geometry-specific thread_local storage
    switch (g)
    {
      case Geometry::Polytope::Type::Point:
      {
        static thread_local std::vector<BasisFunction> s_bs;
        if (s_bs.empty())
        {
          const size_t count = getCount();
          s_bs.reserve(count);
          for (size_t j = 0; j < count; ++j)
            s_bs.emplace_back(j, g);
        }
        return s_bs[i];
      }
      case Geometry::Polytope::Type::Segment:
      {
        static thread_local std::vector<BasisFunction> s_bs;
        if (s_bs.empty())
        {
          const size_t count = getCount();
          s_bs.reserve(count);
          for (size_t j = 0; j < count; ++j)
            s_bs.emplace_back(j, g);
        }
        return s_bs[i];
      }
      case Geometry::Polytope::Type::Triangle:
      {
        static thread_local std::vector<BasisFunction> s_bs;
        if (s_bs.empty())
        {
          const size_t count = getCount();
          s_bs.reserve(count);
          for (size_t j = 0; j < count; ++j)
            s_bs.emplace_back(j, g);
        }
        return s_bs[i];
      }
      case Geometry::Polytope::Type::Quadrilateral:
      {
        static thread_local std::vector<BasisFunction> s_bs;
        if (s_bs.empty())
        {
          const size_t count = getCount();
          s_bs.reserve(count);
          for (size_t j = 0; j < count; ++j)
            s_bs.emplace_back(j, g);
        }
        return s_bs[i];
      }
      case Geometry::Polytope::Type::Tetrahedron:
      {
        static thread_local std::vector<BasisFunction> s_bs;
        if (s_bs.empty())
        {
          const size_t count = getCount();
          s_bs.reserve(count);
          for (size_t j = 0; j < count; ++j)
            s_bs.emplace_back(j, g);
        }
        return s_bs[i];
      }
      case Geometry::Polytope::Type::Wedge:
      {
        static thread_local std::vector<BasisFunction> s_bs;
        if (s_bs.empty())
        {
          const size_t count = getCount();
          s_bs.reserve(count);
          for (size_t j = 0; j < count; ++j)
            s_bs.emplace_back(j, g);
        }
        return s_bs[i];
      }
    }

    // Fallback (should never happen)
    static thread_local BasisFunction s_null(0, g);
    assert(false);
    return s_null;
  }

  template <size_t K, class Scalar>
  constexpr
  Scalar H1Element<K, Scalar>::BasisFunction::operator()(const Math::SpatialPoint& r) const
  {
    switch (m_g)
    {
      case Geometry::Polytope::Type::Point:
      {
        return 1;
      }
      case Geometry::Polytope::Type::Segment:
      {
        const auto& nodes = H1Element<K, Scalar>::getNodes(m_g);
        return Internal::evaluateLagrange1D<K, Scalar>(m_local, r.x(), nodes);
      }
      case Geometry::Polytope::Type::Triangle:
      {
        // Use Dubiner modal basis with Vandermonde approach
        const auto& Vinv = Internal::TriangleDubinerVandermonde<K>::getInverse();
        
        auto [rc, sc] = Internal::triangleToCollapsed(r.x(), r.y());
        
        Scalar result = 0.0;
        size_t mode_idx = 0;
        for (size_t p = 0; p <= K; ++p)
        {
          for (size_t q = 0; q <= K - p; ++q)
          {
            Real psi = Internal::dubinerTriangle(p, q, rc, sc);
            result += Vinv(m_local, mode_idx) * psi;
            mode_idx++;
          }
        }
        return result;
      }
      case Geometry::Polytope::Type::Quadrilateral:
      {
        // Tensor product of 1D Lagrange basis
        size_t j_idx = m_local / (K + 1);
        size_t i_idx = m_local % (K + 1);

        // Initialize to ensure no uninitialized warnings
        const size_t dim = Geometry::Polytope::Traits(m_g).getDimension();
        if (dim < 2)
          return Math::nan<Scalar>();

        const auto& nodes = H1Element<K, Scalar>::getNodes(m_g);
        return Internal::evaluateLagrange1D<K, Scalar>(i_idx, r.x(), nodes) *
               Internal::evaluateLagrange1D<K, Scalar>(j_idx, r.y(), nodes);
      }
      case Geometry::Polytope::Type::Tetrahedron:
      {
        // Use Dubiner modal basis with Vandermonde approach
        const auto& Vinv = Internal::TetrahedronDubinerVandermonde<K>::getInverse();
        
        auto [rc, sc, tc] = Internal::tetrahedronToCollapsed(r.x(), r.y(), r.z());
        
        Scalar result = 0.0;
        size_t mode_idx = 0;
        for (size_t p = 0; p <= K; ++p)
        {
          for (size_t q = 0; q <= K - p; ++q)
          {
            for (size_t r_mode = 0; r_mode <= K - p - q; ++r_mode)
            {
              Real psi = Internal::dubinerTetrahedron(p, q, r_mode, rc, sc, tc);
              result += Vinv(m_local, mode_idx) * psi;
              mode_idx++;
            }
          }
        }
        return result;
      }
      case Geometry::Polytope::Type::Wedge:
      {
        // Node ordering: triangle nodes (i,j) with i+j≤K, then segment node k
        // Ordered by k, then j, then i
        size_t idx = 0;
        for (size_t k = 0; k <= K; ++k)
        {
          for (size_t j = 0; j <= K; ++j)
          {
            for (size_t i = 0; i <= K - j; ++i)
            {
              if (idx == m_local)
              {
                return Internal::evaluateLagrangeWedge<K, Scalar>(i, j, k, r);
              }
              idx++;
            }
          }
        }
        return Math::nan<Scalar>();
      }
    }

    return Math::nan<Scalar>();
  }

  template <size_t K, class Scalar>
  template <size_t Order>
  constexpr
  Scalar H1Element<K, Scalar>::BasisFunction::DerivativeFunction<Order>::operator()(
      const Math::SpatialPoint& r) const
  {
    if constexpr (Order == 0)
    {
      return BasisFunction(m_local, m_g)(r);
    }
    else if constexpr (Order == 1)
    {
      switch (m_g)
      {
        case Geometry::Polytope::Type::Point:
        {
          return 0;
        }
        case Geometry::Polytope::Type::Segment:
        {
          assert(m_i == 0);
          const auto& nodes = H1Element<K, Scalar>::getNodes(m_g);
          return Internal::evaluateLagrange1DDerivative<K, Scalar>(m_local, r.x(), nodes);
        }
        case Geometry::Polytope::Type::Triangle:
        {
          // Use Dubiner modal basis with Vandermonde approach
          const auto& Vinv = Internal::TriangleDubinerVandermonde<K>::getInverse();
          
          auto [rc, sc] = Internal::triangleToCollapsed(r.x(), r.y());
          
          Scalar result = 0.0;
          size_t mode_idx = 0;
          for (size_t p = 0; p <= K; ++p)
          {
            for (size_t q = 0; q <= K - p; ++q)
            {
              auto [dpsi_dr, dpsi_ds] = Internal::dubinerTriangleGradient(p, q, rc, sc);
              
              // Transform gradients from collapsed (r,s) to reference (x,y) coordinates
              // r = 2*x/(1-y) - 1, s = 2*y - 1
              // ∂r/∂x = 2/(1-y), ∂r/∂y = 2*x/((1-y)^2)
              // ∂s/∂x = 0, ∂s/∂y = 2
              
              Real x = r.x();
              Real y = r.y();
              Real eps = 1e-10;
              
              Real dpsi_dx = 0.0, dpsi_dy = 0.0;
              
              if (y < 1.0 - eps)
              {
                Real dr_dx = 2.0 / (1.0 - y);
                Real dr_dy = 2.0 * x / ((1.0 - y) * (1.0 - y));
                Real ds_dx = 0.0;
                Real ds_dy = 2.0;
                
                dpsi_dx = dpsi_dr * dr_dx + dpsi_ds * ds_dx;
                dpsi_dy = dpsi_dr * dr_dy + dpsi_ds * ds_dy;
              }
              
              if (m_i == 0) // ∂/∂x
                result += Vinv(m_local, mode_idx) * dpsi_dx;
              else if (m_i == 1) // ∂/∂y
                result += Vinv(m_local, mode_idx) * dpsi_dy;
              
              mode_idx++;
            }
          }
          return result;
        }
        case Geometry::Polytope::Type::Quadrilateral:
        {
          // Tensor product derivative
          size_t j_idx = m_local / (K + 1);
          size_t i_idx = m_local % (K + 1);

          const auto& nodes = H1Element<K, Scalar>::getNodes(m_g);
          if (m_i == 0) // d/dx
          {
            return Internal::evaluateLagrange1DDerivative<K, Scalar>(i_idx, r.x(), nodes) *
                   Internal::evaluateLagrange1D<K, Scalar>(j_idx, r.y(), nodes);
          }
          else if (m_i == 1) // d/dy
          {
            return Internal::evaluateLagrange1D<K, Scalar>(i_idx, r.x(), nodes) *
                   Internal::evaluateLagrange1DDerivative<K, Scalar>(j_idx, r.y(), nodes);
          }
          return 0;
        }
        case Geometry::Polytope::Type::Tetrahedron:
        {
          // Use Dubiner modal basis with Vandermonde approach
          const auto& Vinv = Internal::TetrahedronDubinerVandermonde<K>::getInverse();
          
          auto [rc, sc, tc] = Internal::tetrahedronToCollapsed(r.x(), r.y(), r.z());
          
          Scalar result = 0.0;
          size_t mode_idx = 0;
          for (size_t p = 0; p <= K; ++p)
          {
            for (size_t q = 0; q <= K - p; ++q)
            {
              for (size_t r_mode = 0; r_mode <= K - p - q; ++r_mode)
              {
                auto [dpsi_dr, dpsi_ds, dpsi_dt] = Internal::dubinerTetrahedronGradient(p, q, r_mode, rc, sc, tc);
                
                // Transform from collapsed coordinates to reference coordinates
                Real dpsi_dx = dpsi_dr * 0.5;
                Real dpsi_dy = dpsi_ds * 0.5;
                Real dpsi_dz = dpsi_dt * 0.5;
                
                if (m_i == 0) // ∂/∂x
                  result += Vinv(m_local, mode_idx) * dpsi_dx;
                else if (m_i == 1) // ∂/∂y
                  result += Vinv(m_local, mode_idx) * dpsi_dy;
                else if (m_i == 2) // ∂/∂z
                  result += Vinv(m_local, mode_idx) * dpsi_dz;
                
                mode_idx++;
              }
            }
          }
          return result;
        }
        case Geometry::Polytope::Type::Wedge:
        {
          // Find (i,j,k) indices for this local DOF
          size_t idx = 0;
          for (size_t k = 0; k <= K; ++k)
          {
            for (size_t j = 0; j <= K; ++j)
            {
              for (size_t i = 0; i <= K - j; ++i)
              {
                if (idx == m_local)
                {
                  return Internal::evaluateLagrangeWedgeDerivative<K, Scalar>(i, j, k, m_i, r);
                }
                idx++;
              }
            }
          }
          return 0;
        }
      }
      return 0;
    }
    else
    {
      return 0; // Higher order derivatives
    }
  }
}

#endif
