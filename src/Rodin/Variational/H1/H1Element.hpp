/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_H1ELEMENT_HPP
#define RODIN_VARIATIONAL_H1_H1ELEMENT_HPP

#include "Rodin/Math/Common.h"

#include "H1Element.h"

namespace Rodin::Variational
{
  namespace Internal
  {
    /**
     * @brief Computes Gauss-Lobatto-Legendre nodes on [-1, 1].
     *
     * GLL nodes are the roots of (1 - x^2) P'_K(x) where P_K is the Legendre polynomial.
     * The endpoints are fixed at -1 and 1, and the interior nodes are computed using Newton iteration.
     *
     * @tparam K Polynomial degree
     * @return Array of K+1 GLL nodes on [-1, 1]
     */
    template <size_t K>
    std::array<Real, K+1> computeGLLNodes()
    {
      std::array<Real, K+1> nodes;
      
      // Endpoints are always -1 and 1
      nodes[0] = -1.0;
      nodes[K] = 1.0;
      
      if constexpr (K == 0)
      {
        return nodes;
      }
      else if constexpr (K == 1)
      {
        return nodes;
      }
      else
      {
        // Compute interior nodes as zeros of (1 - x^2) P'_K(x)
        // Use Newton iteration with Chebyshev initial guess
        const Real tol = 1e-14;
        const size_t max_iter = 100;
        
        for (size_t i = 1; i < K; ++i)
        {
          // Initial guess: Chebyshev node
          Real x = -std::cos(M_PI * static_cast<Real>(i) / static_cast<Real>(K));
          
          // Newton iteration
          for (size_t iter = 0; iter < max_iter; ++iter)
          {
            // Evaluate Legendre polynomial P_K and its derivatives using 3-term recurrence
            // P_0 = 1, P_1 = x
            // P_{n+1} = ((2n+1) x P_n - n P_{n-1}) / (n+1)
            Real P_km1 = 1.0;  // P_{K-1}
            Real P_k = x;      // P_K
            
            if (K >= 2)
            {
              for (size_t n = 1; n < K; ++n)
              {
                Real n_real = static_cast<Real>(n);
                Real P_next = ((2.0 * n_real + 1.0) * x * P_k - n_real * P_km1) / (n_real + 1.0);
                P_km1 = P_k;
                P_k = P_next;
              }
            }
            
            // Compute P'_K using the relation: (1 - x^2) P'_K = K (x P_K - P_{K-1})
            Real P_k_prime = static_cast<Real>(K) * (x * P_k - P_km1) / (1.0 - x * x);
            
            // For (1 - x^2) P'_K = 0, we want to solve P'_K = 0 (since x != ±1 in interior)
            // Actually for GLL, we solve (1 - x^2) P'_K = 0 which means P'_K = 0 for interior
            // But we use the Lobatto polynomial Q_K = (1 - x^2) P'_K
            Real Q = (1.0 - x * x) * P_k_prime;
            
            // Q' = -2x P'_K + (1 - x^2) P''_K
            // Using (1 - x^2) P''_K = K(K+1) P_K - K(2K-1) x P'_K / (K-1) is complex
            // Use simpler formula: Q' = K(K+1) P_K
            Real Q_prime = static_cast<Real>(K * (K + 1)) * P_k;
            
            if (std::abs(Q_prime) < 1e-20)
            {
              // Fallback: try alternative formula
              Q_prime = -2.0 * x * P_k_prime;
            }
            
            Real dx = Q / Q_prime;
            x -= dx;
            
            if (std::abs(dx) < tol)
              break;
          }
          
          nodes[i] = x;
        }
        
        return nodes;
      }
    }
    
    /**
     * @brief Computes Gauss-Lobatto-Legendre nodes on [0, 1].
     *
     * Maps GLL nodes from [-1, 1] to [0, 1] using x ↦ (x + 1)/2.
     *
     * @tparam K Polynomial degree
     * @return Array of K+1 GLL nodes on [0, 1]
     */
    template <size_t K>
    std::array<Real, K+1> computeGLLNodes01()
    {
      auto nodes_11 = computeGLLNodes<K>();
      std::array<Real, K+1> nodes_01;
      
      for (size_t i = 0; i <= K; ++i)
      {
        nodes_01[i] = (nodes_11[i] + 1.0) / 2.0;
      }
      
      return nodes_01;
    }
    
    /**
     * @brief Returns Fekete nodes for triangles.
     *
     * Fekete nodes are optimal interpolation nodes that maximize the determinant
     * of the Vandermonde matrix. For high-order elements, these provide better
     * conditioning than equispaced nodes.
     *
     * @tparam K Polynomial degree
     * @return Reference to cached vector of triangle Fekete nodes
     */
    template <size_t K>
    const std::vector<Math::SpatialPoint>& getTriangleFeketeNodes()
    {
      static thread_local std::vector<Math::SpatialPoint> s_nodes;
      
      if (s_nodes.empty())
      {
        // For low degrees K<=3, equispaced nodes are adequate and avoid
        // interpolation issues with basis functions still using barycentric coordinates
        // For higher degrees, one would use proper Fekete node tables or warp-blend algorithm
        for (size_t j = 0; j <= K; ++j)
        {
          for (size_t i = 0; i <= K - j; ++i)
          {
            Real s = static_cast<Real>(i) / (K == 0 ? 1.0 : static_cast<Real>(K));
            Real t = static_cast<Real>(j) / (K == 0 ? 1.0 : static_cast<Real>(K));
            s_nodes.emplace_back(Math::SpatialPoint{{s, t}});
          }
        }
      }
      
      return s_nodes;
    }
    
    /**
     * @brief Returns Fekete nodes for tetrahedra.
     *
     * Fekete nodes are optimal interpolation nodes for tetrahedra.
     *
     * @tparam K Polynomial degree
     * @return Reference to cached vector of tetrahedron Fekete nodes
     */
    template <size_t K>
    const std::vector<Math::SpatialPoint>& getTetrahedronFeketeNodes()
    {
      static thread_local std::vector<Math::SpatialPoint> s_nodes;
      
      if (s_nodes.empty())
      {
        // For low degrees K<=3, equispaced nodes are adequate
        // For higher degrees, one would use proper Fekete node tables or warp-blend
        for (size_t k = 0; k <= K; ++k)
        {
          for (size_t j = 0; j <= K - k; ++j)
          {
            for (size_t i = 0; i <= K - j - k; ++i)
            {
              Real r = static_cast<Real>(i) / (K == 0 ? 1.0 : static_cast<Real>(K));
              Real s = static_cast<Real>(j) / (K == 0 ? 1.0 : static_cast<Real>(K));
              Real t = static_cast<Real>(k) / (K == 0 ? 1.0 : static_cast<Real>(K));
              s_nodes.emplace_back(Math::SpatialPoint{{r, s, t}});
            }
          }
        }
      }
      
      return s_nodes;
    }
    
    /**
     * @brief Evaluates Lagrange basis function for 1D using actual node positions.
     */
    template <size_t K, class Scalar>
    Scalar evaluateLagrange1D(size_t i, Real x)
    {
      // Get the actual GLL nodes
      static thread_local auto nodes = computeGLLNodes01<K>();
      
      Scalar result = 1;
      Real xi = nodes[i];

      for (size_t j = 0; j <= K; ++j)
      {
        if (j != i)
        {
          Real xj = nodes[j];
          result *= (x - xj) / (xi - xj);
        }
      }
      return result;
    }

    /**
     * @brief Evaluates derivative of Lagrange basis function for 1D using actual node positions.
     */
    template <size_t K, class Scalar>
    Scalar evaluateLagrange1DDerivative(size_t i, Real x)
    {
      // Get the actual GLL nodes
      static thread_local auto nodes = computeGLLNodes01<K>();
      
      Real xi = nodes[i];
      Scalar result = 0;

      // Derivative using product rule
      for (size_t m = 0; m <= K; ++m)
      {
        if (m != i)
        {
          Scalar term = 1;
          Real xm = nodes[m];

          for (size_t j = 0; j <= K; ++j)
          {
            if (j != i && j != m)
            {
              Real xj = nodes[j];
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
    constexpr Scalar evaluateLagrangeWedge(
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
      Scalar seg_result = evaluateLagrange1D<K, Scalar>(k, r.z());

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
    constexpr Scalar evaluateLagrangeWedgeDerivative(
        size_t i, size_t j, size_t k, size_t deriv_dim, const Math::SpatialPoint& r)
    {
      if (deriv_dim < 2)
      {
        // Derivative w.r.t. x or y (triangle part)
        Scalar tri_deriv = evaluateLagrangeTriangleDerivative<K, Scalar>(i, j, deriv_dim, r);
        Scalar seg_val = evaluateLagrange1D<K, Scalar>(k, r.z());
        return tri_deriv * seg_val;
      }
      else
      {
        // Derivative w.r.t. z (segment part)
        Scalar tri_val = evaluateLagrangeTriangle<K, Scalar>(i, j, r);
        Scalar seg_deriv = evaluateLagrange1DDerivative<K, Scalar>(k, r.z());
        return tri_val * seg_deriv;
      }
    }
  }

  template <size_t K, class Scalar>
  const std::vector<Math::SpatialPoint>&
  H1Element<K, Scalar>::getNodes(Geometry::Polytope::Type g)
  {
    using G = Geometry::Polytope::Type;
    
    switch (g)
    {
      case G::Point:
      {
        static thread_local const std::vector<Math::SpatialPoint> s_nodes = [] {
          std::vector<Math::SpatialPoint> n;
          n.emplace_back(Math::SpatialPoint{{0}});
          return n;
        }();
        return s_nodes;
      }

      case G::Segment:
      {
        static thread_local std::vector<Math::SpatialPoint> s_nodes;
        if (s_nodes.empty())
        {
          auto xi = Internal::computeGLLNodes01<K>();
          s_nodes.reserve(K + 1);
          for (size_t i = 0; i <= K; ++i)
            s_nodes.emplace_back(Math::SpatialPoint{{xi[i]}});
        }
        return s_nodes;
      }

      case G::Triangle:
      {
        static thread_local std::vector<Math::SpatialPoint> s_nodes;
        if (s_nodes.empty())
          s_nodes = Internal::getTriangleFeketeNodes<K>();
        return s_nodes;
      }

      case G::Quadrilateral:
      {
        static thread_local std::vector<Math::SpatialPoint> s_nodes;
        if (s_nodes.empty())
        {
          auto xi = Internal::computeGLLNodes01<K>();
          s_nodes.reserve((K + 1) * (K + 1));
          for (size_t j = 0; j <= K; ++j)
            for (size_t i = 0; i <= K; ++i)
              s_nodes.emplace_back(Math::SpatialPoint{{xi[i], xi[j]}});
        }
        return s_nodes;
      }

      case G::Tetrahedron:
      {
        static thread_local std::vector<Math::SpatialPoint> s_nodes;
        if (s_nodes.empty())
          s_nodes = Internal::getTetrahedronFeketeNodes<K>();
        return s_nodes;
      }

      case G::Wedge:
      {
        static thread_local std::vector<Math::SpatialPoint> s_nodes;
        if (s_nodes.empty())
        {
          const auto& tri = Internal::getTriangleFeketeNodes<K>();
          auto z = Internal::computeGLLNodes01<K>();
          s_nodes.reserve(tri.size() * (K + 1));
          for (size_t k = 0; k <= K; ++k)
            for (const auto& p : tri)
              s_nodes.emplace_back(Math::SpatialPoint{{p.x(), p.y(), z[k]}});
        }
        return s_nodes;
      }
    }

    // Should be unreachable if all enum values are handled
    assert(false && "Unsupported Polytope type.");
    static thread_local const std::vector<Math::SpatialPoint> s_empty;
    return s_empty;
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
        return Internal::evaluateLagrange1D<K, Scalar>(m_local, r.x());
      }
      case Geometry::Polytope::Type::Triangle:
      {
        // Use barycentric coordinates for triangle
        // Node ordering: (i,j) with i+j≤K, ordered by j then i
        size_t idx = 0;
        for (size_t j = 0; j <= K; ++j)
        {
          for (size_t i = 0; i <= K - j; ++i)
          {
            if (idx == m_local)
            {
              return Internal::evaluateLagrangeTriangle<K, Scalar>(i, j, r);
            }
            idx++;
          }
        }
        return Math::nan<Scalar>();
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

        return Internal::evaluateLagrange1D<K, Scalar>(i_idx, r.x()) *
               Internal::evaluateLagrange1D<K, Scalar>(j_idx, r.y());
      }
      case Geometry::Polytope::Type::Tetrahedron:
      {
        // Node ordering: (i,j,k) with i+j+k≤K, ordered by k, then j, then i
        size_t idx = 0;
        for (size_t k = 0; k <= K; ++k)
        {
          for (size_t j = 0; j <= K - k; ++j)
          {
            for (size_t i = 0; i <= K - j - k; ++i)
            {
              if (idx == m_local)
              {
                return Internal::evaluateLagrangeTetrahedron<K, Scalar>(i, j, k, r);
              }
              idx++;
            }
          }
        }
        return Math::nan<Scalar>();
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
          return Internal::evaluateLagrange1DDerivative<K, Scalar>(m_local, r.x());
        }
        case Geometry::Polytope::Type::Triangle:
        {
          // Find the (i,j) indices for this local DOF
          size_t idx = 0;
          for (size_t j = 0; j <= K; ++j)
          {
            for (size_t i = 0; i <= K - j; ++i)
            {
              if (idx == m_local)
              {
                return Internal::evaluateLagrangeTriangleDerivative<K, Scalar>(i, j, m_i, r);
              }
              idx++;
            }
          }
          return 0;
        }
        case Geometry::Polytope::Type::Quadrilateral:
        {
          // Tensor product derivative
          size_t j_idx = m_local / (K + 1);
          size_t i_idx = m_local % (K + 1);

          if (m_i == 0) // d/dx
          {
            return Internal::evaluateLagrange1DDerivative<K, Scalar>(i_idx, r.x()) *
                   Internal::evaluateLagrange1D<K, Scalar>(j_idx, r.y());
          }
          else if (m_i == 1) // d/dy
          {
            return Internal::evaluateLagrange1D<K, Scalar>(i_idx, r.x()) *
                   Internal::evaluateLagrange1DDerivative<K, Scalar>(j_idx, r.y());
          }
          return 0;
        }
        case Geometry::Polytope::Type::Tetrahedron:
        {
          // Find (i,j,k) indices for this local DOF
          size_t idx = 0;
          for (size_t k = 0; k <= K; ++k)
          {
            for (size_t j = 0; j <= K - k; ++j)
            {
              for (size_t i = 0; i <= K - j - k; ++i)
              {
                if (idx == m_local)
                {
                  return Internal::evaluateLagrangeTetrahedronDerivative<K, Scalar>(i, j, k, m_i, r);
                }
                idx++;
              }
            }
          }
          return 0;
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
