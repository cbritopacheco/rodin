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
     * @brief Evaluates Lagrange basis function for 1D.
     */
    template <size_t K, class Scalar>
    constexpr Scalar evaluateLagrange1D(size_t i, Real x)
    {
      Scalar result = 1;
      Real xi = static_cast<Real>(i) / static_cast<Real>(K);

      for (size_t j = 0; j <= K; ++j)
      {
        if (j != i)
        {
          Real xj = static_cast<Real>(j) / static_cast<Real>(K);
          result *= (x - xj) / (xi - xj);
        }
      }
      return result;
    }

    /**
     * @brief Evaluates derivative of Lagrange basis function for 1D.
     */
    template <size_t K, class Scalar>
    constexpr Scalar evaluateLagrange1DDerivative(size_t i, Real x)
    {
      Real xi = static_cast<Real>(i) / static_cast<Real>(K);
      Scalar result = 0;

      // Derivative using product rule
      for (size_t m = 0; m <= K; ++m)
      {
        if (m != i)
        {
          Scalar term = 1;
          Real xm = static_cast<Real>(m) / static_cast<Real>(K);

          for (size_t j = 0; j <= K; ++j)
          {
            if (j != i && j != m)
            {
              Real xj = static_cast<Real>(j) / static_cast<Real>(K);
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
  constexpr
  const Math::SpatialPoint& H1Element<K, Scalar>::getNode(size_t i) const
  {
    return getLagrangeNodes(this->getGeometry())[i];
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
