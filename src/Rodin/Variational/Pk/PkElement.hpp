/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_PK_PKELEMENT_HPP
#define RODIN_VARIATIONAL_PK_PKELEMENT_HPP

#include "Rodin/Math/Common.h"

#include "PkElement.h"

namespace Rodin::Variational
{
  namespace Internal
  {
    /**
     * @brief Computes the number of DOFs for a Pk element on a given geometry.
     */
    constexpr size_t getPkDofCount(size_t k, Geometry::Polytope::Type g)
    {
      switch (g)
      {
        case Geometry::Polytope::Type::Point:
          return 1;
        case Geometry::Polytope::Type::Segment:
          return k + 1;
        case Geometry::Polytope::Type::Triangle:
          return (k + 1) * (k + 2) / 2;
        case Geometry::Polytope::Type::Quadrilateral:
          return (k + 1) * (k + 1);
        case Geometry::Polytope::Type::Tetrahedron:
          return (k + 1) * (k + 2) * (k + 3) / 6;
        case Geometry::Polytope::Type::Wedge:
          return (k + 1) * (k + 1) * (k + 2) / 2;
      }
      return 0;
    }

    /**
     * @brief Generates Lagrange nodes for degree k on reference element.
     */
    template <size_t K>
    constexpr void generateLagrangeNodes(
        std::vector<Math::SpatialPoint>& nodes,
        Geometry::Polytope::Type g)
    {
      nodes.clear();
      
      switch (g)
      {
        case Geometry::Polytope::Type::Point:
        {
          nodes.push_back(Math::SpatialPoint{{0}});
          break;
        }
        case Geometry::Polytope::Type::Segment:
        {
          // Uniformly spaced nodes on [0, 1]
          for (size_t i = 0; i <= K; ++i)
          {
            Real t = static_cast<Real>(i) / static_cast<Real>(K);
            nodes.push_back(Math::SpatialPoint{{t}});
          }
          break;
        }
        case Geometry::Polytope::Type::Triangle:
        {
          // Nodes on reference triangle vertices (0,0), (1,0), (0,1)
          for (size_t j = 0; j <= K; ++j)
          {
            for (size_t i = 0; i <= K - j; ++i)
            {
              Real s = static_cast<Real>(i) / static_cast<Real>(K);
              Real t = static_cast<Real>(j) / static_cast<Real>(K);
              nodes.push_back(Math::SpatialPoint{{s, t}});
            }
          }
          break;
        }
        case Geometry::Polytope::Type::Quadrilateral:
        {
          // Tensor product nodes on [0,1]^2
          for (size_t j = 0; j <= K; ++j)
          {
            for (size_t i = 0; i <= K; ++i)
            {
              Real s = static_cast<Real>(i) / static_cast<Real>(K);
              Real t = static_cast<Real>(j) / static_cast<Real>(K);
              nodes.push_back(Math::SpatialPoint{{s, t}});
            }
          }
          break;
        }
        case Geometry::Polytope::Type::Tetrahedron:
        {
          // Nodes on reference tetrahedron vertices (0,0,0), (1,0,0), (0,1,0), (0,0,1)
          for (size_t k = 0; k <= K; ++k)
          {
            for (size_t j = 0; j <= K - k; ++j)
            {
              for (size_t i = 0; i <= K - j - k; ++i)
              {
                Real r = static_cast<Real>(i) / static_cast<Real>(K);
                Real s = static_cast<Real>(j) / static_cast<Real>(K);
                Real t = static_cast<Real>(k) / static_cast<Real>(K);
                nodes.push_back(Math::SpatialPoint{{r, s, t}});
              }
            }
          }
          break;
        }
        case Geometry::Polytope::Type::Wedge:
        {
          // Tensor product of triangle and segment
          for (size_t k = 0; k <= K; ++k)
          {
            for (size_t j = 0; j <= K; ++j)
            {
              for (size_t i = 0; i <= K - j; ++i)
              {
                Real r = static_cast<Real>(i) / static_cast<Real>(K);
                Real s = static_cast<Real>(j) / static_cast<Real>(K);
                Real t = static_cast<Real>(k) / static_cast<Real>(K);
                nodes.push_back(Math::SpatialPoint{{r, s, t}});
              }
            }
          }
          break;
        }
      }
    }

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
  }

  template <size_t K, class Scalar>
  constexpr
  size_t PkElement<K, Scalar>::getCount() const
  {
    return Internal::getPkDofCount(K, this->getGeometry());
  }

  template <size_t K, class Scalar>
  constexpr
  void PkElement<K, Scalar>::buildNodes()
  {
    Internal::generateLagrangeNodes<K>(m_nodes, this->getGeometry());
  }

  template <size_t K, class Scalar>
  const typename PkElement<K, Scalar>::LinearForm&
  PkElement<K, Scalar>::getLinearForm(size_t i) const
  {
    const Geometry::Polytope::Type g = this->getGeometry();
    static thread_local std::vector<LinearForm> s_lfs;
    
    // Thread-safe lazy initialization
    if (s_lfs.empty())
    {
      const size_t count = getCount();
      s_lfs.reserve(count);
      for (size_t j = 0; j < count; ++j)
        s_lfs.emplace_back(j, g);
    }
    
    return s_lfs[i];
  }

  template <size_t K, class Scalar>
  const typename PkElement<K, Scalar>::BasisFunction&
  PkElement<K, Scalar>::getBasis(size_t i) const
  {
    const Geometry::Polytope::Type g = this->getGeometry();
    static thread_local std::vector<BasisFunction> s_bs;
    
    // Thread-safe lazy initialization
    if (s_bs.empty())
    {
      const size_t count = getCount();
      s_bs.reserve(count);
      for (size_t j = 0; j < count; ++j)
        s_bs.emplace_back(j, g);
    }
    
    return s_bs[i];
  }

  template <size_t K, class Scalar>
  constexpr
  Scalar PkElement<K, Scalar>::BasisFunction::operator()(const Math::SpatialPoint& r) const
  {
    // Get nodes for this element
    PkElement<K, Scalar> elem(m_g);
    const auto& nodes = elem.m_nodes;
    
    if (m_local >= nodes.size())
      return Math::nan<Scalar>();
    
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
        // For simplicity, we use the tensor product approach
        // This is a simplified implementation
        size_t idx = 0;
        for (size_t j = 0; j <= K; ++j)
        {
          for (size_t i = 0; i <= K - j; ++i)
          {
            if (idx == m_local)
            {
              // Lagrange basis on triangle using area coordinates
              Real lambda0 = 1.0 - r.x() - r.y();
              Real lambda1 = r.x();
              Real lambda2 = r.y();
              
              // This is a placeholder - proper implementation would use
              // Silvester-Lagrange polynomials on triangles
              Scalar result = 1;
              const auto& node = nodes[m_local];
              Real node_lambda0 = 1.0 - node.x() - node.y();
              Real node_lambda1 = node.x();
              Real node_lambda2 = node.y();
              
              // Simplified Lagrange interpolation
              // In practice, this needs proper barycentric Lagrange polynomials
              Real tol = 1e-10;
              if (std::abs(lambda0 - node_lambda0) < tol &&
                  std::abs(lambda1 - node_lambda1) < tol &&
                  std::abs(lambda2 - node_lambda2) < tol)
                return 1;
              
              // For K=1 (P1), use standard linear basis
              if constexpr (K == 1)
              {
                if (i == 0 && j == 0)
                  return lambda0;
                else if (i == 1 && j == 0)
                  return lambda1;
                else if (i == 0 && j == 1)
                  return lambda2;
              }
              
              return 0; // Placeholder
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
        
        return Internal::evaluateLagrange1D<K, Scalar>(i_idx, r.x()) *
               Internal::evaluateLagrange1D<K, Scalar>(j_idx, r.y());
      }
      case Geometry::Polytope::Type::Tetrahedron:
      case Geometry::Polytope::Type::Wedge:
      {
        // Placeholder for 3D elements
        // Proper implementation would use 3D Lagrange polynomials
        return 0;
      }
    }
    
    return Math::nan<Scalar>();
  }

  template <size_t K, class Scalar>
  template <size_t Order>
  constexpr
  Scalar PkElement<K, Scalar>::BasisFunction::DerivativeFunction<Order>::operator()(
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
          // For K=1 (P1 element), use analytical derivatives
          if constexpr (K == 1)
          {
            size_t idx = 0;
            for (size_t j = 0; j <= K; ++j)
            {
              for (size_t i = 0; i <= K - j; ++i)
              {
                if (idx == m_local)
                {
                  if (i == 0 && j == 0) // lambda0 = 1 - x - y
                  {
                    return (m_i == 0) ? -1 : -1;
                  }
                  else if (i == 1 && j == 0) // lambda1 = x
                  {
                    return (m_i == 0) ? 1 : 0;
                  }
                  else if (i == 0 && j == 1) // lambda2 = y
                  {
                    return (m_i == 0) ? 0 : 1;
                  }
                }
                idx++;
              }
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
        case Geometry::Polytope::Type::Wedge:
        {
          return 0; // Placeholder
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
