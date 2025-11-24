/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_LAGRANGEBASIS_H
#define RODIN_VARIATIONAL_H1_LAGRANGEBASIS_H

#include <array>
#include <cstddef>
#include <functional>

#include "Rodin/Math/Vector.h"
#include "GLL.h" // GLL<K>, GLL01<K>

namespace Rodin::Variational
{
  //==========================================================================
  // Generic 1D Lagrange on arbitrary nodes (keeps nodes)
  //==========================================================================

  template <size_t K>
  class LagrangeBasis1D
  {
    public:
      static constexpr
      Real getBasis(size_t i, Real x, const std::array<Real, K + 1>& nodes)
      {
        Real result = 1;
        const Real xi = nodes[i];

        for (size_t j = 0; j <= K; ++j)
        {
          if (j != i)
          {
            const Real xj = nodes[j];
            result *= (x - xj) / (xi - xj);
          }
        }
        return result;
      }

      static constexpr
      Real getDerivative(size_t i, Real x, const std::array<Real, K + 1>& nodes)
      {
        const Real xi = nodes[i];
        Real result = 0;

        // Derivative using product rule
        for (size_t m = 0; m <= K; ++m)
        {
          if (m != i)
          {
            Real term = 1;
            const Real xm = nodes[m];

            for (size_t j = 0; j <= K; ++j)
            {
              if (j != i && j != m)
              {
                const Real xj = nodes[j];
                term *= (x - xj) / (xi - xj);
              }
            }
            term /= (xi - xm);
            result += term;
          }
        }
        return result;
      }
  };

  //==========================================================================
  // Point basis (reference 0D element)
  //==========================================================================

  template <size_t K>
  class LagrangeBasisPoint
  {
    public:
      static constexpr Real getBasis()
      {
        return 1.0;
      }

      static constexpr Real getDerivative()
      {
        return 0.0;
      }
  };

  //==========================================================================
  // Segment basis on canonical [0,1] with GLL01<K> nodes
  //==========================================================================

  template <size_t K>
  class LagrangeBasisSegment
  {
    public:
      // Node index i, 0 <= i <= K, x in [0,1]
      static constexpr Real getBasis(size_t i, Real x)
      {
        const auto& nodes = GLL01<K>::getNodes();

        Real result = 1;
        const Real xi = nodes[i];

        for (size_t j = 0; j <= K; ++j)
        {
          if (j != i)
          {
            const Real xj = nodes[j];
            result *= (x - xj) / (xi - xj);
          }
        }
        return result;
      }

      static constexpr Real getDerivative(size_t i, Real x)
      {
        const auto& nodes = GLL01<K>::getNodes();

        const Real xi = nodes[i];
        Real result = 0;

        for (size_t m = 0; m <= K; ++m)
        {
          if (m != i)
          {
            Real term = 1;
            const Real xm = nodes[m];

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
  };

  //==========================================================================
  // Triangle basis (Pk on reference triangle (0,0),(1,0),(0,1))
  //==========================================================================

  template <size_t K>
  class LagrangeBasisTriangle
  {
    public:
      // Node identified by (i,j) with i+j <= K.
      static constexpr Real getBasis(size_t i, size_t j, Real x, Real y)
      {
        // Barycentric coordinates: λ0 = 1 - x - y, λ1 = x, λ2 = y
        const Real lambda[3] = { 1.0 - x - y, x, y };

        // Node (i,j) ↔ barycentric ( (K-i-j)/K, i/K, j/K )
        const size_t indices[3] = { K - i - j, i, j };

        Real result = 1;
        for (size_t dim = 0; dim < 3; ++dim)
        {
          const size_t n = indices[dim];
          if (n == 0)
            continue;

          Real L_n = 1;
          for (size_t m = 0; m < n; ++m)
          {
            L_n *= (static_cast<Real>(K) * lambda[dim] - static_cast<Real>(m));
            L_n /= static_cast<Real>(m + 1);
          }
          result *= L_n;
        }

        return result;
      }

      // deriv_dim = 0 -> ∂/∂x, deriv_dim = 1 -> ∂/∂y
      static constexpr
      Real getDerivative(size_t i, size_t j, size_t deriv_dim, Real x, Real y)
      {
        const Real lambda[3] = { 1.0 - x - y, x, y };

        // dλ0/dx = -1, dλ0/dy = -1
        // dλ1/dx =  1, dλ1/dy =  0
        // dλ2/dx =  0, dλ2/dy =  1
        const Real dlambda[3][2] = { {-1, -1}, { 1, 0 }, { 0, 1 } };

        const size_t indices[3] = { K - i - j, i, j };

        Real result = 0;

        // Chain rule over λ_d
        for (size_t d = 0; d < 3; ++d)
        {
          Real term = dlambda[d][deriv_dim];

          for (size_t dim = 0; dim < 3; ++dim)
          {
            const size_t n = indices[dim];

            if (dim == d)
            {
              if (n == 0)
              {
                term = Real(0);
                break;
              }

              Real dL_n = 0;
              for (size_t p = 0; p < n; ++p)
              {
                Real prod = static_cast<Real>(K);
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
              if (n == 0)
                continue;

              Real L_n = 1;
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
  };

  //==========================================================================
  // Tetrahedron basis (Pk on reference tetra (0,0,0),(1,0,0),(0,1,0),(0,0,1))
  //==========================================================================

  template <size_t K>
  class LagrangeBasisTetrahedron
  {
    public:
      // Node (i,j,k) with i+j+k <= K.
      static constexpr Real getBasis(
        size_t i, size_t j, size_t k, Real x, Real y, Real z)
      {
        // Barycentric: λ0 = 1-x-y-z, λ1 = x, λ2 = y, λ3 = z
        const Real lambda[4] = { 1.0 - x - y - z, x, y, z };

        const size_t indices[4] = { K - i - j - k, i, j, k };

        Real result = 1;
        for (size_t dim = 0; dim < 4; ++dim)
        {
          const size_t n = indices[dim];
          if (n == 0)
            continue;

          Real L_n = 1;
          for (size_t m = 0; m < n; ++m)
          {
            L_n *= (static_cast<Real>(K) * lambda[dim] - static_cast<Real>(m));
            L_n /= static_cast<Real>(m + 1);
          }
          result *= L_n;
        }

        return result;
      }

      // deriv_dim = 0 -> ∂/∂x, 1 -> ∂/∂y, 2 -> ∂/∂z
      static constexpr Real getDerivative(
        size_t i, size_t j, size_t k, size_t deriv_dim,
        Real x, Real y, Real z)
      {
        const Real lambda[4] = { 1.0 - x - y - z, x, y, z };

        // dλ0 = (-1,-1,-1)
        // dλ1 = ( 1, 0, 0)
        // dλ2 = ( 0, 1, 0)
        // dλ3 = ( 0, 0, 1)
        const Real dlambda[4][3] = {
          { -1, -1, -1 },
          {  1,  0,  0 },
          {  0,  1,  0 },
          {  0,  0,  1 }
        };

        size_t indices[4] = { K - i - j - k, i, j, k };

        Real result = 0;

        for (size_t d = 0; d < 4; ++d)
        {
          Real term = dlambda[d][deriv_dim];

          for (size_t dim = 0; dim < 4; ++dim)
          {
            const size_t n = indices[dim];

            if (dim == d)
            {
              if (n == 0)
              {
                term = 0;
                break;
              }

              Real dL_n = 0;
              for (size_t p = 0; p < n; ++p)
              {
                Real prod = static_cast<Real>(K);
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
              if (n == 0)
                continue;

              Real L_n = 1;
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
  };

  //==========================================================================
  // Quadrilateral basis (tensor product on [0,1]×[0,1], GLL01<K> in each dir)
  //==========================================================================

  template <size_t K>
  class LagrangeBasisQuadrilateral
  {
    public:
      // Node (i,j), 0 ≤ i,j ≤ K, φ_{i,j}(x,y) = L_i^K(x) L_j^K(y)
      static constexpr Real getBasis(size_t i, size_t j, Real x, Real y)
      {
        const auto& nodes = GLL01<K>::getNodes();

        // L_i(x)
        Real Lix = 1;
        const Real xi = nodes[i];
        for (size_t m = 0; m <= K; ++m)
        {
          if (m != i)
          {
            const Real xm = nodes[m];
            Lix *= (x - xm) / (xi - xm);
          }
        }

        // L_j(y)
        Real Ljy = 1;
        const Real yj = nodes[j];
        for (size_t n = 0; n <= K; ++n)
        {
          if (n != j)
          {
            const Real yn = nodes[n];
            Ljy *= (y - yn) / (yj - yn);
          }
        }

        return Lix * Ljy;
      }

      // deriv_dim = 0 -> ∂/∂x, deriv_dim = 1 -> ∂/∂y
      static constexpr Real getDerivative(
          size_t i, size_t j, size_t deriv_dim, Real x, Real y)
      {
        const auto& nodes = GLL01<K>::getNodes();

        if (deriv_dim == 0)
        {
          // ∂/∂x L_i(x) * L_j(y)
          // dL_i/dx
          const Real xi = nodes[i];
          Real dLix = 0;

          for (size_t m = 0; m <= K; ++m)
          {
            if (m != i)
            {
              const Real xm = nodes[m];
              Real term = 1;

              for (size_t q = 0; q <= K; ++q)
              {
                if (q != i && q != m)
                {
                  const Real xq = nodes[q];
                  term *= (x - xq) / (xi - xq);
                }
              }
              term /= (xi - xm);
              dLix += term;
            }
          }

          // L_j(y)
          Real Ljy = 1;
          const Real yj = nodes[j];
          for (size_t n = 0; n <= K; ++n)
          {
            if (n != j)
            {
              const Real yn = nodes[n];
              Ljy *= (y - yn) / (yj - yn);
            }
          }

          return dLix * Ljy;
        }
        else
        {
          // ∂/∂y L_i(x) * dL_j/dy
          // L_i(x)
          Real Lix = 1;
          const Real xi = nodes[i];
          for (size_t m = 0; m <= K; ++m)
          {
            if (m != i)
            {
              const Real xm = nodes[m];
              Lix *= (x - xm) / (xi - xm);
            }
          }

          // dL_j/dy
          const Real yj = nodes[j];
          Real dLjy = 0;

          for (size_t n = 0; n <= K; ++n)
          {
            if (n != j)
            {
              const Real yn = nodes[n];
              Real term = 1;

              for (size_t q = 0; q <= K; ++q)
              {
                if (q != j && q != n)
                {
                  const Real yq = nodes[q];
                  term *= (y - yq) / (yj - yq);
                }
              }
              term /= (yj - yn);
              dLjy += term;
            }
          }

          return Lix * dLjy;
        }
      }
  };

  //==========================================================================
  // Wedge basis (triangle × segment, [0,1] in z, triangle as above)
  //==========================================================================

  template <size_t K>
  class LagrangeBasisWedge
  {
    public:
      // Node (i,j,k): triangle indices (i,j) with i+j ≤ K, segment index k.
      static constexpr Real getBasis(
        size_t i, size_t j, size_t k,
        Real x, Real y, Real z)
      {
        // Triangle part
        const Real lambda[3] = { 1.0 - x - y, x, y };

        const size_t tri_indices[3] = { K - i - j, i, j };

        Real tri_val = 1;
        for (size_t dim = 0; dim < 3; ++dim)
        {
          const size_t n = tri_indices[dim];
          if (n == 0)
            continue;

          Real L_n = 1;
          for (size_t m = 0; m < n; ++m)
          {
            L_n *= (static_cast<Real>(K) * lambda[dim] - static_cast<Real>(m));
            L_n /= static_cast<Real>(m + 1);
          }
          tri_val *= L_n;
        }

        // Segment part in z on GLL01<K>
        const auto& nodes = GLL01<K>::getNodes();
        Real seg_val = 1;
        const Real zk = nodes[k];
        for (size_t s = 0; s <= K; ++s)
        {
          if (s != k)
          {
            const Real zs = nodes[s];
            seg_val *= (z - zs) / (zk - zs);
          }
        }

        return tri_val * seg_val;
      }

      // deriv_dim = 0 -> ∂/∂x, 1 -> ∂/∂y, 2 -> ∂/∂z
      static constexpr Real getDerivative(
        size_t i, size_t j, size_t k, size_t deriv_dim,
        Real x, Real y, Real z)
      {
        const auto& seg_nodes = GLL01<K>::getNodes();

        if (deriv_dim < 2)
        {
          // d/dx or d/dy on triangle part, segment untouched
          // Triangle derivative
          const Real lambda[3] = { 1.0 - x - y, x, y };
          const Real dlambda[3][2] = { {-1, -1}, { 1, 0 }, { 0, 1 } };
          const size_t tri_indices[3] = { K - i - j, i, j };

          Real tri_deriv = 0;

          for (size_t d = 0; d < 3; ++d)
          {
            Real term = dlambda[d][deriv_dim];

            for (size_t dim = 0; dim < 3; ++dim)
            {
              const size_t n = tri_indices[dim];

              if (dim == d)
              {
                if (n == 0)
                {
                  term = 0;
                  break;
                }

                Real dL_n = 0;
                for (size_t p = 0; p < n; ++p)
                {
                  Real prod = static_cast<Real>(K);
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
                if (n == 0)
                  continue;

                Real L_n = 1;
                for (size_t m = 0; m < n; ++m)
                {
                  L_n *= (static_cast<Real>(K) * lambda[dim] - static_cast<Real>(m));
                  L_n /= static_cast<Real>(m + 1);
                }
                term *= L_n;
              }
            }

            tri_deriv += term;
          }

          // Segment value
          Real seg_val = 1;
          const Real zk = seg_nodes[k];
          for (size_t s = 0; s <= K; ++s)
          {
            if (s != k)
            {
              const Real zs = seg_nodes[s];
              seg_val *= (z - zs) / (zk - zs);
            }
          }

          return tri_deriv * seg_val;
        }
        else
        {
          // d/dz on segment part, triangle unchanged
          // Triangle value
          const Real lambda[3] = { 1.0 - x - y, x, y };

          size_t tri_indices[3] = { K - i - j, i, j };

          Real tri_val = 1;
          for (size_t dim = 0; dim < 3; ++dim)
          {
            size_t n = tri_indices[dim];
            if (n == 0)
              continue;

            Real L_n = 1;
            for (size_t m = 0; m < n; ++m)
            {
              L_n *= (static_cast<Real>(K) * lambda[dim] - static_cast<Real>(m));
              L_n /= static_cast<Real>(m + 1);
            }
            tri_val *= L_n;
          }

          // Segment derivative
          const Real zk = seg_nodes[k];
          Real dseg = 0;

          for (size_t s = 0; s <= K; ++s)
          {
            if (s != k)
            {
              const Real zs = seg_nodes[s];
              Real term = 1;

              for (size_t q = 0; q <= K; ++q)
              {
                if (q != k && q != s)
                {
                  const Real zq = seg_nodes[q];
                  term *= (z - zq) / (zk - zq);
                }
              }
              term /= (zk - zs);
              dseg += term;
            }
          }

          return tri_val * dseg;
        }
      }
  };
}

#endif
