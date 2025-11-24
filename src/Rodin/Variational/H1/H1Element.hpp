/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_H1ELEMENT_HPP
#define RODIN_VARIATIONAL_H1_H1ELEMENT_HPP

#include "Rodin/Math/Common.h"

#include "H1Element.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace Rodin::Variational
{
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
            result += Vinv(mode_idx, m_local) * psi;
            mode_idx++;
          }
        }
        return result;
      }
      case Geometry::Polytope::Type::Quadrilateral:
      {
        // Tensor product of 1D Lagrange basis on GLL nodes in [0,1]
        size_t j_idx = m_local / (K + 1);
        size_t i_idx = m_local % (K + 1);

        // 1D GLL abscissas on [0,1]
        static thread_local std::vector<Math::SpatialPoint> s_nodes1D;
        if (s_nodes1D.empty())
        {
          auto xi = Internal::getGLLNodes01<K>();
          s_nodes1D.reserve(K + 1);
          for (size_t i = 0; i <= K; ++i)
            s_nodes1D.emplace_back(Math::SpatialPoint{{xi[i]}});
        }

        const auto& nodes1D = s_nodes1D;
        return Internal::evaluateLagrange1D<K, Scalar>(i_idx, r.x(), nodes1D) *
               Internal::evaluateLagrange1D<K, Scalar>(j_idx, r.y(), nodes1D);
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
              result += Vinv(mode_idx, m_local) * psi;
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
                result += Vinv(mode_idx, m_local) * dpsi_dx;
              else if (m_i == 1) // ∂/∂y
                result += Vinv(mode_idx, m_local) * dpsi_dy;

              mode_idx++;
            }
          }
          return result;
        }
        case Geometry::Polytope::Type::Quadrilateral:
        {
          // Tensor product derivative with clean 1D GLL nodes
          size_t j_idx = m_local / (K + 1);
          size_t i_idx = m_local % (K + 1);

          static thread_local std::vector<Math::SpatialPoint> s_nodes1D;
          if (s_nodes1D.empty())
          {
            auto xi = Internal::getGLLNodes01<K>();
            s_nodes1D.reserve(K + 1);
            for (size_t i = 0; i <= K; ++i)
              s_nodes1D.emplace_back(Math::SpatialPoint{{xi[i]}});
          }
          const auto& nodes1D = s_nodes1D;

          if (m_i == 0) // d/dx
          {
            return Internal::evaluateLagrange1DDerivative<K, Scalar>(i_idx, r.x(), nodes1D) *
                   Internal::evaluateLagrange1D<K, Scalar>(j_idx, r.y(), nodes1D);
          }
          else if (m_i == 1) // d/dy
          {
            return Internal::evaluateLagrange1D<K, Scalar>(i_idx, r.x(), nodes1D) *
                   Internal::evaluateLagrange1DDerivative<K, Scalar>(j_idx, r.y(), nodes1D);
          }
          return 0;
        }
        case Geometry::Polytope::Type::Tetrahedron:
        {
          const auto& Vinv = Internal::TetrahedronDubinerVandermonde<K>::getInverse();

          auto [ac, bc, cc] = Internal::tetrahedronToCollapsed(r.x(), r.y(), r.z());

          Scalar result = 0.0;
          size_t mode_idx = 0;
          for (size_t p = 0; p <= K; ++p)
          {
            for (size_t q = 0; q <= K - p; ++q)
            {
              for (size_t r_mode = 0; r_mode <= K - p - q; ++r_mode)
              {
                auto [dpsi_dx, dpsi_dy, dpsi_dz] =
                    Internal::dubinerTetrahedronGradient(p, q, r_mode, ac, bc, cc);

                if (m_i == 0)      // ∂/∂x
                  result += Vinv(mode_idx, m_local) * dpsi_dx;
                else if (m_i == 1) // ∂/∂y
                  result += Vinv(mode_idx, m_local) * dpsi_dy;
                else if (m_i == 2) // ∂/∂z
                  result += Vinv(mode_idx, m_local) * dpsi_dz;

                ++mode_idx;
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
