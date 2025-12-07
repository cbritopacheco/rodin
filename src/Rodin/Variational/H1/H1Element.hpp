/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_H1ELEMENT_HPP
#define RODIN_VARIATIONAL_H1_H1ELEMENT_HPP

/**
 * @file
 * @brief Template implementation for H1Element class.
 *
 * This file contains the template member function implementations for
 * H1Element, including:
 * - Basis function evaluation using modal Dubiner basis
 * - First-order derivative computation via chain rule
 * - Vandermonde matrix inversion for nodal-to-modal conversion
 *
 * The implementation uses:
 * - **Segments**: Classical Lagrange interpolation with GLL nodes
 * - **Triangles**: Dubiner modal basis with Fekete nodes
 * - **Quadrilaterals**: Tensor product of 1D GLL Lagrange bases
 * - **Tetrahedra**: 3D Dubiner modal basis with Fekete nodes
 * - **Wedges**: Product of triangle Fekete basis and 1D GLL basis
 */

#include "Rodin/Math/Common.h"

#include "H1Element.h"

#include "LagrangeBasis.h"
#include "Dubiner.h"
#include "Fekete.h"
#include "WarpBlend.h"
#include "GLL.h"
#include "LegendrePolynomial.h"

namespace Rodin::Variational
{
  template <size_t K, class Scalar>
  constexpr
  const Math::SpatialPoint& H1Element<K, Scalar>::getNode(size_t i) const
  {
    return this->getNodes(this->getGeometry())[i];
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
      case Geometry::Polytope::Type::Hexahedron:
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
      case Geometry::Polytope::Type::Hexahedron:
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
  Scalar H1Element<K, Scalar>::BasisFunction::operator()(
      const Math::SpatialPoint& r) const
  {
    switch (m_g)
    {
      case Geometry::Polytope::Type::Point:
      {
        return Scalar(1);
      }
      case Geometry::Polytope::Type::Segment:
      {
        // Canonical [0,1] segment with GLL01 nodes
        return LagrangeBasisSegment<K>::getBasis(m_local, r.value());
      }
      case Geometry::Polytope::Type::Triangle:
      {
        // Modal Dubiner basis + Vandermonde on Fekete nodes
        const auto& inverse = VandermondeTriangle<K>::getInverse();

        Real rc, sc;
        DubinerTriangle<K>::getCollapsed(rc, sc, r.x(), r.y());

        Scalar result = Scalar(0);
        size_t mode_idx = 0;

        Rodin::Utility::ForIndex<K + 1>(
            [&](auto p_idx)
            {
              constexpr size_t P = p_idx.value;
              Rodin::Utility::ForIndex<K + 1 - P>(
                  [&](auto q_idx)
                  {
                    constexpr size_t Q = q_idx.value;
                    Real psi;
                    DubinerTriangle<K>::template getBasis<P, Q>(psi, rc, sc);
                    result += inverse(mode_idx, m_local) * psi;
                    ++mode_idx;
                  });
            });

        return result;
      }
      case Geometry::Polytope::Type::Quadrilateral:
      {
        // Tensor product GLL01 × GLL01 (Lagrange)
        const size_t j_idx = m_local / (K + 1);
        const size_t i_idx = m_local % (K + 1);
        return LagrangeBasisQuadrilateral<K>::getBasis(
            i_idx, j_idx, r.x(), r.y());
      }
      case Geometry::Polytope::Type::Tetrahedron:
      {
        // Modal Dubiner basis + Vandermonde on tetra Fekete nodes
        const auto& inverse = VandermondeTetrahedron<K>::getInverse();

        Real ac, bc, cc;
        DubinerTetrahedron<K>::getCollapsed(
            ac, bc, cc, r.x(), r.y(), r.z());

        Real result = 0;
        size_t mode_idx = 0;

        Rodin::Utility::ForIndex<K + 1>(
            [&](auto p_idx)
            {
              constexpr size_t P = p_idx.value;
              Rodin::Utility::ForIndex<K + 1 - P>(
                  [&](auto q_idx)
                  {
                    constexpr size_t Q = q_idx.value;
                    Rodin::Utility::ForIndex<K + 1 - P - Q>(
                        [&](auto r_idx)
                        {
                          constexpr size_t R = r_idx.value;
                          Real psi;
                          DubinerTetrahedron<K>::template getBasis<P, Q, R>(
                              psi, ac, bc, cc);
                          result += inverse(mode_idx, m_local) * psi;
                          ++mode_idx;
                        });
                  });
            });

        return result;
      }
      case Geometry::Polytope::Type::Wedge:
      {
        // Number of triangle DOFs (Fekete triangle)
        constexpr size_t ntri = FeketeTriangle<K>::Count; // same as FeketeTriangle<K>::Count

        const size_t k     = m_local / ntri; // z-index, 0..K
        const size_t alpha = m_local % ntri; // triangle node index

        // --- triangle factor: nodal basis on Fekete nodes ---
        const auto& Vinv = VandermondeTriangle<K>::getInverse();

        Real rc, sc;
        DubinerTriangle<K>::getCollapsed(rc, sc, r.x(), r.y());

        Scalar tri_val = Scalar(0);
        size_t mode_idx = 0;

        Rodin::Utility::ForIndex<K + 1>(
            [&](auto p_idx)
            {
              constexpr size_t P = p_idx.value;
              Rodin::Utility::ForIndex<K + 1 - P>(
                  [&](auto q_idx)
                  {
                    constexpr size_t Q = q_idx.value;
                    Real psi;
                    DubinerTriangle<K>::template getBasis<P, Q>(psi, rc, sc);
                    tri_val += Vinv(mode_idx, alpha) * psi;
                    ++mode_idx;
                  });
            });

        // --- segment factor in z (GLL01, Lagrange) ---
        const Real z = r(2);
        const Real seg_val = LagrangeBasisSegment<K>::getBasis(k, z);

        return tri_val * seg_val;
      }
      case Geometry::Polytope::Type::Hexahedron:
      {
        // Tensor product GLL01 × GLL01 × GLL01
        // Node ordering matches getNodes(): i + (K+1)*(j + (K+1)*k)
        const size_t n1 = K + 1;
        const size_t k  = m_local / (n1 * n1);
        const size_t r2 = m_local % (n1 * n1);
        const size_t j  = r2 / n1;
        const size_t i  = r2 % n1;

        const Real x = r.x();
        const Real y = r.y();
        const Real z = r.z();

        const Real lx = LagrangeBasisSegment<K>::getBasis(i, x);
        const Real ly = LagrangeBasisSegment<K>::getBasis(j, y);
        const Real lz = LagrangeBasisSegment<K>::getBasis(k, z);

        return static_cast<Scalar>(lx * ly * lz);
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
          return Scalar(0);
        }

        case Geometry::Polytope::Type::Segment:
        {
          // ∂φ_i / ∂x on [0,1] with GLL01 nodes
          assert(m_i == 0);
          return LagrangeBasisSegment<K>::getDerivative(
              m_local, r.value());
        }

        case Geometry::Polytope::Type::Triangle:
        {
          // Dubiner modal gradients + chain rule (r,s) → (x,y)
          const auto& Vinv = VandermondeTriangle<K>::getInverse();

          Scalar rc, sc;
          DubinerTriangle<K>::getCollapsed(rc, sc, r.x(), r.y());

          const Scalar x = r.x();
          const Scalar y = r.y();
          const Scalar eps = RODIN_VARIATIONAL_H1ELEMENT_TOLERANCE;

          Scalar result = Scalar(0);
          size_t mode_idx = 0;

          Rodin::Utility::ForIndex<K + 1>(
              [&](auto p_idx)
              {
                constexpr size_t P = p_idx.value;
                Rodin::Utility::ForIndex<K + 1 - P>(
                    [&](auto q_idx)
                    {
                      constexpr size_t Q = q_idx.value;

                      Scalar dpsi_dr = Scalar(0), dpsi_ds = Scalar(0);
                      DubinerTriangle<K>::template getGradient<P, Q>(
                          dpsi_dr, dpsi_ds, rc, sc);

                      Scalar dpsi_dx = Scalar(0), dpsi_dy = Scalar(0);

                      // r = 2x/(1-y) - 1, s = 2y - 1
                      if (Math::abs(Scalar(1) - y) > eps)
                      {
                        const Scalar denom = Scalar(1) - y;
                        const Scalar dr_dx = Scalar(2) / denom;
                        const Scalar dr_dy =
                            Scalar(2) * x / (denom * denom);
                        const Scalar ds_dx = Scalar(0);
                        const Scalar ds_dy = Scalar(2);

                        dpsi_dx = dpsi_dr * dr_dx + dpsi_ds * ds_dx;
                        dpsi_dy = dpsi_dr * dr_dy + dpsi_ds * ds_dy;
                      }

                      if (m_i == 0)      // ∂/∂x
                        result += Vinv(mode_idx, m_local) * dpsi_dx;
                      else if (m_i == 1) // ∂/∂y
                        result += Vinv(mode_idx, m_local) * dpsi_dy;

                      ++mode_idx;
                    });
              });

          return result;
        }

        case Geometry::Polytope::Type::Quadrilateral:
        {
          // Tensor product Lagrange on GLL01 × GLL01
          const size_t j_idx = m_local / (K + 1);
          const size_t i_idx = m_local % (K + 1);

          if (m_i == 0) // ∂/∂x
          {
            return LagrangeBasisQuadrilateral<K>::getDerivative(
                       i_idx, j_idx, 0, r.x(), r.y());
          }
          else if (m_i == 1) // ∂/∂y
          {
            return LagrangeBasisQuadrilateral<K>::getDerivative(
                       i_idx, j_idx, 1, r.x(), r.y());
          }
          return Scalar(0);
        }

        case Geometry::Polytope::Type::Tetrahedron:
        {
          // Dubiner modal gradients + chain rule (a,b,c) → (x,y,z)
          const auto& Vinv = VandermondeTetrahedron<K>::getInverse();

          Scalar ac, bc, cc;
          DubinerTetrahedron<K>::getCollapsed(
              ac, bc, cc, r.x(), r.y(), r.z());

          const Scalar x = r.x();
          const Scalar y = r.y();
          const Scalar z = r.z();
          const Scalar eps = RODIN_VARIATIONAL_H1ELEMENT_TOLERANCE;

          Scalar result = Scalar(0);
          size_t mode_idx = 0;

          Rodin::Utility::ForIndex<K + 1>(
              [&](auto p_idx)
              {
                constexpr size_t P = p_idx.value;
                Rodin::Utility::ForIndex<K + 1 - P>(
                    [&](auto q_idx)
                    {
                      constexpr size_t Q = q_idx.value;
                      Rodin::Utility::ForIndex<K + 1 - P - Q>(
                          [&](auto r_idx)
                          {
                            constexpr size_t R = r_idx.value;

                            Scalar dpsi_da = Scalar(0);
                            Scalar dpsi_db = Scalar(0);
                            Scalar dpsi_dc = Scalar(0);
                            DubinerTetrahedron<K>::template getGradient<P, Q, R>(
                                dpsi_da, dpsi_db, dpsi_dc, ac, bc, cc);

                            Scalar dpsi_dx = Scalar(0);
                            Scalar dpsi_dy = Scalar(0);
                            Scalar dpsi_dz = Scalar(0);

                            const Scalar denom2 = Scalar(1) - z;       // 1 - z
                            const Scalar denom3 = Scalar(1) - y - z;   // 1 - y - z

                            if (Math::abs(denom2) > eps &&
                                Math::abs(denom3) > eps)
                            {
                              // a = 2x / (1 - y - z) - 1
                              const Scalar da_dx = Scalar(2) / denom3;
                              const Scalar da_dy =
                                  Scalar(2) * x / (denom3 * denom3);
                              const Scalar da_dz = da_dy;

                              // b = 2y / (1 - z) - 1
                              const Scalar db_dx = Scalar(0);
                              const Scalar db_dy = Scalar(2) / denom2;
                              const Scalar db_dz =
                                  Scalar(2) * y / (denom2 * denom2);

                              // c = 2z - 1
                              const Scalar dc_dx = Scalar(0);
                              const Scalar dc_dy = Scalar(0);
                              const Scalar dc_dz = Scalar(2);

                              dpsi_dx = dpsi_da * da_dx
                                      + dpsi_db * db_dx
                                      + dpsi_dc * dc_dx;

                              dpsi_dy = dpsi_da * da_dy
                                      + dpsi_db * db_dy
                                      + dpsi_dc * dc_dy;

                              dpsi_dz = dpsi_da * da_dz
                                      + dpsi_db * db_dz
                                      + dpsi_dc * dc_dz;
                            }

                            if (m_i == 0)      // ∂/∂x
                              result += Vinv(mode_idx, m_local) * dpsi_dx;
                            else if (m_i == 1) // ∂/∂y
                              result += Vinv(mode_idx, m_local) * dpsi_dy;
                            else if (m_i == 2) // ∂/∂z
                              result += Vinv(mode_idx, m_local) * dpsi_dz;

                            ++mode_idx;
                          });
                    });
              });

          return result;
        }

        case Geometry::Polytope::Type::Wedge:
        {
          constexpr size_t ntri = FeketeTriangle<K>::Count;

          const size_t k = m_local / ntri;
          const size_t alpha = m_local % ntri;

          const Real z = r(2);

          if (m_i < 2) // ∂/∂x or ∂/∂y
          {
            // --- triangle gradient (same as Triangle case, but index = alpha) ---
            const auto& Vinv = VandermondeTriangle<K>::getInverse();

            Scalar rc, sc;
            DubinerTriangle<K>::getCollapsed(rc, sc, r.x(), r.y());

            const Scalar x   = r.x();
            const Scalar y   = r.y();
            const Scalar eps = RODIN_VARIATIONAL_H1ELEMENT_TOLERANCE;

            Scalar tri_deriv = Scalar(0);
            size_t mode_idx  = 0;

            Rodin::Utility::ForIndex<K + 1>(
                [&](auto p_idx)
                {
                  constexpr size_t P = p_idx.value;
                  Rodin::Utility::ForIndex<K + 1 - P>(
                      [&](auto q_idx)
                      {
                        constexpr size_t Q = q_idx.value;

                        Scalar dpsi_dr = Scalar(0), dpsi_ds = Scalar(0);
                        DubinerTriangle<K>::template getGradient<P, Q>(
                            dpsi_dr, dpsi_ds, rc, sc);

                        Scalar dpsi_dx = Scalar(0), dpsi_dy = Scalar(0);

                        // r = 2x/(1-y) - 1, s = 2y - 1
                        if (Math::abs(Scalar(1) - y) > eps)
                        {
                          const Scalar denom = Scalar(1) - y;
                          const Scalar dr_dx = Scalar(2) / denom;
                          const Scalar dr_dy =
                              Scalar(2) * x / (denom * denom);
                          const Scalar ds_dx = Scalar(0);
                          const Scalar ds_dy = Scalar(2);

                          dpsi_dx = dpsi_dr * dr_dx + dpsi_ds * ds_dx;
                          dpsi_dy = dpsi_dr * dr_dy + dpsi_ds * ds_dy;
                        }

                        if (m_i == 0)      // ∂/∂x
                          tri_deriv += Vinv(mode_idx, alpha) * dpsi_dx;
                        else if (m_i == 1) // ∂/∂y
                          tri_deriv += Vinv(mode_idx, alpha) * dpsi_dy;

                        ++mode_idx;
                      });
                });

            // --- segment value in z ---
            const Real seg_val = LagrangeBasisSegment<K>::getBasis(k, z);
            return tri_deriv * seg_val;
          }
          else // m_i == 2 → ∂/∂z
          {
            // --- triangle value (same as in BasisFunction wedge case) ---
            const auto& Vinv = VandermondeTriangle<K>::getInverse();

            Real rc, sc;
            DubinerTriangle<K>::getCollapsed(rc, sc, r.x(), r.y());

            Scalar tri_val = Scalar(0);
            size_t mode_idx = 0;

            Rodin::Utility::ForIndex<K + 1>(
                [&](auto p_idx)
                {
                  constexpr size_t P = p_idx.value;
                  Rodin::Utility::ForIndex<K + 1 - P>(
                      [&](auto q_idx)
                      {
                        constexpr size_t Q = q_idx.value;
                        Real psi;
                        DubinerTriangle<K>::template getBasis<P, Q>(psi, rc, sc);
                        tri_val += Vinv(mode_idx, alpha) * psi;
                        ++mode_idx;
                      });
                });

            // --- 1D derivative in z ---
            const Real dseg = LagrangeBasisSegment<K>::getDerivative(k, z);
            return tri_val * dseg;
          }
        }
        case Geometry::Polytope::Type::Hexahedron:
        {
          // Tensor product derivative on [0,1]^3 with GLL nodes.
          // Node ordering: i + (K+1)*(j + (K+1)*k)
          const size_t n1 = K + 1;
          const size_t k  = m_local / (n1 * n1);
          const size_t r2 = m_local % (n1 * n1);
          const size_t j  = r2 / n1;
          const size_t i  = r2 % n1;

          const Real x = r.x();
          const Real y = r.y();
          const Real z = r.z();

          const Real lx = LagrangeBasisSegment<K>::getBasis(i, x);
          const Real ly = LagrangeBasisSegment<K>::getBasis(j, y);
          const Real lz = LagrangeBasisSegment<K>::getBasis(k, z);

          const Real dlx = LagrangeBasisSegment<K>::getDerivative(i, x);
          const Real dly = LagrangeBasisSegment<K>::getDerivative(j, y);
          const Real dlz = LagrangeBasisSegment<K>::getDerivative(k, z);

          Real val = 0;
          if (m_i == 0)      // ∂/∂x
            val = dlx * ly * lz;
          else if (m_i == 1) // ∂/∂y
            val = lx * dly * lz;
          else if (m_i == 2) // ∂/∂z
            val = lx * ly * dlz;

          return static_cast<Scalar>(val);
        }
      }

      return Scalar(0);
    }
    else
    {
      // Higher-order derivatives not implemented
      return Scalar(0);
    }
  }
}

#endif
