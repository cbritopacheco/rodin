/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file P1Element.hpp
 * @brief Template implementation of P1 basis functions and derivatives.
 *
 * This file provides the explicit implementations of P1 basis functions
 * and their first derivatives for each supported polytope geometry.
 *
 * ## Basis Functions by Geometry
 *
 * **Point** (0D): @f$ \phi_0 = 1 @f$
 *
 * **Segment** (1D): Linear interpolation on @f$ [0,1] @f$
 * - @f$ \phi_0(x) = 1 - x @f$
 * - @f$ \phi_1(x) = x @f$
 *
 * **Triangle** (2D): Barycentric coordinates on reference triangle
 * - @f$ \phi_0(x,y) = 1 - x - y @f$
 * - @f$ \phi_1(x,y) = x @f$
 * - @f$ \phi_2(x,y) = y @f$
 *
 * **Quadrilateral** (2D): Bilinear functions on @f$ [0,1]^2 @f$
 * - @f$ \phi_0 = (1-x)(1-y) @f$
 * - @f$ \phi_1 = x(1-y) @f$
 * - @f$ \phi_2 = xy @f$
 * - @f$ \phi_3 = (1-x)y @f$
 *
 * **Tetrahedron** (3D): Barycentric coordinates
 * - @f$ \phi_0 = 1 - x - y - z @f$
 * - @f$ \phi_1 = x @f$
 * - @f$ \phi_2 = y @f$
 * - @f$ \phi_3 = z @f$
 *
 * **Wedge** (3D): Product of triangle and segment basis
 *
 * All derivatives are computed analytically as the partial derivatives
 * of these basis functions.
 */
#ifndef RODIN_VARIATIONAL_P1_P1ELEMENT_HPP
#define RODIN_VARIATIONAL_P1_P1ELEMENT_HPP

#include "Rodin/Math/Common.h"

#include "P1Element.h"

namespace Rodin::Variational
{
  template <class Scalar>
  constexpr
  Scalar P1Element<Scalar>::BasisFunction::operator()(const Math::SpatialPoint& r) const
  {
    switch (m_g)
    {
      case Geometry::Polytope::Type::Point:
      {
        return 1;
      }
      case Geometry::Polytope::Type::Segment:
      {
        switch (m_local)
        {
          case 0:
          {
            return 1 - r.x();
          }
          case 1:
          {
            return r.x();
          }
          default: [[unlikely]]
          {
            assert(false);
            return Math::nan<Scalar>();
          }
        }
      }
      case Geometry::Polytope::Type::Triangle:
      {
        switch (m_local)
        {
          case 0:
          {
            return -r.x() - r.y() + 1;
          }
          case 1:
          {
            return r.x();
          }
          case 2:
          {
            return r.y();
          }
          default: [[unlikely]]
          {
            assert(false);
            return Math::nan<Scalar>();
          }
        }
      }
      case Geometry::Polytope::Type::Quadrilateral:
      {
        switch (m_local)
        {
          case 0:
          {
            const auto& x = r.x();
            const auto& y = r.y();
            return x * y - x - y + 1;
          }
          case 1:
          {
            return r.x() * (1 - r.y());
          }
          case 2:
          {
            return r.x() * r.y();
          }
          case 3:
          {
            return r.y() * (1 - r.x());
          }
          default: [[unlikely]]
          {
            assert(false);
            return Math::nan<Scalar>();
          }
        }
      }
      case Geometry::Polytope::Type::Tetrahedron:
      {
        switch (m_local)
        {
          case 0:
          {
            return - r.x() - r.y() - r.z() + 1;
          }
          case 1:
          {
            return r.x();
          }
          case 2:
          {
            return r.y();
          }
          case 3:
          {
            return r.z();
          }
          default: [[unlikely]]
          {
            assert(false);
            return Math::nan<Scalar>();
          }
        }
      }
      case Geometry::Polytope::Type::Wedge:
      {
        switch (m_local)
        {
          case 0:
          {
            return r.x() * r.z() - r.x()  + r.y() * r.z() - r.y() - r.z() + 1;
          }
          case 1:
          {
            return r.x() * (1 - r.z());
          }
          case 2:
          {
            return r.y() * (1 - r.z());
          }
          case 3:
          {
            return r.z() * (1 - r.x() - r.y());
          }
          case 4:
          {
            return r.x() * r.z();
          }
          case 5:
          {
            return r.y() * r.z();
          }
          default: [[unlikely]]
          {
            assert(false);
            return Math::nan<Scalar>();
          }
        }
      }
    }
    assert(false);
    return Math::nan<Scalar>();
  }

  template <class Scalar>
  template <size_t Order>
  constexpr
  Scalar P1Element<Scalar>::BasisFunction::DerivativeFunction<Order>::operator()(
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
          assert(m_local == 0);
          return 0;
        }
        case Geometry::Polytope::Type::Segment:
        {
          switch (m_local)
          {
            case 0:
            {
              return -1;
            }
            case 1:
            {
              return 1;
            }
            default: [[unlikely]]
            {
              assert(false);
              return Math::nan<Scalar>();
            }
          }
        }
        case Geometry::Polytope::Type::Triangle:
        {
          switch (m_local)
          {
            case 0:
            {
              if (m_i == 0)
              {
                return -1;
              }
              else if (m_i == 1)
              {
                return -1;
              }
              else [[unlikely]]
              {
                assert(false);
                return Math::nan<Scalar>();
              }
            }
            case 1:
            {
              if (m_i == 0)
              {
                return 1;
              }
              else if (m_i == 1)
              {
                return 0;
              }
              else [[unlikely]]
              {
                assert(false);
                return Math::nan<Scalar>();
              }
            }
            case 2:
            {
              if (m_i == 0)
              {
                return 0;
              }
              else if (m_i == 1)
              {
                return 1;
              }
              else [[unlikely]]
              {
                assert(false);
                return Math::nan<Scalar>();
              }
            }
            default: [[unlikely]]
            {
              assert(false);
              return Math::nan<Scalar>();
            }
          }
        }
        case Geometry::Polytope::Type::Quadrilateral:
        {
          switch (m_local)
          {
            case 0:
            {
              if (m_i == 0)
              {
                return r.y() - 1;
              }
              else if (m_i == 1)
              {
                return r.x() - 1;
              }
              else [[unlikely]]
              {
                assert(false);
                return Math::nan<Scalar>();
              }
            }
            case 1:
            {
              if (m_i == 0)
              {
                return 1 - r.y();
              }
              else if (m_i == 1)
              {
                return -r.x();
              }
              else [[unlikely]]
              {
                assert(false);
                return Math::nan<Scalar>();
              }
            }
            case 2:
            {
              if (m_i == 0)
              {
                return r.y();
              }
              else if (m_i == 1)
              {
                return r.x();
              }
              else [[unlikely]]
              {
                assert(false);
                return Math::nan<Scalar>();
              }
            }
            case 3:
            {
              if (m_i == 0)
              {
                return -r.y();
              }
              else if (m_i == 1)
              {
                return 1 - r.x();
              }
              else [[unlikely]]
              {
                assert(false);
                return Math::nan<Scalar>();
              }
            }
            default:
            {
              assert(false);
              return Math::nan<Scalar>();
            }
          }
        }
        case Geometry::Polytope::Type::Tetrahedron:
        {
          switch (m_local)
          {
            case 0:
            {
              if (m_i == 0)
              {
                return -1;
              }
              else if (m_i == 1)
              {
                return -1;
              }
              else if (m_i == 2)
              {
                return -1;
              }
              else [[unlikely]]
              {
                assert(false);
                return Math::nan<Scalar>();
              }
            }
            case 1:
            {
              if (m_i == 0)
              {
                return 1;
              }
              else if (m_i == 1)
              {
                return 0;
              }
              else if (m_i == 2)
              {
                return 0;
              }
              else [[unlikely]]
              {
                assert(false);
                return Math::nan<Scalar>();
              }
            }
            case 2:
            {
              if (m_i == 0)
              {
                return 0;
              }
              else if (m_i == 1)
              {
                return 1;
              }
              else if (m_i == 2)
              {
                return 0;
              }
              else [[unlikely]]
              {
                assert(false);
                return Math::nan<Scalar>();
              }
            }
            case 3:
            {
              if (m_i == 0)
              {
                return 0;
              }
              else if (m_i == 1)
              {
                return 0;
              }
              else if (m_i == 2)
              {
                return 1;
              }
              else [[unlikely]]
              {
                assert(false);
                return Math::nan<Scalar>();
              }
            }
            default: [[unlikely]]
            {
              assert(false);
              return Math::nan<Scalar>();
            }
          }
        }
        case Geometry::Polytope::Type::Wedge:
        {
          switch (m_local)
          {
            case 0:
            {
              if (m_i == 0)
              {
                return r.z() - 1;
              }
              else if (m_i == 1)
              {
                return r.z() - 1;
              }
              else if (m_i == 2)
              {
                return r.x() + r.y() - 1;
              }
              else [[unlikely]]
              {
                assert(false);
                return Math::nan<Scalar>();
              }
            }
            case 1:
            {
              if (m_i == 0)
              {
                return 1 - r.z();
              }
              else if (m_i == 1)
              {
                return 0;
              }
              else if (m_i == 2)
              {
                return -r.x();
              }
              else [[unlikely]]
              {
                assert(false);
                return Math::nan<Scalar>();
              }
            }
            case 2:
            {
              if (m_i == 0)
              {
                return 0;
              }
              else if (m_i == 1)
              {
                return 1 - r.z();
              }
              else if (m_i == 2)
              {
                return -r.y();
              }
              else [[unlikely]]
              {
                assert(false);
                return Math::nan<Scalar>();
              }
            }
            case 3:
            {
              if (m_i == 0)
              {
                return -r.z();
              }
              else if (m_i == 1)
              {
                return -r.z();
              }
              else if (m_i == 2)
              {
                return 1 - r.x() - r.y();
              }
              else [[unlikely]]
              {
                assert(false);
                return Math::nan<Scalar>();
              }
            }
            case 4:
            {
              if (m_i == 0)
              {
                return r.z();
              }
              else if (m_i == 1)
              {
                return 0;
              }
              else if (m_i == 2)
              {
                return r.x();
              }
              else [[unlikely]]
              {
                assert(false);
                return Math::nan<Scalar>();
              }
            }
            case 5:
            {
              if (m_i == 0)
              {
                return 0;
              }
              else if (m_i == 1)
              {
                return r.z();
              }
              else if (m_i == 2)
              {
                return r.y();
              }
              else [[unlikely]]
              {
                assert(false);
                return Math::nan<Scalar>();
              }
            }
            default: [[unlikely]]
            {
              assert(false);
              return Math::nan<Scalar>();
            }
          }
        }
      }
      assert(false);
      return Math::nan<Scalar>();
    }
    else
    {
      return 0;
    }
  }
}

#endif

