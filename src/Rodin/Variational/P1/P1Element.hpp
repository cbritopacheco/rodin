/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P1_P1ELEMENT_HPP
#define RODIN_VARIATIONAL_P1_P1ELEMENT_HPP

#include "P1Element.h"

namespace Rodin::Variational
{
  template <size_t Order>
  constexpr
  Real P1Element<Real>::BasisFunction::DerivativeFunction<Order>::operator()(const Math::SpatialVector<Real>& r) const
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
            default:
            {
              assert(false);
              return Math::nan<Real>();
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
              else
              {
                assert(false);
                return Math::nan<Real>();
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
              else
              {
                assert(false);
                return Math::nan<Real>();
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
              else
              {
                assert(false);
                return Math::nan<Real>();
              }
            }
            default:
            {
              assert(false);
              return Math::nan<Real>();
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
              else
              {
                assert(false);
                return Math::nan<Real>();
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
              else
              {
                assert(false);
                return Math::nan<Real>();
              }
            }
            case 2:
            {
              if (m_i == 0)
              {
                return -r.y();
              }
              else if (m_i == 1)
              {
                return 1 - r.x();
              }
              else
              {
                assert(false);
                return Math::nan<Real>();
              }
            }
            case 3:
            {
              if (m_i == 0)
              {
                return r.y();
              }
              else if (m_i == 1)
              {
                return r.x();
              }
              else
              {
                assert(false);
                return Math::nan<Real>();
              }
            }
            default:
            {
              assert(false);
              return Math::nan<Real>();
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
              else
              {
                assert(false);
                return Math::nan<Real>();
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
              else
              {
                assert(false);
                return Math::nan<Real>();
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
              else
              {
                assert(false);
                return Math::nan<Real>();
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
              else
              {
                assert(false);
                return Math::nan<Real>();
              }
            }
            default:
            {
              assert(false);
              return Math::nan<Real>();
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
              else
              {
                assert(false);
                return Math::nan<Real>();
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
              else
              {
                assert(false);
                return Math::nan<Real>();
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
              else
              {
                assert(false);
                return Math::nan<Real>();
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
              else
              {
                assert(false);
                return Math::nan<Real>();
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
              else
              {
                assert(false);
                return Math::nan<Real>();
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
              else
              {
                assert(false);
                return Math::nan<Real>();
              }
            }
            default:
            {
              assert(false);
              return Math::nan<Real>();
            }
          }
        }
      }
      assert(false);
      return Math::nan<Real>();
    }
    else
    {
      return 0;
    }
  }

  constexpr
  Real RealP1Element::BasisFunction::operator()(const Math::SpatialVector<Real>& r) const
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
          default:
          {
            assert(false);
            return Math::nan<Real>();
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
          default:
          {
            assert(false);
            return Math::nan<Real>();
          }
        }
      }
      case Geometry::Polytope::Type::Quadrilateral:
      {
        switch (m_local)
        {
          case 0:
          {
            const Real x = r.x();
            const Real y = r.y();
            return x * y - x - y + 1;
          }
          case 1:
          {
            return r.x() * (1 - r.y());
          }
          case 2:
          {
            return r.y() * (1 - r.x());
          }
          case 3:
          {
            return r.x() * r.y();
          }
          default:
          {
            assert(false);
            return Math::nan<Real>();
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
          default:
          {
            assert(false);
            return Math::nan<Real>();
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
          default:
          {
            assert(false);
            return Math::nan<Real>();
          }
        }
      }
    }
    assert(false);
    return Math::nan<Real>();
  }
}

#endif

