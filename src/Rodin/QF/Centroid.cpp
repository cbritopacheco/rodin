/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Implementation of the Centroid quadrature formula.
 */

#include "Centroid.h"

namespace Rodin::QF
{
  const Math::SpatialVector<Real>& Centroid::getPoint(size_t i) const
  {
    assert(i == 0);
    switch (this->getGeometry())
    {
      case Geometry::Polytope::Type::Point:
      {
        static thread_local const Math::SpatialPoint s_point{};
        return s_point;
      }
      case Geometry::Polytope::Type::Segment:
      {
        static thread_local const Math::SpatialPoint s_point{ 0.5 };
        return s_point;
      }
      case Geometry::Polytope::Type::Triangle:
      {
        static thread_local const Math::SpatialPoint s_point{{ Real(1) / Real(3), Real(1) / Real(3) }};
        return s_point;
      }
      case Geometry::Polytope::Type::Quadrilateral:
      {
        static thread_local const Math::SpatialPoint s_point{{ 0.5, 0.5 }};
        return s_point;
      }
      case Geometry::Polytope::Type::Tetrahedron:
      {
        static thread_local const Math::SpatialPoint s_point{{ 0.25, 0.25, 0.25 }};
        return s_point;
      }
      case Geometry::Polytope::Type::Wedge:
      {
        static thread_local const Math::SpatialPoint s_point{{ Real(1) / Real(3), Real(1) / Real(3), 0.5 }};
        return s_point;
      }
      case Geometry::Polytope::Type::Hexahedron:
      {
        static thread_local const Math::SpatialPoint s_point{{ 0.5, 0.5, 0.5 }};
        return s_point;
      }
    }
    assert(false);
    static thread_local Math::SpatialPoint s_null;
    return s_null;
  }

  Real Centroid::getWeight(size_t i) const
  {
    assert(i == 0);
    switch (this->getGeometry())
    {
      case Geometry::Polytope::Type::Point:
        return 1;
      case Geometry::Polytope::Type::Segment:
        return 1;
      case Geometry::Polytope::Type::Triangle:
        return 0.5;
      case Geometry::Polytope::Type::Quadrilateral:
        return 1;
      case Geometry::Polytope::Type::Tetrahedron:
        return Real(1.0) / Real(6);
      case Geometry::Polytope::Type::Wedge:
        return 0.5;
      case Geometry::Polytope::Type::Hexahedron:
        return 1;
    }
    assert(false);
    return std::numeric_limits<Real>::quiet_NaN();
  }
}
