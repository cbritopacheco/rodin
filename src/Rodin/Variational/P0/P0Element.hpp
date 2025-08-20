/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P0_P0ELEMENT_HPP
#define RODIN_VARIATIONAL_P0_P0ELEMENT_HPP

#include "P0Element.h"

namespace Rodin::Variational
{
  template <class Scalar>
  const Geometry::GeometryIndexed<Math::SpatialVector<Real>>
  P0Element<Scalar>::s_nodes =
  {
    { Geometry::Polytope::Type::Point,
      Math::SpatialVector<Real>{{ 0 }}
    },
    { Geometry::Polytope::Type::Segment,
      Math::SpatialVector<Real>{{ 0.5 }}
    },
    { Geometry::Polytope::Type::Triangle,
      Math::SpatialVector<Real>{{ Real(1) / Real(3), Real(1) / Real(3) }}
    },
    { Geometry::Polytope::Type::Quadrilateral,
      Math::SpatialVector<Real>{{ 0.5, 0.5 }}
    },
    { Geometry::Polytope::Type::Tetrahedron,
      Math::SpatialVector<Real>{{ 0.25, 0.25, 0.25 }}
    },
    { Geometry::Polytope::Type::Wedge,
      Math::SpatialVector<Real>{{ Real(1) / Real(3), Real(1) / Real(3), 0.5 }}
    }
  };
}

#endif

