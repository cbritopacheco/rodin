/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_BARYCENTRICSUBDIVISION_H
#define RODIN_GEOMETRY_BARYCENTRICSUBDIVISION_H

#include "Mesh.h"

namespace Rodin::Geometry
{
  class BarycentricSubdivision
  {
    public:
      BarycentricSubdivision();

      void refine(Mesh<Context::Local>& mesh);

    private:
  };
}

#endif

