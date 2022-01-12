/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_DISTANCER2D_H
#define RODIN_EXTERNAL_MMG_DISTANCER2D_H

#include "ForwardDecls.h"

#include "MshdistProcess.h"

namespace Rodin::External::MMG
{
  class Distancer2D
  {
    public:
      ScalarSolution2D<> distance(Mesh2D& box, Mesh2D& contour) const;

    private:
      MshdistProcess m_mshdist;
  };
}

#endif
