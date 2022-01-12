/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_REDISTANCER2D_H
#define RODIN_EXTERNAL_MMG_REDISTANCER2D_H

#include "ForwardDecls.h"

#include "MshdistProcess.h"

namespace Rodin::External::MMG
{
  class Redistancer2D
  {
    public:
      void redistance(ScalarSolution2D<>& sol) const;

    private:
      MshdistProcess m_mshdist;
  };
}

#endif

