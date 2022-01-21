/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_REDISTANCER2D_H
#define RODIN_EXTERNAL_MMG_REDISTANCER2D_H

#include "Configure.h"
#include "ForwardDecls.h"
#include "ISCDProcess.h"

namespace Rodin::External::MMG
{
  /**
   * @brief Class which performs the "redistancing" of a level set function.
   */
  class Redistancer2D
  {
    public:
      Redistancer2D()
        : m_mshdist(MSHDIST_EXECUTABLE)
      {}

      /**
       * @brief Redistances a level set function.
       *
       * Given a level set function defined on the vertices of a bounding box
       * @f$ D @f$, this method will regenerate the signed distance function
       * associated to the the subdomain @f$ \Omega \subset D @f$ which the
       * level set function represents.
       *
       * @param[in,out] sol Level set function
       */
      void redistance(ScalarSolution2D& sol) const;

    private:
      ISCDProcess m_mshdist;
  };
}

#endif

