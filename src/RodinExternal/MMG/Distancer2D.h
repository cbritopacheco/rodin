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
  /**
   * @brief Class which performs the "distancing" of a domain.
   */
  class Distancer2D
  {
    public:
      /**
       * @brief Computes a signed distance function from a given bounding box
       * and contour.
       *
       * This function generates the signed distance function @f$ d_\Omega @f$
       * to a domain @f$ \Omega @f$ supplied by means of a mesh of its boundary
       * @f$ \partial \Omega @f$, at the vertices of a computational mesh of a
       * bounding box @f$ D @f$.
       *
       * @param[in] box Bounding box @f$ D @f$ containing the contour.
       * @param[in] contour Orientable contour @f$ \partial \Omega @f$ to distance.
       * @note The contour mesh is allowed to contain a volume part, in which
       * case only the edge (2D) or triangle (3D) information will be retained.
       */
      ScalarSolution2D<> distance(Mesh2D& box, Mesh2D& contour) const;

    private:
      MshdistProcess m_mshdist;
  };
}

#endif
