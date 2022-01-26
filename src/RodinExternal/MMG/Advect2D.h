/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_ADVECT2D_H
#define RODIN_EXTERNAL_MMG_ADVECT2D_H

#include <vector>
#include <optional>

#include "Common.h"
#include "Configure.h"
#include "ForwardDecls.h"
#include "ISCDProcess.h"

namespace Rodin::External::MMG
{
   /**
    * @brief Advection of a level set function by a velocity field.
    *
    * Solves the advection equation
    * @f[
    * \left\{
    *    \begin{aligned}
    *       \dfrac{\partial u}{\partial t} + v(x) \cdot \nabla u (t, x) &= 0
    *          && \text{ in } \Omega \times (0, + \infty) \\
    *       u(x, 0) &= u_0(x) && \text{ on } \Omega \times \{ t = 0 \}
    *    \end{aligned}
    * \right.
    * @f]
    * where @f$ u_0 : \mathbb{R}^d \rightarrow \mathbb{R} @f$ is known.
    */
  class Advect2D
  {
    public:
      Advect2D(ScalarSolution2D& ls, VectorSolution2D& disp)
        : m_ls(ls),
          m_disp(disp),
          m_advect(ADVECTION_EXECUTABLE)
      {}

      void step(double dt)
      {}

    private:
      ScalarSolution2D& m_ls;
      VectorSolution2D& m_disp;
      ISCDProcess m_advect;
  };
}
#endif
