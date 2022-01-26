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
      /**
       * @brief Constructs an Advect2D object
       * @param[in, out] ls Function to advect
       * @param[in] disp Displacement velocity field
       */
      Advect2D(ScalarSolution2D& ls, VectorSolution2D& disp);

      /**
       * @brief Specifies whether to enable or disable truncation of time in the CFL condition.
       *
       * If true, avoids truncation of the time period for advection due to
       * the CFL condition. Otherwise, no action is taken.
       * By default it is set to false.
       */
      Advect2D& avoidTimeTruncation(bool cfl);

      /**
       * @brief Specifies whether to extrapolate characteristic lines or not.
       *
       * If true, characteristic lines are extrapolated outside the domain when
       * the input velocity field causes them to do so. Otherwise, no action is
       * taken.
       * By default it is set to true.
       */
      Advect2D& enableExtrapolation(bool ex = true);

      /**
       * @brief Advances the level set function by the time step `dt`.
       * @param[in] dt Advection timestep
       */
      void step(double dt);

    private:
      double m_t;
      bool m_cfl, m_ex;
      ScalarSolution2D& m_ls;
      VectorSolution2D& m_disp;
      ISCDProcess m_advect;
  };
}
#endif
