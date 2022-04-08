/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_ADVECTS_H
#define RODIN_EXTERNAL_MMG_ADVECTS_H

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
    *
    * @note Unlike Advect2D and Advect3D, time truncation is always avoided.
    */
  class AdvectS
  {
    public:
      /**
       * @brief Constructs an AdvectS object
       * @param[in, out] ls Function to advect
       * @param[in] disp Displacement velocity field
       */
      AdvectS(ScalarSolutionS& ls, VectorSolutionS& disp);

      /**
       * @brief Specifies whether to extrapolate characteristic lines or not.
       *
       * If true, characteristic lines are extrapolated outside the domain when
       * the input velocity field causes them to do so. Otherwise, no action is
       * taken.
       *
       * By default it is set to true.
       */
      AdvectS& enableExtrapolation(bool ex = true);

      /**
       * @brief Advances the level set function by the time step `dt`.
       * @param[in] dt Advection timestep
       */
      void step(double dt);

    private:
      double m_t;
      bool m_ex;
      ScalarSolutionS& m_ls;
      VectorSolutionS& m_disp;
      ISCDProcess m_advect;
  };
}
#endif

