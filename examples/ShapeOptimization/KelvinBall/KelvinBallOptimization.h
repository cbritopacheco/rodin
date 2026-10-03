/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef KELVIN_BALL_OPTIMIZATION_H
#define KELVIN_BALL_OPTIMIZATION_H

namespace KelvinBall
{
  /**
   * @brief Constrained shape optimization of the Kelvin ball.
   *
   * Each iteration follows the mathematical construction in order: chamber
   * discretization, fluid extraction, Stokes state solution, resistance
   * evaluation, shape differentiation, Hilbert identification, Feppon
   * null-space projection, sewn output, level-set advection, and MMG
   * reconstruction. The initial sphere and every accepted iterate are sewn
   * into a complete 24-copy design before the level set is advanced.
   *
   * Architecture:
   * - configuration fixes the chamber, discretization, and descent scales;
   * - state evaluation computes the translational and rotational families;
   * - descent construction identifies and projects the two shape gradients;
   * - output transports chamber fields to the complete reconstructed design;
   * - evolution advects the chamber level set and reconstructs its interface.
   */
  class KelvinBallOptimization
  {
    public:
      KelvinBallOptimization(int argc, char** argv);

      int run();

    private:
      class Implementation;

      int m_argc;
      char** m_argv;
  };
}

#endif
