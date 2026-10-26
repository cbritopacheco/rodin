/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ATRIUM_SOLID_SOLID_H
#define RODIN_ATRIUM_SOLID_SOLID_H

#include <memory>
#include "ProblemData.h"

namespace Rodin::Examples::Heart::AtriumContraction
{
  /**
   * @brief Transient active atrial wall on a local wedge mesh.
   *
   * ## Architecture
   *
   * The component retains its reference mesh, prescribed fields, finite-element
   * spaces, active material, nonlinear problem, and solver between time steps.
   * Newmark updates advance displacement, velocity, and acceleration. Output
   * observes the retained fields. This component uses the local Eigen backend;
   * its mesh and interface are independent of the MPI fluid component.
   */
  class Solid
  {
    public:
      using Parameters = Options;

      /**
       * @brief Constructs the wall component and prescribed material inputs.
       * @param parameters Validated wall, material, wave, and time configuration.
       */
      explicit Solid(const Parameters& parameters);
      ~Solid();
      Solid(const Solid&) = delete;
      Solid& operator=(const Solid&) = delete;

      /**
       * @brief Computes the initial acceleration and writes the initial state.
       * @returns This component.
       */
      Solid& initialize();

      /**
       * @brief Advances one Newmark time step and writes the state.
       * @note A failed step cannot be retried on the same component.
       */
      void step();

      /// @brief Finalizes output; no further time steps may be taken.
      void finish();

      /// @brief Advances the configured number of intervals and finalizes output.
      void run();

      /// @brief Current physical time in seconds.
      Real getTime() const noexcept;

      /// @brief Configured time increment in seconds.
      Real getTimeStep() const noexcept;

      /// @brief Number of successfully completed time steps.
      size_t getStepCount() const noexcept;

      /// @brief Returns the wall configuration.
      const Parameters& getParameters() const noexcept;

      /// @brief Returns the reference wall mesh, with coordinates in mm.
      const Geometry::LocalMesh& getMesh() const noexcept;

      /// @brief Returns the current displacement in mm.
      const VectorField& getDisplacement() const noexcept;

      /// @brief Returns the current velocity in mm/s.
      const VectorField& getVelocity() const noexcept;

      /// @brief Returns the current acceleration in mm/s².
      const VectorField& getAcceleration() const noexcept;

    private:
      struct Impl;
      std::unique_ptr<Impl> m_impl;
  };
}

#endif
