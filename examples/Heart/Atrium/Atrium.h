/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ATRIUM_ATRIUM_H
#define RODIN_ATRIUM_ATRIUM_H

#include <memory>
#include "Solid/Solid.h"
#ifdef RODIN_ATRIUM_HAS_FLUID
#include "Fluid.h"
#endif

namespace Rodin::Examples::Heart
{
  /**
   * @brief Coordinates independent atrial solid and fluid time steps.
   *
   * Independent mode advances both components without exchanging interface
   * data. Each component owns its mesh, problem, solver, and output.
   */
  class Atrium
  {
    public:
      using Solid = AtriumContraction::Solid;
#ifdef RODIN_ATRIUM_HAS_FLUID
      using Fluid = AtriumFluid;
#endif
      enum class Mode
      {
        SolidOnly,
        FluidOnly,
        Independent
      };

      struct Parameters
      {
        Mode mode = Mode::SolidOnly;
        Solid::Parameters solid;
#ifdef RODIN_ATRIUM_HAS_FLUID
        Fluid::Parameters fluid;
#endif
      };

      /// @brief Constructs a solid-only model without an MPI context.
      explicit Atrium(const Parameters& parameters);
#ifdef RODIN_ATRIUM_HAS_FLUID
      /// @brief Constructs the selected model using an MPI context for the fluid.
      Atrium(const Context::MPI& context, const Parameters& parameters);
#endif
      ~Atrium();
      Atrium(const Atrium&) = delete;
      Atrium& operator=(const Atrium&) = delete;

      /// @brief Initializes the selected components and their outputs.
      Atrium& initialize();

      /// @brief Advances each selected component by one time step.
      void step();

      /// @brief Finalizes the selected outputs; no further steps are allowed.
      void finish();

      /// @brief Returns the number of completed coordinated steps.
      size_t getStepCount() const noexcept;

      /// @brief Returns the current physical time in seconds.
      Real getTime() const noexcept;

      /// @brief Returns the configured solid, or null when it is disabled.
      const Solid* getSolid() const noexcept;
#ifdef RODIN_ATRIUM_HAS_FLUID
      /// @brief Returns the configured fluid, or null when it is disabled.
      const Fluid* getFluid() const noexcept;
#endif

    private:
      struct Impl;
      std::unique_ptr<Impl> m_impl;
  };
}

#endif
