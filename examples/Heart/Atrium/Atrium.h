/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ATRIUM_ATRIUM_H
#define RODIN_ATRIUM_ATRIUM_H

#include <memory>
#include <fstream>
#include "CavityVolume.h"
#include "Parameters.h"
#include "Solid/Solid.h"
#include "Fluid/Fluid.h"
#include "FSI/FSI.h"

namespace Rodin::Examples::Heart::Atrium
{
  /**
   * @brief Advances atrial solid, fluid, or contraction-driven coupled dynamics.
   *
   * FSI mode exchanges motion and traction on the shared material interface.
   * Each component retains its mesh, problem, solver, and physical history.
   */
  class Atrium
  {
    public:
      using MeshType = Geometry::Mesh<Context::MPI>;
      using Solid = ::Rodin::Examples::Heart::Atrium::Solid;
      using Fluid = ::Rodin::Examples::Heart::Atrium::Fluid;
      using Mode = ::Rodin::Examples::Heart::Atrium::Mode;
      using Parameters = ::Rodin::Examples::Heart::Atrium::Parameters;

      /// @brief Constructs the selected model using the shared MPI context.
      Atrium(const Context::MPI& context, const Parameters& parameters);
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

      /// @brief Returns the combined distributed parent mesh in source length units.
      const MeshType& getMesh() const noexcept;

      /// @brief Returns the configured solid, or null when it is disabled.
      const Solid* getSolid() const noexcept;

      /// @brief Returns the configured fluid, or null when it is disabled.
      const Fluid* getFluid() const noexcept;

      /// @brief Returns the latest solid-enclosed cavity volume in mL, or empty when disabled.
      Optional<Real> getSolidCavityVolume() const noexcept;

      /// @brief Returns the latest fluid-domain volume in mL, or empty when disabled.
      Optional<Real> getFluidDomainVolume() const noexcept;

    private:
      /// @brief Collectively publishes initial and stepwise run diagnostics.
      void writeDiagnostics();

      /**
       * @brief Collectively loads and Scotch-partitions the combined source once.
       * @param context Communicator used by both component submeshes.
       * @param parameters Common source path and material labels.
       * @returns Reconciled parent mesh in source length units.
       */
      static MeshType makeMesh(const Context::MPI& context, const Parameters& parameters);

      // Source geometry uses mm; fluid component geometry uses m.
      static constexpr Real CubicMillimetresToMillilitres = 1e-3;
      static constexpr Real CubicMetresToMillilitres = 1e6;
      Parameters m_parameters;
      MeshType m_mesh;
      CavityVolume m_cavity;
      std::string m_csvPath;
      std::ofstream m_csvOutput;
      Optional<Real> m_solidCavityVolume;
      Optional<Real> m_fluidDomainVolume;
      Optional<Solid> m_solid;
      Optional<Fluid> m_fluid;
      Optional<FSI> m_fsi;
      size_t m_steps = 0;
      bool m_initialized = false;
      bool m_finished = false;
      bool m_failed = false;
  };
}

#endif
