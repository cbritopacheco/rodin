/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ATRIUM_PARAMETERS_H
#define RODIN_ATRIUM_PARAMETERS_H

#include "RunParameters.h"
#include "Solid/Solid.h"
#include "Fluid/Fluid.h"
#include "FSI/FSI.h"

namespace Rodin::Examples::Heart::Atrium
{
  /// @brief Components advanced by the atrium coordinator.
  enum class Mode
  {
    Solid,
    Fluid,
    FSI
  };

  /**
   * @brief Authoritative configuration of a complete atrium simulation.
   *
   * RunParameters supplies the common time and output configuration. The mesh
   * is loaded and partitioned once; physical and numerical parameters are grouped
   * by component. Copies passed to components have the same run configuration.
   */
  struct Parameters : RunParameters
  {
    Mode mode = Mode::Fluid;
    std::string meshPath =
      "resources/examples/Heart/Atrium/AtriumCoarse_FluidSolid_1mm.mesh";
    Geometry::Attribute solidMaterial = 4;
    Geometry::Attribute fluidMaterial = 1;
    Solid::Parameters solid;
    Fluid::Parameters fluid;
    FSI::Parameters fsi;

    /**
     * @brief Reads common and component options from the PETSc database collectively.
     * @returns Validated simulation parameters.
     */
    static Parameters fromOptions();

    /**
     * @brief Checks common parameters before any mesh is loaded.
     * @returns A validated copy of the configuration.
     */
    Parameters validated() const;
  };
}

#endif
