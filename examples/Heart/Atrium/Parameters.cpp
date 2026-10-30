/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Parameters.h"

#include <cmath>
#include <limits>
#include <type_traits>
#include <Rodin/Alert.h>

namespace Rodin::Examples::Heart::Atrium
{
  Parameters Parameters::fromOptions()
  {
    Parameters out;
    PetscErrorCode ierr = PETSC_SUCCESS;
    const PetscErrorCode optionError = [&]() -> PetscErrorCode {
      PetscFunctionBeginUser;
      PetscOptionsBegin(PETSC_COMM_WORLD, nullptr,
        "Atrium: common run and component parameters", nullptr);
      const auto real = [&](const char* name, const char* description, Real& value) {
        ierr = PetscOptionsReal(name, description, nullptr, value, &value, nullptr);
        assert(ierr == PETSC_SUCCESS);
      };
      const auto integer = [&](const char* name, const char* description, auto& value) {
        PetscInt parsed = static_cast<PetscInt>(value);
        ierr = PetscOptionsInt(name, description, nullptr, parsed, &parsed, nullptr);
        assert(ierr == PETSC_SUCCESS);
        if (parsed < 0 || parsed > std::numeric_limits<int>::max())
        {
          Alert::Exception() << name << " must be a nonnegative integer in range"
                             << Alert::Raise;
        }
        value = static_cast<std::remove_reference_t<decltype(value)>>(parsed);
      };
      const auto boolean = [&](const char* name, const char* description, bool& value) {
        PetscBool parsed = value ? PETSC_TRUE : PETSC_FALSE;
        ierr = PetscOptionsBool(name, description, nullptr, parsed, &parsed, nullptr);
        assert(ierr == PETSC_SUCCESS);
        value = parsed == PETSC_TRUE;
      };
      char mode[16] = "fluid";
      ierr = PetscOptionsString("-atrium_mode", "solid, fluid, or fsi", nullptr, mode,
        mode, sizeof(mode), nullptr);
      assert(ierr == PETSC_SUCCESS);
      if (std::string(mode) == "solid")
        out.mode = Mode::Solid;
      else if (std::string(mode) == "fluid")
        out.mode = Mode::Fluid;
      else if (std::string(mode) == "fsi")
        out.mode = Mode::FSI;
      else
        Alert::Exception() << "Unknown atrium mode: " << mode << Alert::Raise;
      char path[PETSC_MAX_PATH_LEN];
      ierr = PetscOptionsString("-atrium_mesh", "Combined MEDIT mesh in mm", nullptr,
        out.meshPath.c_str(), path, sizeof(path), nullptr);
      assert(ierr == PETSC_SUCCESS);
      out.meshPath = path;
      ierr = PetscOptionsString("-atrium_output", "Run output directory", nullptr,
        out.outputDirectory.c_str(), path, sizeof(path), nullptr);
      assert(ierr == PETSC_SUCCESS);
      out.outputDirectory = path;
      integer("-atrium_solid_material", "Solid volume material label", out.solidMaterial);
      integer("-atrium_fluid_material", "Fluid volume material label", out.fluidMaterial);
      real("-atrium_dt", "Common time step (s)", out.timeStep);
      integer(
        "-atrium_steps", "Physical time intervals, excluding initial state", out.steps);
      integer("-atrium_output_every",
        "Field export interval; zero keeps endpoints and fluid cycle boundaries",
        out.outputEvery);
      real("-atrium_solid_activation", "Electrical activation peak (1/s)",
        out.solid.activation);
      real("-atrium_solid_relaxation", "Resting electrical drive magnitude (1/s)",
        out.solid.relaxation);
      real("-atrium_solid_speed", "Plane-wave speed (mm/s)", out.solid.speed);
      real("-atrium_solid_width", "Activation band width (mm)", out.solid.width);
      real("-atrium_solid_matrix_scale", "matrix scale (kPa)",
        out.solid.passive.matrixScale);
      real("-atrium_solid_matrix_exponent", "matrix exponent (dimensionless)",
        out.solid.passive.matrixExponent);
      real(
        "-atrium_solid_fiber_scale", "fiber scale (kPa)", out.solid.passive.fiberScale);
      real("-atrium_solid_fiber_exponent", "fiber exponent (dimensionless)",
        out.solid.passive.fiberExponent);
      real("-atrium_solid_volume_penalty", "volume penalty (kPa)",
        out.solid.passive.volumePenalty);
      real(
        "-atrium_solid_viscosity", "Kelvin-Voigt viscosity (kPa s)", out.solid.viscosity);
      real("-atrium_solid_active_damping", "Contractile damping (kPa s)",
        out.solid.activeDamping);
      real("-atrium_solid_cross_bridge_stiffness", "cross bridge stiffness (kPa)",
        out.solid.crossBridgeStiffness);
      real("-atrium_solid_contractility", "contractility (kPa)", out.solid.contractility);
      real("-atrium_solid_destruction_rate", "destruction rate (dimensionless)",
        out.solid.destructionRate);
      real("-atrium_solid_density", "Wall reference density (kg/m^3)", out.solid.density);
      real("-atrium_solid_active_stiffness", "active stiffness (kPa)",
        out.solid.activeStiffness);
      real("-atrium_solid_apex_radius", "Appendage apex patch radius (mm)",
        out.solid.apexRadius);
      integer("-atrium_solid_base_label", "Basal boundary label", out.solid.baseLabel);
      boolean("-atrium_solid_cache_constitutive",
        "Cache constitutive evaluations within a Newton state",
        out.solid.cacheConstitutive);
      boolean("-atrium_solid_fibers_only",
        "Export prescribed fibers without advancing mechanics", out.solid.fibersOnly);
      real(
        "-atrium_fluid_mesh_scale", "Source coordinates to metres", out.fluid.meshScale);
      real("-atrium_fluid_pin", "Mean inlet pressure (Pa)", out.fluid.inletPressureMean);
      real("-atrium_fluid_pin_amplitude", "Inlet pressure amplitude (Pa)",
        out.fluid.inletPressureAmplitude);
      real("-atrium_fluid_outlet_pressure",
        "Outlet pressure without ventricular model (Pa)",
        out.fluid.prescribedOutletPressure);
      real("-atrium_fluid_outlet_resistance", "Outlet resistance (Pa s/m^3)",
        out.fluid.prescribedOutletResistance);
      real(
        "-atrium_fluid_vms_scale", "Orthogonal subscale multiplier", out.fluid.vmsScale);
      real("-atrium_fluid_graddiv_scale", "Grad-div multiplier", out.fluid.gradDivScale);
      real("-atrium_fluid_pspg_scale", "PSPG multiplier", out.fluid.pspgScale);
      real("-atrium_fluid_period", "Pressure and statistics cycle period (s)",
        out.fluid.period);
      real("-atrium_fluid_max_velocity", "Numerical divergence guard on speed (m/s)",
        out.fluid.maxVelocity);
      integer("-atrium_fluid_lag_preconditioner",
        "Numeric preconditioner refresh interval in physical steps",
        out.fluid.preconditionerLag);
      integer("-atrium_fluid_flow_cycles",
        "Cycles before optional species transport begins", out.fluid.flowCycles);
      boolean("-atrium_fluid_lv_model", "Enable 0D ventricular driver (fluid mode only)",
        out.fluid.useVentricularModel);
      boolean("-atrium_fluid_vms", "Enable orthogonal subscales (fluid mode only)",
        out.fluid.useVMS);
      boolean("-atrium_fluid_kinetics", "Enable coagulation transport (fluid mode only)",
        out.fluid.solveKinetics);
      boolean("-atrium_fsi_newton", "Use analytic matrix-free interface Newton",
        out.fsi.newton);
      real("-atrium_fsi_atol", "Absolute nodal interface residual tolerance (kPa)",
        out.fsi.absoluteTolerance);
      real("-atrium_fsi_solid_atol",
        "Absolute solid SNES tolerance during interface Newton",
        out.fsi.solidAbsoluteTolerance);
      boolean("-atrium_fsi_check_tangent",
        "Check analytic interface tangents with Taylor remainders", out.fsi.checkTangent);
      integer("-atrium_fsi_iterations", "Maximum Robin-Robin coupling iterations",
        out.fsi.maximumIterations);
      real("-atrium_fsi_robin",
        "Reference-area Robin impedance (Pa s/m); zero selects automatic scale",
        out.fsi.robinCoefficient);
      real("-atrium_fsi_rtol", "Relative solid displacement-change tolerance",
        out.fsi.relativeTolerance);
      PetscReal apex[3] = {out.solid.appendageApex(0), out.solid.appendageApex(1),
        out.solid.appendageApex(2)};
      PetscInt count = 3;
      PetscBool set = PETSC_FALSE;
      ierr = PetscOptionsRealArray("-atrium_solid_apex",
        "Appendage apex x,y,z (comma-separated, mm)", nullptr, apex, &count, &set);
      assert(ierr == PETSC_SUCCESS);
      if (set && count != 3)
      {
        Alert::Exception()
          << "-atrium_solid_apex requires three comma-separated coordinates"
          << Alert::Raise;
      }
      for (size_t i = 0; i < 3; ++i)
        out.solid.appendageApex(i) = apex[i];
      PetscOptionsEnd();
      PetscFunctionReturn(PETSC_SUCCESS);
    }();
    if (optionError != PETSC_SUCCESS)
    {
      Alert::Exception() << "Cannot parse atrium options (PETSc error " << optionError
                         << ")" << Alert::Raise;
    }
    PetscInt count = 0;
    char **names = nullptr, **values = nullptr;
    ierr = PetscOptionsLeftGet(nullptr, &count, &names, &values);
    assert(ierr == PETSC_SUCCESS);
    Optional<std::string> unknown;
    for (PetscInt i = 0; i < count; ++i)
    {
      const StringView name(names[i]);
      // Projection, coordinate, and ALE solvers consume these options later.
      const bool solver = name.starts_with("-atrium_sproj_") ||
        name.starts_with("-atrium_vproj_") || name.starts_with("-atrium_wss_") ||
        name.starts_with("-atrium_transmural_") || name.starts_with("-atrium_fibers_") ||
        name.starts_with("-atrium_ale_") || name.starts_with("-atrium_fsi_interface_");
      if (name.starts_with("-atrium_") && !solver && !unknown)
        unknown = std::string(name);
    }
    ierr = PetscOptionsLeftRestore(nullptr, &count, &names, &values);
    assert(ierr == PETSC_SUCCESS);
    if (unknown)
    {
      Alert::Exception() << "Unknown atrium option " << *unknown << "; use -help"
                         << Alert::Raise;
    }
    (void)ierr;
    return out.validated();
  }

  Parameters Parameters::validated() const
  {
    if ((mode != Mode::Solid && mode != Mode::Fluid && mode != Mode::FSI) ||
      !std::isfinite(timeStep) || timeStep <= 0 ||
      steps > std::numeric_limits<int>::max() || outputDirectory.empty() ||
      meshPath.empty())
      Alert::Exception() << "Invalid common atrium run parameters" << Alert::Raise;
    if (mode == Mode::FSI && solid.fibersOnly)
    {
      Alert::Exception() << "Fiber inspection is available only in solid mode"
                         << Alert::Raise;
    }
    if (mode != Mode::Fluid)
      validateParameters(solid);
    if (mode != Mode::Solid &&
      (!std::isfinite(fluid.maxVelocity) || fluid.maxVelocity <= 0))
    {
      Alert::Exception() << "-atrium_fluid_max_velocity must be finite and positive"
                         << Alert::Raise;
    }
    return *this;
  }
}
