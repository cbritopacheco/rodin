/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <algorithm>
#include <cassert>
#include <exception>
#include <iostream>
#include <limits>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <Rodin/PETSc.h>
#include <Rodin/MPI.h>
#include "Atrium.h"

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  const auto setPETScDefault = [](const char* name, const char* value) {
    PetscBool set = PETSC_FALSE;
    PetscErrorCode ierr = PetscOptionsHasName(PETSC_NULLPTR, PETSC_NULLPTR, name, &set);
    if (ierr == PETSC_SUCCESS && !set)
      ierr = PetscOptionsSetValue(PETSC_NULLPTR, name, value);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
  };

  setPETScDefault("-ksp_type", "gmres");
  setPETScDefault("-pc_type", "asm");
  setPETScDefault("-pc_asm_overlap", "2");
  setPETScDefault("-sub_pc_type", "ilu");
  setPETScDefault("-sub_pc_factor_levels", "1");
  setPETScDefault("-ksp_gmres_restart", "100");
  setPETScDefault("-ksp_rtol", "1e-10");
  setPETScDefault("-ksp_diagonal_scale", "true");
  setPETScDefault("-ksp_diagonal_scale_fix", "true");

  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world(PETSC_COMM_WORLD, boost::mpi::comm_attach);
  Rodin::Context::MPI context(env, world);

  try
  {
    int status = 0;

    {
      Rodin::Examples::Heart::Atrium::Atrium::Parameters parameters;
      Rodin::Examples::Heart::Atrium::Fluid::Config cfg;

      char buffer[512];
      PetscBool got = PETSC_FALSE;
      PetscOptionsGetString(
        PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_mesh", buffer, sizeof(buffer), &got);
      if (got)
        parameters.meshPath = buffer;

      got = PETSC_FALSE;
      PetscOptionsGetString(
        PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_output", buffer, sizeof(buffer), &got);
      if (got)
      {
        cfg.xdmfBasename = buffer;
        cfg.csvPath = std::string(buffer) + ".csv";
      }
      PetscReal real = 0.0;
      got = PETSC_FALSE;
      PetscOptionsGetReal(
        PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_mesh_scale", &real, &got);
      if (got)
        cfg.meshScale = real;

      got = PETSC_FALSE;
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_dt", &real, &got);
      if (got)
        cfg.dt = real;

      got = PETSC_FALSE;
      PetscOptionsGetReal(
        PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_pin_amplitude", &real, &got);
      if (got)
        cfg.inletPressureAmplitude = real;

      got = PETSC_FALSE;
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_pin", &real, &got);
      if (got)
        cfg.inletPressureMean = real;

      PetscInt integer = 0;
      got = PETSC_FALSE;
      PetscOptionsGetInt(
        PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_fluid_lag_preconditioner", &integer, &got);
      if (got)
      {
        if (integer < 1 || integer > std::numeric_limits<int>::max())
        {
          Rodin::Alert::Exception()
            << "Fluid preconditioner lag must be a positive integer"
            << Rodin::Alert::Raise;
        }
        cfg.preconditionerLag = static_cast<int>(integer);
      }
      got = PETSC_FALSE;
      PetscOptionsGetInt(
        PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_flow_cycles", &integer, &got);
      if (got)
        cfg.flowCycles = static_cast<int>(integer);

      got = PETSC_FALSE;
      PetscOptionsGetInt(
        PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_species_cycles", &integer, &got);
      if (got)
        cfg.speciesCycles = static_cast<int>(integer);

      got = PETSC_FALSE;
      PetscOptionsGetInt(
        PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_output_every", &integer, &got);
      if (got)
        cfg.outputEvery = static_cast<int>(integer);

      got = PETSC_FALSE;
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_vms_scale", &real, &got);
      if (got)
        cfg.vmsScale = real;

      got = PETSC_FALSE;
      PetscOptionsGetReal(
        PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_graddiv_scale", &real, &got);
      if (got)
        cfg.gradDivScale = real;

      got = PETSC_FALSE;
      PetscOptionsGetReal(
        PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_pspg_scale", &real, &got);
      if (got)
        cfg.pspgScale = real;

      got = PETSC_FALSE;
      PetscOptionsGetReal(
        PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_outlet_pressure", &real, &got);
      if (got)
        cfg.prescribedOutletPressure = real;

      got = PETSC_FALSE;
      PetscOptionsGetReal(
        PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_outlet_resistance", &real, &got);
      if (got)
        cfg.prescribedOutletResistance = real;

      PetscBool flag = PETSC_FALSE;
      got = PETSC_FALSE;
      PetscOptionsGetBool(PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_lv_model", &flag, &got);
      if (got)
        cfg.useVentricularModel = (flag == PETSC_TRUE);

      got = PETSC_FALSE;
      PetscOptionsGetBool(PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_vms", &flag, &got);
      if (got)
        cfg.useVMS = (flag == PETSC_TRUE);

      got = PETSC_FALSE;
      PetscOptionsGetBool(PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_kinetics", &flag, &got);
      if (got)
        cfg.solveKinetics = (flag == PETSC_TRUE);

      using Atrium = Rodin::Examples::Heart::Atrium::Atrium;
      parameters.mode = Atrium::Mode::Fluid;
      parameters.fluid = cfg;
      parameters.solid = Rodin::Examples::Heart::Atrium::readOptions(
        argc, argv, parameters.meshPath);
      got = PETSC_FALSE;
      PetscOptionsGetString(
        PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_mode", buffer, sizeof(buffer), &got);
      if (got)
      {
        const std::string mode(buffer);
        if (mode == "solid")
          parameters.mode = Atrium::Mode::Solid;
        else if (mode == "fluid")
          parameters.mode = Atrium::Mode::Fluid;
        else if (mode == "fsi")
          parameters.mode = Atrium::Mode::FSI;
        else
        {
          Rodin::Alert::Exception()
            << "Unknown atrium mode " << mode << Rodin::Alert::Raise;
        }
      }
      if (parameters.mode == Atrium::Mode::FSI && parameters.solid.timeStep == 0)
        parameters.solid.timeStep = cfg.dt;
      for (int i = 1; i < argc; ++i)
      {
        if (std::string(argv[i]).starts_with("--"))
          PetscOptionsClearValue(PETSC_NULLPTR, argv[i]);
      }
      size_t count = parameters.mode == Atrium::Mode::Fluid
        ? static_cast<size_t>((cfg.flowCycles + cfg.speciesCycles) *
            static_cast<int>(cfg.period / cfg.dt + 0.5))
        : parameters.solid.steps;
      if (parameters.mode != Atrium::Mode::Fluid && parameters.solid.stopAfter)
        count = std::min(count, parameters.solid.stopAfter);
      integer = static_cast<PetscInt>(count);
      PetscOptionsGetInt(
        PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_steps", &integer, nullptr);
      if (integer < 0)
      {
        Rodin::Alert::Exception()
          << "Atrium step count must be nonnegative" << Rodin::Alert::Raise;
      }
      Atrium simulation(context, parameters);
      simulation.initialize();
      if (!(parameters.mode == Atrium::Mode::Solid && parameters.solid.fibersOnly))
      {
        for (PetscInt step = 0; step < integer; ++step)
          simulation.step();
      }
      simulation.finish();
    }

    PetscFinalize();
    return status;
  }
  catch (const std::exception& e)
  {
    std::cerr << "Fatal error: " << e.what() << "\n";
    PetscFinalize();
    return 1;
  }
}
