/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cassert>
#include <exception>
#include <iostream>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <Rodin/PETSc.h>
#include <Rodin/MPI.h>
#include "Fluid.h"

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

  setPETScDefault("-ksp_type", "preonly");
  setPETScDefault("-pc_type", "lu");
  setPETScDefault("-pc_factor_mat_solver_type", "mumps");
  setPETScDefault("-mat_mumps_icntl_20", "0");
  setPETScDefault("-mat_mumps_icntl_21", "0");

  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world(PETSC_COMM_WORLD, boost::mpi::comm_attach);
  Rodin::Context::MPI context(env, world);

  try
  {
    int status = 0;

    {
      Rodin::Examples::Heart::AtriumFluid::Config cfg;

      char buffer[512];
      PetscBool got = PETSC_FALSE;
      PetscOptionsGetString(
        PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_mesh", buffer, sizeof(buffer), &got);
      if (got)
        cfg.meshPath = buffer;

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
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_pin", &real, &got);
      if (got)
        cfg.inletPressureMean = real;

      PetscInt integer = 0;
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

      Rodin::Examples::Heart::AtriumFluid simulation(context, cfg);
      status = simulation.initialize().run();
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
