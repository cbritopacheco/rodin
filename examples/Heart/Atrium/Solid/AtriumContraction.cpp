/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/// @brief Executable for the transient active atrial solid component.
#include <algorithm>
#include <cassert>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <exception>
#include <iostream>
#include "../Atrium.h"

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);
  int status = 0;
  {
    boost::mpi::environment environment(argc, argv);
    boost::mpi::communicator world(PETSC_COMM_WORLD, boost::mpi::comm_attach);
    Rodin::Context::MPI context(environment, world);
    try
    {
      namespace AtriumSolid = Rodin::Examples::Heart::AtriumContraction;
      const auto setDefault = [](const char* name, const char* value) {
        PetscBool set = PETSC_FALSE;
        PetscErrorCode ierr =
          PetscOptionsHasName(PETSC_NULLPTR, PETSC_NULLPTR, name, &set);
        assert(ierr == PETSC_SUCCESS);
        if (!set)
        {
          ierr = PetscOptionsSetValue(PETSC_NULLPTR, name, value);
          assert(ierr == PETSC_SUCCESS);
        }
        (void)ierr;
      };
      setDefault("-ksp_type", "gmres");
      setDefault("-pc_type", "asm");
      setDefault("-pc_asm_overlap", "2");
      setDefault("-sub_pc_type", "ilu");
      setDefault("-sub_pc_factor_levels", "1");
      setDefault("-ksp_gmres_restart", "100");
      setDefault("-ksp_rtol", "1e-10");
      setDefault("-ksp_diagonal_scale", "true");
      setDefault("-ksp_diagonal_scale_fix", "true");
      setDefault("-snes_linesearch_max_it", "60");
      using Atrium = Rodin::Examples::Heart::Atrium;
      Atrium::Parameters parameters;
      parameters.mode = Atrium::Mode::SolidOnly;
      parameters.solid = AtriumSolid::readOptions(argc, argv, parameters.meshPath);
      for (int i = 1; i < argc; ++i)
      {
        if (std::string(argv[i]).starts_with("--"))
        {
          const auto ierr = PetscOptionsClearValue(PETSC_NULLPTR, argv[i]);
          assert(ierr == PETSC_SUCCESS);
          (void)ierr;
        }
      }
      Atrium atrium(context, parameters);
      atrium.initialize();
      const size_t count = parameters.solid.stopAfter
        ? std::min(parameters.solid.steps, parameters.solid.stopAfter)
        : parameters.solid.steps;
      if (!parameters.solid.fibersOnly)
      {
        for (size_t step = 0; step < count; ++step)
          atrium.step();
      }
      atrium.finish();
    }
    catch (const std::exception& error)
    {
      if (world.rank() == 0)
        std::cerr << error.what() << std::endl;
      status = 1;
    }
  }
  PetscFinalize();
  return status;
}
