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
      using Atrium = Rodin::Examples::Heart::Atrium::Atrium;
      for (int i = 1; i < argc; ++i)
      {
        if (std::string(argv[i]).starts_with("--"))
        {
          Rodin::Alert::Exception()
            << "Unknown option " << argv[i]
            << "; use -help for the unified -atrium_* options" << Rodin::Alert::Raise;
        }
      }
      const auto parameters = Atrium::Parameters::fromOptions();
      PetscBool help = PETSC_FALSE;
      const auto ierr = PetscOptionsHasName(nullptr, nullptr, "-help", &help);
      assert(ierr == PETSC_SUCCESS);
      if (!help)
      {
        Atrium simulation(context, parameters);
        simulation.initialize();
        if (!(parameters.mode == Atrium::Mode::Solid && parameters.solid.fibersOnly))
        {
          for (size_t step = 0; step < parameters.steps; ++step)
            simulation.step();
        }
        simulation.finish();
      }
    }

    PetscFinalize();
    return status;
  }
  catch (const std::exception& e)
  {
    if (world.rank() == 0)
      std::cerr << "Fatal error: " << e.what() << "\n";
    PetscFinalize();
    return 1;
  }
}
