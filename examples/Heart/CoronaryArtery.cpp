/**
 * @file CoronaryArtery.cpp
 * @brief Driver for the coupled LV-0D / coronary 3D Newtonian flow example.
 *
 * This executable configures and runs `CoupledLV0DCoronary3D` with default
 * paths for the example mesh and output files:
 * - Mesh: `../resources/examples/Heart/CoronaryArtery_Fluid.medit.mesh`
 * - XDMF basename: `CoronaryArtery`
 * - CSV output: `CoronaryArtery.csv`
 *
 * VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
 * > mpirun -n 8 ./examples/Heart/CoronaryArtery \
 * >   -snes_atol 1e-8 \
 * >   -snes_rtol 1e-8 \
 * >   -snes_stol 1e-10 \
 * >   -mat_mumps_icntl_7 7 \
 * >   -ksp_converged_reason -snes_monitor -snes_converged_reason -ksp_monitor
 */
#include <cassert>
#include <exception>
#include <iostream>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <petscsys.h>

#include <Rodin/MPI.h>

#include "CoronaryArtery/CoupledLV0DCoronary3D.h"

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  const auto setPETScDefault =
    [](const char* name, const char* value)
    {
      PetscBool set = PETSC_FALSE;
      PetscErrorCode ierr = PetscOptionsHasName(PETSC_NULLPTR, PETSC_NULLPTR, name, &set);
      if (ierr == PETSC_SUCCESS && !set)
        ierr = PetscOptionsSetValue(PETSC_NULLPTR, name, value);
      assert(ierr == PETSC_SUCCESS);
      (void) ierr;
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
      Rodin::Examples::Heart::CoupledLV0DCoronary3D::Config cfg;
      cfg.meshPath = "../resources/examples/Heart/CoronaryArtery_Fluid.medit.mesh";
      cfg.xdmfBasename = "CoronaryArtery";
      cfg.csvPath = "CoronaryArtery.csv";

      Rodin::Examples::Heart::CoupledLV0DCoronary3D simulation(context, cfg);
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
