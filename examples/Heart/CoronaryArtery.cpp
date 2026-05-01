/**
 * @file CoronaryArtery.cpp
 * @brief Driver for the coupled LV-0D / coronary 3D flow example.
 *
 * This executable configures and runs `CoupledLV0DCoronary3D` with default
 * paths for the example mesh and output files:
 * - Mesh: `../resources/examples/Heart/CoronaryArtery_Fluid.medit.mesh`
 * - XDMF basename: `CoronaryArtery`
 * - CSV output: `CoronaryArtery.csv`
 *
 * The driver accepts the following Rodin-specific PETSc option:
 * - `-coronary_flow_mode <newton|oseen>`
 *   Selects the 3D Navier-Stokes linearization strategy. The default is
 *   `newton`, which assembles the full Newton Jacobian for the convective
 *   term and Carreau-Yasuda viscosity and solves it with PETSc SNES. `oseen`
 *   assembles one lagged linear Oseen/Picard system per time step and solves
 *   it directly with PETSc KSP, without SNES nonlinear iterations.
 *
 * Unless the user overrides them on the command line, the executable installs
 * the following PETSc defaults:
 * - `-ksp_type preonly`
 * - `-pc_type lu`
 * - `-pc_factor_mat_solver_type mumps`
 * - `-mat_mumps_icntl_20 0`
 * - `-mat_mumps_icntl_21 0`
 *
 * The simulation defaults inherited from `CoupledLV0DCoronary3D::Config`
 * include `dt = 1e-3 s`, `nsteps = 2550`, `rho = 1060 kg/m^3`,
 * `eps = 1e-12`, `meshScale = 1e-3`, wall attribute `2`, inlet attribute
 * `3`, outlet attributes `4..9`, and default RCR parameters
 * `(Rp, C, Rd, pv0, pd0, par0) = (5e8, 5e-9, 1e9, 400, 10500, 11000)`.
 * The default Carreau-Yasuda blood model is
 * `(mu0, muInf, lambda, n, yasuda, gammaReg) =
 * (0.04868, 0.003605, 3.39, 0.198, 1.235, 1e-3)`.
 * The non-Newtonian outlet update uses proximal vessel
 * `(radius, length) = (0.004, 0.015)` and distal vessel
 * `(radius, length) = (0.0004, 0.002)`, with root-solve tolerances and
 * bracketing limits stored in `Config::outletFlowLaw`.
 * The 0D LV model defaults, initial conditions, activation waveform, and
 * atrial pressure waveform are stored in `Config::lv`, `Config::activation`,
 * and `Config::atrialPressure`.
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
#include <algorithm>
#include <cctype>
#include <exception>
#include <iostream>
#include <string>

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

      char flowMode[32] = {};
      PetscBool flowModeSet = PETSC_FALSE;
      PetscErrorCode ierr = PetscOptionsGetString(
          PETSC_NULLPTR,
          PETSC_NULLPTR,
          "-coronary_flow_mode",
          flowMode,
          sizeof(flowMode),
          &flowModeSet);
      assert(ierr == PETSC_SUCCESS);
      (void) ierr;

      if (flowModeSet)
      {
        std::string mode(flowMode);
        std::transform(
            mode.begin(), mode.end(), mode.begin(),
            [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

        if (mode == "newton")
        {
          cfg.flowMode =
            Rodin::Examples::Heart::CoupledLV0DCoronary3D::FlowMode::Newton;
        }
        else if (mode == "oseen")
        {
          cfg.flowMode =
            Rodin::Examples::Heart::CoupledLV0DCoronary3D::FlowMode::Oseen;
        }
        else
        {
          throw std::runtime_error(
              "Invalid -coronary_flow_mode. Expected newton or oseen.");
        }
      }

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
