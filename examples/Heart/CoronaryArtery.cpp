/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file CoronaryArtery.cpp
 * @brief Driver of the coupled 0D left ventricle / 3D coronary flow example.
 *
 * Configures and runs CoupledLV0DCoronary3D. The mesh, the output basename
 * and every physical default are taken from CoupledLV0DCoronary3D::Config and
 * may be overridden through the PETSc option database:
 *
 * - `-coronary_flow_mode <newton|oseen>`
 *   3D linearization. `oseen` (default) assembles one lagged linear
 *   Oseen/Picard system per time step and solves it with PETSc KSP; `newton`
 *   assembles the full Jacobian of the convective term and of the
 *   Carreau-Yasuda viscosity and solves it with PETSc SNES.
 * - `-coronary_dt <seconds>`
 *   Nominal physical time step.
 * - `-coronary_nsteps <count>`
 *   Number of accepted coupled time steps.
 * - `-coronary_time_adaptivity_reduction_factor <value>` in (0, 1)
 *   Factor applied to the 3D time step after a failed KSP/SNES solve.
 * - `-coronary_time_adaptivity_max_levels <count>`
 *   Maximum number of reductions attempted for one accepted step.
 * - `-coronary_outlet_backflow_stabilization <value>`
 *   Multiplies the outlet damping term <0.5 rho max(-(u_old.n), 0) u, v>.
 *   0 disables it.
 * - `-coronary_inlet_backflow_stabilization <value>`
 *   Multiplies the inlet damping term <0.5 rho max(u_old.n, 0) u, v>. The
 *   intended pressure-driven inflow has u_old.n < 0 and is not damped.
 * - `-coronary_operating_pra <Pa>`
 *   Right atrial pressure seen by the running outlets, with the calibration
 *   left at the healthy baseline (venous-hypertension scenarios).
 * - `-coronary_alpha_im <value>` in [0, 1]
 *   Intramyocardial transmission fraction, p_im = alpha p_LV.
 * - `-coronary_compliance_total <m^3/Pa>`
 *   Total microvascular compliance, split across outlets by the Murray weight.
 * - `-coronary_constant_outlet <bool>`
 *   Freezes the reduced outlet closure at its high-shear plateau, giving a
 *   constant outlet resistance while the 3D field stays non-Newtonian.
 * - `-coronary_output_prefix <dir>`
 *   Writes <dir>/CoronaryArtery.{xdmf,csv}.
 *
 * Unless overridden on the command line, the executable installs the PETSc
 * defaults `-ksp_type preonly`, `-pc_type lu`,
 * `-pc_factor_mat_solver_type mumps`, `-mat_mumps_icntl_20 0` and
 * `-mat_mumps_icntl_21 0`.
 *
 * The 3D solve uses local time-step adaptivity: on a failed KSP/SNES solve
 * the 3D state is restored and retried with a reduced dt, and the step grows
 * back to the nominal value once solves are accepted again. Failures of the
 * 0D Newton are not retried by this mechanism.
 *
 * Example:
 *
 *   VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
 *   mpirun -n 8 ./examples/Heart/CoronaryArtery \
 *     -snes_atol 1e-8 -snes_rtol 1e-8 -snes_stol 1e-10 \
 *     -mat_mumps_icntl_7 7 \
 *     -ksp_converged_reason -snes_monitor -snes_converged_reason -ksp_monitor
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
      Rodin::Examples::Heart::CoupledLV0DCoronary3D::Config cfg;
      cfg.meshPath = "../resources/examples/Heart/coronaria3d.mesh";
      cfg.xdmfBasename = "hyp2/CoronaryArtery";
      cfg.csvPath = "hyp2/CoronaryArtery.csv";

      char flowMode[32] = {};
      PetscBool flowModeSet = PETSC_FALSE;
      PetscErrorCode ierr = PetscOptionsGetString(PETSC_NULLPTR, PETSC_NULLPTR,
        "-coronary_flow_mode", flowMode, sizeof(flowMode), &flowModeSet);
      assert(ierr == PETSC_SUCCESS);
      (void)ierr;

      if (flowModeSet)
      {
        std::string mode(flowMode);
        std::transform(mode.begin(), mode.end(), mode.begin(),
          [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

        if (mode == "newton")
        {
          cfg.flowMode = Rodin::Examples::Heart::CoupledLV0DCoronary3D::FlowMode::Newton;
        }
        else if (mode == "oseen")
        {
          cfg.flowMode = Rodin::Examples::Heart::CoupledLV0DCoronary3D::FlowMode::Oseen;
        }
        else
        {
          throw std::runtime_error(
            "Invalid -coronary_flow_mode. Expected newton or oseen.");
        }
      }

      PetscReal backflowStabilization = cfg.outletBackflowStabilization;
      PetscBool backflowStabilizationSet = PETSC_FALSE;
      ierr = PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR,
        "-coronary_outlet_backflow_stabilization", &backflowStabilization,
        &backflowStabilizationSet);
      assert(ierr == PETSC_SUCCESS);
      if (backflowStabilizationSet)
        cfg.outletBackflowStabilization = backflowStabilization;

      PetscReal inletBackflowStabilization = cfg.inletBackflowStabilization;
      PetscBool inletBackflowStabilizationSet = PETSC_FALSE;
      ierr = PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR,
        "-coronary_inlet_backflow_stabilization", &inletBackflowStabilization,
        &inletBackflowStabilizationSet);
      assert(ierr == PETSC_SUCCESS);
      if (inletBackflowStabilizationSet)
        cfg.inletBackflowStabilization = inletBackflowStabilization;

      PetscReal dt = cfg.dt;
      PetscBool dtSet = PETSC_FALSE;
      ierr =
        PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_dt", &dt, &dtSet);
      assert(ierr == PETSC_SUCCESS);
      if (dtSet)
      {
        if (dt <= 0)
          throw std::runtime_error("-coronary_dt must be positive.");
        cfg.dt = dt;
      }

      PetscInt nsteps = static_cast<PetscInt>(cfg.nsteps);
      PetscBool nstepsSet = PETSC_FALSE;
      ierr = PetscOptionsGetInt(
        PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_nsteps", &nsteps, &nstepsSet);
      assert(ierr == PETSC_SUCCESS);
      if (nstepsSet)
      {
        if (nsteps < 0)
          throw std::runtime_error("-coronary_nsteps must be nonnegative.");
        cfg.nsteps = static_cast<size_t>(nsteps);
      }

      PetscReal reductionFactor = cfg.timeAdaptivityReductionFactor;
      PetscBool reductionFactorSet = PETSC_FALSE;
      ierr = PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR,
        "-coronary_time_adaptivity_reduction_factor", &reductionFactor,
        &reductionFactorSet);
      assert(ierr == PETSC_SUCCESS);
      if (reductionFactorSet)
      {
        if (reductionFactor <= 0 || reductionFactor >= 1)
        {
          throw std::runtime_error(
            "-coronary_time_adaptivity_reduction_factor must be in (0, 1).");
        }
        cfg.timeAdaptivityReductionFactor = reductionFactor;
      }

      PetscInt maxAdaptivityLevels = cfg.timeAdaptivityMaxLevels;
      PetscBool maxAdaptivityLevelsSet = PETSC_FALSE;
      ierr = PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR,
        "-coronary_time_adaptivity_max_levels", &maxAdaptivityLevels,
        &maxAdaptivityLevelsSet);
      assert(ierr == PETSC_SUCCESS);
      if (maxAdaptivityLevelsSet)
      {
        if (maxAdaptivityLevels < 0)
        {
          throw std::runtime_error(
            "-coronary_time_adaptivity_max_levels must be nonnegative.");
        }
        cfg.timeAdaptivityMaxLevels = maxAdaptivityLevels;
      }

      // ---- Scenario parameters -------------------------------------------
      // Venous hypertension: runtime drainage pressure of the outlets, with
      // the calibration left at the healthy baseline (frozen bed geometry).
      PetscReal operatingPra = cfg.operatingRightAtrialPressure;
      PetscBool operatingPraSet = PETSC_FALSE;
      ierr = PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR,
        "-coronary_operating_pra", &operatingPra, &operatingPraSet);
      assert(ierr == PETSC_SUCCESS);
      if (operatingPraSet)
      {
        if (operatingPra <= 0)
          throw std::runtime_error("-coronary_operating_pra must be positive (Pa).");
        cfg.operatingRightAtrialPressure = operatingPra;
      }

      // Pressure overload or subendocardial territory: transmission fraction
      // p_im = alpha p_LV. Baseline 0.7, subendocardium up to ~0.9.
      PetscReal alphaIm = cfg.intramyocardialFraction;
      PetscBool alphaImSet = PETSC_FALSE;
      ierr = PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR,
        "-coronary_alpha_im", &alphaIm, &alphaImSet);
      assert(ierr == PETSC_SUCCESS);
      if (alphaImSet)
      {
        if (alphaIm < 0 || alphaIm > 1)
          throw std::runtime_error("-coronary_alpha_im must be in [0, 1].");
        cfg.intramyocardialFraction = alphaIm;
      }

      // Total microvascular compliance (m^3/Pa). Moves the retrograde spike
      // and hardly the mean flow; identified range 1e-10 to 3e-9.
      PetscReal complianceTotal = cfg.coronaryComplianceTotal;
      PetscBool complianceTotalSet = PETSC_FALSE;
      ierr = PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR,
        "-coronary_compliance_total", &complianceTotal, &complianceTotalSet);
      assert(ierr == PETSC_SUCCESS);
      if (complianceTotalSet)
      {
        if (complianceTotal <= 0)
        {
          throw std::runtime_error(
            "-coronary_compliance_total must be positive (m^3/Pa).");
        }
        cfg.coronaryComplianceTotal = complianceTotal;
      }

      // Constant-resistance outlet closure (the mu_inf comparison run): the
      // reduced outlets lose their shear dependence while the 3D field is left
      // untouched, so the pair of runs isolates the closure.
      PetscBool constantOutlet = PETSC_FALSE;
      PetscBool constantOutletSet = PETSC_FALSE;
      ierr = PetscOptionsGetBool(PETSC_NULLPTR, PETSC_NULLPTR,
        "-coronary_constant_outlet", &constantOutlet, &constantOutletSet);
      assert(ierr == PETSC_SUCCESS);
      if (constantOutletSet)
        cfg.constantOutletResistance = (constantOutlet == PETSC_TRUE);

      // Per-scenario output directory: <prefix>/CoronaryArtery.{xdmf,csv}.
      char outPrefix[256] = {};
      PetscBool outPrefixSet = PETSC_FALSE;
      ierr = PetscOptionsGetString(PETSC_NULLPTR, PETSC_NULLPTR,
        "-coronary_output_prefix", outPrefix, sizeof(outPrefix), &outPrefixSet);
      assert(ierr == PETSC_SUCCESS);
      if (outPrefixSet)
      {
        cfg.xdmfBasename = std::string(outPrefix) + "/CoronaryArtery";
        cfg.csvPath = std::string(outPrefix) + "/CoronaryArtery.csv";
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
