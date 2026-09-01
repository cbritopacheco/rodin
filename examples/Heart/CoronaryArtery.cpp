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
 *   `oseen`, which assembles one lagged linear Oseen/Picard system per time
 *   step and solves it directly with PETSc KSP, without SNES nonlinear
 *   iterations. `newton` assembles the full Newton Jacobian for the convective
 *   term and Carreau-Yasuda viscosity and solves it with PETSc SNES.
 * - `-coronary_outlet_backflow_stabilization <value>`
 *   Multiplies the outlet backflow damping term
 *   `<0.5 rho max(-(u_old.n), 0) u, v>`. The default is `1`; use `0` to
 *   disable this stabilization for diagnostics.
 * - `-coronary_inlet_backflow_stabilization <value>`
 *   Multiplies the inlet reversed-flow damping term
 *   `<0.5 rho max(u_old.n, 0) u, v>`. The default is `1`; the intended
 *   pressure-driven inflow has `u_old.n < 0` and is not damped by this term.
 * - `-coronary_dt <seconds>`
 *   Sets the nominal physical time step. The default is `1e-3`.
 * - `-coronary_nsteps <count>`
 *   Sets the number of accepted coupled time steps. The default is `2550`.
 * - `-coronary_time_adaptivity_reduction_factor <value>`
 *   Sets the factor applied to the 3D solver time step after a failed KSP/SNES
 *   solve. The default is `0.5`.
 * - `-coronary_time_adaptivity_max_levels <count>`
 *   Sets the maximum number of reductions attempted for one accepted time
 *   step. The default is `8`.
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
 * `eps = 1e-12`, `meshScale = 1e-2`, inlet/outlet backflow stabilization `1`,
 * wall attribute `2`, inlet attribute `3`, outlet attributes `4..9`.
 * Each outlet is a Starling resistor in an intramyocardial bed with a single
 * state, the microvascular transmural pressure `p_tm = p_c - p_im`; its three
 * lumped constants `(R_a, R_v, C)` are produced by the calibration from
 * `lcaTargetFlow`, `venularPressureFraction` and `coronaryComplianceTotal`.
 * The arteriolar resistance is assembled implicitly on the outlet boundary,
 * so the 3D-0D coupling is unconditionally stable in `dt`. See
 * `CoronaryArtery/RCR_formulacion_minima.tex`.
 * The 3D solve has local time-step adaptivity enabled by default: if the
 * PETSc KSP/SNES solve fails, the 3D flow state is restored and retried with
 * solver `dt *= 0.5`, down to 8 reductions. After accepted reduced solves,
 * the next step grows by `1 / 0.5` until the original `dt` is recovered. 0D
 * Newton failures are not retried by this mechanism.
 * The default Carreau-Yasuda blood model is
 * `(mu0, muInf, lambda, n, yasuda, gammaReg) =
 * (0.04868, 0.003605, 3.39, 0.198, 1.235, 1e-3)`.
 * The non-Newtonian outlet update evaluates the WRMS closure through the
 * universal apparent-viscosity table `mu_ap(tau_w) = tau_w^4 / (4 I(tau_w))`,
 * built once at calibration; its bounds, node count and the outlet Newton
 * tolerances are stored in `Config::outletFlowLaw`. The rheological operating
 * point of each limb is set by `arteriolarShearRate` and `venularShearRate`,
 * and the modulation is normalized by `newtonianCalibrationViscosity`, never
 * by the running rheology, so a change of blood properties moves the flow.
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
      cfg.meshPath = "../resources/examples/Heart/salida_3d.mesh";
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

      // ---- Pathological-scenario knobs -----------------------------------
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

      // Pressure-overload / subendocardial territory: transmission fraction
      // p_im = alpha p_LV. Baseline 0.7; subendocardium up to ~0.9.
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

      // Sensitivity knob: total microvascular compliance (m^3/Pa). Changes the
      // retrograde spike, hardly the mean flow; identified range 1e-10..3e-9.
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
      // reduced outlets lose their shear dependence while the 3D Carreau-Yasuda
      // field is left untouched, so the pair of runs isolates the closure.
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
