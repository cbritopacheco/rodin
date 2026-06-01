// CoupledLV0DCoronary3D.h
#ifndef EXAMPLES_HEART_CORONARYARTERY_COUPLEDLV0DCORONARY3D_H
#define EXAMPLES_HEART_CORONARYARTERY_COUPLEDLV0DCORONARY3D_H

#include <array>
#include <fstream>
#include <map>
#include <string>

#include "Rodin/Heart/CCMLC2014.h"
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/MPI.h>
#include <Rodin/PETSc.h>
#include <Rodin/Solver.h>
#include <Rodin/Types.h>
#include <Rodin/Variational.h>

#include "CoronaryArteryTiming.h"
#include "VMSConvectionIntegrator.h"

namespace Rodin::Examples::Heart {
class CoupledLV0DCoronary3D {
public:
  using Real = Rodin::Real;
  using Model = Rodin::Heart::CCMLC2014T<>;
  using Attribute = Rodin::Geometry::Attribute;
  using MeshType = Rodin::Geometry::Mesh<Rodin::Context::MPI>;

  using VelocityFESType =
      Rodin::Variational::H1<2, Rodin::Math::SpatialVector<Real>, MeshType>;

  using PressureFESType = Rodin::Variational::H1<1, Real, MeshType>;

  using VelocityGridFunctionType =
      Rodin::PETSc::Variational::GridFunction<VelocityFESType>;

  using PressureGridFunctionType =
      Rodin::PETSc::Variational::GridFunction<PressureFESType>;

  using VelocityTrialFunctionType =
      Rodin::PETSc::Variational::TrialFunction<VelocityGridFunctionType,
                                               VelocityFESType>;

  using VelocityTestFunctionType =
      Rodin::PETSc::Variational::TestFunction<VelocityFESType>;

  using PressureTrialFunctionType =
      Rodin::PETSc::Variational::TrialFunction<PressureGridFunctionType,
                                               PressureFESType>;

  using PressureTestFunctionType =
      Rodin::PETSc::Variational::TestFunction<PressureFESType>;

  using FluxLinearFormType =
      Rodin::Variational::LinearForm<PressureFESType, ::Vec>;

  using LinearSystemType = Rodin::PETSc::Math::LinearSystem;

  using FlowProblemType = Rodin::Variational::Problem<
      LinearSystemType, VelocityTrialFunctionType, PressureTrialFunctionType,
      VelocityTestFunctionType, PressureTestFunctionType>;

  using ScalarProjectionProblemType =
      Rodin::Variational::Problem<LinearSystemType, PressureTrialFunctionType,
                                  PressureTestFunctionType>;

  using TauFESType = Rodin::Variational::P1<Real, MeshType>;

  using TauGridFunctionType =
      Rodin::PETSc::Variational::GridFunction<TauFESType>;

  using TauTrialFunctionType =
      Rodin::PETSc::Variational::TrialFunction<TauGridFunctionType, TauFESType>;

  using TauTestFunctionType =
      Rodin::PETSc::Variational::TestFunction<TauFESType>;

  using TauProblemType =
      Rodin::Variational::Problem<LinearSystemType, TauTrialFunctionType,
                                  TauTestFunctionType>;

  using VMSFESType =
      Rodin::Variational::H1<2, Rodin::Math::SpatialVector<Real>, MeshType>;

  using VMSGridFunctionType =
      Rodin::PETSc::Variational::GridFunction<VMSFESType>;

  using VMSTrialFunctionType =
      Rodin::PETSc::Variational::TrialFunction<VMSGridFunctionType, VMSFESType>;

  using VMSTestFunctionType =
      Rodin::PETSc::Variational::TestFunction<VMSFESType>;

  using VMSProblemType =
      Rodin::Variational::Problem<LinearSystemType, VMSTrialFunctionType,
                                  VMSTestFunctionType>;

  struct RCR {
    /// @brief Proximal resistance.
    Real Rp = 0.0;
    /// @brief Compliance.
    Real C = 0.0;
    /// @brief Distal resistance.
    Real Rd = 0.0;
    /// @brief Distal pressure state.
    Real pd = 0.0;
    /// @brief Capacitor pressure state.
    Real pc = 0.0;
    /// @brief Outlet pressure applied to the 3D model.
    Real pout = 0.0;
    /// @brief Distal branch flow leaving the capacitor.
    Real qd = 0.0;
  };

  /**
   * @brief Carreau-Yasuda blood viscosity parameters.
   */
  struct CarreauYasuda {
    /// @brief Low-shear viscosity.
    Real mu0 = 0.104;
    /// @brief Infinite-shear viscosity.
    Real muInf = 0.00536;
    /// @brief Relaxation time.
    Real lambda = 11.1048;
    /// @brief Power-law index.
    Real n = 0.1502;
    /// @brief Yasuda transition exponent.
    Real yasuda = 0.8625;
    /// @brief Shear-rate regularization used in the 3D viscosity.
    Real gammaRegularization = 1.0e-3;
  };

  /**
   * @brief Geometry and nonlinear-solve parameters for the coronary
   *        outlet flow law used in the non-Newtonian RCR update.
   */
  struct OutletFlowLaw {
    /// @brief Proximal surrogate vessel radius.
    Real proximalRadius = 6.e-4;
    /// @brief Proximal surrogate vessel length.
    Real proximalLength = 0.00075;
    /// @brief Distal surrogate vessel radius.
    Real distalRadius = 1e-4;
    /// @brief Distal surrogate vessel length.
    Real distalLength = 0.0025;
    /// @brief Pressure-drop threshold for the Poiseuille fallback.
    Real pressureDropTolerance = 1.0e-12;
    /// @brief Minimum shear-rate bracket.
    Real minShearRate = 1.0e-8;
    /// @brief Number of RK4 substeps for the WRMS flow integral.
    int integralSteps = 100;
    /// @brief Maximum bracketing expansions for outlet scalar solves.
    int maxBracketIterations = 100;
    /// @brief Wall shear root solver absolute tolerance.
    Real shearAbsoluteTolerance = 1.0e-12;
    /// @brief Wall shear root solver relative tolerance.
    Real shearRelativeTolerance = 1.0e-10;
    /// @brief Wall shear root solver step tolerance.
    Real shearStepTolerance = 1.0e-12;
    /// @brief Wall shear root solver maximum iterations.
    int shearMaxIterations = 50;
    /// @brief Flow inversion root solver absolute tolerance.
    Real flowAbsoluteTolerance = 1.0e-10;
    /// @brief Flow inversion root solver relative tolerance.
    Real flowRelativeTolerance = 1.0e-9;
    /// @brief Flow inversion root solver step tolerance.
    Real flowStepTolerance = 1.0e-12;
    /// @brief Flow inversion root solver maximum iterations.
    int flowMaxIterations = 50;
    /// @brief Flow magnitude treated as zero in pressure-drop inversion.
    Real zeroFlowTolerance = 1.0e-16;
    /// @brief Minimum pressure-drop bracket.
    Real pressureDropBracketMin = 1.0;
    /// @brief Distal capacitor bracket pressure pad.
    Real distalPressureBracketPad = 1000.0;
  };

  /**
   * @brief Piecewise-linear LV activation waveform parameters.
   */
  struct Activation {
    /// @brief Period of the prescribed cardiac cycle.
    Real period = 0.85;
    /// @brief Activation ramp start.
    Real tRampStart = 0.15 - 0.05;
    /// @brief Activation ramp end.
    Real tRampEnd = 0.21 - 0.05;
    /// @brief Activation plateau end.
    Real tPlateauEnd = 0.36 - 0.05;
    /// @brief Relaxation ramp end.
    Real tRelaxEnd = 0.45 - 0.05;
    /// @brief Negative activation plateau end.
    Real tNegativeEnd = 0.6 - 0.05;
    /// @brief Positive activation plateau value.
    Real positiveValue = 35.0;
    /// @brief Negative activation plateau value.
    Real negativeValue = -20.0;
  };

  /**
   * @brief Piecewise-linear atrial pressure waveform parameters.
   */
  struct AtrialPressure {
    /// @brief Period of the prescribed cardiac cycle.
    Real period = 0.85;
    /// @brief Baseline atrial pressure.
    Real minValue = 500.0;
    /// @brief First plateau atrial pressure.
    Real maxValue = 1000.0;
    /// @brief Second plateau atrial pressure.
    Real secondThreshold = 1250.0;
    /// @brief End of first down-ramp.
    Real t1 = 0.02;
    /// @brief End of first plateau.
    Real t2 = 0.15;
    /// @brief End of first return ramp.
    Real t3 = 0.17;
    /// @brief End of second up-ramp.
    Real t4 = 0.56;
    /// @brief End of second plateau.
    Real t5 = 0.62;
    /// @brief End of cycle return ramp.
    Real t6 = 0.85;
  };

  /**
   * @brief Parameters passed to the 0D LV model.
   */
  struct LVModel {
    /// @brief 0D fluid density.
    Real rho = 1.0e3;
    /// @brief Reference radius.
    Real R0 = 2.45e-2;
    /// @brief Reference wall thickness.
    Real d0 = 1.4e-2;
    /// @brief Passive elastic stiffness.

    Real Es = 3.0e6;
    /// @brief Viscous parameter.
    Real mu = 70.0;
    /// @brief Viscous parameter.
    Real eta = 70.0;
    /// @brief Active stress gain.
    Real alpha = 1.5;
    /// @brief Load-dependent relaxation time scale.
    Real alphaR = 0.12;
    /// @brief Active stiffness scale.
    Real k0 = 1.0e5;
    /// @brief Active stress scale.
    Real sigma0 = 1.25e5;
    /// @brief Proximal arterial resistance.
    Real Rp = 5e7;
    /// @brief Proximal arterial compliance.
    Real Cp = 5e-9;
    /// @brief Distal arterial resistance.
    Real Rd = 1.0e8;
    /// @brief Distal arterial compliance.
    Real Cd = 5.0e-10;
    /// @brief Atrial valve coefficient.
    Real Kat = 4.0e-6;
    /// @brief Peripheral valve coefficient.
    Real Kp = 5.0e-10;
    /// @brief Arterial valve coefficient.
    Real Kar = 2.e-7;
    /// @brief LV cavity capacity.
    Real cavityCapacity = 5.0e-12;
    /// @brief Local 0D Newton absolute tolerance.
    Real localTolerance = 1.0e-12;
    /// @brief Local 0D Newton maximum iterations.
    int localMaxIterations = 50;
    /// @brief Local 0D Newton damping.
    Real localDamping = 1.0;
    /// @brief Absolute-value regularization.
    Real absRegularization = 1.0e-14;
    /// @brief Initial fiber deformation.
    Real initFibDef = 0.0;
    /// @brief Initial active stiffness.
    Real initActiveStiffness = 0.0;
    /// @brief Initial active stress.
    Real initActiveStress = 0.0;
    /// @brief Low-fiber-deformation target for load-dependent relaxation.
    Real relaxationM0Low = 1.6;
    /// @brief High-fiber-deformation target for load-dependent relaxation.
    Real relaxationM0High = 1.0;
    /// @brief Fiber deformation at which m0 reaches relaxationM0Low.
    Real relaxationM0LowEc = 0.0;
    /// @brief Fiber deformation at which m0 reaches relaxationM0High.
    Real relaxationM0HighEc = 2.0;
    /// @brief Systemic venous pressure callback value.
    Real systemicVenousPressure = 1.0e3;
    /// @brief Passive energy parameter mu1.
    Real passiveMu1 = 0.0;
    /// @brief Passive energy parameter mu2.
    Real passiveMu2 = 0.0;
    /// @brief Passive energy parameter C0.
    Real passiveC0 = 1.9e3;
    /// @brief Passive energy parameter C1.
    Real passiveC1 = 1.1e-1;
    /// @brief Passive energy parameter C2.
    Real passiveC2 = 1.9e3;
    /// @brief Passive energy parameter C3.
    Real passiveC3 = 1.1e-1;
    /// @brief 0D model maximum Newton iterations.
    int maxIterations = 200;
    /// @brief 0D model absolute tolerance.
    Real absoluteTolerance = 1.0e-8;
    /// @brief 0D model relative tolerance.
    Real relativeTolerance = 1.0e-8;
    /// @brief 0D model step tolerance.
    Real stepTolerance = 1.0e-10;
    /// @brief 0D model damping factor.
    Real dampingFactor = 1.0;
    /// @brief Initial LV fiber deformation state.
    Real initialY = 0.0;
    /// @brief Initial LV velocity state.
    Real initialV = 0.0;
    /// @brief Offset applied to atrial pressure to initialize pv.
    Real initialPvOffset = -100.0;
    /// @brief Initial arterial pressure.
    Real initialPar = 11000.0;
    /// @brief Initial distal pressure.
    Real initialPd = 10000.0;
    /// @brief Low-shear viscosity.
    Real mu_0 = 0.104;
    /// @brief Infinite-shear viscosity.
    Real mu_Inf = 0.00536;
    /// @brief Relaxation time.
    Real lambda = 11.1048;
    /// @brief Power-law index.
    Real n = 0.1502;
    /// @brief Yasuda transition exponent.
    Real yasuda = 0.8625;
    /// @brief Proximal surrogate vessel radius.
    Real proximalRadius = 0.0125;
    /// @brief Proximal surrogate vessel length.
    Real proximalLength = 0.4;
    /// @brief Distal surrogate vessel radius.
    Real distalRadius = 0.005;
    /// @brief Distal surrogate vessel length.
    Real distalLength = 0.2;
  };

  /**
   * @brief Linearization strategy for the 3D coronary flow solve.
   */
  enum class FlowMode {
    /**
     * @brief Full Newton linearization of the nonlinear 3D flow residual.
     *
     * The convective term is differentiated with respect to the current
     * velocity iterate, and the Carreau-Yasuda viscosity includes its
     * directional derivative in the Jacobian.
     */
    Newton,

    /**
     * @brief Oseen/Picard linearization with lagged transport coefficients.
     *
     * The advecting velocity, Temam divergence factor, and
     * Carreau-Yasuda viscosity are evaluated using the previous time-step
     * velocity. This gives a lagged linearized 3D flow system that is
     * solved directly with PETSc KSP, without PETSc SNES nonlinear
     * iterations.
     */
    Oseen
  };

  struct Config {
    /// @brief Input coronary fluid mesh path.
    std::string meshPath = "CoronaryArtery.mesh";
    /// @brief Basename for XDMF and related output files.
    std::string xdmfBasename = "CoronaryArtery";
    /// @brief CSV diagnostics output path.
    std::string csvPath = "CoronaryArtery.csv";

    /// @brief No-slip wall boundary attribute.
    Attribute wall = 2;
    /// @brief Inlet boundary attribute.
    Attribute inlet = 4;
    /// @brief Outlet boundary attributes, in the same order used by RCR data.
    std::array<Attribute, 6> outlets{{7, 8, 9, 10, 14,15}};

    /// @brief Mesh coordinate scale applied after partitioning.
    Real meshScale = 1.0e-3;
    /// @brief Pressure stabilization parameter.
    Real eps = 1.0e-12;
    /// @brief 3D blood density.
    Real rho = 1060.0;
    /// @brief Inlet reversed-flow damping multiplier. Set to 0 to disable.
    Real inletBackflowStabilization = 1.0;
    /// @brief Inlet normal impedance coefficient in Pa s / m. Set to 0 to
    /// disable.
    /// @details Defaults to defaultRCR.Rp times the scaled inlet area.
    Real inletImpedance = 5.e1;
    /// @brief Outlet backflow damping multiplier. Set to 0 to disable.
    Real outletBackflowStabilization = 1.0;

    /// @brief Time-step size.
    Real dt = 1.0e-3;
    /// @brief Number of time steps.
    size_t nsteps = 3 * static_cast<int>(0.85 / 1.0e-3);
    /// @brief Factor applied to dt when the 3D KSP/SNES solve fails.
    Real timeAdaptivityReductionFactor = 0.5;
    /// @brief Maximum number of successive dt reductions per accepted step.
    int timeAdaptivityMaxLevels = 8;

    /// @brief 3D coronary flow linearization mode. Defaults to Oseen/Picard.
    FlowMode flowMode = FlowMode::Oseen;
    /// @brief Blood viscosity model shared by 3D flow and outlet laws.
    CarreauYasuda viscosity;
    /// @brief Non-Newtonian outlet flow-law parameters.
    OutletFlowLaw outletFlowLaw;
    /// @brief Prescribed LV activation waveform parameters.
    Activation activation;
    /// @brief Prescribed atrial pressure waveform parameters.
    AtrialPressure atrialPressure;
    /// @brief 0D LV model parameters and initial conditions.
    LVModel lv;
    /// @brief Default RCR parameters copied to every outlet at startup.
    RCR defaultRCR{5.0e8, 5.0e-11, 1.0e9, 500.0, 10400.0, 10800.0};

    Real inletTangentialDamping = 1e3;
    Real inletVelocityDamping = 0.0;
  };

  explicit CoupledLV0DCoronary3D(const Rodin::Context::MPI &context);
  CoupledLV0DCoronary3D(const Rodin::Context::MPI &context, const Config &cfg);
  ~CoupledLV0DCoronary3D();

  CoupledLV0DCoronary3D(const CoupledLV0DCoronary3D &) = delete;
  CoupledLV0DCoronary3D &operator=(const CoupledLV0DCoronary3D &) = delete;

  int run();
  CoupledLV0DCoronary3D &initialize();

  Config &getConfig() noexcept { return m_cfg; }
  const Config &getConfig() const noexcept { return m_cfg; }

  Model &getModel() noexcept { return m_model; }
  const Model &getModel() const noexcept { return m_model; }

private:
  struct StepData {
    Real t = 0.0;
    Real pat = 0.0;
    Real psv = 0.0;

    Real y = 0.0;
    Real v = 0.0;
    Real radius = 0.0;
    Real lvVolume = 0.0;
    Real lvFlow = 0.0;
    Real pv = 0.0;
    Real par = 0.0;
    Real pd = 0.0;

    Real ec = 0.0;
    Real gamma = 0.0;
    Real beta = 0.0;
    Real w = 1.0;
    Real kc = 0.0;
    Real tauc = 0.0;

    Real qIn = 0.0;
    Real qOutSum = 0.0;
    Real qDistalSum = 0.0;
    Real qCapChargingSum = 0.0;
    Real flowBalance = 0.0;

    std::map<Attribute, Real> qOut;
    std::map<Attribute, Real> qDistal;
    std::map<Attribute, Real> pc;
    std::map<Attribute, Real> pOut;
  };

  static Model::Input makeInput(const Config &cfg);
  static MeshType makeMesh(const Rodin::Context::MPI &context,
                           const Config &cfg);

  static void updateRCR(const Model &model, RCR &bc, Real Q, Real dt);
  static void updateRCRNonNew(const Config &cfg, const Model &model, RCR &bc,
                              Real Q, Real dt);
  static Real periodic_activation(const Activation &cfg, Real t);
  static Real atrial_pressure(const AtrialPressure &cfg, Real t);

  void setupModel();
  void setupMeshAndSpaces();
  void setupDiagnostics();
  void printInitialState() const;
  bool isRoot() const;

  bool advance0D();
  void solveStatic();
  bool solve3D();
  void computeFluxes();
  void updateHistory();
  void writeOutputs();
  void writeCSVHeader();
  void writeCSVRow();
  StepData collectStepData() const;
  void printStepTiming(int step) const;

  Config m_cfg;
  Model::Input m_input;
  Model m_model;

  MeshType m_mesh;
  Rodin::IO::XDMF m_xdmf;

  VelocityFESType m_uh;
  PressureFESType m_ph;
  VMSFESType m_uph;
  TauFESType m_tauh;

  VelocityTrialFunctionType m_u;
  PressureTrialFunctionType m_p;
  PressureTrialFunctionType m_mu;

  VelocityTestFunctionType m_v;
  PressureTestFunctionType m_q;
  PressureTestFunctionType m_r;

  VelocityGridFunctionType m_uOld;
  PressureGridFunctionType m_pOld;
  PressureGridFunctionType m_one;
  PressureTestFunctionType m_qFlux;

  FluxLinearFormType m_flux;

  VMSTrialFunctionType m_sub;
  VMSGridFunctionType m_subOld;
  VMSTrialFunctionType m_up;
  VMSTestFunctionType m_vp;

  TauGridFunctionType m_tauOld;
  TauTrialFunctionType m_tau;
  TauTestFunctionType m_t;

  FlowProblemType m_flow;
  VMSProblemType m_l2ConvU;
  VMSProblemType m_subProjection;
  TauProblemType m_tauProjection;
  Rodin::Solver::KSP m_flowKSP;
  Rodin::Solver::SNES m_flowSolver;
  ScalarProjectionProblemType m_viscosityProjection;
  Rodin::Solver::KSP m_viscosityProjectionKSP;
  Rodin::Solver::KSP m_l2ConvUSolver;
  Rodin::Solver::KSP m_subProjectionSolver;
  Rodin::Solver::KSP m_tauProjectionSolver;

  std::map<Attribute, RCR> m_wk;

  mutable StepData m_stepData;
  std::ofstream m_csv;
  StepTiming m_stepTiming;
  bool m_flowFieldSplitsSet = false;
  bool m_initialized = false;
};
} // namespace Rodin::Examples::Heart

#endif
