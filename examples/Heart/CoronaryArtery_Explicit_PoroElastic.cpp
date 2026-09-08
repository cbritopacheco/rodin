/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * Explicit (staggered) sequential PETSc coronary ALE FSI, P1/P1 fluid and P1
 * solid, coupled in a CLOSED LOOP to the 0D poroelastic left ventricle
 * (Rodin/Heart/PoroelasticSphere).
 *
 * Closed loop (all exchanges lagged one step, staggered):
 *   0D -> 3D : inlet pressure p_in = p_ar (proximal Windkessel node);
 *              tissue pressure p_f (interstitial) as reference and drainage
 *              pressure of every outlet compartment.
 *   3D -> 0D : the inlet flux Q_in is drawn from the proximal Windkessel
 *              node (Cp dp_ar/dt + q_p + Q_in = 0); the compartment drainage
 *              sum_i q_v,i enters the porosity balance
 *              V_w0 dPhi/dt = Q_lumped + sum_i q_v,i - gamma_ven (p_f - p_sv).
 *   Outlet i : 3D --[R_a Phi_a, implicit]--> (p_c, C) --[R_v Phi_v]--> p_f,
 *              p_tm = p_c - p_f, Starling throat shut for p_tm <= 0.
 *   The residual conductance gamma_ar of the 0D input represents the
 *   myocardium NOT perfused by the resolved tree.
 *
 * Fluid: Carreau-Yasuda, conservative BDF1 ALE, orthogonal-subscale VMS
 * (convective + grad-div + PSPG), all stabilization parameters evaluated
 * pointwise at the quadrature point from the local lagged viscosity and the
 * element diameter.  Solid: graded Yeoh, total Lagrangian, Newmark.  Coupling:
 * Robin-Robin loose coupling (Burman, Durst, Fernandez, Guzman & Ruz 2025).
 *
 * Wall traction and wall shear stress are evaluated pointwise from the L2
 * recovered velocity gradient G: D = (G + G^T)/2, gamma = sqrt(2 D:D),
 * mu = mu_CY(gamma), t = p n - 2 mu D n, tau_w = tangential part of 2 mu D n.
 * The same recovered gradient feeds the Robin datum, the load transfer and
 * the WSS output, so the three cannot drift apart.
 *
 * Attributes: fluid FSI 2, inlet 4, outlets 7,8,9,10,14,15; solid FSI 1,
 * inlet 150, outlets 151..156; 99 = clamped ring band; 101 = heart contact.
 */

#include "Rodin/Solid/Integrators/InternalVirtualWorkResidual.h"
#include "Rodin/Solid/Integrators/InternalVirtualWorkTangent.h"
#include "Rodin/Variational/BoundaryIntegral.h"
#include "Rodin/Variational/ForwardDecls.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <petscsys.h>
#include <petscvec.h>

#include <Rodin/Alert.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/PETSc.h>
#include <Rodin/Solid.h>
#include <Rodin/Solver.h>
#include <Rodin/Variational.h>

#include "Rodin/Heart/PoroelasticSphere.h"

#include "CoronaryArtery/GradedYeoh.h"
#include "CoronaryArtery/VMSConvectionIntegrator.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Math;
using namespace Rodin::Solver;
using namespace Rodin::Variational;
using namespace Rodin::Heart;

namespace
{
  using Model = Rodin::Heart::PoroelasticSphereT<>;
  using MeshType = Geometry::Mesh<Context::Local>;

  struct BoundaryFluid
  {
      static constexpr Attribute FSI = 2;
      static constexpr Attribute Inlet = 4;
      static constexpr std::array<Attribute, 6> Outlets{{7, 8, 9, 10, 14, 15}};
      static constexpr Attribute FSIRing = 99;
  };

  struct BoundarySolid
  {
      static constexpr Attribute FSI = 1;
      static constexpr Attribute Inlet = 150;
      static constexpr std::array<Attribute, 6> Outlets{{151, 152, 153, 154, 155, 156}};
      static constexpr Attribute FSIRing = 99;
      static constexpr std::array<Attribute, 1> Contact{{101}};
      static constexpr std::array<Attribute, 2> Outer{{2, 101}};
  };

  /// One R-mu compartment per outlet, embedded in the poroelastic tissue.
  struct Outlet
  {
      Real ptm = 0.0;   ///< Transmural pressure p_tm = p_c - p_f (the state).
      Real pc = 0.0;    ///< Compartment pressure p_c = p_tm + p_f.
      Real pout = 0.0;  ///< Neumann datum of the 3D outlet (= p_c; R_a Phi_a Q is implicit).
      Real qv = 0.0;    ///< Drainage into the interstitium.
      Real pf = 0.0;    ///< Tissue (interstitial) pressure seen by the compartment.
      Real Ra = 0.0;    ///< Arteriolar resistance at mu_N.
      Real Rv = 0.0;    ///< Venular resistance at mu_N.
      Real C = 0.0;     ///< Compartment compliance.
      Real area = 0.0;  ///< Outlet area (implicit term R_a Phi_a A (u.n)(v.n)).
      Real q0 = 0.0;    ///< Calibrated branch flow (reference of Phi(q)).
      Real gammaA = 0.0; ///< Resting wall shear rate, arteriolar limb.
      Real gammaV = 0.0; ///< Resting wall shear rate, venular limb.
      Real phiA = 1.0;  ///< mu_ap/mu_N, arteriolar limb.
      Real phiV = 1.0;  ///< mu_ap/mu_N, venular limb.
      Real vol = 0.0;   ///< Stored volume C max(p_tm, 0).
  };

  struct CarreauYasuda
  {
      Real mu0 = 0.301;
      Real muInf = 0.0055;
      Real lambda = 16.15;
      Real n = 0.21;
      Real yasuda = 0.77;
      Real gammaRegularization = 1.0e-3;

      Real operator()(Real shearRate) const
      {
        return muInf +
          (mu0 - muInf) *
          std::pow(1.0 + std::pow(lambda * shearRate, yasuda), (n - 1.0) / yasuda);
      }
  };

  enum class RheologyModel
  {
      CarreauYasuda,
      Quemada
  };

  struct Quemada
  {
      Real plasmaViscosity = 0.0017963;
      Real hematocrit = 0.45;
      Real k0 = 3.7;      // derived from phi if <= 0
      Real kInf = 1.66;   // derived from phi if <= 0
      Real gammaC = 2.29; // derived from phi if <= 0
  };

  struct OutletFlowLaw
  {
      Real tableTauMin = 1.0e-6;
      Real tableTauMax = 1.0e4;
      int tableNodes = 241;
      int integralSteps = 2000;
      Real shearStepTolerance = 1.0e-12;
      int shearMaxIterations = 200;
      Real outletStepTolerance = 1.0e-9;
      int outletMaxIterations = 50;
      Real zeroFlowTolerance = 1.0e-16;
  };

  /// Universal WRMS apparent viscosity mu_ap(gammadot_nom), log-log, clamped.
  struct WRMSTable
  {
      std::vector<Real> logGamma;
      std::vector<Real> logMu;

      Real operator()(Real gamma) const
      {
        if (logGamma.size() < 2)
          return std::exp(logMu.empty() ? 0.0 : logMu.front());

        const Real lg = std::log(std::max(gamma, 1e-300));

        if (lg <= logGamma.front())
          return std::exp(logMu.front());
        if (lg >= logGamma.back())
          return std::exp(logMu.back());

        const auto it = std::upper_bound(logGamma.begin(), logGamma.end(), lg);
        const std::size_t i = static_cast<std::size_t>(it - logGamma.begin()) - 1;
        const Real w = (lg - logGamma[i]) / (logGamma[i + 1] - logGamma[i]);

        return std::exp(logMu[i] + w * (logMu[i + 1] - logMu[i]));
      }
  };

  struct Config
  {
      std::string resultsDir = "results_poroelastic";
      std::string xdmfBasename = "results_poroelastic/CoronaryArtery_Explicit_PoroElastic";
      std::string csvPath = "results_poroelastic/CoronaryArtery_Explicit_PoroElastic.csv";
      Real meshScale = 1.0e-3;

      Real cyclePeriod = 0.85;
      Real dt = 1.0e-3;
      size_t nsteps = 4 * static_cast<size_t>(0.85 / 1.0e-3);

      /// 0D cycles integrated alone before the coupling; the resolved tree is
      /// replaced by the lumped conductance gammaLca = Q_lca/(p_ar - p_f) so
      /// the hand-off starts from a periodic 0D state.  Time is reset to 0.
      size_t zeroDWarmupCycles = 2;
      /// Tissue pressure used for the lumped-tree estimate before any 0D
      /// step has produced p_f (only relevant with zeroDWarmupCycles == 0).
      Real restingTissuePressure = 4.6e3;

      // Outlet calibration: areas measured on the mesh, Murray split
      // Q_i ~ r_i^3, resting budget dP = p_ar(0) - p_f(0):
      //   R_v,i = f_v dP/Q_i,  R_a,i = (1 - f_v) dP/Q_i,  C_i = C_tot w_i.
      bool autoCalibrateOutlets = true;
      Real lcaTargetFlow = 1.5e-6; // m^3/s
      Real newtonianCalibrationViscosity = 0.0035;
      Real coronaryComplianceTotal = 4e-10; // m^3/Pa
      Real venularPressureFraction = 0.13;

      // Morphometric operating point (r, v) of each limb: gamma_0 = 4v/r.
      Real arteriolarRadius = 25.0e-6;
      Real arteriolarVelocity = 5.0e-3;
      Real venularRadius = 30.0e-6;
      Real venularVelocity = 3.0e-3;
      Real referenceTransitTime = 1.5;

      RheologyModel rheologyModel = RheologyModel::CarreauYasuda;
      Quemada quemada;

      Real fluidDensity = 1060.0;
      Real inletBackflowStabilization = 1.0;
      Real outletBackflowStabilization = 1.0;
      CarreauYasuda viscosity;
      OutletFlowLaw outletFlowLaw;

      Real inletImpedance = 1e3;
      Real inletTangentialDamping = 1e3;
      Real outletResistanceScale = 1.0;

      size_t prestressSteps = 50;
      Real prestressFraction = 0.975;
      size_t prestressRampSteps = 10;

      Real vmsScale = 1.0;
      Real gradDivScale = 1.0;
      Real pgpScale = 1.0;

      // ALE harmonic lift weight w = (aleRefSize/h_K)^aleStiffPower.
      Real aleStiffPower = 0.75;
      Real aleRefSize = 5.0e-4;

      bool meshConsistentInterfaceVelocity = true;
      bool subtractHeartHandoffOffset = true;

      Real solidDensity = 1060.0;
      Real solidViscosity = 4.e3;
      Real newmarkBeta = 0.25;
      Real newmarkGamma = 0.5;

      size_t couplingIterations = 1;
      Real couplingTolerance = 1.0e-6;

      Real robinAlpha = 0.0;
      Real robinGamma = 1.0;

      Real heartDisplacementPenalty = 1.e7;
      Real heartDisplacementScale = 1.0;

      // Viscoelastic tethering of the outlet annuli (Pa/m, Pa s/m).
      Real aViscCondition = 1.0e5;
      Real bViscCondition = 1.0e3;

      // Transmural Yeoh multipliers (xi < 1/3 intima, middle media, > 2/3
      // adventitia; smoothstep half-width gradeTransitionWidth).
      Real gradeIntima = 0.5;
      Real gradeMedia = 1.0;
      Real gradeAdventitia = 1.5;
      Real gradeTransitionWidth = 0.08;
  };

  static Real smoothstep(Real s)
  {
    s = s < 0.0 ? 0.0 : (s > 1.0 ? 1.0 : s);
    return s * s * (3.0 - 2.0 * s);
  }

  static Real periodic_activation(Real t)
  {
    const Real T = 0.85;
    const Real tau = t - T * std::floor(t / T);

    const Real tRampStart = 0.15, tRampEnd = 0.21, tPlateauEnd = 0.36;
    const Real tRelaxEnd = 0.45, tNegativeEnd = 0.6;
    const Real positiveValue = 35.0, negativeValue = -20.0;

    if (tau < tRampStart)
      return 0.0;
    if (tau < tRampEnd)
      return positiveValue * smoothstep((tau - tRampStart) / (tRampEnd - tRampStart));
    if (tau < tPlateauEnd)
      return positiveValue;
    if (tau < tRelaxEnd)
      return positiveValue +
        (negativeValue - positiveValue) *
        smoothstep((tau - tPlateauEnd) / (tRelaxEnd - tPlateauEnd));
    if (tau < tNegativeEnd)
      return negativeValue;
    return negativeValue * (1.0 - smoothstep((tau - tNegativeEnd) / (T - tNegativeEnd)));
  }

  static Real load_dependent_relaxation_m0(Real ec)
  {
    const Real lowEc = 0.0, highEc = 2.0, lowValue = 1.6, highValue = 1.0;
    if (ec <= lowEc)
      return lowValue;
    if (ec >= highEc)
      return highValue;
    const Real s = (ec - lowEc) / (highEc - lowEc);
    return (1.0 - s) * lowValue + s * highValue;
  }

  static Real load_dependent_relaxation_dm0(Real ec)
  {
    const Real lowEc = 0.0, highEc = 2.0, lowValue = 1.6, highValue = 1.0;
    if (ec <= lowEc || ec >= highEc)
      return 0.0;
    return (highValue - lowValue) / (highEc - lowEc);
  }

  static Real atrial_pressure(Real t)
  {
    const Real T = 0.85;
    const Real tau = t - T * std::floor(t / T);
    const Real minValue = 500.0, maxValue = 1000.0, secondThreshold = 1250.0;
    const Real t1 = 0.02, t2 = 0.15, t3 = 0.17, t4 = 0.56, t5 = 0.62, t6 = 0.85;

    auto ramp = [](Real a, Real b, Real s) { return a + (b - a) * smoothstep(s); };

    if (tau < t1)
      return ramp(minValue, maxValue, tau / t1);
    if (tau < t2)
      return maxValue;
    if (tau < t3)
      return ramp(maxValue, minValue, (tau - t2) / (t3 - t2));
    if (tau < t4)
      return ramp(minValue, secondThreshold, (tau - t3) / (t4 - t3));
    if (tau < t5)
      return secondThreshold;
    if (tau < t6)
      return ramp(secondThreshold, minValue, (tau - t5) / (t6 - t5));
    return minValue;
  }

  /// 0D poroelastic LV input.  The two external sources read the lagged 3D
  /// fluxes through the references (staggered coupling: data at t_{n+1},
  /// absent from the 0D Jacobian).
  static Model::Input makeModelInput(const Real& qArterialLag, const Real& qTissueLag)
  {
    Model::Input in;

    in.R0 = 2.4e-2;
    in.d0 = 1.45e-2;
    in.phi0 = 0.1;

    in.Es = 3.0e6;
    in.mu = 70.0;
    in.eta = 70.0;
    in.alpha = 1.5;
    in.alphaR = 0.12;
    in.k0 = 1.0e5;
    in.sigma0 = 1.25e5;

    // Perfusion: storage modulus, residual lumped arterial conductance (the
    // territory not resolved by the 3D tree) and venous conductance.
    in.KPhi = 2.0e5;
    in.gammaAr = 4.7e-10;
    in.gammaVen = 7.0e-10;

    in.Rp = 5.0e7;
    in.Cp = 6e-9;
    in.Rd = 1.0e8;
    in.Cd = 1.0e-9;

    in.mu_0 = 5.35;
    in.mu_Inf = 0.0033;
    in.lambda = 14.445;
    in.n = 0.8;
    in.m = 0.003;
    in.yasuda = 0.62;
    in.mu_plasma = 0.0032704;
    in.k_0 = 3.5678;
    in.gamma_c = 10.2754;
    in.k_Inf = 1.5352;
    in.proximalRadius = 0.0125;
    in.proximalLength = 0.35;
    in.distalRadius = 0.002;
    in.distalLength = 0.55;
    in.windkesselRheology = Rodin::Heart::CCMLC2014::Model::WindkesselRheology::Cross;

    in.Kat = 6.0e-7;
    in.Kp = 5.0e-11;
    in.Kar = 1.0e-7;
    in.cavityCapacity = 5.0e-12;

    in.absRegularization = 1e-14;
    in.wallQuadraturePoints = 8;

    in.initFibDef = 0.0;
    in.initActiveStiffness = 0.0;
    in.initActiveStress = 0.0;

    in.pSv = [](Real) { return 1.0e3; };
    in.pAt = atrial_pressure;
    in.u = periodic_activation;
    in.m0 = load_dependent_relaxation_m0;
    in.dm0 = load_dependent_relaxation_dm0;

    in.qArterialExternal = [&qArterialLag](Real) { return qArterialLag; };
    in.qPerfusionExternal = [&qTissueLag](Real) { return qTissueLag; };

    using PassiveEnergy = std::decay_t<decltype(in.passiveEnergy)>;
    typename PassiveEnergy::Parameters hp;
    hp.mu1 = 0.0;
    hp.mu2 = 0.0;
    hp.C0 = 1.9e3;
    hp.C1 = 1.1e-1;
    hp.C2 = 1.9e3;
    hp.C3 = 1.1e-1;
    in.passiveEnergy = PassiveEnergy(hp);

    return in;
  }

  static void initializeModel(Model& model, const Model::Input& in)
  {
    model.setMaxIterations(200)
      .setAbsoluteTolerance(1.0e-8)
      .setRelativeTolerance(1.0e-8)
      .setStepTolerance(1.0e-10)
      .setDampingFactor(1.0);

    Model::State s0;
    s0.t = 0.0;
    s0.y = 0.0;
    s0.phi = in.phi0;
    s0.pv = in.pAt(0.0) - 100.0;
    s0.par = 11000.0;
    s0.pd = 10000.0;
    s0.ec = in.initFibDef;
    s0.gamma = std::sqrt(std::max<Real>(in.initActiveStiffness, 0.0));
    s0.beta = (s0.gamma > 0.0) ? (in.initActiveStress / s0.gamma) : 0.0;
    s0.kc = s0.gamma * s0.gamma;
    s0.tauc = s0.gamma * s0.beta;
    s0.w = in.m0(s0.ec);

    model.initialize(s0);
  }

  /// Epicardial radial displacement r_out - R_out of the uniform-J shell:
  /// r_out^3 = (R_in + y)^3 + J (R_out^3 - R_in^3), J = 1 - phi0 + Phi.
  static Real epicardialDisplacement(const Model::Input& in, const Model::State& s)
  {
    const Real Rin = in.R0;
    const Real Rout = in.R0 + in.d0;
    const Real rin = Rin + s.y;
    const Real J = 1.0 - in.phi0 + s.phi;
    const Real rout3 = rin * rin * rin + J * (Rout * Rout * Rout - Rin * Rin * Rin);
    return std::cbrt(std::max<Real>(rout3, 0.0)) - Rout;
  }

  static MeshType makeMesh(const Config& cfg, const std::string& meshPath)
  {
    MeshType mesh;
    mesh.load(meshPath, IO::FileFormat::MEDIT);

    if (mesh.getSpaceDimension() != 3)
      throw std::runtime_error("Expected a 3D coronary FSI mesh.");

    mesh.scale(cfg.meshScale);

    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);
    mesh.getConnectivity().compute(D, 0);
    mesh.getConnectivity().compute(D, D - 1);
    mesh.getConnectivity().compute(D - 1, D);
    mesh.getConnectivity().compute(D - 1, 0);
    mesh.getConnectivity().compute(D - 1, 1);
    mesh.getConnectivity().compute(1, 0);

    return mesh;
  }

  template <class MeshT>
  static void saveReferenceVertices(
    const MeshT& mesh, std::vector<Math::SpatialPoint>& vertices)
  {
    vertices.resize(mesh.getVertexCount());
    for (auto it = mesh.getVertex(); it; ++it)
      vertices[it->getIndex()] = mesh.getVertexCoordinates(it->getIndex());
  }

  template <class MeshT>
  static void restoreMeshToReference(
    MeshT& mesh, const std::vector<Math::SpatialPoint>& referenceVertices)
  {
    assert(mesh.getVertexCount() == referenceVertices.size());
    for (auto it = mesh.getVertex(); it; ++it)
    {
      const Index vertex = it->getIndex();
      mesh.setVertexCoordinates(vertex, referenceVertices[vertex]);
    }
    mesh.flush();
  }

  template <class MeshT, class FESType, class GridFunctionType>
  static void moveMeshWithVertexDisplacement(MeshT& mesh,
    const std::vector<Math::SpatialPoint>& referenceVertices,
    const FESType& displacementFES, const GridFunctionType& displacement)
  {
    assert(mesh.getVertexCount() == referenceVertices.size());
    const size_t dim = mesh.getSpaceDimension();

    for (auto it = mesh.getVertex(); it; ++it)
    {
      const Index vertex = it->getIndex();
      auto x = referenceVertices[vertex];

      for (Index c = 0; c < static_cast<Index>(dim); ++c)
        x(c) += displacement[displacementFES.getGlobalIndex({0, vertex}, c)];

      mesh.setVertexCoordinates(vertex, x);
    }
    mesh.flush();
  }

  /// Element diameter (longest edge) of every cell of the CURRENT mesh
  /// configuration; the stabilization length h_K of the VMS parameters.
  static void computeCellDiameters(const MeshType& mesh, std::vector<Real>& h)
  {
    h.assign(mesh.getCellCount(), 0.0);
    std::vector<Index> verts;
    for (auto it = mesh.getCell(); it; ++it)
    {
      verts.clear();
      for (const auto& v : it->getVertices())
        verts.push_back(v);

      Real hmax = 0.0;
      for (std::size_t a = 0; a < verts.size(); ++a)
        for (std::size_t b = a + 1; b < verts.size(); ++b)
        {
          const Real d = Real(
            (mesh.getVertexCoordinates(verts[a]) - mesh.getVertexCoordinates(verts[b])).norm());
          hmax = std::max(hmax, d);
        }

      h[it->getIndex()] = hmax;
    }
  }

  struct InterfaceSegment
  {
      Index fluid;
      Index solid;
  };

  struct InterfaceMap
  {
      std::unordered_map<Index, Index> fluidToSolid;
      std::unordered_map<Index, Index> solidToFluid;
      std::vector<InterfaceSegment> segments;
  };

  static Math::SpatialPoint centroid(const MeshType& mesh, const Polytope& polytope)
  {
    Math::SpatialPoint c(mesh.getSpaceDimension());
    c.setZero();

    const auto& vertices = polytope.getVertices();
    for (const auto& v : vertices)
      c += mesh.getVertexCoordinates(v);

    c /= static_cast<Real>(vertices.size());
    return c;
  }

  /// Relabel to `ring` every FSI face touching a vertex of the given caps.
  static std::size_t tagFSIRingBand(MeshType& mesh, Attribute fsi, Attribute ring,
    const std::vector<Attribute>& caps)
  {
    const std::size_t faceDim = mesh.getDimension() - 1;

    const auto isCap = [&](const Optional<Attribute>& a) {
      if (!a)
        return false;
      for (const Attribute c : caps)
        if (*a == c)
          return true;
      return false;
    };

    std::unordered_set<Index> capVertices;
    for (auto it = mesh.getBoundary(); it; ++it)
    {
      if (!isCap(it->getAttribute()))
        continue;
      for (const auto& v : it->getVertices())
        capVertices.insert(v);
    }

    std::vector<Index> toRelabel;
    for (auto it = mesh.getBoundary(); it; ++it)
    {
      if (it->getAttribute() != fsi)
        continue;
      for (const auto& v : it->getVertices())
        if (capVertices.count(v))
        {
          toRelabel.push_back(it->getIndex());
          break;
        }
    }

    for (const Index f : toRelabel)
      mesh.setAttribute({faceDim, f}, ring);

    return toRelabel.size();
  }

  /// Fluid/solid FSI face pairing by centroid proximity (conforming interface).
  static InterfaceMap buildInterfaceMap(
    const MeshType& fluidReferenceMesh, const MeshType& solidReferenceMesh)
  {
    InterfaceMap map;

    std::vector<std::pair<Math::SpatialPoint, Index>> solidFaces;
    Math::SpatialPoint lo, hi;
    bool haveBox = false;
    for (auto it = solidReferenceMesh.getBoundary(); it; ++it)
    {
      if (it->getAttribute() != BoundarySolid::FSI)
        continue;
      Math::SpatialPoint c = centroid(solidReferenceMesh, *it);
      solidFaces.emplace_back(c, it->getIndex());
      if (!haveBox)
      {
        lo = c;
        hi = c;
        haveBox = true;
      }
      else
      {
        for (Index i = 0; i < static_cast<Index>(c.size()); ++i)
        {
          lo(i) = std::min(lo(i), c(i));
          hi(i) = std::max(hi(i), c(i));
        }
      }
    }

    const Real diag = haveBox ? Real((hi - lo).norm()) : Real(0);
    const Real tol = std::max(Real(1e-12), Real(1e-6) * diag);

    auto cellOf = [&](const Math::SpatialPoint& c) {
      std::array<long long, 3> k{0, 0, 0};
      const Index n = std::min<Index>(3, static_cast<Index>(c.size()));
      for (Index i = 0; i < n; ++i)
        k[i] = static_cast<long long>(std::floor(c(i) / tol));
      return k;
    };
    std::map<std::array<long long, 3>, std::vector<std::size_t>> grid;
    for (std::size_t s = 0; s < solidFaces.size(); ++s)
      grid[cellOf(solidFaces[s].first)].push_back(s);

    std::size_t fluidCount = 0;
    std::size_t unmatched = 0;
    Real worst = 0.0;
    for (auto it = fluidReferenceMesh.getBoundary(); it; ++it)
    {
      if (it->getAttribute() != BoundaryFluid::FSI)
        continue;
      ++fluidCount;

      const Math::SpatialPoint c = centroid(fluidReferenceMesh, *it);
      const auto base = cellOf(c);

      Real best = std::numeric_limits<Real>::max();
      Index bestSolid = 0;
      bool found = false;
      for (long long dx = -1; dx <= 1; ++dx)
        for (long long dy = -1; dy <= 1; ++dy)
          for (long long dz = -1; dz <= 1; ++dz)
          {
            const std::array<long long, 3> key{base[0] + dx, base[1] + dy, base[2] + dz};
            const auto g = grid.find(key);
            if (g == grid.end())
              continue;
            for (const std::size_t s : g->second)
            {
              const Real d = Real((solidFaces[s].first - c).norm());
              if (d < best)
              {
                best = d;
                bestSolid = solidFaces[s].second;
                found = true;
              }
            }
          }

      if (!found || best > tol)
      {
        ++unmatched;
        if (found)
          worst = std::max(worst, best);
        continue;
      }

      const Index fluidFace = it->getIndex();
      map.fluidToSolid.emplace(fluidFace, bestSolid);
      map.solidToFluid.emplace(bestSolid, fluidFace);
      map.segments.push_back({fluidFace, bestSolid});
    }

    if (unmatched > 0)
    {
      std::ostringstream os;
      os << "buildInterfaceMap: " << unmatched << " of " << fluidCount
         << " fluid FSI face(s) had no solid twin within tol=" << tol
         << " (nearest distance up to " << worst
         << "; solid FSI faces=" << solidFaces.size() << ").";
      throw std::runtime_error(os.str());
    }

    return map;
  }

  static Point forwardFluidPointToSolid(
    const Point& p, const MeshType& solidReferenceMesh, const InterfaceMap& map)
  {
    const auto found = map.fluidToSolid.find(p.getPolytope().getIndex());
    if (found == map.fluidToSolid.end())
      throw std::runtime_error("Fluid point is not on a mapped FSI face.");

    auto solidFace = solidReferenceMesh.getFace(found->second);
    const Math::SpatialPoint pc = p.getPhysicalCoordinates();
    Math::SpatialPoint rc;
    solidFace->getTransformation().inverse(rc, pc);
    return Point(*solidFace, rc, pc);
  }

  static Point forwardSolidPointToFluid(
    const Point& p, const MeshType& fluidMesh, const InterfaceMap& map)
  {
    const auto found = map.solidToFluid.find(p.getPolytope().getIndex());
    if (found == map.solidToFluid.end())
      throw std::runtime_error("Solid point is not on a mapped FSI face.");

    auto fluidFace = fluidMesh.getFace(found->second);
    const Math::SpatialPoint pc = p.getPhysicalCoordinates();
    Math::SpatialPoint rc;
    fluidFace->getTransformation().inverse(rc, pc);
    return Point(*fluidFace, rc, pc);
  }

  static bool isFiniteVec(const Math::SpatialVector<Real>& x)
  {
    for (Index i = 0; i < x.size(); ++i)
      if (!std::isfinite(x(i)))
        return false;
    return true;
  }

  static void setPETScDefault(const char* key, const char* value)
  {
    PetscBool set = PETSC_FALSE;
    PetscOptionsHasName(PETSC_NULLPTR, PETSC_NULLPTR, key, &set);
    if (!set)
      PetscOptionsSetValue(PETSC_NULLPTR, key, value);
  }

  static void readOptions(Config& cfg)
  {
    auto optReal = [](const char* key, Real& value, Real lo, Real hi) {
      PetscReal v = value;
      PetscBool set = PETSC_FALSE;
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, key, &v, &set);
      if (set)
        value = std::clamp<Real>(v, lo, hi);
    };
    auto optInt = [](const char* key, size_t& value, PetscInt lo) {
      PetscInt v = static_cast<PetscInt>(value);
      PetscBool set = PETSC_FALSE;
      PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, key, &v, &set);
      if (set)
        value = static_cast<size_t>(std::max<PetscInt>(lo, v));
    };
    auto optBool = [](const char* key, bool& value) {
      PetscBool v = value ? PETSC_TRUE : PETSC_FALSE;
      PetscBool set = PETSC_FALSE;
      PetscOptionsGetBool(PETSC_NULLPTR, PETSC_NULLPTR, key, &v, &set);
      if (set)
        value = (v == PETSC_TRUE);
    };

    const Real inf = std::numeric_limits<Real>::infinity();

    optReal("-coronary_dt", cfg.dt, 1.0e-9, inf);
    optInt("-coronary_nsteps", cfg.nsteps, 0);
    optInt("-coronary_warmup_cycles", cfg.zeroDWarmupCycles, 0);
    optReal("-coronary_resting_tissue_pressure", cfg.restingTissuePressure, 0.0, inf);
    optInt("-coronary_coupling_iterations", cfg.couplingIterations, 1);
    optReal("-coronary_coupling_tolerance", cfg.couplingTolerance, 0.0, inf);
    optInt("-coronary_prestress_steps", cfg.prestressSteps, 0);
    optReal("-coronary_prestress_fraction", cfg.prestressFraction, 0.0, 1.0);
    optInt("-coronary_prestress_ramp_steps", cfg.prestressRampSteps, 0);
    optReal("-coronary_vms_scale", cfg.vmsScale, 0.0, inf);
    optReal("-coronary_graddiv_scale", cfg.gradDivScale, 0.0, inf);
    optReal("-coronary_pgp_scale", cfg.pgpScale, 0.0, inf);
    optReal("-coronary_inlet_impedance", cfg.inletImpedance, 0.0, inf);
    optReal("-coronary_heart_disp_penalty", cfg.heartDisplacementPenalty, 0.0, inf);
    optReal("-coronary_heart_disp_scale", cfg.heartDisplacementScale, -inf, inf);
    optBool("-coronary_mesh_consistent_interface_vel", cfg.meshConsistentInterfaceVelocity);
    optBool("-coronary_subtract_heart_offset", cfg.subtractHeartHandoffOffset);

    // Newmark: gamma >= 1/2, beta = (gamma + 1/2)^2 / 4 unless overridden.
    {
      PetscReal g = cfg.newmarkGamma;
      PetscBool set = PETSC_FALSE;
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_newmark_gamma", &g, &set);
      if (set)
      {
        cfg.newmarkGamma = std::max<Real>(0.5, g);
        const Real s = cfg.newmarkGamma + 0.5;
        cfg.newmarkBeta = 0.25 * s * s;
      }
      optReal("-coronary_newmark_beta", cfg.newmarkBeta, 1.0e-6, inf);
    }

    optReal("-coronary_solid_viscosity", cfg.solidViscosity, 0.0, inf);
    optReal("-coronary_ale_stiff_power", cfg.aleStiffPower, 0.0, inf);
    optReal("-coronary_ale_ref_size", cfg.aleRefSize, 1.0e-30, inf);
    optReal("-coronary_outlet_resistance_scale", cfg.outletResistanceScale, 0.0, inf);
    optBool("-coronary_auto_calibrate", cfg.autoCalibrateOutlets);
    optReal("-coronary_lca_flow", cfg.lcaTargetFlow, 1e-12, inf);
    optReal("-coronary_robin_alpha", cfg.robinAlpha, -inf, inf);
    optReal("-coronary_robin_gamma", cfg.robinGamma, -inf, inf);
    optReal("-coronary_grade_intima", cfg.gradeIntima, 0.0, inf);
    optReal("-coronary_grade_media", cfg.gradeMedia, 0.0, inf);
    optReal("-coronary_grade_adventitia", cfg.gradeAdventitia, 0.0, inf);
  }

  // WRMS closure for a tube carrying a generalized-Newtonian fluid:
  //   Q = pi R^3 I(tau_w)/tau_w^3,  I = int_0^{tau_w} tau^2 gd dtau,
  //   mu_ap(tau_w) = tau_w^4/(4 I),  gd_nom = 4Q/(pi R^3) = tau_w/mu_ap,
  // independent of R and L, so it is tabulated once and shared by all outlets.
  static WRMSTable buildWRMSTable(const Config& cfg, const CarreauYasuda& visc)
  {
    const auto& law = cfg.outletFlowLaw;

    const Real mu0 = visc.mu0;
    const Real muInf = visc.muInf;
    const Real lambda = visc.lambda;
    const Real n = visc.n;
    const Real yasuda = visc.yasuda;
    const Real delta = mu0 - muInf;

    const bool quemada = (cfg.rheologyModel == RheologyModel::Quemada);
    const auto& qp = cfg.quemada;

    const Real phi = std::clamp<Real>(qp.hematocrit, 0.0, 0.75);
    const Real p2 = phi * phi;
    const Real p3 = p2 * phi;

    // Cokelet correlations when the Quemada coefficients are left non-positive.
    const Real qk0 = (qp.k0 > 0.0)
      ? qp.k0 : std::exp(3.874 - 10.41 * phi + 13.8 * p2 - 6.738 * p3);
    const Real qkInf = (qp.kInf > 0.0)
      ? qp.kInf : std::exp(1.3435 - 2.803 * phi + 2.711 * p2 - 0.6479 * p3);
    const Real qgc = (qp.gammaC > 0.0)
      ? qp.gammaC : std::exp(-6.1508 + 27.923 * phi - 25.6 * p2 + 3.697 * p3);

    auto muQ = [&](Real g) -> Real {
      const Real s = std::sqrt(std::max<Real>(g, 0.0) / qgc);
      const Real k = (qk0 + qkInf * s) / (1.0 + s);
      const Real b = std::max<Real>(1.0 - 0.5 * k * phi, 1e-3);
      return qp.plasmaViscosity / (b * b);
    };

    auto mu = [&](Real g) -> Real {
      if (quemada)
        return muQ(g);
      return muInf +
        delta * std::pow(1.0 + std::pow(lambda * g, yasuda), (n - 1.0) / yasuda);
    };

    auto dmu = [&](Real g) -> Real {
      if (quemada)
      {
        const Real h = 1e-6 * std::max<Real>(g, 1e-12);
        return (muQ(g + h) - muQ(std::max<Real>(g - h, 0.0))) / (2.0 * h);
      }
      const Real base = 1.0 + std::pow(lambda * g, yasuda);
      return delta * (n - 1.0) * std::pow(base, (n - 1.0 - yasuda) / yasuda) *
        std::pow(lambda, yasuda) * std::pow(g, yasuda - 1.0);
    };

    // gd(tau_w): tau = mu(gd) gd is strictly increasing; safeguarded Newton
    // on log gd.
    auto shearAt = [&](Real tauW) -> Real {
      Real lo = std::log(1e-14);
      Real hi = std::log(std::max<Real>(tauW / muInf, 1e-12)) + 1.0;
      Real s = 0.5 * (lo + hi);

      for (int it = 0; it < law.shearMaxIterations; ++it)
      {
        const Real g = std::exp(s);
        const Real f = mu(g) * g - tauW;

        if (f < 0.0)
          lo = s;
        else
          hi = s;

        const Real df = (mu(g) + g * dmu(g)) * g;
        Real sNext = (std::abs(df) > 0.0) ? s - f / df : 0.5 * (lo + hi);

        if (!(sNext > lo && sNext < hi))
          sNext = 0.5 * (lo + hi);

        const Real step = std::abs(sNext - s);
        s = sNext;

        if (step < law.shearStepTolerance)
          break;
      }

      return std::exp(s);
    };

    // I = int_0^{gd_w} gd^3 mu^2 (mu + gd mu') dgd (tau = mu gd), Simpson.
    auto rheologicalIntegral = [&](Real gammaW) -> Real {
      const int m = 2 * (law.integralSteps / 2);
      const Real h = gammaW / static_cast<Real>(m);

      auto f = [&](Real g) -> Real {
        if (g <= 0.0)
          return 0.0;
        const Real mg = mu(g);
        return g * g * g * mg * mg * (mg + g * dmu(g));
      };

      Real sum = f(0.0) + f(gammaW);
      for (int i = 1; i < m; ++i)
        sum += ((i % 2 == 1) ? 4.0 : 2.0) * f(static_cast<Real>(i) * h);

      return sum * h / 3.0;
    };

    WRMSTable table;

    const int nodes = std::max(2, law.tableNodes);
    const Real logTauMin = std::log(law.tableTauMin);
    const Real logTauMax = std::log(law.tableTauMax);

    table.logGamma.reserve(static_cast<std::size_t>(nodes));
    table.logMu.reserve(static_cast<std::size_t>(nodes));

    for (int i = 0; i < nodes; ++i)
    {
      const Real tauW = std::exp(logTauMin +
        (logTauMax - logTauMin) * static_cast<Real>(i) / static_cast<Real>(nodes - 1));

      const Real gammaW = shearAt(tauW);
      const Real I = rheologicalIntegral(gammaW);
      if (!(I > 0.0) || !std::isfinite(I))
        continue;

      const Real muAp = std::pow(tauW, 4.0) / (4.0 * I);
      if (!(muAp > 0.0) || !std::isfinite(muAp))
        continue;

      const Real gammaNom = tauW / muAp;
      if (!table.logGamma.empty() && std::log(gammaNom) <= table.logGamma.back())
        continue;

      table.logGamma.push_back(std::log(gammaNom));
      table.logMu.push_back(std::log(muAp));
    }

    if (table.logGamma.size() < 2)
    {
      table.logGamma = {std::log(1e-6), std::log(1e6)};
      table.logMu = {std::log(mu0), std::log(mu0)};
    }

    return table;
  }

  // Outlet compartment embedded in the poroelastic tissue at pressure p_f:
  //   3D --[R_a Phi_a, implicit in the 3D form]--> (p_c, C) --[R_v Phi_v]--> p_f
  //   C d(p_tm)/dt = Q - q_v,    p_tm = p_c - p_f,
  //   q_v = p_tm / (R_v Phi_v(|q_v|))  if p_tm > 0,   q_v = 0  otherwise
  // (Starling throat: external and downstream pressures both equal p_f).
  // Implicit Euler + scalar Newton on p_tm; dq_v/dp_tm >= 0 so R' > 0.
  // Exports p_out = p_c (the resistive part R_a Phi_a Q is assembled
  // implicitly on the 3D outlet) and the rheological factors Phi_a, Phi_v.
  static void updateOutlet0D(const Config& cfg, const WRMSTable& wrms, Real pf,
    Outlet& bc, Real Q, Real dt)
  {
    const auto& law = cfg.outletFlowLaw;
    const Real ptmOld = bc.ptm;

    const Real Rv = std::max<Real>(bc.Rv, 1e-300);
    const Real C = std::max<Real>(bc.C, 1e-300);
    const Real q0 = std::max<Real>(std::abs(bc.q0), 1e-300);
    const Real muN = std::max<Real>(cfg.newtonianCalibrationViscosity, 1e-300);

    auto viscosityFactorV = [&](Real q) -> Real {
      const Real aq = std::abs(q);
      const Real g = (aq < law.zeroFlowTolerance) ? law.zeroFlowTolerance
                                                  : bc.gammaV * aq / q0;
      return wrms(g) / muN;
    };

    Real ptm = ptmOld;
    Real qv = bc.qv;
    Real phiV = viscosityFactorV(qv);
    bool converged = false;

    for (int it = 0; it < law.outletMaxIterations; ++it)
    {
      phiV = viscosityFactorV(qv);
      const Real Gv = 1.0 / (Rv * phiV);
      const bool open = ptm > 0.0;

      qv = open ? ptm * Gv : 0.0;

      const Real R = C * (ptm - ptmOld) / dt - Q + qv;
      const Real J = C / dt + (open ? Gv : 0.0);

      const Real d = -R / J;
      ptm += d;

      if (std::abs(d) < law.outletStepTolerance * (1.0 + std::abs(ptm)))
      {
        converged = true;
        break;
      }
    }

    if (!converged || !std::isfinite(ptm))
    {
      std::cerr << "Warning: coronary outlet solve did not converge; keeping "
                << "the previous state.\n";
      ptm = ptmOld;
    }

    phiV = viscosityFactorV(qv);
    qv = (ptm > 0.0) ? ptm / (Rv * phiV) : 0.0;

    bc.ptm = ptm;
    bc.pf = pf;
    bc.pc = ptm + pf;
    bc.qv = qv;
    bc.vol = C * std::max<Real>(ptm, 0.0);
    bc.phiV = phiV;

    const Real aQ = std::abs(Q);
    const Real gammaA = (aQ < law.zeroFlowTolerance) ? law.zeroFlowTolerance
                                                     : bc.gammaA * aQ / q0;
    bc.phiA = wrms(gammaA) / muN;

    bc.pout = bc.pc;
  }

} // namespace

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  setPETScDefault("-ksp_type", "preonly");
  setPETScDefault("-pc_type", "lu");
  setPETScDefault("-pc_factor_mat_solver_type", "mumps");

  // Mass-matrix projections (VMS subscales, gradient recovery) are SPD.
  setPETScDefault("-coronary_mass_ksp_type", "cg");
  setPETScDefault("-coronary_mass_pc_type", "jacobi");

  const bool isRoot = true;

  try
  {
    Config cfg;
    readOptions(cfg);

    if (!cfg.resultsDir.empty())
    {
      std::error_code ec;
      std::filesystem::create_directories(cfg.resultsDir, ec);
      if (ec && isRoot)
        std::cerr << "Warning: could not create results directory '" << cfg.resultsDir
                  << "': " << ec.message() << '\n';
    }

    // ------------------------------------------------------------------
    // 0D poroelastic LV.  The lagged 3D fluxes enter through these two
    // scalars (flow drawn from p_ar; flow delivered into the interstitium).
    // ------------------------------------------------------------------
    Real qArterialLag = 0.0;
    Real qTissueLag = 0.0;
    Model::Input modelInput = makeModelInput(qArterialLag, qTissueLag);
    Model model(modelInput);
    initializeModel(model, modelInput);

    const Real dt = cfg.dt;
    const size_t stepsPerCycle =
      static_cast<size_t>(std::llround(cfg.cyclePeriod / dt));

    // 0D warm-up with the resolved tree replaced by a lumped conductance
    // gammaLca = Q_lca / (p_ar(0) - p_f,rest); the closed loop then hands
    // off from a periodic state.  Times are shifted back by whole cycles so
    // the coupled run starts at t = 0 with the same activation phase.
    const Real gammaLca = cfg.lcaTargetFlow /
      std::max<Real>(model.getState().par - cfg.restingTissuePressure, 1.0);
    if (cfg.zeroDWarmupCycles > 0)
    {
      const size_t warmSteps = cfg.zeroDWarmupCycles * stepsPerCycle;
      for (size_t k = 0; k < warmSteps; ++k)
      {
        const auto rep = model.step(dt);
        if (!rep.converged)
          throw std::runtime_error("0D warm-up failed at step " + std::to_string(k));
        const auto& s = model.getState();
        qArterialLag = gammaLca * (s.par - s.pf);
        qTissueLag = qArterialLag;
      }

      auto st = model.getState();
      auto hist = model.getHistory();
      const Real tShift = st.t;
      st.t -= tShift;
      hist.n.t -= tShift;
      hist.nm1.t -= tShift;
      if (hist.nm2)
        hist.nm2->t -= tShift;
      model.restore(st, hist, model.getUnknowns(), model.getReport());

      if (isRoot)
        Alert::Info() << "[0D] warm-up: " << cfg.zeroDWarmupCycles << " cycle(s), "
                      << "p_ar=" << st.par << " Pa  p_f=" << st.pf << " Pa  phi="
                      << st.phi << "  y=" << st.y << " m  (gammaLca=" << gammaLca
                      << " m^3/(Pa s))" << Alert::Raise;
    }

    const Real par0 = model.getState().par;
    const Real pf0 = (cfg.zeroDWarmupCycles > 0) ? model.getState().pf
                                                 : cfg.restingTissuePressure;

    const std::string fluidMesh = "../resources/examples/Heart/coronaria_estenosis_80.mesh";
    MeshType meshFluid = makeMesh(cfg, fluidMesh);
    const size_t dimFluid = meshFluid.getSpaceDimension();

    const std::string solidMesh = "../resources/examples/Heart/coronaria_80_prismatica.mesh";
    MeshType meshSolid = makeMesh(cfg, solidMesh);
    const size_t dimSolid = meshSolid.getSpaceDimension();

    const size_t dim = dimFluid;

    // Ring band only around the inlet; the outlet ends are held by the
    // viscoelastic tethering springs.
    const std::size_t fluidRingFaces = tagFSIRingBand(
      meshFluid, BoundaryFluid::FSI, BoundaryFluid::FSIRing, {BoundaryFluid::Inlet});
    const std::size_t solidRingFaces = tagFSIRingBand(
      meshSolid, BoundarySolid::FSI, BoundarySolid::FSIRing, {BoundarySolid::Inlet});
    if (isRoot)
      std::cout << "FSI ring band: fluid=" << fluidRingFaces
                << " face(s), solid=" << solidRingFaces << " face(s)\n";

    std::vector<Math::SpatialPoint> referenceVertices;
    saveReferenceVertices(meshFluid, referenceVertices);

    const InterfaceMap interfaceMap = buildInterfaceMap(meshFluid, meshSolid);

    using VelocityFES = H1<1, Math::SpatialVector<Real>, MeshType>;
    using PressureFES = H1<1, Real, MeshType>;
    using DisplacementFluidFES = H1<1, Math::SpatialVector<Real>, MeshType>;
    using DisplacementFES = H1<1, Math::SpatialVector<Real>, MeshType>;
    using LaplacianFES = H1<1, Real, MeshType>;

    VelocityFES uh(std::integral_constant<size_t, 1>{}, meshFluid, dimFluid);
    PressureFES ph(std::integral_constant<size_t, 1>{}, meshFluid);
    DisplacementFluidFES dfh(std::integral_constant<size_t, 1>{}, meshFluid, dimFluid);
    DisplacementFES dh(std::integral_constant<size_t, 1>{}, meshSolid, dimSolid);
    LaplacianFES lh(std::integral_constant<size_t, 1>{}, meshSolid);

    PETSc::Variational::TrialFunction u(uh);
    PETSc::Variational::TrialFunction p(ph);
    PETSc::Variational::TestFunction v(uh);
    PETSc::Variational::TestFunction q(ph);

    PETSc::Variational::TrialFunction d(dh);
    PETSc::Variational::TestFunction w(dh);

    u.setName("u");
    p.setName("p");
    d.setName("disp");

    PETSc::Variational::GridFunction uOld(uh);
    PETSc::Variational::GridFunction pOld(ph);
    PETSc::Variational::GridFunction one(ph);
    PETSc::Variational::GridFunction shearWall(uh);
    PETSc::Variational::GridFunction gradRec0(uh);
    PETSc::Variational::GridFunction gradRec1(uh);
    PETSc::Variational::GridFunction gradRec2(uh);

    PETSc::Variational::GridFunction aleDisp(uh);
    PETSc::Variational::GridFunction aleDispOld(uh);
    PETSc::Variational::GridFunction meshVelocity(uh);

    PETSc::Variational::GridFunction dState(dh);
    PETSc::Variational::GridFunction dOld(dh);
    PETSc::Variational::GridFunction dIter(dh);
    PETSc::Variational::GridFunction etaState(dh);
    PETSc::Variational::GridFunction dPred(dh);
    PETSc::Variational::GridFunction vPred(dh);
    PETSc::Variational::GridFunction solidVelocity(dh);
    PETSc::Variational::GridFunction solidAcceleration(dh);
    PETSc::Variational::GridFunction solidVelocityOld(dh);
    PETSc::Variational::GridFunction solidAccelerationOld(dh);

    PETSc::Variational::GridFunction fluidTraction(dfh);
    PETSc::Variational::GridFunction tractionTransfer(dfh);
    PETSc::Variational::GridFunction uWall(dfh);

    PETSc::Variational::TrialFunction l(lh);
    PETSc::Variational::TestFunction t(lh);

    // Transmural coordinate xi (0 lumen, 1 adventitia), harmonic.
    PETSc::Variational::TrialFunction xiT(lh);
    PETSc::Variational::TestFunction xiTest(lh);
    PETSc::Variational::GridFunction xi(lh);
    xi = 0.5;

    auto zero = VectorFunction(dim, [&](const Point&) {
      Math::SpatialVector<Real> value(dim);
      value.setZero();
      return value;
    });

    uOld = zero;
    shearWall = zero;
    gradRec0 = zero;
    gradRec1 = zero;
    gradRec2 = zero;
    pOld = 0.0;
    one = 1.0;
    dState = zero;
    dOld = zero;
    dIter = zero;
    etaState = zero;
    dPred = zero;
    vPred = zero;
    solidVelocity = zero;
    solidAcceleration = zero;
    solidVelocityOld = zero;
    solidAccelerationOld = zero;
    aleDisp = zero;
    aleDispOld = zero;
    meshVelocity = zero;
    fluidTraction = zero;
    tractionTransfer = zero;
    uWall = zero;

    uOld.setName("FluidVelocity");
    shearWall.setName("shearStress");
    pOld.setName("FluidPressure");
    dState.setName("Displacement");
    solidVelocity.setName("SolidVelocity");
    meshVelocity.setName("ALEMeshVelocity");
    fluidTraction.setName("FluidTraction");
    aleDisp.setName("ALEDisp");

    IO::XDMF xdmf_fluid(cfg.xdmfBasename + "fluid");
    xdmf_fluid.setMesh(meshFluid);
    xdmf_fluid.add("FluidVelocity", uOld);
    xdmf_fluid.add("FluidPressure", pOld);
    xdmf_fluid.add("ALEMeshVelocity", meshVelocity);
    xdmf_fluid.add("FluidTraction", fluidTraction);
    xdmf_fluid.add("ALEDisp", aleDisp);
    xdmf_fluid.add("shearStress", shearWall);
    xdmf_fluid.write(0.0).flush();

    IO::XDMF xdmf_solid(cfg.xdmfBasename + "solid");
    xdmf_solid.setMesh(meshSolid);
    xdmf_solid.add("Displacement", dState);
    xdmf_solid.add("SolidVelocity", solidVelocity);
    xdmf_solid.write(0.0).flush();

    IO::XDMF xdmf_laplacian(cfg.xdmfBasename + "laplacian");
    xdmf_laplacian.setMesh(meshSolid);
    xdmf_laplacian.add("laplacian", l.getSolution());

    std::ofstream csv;
    if (isRoot)
    {
      csv.open(cfg.csvPath);
      csv << std::setprecision(12);
      csv << "t,lv_y,lv_y_epi,lv_phi,lv_pv,lv_par,lv_pd,lv_pf,lv_lambda_bar"
          << ",q_in,q_out_total,q_tissue,q_lumped";
      for (const Attribute outlet : BoundaryFluid::Outlets)
        csv << ",q_out_" << outlet << ",p_out_" << outlet;
      csv << ",ptm_mean,phi_a_mean,phi_v_mean,stored_volume";
      csv << ",e_interface,coupling_rel\n";
    }

    std::map<Attribute, Outlet> wk;
    for (const Attribute outlet : BoundaryFluid::Outlets)
      wk.emplace(outlet, Outlet{});

    const WRMSTable wrms = buildWRMSTable(cfg, cfg.viscosity);

    const Real PI = std::numbers::pi_v<Real>;
    std::map<Attribute, Real> rEq;
    Real sumR3 = 0.0;
    {
      PETSc::Variational::TestFunction qCal(ph);
      LinearForm<PressureFES, ::Vec> calForm(qCal);
      for (const Attribute tag : BoundaryFluid::Outlets)
      {
        calForm = BoundaryIntegral(one, qCal).over(tag);
        calForm.assemble();
        const Real area = std::max<Real>(calForm(one), 1e-12);
        wk.at(tag).area = area;
        rEq[tag] = std::sqrt(area / PI);
        sumR3 += rEq[tag] * rEq[tag] * rEq[tag];
      }
    }

    // Outlet calibration on the resting budget dP = p_ar(0) - p_f(0):
    //   R_v,i = f_v dP/Q_i, R_a,i = (1 - f_v) dP/Q_i, C_i = C_tot w_i,
    //   p_tm,i(0) = f_v dP Phi_v0 (steady state of the calibrated network).
    if (cfg.autoCalibrateOutlets)
    {
      const Real dP = std::max<Real>(par0 - pf0, 1.0);
      const Real fv = cfg.venularPressureFraction;
      const Real dPv = std::max<Real>(fv * dP, 1.0);
      const Real dPa = std::max<Real>((1.0 - fv) * dP, 1.0);
      const Real muN = std::max<Real>(cfg.newtonianCalibrationViscosity, 1e-300);

      // Morphometric operating point (r, v) per limb: g_0 = 4v/r,
      // L = r dP_share/(2 mu_N g_0), T = L/v.
      const Real ra = std::max<Real>(cfg.arteriolarRadius, 1e-12);
      const Real va = std::max<Real>(cfg.arteriolarVelocity, 1e-12);
      const Real rv = std::max<Real>(cfg.venularRadius, 1e-12);
      const Real vv = std::max<Real>(cfg.venularVelocity, 1e-12);

      const Real gammaA0 = 4.0 * va / ra;
      const Real gammaV0 = 4.0 * vv / rv;
      const Real Ta = ra * dPa / (2.0 * muN * gammaA0) / va;
      const Real Tv = rv * dPv / (2.0 * muN * gammaV0) / vv;

      const Real phiA0 = wrms(gammaA0) / muN;
      const Real phiV0 = wrms(gammaV0) / muN;
      const Real ptmRest = dPv * phiV0;

      for (const Attribute tag : BoundaryFluid::Outlets)
      {
        const Real wgt = (rEq[tag] * rEq[tag] * rEq[tag]) / std::max(sumR3, 1e-30);
        const Real Qi = std::max<Real>(cfg.lcaTargetFlow * wgt, 1e-12);

        auto& bc = wk.at(tag);
        bc.q0 = Qi;
        bc.Ra = dPa / Qi;
        bc.Rv = dPv / Qi;
        bc.C = std::max<Real>(cfg.coronaryComplianceTotal * wgt, 1e-300);
        bc.gammaA = gammaA0;
        bc.gammaV = gammaV0;
        bc.ptm = ptmRest;
        bc.pf = pf0;
        bc.pc = ptmRest + pf0;
        bc.pout = bc.pc;
        bc.qv = Qi;
        bc.vol = bc.C * std::max<Real>(ptmRest, 0.0);
        bc.phiA = phiA0;
        bc.phiV = phiV0;

        if (isRoot)
          Alert::Info() << "  [calib] outlet " << tag << "  A=" << bc.area << " m^2"
                        << "  Q=" << (Qi * 6.0e7) << " mL/min" << "  Ra=" << bc.Ra
                        << "  Rv=" << bc.Rv << " Pa s/m^3" << "  C=" << bc.C
                        << " m^3/Pa" << "  tau=C*Rv=" << (bc.C * bc.Rv) << " s"
                        << "  ptm0=" << ptmRest << "  pc0=" << bc.pc
                        << " Pa (par=" << par0 << ", pf=" << pf0 << ")" << Alert::Raise;
      }

      if (isRoot)
      {
        Alert::Info() << "  [calib] WRMS nodes=" << wrms.logGamma.size() << "  rheology="
                      << (cfg.rheologyModel == RheologyModel::Quemada ? "Quemada"
                                                                       : "Carreau-Yasuda")
                      << "  Phi_a0=" << phiA0 << "  Phi_v0=" << phiV0
                      << "  |  T_a=" << Ta << " s  T_v=" << Tv << " s  (reference "
                      << cfg.referenceTransitTime << " s)" << Alert::Raise;

        const Real Tref = cfg.referenceTransitTime;
        if (Ta < 0.5 * Tref || Ta > 2.0 * Tref || Tv < 0.5 * Tref || Tv > 2.0 * Tref)
          Alert::Warning()
            << "  [calib] derived transit time departs more than 2x from the "
            << "reference: radius, velocity and pressure split are not consistent "
            << "with a single effective segment per limb." << Alert::Raise;
      }
    }
    else
    {
      for (const Attribute tag : BoundaryFluid::Outlets)
      {
        auto& bc = wk.at(tag);
        const Real Qi = cfg.lcaTargetFlow / static_cast<Real>(BoundaryFluid::Outlets.size());
        bc.q0 = Qi;
        bc.Ra = 4.5e9;
        bc.Rv = 6.8e8;
        bc.C = cfg.coronaryComplianceTotal / static_cast<Real>(BoundaryFluid::Outlets.size());
        bc.ptm = 1400.0;
        bc.pf = pf0;
        bc.pc = bc.ptm + pf0;
        bc.pout = bc.pc;
        bc.qv = Qi;
        bc.gammaA = 4.0 * cfg.arteriolarVelocity / cfg.arteriolarRadius;
        bc.gammaV = 4.0 * cfg.venularVelocity / cfg.venularRadius;
        bc.vol = bc.C * std::max<Real>(bc.ptm, 0.0);
      }
    }

    Real pinValue = par0;
    std::map<Attribute, Real> outletPressureValue;
    for (const Attribute outlet : BoundaryFluid::Outlets)
      outletPressureValue[outlet] = wk.at(outlet).pout;

    auto pin = RealFunction([&](const Point&) { return pinValue; });
    auto pout0 = RealFunction(
      [&](const Point&) { return outletPressureValue[BoundaryFluid::Outlets[0]]; });
    auto pout1 = RealFunction(
      [&](const Point&) { return outletPressureValue[BoundaryFluid::Outlets[1]]; });
    auto pout2 = RealFunction(
      [&](const Point&) { return outletPressureValue[BoundaryFluid::Outlets[2]]; });
    auto pout3 = RealFunction(
      [&](const Point&) { return outletPressureValue[BoundaryFluid::Outlets[3]]; });
    auto pout4 = RealFunction(
      [&](const Point&) { return outletPressureValue[BoundaryFluid::Outlets[4]]; });
    auto pout5 = RealFunction(
      [&](const Point&) { return outletPressureValue[BoundaryFluid::Outlets[5]]; });

    const Real betaN = cfg.newmarkBeta;
    const Real gammaN = cfg.newmarkGamma;
    const Real solidMass = cfg.solidDensity / (betaN * dt * dt);
    const Real solidVelocityCoeff = gammaN / (betaN * dt);

    const Real yeohC1 = 80000.0;
    const Real yeohC2 = 400000.0;
    const Real yeohC3 = 5000000.0;
    const Real yeohKappa = 12000000.0;

    const Real solidShearEquiv = 2.0 * yeohC1;
    const Real solidYoungEquiv =
      9.0 * yeohKappa * solidShearEquiv / (3.0 * yeohKappa + solidShearEquiv);

    // Robin parameter alpha = gamma sqrt(rho_s E_eq) unless overridden.
    const Real robinAlpha = (cfg.robinAlpha > 0.0)
      ? cfg.robinAlpha
      : cfg.robinGamma * std::sqrt(cfg.solidDensity * solidYoungEquiv);
    if (isRoot)
      Alert::Info() << "Robin parameter alpha = " << robinAlpha
                    << "  (gamma * sqrt(rho_s E_eq), gamma = " << cfg.robinGamma
                    << ", E_eq = " << solidYoungEquiv << " Pa)" << Alert::Raise;
    const Real robinVelocityCoeff = robinAlpha * solidVelocityCoeff;

    const auto& cy = cfg.viscosity;
    const Real gammaReg = cy.gammaRegularization;
    const Real deltaMu = cy.mu0 - cy.muInf;

    const auto normalFluid = BoundaryNormal(meshFluid);

    // Graded Yeoh: m(xi) blends intima/media/adventitia with smoothsteps at
    // xi = 1/3, 2/3; the energy is linear in (c1, c2, c3, kappa), so the
    // tangent stays consistent.
    auto wallGrade = [&](const Point& pt) -> Real {
      const Real s = std::clamp<Real>(xi.getValue(pt), 0.0, 1.0);
      const Real hw = std::max<Real>(cfg.gradeTransitionWidth, 1.0e-6);
      const Real f1 = smoothstep((s - (1.0 / 3.0 - hw)) / (2.0 * hw));
      const Real f2 = smoothstep((s - (2.0 / 3.0 - hw)) / (2.0 * hw));
      return cfg.gradeIntima + (cfg.gradeMedia - cfg.gradeIntima) * f1 +
        (cfg.gradeAdventitia - cfg.gradeMedia) * f2;
    };

    Rodin::Examples::Heart::GradedYeoh law(
      Solid::Yeoh(yeohC1, yeohC2, yeohC3, yeohKappa), wallGrade);
    Solid::InternalVirtualWorkTangent solidTangent(law, d, w, dState);
    Solid::InternalVirtualWorkResidual solidInternal(law, w, dState);

    PETSc::Variational::GridFunction uCur(uh);
    PETSc::Variational::GridFunction pCur(ph);
    uCur = zero;
    pCur = 0.0;

    // ------------------------------------------------------------------
    // L2 gradient recovery: gradRec_i = Pi_h[row i of grad u].  A P1
    // gradient is elementwise constant; the recovered nodal field is what
    // the wall traction and the WSS are built from, pointwise.
    // ------------------------------------------------------------------
    PETSc::Variational::TrialFunction gradRecTrial(uh);
    PETSc::Variational::TestFunction gradRecTest(uh);
    const auto jacRow0 = VectorFunction(Component(Jacobian(uCur), 0, 0),
      Component(Jacobian(uCur), 0, 1), Component(Jacobian(uCur), 0, 2));
    const auto jacRow1 = VectorFunction(Component(Jacobian(uCur), 1, 0),
      Component(Jacobian(uCur), 1, 1), Component(Jacobian(uCur), 1, 2));
    const auto jacRow2 = VectorFunction(Component(Jacobian(uCur), 2, 0),
      Component(Jacobian(uCur), 2, 1), Component(Jacobian(uCur), 2, 2));
    Problem gradRecProj0(gradRecTrial, gradRecTest);
    gradRecProj0 = Integral(gradRecTrial, gradRecTest) - Integral(jacRow0, gradRecTest);
    Problem gradRecProj1(gradRecTrial, gradRecTest);
    gradRecProj1 = Integral(gradRecTrial, gradRecTest) - Integral(jacRow1, gradRecTest);
    Problem gradRecProj2(gradRecTrial, gradRecTest);
    gradRecProj2 = Integral(gradRecTrial, gradRecTest) - Integral(jacRow2, gradRecTest);

    auto solveMass = [](auto& problem) {
      problem.assemble();
      Solver::KSP ksp(problem);
      ksp.setPrefix("coronary_mass_");
      ksp.solve();
    };

    auto recoverGradient = [&]() {
      solveMass(gradRecProj0);
      gradRec0.setData(gradRecTrial.getSolution().getData());
      solveMass(gradRecProj1);
      gradRec1.setData(gradRecTrial.getSolution().getData());
      solveMass(gradRecProj2);
      gradRec2.setData(gradRecTrial.getSolution().getData());
    };

    // Pointwise wall stress sample at a boundary point of the fluid mesh:
    //   G = recovered grad u,  D = (G + G^T)/2,  gamma = sqrt(2 D:D),
    //   mu = mu_CY(gamma),  t_visc = 2 mu D n,
    //   traction = p n - t_visc (= -sigma_f n_f),  shear = t_visc - (t_visc.n) n.
    struct WallSample
    {
        Math::SpatialVector<Real> traction;
        Math::SpatialVector<Real> shear;
    };
    auto wallSampleAt = [&](const Point& pt) -> WallSample {
      const Math::SpatialVector<Real> nRaw = normalFluid.getValue(pt);
      const auto g0 = gradRec0(pt);
      const auto g1 = gradRec1(pt);
      const auto g2 = gradRec2(pt);

      Math::SpatialVector<Real> n(3);
      Math::SpatialMatrix<Real> G(3, 3);
      for (std::uint8_t j = 0; j < 3; ++j)
      {
        n(j) = nRaw(j);
        G(0, j) = g0(j);
        G(1, j) = g1(j);
        G(2, j) = g2(j);
      }

      Math::SpatialMatrix<Real> D(3, 3);
      Real dd = 0.0;
      for (std::uint8_t i = 0; i < 3; ++i)
        for (std::uint8_t j = 0; j < 3; ++j)
        {
          D(i, j) = 0.5 * (G(i, j) + G(j, i));
          dd += D(i, j) * D(i, j);
        }

      const Real shearRate = std::sqrt(gammaReg * gammaReg + 2.0 * dd);
      const Real mu = cy(shearRate);
      const Real pw = pCur.getValue(pt);

      WallSample s;
      s.shear.resize(3);
      s.traction.resize(3);
      Real tn = 0.0;
      for (std::uint8_t i = 0; i < 3; ++i)
      {
        Real Dn = 0.0;
        for (std::uint8_t j = 0; j < 3; ++j)
          Dn += D(i, j) * n(j);
        s.shear(i) = 2.0 * mu * Dn;
        s.traction(i) = pw * n(i) - s.shear(i);
        tn += s.shear(i) * n(i);
      }
      for (std::uint8_t i = 0; i < 3; ++i)
        s.shear(i) -= tn * n(i);

      if (!isFiniteVec(s.traction) || !isFiniteVec(s.shear))
      {
        static bool reported = false;
        if (!reported)
        {
          reported = true;
          Alert::Warning() << "wallSampleAt: non-finite wall sample; using zero "
                           << "there (warned once)." << Alert::Raise;
        }
        s.traction.setZero();
        s.shear.setZero();
      }
      return s;
    };

    auto tractionFSI = VectorFunction(dim, [&](const Point& pt) {
      return wallSampleAt(pt).traction;
    });
    auto wallShear = VectorFunction(dim, [&](const Point& pt) {
      return wallSampleAt(pt).shear;
    });

    // Areal stretch A_t/A_0 of a solid FSI face at the iterate dIter.
    auto arealStretchAt = [&](const Point& xs) -> Real {
      Real stretch = 1.0;
      const auto& verts = xs.getPolytope().getVertices();
      if (dim == 3 && verts.size() == 3)
      {
        Math::SpatialPoint X0 = meshSolid.getVertexCoordinates(verts[0]);
        Math::SpatialPoint X1 = meshSolid.getVertexCoordinates(verts[1]);
        Math::SpatialPoint X2 = meshSolid.getVertexCoordinates(verts[2]);

        Math::SpatialPoint x0 = X0, x1 = X1, x2 = X2;
        for (Index c = 0; c < 3; ++c)
        {
          x0(c) += dIter[dh.getGlobalIndex({0, verts[0]}, c)];
          x1(c) += dIter[dh.getGlobalIndex({0, verts[1]}, c)];
          x2(c) += dIter[dh.getGlobalIndex({0, verts[2]}, c)];
        }

        const auto triArea = [](const Math::SpatialPoint& a, const Math::SpatialPoint& b,
                               const Math::SpatialPoint& c) -> Real {
          const Math::SpatialPoint e1 = b - a;
          const Math::SpatialPoint e2 = c - a;
          const Real nx = e1(1) * e2(2) - e1(2) * e2(1);
          const Real ny = e1(2) * e2(0) - e1(0) * e2(2);
          const Real nz = e1(0) * e2(1) - e1(1) * e2(0);
          return 0.5 * std::sqrt(nx * nx + ny * ny + nz * nz);
        };

        const Real A0 = triArea(X0, X1, X2);
        const Real At = triArea(x0, x1, x2);
        if (A0 > 0.0)
          stretch = At / A0;
      }
      return stretch;
    };

    // Fluid traction on the solid, pulled back to the reference area.
    auto fluidStress = VectorFunction(dim, [&](const Point& xs) {
      const Point xf = forwardSolidPointToFluid(xs, meshFluid, interfaceMap);

      Math::SpatialVector<Real> value(dim);
      value.setZero();
      const auto force = tractionTransfer(xf);
      for (Index i = 0; i < static_cast<Index>(dim); ++i)
        value(i) = force(i);

      if (!isFiniteVec(value))
      {
        static bool reported = false;
        if (!reported)
        {
          reported = true;
          Alert::Warning() << "fluidStress: non-finite cross-mesh traction sample; "
                           << "using zero there (warned once)." << Alert::Raise;
        }
        value.setZero();
        return value;
      }

      value *= arealStretchAt(xs);
      return value;
    });

    // Interface velocity of the fluid Robin datum: the BDF1 mesh velocity
    // (GCL-consistent) or the solid Newmark velocity.
    auto interfaceSolidVelocity = VectorFunction(dim, [&](const Point& xf) {
      Math::SpatialVector<Real> value(dim);
      value.setZero();

      if (cfg.meshConsistentInterfaceVelocity)
      {
        const auto wv = meshVelocity(xf);
        for (Index i = 0; i < static_cast<Index>(dim); ++i)
          value(i) = wv(i);
      }
      else
      {
        const Point xs = forwardFluidPointToSolid(xf, meshSolid, interfaceMap);
        const auto us = solidVelocity(xs);
        for (Index i = 0; i < static_cast<Index>(dim); ++i)
          value(i) = us(i);
      }

      if (!isFiniteVec(value))
      {
        static bool reported = false;
        if (!reported)
        {
          reported = true;
          Alert::Warning() << "interfaceSolidVelocity: non-finite sample; using zero "
                           << "there (warned once)." << Alert::Raise;
        }
        value.setZero();
      }
      return value;
    });

    const auto transportLag = uOld - meshVelocity;
    const auto convU = Mult(Jacobian(u), transportLag);
    // Conservative two-mesh BDF1: Temam term + geometric companion,
    // coefficient (1/2) div(u^n) - div(w).
    const auto divGeomTemam = Div(uOld) - 2.0 * Div(meshVelocity);

    const auto duNormal = Dot(u, normalFluid) * normalFluid;
    const auto duTangential = u - duNormal;

    const auto symU = 0.5 * (Jacobian(u) + Transpose(Jacobian(u)));
    const auto symV = 0.5 * (Jacobian(v) + Transpose(Jacobian(v)));
    const auto symLag = 0.5 * (Jacobian(uOld) + Transpose(Jacobian(uOld)));
    const auto shearLag = Sqrt(gammaReg * gammaReg + 2.0 * Dot(symLag, symLag));
    const auto muLag = cy.muInf +
      deltaMu * Pow(1.0 + Pow(cy.lambda * shearLag, cy.yasuda), (cy.n - 1.0) / cy.yasuda);

    const auto outletBeta = Max(-Dot(transportLag, normalFluid), 0.0);
    const auto inletBeta = Max(Dot(transportLag, normalFluid), 0.0);
    const auto outletBackflow =
      0.5 * cfg.outletBackflowStabilization * cfg.fluidDensity * outletBeta;
    const auto inletBackflow =
      0.5 * cfg.inletBackflowStabilization * cfg.fluidDensity * inletBeta;

    using namespace Rodin::Examples::Heart;

    // ALE convecting velocity u^n - w (refreshed every coupling iterate).
    PETSc::Variational::GridFunction uConv(uh);
    uConv = uOld;

    PressureFES tauFes(std::integral_constant<size_t, 1>{}, meshFluid);
    PETSc::Variational::TestFunction vmsScalarTest(tauFes);
    PETSc::Variational::TrialFunction vmsPiTilde(tauFes);

    PETSc::Variational::TrialFunction vmsUp(uh);
    PETSc::Variational::TestFunction vmsVp(uh);
    PETSc::Variational::TrialFunction vmsSub(uh);
    PETSc::Variational::GridFunction vmsSubOld(uh);
    vmsSubOld = zero;

    const auto vmsConvectionTarget = Mult(Jacobian(uConv), uConv);

    // ------------------------------------------------------------------
    // Stabilization parameters, POINTWISE at every quadrature point:
    //   nu(x)  = mu_CY(gamma(u^n)(x)) / rho          (local lagged viscosity)
    //   h_K    = element diameter (longest edge) of the CURRENT configuration
    //   tau_1  = 1 / (c1 nu/h_K^2 + c2 |u^n - w|/h_K),  c1 = 4, c2 = 2 (P1)
    //   tau_K  = vmsScale / (rho/dt + rho/tau_1)           (convective subscale)
    //   tau_C  = gradDivScale rho h_K^2 / (4 tau_1)        (grad-div)
    //   tau_p  = pgpScale tau_1 / rho                      (PSPG)
    // h_K is cached per cell and refreshed after every mesh move.
    // ------------------------------------------------------------------
    std::vector<Real> hCell;
    computeCellDiameters(meshFluid, hCell);

    auto cellSizeAt = [&](const Point& pp) -> Real {
      const auto& poly = pp.getPolytope();
      if (poly.getDimension() == meshFluid.getDimension() &&
        poly.getIndex() < hCell.size())
        return hCell[poly.getIndex()];
      return std::pow(poly.getMeasure(), 1.0 / poly.getDimension());
    };

    auto viscosityAt = [&](const Point& pp) -> Real {
      const auto sym = 0.5 * (Jacobian(uOld) + Transpose(Jacobian(uOld)));
      const Real shear = std::sqrt(gammaReg * gammaReg + 2.0 * Dot(sym, sym).getValue(pp));
      return cy(shear);
    };

    auto tau1At = [&](const Point& pp) -> Real {
      const auto uc = uConv.getValue(pp);
      const Real nu = viscosityAt(pp) / cfg.fluidDensity;
      const Real hK = std::max<Real>(cellSizeAt(pp), 1.0e-30);
      const Real speed = std::sqrt(Math::dot(uc, uc));
      return 1.0 / (4.0 * nu / (hK * hK) + 2.0 * speed / hK);
    };

    auto vmsTauAt = [&](const Point& pp) -> Real {
      return cfg.vmsScale / (cfg.fluidDensity / dt + cfg.fluidDensity / tau1At(pp));
    };
    RealFunction vmsTauFn = [&](const Point& pp) -> Real { return vmsTauAt(pp); };

    auto sqrtTauCAt = [&](const Point& pp) -> Real {
      const Real hK = cellSizeAt(pp);
      return std::sqrt(std::max<Real>(
        0.0, cfg.gradDivScale * cfg.fluidDensity * hK * hK / (4.0 * tau1At(pp))));
    };
    RealFunction sqrtTauCFn = [&](const Point& pp) -> Real { return sqrtTauCAt(pp); };
    RealFunction tauCFn = [&](const Point& pp) -> Real {
      const Real s = sqrtTauCAt(pp);
      return s * s;
    };

    RealFunction vmsTauPFn = [&](const Point& pp) -> Real {
      return cfg.pgpScale * tau1At(pp) / cfg.fluidDensity;
    };

    // Dynamic subscale u'^{n+1} = tau_K rho (u'^n/dt - ((grad u^n)u^n - Pi[...])).
    auto vmsSubUpdate =
      VectorFunction(dim, [&](const Point& pp) -> Math::SpatialVector<Real> {
        const auto conv = vmsConvectionTarget.getValue(pp);
        const auto proj = vmsUp.getSolution().getValue(pp);
        const auto old = vmsSubOld.getValue(pp);
        const Real tau = vmsTauAt(pp);

        Math::SpatialVector<Real> out(dim);
        for (Index c = 0; c < static_cast<Index>(dim); ++c)
          out(c) = tau * cfg.fluidDensity * (1.0 / dt * old(c) - (conv(c) - proj(c)));
        return out;
      });

    Problem vmsL2Conv(vmsUp, vmsVp);
    vmsL2Conv = Integral(vmsUp, vmsVp) - Integral(vmsConvectionTarget, vmsVp);

    Problem vmsSubProj(vmsSub, vmsVp);
    vmsSubProj = Integral(vmsSub, vmsVp) - Integral(vmsSubUpdate, vmsVp);

    // pi~ = Pi[sqrt(tau_C) div u^n]; the same sqrt(tau_C) multiplies div v
    // in the linear term and squares into the implicit coefficient.
    Problem vmsPiTildeProj(vmsPiTilde, vmsScalarTest);
    vmsPiTildeProj = Integral(vmsPiTilde, vmsScalarTest) -
      Integral(sqrtTauCFn * Div(uOld), vmsScalarTest);

    // Implicit outlet resistance R_a Phi_a A (u.n)(v.n): the resistive part
    // of p_out is not lagged, so the 0D-3D coupling is stable for any dt.
    auto outletZAt = [&](size_t i) -> Real {
      const auto& bc = wk.at(BoundaryFluid::Outlets[i]);
      return cfg.outletResistanceScale * bc.Ra * bc.phiA * bc.area;
    };
    auto zFn0 = RealFunction([&](const Point&) { return outletZAt(0); });
    auto zFn1 = RealFunction([&](const Point&) { return outletZAt(1); });
    auto zFn2 = RealFunction([&](const Point&) { return outletZAt(2); });
    auto zFn3 = RealFunction([&](const Point&) { return outletZAt(3); });
    auto zFn4 = RealFunction([&](const Point&) { return outletZAt(4); });
    auto zFn5 = RealFunction([&](const Point&) { return outletZAt(5); });

    // BDF1 mass split: implicit (rho/dt)(u, v)_{n+1} here, the explicit
    // (rho/dt)(u^n, v)_n assembled on the previous configuration (massOld).
    Problem flow(u, p, v, q);
    flow = (cfg.fluidDensity / dt) * Integral(u, v) +
      cfg.fluidDensity * Integral(Dot(convU, v)) +
      0.5 * cfg.fluidDensity * Integral(divGeomTemam * Dot(u, v)) +
      VMSConvectionBilinearIntegrator(u, v, uConv, vmsTauFn, cfg.fluidDensity) -
      VMSConvectionLinearIntegrator(v, vmsSub.getSolution(), uConv, vmsUp.getSolution(),
        vmsTauFn, cfg.fluidDensity, dt) +
      VMSGradDivBilinearIntegrator(u, v, tauCFn) -
      VMSGradDivLinearIntegrator(v, vmsPiTilde.getSolution(), sqrtTauCFn) +
      2.0 * Integral(muLag * symU, symV) - Integral(p, Div(v)) + Integral(Div(u), q) +
      Integral(vmsTauPFn * Grad(p), Grad(q)) +
      BoundaryIntegral(inletBackflow * Dot(u, v)).over(BoundaryFluid::Inlet) +
      BoundaryIntegral(outletBackflow * Dot(u, v))
        .over(BoundaryFluid::Outlets[0], BoundaryFluid::Outlets[1],
          BoundaryFluid::Outlets[2], BoundaryFluid::Outlets[3], BoundaryFluid::Outlets[4],
          BoundaryFluid::Outlets[5]) +
      BoundaryIntegral(pin * Dot(v, normalFluid)).over(BoundaryFluid::Inlet) +
      BoundaryIntegral(pout0 * Dot(v, normalFluid)).over(BoundaryFluid::Outlets[0]) +
      BoundaryIntegral(pout1 * Dot(v, normalFluid)).over(BoundaryFluid::Outlets[1]) +
      BoundaryIntegral(pout2 * Dot(v, normalFluid)).over(BoundaryFluid::Outlets[2]) +
      BoundaryIntegral(pout3 * Dot(v, normalFluid)).over(BoundaryFluid::Outlets[3]) +
      BoundaryIntegral(pout4 * Dot(v, normalFluid)).over(BoundaryFluid::Outlets[4]) +
      BoundaryIntegral(pout5 * Dot(v, normalFluid)).over(BoundaryFluid::Outlets[5]) +
      BoundaryIntegral(zFn0 * Dot(Dot(u, normalFluid) * normalFluid, v))
        .over(BoundaryFluid::Outlets[0]) +
      BoundaryIntegral(zFn1 * Dot(Dot(u, normalFluid) * normalFluid, v))
        .over(BoundaryFluid::Outlets[1]) +
      BoundaryIntegral(zFn2 * Dot(Dot(u, normalFluid) * normalFluid, v))
        .over(BoundaryFluid::Outlets[2]) +
      BoundaryIntegral(zFn3 * Dot(Dot(u, normalFluid) * normalFluid, v))
        .over(BoundaryFluid::Outlets[3]) +
      BoundaryIntegral(zFn4 * Dot(Dot(u, normalFluid) * normalFluid, v))
        .over(BoundaryFluid::Outlets[4]) +
      BoundaryIntegral(zFn5 * Dot(Dot(u, normalFluid) * normalFluid, v))
        .over(BoundaryFluid::Outlets[5]) +
      cfg.inletImpedance *
        BoundaryIntegral(Dot(Dot(u, normalFluid) * normalFluid, v))
          .over(BoundaryFluid::Inlet) +
      cfg.inletTangentialDamping *
        BoundaryIntegral(Dot(duTangential, v)).over(BoundaryFluid::Inlet)
      // Robin-Robin (fluid side): sigma_f n + alpha u = alpha d_dot_s + lambda^{k-1}.
      + robinAlpha * BoundaryIntegral(u, v).over(BoundaryFluid::FSI) -
      robinAlpha * BoundaryIntegral(interfaceSolidVelocity, v).over(BoundaryFluid::FSI) +
      BoundaryIntegral(tractionFSI, v).over(BoundaryFluid::FSI)
      // Interface convective stabilization -(rho/2)((u^n - w).n)(u.v).
      - 0.5 * cfg.fluidDensity *
        BoundaryIntegral(Dot(transportLag, normalFluid) * Dot(u, v))
          .over(BoundaryFluid::FSI) +
      DirichletBC(u, zero).on(BoundaryFluid::FSIRing);

    PETSc::Variational::TestFunction vMass(uh);
    LinearForm<VelocityFES, ::Vec> massOld(vMass);
    massOld = (cfg.fluidDensity / dt) * Integral(uOld, vMass);

    auto interfaceSolidDisplacement = VectorFunction(dim, [&](const Point& xf) {
      const Point xs = forwardFluidPointToSolid(xf, meshSolid, interfaceMap);

      Math::SpatialVector<Real> value(dim);
      value.setZero();

      const auto ds = dIter(xs);
      for (Index i = 0; i < static_cast<Index>(dim); ++i)
        value(i) = ds(i);

      if (!isFiniteVec(value))
      {
        static bool reported = false;
        if (!reported)
        {
          reported = true;
          Alert::Warning() << "interfaceSolidDisplacement: non-finite cross-mesh "
                           << "sample; using zero there (warned once)." << Alert::Raise;
        }
        value.setZero();
      }
      return value;
    });

    // Harmonic ALE lift on the reference configuration, element-size stiffened.
    PETSc::Variational::TrialFunction dMove(uh);
    PETSc::Variational::TestFunction vMove(uh);

    RealFunction aleStiffFn = [&](const Point& pp) -> Real {
      if (cfg.aleStiffPower <= 0.0)
        return 1.0;
      const Real hK =
        std::pow(pp.getPolytope().getMeasure(), 1.0 / pp.getPolytope().getDimension());
      return std::pow(cfg.aleRefSize / std::max(hK, 1.0e-30), cfg.aleStiffPower);
    };

    Problem ale(dMove, vMove);
    ale = Integral(aleStiffFn * Jacobian(dMove), Jacobian(vMove)) +
      DirichletBC(dMove, interfaceSolidDisplacement).on(BoundaryFluid::FSI) +
      DirichletBC(dMove, zero).on(BoundaryFluid::Inlet, BoundaryFluid::FSIRing);

    // Solid Robin datum (per reference area):
    //   rVC (dState - dPred) + alpha vPred - alpha u_f^{k-1}.
    auto robinInterfaceData = VectorFunction(dim, [&](const Point& xs) {
      const Point xf = forwardSolidPointToFluid(xs, meshFluid, interfaceMap);

      Math::SpatialVector<Real> value(dim);
      value.setZero();

      Math::SpatialVector<Real> uf(dim);
      uf.setZero();
      {
        const auto ufRaw = uWall(xf);
        for (Index i = 0; i < static_cast<Index>(dim); ++i)
          uf(i) = ufRaw(i);
        if (!isFiniteVec(uf))
        {
          static bool reported = false;
          if (!reported)
          {
            reported = true;
            Alert::Warning() << "robinInterfaceData: non-finite cross-mesh velocity "
                             << "sample; using zero there (warned once)." << Alert::Raise;
          }
          uf.setZero();
        }
      }
      const auto vp = vPred(xs);
      const auto dS = dState(xs);
      const auto dP = dPred(xs);
      for (Index i = 0; i < static_cast<Index>(dim); ++i)
        value(i) = robinVelocityCoeff * (dS(i) - dP(i)) + robinAlpha * vp(i) -
          robinAlpha * uf(i);
      return value;
    });

    // 0D heart motion: weak d.n = disp_0D on the contact patch, weighted by
    // the heart-contact laplacian l.
    const auto normalSolid = BoundaryNormal(meshSolid);
    Real disp0D = 0.0;
    Real disp0DOffset = 0.0;
    auto disp0DFn = RealFunction([&](const Point&) { return -disp0D; });
    const Real heartK = cfg.heartDisplacementPenalty;

    const Real solidViscImpl = cfg.solidViscosity * solidVelocityCoeff;

    const Real a = cfg.aViscCondition;
    const Real b = cfg.bViscCondition;
    const Real aVel = b * gammaN / (betaN * dt);

    Problem solid(d, w);
    solid = solidMass * Integral(d, w) + solidTangent + solidMass * Integral(dState, w) -
      solidMass * Integral(dPred, w) + solidInternal
      // Kelvin-Voigt damping with the Newmark velocity vPred + rVC (dState - dPred).
      + solidViscImpl * Integral(Jacobian(d), Jacobian(w)) +
      solidViscImpl * Integral(Jacobian(dState), Jacobian(w)) -
      solidViscImpl * Integral(Jacobian(dPred), Jacobian(w)) +
      cfg.solidViscosity * Integral(Jacobian(vPred), Jacobian(w)) +
      DirichletBC(d, zero).on(BoundarySolid::Inlet) +
      DirichletBC(d, zero).on(BoundarySolid::FSIRing) +
      a *
        BoundaryIntegral(d, w).over(BoundarySolid::Outlets[0], BoundarySolid::Outlets[1],
          BoundarySolid::Outlets[2], BoundarySolid::Outlets[3], BoundarySolid::Outlets[4],
          BoundarySolid::Outlets[5]) +
      aVel *
        BoundaryIntegral(d, w).over(BoundarySolid::Outlets[0], BoundarySolid::Outlets[1],
          BoundarySolid::Outlets[2], BoundarySolid::Outlets[3], BoundarySolid::Outlets[4],
          BoundarySolid::Outlets[5]) +
      a *
        BoundaryIntegral(dState, w).over(BoundarySolid::Outlets[0],
          BoundarySolid::Outlets[1], BoundarySolid::Outlets[2], BoundarySolid::Outlets[3],
          BoundarySolid::Outlets[4], BoundarySolid::Outlets[5]) +
      aVel *
        BoundaryIntegral(dState, w).over(BoundarySolid::Outlets[0],
          BoundarySolid::Outlets[1], BoundarySolid::Outlets[2], BoundarySolid::Outlets[3],
          BoundarySolid::Outlets[4], BoundarySolid::Outlets[5]) -
      aVel *
        BoundaryIntegral(dPred, w).over(BoundarySolid::Outlets[0],
          BoundarySolid::Outlets[1], BoundarySolid::Outlets[2], BoundarySolid::Outlets[3],
          BoundarySolid::Outlets[4], BoundarySolid::Outlets[5]) +
      b *
        BoundaryIntegral(vPred, w).over(BoundarySolid::Outlets[0],
          BoundarySolid::Outlets[1], BoundarySolid::Outlets[2], BoundarySolid::Outlets[3],
          BoundarySolid::Outlets[4], BoundarySolid::Outlets[5]) +
      heartK *
        BoundaryIntegral(Dot(d, normalSolid), Dot(w, normalSolid))
          .over(BoundarySolid::Contact[0]) +
      heartK *
        BoundaryIntegral(Dot(dState, normalSolid), Dot(w, normalSolid))
          .over(BoundarySolid::Contact[0]) -
      heartK *
        BoundaryIntegral(disp0DFn * l.getSolution(), Dot(w, normalSolid))
          .over(BoundarySolid::Contact[0])
      // Robin-Robin (solid side): sigma_s n_s + alpha d_dot = alpha u_f^lag + t_f^lag.
      + robinVelocityCoeff * BoundaryIntegral(Dot(d, w)).over(BoundarySolid::FSI) +
      BoundaryIntegral(robinInterfaceData, w).over(BoundarySolid::FSI) -
      BoundaryIntegral(fluidStress, w).over(BoundarySolid::FSI);

    solid.assemble();
    Solver::KSP kspSolid(solid);
    Solver::SNES snes(kspSolid);
    snes.setTolerances(1.0e-10, 1.0e-8, 1.0e-10, 50, 10000);
    snes.setStateUpdate([&](const PETSc::Math::Vector& state) {
      etaState.setData(state, 0);
      dState = dOld;
      dState += etaState;
    });

    // Heart-contact weight: harmonic between the inlet (0) and the outlets.
    Problem laplacian(l, t);
    laplacian = Integral(Grad(l), Grad(t)) +
      DirichletBC(l, RealFunction(0.0)).on(BoundarySolid::Inlet) +
      DirichletBC(l, RealFunction(0.95)).on(BoundarySolid::Outlets[3]) +
      DirichletBC(l, RealFunction(0.85)).on(BoundarySolid::Outlets[4]) +
      DirichletBC(l, RealFunction(0.75))
        .on(BoundarySolid::Outlets[0], BoundarySolid::Outlets[1],
          BoundarySolid::Outlets[2]) +
      DirichletBC(l, RealFunction(0.45)).on(BoundarySolid::Outlets[5]);

    if (isRoot)
    {
      PETSc::Variational::GridFunction oneSolid(lh);
      oneSolid = 1.0;
      LinearForm<LaplacianFES, ::Vec> solidArea(t);
      auto faceArea = [&](Attribute tag) -> Real {
        solidArea = BoundaryIntegral(oneSolid, t).over(tag);
        solidArea.assemble();
        return solidArea(oneSolid);
      };
      Alert::Info() << "  [solid-bdr] Inlet(" << BoundarySolid::Inlet
                    << ") area=" << faceArea(BoundarySolid::Inlet) << "  FSI("
                    << BoundarySolid::FSI << ") area=" << faceArea(BoundarySolid::FSI)
                    << "  FSIRing(" << BoundarySolid::FSIRing
                    << ") area=" << faceArea(BoundarySolid::FSIRing) << Alert::Raise;
      for (size_t i = 0; i < BoundarySolid::Outlets.size(); ++i)
        Alert::Info() << "  [solid-bdr] Outlet(" << BoundarySolid::Outlets[i]
                      << ") area=" << faceArea(BoundarySolid::Outlets[i]) << Alert::Raise;
    }

    laplacian.assemble();
    Solver::KSP(laplacian).solve();

    // xi must be solved before the prestress: the graded law reads it at
    // every quadrature point.
    Problem thickness(xiT, xiTest);
    thickness = Integral(Grad(xiT), Grad(xiTest)) +
      DirichletBC(xiT, RealFunction(0.0)).on(BoundarySolid::FSI) +
      DirichletBC(xiT, RealFunction(1.0))
        .on(BoundarySolid::Outer[0], BoundarySolid::Outer[1]);
    thickness.assemble();
    Solver::KSP(thickness).solve();
    xi.setData(xiT.getSolution().getData());
    if (isRoot)
      Alert::Info() << "  [grade] transmural xi solved; multipliers"
                    << "  intima=" << cfg.gradeIntima << "  media=" << cfg.gradeMedia
                    << "  adventitia=" << cfg.gradeAdventitia
                    << "  (transition half-width " << cfg.gradeTransitionWidth << ")"
                    << Alert::Raise;

    xdmf_laplacian.add("xi", xi);
    xdmf_laplacian.write().flush();
    xdmf_laplacian.close();

    if (cfg.prestressSteps > 0)
    {
      // Static follower-pressure prestress to prestressFraction * p_ar(0);
      // the dynamic loop ramps the remainder over prestressRampSteps.
      const Real p0 = cfg.prestressFraction * par0;
      Real prestressPressure = 0.0;

      PETSc::Variational::TrialFunction dPre(dh);
      PETSc::Variational::TestFunction wPre(dh);
      Solid::InternalVirtualWorkTangent preTangent(law, dPre, wPre, dState);
      Solid::InternalVirtualWorkResidual preInternal(law, wPre, dState);

      Solid::FollowerPressureForce preLoad(prestressPressure, wPre, dState);
      preLoad.over(BoundarySolid::FSI);
      Solid::FollowerPressureTangent preLoadK(prestressPressure, dPre, wPre, dState);
      preLoadK.over(BoundarySolid::FSI);

      Problem prestress(dPre, wPre);
      prestress = preTangent + preInternal + preLoadK + preLoad +
        a *
          BoundaryIntegral(dPre, wPre).over(BoundarySolid::Outlets[0],
            BoundarySolid::Outlets[1], BoundarySolid::Outlets[2],
            BoundarySolid::Outlets[3], BoundarySolid::Outlets[4],
            BoundarySolid::Outlets[5]) +
        a *
          BoundaryIntegral(dState, wPre).over(BoundarySolid::Outlets[0],
            BoundarySolid::Outlets[1], BoundarySolid::Outlets[2],
            BoundarySolid::Outlets[3], BoundarySolid::Outlets[4],
            BoundarySolid::Outlets[5]) +
        DirichletBC(dPre, zero).on(BoundarySolid::Inlet, BoundarySolid::FSIRing);

      prestress.assemble();
      Solver::KSP kspPre(prestress);
      Solver::SNES snesPre(kspPre);
      snesPre.setTolerances(1.0e-10, 1.0e-8, 1.0e-10, 50, 10000);
      snesPre.setStateUpdate([&](const PETSc::Math::Vector& state) {
        etaState.setData(state, 0);
        dState = dOld;
        dState += etaState;
      });

      for (size_t k = 1; k <= cfg.prestressSteps; ++k)
      {
        prestressPressure =
          (static_cast<Real>(k) / static_cast<Real>(cfg.prestressSteps)) * p0;
        snesPre.solve();
        if (!snesPre.converged())
        {
          if (isRoot)
            std::cerr << "Prestress SNES failed at increment " << k << " / "
                      << cfg.prestressSteps << "; continuing with the last "
                      << "converged (partial) prestress state.\n";
          break;
        }
      }
      if (isRoot)
        Alert::Info() << "Prestressed wall to " << prestressPressure << " Pa in "
                      << cfg.prestressSteps << " increment(s)" << Alert::Raise;

      dOld.setData(dState.getData());
      dIter.setData(dState.getData());
      restoreMeshToReference(meshFluid, referenceVertices);
      ale.assemble();
      Solver::KSP(ale).solve();
      aleDisp.setData(dMove.getSolution().getData());
      aleDispOld.setData(aleDisp.getData());
      moveMeshWithVertexDisplacement(meshFluid, referenceVertices, uh, aleDisp);
      computeCellDiameters(meshFluid, hCell);

      // Seed the fluid pressure at p0 so the transferred traction balances
      // the prestressed wall at the hand-off (recovered gradient is zero).
      pOld = p0;
      pCur = p0;
      tractionTransfer.project(Region::Faces, tractionFSI, BoundaryFluid::FSI);
    }

    using FluxLinearForm = LinearForm<PressureFES, ::Vec>;
    PETSc::Variational::TestFunction qFlux(ph);
    FluxLinearForm flux(qFlux);

    Real qIn = 0.0;
    std::map<Attribute, Real> qOut;
    for (const Attribute outlet : BoundaryFluid::Outlets)
      qOut[outlet] = 0.0;

    // Closed-loop hand-off: until the first 3D fluxes exist, the 0D model
    // keeps receiving the calibrated resting flows.
    qArterialLag = cfg.lcaTargetFlow;
    qTissueLag = 0.0;
    for (const auto& [tag, bc] : wk)
      qTissueLag += bc.qv;

    // ==================================================================
    // Time loop.  Step n -> n+1:
    //   (1) 0D poroelastic LV advanced with the lagged 3D fluxes;
    //   (2) p_in = p_ar^{n+1}, p_out,i = p_tm,i^n + p_f^{n+1};
    //   (3) coupling iterates: solid SNES -> ALE lift -> fluid Oseen;
    //   (4) fluxes, gradient recovery, WSS, outlet compartments with p_f^{n+1},
    //       lagged fluxes for the next 0D step, commit.
    // ==================================================================
    for (size_t step = 1; step <= cfg.nsteps; ++step)
    {
      const auto rep = model.step(dt);
      if (!rep.converged)
      {
        if (isRoot)
          std::cerr << "0D model failed at step " << step << '\n';
        break;
      }

      const auto& s = model.getState();
      const Real yEpi = epicardialDisplacement(modelInput, s);

      if (cfg.subtractHeartHandoffOffset)
      {
        if (step == 1)
          disp0DOffset = yEpi;
        disp0D = cfg.heartDisplacementScale * (yEpi - disp0DOffset);
      }
      else
      {
        disp0D = cfg.heartDisplacementScale * yEpi;
      }

      Real loadRamp = 1.0;
      if (cfg.prestressFraction < 1.0 && cfg.prestressRampSteps > 0)
      {
        const Real sRamp = std::min(Real(1.0),
          static_cast<Real>(step - 1) / static_cast<Real>(cfg.prestressRampSteps));
        loadRamp = cfg.prestressFraction + (1.0 - cfg.prestressFraction) * smoothstep(sRamp);
      }

      pinValue = loadRamp * s.par;
      for (const auto& [tag, bc] : wk)
        outletPressureValue[tag] = loadRamp * (bc.ptm + s.pf);

      // Newmark predictors.
      dPred = dOld;
      auto tmp = solidVelocityOld;
      tmp *= dt;
      dPred += tmp;
      tmp = solidAccelerationOld;
      tmp *= dt * dt * (0.5 - betaN);
      dPred += tmp;

      vPred = solidVelocityOld;
      tmp = solidAccelerationOld;
      tmp *= dt * (1.0 - gammaN);
      vPred += tmp;

      dIter = dPred;

      if (isRoot)
      {
        Alert::Info() << "Coronary explicit ALE FSI step " << step << " / " << cfg.nsteps
                      << "  (coupling iterations: " << cfg.couplingIterations << ")"
                      << Alert::Raise;
        Alert::Info() << "  [0D] t=" << s.t << " s  p_ar=" << s.par << "  p_v=" << s.pv
                      << "  p_f=" << s.pf << " Pa  phi=" << s.phi
                      << "  y_epi=" << yEpi << " m  disp0D=" << disp0D << " m"
                      << "  Q_in^lag=" << (qArterialLag * 6.0e7) << " mL/min"
                      << "  Q_tissue^lag=" << (qTissueLag * 6.0e7) << " mL/min"
                      << Alert::Raise;
      }

      Real qOutSum = 0.0;
      bool stepFailed = false;
      Real lastRel = 0.0;
      size_t couplesDone = 0;

      // Omega^n and u^n are fixed within the step: assemble massOld once.
      moveMeshWithVertexDisplacement(meshFluid, referenceVertices, uh, aleDispOld);
      massOld.assemble();

      for (size_t couple = 1; couple <= cfg.couplingIterations; ++couple)
      {
        snes.solve();
        if (!snes.converged())
        {
          if (isRoot)
            std::cerr << "Solid SNES failed at step " << step << " (coupling iterate "
                      << couple << ") after " << snes.getIterationNumber()
                      << " iterations.\n";
          stepFailed = true;
          break;
        }

        auto delta = dState;
        delta -= dIter;
        PetscReal deltaNorm = 0.0;
        PetscReal stateNorm = 0.0;
        VecNorm(delta.getData(), NORM_2, &deltaNorm);
        VecNorm(dState.getData(), NORM_2, &stateNorm);
        const Real rel = (stateNorm > 0.0)
          ? (static_cast<Real>(deltaNorm) / static_cast<Real>(stateNorm))
          : static_cast<Real>(deltaNorm);
        lastRel = rel;
        couplesDone = couple;
        dIter = dState;

        solidAcceleration = dIter;
        solidAcceleration -= dPred;
        solidAcceleration *= 1.0 / (betaN * dt * dt);
        solidVelocity = vPred;
        tmp = solidAcceleration;
        tmp *= gammaN * dt;
        solidVelocity += tmp;

        if (isRoot)
          Alert::Info() << "  coupling iterate " << couple << " / " << cfg.couplingIterations
                        << "  relative interface change = " << rel
                        << "  |d - dPrev| = " << deltaNorm << Alert::Raise;

        if (couple > 1 && rel < cfg.couplingTolerance)
          break;

        restoreMeshToReference(meshFluid, referenceVertices);
        ale.assemble();
        Solver::KSP(ale).solve();
        aleDisp.setData(dMove.getSolution().getData());

        {
          PetscReal aleNorm = 0.0;
          VecNorm(aleDisp.getData(), NORM_INFINITY, &aleNorm);
          if (!std::isfinite(static_cast<Real>(aleNorm)))
          {
            if (isRoot)
              std::cerr << "ALE lift non-finite at step " << step << " (coupling iterate "
                        << couple << "): the fluid mesh is tangling.\n";
            stepFailed = true;
            break;
          }
        }

        meshVelocity = aleDisp;
        meshVelocity -= aleDispOld;
        meshVelocity *= 1.0 / dt;

        moveMeshWithVertexDisplacement(meshFluid, referenceVertices, uh, aleDisp);
        computeCellDiameters(meshFluid, hCell);

        uConv = uOld;
        uConv -= meshVelocity;
        solveMass(vmsL2Conv);
        solveMass(vmsSubProj);
        solveMass(vmsPiTildeProj);

        flow.assemble().setFieldSplits();

        {
          ::Vec bvec = flow.getLinearSystem().getVector();
          const PetscInt vOff = static_cast<PetscInt>(flow.getTestOffsets()[0]);
          const ::Vec& mOld = massOld.getVector();
          PetscInt lo = 0, hi = 0;
          VecGetOwnershipRange(mOld, &lo, &hi);
          const PetscScalar* arr = nullptr;
          VecGetArrayRead(mOld, &arr);
          for (PetscInt i = lo; i < hi; ++i)
            if (arr[i - lo] != PetscScalar(0))
              VecSetValue(bvec, vOff + i, arr[i - lo], ADD_VALUES);
          VecRestoreArrayRead(mOld, &arr);
          VecAssemblyBegin(bvec);
          VecAssemblyEnd(bvec);
        }

        Solver::KSP(flow).solve();
        uCur.setData(u.getSolution().getData());
        pCur.setData(p.getSolution().getData());

        // Recovered gradient of the new iterate; the traction read by the
        // next solid solve and by the fluid Robin datum is built from it.
        recoverGradient();
        fluidTraction.project(Region::Faces, tractionFSI, BoundaryFluid::FSI);
        tractionTransfer.project(Region::Faces, tractionFSI, BoundaryFluid::FSI);
        uWall.project(Region::Faces, uCur, BoundaryFluid::FSI);
      }

      if (stepFailed)
        break;

      // Nodal WSS: area-weighted (lumped) boundary average of the pointwise
      // wall shear tau_w(x) on the FSI surface.
      {
        PETSc::Variational::TestFunction wssTest(uh);
        const auto onesVec = VectorFunction(dim, [&](const Point&) {
          Math::SpatialVector<Real> o(dim);
          for (Index c = 0; c < static_cast<Index>(dim); ++c)
            o(c) = 1.0;
          return o;
        });
        LinearForm<VelocityFES, ::Vec> wssLoad(wssTest);
        wssLoad = BoundaryIntegral(wallShear, wssTest).over(BoundaryFluid::FSI);
        wssLoad.assemble();
        LinearForm<VelocityFES, ::Vec> wssArea(wssTest);
        wssArea = BoundaryIntegral(onesVec, wssTest).over(BoundaryFluid::FSI);
        wssArea.assemble();

        ::Vec bvec = wssLoad.getVector();
        ::Vec mvec = wssArea.getVector();
        ::Vec svec = shearWall.getData();
        PetscInt nb = 0, nm = 0, ns = 0;
        VecGetLocalSize(bvec, &nb);
        VecGetLocalSize(mvec, &nm);
        VecGetLocalSize(svec, &ns);
        const PetscInt n = std::min(ns, std::min(nb, nm));
        const PetscScalar *barr = nullptr, *marr = nullptr;
        PetscScalar* sarr = nullptr;
        VecGetArrayRead(bvec, &barr);
        VecGetArrayRead(mvec, &marr);
        VecGetArray(svec, &sarr);
        for (PetscInt i = 0; i < n; ++i)
          sarr[i] = (std::abs(marr[i]) > 1.0e-30) ? (barr[i] / marr[i]) : PetscScalar(0);
        VecRestoreArray(svec, &sarr);
        VecRestoreArrayRead(mvec, &marr);
        VecRestoreArrayRead(bvec, &barr);
      }

      dState = dIter;
      solidAcceleration = dState;
      solidAcceleration -= dPred;
      solidAcceleration *= 1.0 / (betaN * dt * dt);
      solidVelocity = vPred;
      tmp = solidAcceleration;
      tmp *= gammaN * dt;
      solidVelocity += tmp;

      flux = BoundaryIntegral(Dot(uCur, normalFluid), qFlux).over(BoundaryFluid::Inlet);
      flux.assemble();
      qIn = flux(one);
      qOutSum = 0.0;
      for (const Attribute outlet : BoundaryFluid::Outlets)
      {
        flux = BoundaryIntegral(Dot(uCur, normalFluid), qFlux).over(outlet);
        flux.assemble();
        const Real qo = flux(one);
        qOut[outlet] = qo;
        qOutSum += qo;
      }

      // Interface power E = int (sigma_f n_f).(u_f - d_dot_s) on Gamma_FSI.
      auto sigmaFn = VectorFunction(dim, [&](const Point& pt) {
        Math::SpatialVector<Real> tr = wallSampleAt(pt).traction;
        for (Index i = 0; i < static_cast<Index>(dim); ++i)
          tr(i) = -tr(i);
        return tr;
      });
      flux = BoundaryIntegral(Dot(sigmaFn, uCur), qFlux).over(BoundaryFluid::FSI);
      flux.assemble();
      const Real ePowerFluid = flux(one);
      flux = BoundaryIntegral(Dot(sigmaFn, interfaceSolidVelocity), qFlux)
               .over(BoundaryFluid::FSI);
      flux.assemble();
      const Real ePowerSolid = flux(one);
      const Real eInterface = ePowerFluid - ePowerSolid;

      Real slipRms = 0.0;
      {
        const auto slipVec = uCur - interfaceSolidVelocity;
        flux = BoundaryIntegral(Dot(slipVec, slipVec), qFlux).over(BoundaryFluid::FSI);
        flux.assemble();
        const Real slipSq = flux(one);
        flux = BoundaryIntegral(RealFunction(1.0), qFlux).over(BoundaryFluid::FSI);
        flux.assemble();
        const Real fsiArea = flux(one);
        slipRms = (fsiArea > 0.0) ? std::sqrt(std::max(Real(0.0), slipSq / fsiArea)) : 0.0;
      }

      {
        PetscReal pMin = 0.0, pMax = 0.0;
        VecMin(pCur.getData(), PETSC_NULLPTR, &pMin);
        VecMax(pCur.getData(), PETSC_NULLPTR, &pMax);
        if (isRoot)
          Alert::Info() << "  [diag] p in [" << pMin << ", " << pMax << "] Pa"
                        << "  mass(qIn+qOut) = " << (qIn + qOutSum)
                        << "  E_iface = " << eInterface << " W"
                        << "  slip(RMS|u_f-u_s|) = " << slipRms << " m/s"
                        << "  | coupling: iters = " << couplesDone << "/"
                        << cfg.couplingIterations << "  interface change = " << lastRel
                        << Alert::Raise;
      }

      // Outlet compartments with the tissue pressure p_f^{n+1}; the lagged
      // fluxes handed to the next 0D step close the loop:
      //   Q_ar^ext = -q_in (inflow drawn from p_ar), Q_Phi^ext = sum_i q_v,i.
      for (const Attribute outlet : BoundaryFluid::Outlets)
        updateOutlet0D(cfg, wrms, s.pf, wk[outlet], qOut[outlet], dt);

      qArterialLag = -qIn;
      qTissueLag = 0.0;
      for (const auto& [tag, bc] : wk)
        qTissueLag += bc.qv;

      uOld.setData(uCur.getData());
      pOld.setData(pCur.getData());
      dOld.setData(dState.getData());
      solidVelocityOld.setData(solidVelocity.getData());
      solidAccelerationOld.setData(solidAcceleration.getData());
      vmsSubOld.setData(vmsSub.getSolution().getData());
      aleDispOld.setData(aleDisp.getData());

      xdmf_fluid.write(s.t).flush();
      xdmf_solid.write(s.t).flush();

      if (isRoot && csv)
      {
        const Real qLumped = modelInput.gammaAr * (s.par - s.pf);
        csv << s.t << ',' << s.y << ',' << yEpi << ',' << s.phi << ',' << s.pv << ','
            << s.par << ',' << s.pd << ',' << s.pf << ',' << s.lambdaBar << ',' << qIn
            << ',' << qOutSum << ',' << qTissueLag << ',' << qLumped;
        for (const Attribute outlet : BoundaryFluid::Outlets)
          csv << ',' << qOut[outlet] << ',' << wk[outlet].pout;
        {
          Real ptmSum = 0.0, phiASum = 0.0, phiVSum = 0.0, volSum = 0.0;
          for (const auto& [tag, bc] : wk)
          {
            ptmSum += bc.ptm;
            phiASum += bc.phiA;
            phiVSum += bc.phiV;
            volSum += bc.vol;
          }
          const Real nOut = static_cast<Real>(wk.size());
          csv << ',' << (ptmSum / nOut) << ',' << (phiASum / nOut) << ','
              << (phiVSum / nOut) << ',' << volSum;
        }
        csv << ',' << eInterface << ',' << lastRel << '\n';
        csv.flush();
      }
    }

    xdmf_fluid.close();
    xdmf_solid.close();
    if (isRoot && csv)
      csv.close();
  }
  catch (const std::exception& e)
  {
    std::cerr << "CoronaryArtery_Explicit_PoroElastic failed: " << e.what() << '\n';
    PetscFinalize();
    return 1;
  }

  PetscFinalize();
  return 0;
}
