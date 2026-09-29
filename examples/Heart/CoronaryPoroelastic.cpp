/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file CoronaryPoroelastic.cpp
 * @brief Rigid-wall coronary flow coupled in closed loop to the 0D
 * poroelastic left ventricle.
 *
 * Fluid-only counterpart of CoronaryArtery_Explicit_PoroElastic: the same
 * P1/P1 Carreau-Yasuda VMS Navier-Stokes on the lumen, without wall, ALE or
 * heart motion, coupled (staggered, one-step lag) to the healthy-calibrated
 * 0D poroelastic LV of CoronaryPoroelastic/HealthyLV.h.
 *
 *   0D -> 3D : inlet traction -p_ar n; tissue pressure p_f as reference and
 *              drainage pressure of every outlet compartment.
 *   3D -> 0D : the inlet flux is drawn from the proximal Windkessel node; the
 *              compartment drainage sum_i q_v,i enters the porosity balance.
 *   The 0D arterial conductance keeps the share (1 - lcaFraction) of the
 *   myocardium that the resolved tree does not perfuse.
 *
 * Fluid, step n -> n+1 on the fixed mesh (BDF1, Oseen linearization):
 *   rho/dt (u - u^n, v) + rho ((grad u) u^n, v) + rho/2 (div u^n u, v)
 *   + 2 (mu(u^n) eps(u), eps(v)) - (p, div v) + (div u, q)
 *   + VMS convective (tau_K) and grad-div (tau_C) orthogonal subscales
 *   + tau_p (grad p - Pi[grad p^n], grad q)                       (OSS-PSPG)
 *   + backflow stabilization, inlet traction, no-slip on the wall and the
 *     outlet traction -(p_c + R Q^n + R A (u.n - u^n.n)) n, R = R_a Phi_a.
 *
 * The outlet resistance is split implicit-explicit: the pointwise implicit
 * term R A (u.n) keeps the 3D-0D coupling stable for any dt, and the explicit
 * defect R (Q^n - A u^n.n) makes the mean traction exactly p_c + R Q once the
 * flux settles. The pointwise term alone is not a resistance R: the test
 * functions vanish on the no-slip rim of the outlet, so it enforces
 * p - p_c = R A u.n only in a rim-free weighted mean, and on coarse outlets
 * the resolved flux falls 15-20 % below (p - p_c)/R.
 *
 * Attributes: wall 2, inlet 4, outlets 7, 8, 9, 10, 14, 15.
 */

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numbers>

#include <petscsys.h>
#include <petscvec.h>

#include <Rodin/Alert.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/PETSc.h>
#include <Rodin/Solver.h>
#include <Rodin/Variational.h>

#include "CoronaryArtery/VMSConvectionIntegrator.h"
#include "CoronaryPoroelastic/CoronaryPoroelastic.h"
#include "CoronaryPoroelastic/HealthyLV.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Examples::Heart;
using namespace Rodin::Examples::Heart::CoronaryPoroelastic;

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  setPETScDefault("-ksp_type", "preonly");
  setPETScDefault("-pc_type", "lu");
  setPETScDefault("-pc_factor_mat_solver_type", "mumps");
  setPETScDefault("-coronary_mass_ksp_type", "cg");
  setPETScDefault("-coronary_mass_pc_type", "jacobi");

  try
  {
    Config cfg;
    readOptions(cfg);
    std::filesystem::create_directories(cfg.resultsDir);

    // ------------------------------------------------------------------
    // 0D poroelastic LV (healthy calibration).  The resolved tree takes the
    // share lcaFraction of the arterial conductance; the lagged 3D fluxes
    // enter through the two external sources.
    // ------------------------------------------------------------------
    Real qArterialLag = 0.0;
    Real qTissueLag = 0.0;

    HealthyTargets targets;
    HealthyCalibration calibration;
    Model::Input modelInput = makeHealthyInput(targets, &calibration);
    const Real gammaArFull = modelInput.gammaAr;
    modelInput.gammaAr = (1.0 - cfg.lcaFraction) * gammaArFull;
    modelInput.qArterialExternal = [&qArterialLag](Real) { return qArterialLag; };
    modelInput.qPerfusionExternal = [&qTissueLag](Real) { return qTissueLag; };

    Model model(modelInput);
    model.setMaxIterations(200)
      .setAbsoluteTolerance(1.0e-8)
      .setRelativeTolerance(1.0e-8)
      .setStepTolerance(1.0e-10)
      .setDampingFactor(1.0);
    model.initialize(makeHealthyInitialState(modelInput));

    Alert::Info() << "[0D] healthy LV: V_w0=" << calibration.wallVolume * 1e6
                  << " mL  R_p=" << calibration.Rp << "  R_d=" << calibration.Rd
                  << " Pa s/m^3  Murray r0=" << calibration.proximal.radius * 1e3
                  << " mm  L_p=" << calibration.proximal.length
                  << " m  L_d=" << calibration.distal.length
                  << " m  gamma_ar=" << gammaArFull
                  << "  gamma_ven=" << modelInput.gammaVen << Alert::Raise;

    const Real dt = cfg.dt;
    const size_t stepsPerCycle = static_cast<size_t>(std::llround(cfg.cyclePeriod / dt));

    // Warm-up: the resolved tree is the lumped conductance
    // lcaFraction gamma_ar(Phi); the last cycle gives the resting means.
    auto lumpedTree = [&](const Model::State& s) {
      const Real g = cfg.lcaFraction * gammaArFull *
        std::pow(s.phi / modelInput.phi0, modelInput.gammaArExponent);
      return g * (s.par - s.pf);
    };
    Real qLcaMean = 0.0, budgetMean = 0.0;
    {
      const size_t warmSteps = cfg.zeroDWarmupCycles * stepsPerCycle;
      for (size_t k = 0; k < warmSteps; ++k)
      {
        if (!model.step(dt).converged)
          throw std::runtime_error("0D warm-up failed at step " + std::to_string(k));
        const auto& s = model.getState();
        qArterialLag = lumpedTree(s);
        qTissueLag = qArterialLag;
        if (k + stepsPerCycle >= warmSteps)
        {
          qLcaMean += qArterialLag / static_cast<Real>(stepsPerCycle);
          budgetMean += (s.par - s.pf) / static_cast<Real>(stepsPerCycle);
        }
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

      Alert::Info() << "[0D] warm-up " << cfg.zeroDWarmupCycles << " cycle(s): p_ar="
                    << st.par << " Pa  p_f=" << st.pf << " Pa  phi=" << st.phi
                    << "  <Q_LCA>=" << qLcaMean * 6e7 << " mL/min  <p_ar - p_f>="
                    << budgetMean << " Pa" << Alert::Raise;
    }

    // ------------------------------------------------------------------
    // Mesh, spaces, fields.
    // ------------------------------------------------------------------
    MeshType mesh = loadMesh(cfg);
    const size_t dim = mesh.getSpaceDimension();

    using VelocityFES = H1<1, Math::SpatialVector<Real>, MeshType>;
    using PressureFES = H1<1, Real, MeshType>;

    VelocityFES uh(std::integral_constant<size_t, 1>{}, mesh, dim);
    PressureFES ph(std::integral_constant<size_t, 1>{}, mesh);

    PETSc::Variational::TrialFunction u(uh);
    PETSc::Variational::TrialFunction p(ph);
    PETSc::Variational::TestFunction v(uh);
    PETSc::Variational::TestFunction q(ph);

    PETSc::Variational::GridFunction uOld(uh);
    PETSc::Variational::GridFunction pOld(ph);
    PETSc::Variational::GridFunction uCur(uh);
    PETSc::Variational::GridFunction pCur(ph);
    PETSc::Variational::GridFunction one(ph);
    PETSc::Variational::GridFunction shearWall(uh);
    PETSc::Variational::GridFunction gradRec0(uh);
    PETSc::Variational::GridFunction gradRec1(uh);
    PETSc::Variational::GridFunction gradRec2(uh);

    auto zero = VectorFunction(dim, [&](const Point&) {
      Math::SpatialVector<Real> z(dim);
      z.setZero();
      return z;
    });

    uOld = zero;
    uCur = zero;
    shearWall = zero;
    gradRec0 = zero;
    gradRec1 = zero;
    gradRec2 = zero;
    pOld = 0.0;
    pCur = 0.0;
    one = 1.0;

    uOld.setName("FluidVelocity");
    pOld.setName("FluidPressure");
    shearWall.setName("shearStress");

    IO::XDMF xdmf(cfg.resultsDir + "/CoronaryPoroelastic");
    xdmf.setMesh(mesh);
    xdmf.add("FluidVelocity", uOld);
    xdmf.add("FluidPressure", pOld);
    xdmf.add("shearStress", shearWall);
    xdmf.write(0.0).flush();

    // ------------------------------------------------------------------
    // Outlet compartments: Murray split of the warm-up LCA flow.
    // ------------------------------------------------------------------
    const WRMSTable wrms = buildWRMSTable(cfg);
    std::map<Attribute, Outlet> wk;
    {
      PETSc::Variational::TestFunction qCal(ph);
      LinearForm<PressureFES, ::Vec> areaForm(qCal);
      std::map<Attribute, Real> r3;
      Real sumR3 = 0.0;
      for (const Attribute tag : Boundary::Outlets)
      {
        areaForm = BoundaryIntegral(one, qCal).over(tag);
        areaForm.assemble();
        const Real area = std::max<Real>(areaForm(one), 1e-12);
        const Real r = std::sqrt(area / std::numbers::pi_v<Real>);
        wk[tag].area = area;
        r3[tag] = r * r * r;
        sumR3 += r3[tag];
      }

      const Real dP = std::max<Real>(budgetMean, 1.0);
      for (const Attribute tag : Boundary::Outlets)
      {
        const Real w = r3[tag] / sumR3;
        auto& bc = wk[tag];
        calibrateOutlet(cfg, wrms, bc, qLcaMean * w, dP, cfg.coronaryComplianceTotal * w);
        bc.pc = bc.ptm + model.getState().pf;
        Alert::Info() << "  [outlet " << tag << "] A=" << bc.area
                      << " m^2  Q0=" << bc.q0 * 6e7 << " mL/min  Ra=" << bc.Ra
                      << "  Rv=" << bc.Rv << " Pa s/m^3  C=" << bc.C
                      << " m^3/Pa  C Rv=" << bc.C * bc.Rv << " s  ptm0=" << bc.ptm0
                      << " Pa  Phi_a0=" << bc.phiA << "  Phi_v0=" << bc.phiV
                      << Alert::Raise;
      }
    }

    Real pinValue = model.getState().par;
    std::array<Real, 6> poutValue{};
    std::array<Real, 6> zValue{};
    auto refreshOutletData = [&]() {
      for (size_t i = 0; i < Boundary::Outlets.size(); ++i)
      {
        const auto& bc = wk.at(Boundary::Outlets[i]);
        const Real R = cfg.outletResistanceScale * bc.Ra * bc.phiA;
        poutValue[i] = bc.pc + R * bc.q;
        zValue[i] = R * bc.area;
      }
    };
    refreshOutletData();

    auto pin = RealFunction([&](const Point&) { return pinValue; });
    auto pout0 = RealFunction([&](const Point&) { return poutValue[0]; });
    auto pout1 = RealFunction([&](const Point&) { return poutValue[1]; });
    auto pout2 = RealFunction([&](const Point&) { return poutValue[2]; });
    auto pout3 = RealFunction([&](const Point&) { return poutValue[3]; });
    auto pout4 = RealFunction([&](const Point&) { return poutValue[4]; });
    auto pout5 = RealFunction([&](const Point&) { return poutValue[5]; });
    auto zFn0 = RealFunction([&](const Point&) { return zValue[0]; });
    auto zFn1 = RealFunction([&](const Point&) { return zValue[1]; });
    auto zFn2 = RealFunction([&](const Point&) { return zValue[2]; });
    auto zFn3 = RealFunction([&](const Point&) { return zValue[3]; });
    auto zFn4 = RealFunction([&](const Point&) { return zValue[4]; });
    auto zFn5 = RealFunction([&](const Point&) { return zValue[5]; });

    // ------------------------------------------------------------------
    // Stabilization parameters, pointwise at the quadrature point:
    //   nu = mu_CY(gamma(u^n))/rho,  h_K = longest edge,
    //   tau_1 = 1/(4 nu/h_K^2 + 2|u^n|/h_K),
    //   tau_K = vmsScale/(rho/dt + rho/tau_1),  tau_C = rho h_K^2/(4 tau_1),
    //   tau_p = pspgScale tau_1/rho.
    // ------------------------------------------------------------------
    const auto& cy = cfg.viscosity;
    const Real rho = cfg.fluidDensity;
    const Real gammaReg = cy.gammaRegularization;
    const auto normal = BoundaryNormal(mesh);

    std::vector<Real> hCell;
    computeCellDiameters(mesh, hCell);

    // Lagged viscosity mu_CY(gamma(u^n)), cellwise constant for P1.
    std::vector<Real> muCell;
    computeCellViscosity(mesh, uh, uOld, cy, muCell);
    auto muAt = [&](const Point& pp) -> Real {
      const auto& poly = pp.getPolytope();
      if (poly.getDimension() == mesh.getDimension())
        return muCell[poly.getIndex()];
      const auto sym = 0.5 * (Jacobian(uOld) + Transpose(Jacobian(uOld)));
      return cy(std::sqrt(gammaReg * gammaReg + 2.0 * Dot(sym, sym).getValue(pp)));
    };
    RealFunction muLag = [&](const Point& pp) -> Real { return muAt(pp); };

    auto hAt = [&](const Point& pp) -> Real {
      const auto& poly = pp.getPolytope();
      return (poly.getDimension() == mesh.getDimension())
        ? hCell[poly.getIndex()]
        : std::pow(poly.getMeasure(), 1.0 / poly.getDimension());
    };
    auto tau1At = [&](const Point& pp) -> Real {
      const Real hK = std::max<Real>(hAt(pp), 1e-30);
      const Real nu = muAt(pp) / rho;
      const auto uc = uOld.getValue(pp);
      const Real speed = std::sqrt(Math::dot(uc, uc));
      return 1.0 / (4.0 * nu / (hK * hK) + 2.0 * speed / hK);
    };

    RealFunction vmsTauFn = [&](const Point& pp) -> Real {
      return cfg.vmsScale / (rho / dt + rho / tau1At(pp));
    };
    RealFunction sqrtTauCFn = [&](const Point& pp) -> Real {
      const Real hK = hAt(pp);
      return std::sqrt(cfg.gradDivScale * rho * hK * hK / (4.0 * tau1At(pp)));
    };
    RealFunction tauCFn = [&](const Point& pp) -> Real {
      const Real s = sqrtTauCFn(pp);
      return s * s;
    };
    RealFunction tauPFn = [&](const Point& pp) -> Real {
      return cfg.pspgScale * tau1At(pp) / rho;
    };
    const Real orthogonalPSPG = cfg.pspgOrthogonal ? 1.0 : 0.0;

    // Orthogonal-subscale projections (lagged, mass matrix assembled once):
    // Pi[(grad u^n) u^n], dynamic subscale u', Pi[sqrt(tau_C) div u^n],
    // Pi[grad p^n].
    PETSc::Variational::TrialFunction projTrial(uh);
    PETSc::Variational::TestFunction projTest(uh);
    PETSc::Variational::TrialFunction projScalarTrial(ph);
    PETSc::Variational::TestFunction projScalarTest(ph);
    MassProjection projectVector(projTrial, projTest);
    MassProjection projectScalar(projScalarTrial, projScalarTest);

    PETSc::Variational::GridFunction convProj(uh);
    PETSc::Variational::GridFunction vmsSub(uh);
    PETSc::Variational::GridFunction vmsSubOld(uh);
    PETSc::Variational::GridFunction piTilde(ph);
    PETSc::Variational::GridFunction gradPProj(uh);
    convProj = zero;
    vmsSub = zero;
    vmsSubOld = zero;
    piTilde = 0.0;
    gradPProj = zero;

    const auto convTarget = Mult(Jacobian(uOld), uOld);
    auto vmsSubUpdate =
      VectorFunction(dim, [&](const Point& pp) -> Math::SpatialVector<Real> {
        const auto conv = convTarget.getValue(pp);
        const auto proj = convProj.getValue(pp);
        const auto old = vmsSubOld.getValue(pp);
        const Real tau = vmsTauFn(pp);
        Math::SpatialVector<Real> out(dim);
        for (Index c = 0; c < static_cast<Index>(dim); ++c)
          out(c) = tau * rho * (old(c) / dt - (conv(c) - proj(c)));
        return out;
      });

    LinearForm<VelocityFES, ::Vec> convRhs(projTest);
    convRhs = Integral(convTarget, projTest);
    LinearForm<VelocityFES, ::Vec> subRhs(projTest);
    subRhs = Integral(vmsSubUpdate, projTest);
    LinearForm<PressureFES, ::Vec> piTildeRhs(projScalarTest);
    piTildeRhs = Integral(sqrtTauCFn * Div(uOld), projScalarTest);
    LinearForm<VelocityFES, ::Vec> gradPRhs(projTest);
    gradPRhs = Integral(Grad(pOld), projTest);

    // ------------------------------------------------------------------
    // Flow problem.
    // ------------------------------------------------------------------
    const auto symU = 0.5 * (Jacobian(u) + Transpose(Jacobian(u)));
    const auto symV = 0.5 * (Jacobian(v) + Transpose(Jacobian(v)));
    const auto convU = Mult(Jacobian(u), uOld);
    const auto backflowOut = 0.5 * cfg.backflowStabilization * rho *
      Max(-Dot(uOld, normal), 0.0);
    const auto backflowIn = 0.5 * cfg.backflowStabilization * rho *
      Max(Dot(uOld, normal), 0.0);
    const auto uTangential = u - Dot(u, normal) * normal;
    const auto& O = Boundary::Outlets;

    Problem flow(u, p, v, q);
    flow = (rho / dt) * Integral(u, v) - (rho / dt) * Integral(uOld, v) +
      rho * Integral(Dot(convU, v)) + 0.5 * rho * Integral(Div(uOld) * Dot(u, v)) +
      VMSConvectionBilinearIntegrator(u, v, uOld, vmsTauFn, rho) -
      VMSConvectionLinearIntegrator(v, vmsSub, uOld, convProj, vmsTauFn, rho, dt) +
      VMSGradDivBilinearIntegrator(u, v, tauCFn) -
      VMSGradDivLinearIntegrator(v, piTilde, sqrtTauCFn) +
      2.0 * Integral(muLag * symU, symV) - Integral(p, Div(v)) + Integral(Div(u), q) +
      Integral(tauPFn * Grad(p), Grad(q)) -
      orthogonalPSPG * Integral(tauPFn * gradPProj, Grad(q)) +
      BoundaryIntegral(backflowIn * Dot(u, v)).over(Boundary::Inlet) +
      BoundaryIntegral(backflowOut * Dot(u, v)).over(O[0], O[1], O[2], O[3], O[4], O[5]) +
      BoundaryIntegral(pin * Dot(v, normal)).over(Boundary::Inlet) +
      BoundaryIntegral(pout0 * Dot(v, normal)).over(O[0]) +
      BoundaryIntegral(pout1 * Dot(v, normal)).over(O[1]) +
      BoundaryIntegral(pout2 * Dot(v, normal)).over(O[2]) +
      BoundaryIntegral(pout3 * Dot(v, normal)).over(O[3]) +
      BoundaryIntegral(pout4 * Dot(v, normal)).over(O[4]) +
      BoundaryIntegral(pout5 * Dot(v, normal)).over(O[5]) +
      BoundaryIntegral(zFn0 * Dot(Dot(u, normal) * normal, v)).over(O[0]) -
      BoundaryIntegral(zFn0 * Dot(Dot(uOld, normal) * normal, v)).over(O[0]) +
      BoundaryIntegral(zFn1 * Dot(Dot(u, normal) * normal, v)).over(O[1]) -
      BoundaryIntegral(zFn1 * Dot(Dot(uOld, normal) * normal, v)).over(O[1]) +
      BoundaryIntegral(zFn2 * Dot(Dot(u, normal) * normal, v)).over(O[2]) -
      BoundaryIntegral(zFn2 * Dot(Dot(uOld, normal) * normal, v)).over(O[2]) +
      BoundaryIntegral(zFn3 * Dot(Dot(u, normal) * normal, v)).over(O[3]) -
      BoundaryIntegral(zFn3 * Dot(Dot(uOld, normal) * normal, v)).over(O[3]) +
      BoundaryIntegral(zFn4 * Dot(Dot(u, normal) * normal, v)).over(O[4]) -
      BoundaryIntegral(zFn4 * Dot(Dot(uOld, normal) * normal, v)).over(O[4]) +
      BoundaryIntegral(zFn5 * Dot(Dot(u, normal) * normal, v)).over(O[5]) -
      BoundaryIntegral(zFn5 * Dot(Dot(uOld, normal) * normal, v)).over(O[5]) +
      cfg.inletImpedance *
        BoundaryIntegral(Dot(Dot(u, normal) * normal, v)).over(Boundary::Inlet) +
      cfg.inletTangentialDamping *
        BoundaryIntegral(Dot(uTangential, v)).over(Boundary::Inlet) +
      DirichletBC(u, zero).on(Boundary::Wall);

    // ------------------------------------------------------------------
    // Pointwise WSS from the L2-recovered velocity gradient:
    //   D = (G + G^T)/2, gamma = sqrt(2 D:D), tau_w = (I - n n) 2 mu D n.
    // ------------------------------------------------------------------
    const auto jacRow = [&](int i) {
      return VectorFunction(Component(Jacobian(uCur), i, 0),
        Component(Jacobian(uCur), i, 1), Component(Jacobian(uCur), i, 2));
    };
    const auto row0 = jacRow(0);
    const auto row1 = jacRow(1);
    const auto row2 = jacRow(2);
    LinearForm<VelocityFES, ::Vec> row0Rhs(projTest);
    row0Rhs = Integral(row0, projTest);
    LinearForm<VelocityFES, ::Vec> row1Rhs(projTest);
    row1Rhs = Integral(row1, projTest);
    LinearForm<VelocityFES, ::Vec> row2Rhs(projTest);
    row2Rhs = Integral(row2, projTest);

    auto wallShear = VectorFunction(dim, [&](const Point& pt) {
      const Math::SpatialVector<Real> n = normal.getValue(pt);
      const auto g0 = gradRec0(pt);
      const auto g1 = gradRec1(pt);
      const auto g2 = gradRec2(pt);
      Math::SpatialMatrix<Real> G(3, 3);
      for (std::uint8_t j = 0; j < 3; ++j)
      {
        G(0, j) = g0(j);
        G(1, j) = g1(j);
        G(2, j) = g2(j);
      }
      Real dd = 0.0;
      Math::SpatialMatrix<Real> D(3, 3);
      for (std::uint8_t i = 0; i < 3; ++i)
        for (std::uint8_t j = 0; j < 3; ++j)
        {
          D(i, j) = 0.5 * (G(i, j) + G(j, i));
          dd += D(i, j) * D(i, j);
        }
      const Real mu = cy(std::sqrt(gammaReg * gammaReg + 2.0 * dd));
      Math::SpatialVector<Real> t(3);
      Real tn = 0.0;
      for (std::uint8_t i = 0; i < 3; ++i)
      {
        Real Dn = 0.0;
        for (std::uint8_t j = 0; j < 3; ++j)
          Dn += D(i, j) * n(j);
        t(i) = 2.0 * mu * Dn;
        tn += t(i) * n(i);
      }
      for (std::uint8_t i = 0; i < 3; ++i)
        t(i) -= tn * n(i);
      return t;
    });

    auto computeWSS = [&]() {
      projectVector(row0Rhs, gradRec0);
      projectVector(row1Rhs, gradRec1);
      projectVector(row2Rhs, gradRec2);

      // Area-weighted (lumped) nodal average on the wall.
      PETSc::Variational::TestFunction wssTest(uh);
      const auto ones = VectorFunction(dim, [&](const Point&) {
        Math::SpatialVector<Real> o(dim);
        for (Index c = 0; c < static_cast<Index>(dim); ++c)
          o(c) = 1.0;
        return o;
      });
      LinearForm<VelocityFES, ::Vec> load(wssTest);
      load = BoundaryIntegral(wallShear, wssTest).over(Boundary::Wall);
      load.assemble();
      LinearForm<VelocityFES, ::Vec> area(wssTest);
      area = BoundaryIntegral(ones, wssTest).over(Boundary::Wall);
      area.assemble();

      PetscInt n = 0;
      VecGetLocalSize(shearWall.getData(), &n);
      const PetscScalar *b = nullptr, *m = nullptr;
      PetscScalar* s = nullptr;
      VecGetArrayRead(load.getVector(), &b);
      VecGetArrayRead(area.getVector(), &m);
      VecGetArray(shearWall.getData(), &s);
      for (PetscInt i = 0; i < n; ++i)
        s[i] = (std::abs(m[i]) > 1.0e-30) ? b[i] / m[i] : PetscScalar(0);
      VecRestoreArray(shearWall.getData(), &s);
      VecRestoreArrayRead(area.getVector(), &m);
      VecRestoreArrayRead(load.getVector(), &b);
    };

    PETSc::Variational::TestFunction qFlux(ph);
    LinearForm<PressureFES, ::Vec> flux(qFlux);
    auto boundaryFlux = [&](Attribute tag) {
      flux = BoundaryIntegral(Dot(uCur, normal), qFlux).over(tag);
      flux.assemble();
      return flux(one);
    };

    std::ofstream csv(cfg.resultsDir + "/CoronaryPoroelastic.csv");
    csv << std::setprecision(12);
    csv << "t,lv_y,lv_phi,lv_pv,lv_par,lv_pd,lv_pf,lv_V,q_in,q_out_total,q_tissue"
        << ",q_lumped,q_ven";
    for (const Attribute tag : Boundary::Outlets)
      csv << ",q_out_" << tag << ",p_c_" << tag << ",ptm_" << tag;
    csv << ",phi_a_mean,phi_v_mean,stored_volume,p_min,p_max\n";

    // Hand-off: the 0D keeps the lumped-tree fluxes until the first 3D step.
    qTissueLag = 0.0;
    for (const auto& [tag, bc] : wk)
      qTissueLag += bc.qv;

    // ==================================================================
    // Time loop, step n -> n+1:
    //   (1) 0D advanced with the lagged 3D fluxes;
    //   (2) p_in = p_ar^{n+1}, p_out,i = p_tm,i^n + p_f^{n+1} + R_i Q_i^n;
    //   (3) lagged projections, one Oseen solve;
    //   (4) fluxes, outlet compartments with p_f^{n+1}, lagged fluxes.
    // ==================================================================
    for (size_t step = 1; step <= cfg.nsteps; ++step)
    {
      if (!model.step(dt).converged)
      {
        std::cerr << "0D model failed at step " << step << '\n';
        break;
      }
      const auto& s = model.getState();

      computeCellViscosity(mesh, uh, uOld, cy, muCell);
      pinValue = s.par;
      for (auto& [tag, bc] : wk)
        bc.pc = bc.ptm + s.pf;
      refreshOutletData();

      projectVector(convRhs, convProj);
      projectVector(subRhs, vmsSub);
      projectScalar(piTildeRhs, piTilde);
      if (cfg.pspgOrthogonal)
        projectVector(gradPRhs, gradPProj);

      flow.assemble().setFieldSplits();
      Solver::KSP(flow).solve();
      uCur.setData(u.getSolution().getData());
      pCur.setData(p.getSolution().getData());

      const Real qIn = boundaryFlux(Boundary::Inlet);
      std::map<Attribute, Real> qOut;
      Real qOutSum = 0.0;
      for (const Attribute tag : Boundary::Outlets)
      {
        qOut[tag] = boundaryFlux(tag);
        qOutSum += qOut[tag];
      }

      for (const Attribute tag : Boundary::Outlets)
        updateOutlet(cfg, wrms, s.pf, wk[tag], qOut[tag], dt);

      qArterialLag = -qIn;
      qTissueLag = 0.0;
      for (const auto& [tag, bc] : wk)
        qTissueLag += bc.qv;

      PetscReal pMin = 0.0, pMax = 0.0;
      VecMin(pCur.getData(), PETSC_NULLPTR, &pMin);
      VecMax(pCur.getData(), PETSC_NULLPTR, &pMax);

      Alert::Info() << "step " << step << "/" << cfg.nsteps << "  t=" << s.t
                    << "  p_ar=" << s.par << "  p_f=" << s.pf << " Pa  Q_in="
                    << -qIn * 6e7 << "  Q_out=" << qOutSum * 6e7
                    << "  q_v=" << qTissueLag * 6e7 << " mL/min  p in [" << pMin << ", "
                    << pMax << "]" << Alert::Raise;

      uOld.setData(uCur.getData());
      pOld.setData(pCur.getData());
      vmsSubOld.setData(vmsSub.getData());

      if (step % cfg.outputInterval == 0)
      {
        computeWSS();
        xdmf.write(s.t).flush();
      }

      const Real R = modelInput.R0 + s.y;
      const Real V = 4.0 / 3.0 * std::numbers::pi_v<Real> * R * R * R;
      csv << s.t << ',' << s.y << ',' << s.phi << ',' << s.pv << ',' << s.par << ','
          << s.pd << ',' << s.pf << ',' << V << ',' << qIn << ',' << qOutSum << ','
          << qTissueLag << ',' << s.qPerfusionIn << ',' << s.qPerfusionOut;
      Real phiA = 0.0, phiV = 0.0, vol = 0.0;
      for (const Attribute tag : Boundary::Outlets)
      {
        const auto& bc = wk[tag];
        csv << ',' << qOut[tag] << ',' << bc.pc << ',' << bc.ptm;
        phiA += bc.phiA / static_cast<Real>(wk.size());
        phiV += bc.phiV / static_cast<Real>(wk.size());
        vol += bc.vol;
      }
      csv << ',' << phiA << ',' << phiV << ',' << vol << ',' << pMin << ',' << pMax
          << '\n';
      csv.flush();
    }

    xdmf.close();
  }
  catch (const std::exception& e)
  {
    std::cerr << "CoronaryPoroelastic failed: " << e.what() << '\n';
    PetscFinalize();
    return 1;
  }

  PetscFinalize();
  return 0;
}
