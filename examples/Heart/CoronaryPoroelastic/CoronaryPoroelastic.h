/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef EXAMPLES_HEART_CORONARYPOROELASTIC_CORONARYPOROELASTIC_H
#define EXAMPLES_HEART_CORONARYPOROELASTIC_CORONARYPOROELASTIC_H

/**
 * @file CoronaryPoroelastic.h
 * @brief Configuration, outlet compartments and helpers of the rigid-wall
 * coronary flow coupled in closed loop to the 0D poroelastic left ventricle.
 *
 * Outlet i is an intramyocardial compartment embedded in the tissue at the
 * interstitial pressure @f$ p_f @f$ of the 0D model:
 * @f[
 *   \text{3D} \xrightarrow{R_a \Phi_a\ \text{(implicit)}} (p_c, C)
 *   \xrightarrow{R_v(p_{tm}) \Phi_v} p_f, \qquad
 *   C \dot p_{tm} = Q - q_v, \qquad p_{tm} = p_c - p_f .
 * @f]
 * The venular limb is a collapsible bed of fixed vessel count and length,
 * @f$ R_v \propto V^{-\kappa_v} @f$ with @f$ V = C p_{tm} @f$, so
 * @f[
 *   q_v = \frac{p_{tm}}{R_v \Phi_v}\Big(\frac{p_{tm}}{p_{tm,0}}\Big)^{\kappa_v},
 *   \quad p_{tm} > 0, \qquad q_v = 0, \quad p_{tm} \le 0,
 * @f]
 * (@f$ \kappa_v = 2 @f$: Poiseuille, @f$ \kappa_v = 0 @f$: rigid Starling
 * throat). The vessels close continuously when the tissue pressure reaches
 * the compartment pressure in systole. @f$ \Phi_{a,v} = \mu_{ap}/\mu_N @f$ are
 * the WRMS apparent-viscosity factors of the Carreau-Yasuda blood at the
 * limb's nominal shear rate.
 */

#include <algorithm>
#include <cassert>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <petscksp.h>
#include <petscsys.h>

#include <Eigen/Dense>

#include <Rodin/Geometry.h>
#include <Rodin/IO/MEDIT.h>
#include <Rodin/Types.h>
#include <Rodin/PETSc.h>
#include <Rodin/Variational.h>

namespace Rodin::Examples::Heart::CoronaryPoroelastic
{
  using MeshType = Geometry::Mesh<Context::Local>;
  using Geometry::Attribute;

  /// Boundary attributes of the coronary lumen mesh.
  /// Boundary attributes of LCA_volumen_v22_smooth.mesh (left main ostium,
  /// four distal outlets), as in CoupledLV0DCoronary3D.
  struct Boundary
  {
      static constexpr Attribute Wall = 2;
      static constexpr Attribute Inlet = 40;
      static constexpr std::array<Attribute, 4> Outlets{{36, 37, 38, 39}};
  };

  /// Carreau-Yasuda blood viscosity.
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

      Real derivative(Real g) const
      {
        const Real base = 1.0 + std::pow(lambda * g, yasuda);
        return (mu0 - muInf) * (n - 1.0) * std::pow(base, (n - 1.0 - yasuda) / yasuda) *
          std::pow(lambda, yasuda) * std::pow(g, yasuda - 1.0);
      }
  };

  /// One intramyocardial compartment per outlet.
  struct Outlet
  {
      Real ptm = 0.0;   ///< Transmural pressure p_c - p_f (the state).
      Real ptm0 = 0.0;  ///< Resting transmural pressure (reference of R_v(V)).
      Real pc = 0.0;    ///< Compartment pressure.
      Real qv = 0.0;    ///< Venular drainage into the tissue.
      Real q = 0.0;     ///< Last 3D outlet flux Q^n.
      Real Ra = 0.0;    ///< Arteriolar resistance at mu_N.
      Real Rv = 0.0;    ///< Venular resistance at mu_N and p_tm0.
      Real C = 0.0;     ///< Compliance.
      Real area = 0.0;  ///< Outlet area.
      Real q0 = 0.0;    ///< Calibrated resting flow.
      Real gammaA = 0.0; ///< Resting wall shear rate, arteriolar limb.
      Real gammaV = 0.0; ///< Resting wall shear rate, venular limb.
      Real phiA = 1.0;  ///< mu_ap / mu_N, arteriolar limb.
      Real phiV = 1.0;  ///< mu_ap / mu_N, venular limb.
      Real vol = 0.0;   ///< Stored volume C max(p_tm, 0).
  };

  /// Universal WRMS apparent viscosity mu_ap(gammadot_nom), log-log table.
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
      std::string meshPath =
        "../resources/examples/Heart/CoronaryArtery/LCA_volumen_v22_smooth.mesh";
      std::string resultsDir = "results_coronary_poroelastic";
      Real meshScale = 1.0e-3;

      Real cyclePeriod = 0.85;
      Real dt = 1.0e-3;
      size_t nsteps = 3 * static_cast<size_t>(0.85 / 1.0e-3);
      size_t outputInterval = 10; ///< XDMF + WSS every k steps (CSV every step).

      /// 0D cycles alone before the hand-off, the resolved tree lumped into
      /// the fraction lcaFraction of the 0D arterial conductance.
      size_t zeroDWarmupCycles = 4;
      /// Share of the LV perfusion carried by the resolved tree.
      Real lcaFraction = 0.65;

      // Outlets: Murray split Q_i ~ r_i^3 of the warm-up mean LCA flow on the
      // mean budget dP = <p_ar - p_f>: R_v = f_v dP/Q_i, R_a = (1 - f_v) dP/Q_i.
      Real newtonianCalibrationViscosity = 0.0035;
      Real coronaryComplianceTotal = 4.0e-10; // m^3/Pa
      Real venularPressureFraction = 0.13;
      Real venularCollapseExponent = 2.0; ///< kappa_v of R_v(V) ~ V^{-kappa_v}.

      // Morphometric operating point (r, v) of each limb: gamma_0 = 4v/r.
      Real arteriolarRadius = 25.0e-6;
      Real arteriolarVelocity = 5.0e-3;
      Real venularRadius = 30.0e-6;
      Real venularVelocity = 3.0e-3;

      Real fluidDensity = 1060.0;
      CarreauYasuda viscosity;
      Real backflowStabilization = 1.0;
      Real inletImpedance = 1.0e3;
      Real inletTangentialDamping = 1.0e3;
      Real outletResistanceScale = 1.0;

      Real vmsScale = 1.0;
      Real gradDivScale = 1.0;
      Real pspgScale = 1.0;
      bool pspgOrthogonal = true; ///< tau_p (grad p - Pi[grad p^n], grad q).

      // WRMS table and outlet Newton.
      Real tableTauMin = 1.0e-6;
      Real tableTauMax = 1.0e4;
      int tableNodes = 241;
      int integralSteps = 2000;
      Real outletStepTolerance = 1.0e-9;
      int outletMaxIterations = 50;
      Real zeroFlowTolerance = 1.0e-16;
  };

  inline void setPETScDefault(const char* key, const char* value)
  {
    PetscBool set = PETSC_FALSE;
    PetscOptionsHasName(PETSC_NULLPTR, PETSC_NULLPTR, key, &set);
    if (!set)
      PetscOptionsSetValue(PETSC_NULLPTR, key, value);
  }

  inline void readOptions(Config& cfg)
  {
    const Real inf = std::numeric_limits<Real>::infinity();
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
    auto optString = [](const char* key, std::string& value) {
      char buf[PETSC_MAX_PATH_LEN];
      PetscBool set = PETSC_FALSE;
      PetscOptionsGetString(PETSC_NULLPTR, PETSC_NULLPTR, key, buf, sizeof(buf), &set);
      if (set)
        value = buf;
    };

    optString("-coronary_mesh", cfg.meshPath);
    optString("-coronary_results", cfg.resultsDir);
    optReal("-coronary_dt", cfg.dt, 1.0e-9, inf);
    optInt("-coronary_nsteps", cfg.nsteps, 0);
    optInt("-coronary_output_interval", cfg.outputInterval, 1);
    optInt("-coronary_warmup_cycles", cfg.zeroDWarmupCycles, 1);
    optReal("-coronary_lca_fraction", cfg.lcaFraction, 0.0, 1.0);
    optReal("-coronary_compliance", cfg.coronaryComplianceTotal, 0.0, inf);
    optReal("-coronary_venular_fraction", cfg.venularPressureFraction, 0.01, 0.99);
    optReal("-coronary_collapse_exponent", cfg.venularCollapseExponent, 0.0, 8.0);
    optReal("-coronary_vms_scale", cfg.vmsScale, 0.0, inf);
    optReal("-coronary_graddiv_scale", cfg.gradDivScale, 0.0, inf);
    optReal("-coronary_pspg_scale", cfg.pspgScale, 0.0, inf);
    optBool("-coronary_pspg_orthogonal", cfg.pspgOrthogonal);
    optReal("-coronary_inlet_impedance", cfg.inletImpedance, 0.0, inf);
    optReal("-coronary_outlet_resistance_scale", cfg.outletResistanceScale, 0.0, inf);
  }

  /**
   * @brief WRMS closure of a tube carrying Carreau-Yasuda blood:
   * @f$ Q = \pi R^3 I(\tau_w)/\tau_w^3 @f$, @f$ I = \int_0^{\tau_w}\tau^2\dot\gamma\,d\tau @f$,
   * @f$ \mu_{ap} = \tau_w^4/(4I) @f$ at @f$ \dot\gamma_{nom} = \tau_w/\mu_{ap} @f$;
   * independent of @f$ R, L @f$, tabulated once.
   */
  inline WRMSTable buildWRMSTable(const Config& cfg)
  {
    const auto& cy = cfg.viscosity;
    auto shearAt = [&](Real tauW) -> Real {
      Real lo = std::log(1e-14);
      Real hi = std::log(std::max<Real>(tauW / cy.muInf, 1e-12)) + 1.0;
      Real s = 0.5 * (lo + hi);
      for (int it = 0; it < 200; ++it)
      {
        const Real g = std::exp(s);
        const Real f = cy(g) * g - tauW;
        (f < 0.0 ? lo : hi) = s;
        const Real df = (cy(g) + g * cy.derivative(g)) * g;
        Real sNext = (std::abs(df) > 0.0) ? s - f / df : 0.5 * (lo + hi);
        if (!(sNext > lo && sNext < hi))
          sNext = 0.5 * (lo + hi);
        const Real step = std::abs(sNext - s);
        s = sNext;
        if (step < 1e-12)
          break;
      }
      return std::exp(s);
    };

    // I = int_0^{g_w} g^3 mu^2 (mu + g mu') dg (tau = mu g), Simpson.
    auto integral = [&](Real gW) -> Real {
      const int m = 2 * (cfg.integralSteps / 2);
      const Real h = gW / static_cast<Real>(m);
      auto f = [&](Real g) -> Real {
        if (g <= 0.0)
          return 0.0;
        const Real mg = cy(g);
        return g * g * g * mg * mg * (mg + g * cy.derivative(g));
      };
      Real sum = f(0.0) + f(gW);
      for (int i = 1; i < m; ++i)
        sum += ((i % 2 == 1) ? 4.0 : 2.0) * f(static_cast<Real>(i) * h);
      return sum * h / 3.0;
    };

    WRMSTable table;
    const int nodes = std::max(2, cfg.tableNodes);
    const Real a = std::log(cfg.tableTauMin), b = std::log(cfg.tableTauMax);
    for (int i = 0; i < nodes; ++i)
    {
      const Real tauW = std::exp(a + (b - a) * i / static_cast<Real>(nodes - 1));
      const Real I = integral(shearAt(tauW));
      if (!(I > 0.0) || !std::isfinite(I))
        continue;
      const Real muAp = std::pow(tauW, 4.0) / (4.0 * I);
      const Real gNom = tauW / muAp;
      if (!(muAp > 0.0) || !std::isfinite(muAp) ||
        (!table.logGamma.empty() && std::log(gNom) <= table.logGamma.back()))
        continue;
      table.logGamma.push_back(std::log(gNom));
      table.logMu.push_back(std::log(muAp));
    }
    if (table.logGamma.size() < 2)
    {
      table.logGamma = {std::log(1e-6), std::log(1e6)};
      table.logMu = {std::log(cy.mu0), std::log(cy.mu0)};
    }
    return table;
  }

  /**
   * @brief Calibrates one outlet at rest: flow Q, budget dP = <p_ar - p_f>.
   *
   * The effective resistances @f$ R_a\Phi_{a,0} @f$, @f$ R_v\Phi_{v,0} @f$
   * (not the Newtonian ones) carry the budget at the resting shear rates:
   * @f$ R_a\Phi_{a,0} = (1-f_v)\,dP/Q @f$, @f$ R_v\Phi_{v,0} = f_v\,dP/Q @f$,
   * so the resting state @f$ p_{tm,0} = f_v\,dP @f$ drains exactly @f$ Q @f$
   * and the arteriolar drop is @f$ (1-f_v)\,dP @f$.
   */
  inline void calibrateOutlet(
    const Config& cfg, const WRMSTable& wrms, Outlet& bc, Real Q, Real dP, Real C)
  {
    const Real fv = cfg.venularPressureFraction;
    const Real muN = cfg.newtonianCalibrationViscosity;
    bc.q0 = Q;
    bc.C = C;
    bc.gammaA = 4.0 * cfg.arteriolarVelocity / cfg.arteriolarRadius;
    bc.gammaV = 4.0 * cfg.venularVelocity / cfg.venularRadius;
    bc.phiA = wrms(bc.gammaA) / muN;
    bc.phiV = wrms(bc.gammaV) / muN;
    bc.Ra = (1.0 - fv) * dP / (Q * bc.phiA);
    bc.Rv = fv * dP / (Q * bc.phiV);
    bc.ptm0 = fv * dP;
    bc.ptm = bc.ptm0;
    bc.qv = Q;
    bc.vol = C * bc.ptm0;
  }

  /**
   * @brief Advances one outlet compartment by implicit Euler (scalar Newton
   * on @f$ p_{tm} @f$, @f$ \Phi_v @f$ lagged within the iteration).
   *
   * @param[in] pf Tissue pressure at @f$ t_{n+1} @f$.
   * @param[in] Q 3D outlet flux at @f$ t_{n+1} @f$.
   */
  inline void updateOutlet(
    const Config& cfg, const WRMSTable& wrms, Real pf, Outlet& bc, Real Q, Real dt)
  {
    const Real kv = cfg.venularCollapseExponent;
    const Real muN = cfg.newtonianCalibrationViscosity;
    const Real ptmOld = bc.ptm;
    const Real ptm0 = std::max<Real>(bc.ptm0, 1e-300);
    const Real q0 = std::max<Real>(std::abs(bc.q0), 1e-300);

    auto phiV = [&](Real q) {
      const Real aq = std::max<Real>(std::abs(q), cfg.zeroFlowTolerance);
      return wrms(bc.gammaV * aq / q0) / muN;
    };
    // Drainage and its derivative at fixed Phi_v.
    auto drain = [&](Real ptm, Real pv) -> std::pair<Real, Real> {
      if (ptm <= 0.0)
        return {0.0, 0.0};
      const Real G = std::pow(ptm / ptm0, kv) / (bc.Rv * pv);
      return {ptm * G, (1.0 + kv) * G};
    };

    Real ptm = ptmOld;
    Real qv = bc.qv;
    bool converged = false;
    for (int it = 0; it < cfg.outletMaxIterations; ++it)
    {
      const Real pv = phiV(qv);
      const auto [q, dq] = drain(ptm, pv);
      qv = q;
      const Real R = bc.C * (ptm - ptmOld) / dt - Q + q;
      const Real J = bc.C / dt + dq;
      const Real d = -R / J;
      ptm += d;
      if (std::abs(d) < cfg.outletStepTolerance * (1.0 + std::abs(ptm)))
      {
        converged = true;
        break;
      }
    }
    if (!converged || !std::isfinite(ptm))
    {
      std::cerr << "Warning: outlet compartment Newton did not converge; keeping "
                << "the previous state.\n";
      ptm = ptmOld;
    }

    bc.phiV = phiV(qv);
    bc.qv = drain(ptm, bc.phiV).first;
    bc.ptm = ptm;
    bc.pc = ptm + pf;
    bc.q = Q;
    bc.vol = bc.C * std::max<Real>(ptm, 0.0);
    const Real aQ = std::max<Real>(std::abs(Q), cfg.zeroFlowTolerance);
    bc.phiA = wrms(bc.gammaA * aQ / q0) / muN;
  }

  /**
   * @brief L2 projection with a mass matrix assembled once (fixed mesh):
   * each call assembles only the right-hand side and solves M x = b with the
   * "coronary_mass_" KSP.
   */
  class MassProjection
  {
    public:
      template <class Trial, class Test>
      MassProjection(Trial& trial, Test& test)
      {
        Variational::Problem mass(trial, test);
        mass = Variational::Integral(trial, test);
        mass.assemble();
        PetscErrorCode ierr = MatDuplicate(
          mass.getLinearSystem().getOperator(), MAT_COPY_VALUES, &m_M);
        assert(ierr == PETSC_SUCCESS);
        ierr = KSPCreate(PETSC_COMM_SELF, &m_ksp);
        assert(ierr == PETSC_SUCCESS);
        ierr = KSPSetOperators(m_ksp, m_M, m_M);
        assert(ierr == PETSC_SUCCESS);
        ierr = KSPSetOptionsPrefix(m_ksp, "coronary_mass_");
        assert(ierr == PETSC_SUCCESS);
        ierr = KSPSetFromOptions(m_ksp);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
      }

      MassProjection(const MassProjection&) = delete;
      MassProjection& operator=(const MassProjection&) = delete;

      ~MassProjection()
      {
        KSPDestroy(&m_ksp);
        MatDestroy(&m_M);
      }

      /// Assembles @p rhs and writes the projection into @p target.
      template <class LinearForm, class GridFunction>
      void operator()(LinearForm& rhs, GridFunction& target)
      {
        rhs.assemble();
        PetscErrorCode ierr = KSPSolve(m_ksp, rhs.getVector(), target.getData());
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
      }

    private:
      ::Mat m_M = PETSC_NULLPTR;
      ::KSP m_ksp = PETSC_NULLPTR;
  };

  /// Longest edge of every cell (stabilization length h_K).
  inline void computeCellDiameters(const MeshType& mesh, std::vector<Real>& h)
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
          hmax = std::max<Real>(hmax,
            (mesh.getVertexCoordinates(verts[a]) - mesh.getVertexCoordinates(verts[b]))
              .norm());
      h[it->getIndex()] = hmax;
    }
  }

  /**
   * @brief Lagged Carreau-Yasuda viscosity of every cell.
   *
   * For P1 velocity the gradient is constant on each tetrahedron,
   * @f$ \nabla u = [u_i - u_0]\,[x_i - x_0]^{-1} @f$, so this is exactly the
   * pointwise @f$ \mu(\dot\gamma(u^n)) @f$, evaluated once per cell and step.
   */
  template <class FES, class GridFunction>
  inline void computeCellViscosity(const MeshType& mesh, const FES& fes,
    const GridFunction& u, const CarreauYasuda& cy, std::vector<Real>& mu)
  {
    mu.assign(mesh.getCellCount(), cy.mu0);
    for (auto it = mesh.getCell(); it; ++it)
    {
      const auto& verts = it->getVertices();
      if (verts.size() != 4)
        continue;
      Eigen::Matrix3d E, U;
      const auto x0 = mesh.getVertexCoordinates(verts[0]);
      Eigen::Vector3d u0;
      for (int c = 0; c < 3; ++c)
        u0(c) = u[fes.getGlobalIndex({0, verts[0]}, c)];
      for (int k = 1; k < 4; ++k)
      {
        const auto xk = mesh.getVertexCoordinates(verts[k]);
        for (int c = 0; c < 3; ++c)
        {
          E(c, k - 1) = xk(c) - x0(c);
          U(c, k - 1) = u[fes.getGlobalIndex({0, verts[k]}, c)] - u0(c);
        }
      }
      const Eigen::Matrix3d G = U * E.inverse();
      const Eigen::Matrix3d D = 0.5 * (G + G.transpose());
      const Real g2 = cy.gammaRegularization * cy.gammaRegularization +
        2.0 * D.cwiseProduct(D).sum();
      mu[it->getIndex()] = cy(std::sqrt(g2));
    }
  }

  inline MeshType loadMesh(const Config& cfg)
  {
    MeshType mesh;
    mesh.load(cfg.meshPath, IO::FileFormat::MEDIT);
    if (mesh.getSpaceDimension() != 3)
      throw std::runtime_error("Expected a 3D coronary lumen mesh.");
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
}

#endif
