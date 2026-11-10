/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>

#include <Rodin/Alert.h>

#include "Fluid.h"
#include "FlowIntegral.h"
#include "../Log.h"
#include "../../CoronaryArtery/CoronaryArteryTiming.h"

namespace Rodin::Examples::Heart::Atrium
{
  using namespace Rodin;
  using namespace Rodin::Geometry;
  using namespace Rodin::Variational;

  namespace
  {
    /// Refreshes the ghost entries after owned nodal values were written.
    void refreshGhosts(::Vec v)
    {
      PetscErrorCode ierr = VecGhostUpdateBegin(v, INSERT_VALUES, SCATTER_FORWARD);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecGhostUpdateEnd(v, INSERT_VALUES, SCATTER_FORWARD);
      assert(ierr == PETSC_SUCCESS);
      (void)ierr;
    }
  }

  bool Fluid::onWallNode(const Point& p) const
  {
    return m_wallArea.getValue(p) > m_wallAreaThreshold;
  }

  Fluid::Real Fluid::fibrinAt(const Point& p) const
  {
    if (m_cfg.kinetics.immobileFibrin)
      return std::max<Real>(m_fibrin.getValue(p), 0.0);
    const Real G0 = m_cfg.kinetics.fibrinogen0;
    return std::clamp<Real>(G0 - m_gCur.getValue(p), 0.0, G0);
  }

  Fluid::Real Fluid::dragAt(const Point& p) const
  {
    if (!m_speciesActive || !m_cfg.brinkman.enabled)
      return 0.0;
    return m_cfg.brinkman(fibrinAt(p));
  }

  Fluid::Real Fluid::crosswind(const Point& p, Real cur, Real prev,
    const Math::SpatialVector<Real>& gradient, Real diffusivity) const
  {
    const Real gn = std::sqrt(Math::dot(gradient, gradient));
    if (gn < 1.0e-14)
      return 0.0;
    const auto uc = (m_uOld - m_meshVelocity).getValue(p);
    const Real residual = (cur - prev) / m_run.timeStep + Math::dot(uc, gradient);
    return std::max<Real>(0.0,
      m_cfg.crosswindC * std::abs(residual) * cellSize(p) / (2.0 * gn) - diffusivity);
  }

  void Fluid::assignTransport(ScalarProblemType& problem, ScalarTrialFunctionType& c,
    ScalarTestFunctionType& v, ScalarGridFunctionType& cur, ScalarGridFunctionType& prev,
    Real D, Real inlet, Real robinScale, Real loadScale,
    const ScalarGridFunctionType* load)
  {
    const Real dt = m_run.timeStep;
    const Real scale = m_cfg.kinetics.stabilizationScale;

    // ALE transport in nonconservative form. Nodal concentrations move with
    // the mesh, so (c - c^n)/dt is the time derivative at fixed mesh
    // coordinate and the convecting velocity is u - w. The operator is the
    // same for every species and preserves constants, hence P + T - A.
    const auto transport = m_uOld - m_meshVelocity;

    // Cell Peclet numbers are of order 1e6: tau follows the transient and
    // advective scales.
    RealFunction tauFn = [this, D, scale, dt](const Point& p) {
      const auto uc = (m_uOld - m_meshVelocity).getValue(p);
      return scale * supgTau(cellSize(p), std::sqrt(Math::dot(uc, uc)), D, 0.0, dt);
    };
    RealFunction kdcFn = [this, D, &cur, &prev](const Point& p) {
      return crosswind(p, cur.getValue(p), prev.getValue(p), Grad(cur).getValue(p), D);
    };
    RealFunction invSpeedSqFn = [this](const Point& p) {
      const auto uc = (m_uOld - m_meshVelocity).getValue(p);
      const Real u2 = Math::dot(uc, uc);
      return (u2 > 1.0e-20) ? 1.0 / u2 : 0.0;
    };

    // Wall exchange without a reaction layer: an implicit Robin sink
    // robinScale E c of the consumed species and a load loadScale E f of the
    // produced one.
    RealFunction robinFn = [this, robinScale](const Point& p) {
      return robinScale * m_activation.getValue(p);
    };
    RealFunction loadFn = [this, loadScale, load](const Point& p) {
      return load ? loadScale * m_activation.getValue(p) * load->getValue(p) : 0.0;
    };

    // Half the discrete divergence of the advecting fluid velocity, for the
    // skew-symmetric form about the inflow state. On the moving domain the
    // balance of int c^2 contains div(u - w) from the advection and div w
    // from the domain, so the factor is div u.
    RealFunction halfDivFn = [this](const Point& p) -> Real {
      return 0.5 * Div(m_uOld).getValue(p);
    };
    RealFunction halfDivInletFn = [this, inlet](const Point& p) -> Real {
      return 0.5 * inlet * Div(m_uOld).getValue(p);
    };

    const auto streamline = Dot(transport, Grad(v));

    problem = (1.0 / dt) * SpeciesIntegral(c, v) - (1.0 / dt) * SpeciesIntegral(cur, v) +
      D * SpeciesIntegral(Grad(c), Grad(v)) +
      SpeciesIntegral(Dot(transport, Grad(c)), v)

      // SUPG on the transient and advective residual.
      + (1.0 / dt) * SpeciesIntegral(tauFn * c, streamline) -
      (1.0 / dt) * SpeciesIntegral(tauFn * cur, streamline) +
      SpeciesIntegral(tauFn * Dot(transport, Grad(c)), streamline)

      // Codina crosswind: SUPG leaves the crosswind direction free.
      + SpeciesIntegral(kdcFn * Grad(c), Grad(v)) -
      SpeciesIntegral(kdcFn * invSpeedSqFn * Dot(transport, Grad(c)), streamline)

      // Inflow composition at the pulmonary veins; the mitral outlet is free.
      + DirichletBC(c, RealFunction(Real(inlet))).on(m_inletSet);
    if (m_cfg.skewTransport)
    {
      // (1/2) int div(u_h) (c - c_in) v. The inflow state remains an exact
      // solution and a deviation from it is advected by an energy-neutral
      // operator, so the discrete divergence of the equal-order velocity
      // cannot amplify it.
      problem += SpeciesIntegral(halfDivFn * c, v);
      problem -= SpeciesIntegral(halfDivInletFn, v);
    }
    if (robinScale > 0)
      problem += BoundaryIntegral(robinFn * c, v).over(m_cfg.wall);
    if (load && loadScale > 0)
      problem -= BoundaryIntegral(loadFn * v).over(m_cfg.wall);
  }

  void Fluid::setupSpecies()
  {
    if (!m_cfg.solveKinetics)
      return;
    const auto& k = m_cfg.kinetics;

    // Wall exchange is a boundary flux only without a reaction layer;
    // otherwise initiation is the volumetric channel of the reaction step.
    // The direct channel j = j_0 E P/P_ref is a Robin sink of P and the same
    // amount as a load of T. With factor Xa it is replaced by the emission of X.
    const bool layer = k.layerThickness > 0.0;
    const Real robin = (layer || k.factorXa) ? 0.0 : k.wallFlux / k.prothrombinRef;
    assignTransport(m_transportP, m_pT, m_vP, m_pCur, m_pPrev, k.diffusivityProtein,
      k.prothrombin0, robin, 0.0, nullptr);
    assignTransport(m_transportT, m_tT, m_vT, m_tCur, m_tPrev, k.diffusivityProtein,
      k.inletThrombin, 0.0, robin, &m_pCur);
    assignTransport(m_transportA, m_aT, m_vA, m_aCur, m_aPrev, k.diffusivityProtein,
      k.antithrombin0, 0.0, 0.0, nullptr);
    assignTransport(m_transportG, m_gT, m_vG, m_gCur, m_gPrev, k.diffusivityFibrinogen,
      k.fibrinogen0, 0.0, 0.0, nullptr);
    if (k.factorXa)
    {
      // The veins carry no factor Xa. Without a layer the endocardium emits
      // -D d_n X = j_X E; with one, the reaction step emits it in the layer.
      assignTransport(m_transportX, m_xT, m_vX, m_xCur, m_xPrev, k.diffusivityProtein,
        0.0, 0.0, layer ? 0.0 : k.resolvedXaFlux(), &m_one);
    }
    if (layer)
    {
      // Screened Poisson extension of the wall activation into the lumen.
      const Real delta = k.layerThickness;
      m_extension = (delta * delta) * Integral(Grad(m_eT), Grad(m_vE)) +
        Integral(m_eT, m_vE) +
        DirichletBC(m_eT, RealFunction([this](const Point& p) -> Real {
          return m_activation.getValue(p);
        }))
          .on(m_cfg.wall);
    }
    const Real threshold = k.inhibitionRate * k.antithrombin0 * k.amplificationHalf /
      (k.amplificationRate * k.prothrombin0 / k.prothrombinRef -
        k.inhibitionRate * k.antithrombin0);
    if (isRoot())
    {
      Log("Fluid/Kinetics") << "k_on A_0=" << k.inhibitionRate * k.antithrombin0
                            << " 1/s activation_threshold=" << threshold
                            << " mol/m^3 layer_thickness=" << k.layerThickness << " m"
                            << Alert::Raise;
      Log("Fluid/Kinetics") << "factor_Xa=" << (k.factorXa ? "on" : "off")
                            << " j_X=" << k.resolvedXaFlux()
                            << " mol/m^2/s kappa=" << k.xaRate << " 1/s k_dX="
                            << k.xaDecay << " 1/s fibrin="
                            << (k.immobileFibrin ? "immobile" : "advected") << " skew="
                            << (m_cfg.skewTransport ? "on" : "off") << Alert::Raise;
    }
  }

  void Fluid::extendActivation()
  {
    const Real delta = m_cfg.kinetics.layerThickness;
    if (delta <= 0.0)
      return;
    m_extension.assemble();
    m_extension.solve(m_extensionKSP);
    m_layer.setData(m_eT.getSolution().getData());

    // The volumetric dose int (j_0/delta) scale E~ dV equals the wall dose
    // j_0 int E dA. The scale is one on a flat wall and below one where
    // layers of opposite walls overlap.
    m_flux = BoundaryIntegral(m_activation, m_qFlux).over(m_cfg.wall);
    m_flux.assemble();
    const Real wallDose = m_flux(m_one);
    m_flux = Integral(m_layer, m_qFlux);
    m_flux.assemble();
    const Real layerDose = m_flux(m_one);
    m_layerScale = layerDose > 0.0 ? delta * wallDose / layerDose : 0.0;
    if (isRoot())
    {
      Log("Fluid/Kinetics") << "wall_activation_area=" << wallDose
                            << " m^2 layer_volume=" << layerDose
                            << " m^3 dose_scale=" << m_layerScale << Alert::Raise;
    }
  }

  void Fluid::solveSpecies()
  {
    Stage stage("Species transport and reaction", isRoot());
    const auto speciesStart = CoronaryClock::now();
    m_speciesActive = true;

    // P is transported first so that a wall load of T reads the new P.
    const auto advance = [](ScalarProblemType& problem, Rodin::Solver::KSP& ksp,
                           ScalarTrialFunctionType& trial, ScalarGridFunctionType& cur,
                           ScalarGridFunctionType& prev) {
      problem.assemble();
      problem.solve(ksp);

      // The crosswind residual of the next step reads c^n and c^{n-1}.
      prev.setData(cur.getData());
      cur.setData(trial.getSolution().getData());
    };
    advance(m_transportP, m_kspP, m_pT, m_pCur, m_pPrev);
    advance(m_transportT, m_kspT, m_tT, m_tCur, m_tPrev);
    advance(m_transportA, m_kspA, m_aT, m_aCur, m_aPrev);
    advance(m_transportG, m_kspG, m_gT, m_gCur, m_gPrev);
    if (m_cfg.kinetics.factorXa)
      advance(m_transportX, m_kspX, m_xT, m_xCur, m_xPrev);
    reactionStep();
    updateCoagulationDiagnostics();
    m_timing.species = secondsSince(speciesStart);
  }

  void Fluid::reactionStep()
  {
    const auto& k = m_cfg.kinetics;
    const Real dt = m_run.timeStep;

    // Layer initiation inside the same conservative nodal step. The direct
    // channel is the rate s/P = (j_0/delta) scale E~/P_ref of P -> T; with
    // factor Xa the layer emits X at (j_X/delta) scale E~ instead.
    const bool hasLayer = k.layerThickness > 0.0;
    const Real layerRate = hasLayer && !k.factorXa
      ? k.wallFlux / k.layerThickness * m_layerScale / k.prothrombinRef
      : 0.0;
    const Real xaRate = hasLayer && k.factorXa
      ? k.resolvedXaFlux() / k.layerThickness * m_layerScale
      : 0.0;

    // The P1 fields share one layout: the owned nodal values are reacted in
    // place and the ghosts refreshed once. Factor Xa and the immobile fibrin
    // are written only when enabled.
    constexpr size_t Count = 6;
    std::array<::Vec, Count> fields{m_pCur.getData(), m_tCur.getData(), m_aCur.getData(),
      m_gCur.getData(), m_xCur.getData(), m_fibrin.getData()};
    std::array<PetscScalar*, Count> values{};
    const PetscScalar* layer = nullptr;
    PetscInt size = 0;
    PetscErrorCode ierr = VecGetLocalSize(fields[0], &size);
    assert(ierr == PETSC_SUCCESS);
    for (size_t i = 0; i < Count; ++i)
    {
      ierr = VecGetArray(fields[i], &values[i]);
      assert(ierr == PETSC_SUCCESS);
    }
    ierr = VecGetArrayRead(m_layer.getData(), &layer);
    assert(ierr == PETSC_SUCCESS);
    for (PetscInt i = 0; i < size; ++i)
    {
      const Real extended = std::max<Real>(PetscRealPart(layer[i]), 0.0);
      const auto next =
        k.advance({PetscRealPart(values[0][i]), PetscRealPart(values[1][i]),
                    PetscRealPart(values[2][i]), PetscRealPart(values[3][i]),
                    k.factorXa ? PetscRealPart(values[4][i]) : Real(0),
                    k.immobileFibrin ? PetscRealPart(values[5][i]) : Real(0)},
          dt, layerRate * extended, xaRate * extended);
      values[0][i] = next.P;
      values[1][i] = next.T;
      values[2][i] = next.A;
      values[3][i] = next.G;
      if (k.factorXa)
        values[4][i] = next.X;
      if (k.immobileFibrin)
        values[5][i] = next.F;
    }
    ierr = VecRestoreArrayRead(m_layer.getData(), &layer);
    assert(ierr == PETSC_SUCCESS);
    for (size_t i = 0; i < Count; ++i)
    {
      ierr = VecRestoreArray(fields[i], &values[i]);
      assert(ierr == PETSC_SUCCESS);
      refreshGhosts(fields[i]);
    }
    (void)ierr;
  }

  void Fluid::updateCoagulationDiagnostics()
  {
    const auto& k = m_cfg.kinetics;
    const Real invariant = k.prothrombin0 - k.antithrombin0;
    const PetscScalar *P = nullptr, *T = nullptr, *A = nullptr, *G = nullptr;
    PetscScalar *F = nullptr, *I = nullptr, *S = nullptr, *R = nullptr;
    PetscInt size = 0;
    PetscErrorCode ierr = VecGetLocalSize(m_pCur.getData(), &size);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGetArrayRead(m_pCur.getData(), &P);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGetArrayRead(m_tCur.getData(), &T);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGetArrayRead(m_aCur.getData(), &A);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGetArrayRead(m_gCur.getData(), &G);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGetArray(m_fibrin.getData(), &F);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGetArray(m_tat.getData(), &I);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGetArray(m_drag.getData(), &S);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGetArray(m_drift.getData(), &R);
    assert(ierr == PETSC_SUCCESS);
    for (PetscInt i = 0; i < size; ++i)
    {
      // The immobile fibrin is a state of the reaction step; otherwise the
      // fibrin is the consumed fibrinogen.
      const Real fibrin = k.immobileFibrin
        ? std::max<Real>(PetscRealPart(F[i]), 0.0)
        : std::clamp<Real>(k.fibrinogen0 - PetscRealPart(G[i]), 0.0, k.fibrinogen0);
      F[i] = fibrin;
      I[i] = k.antithrombin0 - PetscRealPart(A[i]);
      S[i] = m_cfg.brinkman(fibrin);
      R[i] = PetscRealPart(P[i]) + PetscRealPart(T[i]) - PetscRealPart(A[i]) - invariant;
    }
    ierr = VecRestoreArrayRead(m_pCur.getData(), &P);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecRestoreArrayRead(m_tCur.getData(), &T);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecRestoreArrayRead(m_aCur.getData(), &A);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecRestoreArrayRead(m_gCur.getData(), &G);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecRestoreArray(m_fibrin.getData(), &F);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecRestoreArray(m_tat.getData(), &I);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecRestoreArray(m_drag.getData(), &S);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecRestoreArray(m_drift.getData(), &R);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
    for (::Vec v :
      {m_fibrin.getData(), m_tat.getData(), m_drag.getData(), m_drift.getData()})
      refreshGhosts(v);
  }
}
