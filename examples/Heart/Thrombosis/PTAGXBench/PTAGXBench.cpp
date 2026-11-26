/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file PTAGXBench.cpp
 * @brief P-T-A-G + platelets on a channel or strip with a prescribed
 *        Poiseuille velocity (see PTAGXBench.h). Used by the validation tests in
 *        examples/Heart/Thrombosis/Test.
 *
 * Options (prefix -pg_):
 *   mesh, mesh_scale, xdmf, csv, dt, t_end, csv_interval, output_interval,
 *   log_every, inlet, reactive, collagen, tm (attribute lists, comma
 *   separated; "none" for an empty set), activation, H, gamma, gamma2,
 *   t_switch, crosswind, stop_thrombin, wall_flux, k_tm, k_a, theta0,
 *   immobile_fibrin, F_cap, fibrin_monomer, k_pol, xa, j_x, kappa_x, k_dx, platelets, N0, k_p, T_p, S_max, k_w, k_wN, k_agg,
 *   l_agg, D_pl0, D_pl1, P0, A0, G0, supg.
 */
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/collectives.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <Rodin/Alert.h>
#include <Rodin/Configure.h>

#ifdef RODIN_USE_SCOTCH
#include <Rodin/Scotch/MeshPartitioner.h>
#endif

#include "PTAGXBench.h"
#include "CoronaryArtery/ThrombosisModel.h"

namespace Rodin::Examples::Heart::Thrombosis
{
  using namespace Rodin;
  using namespace Rodin::Math;
  using namespace Rodin::Solver;
  using namespace Rodin::Geometry;
  using namespace Rodin::Variational;

  namespace
  {
    constexpr int RootRank = 0;

    void setPrefixedDefault(
      const std::string& prefix, const char* suffix, const char* value)
    {
      const std::string name = "-" + prefix + suffix;
      PetscBool set = PETSC_FALSE;
      PetscErrorCode ierr =
        PetscOptionsHasName(PETSC_NULLPTR, PETSC_NULLPTR, name.c_str(), &set);
      assert(ierr == PETSC_SUCCESS);
      if (!set)
      {
        ierr = PetscOptionsSetValue(PETSC_NULLPTR, name.c_str(), value);
        assert(ierr == PETSC_SUCCESS);
      }
      (void)ierr;
    }
  }

  // ==========================================================================
  // Construction
  // ==========================================================================
  PTAGXBench::PTAGXBench(const Context::MPI& context, const Config& cfg)
    : m_cfg(cfg),
      m_mesh(makeMesh(context, m_cfg)),
      m_xdmf(context.getCommunicator(), m_cfg.xdmfBasename),
      m_vh(std::integral_constant<size_t, 1>{}, m_mesh, m_mesh.getSpaceDimension()),
      m_sh(std::integral_constant<size_t, 1>{}, m_mesh),
      m_u(m_vh),
      m_pT(m_sh), m_tT(m_sh), m_aT(m_sh), m_gT(m_sh), m_nT(m_sh), m_mT(m_sh),
      m_xT(m_sh), m_fmT(m_sh),
      m_vP(m_sh), m_vT(m_sh), m_vA(m_sh), m_vG(m_sh), m_vN(m_sh), m_vM(m_sh),
      m_vX(m_sh), m_vFm(m_sh),
      m_p(m_sh), m_t(m_sh), m_a(m_sh), m_g(m_sh), m_n(m_sh), m_m(m_sh), m_x(m_sh),
      m_fm(m_sh),
      m_pPrev(m_sh), m_tPrev(m_sh), m_aPrev(m_sh), m_gPrev(m_sh), m_nPrev(m_sh),
      m_mPrev(m_sh), m_xPrev(m_sh), m_fmPrev(m_sh),
      m_f(m_sh), m_s(m_sh), m_sSmooth(m_sh),
      m_pNext(m_sh), m_tNext(m_sh), m_aNext(m_sh), m_gNext(m_sh), m_fNext(m_sh),
      m_nNext(m_sh), m_mNext(m_sh), m_sNext(m_sh), m_xNext(m_sh), m_fmNext(m_sh),
      m_lumped(m_sh), m_deposit(m_sh),
      m_one(m_sh), m_work(m_sh), m_stab(m_sh), m_theta(m_sh), m_tat(m_sh), m_surf(m_sh),
      m_q(m_sh), m_form(m_q),
      m_transportP(m_pT, m_vP), m_transportT(m_tT, m_vT), m_transportA(m_aT, m_vA),
      m_transportG(m_gT, m_vG), m_transportN(m_nT, m_vN), m_transportM(m_mT, m_vM),
      m_transportX(m_xT, m_vX), m_transportFm(m_fmT, m_vFm),
      m_kspP(m_transportP), m_kspT(m_transportT), m_kspA(m_transportA),
      m_kspG(m_transportG), m_kspN(m_transportN), m_kspM(m_transportM),
      m_kspX(m_transportX), m_kspFm(m_transportFm),
      m_smoothTrial(m_sh), m_smoothTest(m_sh), m_smooth(m_smoothTrial, m_smoothTest),
      m_smoothKSP(m_smooth)
  {
    // Platelets adhere on the tissue-factor wall and on exposed collagen.
    m_adhesive = m_cfg.reactive;
    for (const auto& a : m_cfg.collagen)
      m_adhesive.insert(a);

    if (isRoot())
    {
      Alert::Info() << "[mesh] dimension=" << m_mesh.getDimension()
                    << " cells=" << m_mesh.getCellCount()
                    << " vertices=" << m_mesh.getVertexCount() << Alert::Raise;
    }
  }

  bool PTAGXBench::isRoot() const
  {
    return m_mesh.getContext().getCommunicator().rank() == RootRank;
  }

  PTAGXBench::MeshType PTAGXBench::makeMesh(const Context::MPI& context, const Config& cfg)
  {
    const auto& comm = context.getCommunicator();

    Rodin::MPI::Sharder sharder(context);
    if (comm.rank() == RootRank)
    {
      Geometry::Mesh<Context::Local> mesh;
      mesh.load(cfg.meshPath, IO::FileFormat::MEDIT);
      const size_t D = mesh.getDimension();
      mesh.getConnectivity().compute(D, D);
      mesh.getConnectivity().compute(D, 0);
      mesh.getConnectivity().compute(D, D - 1);
      mesh.getConnectivity().compute(D - 1, D);
      mesh.getConnectivity().compute(D - 1, 0);

#ifdef RODIN_USE_SCOTCH
      Scotch::Partitioner partitioner(mesh);
#else
      Geometry::BalancedCompactPartitioner partitioner(mesh);
#endif
      partitioner.partition(static_cast<size_t>(comm.size()));
      sharder.shard(partitioner);
      sharder.scatter(RootRank);
    }

    MeshType mesh = sharder.gather(RootRank);
    mesh.scale(cfg.meshScale);

    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);
    mesh.getConnectivity().compute(D, 0);
    mesh.getConnectivity().compute(D, D - 1);
    mesh.getConnectivity().compute(D - 1, D);
    mesh.getConnectivity().compute(D - 1, 0);
    if (D == 3)
    {
      mesh.getConnectivity().compute(D - 1, 1);
      mesh.getConnectivity().compute(1, 0);
      mesh.reconcile(2);
    }
    mesh.reconcile(1);
    return mesh;
  }

  PTAGXBench::Real PTAGXBench::cellSize(const Point& p)
  {
    return std::pow(p.getPolytope().getMeasure(), 1.0 / p.getPolytope().getDimension());
  }

  // ==========================================================================
  // Velocity and pointwise coefficients
  // ==========================================================================
  PTAGXBench::Real PTAGXBench::wallShearRate(Real t) const
  {
    if (m_cfg.switchTime >= 0.0 && t > m_cfg.switchTime)
      return m_cfg.wallShearRate2;
    return m_cfg.wallShearRate;
  }

  void PTAGXBench::updateVelocity(Real t)
  {
    const Real gamma = wallShearRate(t);
    if (gamma == m_gamma)
      return;
    m_gamma = gamma;
    const Real H = m_cfg.channelHeight;
    const size_t dim = m_mesh.getSpaceDimension();
    m_u.project(VectorFunction(dim,
      [gamma, H, dim](const Point& p) -> Math::SpatialVector<Real> {
        Math::SpatialVector<Real> out(dim);
        for (size_t c = 0; c < dim; ++c)
          out(c) = 0.0;
        const Real y = std::clamp<Real>(p.y(), 0.0, H);
        out(0) = gamma * y * (1.0 - y / H);
        return out;
      }));
    updateCellCache();
    if (isRoot())
    {
      Alert::Info() << "[flow] t=" << t << " s: wall shear rate " << gamma
                    << " 1/s, centreline velocity " << (0.25 * gamma * H) << " m/s"
                    << Alert::Raise;
    }
  }

  void PTAGXBench::updateCellCache()
  {
    const size_t D = m_mesh.getDimension();
    const size_t count = m_mesh.getPolytopeCount(D);
    m_cellH.assign(count, 0.0);
    m_cellSpeed2.assign(count, 0.0);
    m_cellShear.assign(count, 0.0);
    Math::SpatialPoint rc(D);
    for (size_t c = 0; c < D; ++c)
      rc(c) = 1.0 / static_cast<Real>(D + 1);
    const auto sym = 0.5 * (Jacobian(m_u) + Transpose(Jacobian(m_u)));
    for (auto it = m_mesh.getCell(); !it.end(); ++it)
    {
      const Point p(*it, rc);
      const Index i = it->getIndex();
      const auto uc = m_u.getValue(p);
      m_cellH[i] = cellSize(p);
      m_cellSpeed2[i] = Math::dot(uc, uc);
      m_cellShear[i] = std::sqrt(std::max<Real>(2.0 * Dot(sym, sym).getValue(p), 0.0));
    }
  }

  Index PTAGXBench::cellOf(const Point& p) const
  {
    assert(p.getPolytope().getDimension() == m_mesh.getDimension());
    return p.getPolytope().getIndex();
  }

  PTAGXBench::Real PTAGXBench::plateletDiffusivityAt(const Point& p) const
  {
    const auto& k = m_cfg.kinetics;
    return k.plateletDiffusivity0 + k.plateletDiffusivity1 * m_cellShear[cellOf(p)];
  }

  PTAGXBench::Real PTAGXBench::crosswind(const Point& p, Real cur, Real prev,
    const Math::SpatialVector<Real>& gradient, Real diffusivity) const
  {
    if (m_cfg.crosswindC <= 0.0)
      return 0.0;
    const Real gn = std::sqrt(Math::dot(gradient, gradient));
    if (gn < 1.0e-14)
      return 0.0;
    const auto uc = m_u.getValue(p);
    const Real h = m_cellH[cellOf(p)];
    const Real residual = (cur - prev) / m_cfg.dt + Math::dot(uc, gradient);
    // |R|/|grad c| is a velocity; capped by |u| + h/dt as in AneurysmPTAG so
    // that a field uniform up to round-off cannot produce an unbounded
    // diffusivity.
    const Real speed = std::min<Real>(std::abs(residual) / gn,
      std::sqrt(Math::dot(uc, uc)) + h / m_cfg.dt);
    return std::max<Real>(0.0, m_cfg.crosswindC * speed * h / 2.0 - diffusivity);
  }

  bool PTAGXBench::onWall(const Point& p, const FlatSet<Attribute>& set) const
  {
    if (set.empty())
      return false;
    const auto a = p.getPolytope().getAttribute();
    return a.has_value() && set.count(*a) > 0;
  }

  void PTAGXBench::setupSurfaceLayer()
  {
    const auto& k = m_cfg.kinetics;
    m_surf = Real(0);
    m_surfScale = 0.0;
    if (k.surfaceAmplification <= 0.0)
      return;
    const Real d = k.surfaceLayer;
    const Real x0 = m_cfg.patch[0], x1 = m_cfg.patch[1];
    m_surf.project(RealFunction([d, x0, x1](const Point& p) -> Real {
      const Real dx = std::max<Real>({ x0 - p.x(), Real(0), p.x() - x1 });
      return std::exp(-std::hypot(dx, std::max<Real>(p.y(), 0.0)) / d);
    }));
    // patch length inside the mesh
    Real xmin = 1e30, xmax = -1e30;
    m_work.project(RealFunction([](const Point& p) -> Real { return p.x(); }));
    xmin = m_work.min();
    xmax = m_work.max();
    const Real length = std::max<Real>(std::min(x1, xmax) - std::max(x0, xmin), 0.0);
    Real vol = integrate(m_surf);
    if (m_cfg.axisymmetric)
    {
      // disc of radius x1: int w r dr dy against |patch|/(2 pi) = x1^2/2
      m_work.project(RealFunction([this](const Point& p) -> Real {
        return std::max<Real>(p.x(), 0.0) * m_surf.getValue(p);
      }));
      vol = integrate(m_work) / (0.5 * x1);
    }
    m_surfScale = (vol > 0.0) ? (k.surfaceAmplification / d) * (d * length / vol) : 0.0;
    Alert::Info() << "[surface] sigma=" << k.surfaceAmplification << " m, delta_s=" << d
                  << " m, patch " << length * 1e6 << " um, chi=" << (d * length / vol)
                  << ", theta at the wall=" << m_surfScale << Alert::Raise;
  }

  PTAGXBench::Real PTAGXBench::surfaceThetaAt(const Point& p) const
  {
    return m_surfScale > 0.0 ? m_surfScale * m_surf.getValue(p) : 0.0;
  }

  PTAGXBench::Real PTAGXBench::adhesionActivatedAt(const Point& p) const
  {
    const auto& k = m_cfg.kinetics;
    Real velocity = 0.0;
    if (onWall(p, m_cfg.collagen))
      velocity = k.adhesionActivated;
    else if (onWall(p, m_cfg.reactive))
      velocity = k.adhesionActivated * m_cfg.activation;
    if (velocity <= 0.0)
      return 0.0;
    return velocity * std::max<Real>(0.0, 1.0 - m_s.getValue(p) / k.plateletMax);
  }

  PTAGXBench::Real PTAGXBench::adhesionRestingAt(const Point& p) const
  {
    const auto& k = m_cfg.kinetics;
    if (k.adhesionResting <= 0.0 || !onWall(p, m_cfg.collagen))
      return 0.0;
    return k.adhesionResting * std::max<Real>(0.0, 1.0 - m_s.getValue(p) / k.plateletMax);
  }

  // ==========================================================================
  // Setup
  // ==========================================================================
  template <class DFn, class RFn, class LFn>
  void PTAGXBench::assignTransport(int species, ScalarProblemType& problem,
    ScalarTrialFunctionType& c, ScalarTestFunctionType& v,
    ScalarGridFunctionType& cur, ScalarGridFunctionType& prev,
    const DFn& diffusivity, Real inlet, const FlatSet<Attribute>& wall,
    const RFn& robin, const LFn& load, Real sink)
  {
    (void)prev;
    m_speciesDiffusivity[species] = diffusivity;
    const Real dt = m_cfg.dt;
    const Real scale = m_cfg.kinetics.stabilizationScale;

    RealFunction dFn = [diffusivity](const Point& p) { return diffusivity(p); };
    // tau_K from the cell cache (h_K and |u| at the centroid): constant per
    // cell, as the SUPG parameter is by construction.
    RealFunction tauFn = [this, diffusivity, scale, dt](const Point& p) {
      const Index i = cellOf(p);
      return scale * supgTau(m_cellH[i], std::sqrt(m_cellSpeed2[i]), diffusivity(p),
                             0.0, dt);
    };
    // Crosswind coefficient and 1/|u|^2 cellwise (centroid), as the Codina
    // operator is defined elementwise; refreshed before each assembly.
    RealFunction kdcFn = [this, species](const Point& p) -> Real {
      return m_cellCrosswind[species][cellOf(p)];
    };
    RealFunction invSpeedSqFn = [this](const Point& p) -> Real {
      const Real u2 = m_cellSpeed2[cellOf(p)];
      return (u2 > 1.0e-20) ? 1.0 / u2 : 0.0;
    };
    RealFunction robinFn = [robin](const Point& p) { return robin(p); };

    // An empty attribute set means "the whole boundary" to the assembler.
    // The wall coefficients vanish off their sets, so the integration is
    // still exact there; the thrombomodulin sink, a constant, must be
    // switched off explicitly.
    const FlatSet<Attribute>& wallSet = wall.empty() ? m_cfg.reactive : wall;
    const Real sinkEff = m_cfg.thrombomodulin.empty() ? 0.0 : sink;
    const FlatSet<Attribute>& sinkSet =
      m_cfg.thrombomodulin.empty() ? m_cfg.reactive : m_cfg.thrombomodulin;

    const Real rdt = 1.0 / dt;
    const auto pv = Dot(m_u, Grad(v));
    // Streamline-upwind Petrov-Galerkin: the transient and convective residual
    // tested against v + tau u.grad v in one integrator (the reaction is in
    // the nodal step and, with P1, the diffusion drops out of the residual).
    // The Galerkin mass and convection stay in their own integrators: folded
    // into the Petrov-Galerkin weight, the non-polynomial tau lowers the
    // quadrature chosen for the mass matrix (a 3% error in a quiescent test
    // against an exact P1 solution).
    const auto residual = rdt * c + Dot(m_u, Grad(c));
    const auto weight = tauFn * pv;

    // Quadrature orders stated explicitly wherever a lambda coefficient hides
    // the degree from inference (Integrator::setOrder): the SUPG right-hand
    // side (c^n u.grad v, degree 2) would otherwise be integrated at order 1
    // while its bilinear counterpart is integrated at order 2, a spurious
    // source in a uniform field; the wall terms carry (1 - S/S_max) and P,
    // degree 3, and must match the deposition form exactly.
    constexpr size_t SUPGOrder = 2;
    constexpr size_t WallOrder = 3;
    auto supgLhs = Integral(residual, weight);
    supgLhs.setOrder(SUPGOrder);
    auto supgRhs = Integral(rdt * cur, weight);
    supgRhs.setOrder(SUPGOrder);
    auto wallRobin = BoundaryIntegral(robinFn * c, v);
    wallRobin.over(wallSet);
    wallRobin.setOrder(WallOrder);
    auto wallLoad = BoundaryIntegral(load * v);
    wallLoad.over(wallSet);
    wallLoad.setOrder(WallOrder);

    const bool quiescent = (m_cfg.wallShearRate == 0.0 && m_cfg.wallShearRate2 == 0.0);
    if (quiescent && m_cfg.axisymmetric)
    {
      // Weak form in cylindrical coordinates (no swirl): dV = r dr dy, dA = r dr
      // on y = 0. Orders stated explicitly: the weight r hides the degree (2 is the
      // highest order rule used in this file for cells, 3 for edges).
      // The scalars 1/dt and k_TM are folded into the weights: a named
      // integrator scaled by a scalar is held by reference by the problem.
      // The weight is r / (100 um): the equation is unchanged, and the matrix
      // entries stay above the zero-pivot tolerance of the factorization.
      constexpr Real rRef = 1.0e-4;
      RealFunction rFn = [](const Point& p) { return std::max<Real>(p.x(), 0.0) / rRef; };
      RealFunction rDtFn = [rdt](const Point& p) {
        return rdt * std::max<Real>(p.x(), 0.0) / rRef;
      };
      RealFunction rdFn = [diffusivity](const Point& p) {
        return std::max<Real>(p.x(), 0.0) / rRef * diffusivity(p);
      };
      RealFunction rRobinFn = [robin, sinkEff](const Point& p) {
        return std::max<Real>(p.x(), 0.0) / rRef * (robin(p) + sinkEff);
      };
      auto mass = Integral(rDtFn * c, v);
      mass.setOrder(SUPGOrder);
      auto massPrev = Integral(rDtFn * cur, v);
      massPrev.setOrder(SUPGOrder);
      auto stiff = Integral(rdFn * Grad(c), Grad(v));
      stiff.setOrder(2);
      // Robin and thrombomodulin sink share one integrator: same wall set.
      if (sinkEff != 0.0 && sinkSet != wallSet)
        throw std::runtime_error("-pg_axisym: thrombomodulin wall must be the reactive wall");
      auto robinR = BoundaryIntegral(rRobinFn * c, v);
      robinR.over(wallSet);
      robinR.setOrder(WallOrder);
      auto loadR = BoundaryIntegral(rFn * load * v);
      loadR.over(wallSet);
      loadR.setOrder(WallOrder);
      problem =
          mass - massPrev
        + stiff
        + robinR - loadR
        + DirichletBC(c, RealFunction(Real(inlet))).on(m_cfg.inlet);
    }
    else if (quiescent)
    {
      problem =
          rdt * Integral(c, v) - rdt * Integral(cur, v)
        + Integral(dFn * Grad(c), Grad(v))
        + wallRobin - wallLoad
        + sinkEff * BoundaryIntegral(c, v).over(sinkSet)
        + DirichletBC(c, RealFunction(Real(inlet))).on(m_cfg.inlet);
    }
    else if (m_cfg.crosswindC > 0.0)
    {
      problem =
          rdt * Integral(c, v) - rdt * Integral(cur, v)
        + Integral(Dot(m_u, Grad(c)), v)
        + supgLhs - supgRhs
        + Integral(dFn * Grad(c), Grad(v))
        + wallRobin - wallLoad
        + sinkEff * BoundaryIntegral(c, v).over(sinkSet)
        // Codina crosswind: SUPG leaves the crosswind direction free.
        + Integral(kdcFn * Grad(c), Grad(v))
        - Integral(kdcFn * invSpeedSqFn * Dot(m_u, Grad(c)), Dot(m_u, Grad(v)))
        + DirichletBC(c, RealFunction(Real(inlet))).on(m_cfg.inlet);
    }
    else
    {
      problem =
          rdt * Integral(c, v) - rdt * Integral(cur, v)
        + Integral(Dot(m_u, Grad(c)), v)
        + supgLhs - supgRhs
        + Integral(dFn * Grad(c), Grad(v))
        + wallRobin - wallLoad
        + sinkEff * BoundaryIntegral(c, v).over(sinkSet)
        + DirichletBC(c, RealFunction(Real(inlet))).on(m_cfg.inlet);
    }
  }

  void PTAGXBench::setupSpecies()
  {
    const auto& k = m_cfg.kinetics;
    const Real E = m_cfg.activation;
    // With the initiation enzyme the wall emits X only; the legacy direct
    // exchange P -> T at the wall is switched off.
    const Real beta = k.factorXa ? 0.0 : k.wallFlux / k.prothrombinRef * E;
    const Real Dp = k.diffusivityProtein;
    const Real Dg = k.diffusivityFibrinogen;

    const auto protein = [Dp](const Point&) { return Dp; };
    const auto fibrinogen = [Dg](const Point&) { return Dg; };
    const auto zero = [](const Point&) { return 0.0; };
    const auto wallP = [beta](const Point&) { return beta; };
    // The wall load of T is the P that the implicit Robin term removed, written
    // with the grid function itself so that both sides of the exchange use the
    // same quadrature (P + T - A then balances to round-off).
    const auto loadT = beta * m_p;
    const auto noLoad = RealFunction(Real(0));

    assignTransport(0, m_transportP, m_pT, m_vP, m_p, m_pPrev, protein, k.prothrombin0,
      m_cfg.reactive, wallP, noLoad, 0.0);
    assignTransport(1, m_transportT, m_tT, m_vT, m_t, m_tPrev, protein, 0.0,
      m_cfg.reactive, zero, loadT, k.thrombomodulinRate);
    assignTransport(2, m_transportA, m_aT, m_vA, m_a, m_aPrev, protein, k.antithrombin0,
      m_cfg.reactive, zero, noLoad, 0.0);
    assignTransport(3, m_transportG, m_gT, m_vG, m_g, m_gPrev, fibrinogen, k.fibrinogen0,
      m_cfg.reactive, zero, noLoad, 0.0);

    if (k.factorXa)
    {
      // Neumann emission of Xa by the tissue-factor wall.
      const auto loadX = RealFunction(Real(k.xaFlux * E));
      assignTransport(6, m_transportX, m_xT, m_vX, m_x, m_xPrev, protein, 0.0,
        m_cfg.reactive, zero, loadX, 0.0);
    }
    if (k.immobileFibrin && k.fibrinMonomer)
    {
      assignTransport(7, m_transportFm, m_fmT, m_vFm, m_fm, m_fmPrev, fibrinogen, 0.0,
        m_cfg.reactive, zero, noLoad, 0.0);
    }

    if (!k.platelets)
      return;

    const auto platelet = [this](const Point& p) { return plateletDiffusivityAt(p); };
    // Same coefficient functions as depositAtWall(), so that what the Robin
    // terms remove is exactly what the wall nodes receive.
    const auto wallM = [this](const Point& p) { return adhesionActivatedAt(p); };
    const auto wallN = [this](const Point& p) { return adhesionRestingAt(p); };
    assignTransport(4, m_transportN, m_nT, m_vN, m_n, m_nPrev, platelet, k.plateletCount,
      m_cfg.collagen, wallN, noLoad, 0.0);
    assignTransport(5, m_transportM, m_mT, m_vM, m_m, m_mPrev, platelet, 0.0,
      m_adhesive, wallM, noLoad, 0.0);

    const Real l = k.aggregationLength;
    if (l > 0.0)
    {
      m_smooth = (l * l) * Integral(Grad(m_smoothTrial), Grad(m_smoothTest))
               + Integral(m_smoothTrial, m_smoothTest)
               - Integral(m_s, m_smoothTest);
    }
  }

  PTAGXBench& PTAGXBench::initialize()
  {
    const auto& k = m_cfg.kinetics;

    m_one = Real(1);
    for (auto* g : { &m_work, &m_stab, &m_theta, &m_tat, &m_deposit, &m_f, &m_s,
                     &m_sSmooth, &m_fNext, &m_sNext, &m_x, &m_xPrev, &m_xNext, &m_fm,
                     &m_fmPrev, &m_fmNext })
      *g = Real(0);
    for (auto* g : { &m_p, &m_pPrev, &m_pNext }) *g = Real(k.prothrombin0);
    for (auto* g : { &m_t, &m_tPrev, &m_tNext }) *g = Real(0);
    for (auto* g : { &m_a, &m_aPrev, &m_aNext }) *g = Real(k.antithrombin0);
    for (auto* g : { &m_g, &m_gPrev, &m_gNext }) *g = Real(k.fibrinogen0);
    const Real N0 = k.platelets ? k.plateletCount : 0.0;
    for (auto* g : { &m_n, &m_nPrev, &m_nNext }) *g = Real(N0);
    for (auto* g : { &m_m, &m_mPrev, &m_mNext }) *g = Real(0);

    const auto configureTransport = [](Rodin::Solver::KSP& ksp, const char* prefix) {
      ksp.setPrefix(prefix);
    };
    configureTransport(m_kspP, "pg_P_");
    configureTransport(m_kspT, "pg_T_");
    configureTransport(m_kspA, "pg_A_");
    configureTransport(m_kspG, "pg_G_");
    configureTransport(m_kspN, "pg_N_");
    configureTransport(m_kspM, "pg_M_");
    configureTransport(m_kspX, "pg_X_");
    configureTransport(m_kspFm, "pg_Fm_");
    setPrefixedDefault("pg_smooth_", "ksp_type", "cg");
    setPrefixedDefault("pg_smooth_", "pc_type", "jacobi");
    setPrefixedDefault("pg_smooth_", "ksp_rtol", "1e-10");
    m_smoothKSP.setPrefix("pg_smooth_");

    m_gamma = -1.0;
    updateVelocity(0.0);

    // Nodal lumped mass, int phi_i: the volume that a wall deposit fills.
    m_form = Integral(m_one, m_q);
    m_form.assemble();
    m_lumped.setData(m_form.getVector());

    setupSpecies();
    setupSurfaceLayer();

    m_xdmf.setMesh(m_mesh);
    m_xdmf.add("velocity", m_u);
    m_xdmf.add("prothrombin", m_p);
    m_xdmf.add("thrombin", m_t);
    m_xdmf.add("antithrombin", m_a);
    m_xdmf.add("fibrinogen", m_g);
    m_xdmf.add("fibrin", m_f);
    if (k.factorXa)
      m_xdmf.add("factorXa", m_x);
    if (k.immobileFibrin && k.fibrinMonomer)
      m_xdmf.add("fibrinMonomer", m_fm);
    if (k.platelets)
    {
      m_xdmf.add("restingPlatelets", m_n);
      m_xdmf.add("activatedPlatelets", m_m);
      m_xdmf.add("depositedPlatelets", m_s);
      m_xdmf.add("stabilizedPlatelets", m_stab);
    }
    m_xdmf.add("membrane", m_theta);

    const Real Tth = k.inhibitionRate * k.antithrombin0 * k.amplificationHalf /
      (k.amplificationRate * k.prothrombin0 / k.prothrombinRef -
       k.inhibitionRate * k.antithrombin0);
    const Real lw = std::sqrt(k.diffusivityProtein / (k.inhibitionRate * k.antithrombin0));
    const Real DaTM = k.thrombomodulinRate * lw / k.diffusivityProtein;
    if (isRoot())
    {
      Alert::Info() << "[model] platelets " << (k.platelets ? "on" : "off")
                    << ", theta_0=" << k.membraneBackground << ", fibrin "
                    << (k.immobileFibrin ? "immobile" : "advected (F = G0 - G)")
                    << (k.immobileFibrin && k.fibrinMonomer ? " via soluble monomer" : "")
                    << ", F_cap/G0=" << (k.fibrinCap > 0.0 ? k.fibrinCap / k.fibrinogen0 : 0.0)
                    << " (0: none), initiation "
                    << (k.factorXa ? "by Xa (j_X=" : "by direct wall flux (j_0=")
                    << (k.factorXa ? k.xaFlux : k.wallFlux) << " mol/m^2/s)"
                    << (k.factorXa ? ", kappa=" : "")
                    << (k.factorXa ? std::to_string(k.xaRate) + " 1/s, k_dX=" +
                                     std::to_string(k.xaDecay) + " 1/s" : std::string())
                    << ", k_TM="
                    << k.thrombomodulinRate << " m/s on "
                    << m_cfg.thrombomodulin.size() << " attribute(s) (Da_TM="
                    << DaTM << "), collagen on " << m_cfg.collagen.size()
                    << " attribute(s), T_th(theta=1)=" << (Tth * 1e6)
                    << " nM, crosswind C=" << m_cfg.crosswindC << Alert::Raise;
    }

    if (isRoot())
    {
      m_csv.open(m_cfg.csvPath);
      m_csv << std::setprecision(12);
      if (!m_csv)
        throw std::runtime_error("Failed to open " + m_cfg.csvPath);
      writeCSVHeader();
    }
    m_initialized = true;
    return *this;
  }

  // ==========================================================================
  // Time stepping
  // ==========================================================================
  void PTAGXBench::updateCrosswind(int species, const ScalarGridFunctionType& cur,
    const ScalarGridFunctionType& prev)
  {
    auto& kdc = m_cellCrosswind[species];
    const size_t D = m_mesh.getDimension();
    kdc.assign(m_mesh.getPolytopeCount(D), 0.0);
    if (m_cfg.crosswindC <= 0.0)
      return;
    Math::SpatialPoint rc(D);
    for (size_t c = 0; c < D; ++c)
      rc(c) = 1.0 / static_cast<Real>(D + 1);
    for (auto it = m_mesh.getCell(); !it.end(); ++it)
    {
      const Point p(*it, rc);
      // Full residual of the step just completed, (c^n - c^{n-1})/dt +
      // u.grad c^n, as in Codina's definition.
      kdc[it->getIndex()] = crosswind(p, cur.getValue(p), prev.getValue(p),
                                      Grad(cur).getValue(p),
                                      m_speciesDiffusivity[species](p));
    }
  }

  void PTAGXBench::solveTransport(int species, ScalarProblemType& problem,
    Rodin::Solver::KSP& ksp, ScalarTrialFunctionType& trial, ScalarGridFunctionType& cur,
    ScalarGridFunctionType& prev)
  {
    using Clock = std::chrono::steady_clock;
    const auto t0 = Clock::now();
    updateCrosswind(species, cur, prev);
    prev.setData(cur.getData());
    problem.assemble();
    const auto t1 = Clock::now();
    problem.solve(ksp);
    cur.setData(trial.getSolution().getData());
    const auto t2 = Clock::now();
    m_tAssembly += std::chrono::duration<Real>(t1 - t0).count();
    m_tTransport += std::chrono::duration<Real>(t2 - t1).count();
  }

  void PTAGXBench::depositAtWall()
  {
    const auto& k = m_cfg.kinetics;
    const Real dt = m_cfg.dt;
    const Real Smax = k.plateletMax;

    if (m_adhesive.empty())
      return;

    // What the implicit Robin terms removed from M and N in this step,
    // tested against each wall basis function. Same coefficient functions
    // and the same grid functions as the Robin terms of the M and N
    // transports, so that the quadrature matches.
    const RealFunction freeM = [this](const Point& p) { return adhesionActivatedAt(p); };
    const RealFunction freeN = [this](const Point& p) { return adhesionRestingAt(p); };
    auto deposition = BoundaryIntegral(freeM * m_m + freeN * m_n, m_q);
    deposition.over(m_adhesive);
    deposition.setOrder(3);
    m_form = deposition;
    m_form.assemble();
    m_deposit.setData(m_form.getVector());

    // S_i += dt b_i / m_i, capped at S_max; the excess goes back to M at the
    // same node so that N + M + S is conserved exactly.
    m_sNext.project(RealFunction([this, dt, Smax](const Point& p) -> Real {
      const Real mass = m_lumped.getValue(p);
      const Real s = m_s.getValue(p);
      if (mass <= 0.0)
        return s;
      return std::min<Real>(Smax, s + dt * m_deposit.getValue(p) / mass);
    }));
    m_work.project(RealFunction([this, dt](const Point& p) -> Real {
      const Real mass = m_lumped.getValue(p);
      if (mass <= 0.0)
        return 0.0;
      return std::max<Real>(0.0,
        m_s.getValue(p) + dt * m_deposit.getValue(p) / mass - m_sNext.getValue(p));
    }));
    m_clipped += integrate(m_work);
    m_mNext.project(RealFunction([this](const Point& p) -> Real {
      return m_m.getValue(p) + m_work.getValue(p);
    }));
    m_m.setData(m_mNext.getData());
    m_s.setData(m_sNext.getData());
  }

  void PTAGXBench::smoothDeposit()
  {
    if (m_cfg.kinetics.aggregationLength <= 0.0)
    {
      m_sSmooth.setData(m_s.getData());
      return;
    }
    m_smooth.assemble();
    m_smooth.solve(m_smoothKSP);
    m_sSmooth.setData(m_smoothTrial.getSolution().getData());
  }

  void PTAGXBench::reactionStep()
  {
    const auto& k = m_cfg.kinetics;
    const Real dt = m_cfg.dt;

    // One evaluation of the nodal step per DOF: the first projection computes
    // and stores the states (and the mass added by clipping, weighted by the
    // lumped mass), the others read them back in the same order (the MPI
    // interpolation visits each owned DOF once, by increasing global index).
    m_stateCache.clear();
    Real clipPTA = 0.0, clipG = 0.0, clipNMS = 0.0;
    const bool xa = k.factorXa;
    const bool monomer = k.immobileFibrin && k.fibrinMonomer;
    m_pNext.project(RealFunction([&](const Point& p) {
      const PlateletKinetics::State s{ m_p.getValue(p), m_t.getValue(p), m_a.getValue(p),
        m_g.getValue(p), m_f.getValue(p), m_n.getValue(p), m_m.getValue(p),
        m_s.getValue(p), xa ? m_x.getValue(p) : 0.0, monomer ? m_fm.getValue(p) : 0.0 };
      const Real w = m_lumped.getValue(p);
      const auto neg = [](Real c) { return std::max<Real>(-c, 0.0); };
      clipPTA += w * (neg(s.P) + neg(s.T) - neg(s.A));
      clipG += w * (neg(s.G) + neg(s.Fm));
      clipNMS += w * (neg(s.N) + neg(s.M) + neg(s.S));
      m_stateCache.push_back(k.advance(s, dt, m_sSmooth.getValue(p), surfaceThetaAt(p)));
      return m_stateCache.back().P;
    }));
    const auto& comm = m_mesh.getContext().getCommunicator();
    m_clipPTA += boost::mpi::all_reduce(comm, clipPTA, std::plus<Real>());
    m_clipG += boost::mpi::all_reduce(comm, clipG, std::plus<Real>());
    if (k.platelets)
      m_clipNMS += boost::mpi::all_reduce(comm, clipNMS, std::plus<Real>());

    const auto replay = [this](Real PlateletKinetics::State::* member) {
      m_cacheIndex = 0;
      return RealFunction([this, member](const Point&) {
        assert(m_cacheIndex < m_stateCache.size());
        return m_stateCache[m_cacheIndex++].*member;
      });
    };
    using St = PlateletKinetics::State;
    m_tNext.project(replay(&St::T));
    m_aNext.project(replay(&St::A));
    m_gNext.project(replay(&St::G));
    m_fNext.project(replay(&St::F));
    if (xa)
      m_xNext.project(replay(&St::X));
    if (monomer)
      m_fmNext.project(replay(&St::Fm));
    if (k.platelets)
    {
      m_nNext.project(replay(&St::N));
      m_mNext.project(replay(&St::M));
      m_sNext.project(replay(&St::S));
    }

    m_p.setData(m_pNext.getData());
    m_t.setData(m_tNext.getData());
    m_a.setData(m_aNext.getData());
    m_g.setData(m_gNext.getData());
    m_f.setData(m_fNext.getData());
    if (xa)
      m_x.setData(m_xNext.getData());
    if (monomer)
      m_fm.setData(m_fmNext.getData());
    if (k.platelets)
    {
      m_n.setData(m_nNext.getData());
      m_m.setData(m_mNext.getData());
      m_s.setData(m_sNext.getData());
    }
  }

  // ==========================================================================
  // Diagnostics
  // ==========================================================================
  PTAGXBench::Real PTAGXBench::integrate(const ScalarGridFunctionType& f)
  {
    m_form = Integral(f, m_q);
    m_form.assemble();
    return m_form(m_one);
  }

  void PTAGXBench::writeCSVHeader()
  {
    m_csv << "t,gamma,intP,intT,intA,intG,intF,intN,intM,intS,intSs,maxT,maxF,maxS,"
             "maxTheta,gelArea,gelHeight,depositArea,depositHeight,drift,minT,"
             "tmLoss,clipped,clipPTA,clipG,clipNMS,intX,intFm,maxX\n";
    m_csv.flush();
  }

  void PTAGXBench::writeCSVRow()
  {
    const auto& k = m_cfg.kinetics;
    const Real KF = k.gelThreshold;
    const Real Smax = k.plateletMax;
    const Real invariant = k.prothrombin0 - k.antithrombin0;

    m_stab.project(RealFunction([this, &k](const Point& p) {
      return k.stabilized(std::max<Real>(m_s.getValue(p), 0.0),
                          std::max<Real>(m_f.getValue(p), 0.0)); }));
    m_theta.project(RealFunction([this, &k](const Point& p) {
      return k.membrane(std::max<Real>(m_m.getValue(p), 0.0),
                        std::max<Real>(m_s.getValue(p), 0.0)); }));

    const Real intP = integrate(m_p), intT = integrate(m_t), intA = integrate(m_a),
               intG = integrate(m_g), intF = integrate(m_f), intN = integrate(m_n),
               intM = integrate(m_m), intS = integrate(m_s), intSs = integrate(m_stab);
    const Real maxT = m_t.max(), maxF = m_f.max(), maxS = m_s.max(),
               maxTheta = m_theta.max();
    const Real intX = k.factorXa ? integrate(m_x) : 0.0;
    const Real intFm = (k.immobileFibrin && k.fibrinMonomer) ? integrate(m_fm) : 0.0;
    const Real maxX = k.factorXa ? m_x.max() : 0.0;

    m_work.project(RealFunction([this, KF](const Point& p) -> Real {
      return m_f.getValue(p) >= KF ? 1.0 : 0.0; }));
    const Real gelArea = integrate(m_work);
    m_work.project(RealFunction([this, KF](const Point& p) -> Real {
      return m_f.getValue(p) >= KF ? p.y() : 0.0; }));
    const Real gelHeight = m_work.max();
    m_work.project(RealFunction([this, Smax](const Point& p) -> Real {
      return m_s.getValue(p) >= 0.1 * Smax ? 1.0 : 0.0; }));
    const Real depositArea = integrate(m_work);
    m_work.project(RealFunction([this, Smax](const Point& p) -> Real {
      return m_s.getValue(p) >= 0.1 * Smax ? p.y() : 0.0; }));
    const Real depositHeight = m_work.max();
    m_work.project(RealFunction([this, invariant](const Point& p) -> Real {
      return std::abs(m_p.getValue(p) + m_t.getValue(p) - m_a.getValue(p) - invariant); }));
    const Real drift = m_work.max();

    if (isRoot())
    {
      m_csv << m_time << ',' << m_gamma << ',' << intP << ',' << intT << ',' << intA << ','
            << intG << ',' << intF << ',' << intN << ',' << intM << ',' << intS << ','
            << intSs << ',' << maxT << ',' << maxF << ',' << maxS << ',' << maxTheta
            << ',' << gelArea << ',' << gelHeight << ',' << depositArea << ','
            << depositHeight << ',' << drift << ',' << m_minTransport << ',' << m_tmLoss
            << ',' << m_clipped << ',' << m_clipPTA << ',' << m_clipG << ',' << m_clipNMS
            << ',' << intX << ',' << intFm << ',' << maxX << '\n';
      m_csv.flush();
    }
  }

  int PTAGXBench::run()
  {
    if (!m_initialized)
      initialize();

    const auto& k = m_cfg.kinetics;
    const Real dt = m_cfg.dt;
    const int steps = static_cast<int>(std::lround(m_cfg.finalTime / dt));
    const int csvEvery = std::max(1, static_cast<int>(std::lround(m_cfg.csvInterval / dt)));
    const int outEvery = (m_cfg.outputInterval > 0.0)
      ? std::max(1, static_cast<int>(std::lround(m_cfg.outputInterval / dt))) : 0;

    if (isRoot())
    {
      Alert::Info() << "[run] " << steps << " steps of " << dt << " s" << Alert::Raise;
    }

    writeCSVRow();
    if (outEvery > 0)
      m_xdmf.write(m_time).flush();

    const bool tmActive = k.thrombomodulinRate > 0.0 && !m_cfg.thrombomodulin.empty();

    for (int step = 0; step < steps; ++step)
    {
      m_step = step;
      m_time += dt;
      updateVelocity(m_time);

      solveTransport(0, m_transportP, m_kspP, m_pT, m_p, m_pPrev);
      solveTransport(1, m_transportT, m_kspT, m_tT, m_t, m_tPrev);
      solveTransport(2, m_transportA, m_kspA, m_aT, m_a, m_aPrev);
      solveTransport(3, m_transportG, m_kspG, m_gT, m_g, m_gPrev);
      if (k.factorXa)
        solveTransport(6, m_transportX, m_kspX, m_xT, m_x, m_xPrev);
      if (k.immobileFibrin && k.fibrinMonomer)
        solveTransport(7, m_transportFm, m_kspFm, m_fmT, m_fm, m_fmPrev);
      m_minTransport = std::min<Real>(m_minTransport, m_t.min());

      if (tmActive)
      {
        m_form = BoundaryIntegral(m_t, m_q).over(m_cfg.thrombomodulin);
        m_form.assemble();
        m_tmLoss += dt * k.thrombomodulinRate * m_form(m_one);
      }

      using Clock = std::chrono::steady_clock;
      if (k.platelets)
      {
        solveTransport(4, m_transportN, m_kspN, m_nT, m_n, m_nPrev);
        solveTransport(5, m_transportM, m_kspM, m_mT, m_m, m_mPrev);
        const auto d0 = Clock::now();
        depositAtWall();
        smoothDeposit();
        m_tDeposit += std::chrono::duration<Real>(Clock::now() - d0).count();
      }

      const auto r0 = Clock::now();
      reactionStep();
      m_tReaction += std::chrono::duration<Real>(Clock::now() - r0).count();

      const auto g0 = Clock::now();
      const bool stop = m_cfg.stopThrombin > 0.0 && m_t.max() >= m_cfg.stopThrombin;
      if ((step + 1) % csvEvery == 0 || stop)
        writeCSVRow();
      m_tDiagnostics += std::chrono::duration<Real>(Clock::now() - g0).count();
      if (outEvery > 0 && (step + 1) % outEvery == 0)
        m_xdmf.write(m_time).flush();
      if (isRoot() && m_cfg.logEvery > 0 && (step + 1) % m_cfg.logEvery == 0)
      {
        Alert::Info() << "[step " << (step + 1) << "/" << steps << "] t=" << m_time
                      << " s  max T=" << (m_t.max() * 1e6) << " nM  max F/G0="
                      << (m_f.max() / k.fibrinogen0) << "  max S/Smax="
                      << (m_s.max() / k.plateletMax) << Alert::Raise;
      }
      if (stop)
      {
        if (isRoot())
        {
          Alert::Info() << "[run] max T reached " << (m_cfg.stopThrombin * 1e6)
                        << " nM at t=" << m_time << " s: stopping" << Alert::Raise;
        }
        break;
      }
    }

    m_xdmf.close();
    if (isRoot())
    {
      Alert::Info() << "[timing] assembly " << m_tAssembly << " s, solve " << m_tTransport
                    << " s, deposition " << m_tDeposit << " s, reaction " << m_tReaction
                    << " s, diagnostics " << m_tDiagnostics << " s" << Alert::Raise;
    }
    return 0;
  }
}

namespace
{
  /// @brief "1,3,5" -> {1, 3, 5}; "none" (or "") -> empty set.
  Rodin::FlatSet<Rodin::Geometry::Attribute> parseAttributes(const std::string& text)
  {
    Rodin::FlatSet<Rodin::Geometry::Attribute> out;
    std::stringstream ss(text);
    std::string item;
    while (std::getline(ss, item, ','))
    {
      if (item.empty() || item == "none" || item == "-")
        continue;
      out.insert(static_cast<Rodin::Geometry::Attribute>(std::stoul(item)));
    }
    return out;
  }
}

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

  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world(PETSC_COMM_WORLD, boost::mpi::comm_attach);
  Rodin::Context::MPI context(env, world);

  try
  {
    int status = 0;
    {
      using Rodin::Examples::Heart::Thrombosis::PTAGXBench;
      PTAGXBench::Config cfg;
      auto& k = cfg.kinetics;

      const auto getString = [](const char* name, std::string& out) {
        char buffer[512];
        PetscBool got = PETSC_FALSE;
        PetscOptionsGetString(PETSC_NULLPTR, PETSC_NULLPTR, name, buffer, sizeof(buffer), &got);
        if (got) out = buffer;
      };
      const auto getReal = [](const char* name, Rodin::Real& out) {
        PetscReal value = 0.0;
        PetscBool got = PETSC_FALSE;
        PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, name, &value, &got);
        if (got) out = value;
      };
      const auto getInt = [](const char* name, int& out) {
        PetscInt value = 0;
        PetscBool got = PETSC_FALSE;
        PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, name, &value, &got);
        if (got) out = static_cast<int>(value);
      };
      const auto getBool = [](const char* name, bool& out) {
        PetscBool value = PETSC_FALSE;
        PetscBool got = PETSC_FALSE;
        PetscOptionsGetBool(PETSC_NULLPTR, PETSC_NULLPTR, name, &value, &got);
        if (got) out = (value == PETSC_TRUE);
      };
      const auto getAttributes = [&getString](const char* name,
                                              Rodin::FlatSet<Rodin::Geometry::Attribute>& out) {
        std::string text;
        getString(name, text);
        if (!text.empty()) out = parseAttributes(text);
      };

      getString("-pg_mesh", cfg.meshPath);
      getReal("-pg_mesh_scale", cfg.meshScale);
      getString("-pg_xdmf", cfg.xdmfBasename);
      getString("-pg_csv", cfg.csvPath);
      getReal("-pg_dt", cfg.dt);
      getReal("-pg_t_end", cfg.finalTime);
      getReal("-pg_csv_interval", cfg.csvInterval);
      getReal("-pg_output_interval", cfg.outputInterval);
      getInt("-pg_log_every", cfg.logEvery);
      getAttributes("-pg_inlet", cfg.inlet);
      getAttributes("-pg_reactive", cfg.reactive);
      getAttributes("-pg_collagen", cfg.collagen);
      getAttributes("-pg_tm", cfg.thrombomodulin);
      getReal("-pg_activation", cfg.activation);
      getReal("-pg_H", cfg.channelHeight);
      getReal("-pg_gamma", cfg.wallShearRate);
      cfg.wallShearRate2 = cfg.wallShearRate;
      getReal("-pg_gamma2", cfg.wallShearRate2);
      getReal("-pg_t_switch", cfg.switchTime);
      getReal("-pg_crosswind", cfg.crosswindC);
      getReal("-pg_stop_thrombin", cfg.stopThrombin);

      getReal("-pg_wall_flux", k.wallFlux);
      getReal("-pg_k_tm", k.thrombomodulinRate);
      getReal("-pg_k_a", k.amplificationRate);
      getReal("-pg_P0", k.prothrombin0);
      getReal("-pg_A0", k.antithrombin0);
      getReal("-pg_G0", k.fibrinogen0);
      getReal("-pg_supg", k.stabilizationScale);
      getBool("-pg_immobile_fibrin", k.immobileFibrin);
      getReal("-pg_F_cap", k.fibrinCap);
      getReal("-pg_phi_i", k.intrinsicFraction);
      getReal("-pg_eps_x", k.xaCofactorFloor);
      getBool("-pg_fibrin_monomer", k.fibrinMonomer);
      getReal("-pg_k_pol", k.polymerizationRate);
      getBool("-pg_xa", k.factorXa);
      getReal("-pg_j_x", k.xaFlux);
      getReal("-pg_kappa_x", k.xaRate);
      getReal("-pg_k_dx", k.xaDecay);
      getBool("-pg_platelets", k.platelets);
      getBool("-pg_contact", k.contactActivation);
      getBool("-pg_RF_linear", k.linearFibrinogen);
      getReal("-pg_N_ref", k.membraneNorm);
      getReal("-pg_K_G", k.conversionHalf);
      if (!k.platelets)
        k.membraneBackground = 1.0;
      getReal("-pg_theta0", k.membraneBackground);
      getReal("-pg_sigma_s", k.surfaceAmplification);
      getReal("-pg_delta_s", k.surfaceLayer);
      getBool("-pg_axisym", cfg.axisymmetric);
      {
        PetscReal xr[2] = {0, 0};
        PetscInt n = 2;
        PetscBool got = PETSC_FALSE;
        PetscOptionsGetRealArray(PETSC_NULLPTR, PETSC_NULLPTR, "-pg_patch", xr, &n, &got);
        if (got && n == 2)
          cfg.patch = {{xr[0], xr[1]}};
      }
      getReal("-pg_N0", k.plateletCount);
      getReal("-pg_k_p", k.activationRate);
      getReal("-pg_T_p", k.activationHalf);
      getReal("-pg_S_max", k.plateletMax);
      getReal("-pg_k_w", k.adhesionActivated);
      getReal("-pg_k_wN", k.adhesionResting);
      getReal("-pg_k_agg", k.aggregationRate);
      getReal("-pg_l_agg", k.aggregationLength);
      getReal("-pg_D_pl0", k.plateletDiffusivity0);
      getReal("-pg_D_pl1", k.plateletDiffusivity1);

      PTAGXBench simulation(context, cfg);
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
