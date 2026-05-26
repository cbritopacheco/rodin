/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @example Three-layer level-set / graph-cut / TMOP pipeline on a translating circle.
 *
 * Layer 1 — Level set (analytic):
 *   phi(x,y,t) = sqrt((x - cx(t))^2 + (y - cy(t))^2) - R
 *   Interface Gamma(t) = { phi = 0 }, phi < 0 = Inside.
 *   Translation: cx(t) = cx0 + vx*t, cy(t) = cy0 + vy*t.
 *
 * Layer 2 — Variational classifier (seeded min s-t cut):
 *   Minimize: sum_K [chi_K D_K^In + (1-chi_K) D_K^Out]
 *           + lambda sum_{F=K|L} sigma_F |chi_K - chi_L|
 *   D_K^In  = int_K max(phi,0)^2 dx  (penalize positive phi inside)
 *   D_K^Out = int_K max(-phi,0)^2 dx (penalize negative phi outside)
 *   sigma_F = |F| * (epsilon + |phi_bar_F| / h_F)  (interface attraction)
 *   Seeds: cells with |phi_K| > 2*h_K are forced Inside or Outside.
 *   Interface skeleton: Gamma_h = { F = K+|K- : chi_K+ != chi_K- }.
 *
 * Layer 3 — High-order geometry fitting + TMOP quality:
 *   minimize_u  alpha E_TMOP(u) + rho*E_fit(u) + eta*||u||^2
 *   E_TMOP = sum_K int mu(A W^{-1}) det(W) dX   (quality + det barrier)
 *   E_fit  = int_{Gamma_h} phi(x+u)^2 ds        (smooth fit to phi=0)
 *   Minimizer: preconditioned gradient descent with Armijo backtracking +
 *              curved determinant validity gate.
 *
 * Pipeline separation:
 *   - phi defines continuous topology (never modified by TMOP).
 *   - Graph-cut classifies cells and builds the skeleton (topology layer).
 *   - TMOP makes the skeleton geometrically accurate (fitting layer).
 *   - TMOP does NOT decide topology. If the skeleton is wrong, report it.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include <Rodin/Adaptation.h>
#include <Rodin/Adaptation/TargetMatrixOptimization/IsoparametricGeometry.h>
#include <Rodin/Adaptation/TargetMatrixOptimization/LevelSetGraphCutClassifier.h>
#include <Rodin/Adaptation/TargetMatrixOptimization/Metrics.h>
#include <Rodin/Assembly/Default.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Adaptation::TargetMatrixOptimization;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  // ── Level-set configuration ──────────────────────────────────────────────

  constexpr Attribute Inside    = 1;
  constexpr Attribute Outside   = 10;
  constexpr Attribute Interface = 30;
  constexpr size_t   QOrder     = 4;

  constexpr Real CX0 = Real(0.35);
  constexpr Real CY0 = Real(0.50);
  constexpr Real VX  = Real(0.30);   // horizontal advection speed
  constexpr Real VY  = Real(0.00);
  constexpr Real R   = Real(0.15);

  using P2Element = RealH1Element<2>;

  inline Math::SpatialPoint circleCenter(Real t)
  {
    return { CX0 + VX * t, CY0 + VY * t };
  }

  Real phiAt(const Math::SpatialPoint& x, Real t)
  {
    const auto c = circleCenter(t);
    return std::hypot(x[0] - c[0], x[1] - c[1]) - R;
  }

  Math::SpatialPoint gradPhiAt(const Math::SpatialPoint& x, Real t)
  {
    const auto c = circleCenter(t);
    const Real dx = x[0] - c[0];
    const Real dy = x[1] - c[1];
    const Real rho = std::hypot(dx, dy);
    if (rho <= Real(1e-14))
      return { Real(1), Real(0) };
    return { dx / rho, dy / rho };
  }

  Math::SpatialMatrix<Real> hessianPhiAt(const Math::SpatialPoint& x, Real t)
  {
    Math::SpatialMatrix<Real> h(2, 2);
    h.setZero();
    const auto c = circleCenter(t);
    const Real dx = x[0] - c[0];
    const Real dy = x[1] - c[1];
    const Real rho2 = dx * dx + dy * dy;
    const Real rho  = std::sqrt(rho2);
    if (rho <= Real(1e-14))
      return h;
    const Real rho3 = rho2 * rho;
    h(0, 0) =  dy * dy / rho3;
    h(0, 1) = -dx * dy / rho3;
    h(1, 0) =  h(0, 1);
    h(1, 1) =  dx * dx / rho3;
    return h;
  }

  // ── Interface edge P2 sync ────────────────────────────────────────────────
  // Sets ParametricTransformation<P2Element> on every Interface-attributed edge
  // by mapping the edge's two endpoints and the midpoint through the adjacent
  // cell's P2 parametric transformation.  Needed so that the AnalyticLevelSetFitTerm
  // integrates on a curved interface rather than a piecewise-linear one.
  void syncInterfaceEdgeTransformations(
      LocalMesh& mesh,
      const P2Element& fe,
      Attribute interfaceAttribute)
  {
    auto& conn = mesh.getConnectivity();
    conn.compute(2, 1);
    conn.compute(1, 2);
    conn.compute(1, 0);

    P2Element feSeg(Polytope::Type::Segment);
    const auto& refSeg = P2Element::getNodes(Polytope::Type::Segment);

    static const std::array<Math::SpatialPoint, 3> triRef = {{
      Math::SpatialPoint{ Real(0), Real(0) },
      Math::SpatialPoint{ Real(1), Real(0) },
      Math::SpatialPoint{ Real(0), Real(1) }
    }};

    const size_t nEdges = conn.getCount(1);
    for (Index e = 0; e < static_cast<Index>(nEdges); ++e)
    {
      const auto attr = mesh.getAttribute(1, e);
      if (!attr || *attr != interfaceAttribute)
        continue;

      const auto& adjCells = conn.getIncidence({1, 2}, e);
      if (adjCells.empty())
        continue;

      const auto& edgeKey = conn.getPolytope(1, e);
      const auto  cellIdx = adjCells[0];
      const auto& cell    = conn.getPolytope(2, cellIdx);

      size_t localA = 3, localB = 3;
      for (size_t i = 0; i < 3; ++i)
      {
        if (cell(i) == edgeKey(0)) localA = i;
        if (cell(i) == edgeKey(1)) localB = i;
      }
      if (localA >= 3 || localB >= 3) continue;

      auto cellIt = mesh.getPolytope(2, cellIdx);
      PointCloud pc(2, refSeg.size());
      for (size_t j = 0; j < refSeg.size(); ++j)
      {
        const Real s = refSeg[j][0];
        const Math::SpatialPoint rc{
          (Real(1) - s) * triRef[localA][0] + s * triRef[localB][0],
          (Real(1) - s) * triRef[localA][1] + s * triRef[localB][1]
        };
        const auto x = Point(*cellIt, rc).getPhysicalCoordinates();
        pc(0, j) = x[0];
        pc(1, j) = x[1];
      }

      mesh.setPolytopeTransformation(
          { size_t(1), e },
          new ParametricTransformation<P2Element>(pc, feSeg));
    }
  }

  // ── Geometry diagnostics ─────────────────────────────────────────────────
  Real edgeMeasure(LocalMesh& mesh, Attribute attr)
  {
    Real m = Real(0);
    const auto& conn = mesh.getConnectivity();
    const auto& qf   = QF::PolytopeQuadratureFormula::get(QOrder, Polytope::Type::Segment);
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto a = mesh.getAttribute(1, e);
      if (!a || *a != attr) continue;
      auto it = mesh.getPolytope(1, e);
      for (size_t q = 0; q < qf.getSize(); ++q)
      {
        Math::SpatialMatrix<Real> J;
        it->getTransformation().jacobian(J, qf.getPoint(q));
        m += qf.getWeight(q) * std::hypot(J(0, 0), J(1, 0));
      }
    }
    return m;
  }

  Real cellMeasure(LocalMesh& mesh)
  {
    Real m = Real(0);
    const auto& conn = mesh.getConnectivity();
    const auto& qf   = QF::PolytopeQuadratureFormula::get(QOrder, Polytope::Type::Triangle);
    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
    {
      if (conn.getGeometry(2, c) != Polytope::Type::Triangle) continue;
      auto it = mesh.getPolytope(2, c);
      for (size_t q = 0; q < qf.getSize(); ++q)
      {
        Math::SpatialMatrix<Real> J;
        it->getTransformation().jacobian(J, qf.getPoint(q));
        m += qf.getWeight(q) * std::abs(J.determinant());
      }
    }
    return m;
  }

  // ── TMOP geometry step ────────────────────────────────────────────────────
  struct GeometryReport
  {
    bool   converged          = false;
    size_t iterations         = 0;
    Index  acceptedSteps      = 0;
    Index  rejectedSteps      = 0;
    Real   energyInitial      = 0;
    Real   energyFinal        = 0;
    Real   fitEnergyInitial   = 0;
    Real   fitEnergyFinal     = 0;
    Real   gradientInitial    = 0;
    Real   gradientFinal      = 0;
    Real   finalAlpha         = 0;
  };

  // Fit-refinement pass: runs AFTER the main TMOP pass (moveMesh already called).
  // Minimizes only the static fit + deviation: no quality term.
  // This corrects the residual phi error left by the quality-fit tradeoff in pass 1.
  Real fitRefinementPass(
      LocalMesh&       mesh,
      const P2Element& fe,
      Real             t,
      Real             h)
  {
    auto phiValue    = [t](const Math::SpatialPoint& x) { return phiAt(x, t); };
    auto phiGradient = [t](const Math::SpatialPoint& x) { return gradPhiAt(x, t); };
    auto phiHessian  = [t](const Math::SpatialPoint& x) { return hessianPhiAt(x, t); };

    // Sync P2 edges before computing ifaceMeasure and building terms.
    syncInterfaceEdgeTransformations(mesh, fe, Interface);

    const Real ifaceMeasure = std::max(edgeMeasure(mesh, Interface), Real(1e-12));

    // Very high fit weight: phi -> 0 is the only goal.
    // Moderate deviation anchor: prevents large shifts (pass 1 already committed most
    // of the movement; pass 2 only corrects the residual ~0.010 phi error).
    AnalyticLevelSetFitTerm fitRefine(
        phiValue, phiGradient, phiHessian,
        Optional<Attribute>(Interface),
        Real(5000));
    fitRefine.setQuadratureOrder(QOrder).setNormalization(ifaceMeasure);

    // devAnchor weight chosen so that a displacement of size ~0.01 (= fit_rms)
    // costs ~fit penalty.  With fit=5000/ifaceMeasure*0.01^2 and dev=20*0.01^2:
    // fit:dev ≈ 250 → fit dominates.
    DeviationTerm devAnchor(Real(20));

    VectorH1<2, LocalMesh> space(std::integral_constant<size_t, 2>{}, mesh, 2);
    GridFunction u2(space);
    u2.getData().setZero();
    TrialFunction p2(space);
    TestFunction  v2(space);

    auto makeResidual = [&]() {
      return fitRefine.residual(u2, v2) + devAnchor.residual(u2, v2);
    };
    auto energy = [&]() {
      return fitRefine.energy(u2) + devAnchor.energy(u2);
    };

    IsoparametricTMOPMinimizerParameters params;
    params.maxIterations        = 30;
    params.gradientTolerance    = Real(1e-10);
    params.stepTolerance        = Real(1e-10);
    params.energyTolerance      = Real(1e-14);
    params.preconditionerLength = std::max(h, Real(1e-8));  // shorter: local corrections
    params.minDetFloor          = Real(5e-4);

    // Use the curved validity gate (no quality target needed for pass 2).
    auto admissible = [&]() {
      return isCurvedMoveValid(mesh, u2, fe, params.minDetFloor);
    };

    minimizeIsoparametricTMOP(
        mesh, fe, u2, p2, v2, makeResidual, energy, admissible, Interface, params);

    // Sync edge transformations for the refined mesh state.
    syncInterfaceEdgeTransformations(mesh, fe, Interface);
    return fitRefine.energy(u2);
  }

  // dt > 0: advection-driven fit (normal velocity target = dt * V_n).
  // dt = 0: fall back to static analytic level-set fit (phi(x+u)^2).
  GeometryReport solveGeometryStep(
      LocalMesh&       mesh,
      const P2Element& fe,
      Real             t,
      Real             dt,
      Real             h,
      Real             qualityWeight,
      Real             fitWeight,
      Real             deviationWeight,
      Real             targetBlend,
      Index            maxIter)
  {
    auto phiValue    = [t](const Math::SpatialPoint& x) { return phiAt(x, t); };
    auto phiGradient = [t](const Math::SpatialPoint& x) { return gradPhiAt(x, t); };
    auto phiHessian  = [t](const Math::SpatialPoint& x) { return hessianPhiAt(x, t); };

    // Normal velocity: V_n(x) = V · (grad phi / |grad phi|).
    // target(x) = dt * V_n(x).
    auto advectionTarget = [t, dt](const Math::SpatialPoint& x) -> Real
    {
      const auto g   = gradPhiAt(x, t);
      const Real gn  = g.norm();
      if (gn <= Real(1e-14))
        return Real(0);
      // V = (VX, VY) = (0.3, 0)
      return dt * (VX * g[0] + VY * g[1]) / gn;
    };

    CurvedQualityTargetJacobian target(mesh, targetBlend);
    ShapeSizeBlendMetric         metric(Real(0.5));

    VectorH1<2, LocalMesh> space(std::integral_constant<size_t, 2>{}, mesh, 2);
    GridFunction u(space);
    u.getData().setZero();
    TrialFunction p(space);
    TestFunction  v(space);

    QualityTerm quality(metric, target, qualityWeight);
    quality.setQuadratureOrder(QOrder);

    DeviationTerm deviation(deviationWeight);

    const Real ifaceMeasure = std::max(edgeMeasure(mesh, Interface), Real(1e-12));

    GeometryReport rep;

    if (dt > Real(0))
    {
      // Pass 1 — combined predictor + corrector:
      //   E = α E_TMOP + ρ_adv/2 ∫(u·n − dt·V_n)² ds   (linear warm-start)
      //                + ρ_fit/2 ∫φ(x+u)² ds            (exact accuracy)
      //                + η/2 ∫|u|² dx

      InterfaceNormalAdvectionFitTerm advFit(
          advectionTarget,
          Optional<Attribute>(Interface),
          fitWeight * Real(0.5));
      advFit
        .setQuadratureOrder(QOrder)
        .setNormalization(ifaceMeasure);

      AnalyticLevelSetFitTerm staticFit(
          phiValue, phiGradient, phiHessian,
          Optional<Attribute>(Interface),
          fitWeight);
      staticFit
        .setQuadratureOrder(QOrder)
        .setNormalization(ifaceMeasure);

      auto makeResidual = [&]()
      {
        return quality.residual(u, v)
             + deviation.residual(u, v)
             + advFit.residual(u, v)
             + staticFit.residual(u, v);
      };
      auto energy = [&]()
      {
        return quality.energy(u)
             + deviation.energy(u)
             + advFit.energy(u)
             + staticFit.energy(u);
      };

      IsoparametricTMOPMinimizerParameters params;
      params.maxIterations         = std::max<Index>(1, maxIter);
      params.gradientTolerance     = Real(1e-10);
      params.stepTolerance         = Real(1e-10);
      params.energyTolerance       = Real(1e-14);
      params.preconditionerLength  = std::max(Real(2) * h, Real(1e-8));
      params.minDetFloor           = Real(5e-4);

      auto admissible = [&]()
      {
        return isTargetAdmissible(
            mesh, u, fe, target, params.minDetFloor,
            Polytope::Type::Triangle,
            static_cast<Index>(QOrder));
      };

      rep.fitEnergyInitial = staticFit.energy(u);

      const auto res = minimizeIsoparametricTMOP(
          mesh, fe, u, p, v, makeResidual, energy, admissible, Interface, params);

      rep.converged        = res.converged;
      rep.iterations       = res.iterations;
      rep.acceptedSteps    = res.acceptedSteps;
      rep.rejectedSteps    = res.rejectedSteps;
      rep.energyInitial    = res.initialEnergy;
      rep.energyFinal      = res.finalEnergy;
      rep.gradientInitial  = res.initialGradientNorm;
      rep.gradientFinal    = res.finalGradientNorm;
      rep.finalAlpha       = res.finalAlpha;
      rep.fitEnergyFinal   = staticFit.energy(u);
    }
    else
    {
      // Static fit: ρ/2 ∫ phi(x+u)² ds
      AnalyticLevelSetFitTerm fit(
          phiValue, phiGradient, phiHessian,
          Optional<Attribute>(Interface),
          fitWeight);
      fit
        .setQuadratureOrder(QOrder)
        .setNormalization(ifaceMeasure);

      auto makeResidual = [&]()
      {
        return quality.residual(u, v)
             + deviation.residual(u, v)
             + fit.residual(u, v);
      };
      auto energy = [&]()
      {
        return quality.energy(u)
             + deviation.energy(u)
             + fit.energy(u);
      };

      IsoparametricTMOPMinimizerParameters params;
      params.maxIterations         = std::max<Index>(1, maxIter);
      params.gradientTolerance     = Real(1e-10);
      params.stepTolerance         = Real(1e-10);
      params.energyTolerance       = Real(1e-14);
      params.preconditionerLength  = std::max(Real(2) * h, Real(1e-8));
      params.minDetFloor           = Real(5e-4);

      auto admissible = [&]()
      {
        return isTargetAdmissible(
            mesh, u, fe, target, params.minDetFloor,
            Polytope::Type::Triangle,
            static_cast<Index>(QOrder));
      };

      rep.fitEnergyInitial = fit.energy(u);

      const auto res = minimizeIsoparametricTMOP(
          mesh, fe, u, p, v, makeResidual, energy, admissible, Interface, params);

      rep.converged        = res.converged;
      rep.iterations       = res.iterations;
      rep.acceptedSteps    = res.acceptedSteps;
      rep.rejectedSteps    = res.rejectedSteps;
      rep.energyInitial    = res.initialEnergy;
      rep.energyFinal      = res.finalEnergy;
      rep.gradientInitial  = res.initialGradientNorm;
      rep.gradientFinal    = res.finalGradientNorm;
      rep.finalAlpha       = res.finalAlpha;
      rep.fitEnergyFinal   = fit.energy(u);
    }

    // Pass 2 — fit-only refinement (corrects the residual phi error left by
    // the quality-fit tradeoff in pass 1; no quality term; small corrections).
    fitRefinementPass(mesh, fe, t, h);

    return rep;
  }

} // anonymous namespace

int main(int argc, char** argv)
{
  // ── Parameters ─────────────────────────────────────────────────────────────
  size_t resolution   = 20;
  Index  steps        = 20;
  Real   fitWeight    = Real(100);   // high fit drive
  Real   qualWeight   = Real(8);     // fit:quality = 12.5:1; quality actually matters
  Real   devWeight    = Real(5);     // h-scaled regularization
  Real   targetBlend  = Real(0.2);   // 80% equilateral target; resists element collapse
  Index  maxIter      = 400;
  Real   gcPerimeter  = Real(1);
  Real   gcEpsilon    = Real(1e-3);

  if (argc > 1) resolution  = static_cast<size_t>(std::max(3, std::atoi(argv[1])));
  if (argc > 2) steps       = static_cast<Index>(std::max(2, std::atoi(argv[2])));
  if (argc > 3) fitWeight   = static_cast<Real>(std::atof(argv[3]));
  if (argc > 4) qualWeight  = static_cast<Real>(std::atof(argv[4]));
  if (argc > 5) devWeight   = static_cast<Real>(std::atof(argv[5]));
  if (argc > 6) targetBlend = static_cast<Real>(std::atof(argv[6]));
  if (argc > 7) maxIter     = static_cast<Index>(std::max(1, std::atoi(argv[7])));
  if (argc > 8) gcPerimeter = static_cast<Real>(std::atof(argv[8]));
  if (argc > 9) gcEpsilon   = static_cast<Real>(std::atof(argv[9]));

  const Real h = Real(1) / static_cast<Real>(resolution - 1);

  // ── Mesh setup ──────────────────────────────────────────────────────────────
  LocalMesh mesh =
    LocalMesh::UniformGrid(Polytope::Type::Triangle, { resolution, resolution });
  mesh.scale(h);
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 2);
  mesh.getConnectivity().compute(1, 0);

  // Upgrade to P2 (isoparametric high-order geometry)
  P2Element fe(Polytope::Type::Triangle);
  upgradeTransformations(mesh, fe);

  // ── Graph-cut classifier options ────────────────────────────────────────────
  LevelSetGraphCutClassifier classifier;
  LevelSetGraphCutClassifier::Options gcOpts;
  gcOpts.insideAttribute    = Inside;
  gcOpts.outsideAttribute   = Outside;
  gcOpts.interfaceAttribute = Interface;
  gcOpts.perimeterWeight    = gcPerimeter;
  gcOpts.attractionEpsilon  = gcEpsilon;
  gcOpts.seedDistanceFactor = Real(2);
  gcOpts.quadratureOrder    = QOrder;

  // ── XDMF output ─────────────────────────────────────────────────────────────
  IO::XDMF xdmf("AdvectedGraphCutCurvedInterface");
  auto grid = xdmf.grid("p2-graphcut");

  std::cout << std::setprecision(12);
  std::cout
    << "step t "
    << "gc_inside gc_outside gc_interface "
    << "gc_seeds_in gc_seeds_out gc_optcells "
    << "gc_inside_comp gc_outside_comp gc_nonmanifold "
    << "gc_phi_rms gc_phi_max "
    << "gc_data_energy gc_perimeter_energy "
    << "tmop_converged tmop_iter tmop_accepted tmop_rejected "
    << "tmop_energy_i tmop_energy_f "
    << "tmop_grad_i tmop_grad_f tmop_alpha "
    << "fit_energy_i fit_energy_f "
    << "curved_fit_rms curved_fit_max curved_qmin curved_min_det "
    << "curved_invalid curved_overlap "
    << "area\n";

  for (Index step = 0; step < steps; ++step)
  {
    const Real t = static_cast<Real>(step)
        / static_cast<Real>(std::max<Index>(1, steps - 1));

    // ── Layer 1: Level set (analytic phi) ─────────────────────────────────
    // phi is evaluated on-demand by the classifier and TMOP terms; no
    // discrete field is stored on the mesh.  For XDMF visualization we
    // project onto P1.

    // ── Layer 2: Variational graph-cut classifier ─────────────────────────
    const auto gcResult = classifier.classify(
        mesh,
        [](const Math::SpatialPoint& x, Real time) { return phiAt(x, time); },
        t,
        gcOpts);

    // Upgrade interface edge transformations to P2 (needed for curved fit)
    syncInterfaceEdgeTransformations(mesh, fe, Interface);

    // ── Layer 3: High-order geometry fitting + TMOP quality ──────────────
    // Step 0: static phi fit (initial placement).
    // Step k>0: advection-driven fit (u·n = dt * V_n) + static corrector.
    // Combined single-pass objective:
    //   α E_TMOP  +  ρ_adv E_adv  +  ρ_fit E_fit  +  η E_dev
    // TMOP quality (qualWeight=8) and fit (fitWeight=100) compete at 12.5:1.
    // This ratio prevents element collapse while maintaining excellent fit.
    const Real dt = (step > 0)
        ? Real(1) / static_cast<Real>(std::max<Index>(1, steps - 1))
        : Real(0);
    const auto geoReport = solveGeometryStep(
        mesh, fe, t, dt, h,
        qualWeight, fitWeight, devWeight, targetBlend, maxIter);

    // Re-classify after geometry update so diagnostics reflect moved mesh
    const auto gcPost = classifier.classify(
        mesh,
        [](const Math::SpatialPoint& x, Real time) { return phiAt(x, time); },
        t,
        gcOpts);
    syncInterfaceEdgeTransformations(mesh, fe, Interface);

    // Curved diagnostics
    const auto phiValue = [t](const Math::SpatialPoint& x) { return phiAt(x, t); };
    const auto curved   = curvedMetrics(mesh, phiValue, Interface);
    const auto tqm      = targetQualityMetrics(
        mesh,
        ShapeSizeBlendMetric(Real(0.5)),
        CurvedQualityTargetJacobian(mesh, targetBlend),
        Polytope::Type::Triangle,
        static_cast<Index>(QOrder));

    // ── XDMF output ───────────────────────────────────────────────────────
    {
      P1<Real, LocalMesh> p1(mesh);
      GridFunction<P1<Real, LocalMesh>, Math::Vector<Real>> phi(p1);
      RealFunction phiFn([t](const Point& p)
      {
        return phiAt(p.getPhysicalCoordinates(), t);
      });
      phi.project(phiFn);

      grid.setMesh(mesh, IO::XDMF::MeshPolicy::Transient);
      grid.clear();
      grid.add("phi",  phi, IO::XDMF::Center::Node);
      xdmf.write(t).flush();
    }

    // ── Console report ────────────────────────────────────────────────────
    std::cout
      << step   << ' ' << t << ' '
      << gcResult.insideCells    << ' '
      << gcResult.outsideCells   << ' '
      << gcResult.interfaceFacets << ' '
      << gcResult.forcedInsideSeeds  << ' '
      << gcResult.forcedOutsideSeeds << ' '
      << gcResult.optimizedCells     << ' '
      << gcResult.connectedInsideComponents  << ' '
      << gcResult.connectedOutsideComponents << ' '
      << (gcResult.hasNonmanifoldInterface ? 1 : 0) << ' '
      << gcResult.interfacePhiRms << ' '
      << gcResult.interfacePhiMax << ' '
      << gcResult.dataEnergy      << ' '
      << gcResult.perimeterEnergy << ' '
      << (geoReport.converged ? 1 : 0) << ' '
      << geoReport.iterations   << ' '
      << geoReport.acceptedSteps << ' '
      << geoReport.rejectedSteps << ' '
      << geoReport.energyInitial << ' '
      << geoReport.energyFinal   << ' '
      << geoReport.gradientInitial << ' '
      << geoReport.gradientFinal   << ' '
      << geoReport.finalAlpha      << ' '
      << geoReport.fitEnergyInitial << ' '
      << geoReport.fitEnergyFinal   << ' '
      << curved.fitRms  << ' '
      << curved.fitMax  << ' '
      << curved.qmin    << ' '
      << curved.minDet  << ' '
      << curved.invalidJacobianSamples << ' '
      << curved.overlapSamples         << ' '
      << cellMeasure(mesh)
      << '\n';

    // Reset P2 for next step: commit moved corners to the linear backbone,
    // then reset midpoints to geometric midpoints.  This prevents accumulated
    // P2-midpoint distortion across steps, which would otherwise drive some
    // element Jacobians to near-zero and collapse the solver.
    syncLinearBackbone(mesh, fe);
    upgradeTransformations(mesh, fe);
  }

  xdmf.close();
  return 0;
}
