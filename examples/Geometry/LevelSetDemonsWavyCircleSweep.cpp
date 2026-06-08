/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//
// Wavy-circle sweep test for a demons-style level-set mesh motion.
//
// At every frame the background grid is classified from the analytic level
// set. The classified skeleton is converted to a signed-distance
// representative psi by FMM. The mesh displacement is then built by repeated
// robust demons updates:
//
//   r_Gamma,k(X) = phi(X + u_k(X)) / |grad(phi)(X + u_k(X))|,
//   X in Gamma_h,
//
// followed by the linearized extension
//
//   int du . v + ell^2 a_el(du, v)
//     + gamma_Gamma int_Gamma (n_phi . du)(n_phi . v)
//     = - gamma_Gamma int_Gamma r_Gamma (n_phi . v),
//
// with homogeneous Dirichlet data on the box boundary. Thus psi only localizes
// the classified skeleton for output; the interface integrals are evaluated on
// the actual classified interior faces. The accepted step uses the sampled
// quadratic admissibility predictor. There is no energy line search and
// therefore no backtracking; failure is reported as a best-effort admissible
// state.
//
#include <Rodin/Adaptation.h>
#include <Rodin/Assembly.h>
#include <Rodin/Distance/Eikonal.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Math.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>
#include <Rodin/Solver/CG.h>
#include <Rodin/Solver/SparseLU.h>
#include <Rodin/Solid.h>
#include <Rodin/Variational.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Adaptation;

namespace
{
  using Vec2 = Math::SpatialVector<Real>;
  using LocalMesh = Geometry::Mesh<Context::Local>;

  Vec2 vec2(Real x = 0, Real y = 0)
  {
    Vec2 out(2);
    out(0) = x;
    out(1) = y;
    return out;
  }

  template <class Displacement>
  void updateMovedMeshFromDisplacement(
      const LocalMesh& mesh,
      LocalMesh& moved,
      const Displacement& u)
  {
    const auto& uFes = u.getFiniteElementSpace();
    const auto& uData = u.getData();
    const Index vn = mesh.getVertexCount();
    for (Index vertex = 0; vertex < vn; ++vertex)
    {
      const Vec2 x = mesh.getVertexCoordinates(vertex);
      const auto& dofs = uFes.getDOFs(0, vertex);
      moved.setVertexCoordinates(
          vertex,
          vec2(x(0) + uData(dofs[0]), x(1) + uData(dofs[1])));
    }
  }

  struct AdmissibilityDistribution
  {
    Real minJ = std::numeric_limits<Real>::infinity();
    Real p10J = 0, p50J = 0, p90J = 0;
    Real maxQ = 0;
    Real p10Q = 0, p50Q = 0, p90Q = 0;
    std::size_t numQuadPoints = 0;
    std::size_t numNearJSafe = 0;   // j_K < 1.5·j_safe
    std::size_t numNearQMax = 0;    // Q_rel > 0.7·Q_max
    std::size_t numInadmissible = 0;
  };

  /// Sampled per-cell admissibility distribution (j, Q_rel) at all
  /// quadrature points of the displacement field. Returns percentiles
  /// (10/50/90) of both margins plus counts of cells "near" each
  /// constraint boundary. Lighter than evaluateLSRAdmissibilitySampled
  /// only conceptually — same dominant cost (cell-loop gradU eval),
  /// but accumulates a distribution rather than min/max.
  template <class Displacement>
  AdmissibilityDistribution
  evaluateAdmissibilityDistribution(
      Displacement& u,
      Real jSafe,
      Real qMax,
      std::size_t qOrder = 2)
  {
    using Variational::IntegrationPoint;
    using Variational::Jacobian;

    AdmissibilityDistribution dist;
    const auto& fes = u.getFiniteElementSpace();
    const auto& mesh = fes.getMesh();
    const std::size_t dim = mesh.getDimension();
    auto gradU = Jacobian(u);

    std::vector<Real> jSamples;
    std::vector<Real> qSamples;
    jSamples.reserve(mesh.getCellCount() * 4);
    qSamples.reserve(mesh.getCellCount() * 4);

    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const auto& cell = *cellIt;
      const auto& fe =
        fes.getFiniteElement(cell.getDimension(), cell.getIndex());
      const auto& qf =
        QF::PolytopeQuadratureFormula::get(
            qOrder, cell.getGeometry());
      const auto& quadrature = cell.getQuadrature(qf);
      for (std::size_t q = 0; q < quadrature.getSize(); ++q)
      {
        const auto& pt = quadrature.getPoint(q);
        const IntegrationPoint ip(pt, &qf, q);
        Math::SpatialMatrix<Real> F =
          Math::SpatialMatrix<Real>::Identity(dim, dim)
          + gradU.getValue(ip);
        const Real j = F.determinant();
        Real qRel = std::numeric_limits<Real>::infinity();
        if (j > Real(0))
          qRel = F.squaredNorm()
            / (static_cast<Real>(dim)
               * std::pow(j, Real(2) / static_cast<Real>(dim)));
        jSamples.push_back(j);
        qSamples.push_back(qRel);
        if (j < dist.minJ) dist.minJ = j;
        if (qRel > dist.maxQ) dist.maxQ = qRel;
        if (j <= jSafe) ++dist.numInadmissible;
        if (j < Real(1.5) * jSafe) ++dist.numNearJSafe;
        if (qRel > Real(0.7) * qMax) ++dist.numNearQMax;
      }
    }
    dist.numQuadPoints = jSamples.size();
    if (!jSamples.empty())
    {
      auto pct = [](std::vector<Real>& v, double p) -> Real {
        const std::size_t k = static_cast<std::size_t>(
            p * static_cast<double>(v.size()));
        std::nth_element(v.begin(), v.begin() + k, v.end());
        return v[k];
      };
      dist.p10J = pct(jSamples, 0.10);
      dist.p50J = pct(jSamples, 0.50);
      dist.p90J = pct(jSamples, 0.90);
      dist.p10Q = pct(qSamples, 0.10);
      dist.p50Q = pct(qSamples, 0.50);
      dist.p90Q = pct(qSamples, 0.90);
    }
    return dist;
  }

  struct WavyCircleLevelSet
  {
    Real cx = Real(0.5);
    Real cy = Real(0.5);
    Real R0 = Real(0.20);
    Real amp = Real(0.05);
    Real k = Real(6);
    Real phase = Real(0);

    Real phi(const Vec2& p) const
    {
      const Real dx = p(0) - cx;
      const Real dy = p(1) - cy;
      const Real r = std::sqrt(dx * dx + dy * dy);
      const Real theta = std::atan2(dy, dx);
      return r - (R0 + amp * std::cos(k * theta + phase));
    }

    Vec2 grad(const Vec2& p) const
    {
      const Real dx = p(0) - cx;
      const Real dy = p(1) - cy;
      const Real r2 = dx * dx + dy * dy;
      const Real r = std::max(std::sqrt(r2), Real(1e-14));
      const Real theta = std::atan2(dy, dx);
      const Real dRdtheta = -amp * k * std::sin(k * theta + phase);
      const Real r2safe = std::max(r2, Real(1e-28));
      return vec2(
          dx / r - dRdtheta * (-dy / r2safe),
          dy / r - dRdtheta * ( dx / r2safe));
    }
  };

  constexpr std::array<std::array<Real, 3>, 3> TriangleBarycentricQuadrature = {{
    {{ Real(2) / 3, Real(1) / 6, Real(1) / 6 }},
    {{ Real(1) / 6, Real(2) / 3, Real(1) / 6 }},
    {{ Real(1) / 6, Real(1) / 6, Real(2) / 3 }}
  }};

  Real applyPhaseMomentMap(Real phi, Real epsilon)
  {
    return std::tanh(phi / epsilon);
  }

  Vec2 interpolateVec(
      const std::array<Vec2, 3>& values,
      const std::array<Real, 3>& bary)
  {
    return bary[0] * values[0] + bary[1] * values[1] + bary[2] * values[2];
  }

  struct CellMomentInfo
  {
    Index index = 0;
    Real area = 0;
    Real moment = 0;
    std::array<Vec2, 3> x;
  };

  template <class PhiFn>
  std::vector<CellMomentInfo> collectCellMomentInfo(
      const LocalMesh& mesh, PhiFn&& phi, Real epsilon)
  {
    std::vector<CellMomentInfo> cells;
    cells.reserve(mesh.getCellCount());
    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const auto& cell = *cellIt;
      const auto& vertices = cell.getVertices();
      if (vertices.size() != 3)
        throw std::runtime_error(
            "LevelSetDemonsWavyCircleSweep expects triangular cells.");

      CellMomentInfo info;
      info.index = cell.getIndex();
      for (std::size_t i = 0; i < 3; ++i)
        info.x[i] = mesh.getVertexCoordinates(vertices[i]);

      const Vec2 e1 = info.x[1] - info.x[0];
      const Vec2 e2 = info.x[2] - info.x[0];
      info.area =
        std::abs(Real(0.5) * (e1(0) * e2(1) - e1(1) * e2(0)));

      Real moment = 0;
      for (const auto& bary : TriangleBarycentricQuadrature)
        moment += applyPhaseMomentMap(phi(interpolateVec(info.x, bary)), epsilon);
      info.moment = moment / TriangleBarycentricQuadrature.size();
      cells.push_back(std::move(info));
    }
    return cells;
  }

  Real facetLength(const LocalMesh& mesh, Index facet)
  {
    const auto face = mesh.getFace(facet);
    const auto& vertices = face->getVertices();
    const Vec2 a = mesh.getVertexCoordinates(vertices[0]);
    const Vec2 b = mesh.getVertexCoordinates(vertices[1]);
    return (b - a).norm();
  }

  void clearXDMFRegionAttributes(LocalMesh& mesh)
  {
    const std::size_t D = mesh.getDimension();
    for (Index i = 0; i < static_cast<Index>(mesh.getCellCount()); ++i)
      mesh.setAttribute({D, i}, Attribute{0});
    for (auto it = mesh.getFace(); it; ++it)
      mesh.setAttribute({D - 1, it->getIndex()}, Attribute{0});
  }

  std::size_t parseSizeTOption(
      int argc, char** argv, const std::string& name, std::size_t fallback)
  {
    const std::string prefix = "--" + name + "=";
    for (int i = 1; i < argc; ++i)
    {
      const std::string arg(argv[i]);
      if (arg.rfind(prefix, 0) == 0)
        return static_cast<std::size_t>(std::stoul(arg.substr(prefix.size())));
    }
    return fallback;
  }

  Real parseRealOption(
      int argc, char** argv, const std::string& name, Real fallback)
  {
    const std::string prefix = "--" + name + "=";
    for (int i = 1; i < argc; ++i)
    {
      const std::string arg(argv[i]);
      if (arg.rfind(prefix, 0) == 0)
        return static_cast<Real>(std::stod(arg.substr(prefix.size())));
    }
    return fallback;
  }

  bool hasFlag(int argc, char** argv, const std::string& name)
  {
    const std::string flag = "--" + name;
    for (int i = 1; i < argc; ++i)
      if (std::string(argv[i]) == flag)
        return true;
    return false;
  }

}

int main(int argc, char** argv)
{
  const std::size_t n =
    parseSizeTOption(argc, argv, "n", 50);
  const std::size_t nFrames =
    parseSizeTOption(argc, argv, "frames", 40);

  const Real orbitR = parseRealOption(argc, argv, "orbitR", Real(0.10));
  const Real amp = parseRealOption(argc, argv, "amp", Real(0.05));
  const Real R0 = parseRealOption(argc, argv, "R0", Real(0.20));
  const Real kLobes = parseRealOption(argc, argv, "lobes", Real(6));

  const Real h = Real(1) / static_cast<Real>(n - 1);
  const Real epsilon = parseRealOption(argc, argv, "classifier-eps", Real(1.25) * h);
  const Real lambdaC = parseRealOption(argc, argv, "classifier-lambda", Real(0.008));
  const std::size_t maxIterations =
    parseSizeTOption(argc, argv, "demon-steps", 25);
  const Real fitTol =
    parseRealOption(argc, argv, "fit-tol", Real(4) * h * h);
  const Real stepTol =
    parseRealOption(argc, argv, "step-tol", Real(1e-4) * h);
  const Real extensionEll =
    parseRealOption(argc, argv, "extension-ell", Real(3) * h);
  const Real extensionLambda =
    parseRealOption(argc, argv, "extension-lambda", Real(0));
  const Real extensionMu =
    parseRealOption(argc, argv, "extension-mu", Real(1));
  // Interface weight γ_Γ defaults to min(cap, 1/h_ref) — same scaling
  // rule as the LSR examples. With γ_Γ = 1 (the previous default) the
  // bulk L² + Helmholtz-elasticity terms dominate the surface
  // coupling and δu collapses to 0 on under-resolved features.
  constexpr Real kInterfaceWeightCap = Real(50);
  const Real interfaceWeight =
    parseRealOption(argc, argv, "interface-weight",
                    std::min(kInterfaceWeightCap, Real(1) / h));

  // Robust loss σ on the interface residual r_Γ = φ/|∇φ|. Selected
  // by --interface-welsch (Welsch) or default Charbonnier. Welsch
  // fully drops topology-mismatched skeleton points (w = exp(-r²/σ²));
  // Charbonnier caps their contribution (w = 1/√(1+r²/σ²)).
  const Real charbonnierSigma =
    parseRealOption(argc, argv, "charbonnier-sigma", Real(5) * h);
  const bool kInterfaceWelsch =
    hasFlag(argc, argv, "interface-welsch");

  // Direction-as-Dirichlet mode (Step A/B/D architecture):
  //   A : compute d(X) closed-form on Γ_h (no solve).
  //   B : solve barriers + Q_shape Newton step subject to
  //       δu|_Γ_h = d (hard Dirichlet trace), δu|_∂Ω = 0.
  //   D : admissibility-clipped step.
  // Drops the InterfaceTangent + InterfaceResidual penalty in favour
  // of hard trace data. Requires --quality-extension.
  const bool kDirectionDirichlet =
    hasFlag(argc, argv, "direction-dirichlet");
  const Real trustRadius =
    parseRealOption(argc, argv, "trust-radius", Real(0.2) * h);
  const Real jMinRatio = parseRealOption(argc, argv, "j-min", Real(1e-8));
  const Real jSafeRatio = parseRealOption(argc, argv, "j-safe", Real(1e-3));

  // Q3 mode: replace L² mass + linear elasticity with the linearised
  // LSR barriers (shape Q_rel barrier + j-barrier + j-volume-tether).
  // This makes the extension Q_shape-anisotropic and cell-quality-aware
  // — matching the user's proposed quality extension.
  //   --quality-extension : enable barrier-based extension.
  //   --gamma-shape       : shape barrier weight (default sqrt(h)).
  //   --qbar-act, --qbar-max : Q_rel barrier activation/asymptote.
  //   --jbar-weight, --jbar-safe-ratio : j barrier knobs.
  //   --jvol-weight       : centred (log j)^2 tether weight.
  const bool kQualityExtension =
    hasFlag(argc, argv, "quality-extension");
  const Real kGammaShape =
    parseRealOption(argc, argv, "gamma-shape", std::sqrt(h));
  const Real kQBarrierWeight =
    parseRealOption(argc, argv, "qbar-weight", Real(0));
  const Real kQBarrierAct =
    parseRealOption(argc, argv, "qbar-act",
                    std::numeric_limits<Real>::infinity());
  const Real kQBarrierMax =
    parseRealOption(argc, argv, "qbar-max",
                    std::numeric_limits<Real>::infinity());
  const Real kJBarrierWeight =
    parseRealOption(argc, argv, "jbar-weight", Real(0));
  const Real kJBarrierSafeRatio =
    parseRealOption(argc, argv, "jbar-safe-ratio", Real(1e-3));
  const Real kJVolTetherWeight =
    parseRealOption(argc, argv, "jvol-weight", Real(0));

  // (Two-step mode was tested and underperformed the single-solve
  // version: decoupling surface direction from barrier projection
  // forces information to flow through an intermediate field with
  // its own smoother, costing fit quality. Single-solve with
  // --quality-extension keeps direction, surface, and barriers
  // mutually consistent.)
  const bool kTwoStep = hasFlag(argc, argv, "two-step");
  const Real kDirectionEll =
    parseRealOption(argc, argv, "direction-ell", Real(2) * h);

  // Volumetric push: add the linearised LSR Gauss–Newton push tangent
  // to the demons single-solve bilinear form. This bypasses the
  // classifier-resolution fit floor inherent to surface-only coupling
  // by reading φ volumetrically through the analytic callback.
  //   E_push(u) = (ρ_S/2) · N · ∫_Ω W(ψ)·(φ(X+u) − ψ)² dX,
  //   W(s) = exp(−s²/(2 δ²)).
  const bool kVolumetricPush =
    hasFlag(argc, argv, "volumetric-push");
  const Real kPushWeight =
    parseRealOption(argc, argv, "push-weight", Real(1));
  const Real kPushDeltaW =
    parseRealOption(argc, argv, "push-delta-w", Real(2.6) * h);
  // σ for Welsch on the push residual. 0 disables (quadratic loss).
  // Default 1.5·h matches the LSR merge example's robust setting.
  const Real kPushLossSigma =
    parseRealOption(argc, argv, "push-loss-sigma", Real(1.5) * h);
  const Real jLineSearchRatio =
    parseRealOption(argc, argv, "j-ls", std::max(jMinRatio, Real(10) * jSafeRatio));
  const Real alphaSafety = parseRealOption(argc, argv, "alpha-safety", Real(0.9));
  const std::size_t qOrder = parseSizeTOption(argc, argv, "quad-order", 2);
  const bool verbose = hasFlag(argc, argv, "verbose");
  const bool trace = hasFlag(argc, argv, "trace");

  // Stall detector. The iteration is declared stalled when the
  // best-so-far fit has not improved by more than `stallTol` (relative
  // to itself) over the last `stallWindow` accepted steps. This
  // catches "iteration is still admissible and ticking but not
  // actually moving the fit" — the usual failure mode for demons at
  // the chord-error floor or at a topology-mismatch trap.
  const std::size_t stallWindow =
    parseSizeTOption(argc, argv, "stall-window", 5);
  const Real stallTol =
    parseRealOption(argc, argv, "stall-tol", Real(1e-3));

  constexpr Attribute interiorAttribute = 1;
  constexpr Attribute exteriorAttribute = 2;
  constexpr Attribute interfaceAttribute = 10;
  constexpr Attribute boundaryAttribute = 20;

  LocalMesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { n, n });
  mesh.scale(h);
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 2);
  mesh.getConnectivity().compute(2, 2);
  mesh.getConnectivity().compute(0, 0);

  for (auto faceIt = mesh.getBoundary(); faceIt; ++faceIt)
    mesh.setAttribute({mesh.getDimension() - 1, faceIt->getIndex()},
                      boundaryAttribute);

  using ScalarP1 = P1<Real, LocalMesh>;
  using ScalarP0 = P0<Real, LocalMesh>;
  using VectorP1 = P1<Math::SpatialVector<Real>, LocalMesh>;

  ScalarP1 p1Fes(mesh);
  ScalarP0 p0Fes(mesh);
  VectorP1 vectorFes(mesh, 2);

  GridFunction psi(p1Fes);
  psi.setName("psi");
  GridFunction phiGf(p1Fes);
  phiGf.setName("phi");
  GridFunction cellLabel(p0Fes);
  cellLabel.setName("cell_label");
  GridFunction phaseMoment(p0Fes);
  phaseMoment.setName("phase_moment");
  GridFunction u(vectorFes);
  u.setName("displacement");
  GridFunction du(vectorFes);
  du.setName("demon_step");

  LocalMesh moved(mesh);
  ScalarP0 p0FesMoved(moved);
  ScalarP1 p1FesMoved(moved);
  GridFunction movedLabel(p0FesMoved);
  movedLabel.setName("cell_label");
  GridFunction phiMoved(p1FesMoved);
  phiMoved.setName("phi_moved");
  GridFunction jMoved(p0FesMoved);
  jMoved.setName("j");
  GridFunction qRelMoved(p0FesMoved);
  qRelMoved.setName("q_rel");

  IO::XDMF xdmf("LevelSetDemonsWavyCircleSweep");
  auto psiGrid = xdmf.grid("psi");
  psiGrid.setMesh(mesh, IO::XDMF::MeshPolicy::Transient);
  psiGrid.add(cellLabel, IO::XDMF::Center::Cell);
  psiGrid.add(phaseMoment, IO::XDMF::Center::Cell);
  psiGrid.add(phiGf, IO::XDMF::Center::Node);
  psiGrid.add(psi, IO::XDMF::Center::Node);
  psiGrid.add(u, IO::XDMF::Center::Node);
  psiGrid.add(du, IO::XDMF::Center::Node);

  auto movedGrid = xdmf.grid("moved");
  movedGrid.setMesh(moved, IO::XDMF::MeshPolicy::Transient);
  movedGrid.add(movedLabel, IO::XDMF::Center::Cell);
  movedGrid.add(jMoved, IO::XDMF::Center::Cell);
  movedGrid.add(qRelMoved, IO::XDMF::Center::Cell);
  movedGrid.add(phiMoved, IO::XDMF::Center::Node);

  std::cout << "Wavy-circle demons sweep on " << n << "x" << n
            << " unit-square mesh, " << nFrames << " frames\n";
  std::cout << "  R0=" << R0 << "  amp=" << amp << "  k=" << kLobes
            << "  orbit R=" << orbitR
            << "  interfaceWeight=" << interfaceWeight
            << "  extensionEll=" << extensionEll << '\n';

  std::size_t framesConverged = 0;
  std::vector<Real> finalFitPerFrame;
  finalFitPerFrame.reserve(nFrames);

  for (std::size_t frame = 0; frame < nFrames; ++frame)
  {
    const Real t = static_cast<Real>(frame) / static_cast<Real>(nFrames);
    const Real angle = Real(2) * Real(M_PI) * t;

    WavyCircleLevelSet levelSet;
    levelSet.cx = Real(0.5) + orbitR * std::cos(angle);
    levelSet.cy = Real(0.5) + orbitR * std::sin(angle);
    levelSet.R0 = R0;
    levelSet.amp = amp;
    levelSet.k = kLobes;
    levelSet.phase = angle;

    std::cout << "\n--- Frame " << std::setw(2) << frame
              << " : c=(" << std::fixed << std::setprecision(4)
              << levelSet.cx << ", " << levelSet.cy << ")"
              << "  phase=" << std::setprecision(3) << angle << " rad\n";

    clearXDMFRegionAttributes(mesh);
    for (auto faceIt = mesh.getBoundary(); faceIt; ++faceIt)
      mesh.setAttribute({mesh.getDimension() - 1, faceIt->getIndex()},
                        boundaryAttribute);

    const auto cellMoments =
      collectCellMomentInfo(mesh,
                            [&](const Vec2& p) { return levelSet.phi(p); },
                            epsilon);

    std::unordered_map<Index, std::size_t> cellToLocal;
    std::vector<Index> localToCell;
    cellToLocal.reserve(cellMoments.size());
    localToCell.reserve(cellMoments.size());
    for (std::size_t local = 0; local < cellMoments.size(); ++local)
    {
      cellToLocal[cellMoments[local].index] = local;
      localToCell.push_back(cellMoments[local].index);
    }

    std::vector<Real> volumes(cellMoments.size());
    std::vector<Real> moments(cellMoments.size());
    for (std::size_t local = 0; local < cellMoments.size(); ++local)
    {
      volumes[local] = cellMoments[local].area;
      moments[local] = cellMoments[local].moment;
    }

    std::vector<MinSTCut::Edge> graphEdges;
    for (auto faceIt = mesh.getFace(); faceIt; ++faceIt)
    {
      const Index facet = faceIt->getIndex();
      const auto& incident =
        mesh.getConnectivity().getIncidence({1, 2}, facet);
      if (incident.size() != 2)
        continue;
      const auto itA = cellToLocal.find(incident[0]);
      const auto itB = cellToLocal.find(incident[1]);
      if (itA == cellToLocal.end() || itB == cellToLocal.end())
        continue;
      graphEdges.push_back({
          static_cast<Index>(itA->second),
          static_cast<Index>(itB->second),
          lambdaC * facetLength(mesh, facet),
          facet});
    }

    const MinSTCut cut;
    const MinSTCut::Result classified =
      cut.classify(volumes, moments, graphEdges);

    std::vector<Index> interfaceFacets;
    interfaceFacets.reserve(classified.cutEdges.size());
    for (const MinSTCut::Edge& edge : classified.cutEdges)
      if (edge.index != MinSTCut::InvalidIndex)
        interfaceFacets.push_back(edge.index);

    for (std::size_t local = 0; local < classified.labels.size(); ++local)
    {
      const Index cellIdx = localToCell[local];
      mesh.setAttribute(
          {mesh.getDimension(), cellIdx},
          classified.labels[local] == MinSTCut::Inside
            ? interiorAttribute
            : exteriorAttribute);
    }
    for (const Index facet : interfaceFacets)
      mesh.setAttribute({mesh.getDimension() - 1, facet}, interfaceAttribute);

    psi = Real(0);
    Distance::Eikonal<ScalarP1, Math::Vector<Real>>(psi)
      .setInterface(interfaceAttribute)
      .setInterior(interiorAttribute)
      .solve()
      .sign();

    phiGf = [&](const Geometry::Point& p) -> Real
    {
      const auto& X = p.getCoordinates();
      return levelSet.phi(vec2(X(0), X(1)));
    };

    RealFunction phi(
        [&](const Geometry::Point& p) -> Real
        {
          const auto& X = p.getPhysicalCoordinates();
          return levelSet.phi(vec2(X(0), X(1)));
        });
    AnalyticVectorFunction gradPhi(
        [&](const Geometry::Point& p) -> Math::SpatialVector<Real>
        {
          const auto& X = p.getPhysicalCoordinates();
          return levelSet.grad(vec2(X(0), X(1)));
        },
        /*dimension=*/2);

    u.getData().setZero();
    du.getData().setZero();

    auto computeInterfaceFit = [&]() -> Real
    {
      Real interfacePhi = 0;
      Real interfaceLen = 0;
      for (const Index facet : interfaceFacets)
      {
        const auto face = mesh.getFace(facet);
        const auto& vs = face->getVertices();
        const Vec2 a = mesh.getVertexCoordinates(vs[0]);
        const Vec2 b = mesh.getVertexCoordinates(vs[1]);
        const auto& fes = u.getFiniteElementSpace();
        const auto& dofsA = fes.getDOFs(0, vs[0]);
        const auto& dofsB = fes.getDOFs(0, vs[1]);
        const Vec2 ua = vec2(u.getData()(dofsA[0]), u.getData()(dofsA[1]));
        const Vec2 ub = vec2(u.getData()(dofsB[0]), u.getData()(dofsB[1]));
        const Real phiA = levelSet.phi(a + ua);
        const Real phiB = levelSet.phi(b + ub);
        const Real len = (b - a).norm();
        interfacePhi += Real(0.5) * (phiA * phiA + phiB * phiB) * len;
        interfaceLen += len;
      }
      return std::sqrt(interfacePhi / std::max(interfaceLen, Real(1e-30)));
    };

    Real interfaceFit = computeInterfaceFit();
    if (verbose)
    {
      std::size_t insideCount = 0;
      for (int lbl : classified.labels) if (lbl == MinSTCut::Inside) ++insideCount;
      std::cout << "    debug: facets=" << interfaceFacets.size()
                << "  inside=" << insideCount
                << "  outside=" << (classified.labels.size() - insideCount)
                << "  fit0=" << interfaceFit << "\n";
    }
    Real bestFit = interfaceFit;
    Math::Vector<Real> bestU = u.getData();
    Real minJ = Real(1);
    Real maxQRel = Real(1);
    Real lastAlpha = Real(0);
    Real maxStep = Real(0);
    Real acceptedStep = Real(0);
    std::size_t iterations = 0;

    // Exit-reason tag + ring buffer of best-fit history for stall test.
    const char* exitReason = "iter-budget";
    std::vector<Real> bestFitHistory;
    bestFitHistory.reserve(stallWindow + 1);
    bestFitHistory.push_back(bestFit);
    Real interfaceLen = 0;
    for (const Index facet : interfaceFacets)
      interfaceLen += facetLength(mesh, facet);
    LSRIntegratorParameters interfaceParams;
    interfaceParams.hRef = h;
    interfaceParams.quadratureOrder = qOrder;
    interfaceParams.interfaceAttribute = interfaceAttribute;
    interfaceParams.interfaceWeight = interfaceWeight;
    interfaceParams.interfaceNormalizer =
      Real(1) / std::max(interfaceLen * h, Real(1e-30));
    interfaceParams.interfaceGradientFloor = Real(1e-12);
    interfaceParams.interfaceLossSigma = charbonnierSigma;
    interfaceParams.interfaceLossKind = kInterfaceWelsch
      ? InterfaceLossKind::Welsch
      : InterfaceLossKind::Charbonnier;

    for (; iterations < maxIterations; ++iterations)
    {
      if (interfaceFit <= fitTol)
      {
        exitReason = "geometry-converged";
        break;
      }

      TrialFunction trial(vectorFes);
      TestFunction test(vectorFes);
      auto zero = VectorFunction{ Zero(), Zero() };
      LSRRegistration interfaceTerm(trial, test, u);

      // Two-step demons: solve direction first, then project onto
      // linearised admissibility via barriers.
      if (kTwoStep)
      {
        // ---- Step 1: cheap surface direction ----
        // L² mass + light Sobolev + surface coupling. No barriers.
        TrialFunction trialDir(vectorFes);
        TestFunction  testDir(vectorFes);
        LSRRegistration directionTerm(trialDir, testDir, u);
        Problem direction(trialDir, testDir);
        direction =
            Integral(trialDir, testDir)
          + (kDirectionEll * kDirectionEll)
              * LinearElasticityIntegral(trialDir, testDir)(
                  Real(0), Real(1))
          + directionTerm.InterfaceTangent(
                phi, gradPhi, interfaceParams)
          + directionTerm.InterfaceResidual(
                phi, gradPhi, interfaceParams)
          + DirichletBC(trialDir, zero).on(boundaryAttribute);
        Solver::CG(direction)
          .setTolerance(Real(1e-10)).setMaxIterations(1000).solve();
        GridFunction dField(vectorFes);
        dField.getData() = trialDir.getSolution().getData();

        // ---- Step 2: linearised-admissibility projection ----
        //   (M + γ_shape·Hess_barrier)·v = M·d − γ_shape·∇barrier
        // The barrier residual + tangent linearise the LSR shape/j
        // barriers around the current u; the mass M pulls v toward d.
        JacobianAdmissibilityBarrierSampled barrier(
            trial, test, u, static_cast<std::size_t>(qOrder));
        BarrierParameters barrierParams;
        barrierParams.jMin = jMinRatio;
        barrierParams.domainMeasure = Real(1);
        barrierParams.qBarrierWeight = kQBarrierWeight;
        barrierParams.qBarrierAct = kQBarrierAct;
        barrierParams.qBarrierMax = kQBarrierMax;
        barrierParams.jBarrierWeight = kJBarrierWeight;
        barrierParams.jBarrierSafeRatio = kJBarrierSafeRatio;
        barrierParams.jVolumeTetherWeight = kJVolTetherWeight;
        RealFunction<Real> shapeWeight(kGammaShape);
        // VectorFunction wrapper so Rodin can dispatch a vector-vector
        // L² coupling on the linear form. dField is the direction.
        auto dCoupling = VectorFunction{
            [&](const Geometry::Point& p) -> Real {
              const auto val = dField.getValue(p);
              return val(0);
            },
            [&](const Geometry::Point& p) -> Real {
              const auto val = dField.getValue(p);
              return val(1);
            }
          };
        Problem projection(trial, test);
        projection =
            Integral(trial, test)
          + barrier.TangentPSDProjected(shapeWeight, barrierParams)
          + barrier.Residual(shapeWeight, barrierParams)
          + Integral(dCoupling, test)  // sign TBD by empirical A/B
          + DirichletBC(trial, zero).on(boundaryAttribute);
        Solver::CG(projection)
          .setTolerance(Real(1e-10)).setMaxIterations(1000).solve();
        du.getData() = trial.getSolution().getData();

        maxStep = du.getData().cwiseAbs().maxCoeff();
        if (verbose)
          std::cout << "    debug iter=" << iterations
                    << "  two-step  |d|=" << dField.getData().norm()
                    << "  |du|=" << du.getData().norm()
                    << "  maxStep=" << maxStep << "\n";
      }
      else
      {
      Problem extension(trial, test);
      // Volumetric push parameters (used only when kVolumetricPush).
      LSRIntegratorParameters pushParams;
      pushParams.rhoS = kPushWeight;
      pushParams.deltaW = kPushDeltaW;
      pushParams.hRef = h;
      pushParams.quadratureOrder = qOrder;
      pushParams.normalizer = Real(1) / (h * h);
      pushParams.lossSigma = kPushLossSigma;
      if (kQualityExtension && kVolumetricPush)
      {
        JacobianAdmissibilityBarrierSampled barrier(
            trial, test, u, static_cast<std::size_t>(qOrder));
        BarrierParameters barrierParams;
        barrierParams.jMin = jMinRatio;
        barrierParams.domainMeasure = Real(1);
        barrierParams.qBarrierWeight = kQBarrierWeight;
        barrierParams.qBarrierAct = kQBarrierAct;
        barrierParams.qBarrierMax = kQBarrierMax;
        barrierParams.jBarrierWeight = kJBarrierWeight;
        barrierParams.jBarrierSafeRatio = kJBarrierSafeRatio;
        barrierParams.jVolumeTetherWeight = kJVolTetherWeight;
        RealFunction<Real> shapeWeight(kGammaShape);
        extension =
            barrier.TangentPSDProjected(shapeWeight, barrierParams)
          + barrier.Residual(shapeWeight, barrierParams)
          + interfaceTerm.InterfaceTangent(phi, gradPhi, interfaceParams)
          + interfaceTerm.InterfaceResidual(phi, gradPhi, interfaceParams)
          + interfaceTerm.Tangent(phi, gradPhi, psi, pushParams)
          + interfaceTerm.Residual(phi, gradPhi, psi, pushParams)
          + DirichletBC(trial, zero).on(boundaryAttribute);
      }
      else if (kQualityExtension && kDirectionDirichlet)
      {
        // Step A: closed-form direction at each skeleton vertex.
        //   d(v) = -ρ'(r_Γ(v + u(v))) · n_φ(v + u(v))
        // Stored in a vector GridFunction with values at skeleton
        // vertices, zero elsewhere — acts as the L² target field.
        GridFunction dField(vectorFes);
        dField.getData().setZero();
        std::set<Index> skeletonVerts;
        for (const Index facet : interfaceFacets)
        {
          const auto face = mesh.getFace(facet);
          for (const Index v : face->getVertices())
            skeletonVerts.insert(v);
        }
        const auto& fes = u.getFiniteElementSpace();
        for (const Index v : skeletonVerts)
        {
          const Vec2 X = mesh.getVertexCoordinates(v);
          const auto& dofs = fes.getDOFs(0, v);
          const Vec2 uX =
            vec2(u.getData()(dofs[0]), u.getData()(dofs[1]));
          const Vec2 y = X + uX;
          const Vec2 g = levelSet.grad(y);
          const Real gn = std::max(g.norm(), Real(1e-12));
          const Real rGamma = levelSet.phi(y) / gn;
          const Real w = kInterfaceWelsch
            ? welschWeight(rGamma, charbonnierSigma)
            : charbonnierWeight(rGamma, charbonnierSigma);
          const Real coef = -rGamma * w / gn;
          dField.getData()(dofs[0]) = coef * g(0);
          dField.getData()(dofs[1]) = coef * g(1);
        }

        // Step B: M-metric soft-penalty Galerkin.
        //   (d*, v)_L² + γ_shape · (d*, v)_M = (d, v)_L² − γ_shape · ∇B(u^k) · v
        //
        // where M = linearised barrier + Q_shape Hessian. One direct
        // linear solve, no inner Newton. The L² coupling against d
        // pulls δu* toward the closed-form Newton direction; the
        // M-Hessian regularises into the linearised admissibility
        // metric; the barrier residual gives a Newton-like restoring
        // force toward the current admissible state.
        JacobianAdmissibilityBarrierSampled barrier(
            trial, test, u, static_cast<std::size_t>(qOrder));
        BarrierParameters barrierParams;
        barrierParams.jMin = jMinRatio;
        barrierParams.domainMeasure = Real(1);
        barrierParams.qBarrierWeight = kQBarrierWeight;
        barrierParams.qBarrierAct = kQBarrierAct;
        barrierParams.qBarrierMax = kQBarrierMax;
        barrierParams.jBarrierWeight = kJBarrierWeight;
        barrierParams.jBarrierSafeRatio = kJBarrierSafeRatio;
        barrierParams.jVolumeTetherWeight = kJVolTetherWeight;
        RealFunction<Real> shapeWeight(kGammaShape);

        auto dFunc = VectorFunction{
            [&dField](const Geometry::Point& p) -> Real {
              const auto val = dField.getValue(p);
              return val(0);
            },
            [&dField](const Geometry::Point& p) -> Real {
              const auto val = dField.getValue(p);
              return val(1);
            }
          };

        extension =
            Integral(trial, test)
          + barrier.TangentPSDProjected(shapeWeight, barrierParams)
          + barrier.Residual(shapeWeight, barrierParams)
          - Integral(dFunc, test)
          + DirichletBC(trial, zero).on(boundaryAttribute);
      }
      else if (kQualityExtension)
      {
        JacobianAdmissibilityBarrierSampled barrier(
            trial, test, u, static_cast<std::size_t>(qOrder));
        BarrierParameters barrierParams;
        barrierParams.jMin = jMinRatio;
        barrierParams.domainMeasure = Real(1);
        barrierParams.qBarrierWeight = kQBarrierWeight;
        barrierParams.qBarrierAct = kQBarrierAct;
        barrierParams.qBarrierMax = kQBarrierMax;
        barrierParams.jBarrierWeight = kJBarrierWeight;
        barrierParams.jBarrierSafeRatio = kJBarrierSafeRatio;
        barrierParams.jVolumeTetherWeight = kJVolTetherWeight;
        RealFunction<Real> shapeWeight(kGammaShape);
        extension =
            barrier.TangentPSDProjected(shapeWeight, barrierParams)
          + barrier.Residual(shapeWeight, barrierParams)
          + interfaceTerm.InterfaceTangent(phi, gradPhi, interfaceParams)
          + interfaceTerm.InterfaceResidual(phi, gradPhi, interfaceParams)
          + DirichletBC(trial, zero).on(boundaryAttribute);
      }
      else if (kVolumetricPush)
      {
        extension =
            Integral(trial, test)
          + extensionEll * extensionEll
              * LinearElasticityIntegral(trial, test)(
                  extensionLambda, extensionMu)
          + interfaceTerm.InterfaceTangent(phi, gradPhi, interfaceParams)
          + interfaceTerm.InterfaceResidual(phi, gradPhi, interfaceParams)
          + interfaceTerm.Tangent(phi, gradPhi, psi, pushParams)
          + interfaceTerm.Residual(phi, gradPhi, psi, pushParams)
          + DirichletBC(trial, zero).on(boundaryAttribute);
      }
      else
      {
        extension =
            Integral(trial, test)
          + extensionEll * extensionEll
              * LinearElasticityIntegral(trial, test)(
                  extensionLambda, extensionMu)
          + interfaceTerm.InterfaceTangent(phi, gradPhi, interfaceParams)
          + interfaceTerm.InterfaceResidual(phi, gradPhi, interfaceParams)
          + DirichletBC(trial, zero).on(boundaryAttribute);
      }

      auto solver = Solver::CG(extension);
      solver.setTolerance(Real(1e-10)).setMaxIterations(1000).solve();
      du.getData() = trial.getSolution().getData();

      maxStep = du.getData().cwiseAbs().maxCoeff();
      if (verbose)
      {
        const auto& ls = extension.getLinearSystem();
        std::cout << "    debug iter=" << iterations
                  << "  |rhs|=" << ls.getVector().norm()
                  << "  |du|=" << du.getData().norm()
                  << "  maxStep=" << maxStep << "\n";
      }
      }  // end else (single-solve path)
      if (!std::isfinite(maxStep) || maxStep <= stepTol)
      {
        lastAlpha = Real(0);
        exitReason = std::isfinite(maxStep)
          ? "step-below-stepTol" : "solve-nonfinite";
        break;
      }

      // Pre-clip du by the trust radius BEFORE the admissibility
      // predictor. The sampled quadratic predictor's coefficients on
      // the raw du can be so steep on a sub-resolved skeleton that
      // even infinitesimal α violates jSafe (predicted.alphaMax → 0).
      // Clipping du first decouples the admissibility test from the
      // magnitude of the bulk solve and recovers a usable α.
      Real preScale = Real(1);
      if (trustRadius > Real(0) && maxStep > Real(0))
        preScale = std::min(Real(1), trustRadius / maxStep);
      du.getData() *= preScale;
      const Real clippedMaxStep = preScale * maxStep;

      // Demons-style: there is no Armijo "α ≤ 1" cap. The admissibility
      // predictor is the only governor of step size. We do cap at the
      // unclipped-du recovery scale (1/preScale) so α cannot
      // structurally overshoot the bulk solve's own preferred step.
      const Real alphaRecover =
        preScale > Real(0) ? Real(1) / preScale : Real(1);
      const auto predicted = predictSampledQuadraticAlpha(
          u, u.getData(), du.getData(),
          jLineSearchRatio, alphaSafety, qOrder);
      Real alpha = std::isfinite(predicted.alphaMax)
        ? std::min(alphaRecover, predicted.alphaMax)
        : alphaRecover;

      if (verbose)
        std::cout << "    debug iter=" << iterations
                  << "  preScale=" << preScale
                  << "  clippedMax=" << clippedMaxStep
                  << "  predictor.alphaMax=" << predicted.alphaMax
                  << "  alpha=" << alpha << "\n";
      // Keep maxStep referring to the clipped direction for downstream
      // step bookkeeping.
      maxStep = clippedMaxStep;

      if (alpha <= Real(0))
      {
        lastAlpha = Real(0);
        exitReason = "predictor-rejects-all-alpha";
        break;
      }

      // Bisection backoff: predictor.alphaMax is a quadratic sampled
      // estimate; the actual evaluateLSRAdmissibilitySampled may reject
      // at the predicted α. Backoff up to a few times before giving up.
      const Math::Vector<Real> previousU = u.getData();
      const std::size_t kAdmissibilityBackoffMax = 6;
      LSRAdmissibilityReport adm{};
      bool stepAccepted = false;
      for (std::size_t backoff = 0; backoff <= kAdmissibilityBackoffMax; ++backoff)
      {
        u.getData() = previousU + alpha * du.getData();
        adm = evaluateLSRAdmissibilitySampled(
            u, u.getData(), jMinRatio, qOrder);
        if (adm.inadmissibleCount == 0 && adm.minJRatio > jLineSearchRatio)
        {
          stepAccepted = true;
          break;
        }
        alpha *= Real(0.5);
        if (alpha * maxStep < stepTol)
          break;
      }
      if (!stepAccepted)
      {
        u.getData() = previousU;
        exitReason = "admissibility-bisection-exhausted";
        break;
      }

      minJ = adm.minJRatio;
      maxQRel = adm.maxQRel;
      interfaceFit = computeInterfaceFit();
      if (!std::isfinite(interfaceFit))
      {
        // NaN/inf is a hard failure — rollback and break.
        u.getData() = previousU;
        interfaceFit = bestFit;
        exitReason = "fit-nonfinite";
        break;
      }

      // Accept any admissible step. Track best-so-far separately;
      // the iteration is allowed to take non-monotone moves and only
      // stops at the iteration budget or a stall criterion below.
      lastAlpha = alpha;
      acceptedStep = alpha * maxStep;
      if (interfaceFit < bestFit)
      {
        bestFit = interfaceFit;
        bestU = u.getData();
      }

      // Per-iter trace (one line per accepted step).
      if (trace)
      {
        // Sample admissibility distribution at quadrature points of u.
        // qMax used as a reference scale for "near boundary" count.
        constexpr Real kQMaxRef = Real(10);
        const auto admDist =
          evaluateAdmissibilityDistribution(
              u, jLineSearchRatio, kQMaxRef, qOrder);
        std::cout << "      it=" << std::setw(3) << iterations
                  << "  fit=" << std::scientific << std::setprecision(3)
                  << interfaceFit
                  << "  α=" << alpha
                  << "  Δu_∞=" << acceptedStep
                  << "\n"
                  << "         adm: j ["
                  << admDist.minJ << " | "
                  << admDist.p10J << ", "
                  << admDist.p50J << ", "
                  << admDist.p90J << "]"
                  << "  near_jsafe=" << admDist.numNearJSafe << "/"
                  << admDist.numQuadPoints
                  << "  inadm=" << admDist.numInadmissible
                  << "\n"
                  << "         adm: Q ["
                  << admDist.p10Q << ", "
                  << admDist.p50Q << ", "
                  << admDist.p90Q << " | "
                  << admDist.maxQ << "]"
                  << "  near_qmax=" << admDist.numNearQMax << "/"
                  << admDist.numQuadPoints
                  << "  predictor.α=" << predicted.alphaMax
                  << '\n';
      }

      // Stall detector: relative improvement of the best-so-far over
      // the last `stallWindow` accepted iterations. The ratio
      // (bestFit_{k-window} - bestFit_k) / bestFit_{k-window} measures
      // how much the floor has come down across the window. When it
      // falls below stallTol the iteration is no longer making
      // measurable progress — declare stall.
      bestFitHistory.push_back(bestFit);
      if (bestFitHistory.size() > stallWindow + 1)
      {
        const Real prior =
          bestFitHistory[bestFitHistory.size() - 1 - stallWindow];
        const Real improvement =
          (prior - bestFit) / std::max(prior, Real(1e-30));
        if (improvement < stallTol)
        {
          exitReason = "stall";
          break;
        }
      }

      if (verbose)
      {
        std::cout << "      demon " << std::setw(2) << iterations + 1
                  << " fit=" << std::scientific << std::setprecision(3)
                  << interfaceFit
                  << " alpha=" << lastAlpha
                  << " step=" << maxStep
                  << " min_j=" << minJ
                  << " max_qrel=" << maxQRel << '\n';
      }
    }

    u.getData() = bestU;
    interfaceFit = bestFit;
    const auto bestAdm = evaluateLSRAdmissibilitySampled(
        u, u.getData(), jMinRatio, qOrder);
    minJ = bestAdm.minJRatio;
    maxQRel = bestAdm.maxQRel;

    const bool converged = interfaceFit <= fitTol;
    if (converged)
      ++framesConverged;
    finalFitPerFrame.push_back(interfaceFit);

    const std::size_t D = mesh.getDimension();
    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const Index cellIdx = cellIt->getIndex();
      const std::size_t local = cellToLocal.at(cellIdx);
      const Index dof = p0Fes.getGlobalIndex({D, cellIdx}, 0);
      cellLabel.getData()(dof) =
        static_cast<Real>(classified.labels[local]);
      phaseMoment.getData()(dof) = cellMoments[local].moment;
    }

    updateMovedMeshFromDisplacement(mesh, moved, u);
    for (Index i = 0; i < static_cast<Index>(mesh.getCellCount()); ++i)
    {
      if (const auto a = mesh.getAttribute(D, i))
        moved.setAttribute({D, i}, *a);
      else
        moved.setAttribute({D, i}, Attribute{0});
    }
    for (auto it = mesh.getFace(); it; ++it)
    {
      const Index idx = it->getIndex();
      if (const auto a = mesh.getAttribute(D - 1, idx))
        moved.setAttribute({D - 1, idx}, *a);
      else
        moved.setAttribute({D - 1, idx}, Attribute{0});
    }

    auto [srcCache, srcToLocal] = precomputeCellGeometry(mesh);
    auto [dstCache, dstToLocal] = precomputeCellGeometry(moved);
    for (auto cellIt = moved.getCell(); cellIt; ++cellIt)
    {
      const Index cellIdx = cellIt->getIndex();
      const std::size_t srcLocal = srcToLocal.at(cellIdx);
      const std::size_t dstLocal = dstToLocal.at(cellIdx);
      const Index dof = p0FesMoved.getGlobalIndex({D, cellIdx}, 0);
      const auto& src = srcCache[srcLocal];
      const auto& dst = dstCache[dstLocal];
      const Real sigDetAu = static_cast<Real>(src.sigmaK) * dst.detAK;
      const Real jK = sigDetAu / src.Jscale;
      jMoved.getData()(dof) = jK;
      const Math::SpatialMatrix<Real> F = dst.A * src.A.inverse();
      qRelMoved.getData()(dof) =
        F.squaredNorm() / (Real(2) * std::max(jK, Real(1e-30)));
      movedLabel.getData()(dof) =
        static_cast<Real>(classified.labels[cellToLocal.at(cellIdx)]);
    }

    phiMoved = [&](const Geometry::Point& p) -> Real
    {
      const auto& X = p.getCoordinates();
      return levelSet.phi(vec2(X(0), X(1)));
    };

    std::cout << "    Demons it=" << iterations
              << "  fit=" << std::scientific << std::setprecision(3)
              << interfaceFit
              << "  alpha=" << lastAlpha
              << "  step=" << acceptedStep
              << "  min_j=" << minJ
              << "  max_qrel=" << maxQRel
              << "  converged=" << (converged ? "yes" : "best-effort")
              << "  exit=" << exitReason
              << '\n';

    xdmf.write(t).flush();
  }

  xdmf.close();

  std::cout << "\nSummary\n";
  std::cout << "  frames converged: " << framesConverged
            << " / " << nFrames << '\n';
  if (!finalFitPerFrame.empty())
  {
    const Real fitMin =
      *std::min_element(finalFitPerFrame.begin(), finalFitPerFrame.end());
    const Real fitMax =
      *std::max_element(finalFitPerFrame.begin(), finalFitPerFrame.end());
    Real fitMean = 0;
    for (Real x : finalFitPerFrame)
      fitMean += x;
    fitMean /= static_cast<Real>(finalFitPerFrame.size());
    std::cout << "  ||phi(X+u)||_RMS  min=" << std::scientific
              << std::setprecision(3) << fitMin
              << "  mean=" << fitMean
              << "  max=" << fitMax << '\n';
  }

  return 0;
}
