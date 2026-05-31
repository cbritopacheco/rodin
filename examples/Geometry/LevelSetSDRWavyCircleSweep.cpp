/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//
// Wavy-circle sweep test for the SDR pipeline.
//
// A non-convex implicit shape — a circle whose radius varies sinusoidally
// in the polar angle (k azimuthal lobes) — is moved around the unit
// square along a circular orbit. The whole classify-then-displace
// pipeline (s-t cut + FMM + SDR Newton) is rerun from scratch at every
// frame, and a transient XDMF stream is emitted with both the
// low-fidelity classification and the high-fidelity moved mesh per
// frame.
//
// Compared with `LevelSetSDRReconstruction`, this example exercises:
//   - A non-circular, non-convex interface that genuinely requires
//     tangential motion of the SDR (the cosine perturbation cannot be
//     absorbed by uniform dilation of the cell map).
//   - Repeated classification + Newton solves at moving interface
//     locations: the classified topology changes between frames and the
//     transient XDMF writer dumps a fresh mesh each frame.
//   - The shard-aware XDMF writers under non-trivial cell-attribute
//     turnover.
//
// FES-independence and Newton-tangent choices follow the same honest
// caveats as the parent example: the barrier + shape closed-form
// algebra is specialised to affine P1 triangles in 2D; SDR is the
// Gauss--Newton flavour (CTAD-selected by passing phi/grad/sLF/params,
// without hess). The example is therefore a 2D triangular prototype.
//
#include <Rodin/Adaptation.h>
#include <Rodin/Assembly.h>
#include <Rodin/Distance/Eikonal.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Math.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>
#include <Rodin/Solver/NewtonSolver.h>
#include <Rodin/Solver/SparseLU.h>
#include <Rodin/Variational.h>

#include "LevelSetSDRWavyCircleSweepDiagnostics.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Adaptation;
namespace Diag = LevelSetSDRWavyCircleSweepSupport;

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

  // -------------------------------------------------------------------------
  // Wavy circle level set.
  //
  //   phi(x, y) = r - R(theta),
  //   r        = || (x, y) - c ||,
  //   theta    = atan2(y - cy, x - cx),
  //   R(theta) = R0 + amp * cos(k * theta + phase).
  //
  // phi is smooth away from the centre. Like the plain circle SDF, it is
  // NOT a strict signed distance (the radial perturbation makes the
  // gradient norm differ from 1) but the SDR pipeline only requires a
  // sufficiently smooth implicit function — phi, grad phi at quadrature
  // points — to operate.
  // -------------------------------------------------------------------------
  struct WavyCircleLevelSet
  {
    Real cx = Real(0.5);
    Real cy = Real(0.5);
    Real R0 = Real(0.18);
    Real amp = Real(0.04);
    Real k = Real(5);
    Real phase = Real(0);

    Real phi(const Vec2& p) const
    {
      const Real dx = p(0) - cx;
      const Real dy = p(1) - cy;
      const Real r = std::sqrt(dx * dx + dy * dy);
      const Real theta = std::atan2(dy, dx);
      const Real R = R0 + amp * std::cos(k * theta + phase);
      return r - R;
    }

    // grad phi = (dx/r, dy/r) - dR/dtheta * grad(theta).
    // grad(theta) = (-dy / r^2, dx / r^2).
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

    // Hess phi = Hess(r) - R'(theta) Hess(theta) - R''(theta) grad(theta) grad(theta)^T.
    //
    //   Hess(r)_ij = (delta_ij - dx_i dx_j / r^2) / r,
    //
    //   grad(theta) = (-dy, dx) / r^2,
    //
    //   Hess(theta)_xx = 2 dx dy / r^4,
    //   Hess(theta)_yy = -2 dx dy / r^4,
    //   Hess(theta)_xy = (dy^2 - dx^2) / r^4.
    //
    // Consumed by the PSD-projected full-Newton tangent.
    Math::SpatialMatrix<Real> hess(const Vec2& p) const
    {
      const Real dx = p(0) - cx;
      const Real dy = p(1) - cy;
      const Real r2 = dx * dx + dy * dy;
      const Real r = std::max(std::sqrt(r2), Real(1e-14));
      const Real r3 = r * r * r;
      const Real r4 = r2 * r2;

      const Real theta = std::atan2(dy, dx);
      const Real ka = k * theta + phase;
      const Real Rpp = -amp * k * std::sin(ka);
      const Real Rpp2 = -amp * k * k * std::cos(ka);

      // Hess(r)
      Math::SpatialMatrix<Real> Hr(2, 2);
      Hr(0, 0) = (Real(1) - dx * dx / r2) / r;   // = dy^2 / r^3
      Hr(1, 1) = (Real(1) - dy * dy / r2) / r;
      Hr(0, 1) = -dx * dy / r3;
      Hr(1, 0) = Hr(0, 1);

      // Hess(theta)
      Math::SpatialMatrix<Real> Ht(2, 2);
      Ht(0, 0) =  Real(2) * dx * dy / r4;
      Ht(1, 1) = -Real(2) * dx * dy / r4;
      Ht(0, 1) = (dy * dy - dx * dx) / r4;
      Ht(1, 0) = Ht(0, 1);

      // grad(theta) outer product
      const Vec2 gth = vec2(-dy / r2, dx / r2);
      Math::SpatialMatrix<Real> GG(2, 2);
      GG(0, 0) = gth(0) * gth(0);
      GG(0, 1) = gth(0) * gth(1);
      GG(1, 0) = GG(0, 1);
      GG(1, 1) = gth(1) * gth(1);

      Math::SpatialMatrix<Real> H(2, 2);
      H = Hr - Rpp * Ht - Rpp2 * GG;
      return H;
    }
  };

  // -------------------------------------------------------------------------
  // Triangle-barycentric quadrature for the phase moments (triangle-only,
  // NOT FES-independent — same caveat as the parent example).
  // -------------------------------------------------------------------------
  constexpr std::array<std::array<Real, 3>, 3> TriangleBarycentricQuadrature = {{
    {{ Real(2) / 3, Real(1) / 6, Real(1) / 6 }},
    {{ Real(1) / 6, Real(2) / 3, Real(1) / 6 }},
    {{ Real(1) / 6, Real(1) / 6, Real(2) / 3 }}
  }};

  Real applyPhaseMomentMap(Real phi, Real epsilon)
  {
    return std::tanh(phi / epsilon);
  }

  Vec2 interpolateVec(const std::array<Vec2, 3>& values,
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
    std::array<Index, 3> vertices = {{ 0, 0, 0 }};
  };

  // Templated on a `phi(Vec2)` callable so the helper is reusable across
  // any level set type.
  template <class PhiFn>
  std::vector<CellMomentInfo> collectCellMomentInfo(
      const LocalMesh& mesh, PhiFn&& phiFn, Real epsilon)
  {
    std::vector<CellMomentInfo> cells;
    cells.reserve(mesh.getCellCount());

    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const auto& cellPolytope = *cellIt;
      const auto& vertices = cellPolytope.getVertices();
      if (vertices.size() != 3)
        throw std::runtime_error(
            "LevelSetSDRWavyCircleSweep expects triangular cells.");

      CellMomentInfo info;
      info.index = cellPolytope.getIndex();
      for (size_t i = 0; i < 3; ++i)
      {
        info.vertices[i] = vertices[i];
        info.x[i] = mesh.getVertexCoordinates(vertices[i]);
      }

      const Vec2 e1 = info.x[1] - info.x[0];
      const Vec2 e2 = info.x[2] - info.x[0];
      info.area =
        std::abs(Real(0.5) * (e1(0) * e2(1) - e1(1) * e2(0)));

      Real moment = 0;
      for (const auto& bary : TriangleBarycentricQuadrature)
      {
        const Vec2 xq = interpolateVec(info.x, bary);
        moment += applyPhaseMomentMap(phiFn(xq), epsilon);
      }
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

  // -------------------------------------------------------------------------
  // Reset cell-region / interface / boundary attributes so that the mesh
  // can be reclassified from scratch at the next frame.
  // -------------------------------------------------------------------------
  void clearXDMFRegionAttributes(LocalMesh& mesh)
  {
    const std::size_t D = mesh.getDimension();
    for (Index i = 0; i < static_cast<Index>(mesh.getCellCount()); ++i)
      mesh.setAttribute({D, i}, Attribute{0});
    for (auto it = mesh.getFace(); it; ++it)
      mesh.setAttribute({D - 1, it->getIndex()}, Attribute{0});
  }

}

int main(int, char**)
{
  // -------------------------------------------------------------------------
  // Step 0: frame schedule, mesh, and constants.
  // -------------------------------------------------------------------------
  constexpr std::size_t n = 20;        ///< 1D mesh resolution (n x n nodes).
  constexpr std::size_t nFrames = 20;  ///< Number of orbit snapshots.
  constexpr Real        orbitR = 0.10; ///< Radius of the centre's orbit.
  constexpr Real        amp    = 0.025;///< Wavy-circle radial amplitude.
  constexpr Real        R0     = 0.20; ///< Wavy-circle nominal radius.
  constexpr Real        kLobes = 6;    ///< Number of azimuthal lobes.

  // With orbitR=0.10 and R0=0.20 the contour reach is at most 0.30 from
  // the orbit-of-orbit centre (0.5, 0.5), so the contour stays in
  // [0.20, 0.80] on each axis — i.e. >= 0.20 from the Dirichlet
  // boundary in every direction. This is the headroom that lets Newton
  // take the first un-damped step (damping=1.0 below) without any cell
  // immediately flipping orientation. See the comments around the
  // NewtonSolver setup for the damping rationale.

  const Real h = Real(1) / static_cast<Real>(n - 1);
  const Real epsilon = 1.25 * h;
  const Real lambdaC = 0.008;
  const Real delta   = 1.75 * h;
  const Real deltaW  = 1.5 * delta;

  // -------------------------------------------------------------------------
  // Tikhonov + harmonic regularization knobs.
  //
  // Adding `eps_T * ∫ du·v + eps_H * ∫ ∇du:∇v` to the Newton tangent makes
  // the local quadratic model
  //
  //   min  ½ ‖R + K du‖² + (eps_T/2) ‖du‖_L²² + (eps_H/2) ‖∇du‖_L²²,
  //
  // which is the classic Tikhonov / harmonic regularisation. Both terms are
  // symmetric positive (semi)definite; together with the shape-energy
  // coercivity on H¹₀, they raise the smallest eigenvalue of the global
  // tangent, so Newton contracts even when the PSD-projected r·∇²φ term
  // is small at the minimum.
  //
  // Choice of magnitude:
  //   - Harmonic: scaling with the SDR per-cell magnitude ρs/M_w ≈ O(1)
  //     gives a regularisation that is mesh-independent in `eps_H`.
  //     At n=30 (h ≈ 0.034) the GN noise floor swamps the iterate
  //     unless eps_H >= 1e-2 or so.
  //   - Tikhonov: the mass matrix has h² scaling per cell, so for
  //     mesh-independent regularisation `eps_T` must scale ~1/h². We
  //     default to 0 here and rely on harmonic.
  //
  // Both knobs are opt-in. Set to 0 to recover the un-regularised
  // PSD-projected Newton iteration that was already used at n=60.
  const Real kHarmonicEps = Real(0);      ///< harmonic / Dirichlet energy weight
  const Real kTikhonovEps = Real(0);      ///< L² Tikhonov weight

  // -------------------------------------------------------------------------
  // Admissibility-first backtracking line search.
  //
  //   while (alpha >= alphaMin)
  //     uTrial = u + alpha * p
  //     adm    = evaluateAdmissibility(uTrial)
  //     if (adm.minJ > jMin && adm.inadmissibleCount == 0): accept
  //     else: alpha *= reduction
  //
  // When disabled, the iteration falls back to the static `kDamping`
  // factor (the original failing baseline).
  const bool   kLineSearchEnabled   = true;
  const Real   kLineSearchAlphaInit = Real(1);
  const Real   kLineSearchReduction = Real(0.5);
  const Real   kLineSearchAlphaMin  = Real(1e-6);
  // Safety margin above j_safe used as the line-search admissibility
  // floor: accept only steps with min_K j_K^u > kLineSearchSafetyMargin
  // × j_safe. The previous default (j_min = 10⁻⁸) only protected against
  // catastrophic flips, allowing iterates to enter the FLOOR-BARRIER
  // active region (j ∈ (j_min, j_safe)) where the barrier Hessian is
  // very stiff and subsequent Newton directions become explosively
  // compressive. Keeping the iterate safely above j_safe avoids this.
  const Real   kLineSearchSafetyMargin = Real(10);

  // ----- Initial guess strategy --------------------------------------------
  //   Cold    : u₀ = 0
  //   Warm    : u₀ = previous frame's converged u (cold fallback on failure)
  //   Hilbert : u₀ = Riesz lift of −δE_SDR(0) in the H¹₀ inner product
  //             (a(u₀, v) = −⟨δE_SDR(0), v⟩  on V₀)
  //
  // All three drop back to a cold-start retry if the first attempt fails
  // the convergence criterion (Cold is its own retry; this is a no-op).
  enum class InitialGuessStrategy { Cold, Warm, Hilbert };
  constexpr InitialGuessStrategy kInitialGuessStrategy =
      InitialGuessStrategy::Hilbert;

  // Hilbert-lift bilinear form choice on V₀:
  //   Harmonic     : a(u,v) = ∫ ∇u : ∇v dx                  (Poincaré)
  //   Elasticity   : a(u,v) = ∫ λ(div u)(div v) + 2μ ε(u):ε(v) dx  (Korn II)
  //   ShapeHessian : a(u,v) = δ²E_shape(0)[u,v]              (matches the
  //                                                          well-posedness
  //                                                          theorem's α)
  enum class HilbertBilinearForm { Harmonic, Elasticity, ShapeHessian };
  constexpr HilbertBilinearForm kHilbertBilinearForm =
      HilbertBilinearForm::Harmonic;
  const Real kHilbertElasticityLambda = Real(0);      ///< pure deviatoric
  const Real kHilbertElasticityMu     = Real(0.5);    ///< 2μ = 1, ≈ harmonic scale

  // ----- Convergence criterion ---------------------------------------------
  //   Residual : classic ||R|| < tol (fails to recognise the GN data-floor).
  //   Geometry : ||phi(X+u)||_RMS < ~h^2 (honest "did we capture Γ?").
  //   Either   : converged if EITHER residual or geometry tolerance holds.
  enum class ConvergenceMode { Residual, Geometry, Either };
  constexpr ConvergenceMode kConvergenceMode = ConvergenceMode::Either;
  const Real        kInterfaceFitTol = Real(2) * h * h;   ///< ~2·h² target
  constexpr Real    kResidualAbsTol  = Real(1e-6);
  constexpr Real    kResidualRelTol  = Real(1e-3);

  constexpr Attribute interiorAttribute  = 1;
  constexpr Attribute exteriorAttribute  = 2;
  constexpr Attribute interfaceAttribute = 10;
  constexpr Attribute boundaryAttribute  = 20;

  // Build the background mesh ONCE; only attributes change per frame.
  LocalMesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { n, n });
  mesh.scale(h);
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 2);
  mesh.getConnectivity().compute(2, 2);
  mesh.getConnectivity().compute(0, 0);

  // Tag the static boundary attribute once; the per-frame loop will
  // re-tag interior / exterior / interface cells from the classifier.
  for (auto faceIt = mesh.getBoundary(); faceIt; ++faceIt)
    mesh.setAttribute({mesh.getDimension() - 1, faceIt->getIndex()},
                      boundaryAttribute);

  // -------------------------------------------------------------------------
  // FE spaces and persistent grid functions used across frames.
  // -------------------------------------------------------------------------
  using ScalarP1 = P1<Real, LocalMesh>;
  using VectorP1 = P1<Math::SpatialVector<Real>, LocalMesh>;
  using ScalarP0 = P0<Real, LocalMesh>;

  ScalarP1 p1Fes(mesh);
  ScalarP0 p0Fes(mesh);
  VectorP1 vectorFes(mesh, 2);

  GridFunction sLF(p1Fes);            sLF.setName("s_LF");
  GridFunction phiGf(p1Fes);          phiGf.setName("phi");
  GridFunction cellLabel(p0Fes);      cellLabel.setName("cell_label");
  GridFunction phaseMoment(p0Fes);    phaseMoment.setName("phase_moment");
  GridFunction sigmaKgf(p0Fes);       sigmaKgf.setName("sigma_K");

  TrialFunction du(vectorFes);
  TestFunction  v(vectorFes);
  GridFunction  u(vectorFes);         u.setName("displacement");

  auto zero = VectorFunction{ Zero(), Zero() };

  // -------------------------------------------------------------------------
  // HF (moved) side. Same combinatorics as `mesh`; only vertex coordinates
  // change per frame, so the FES on `moved` is stable across frames and
  // we can register persistent GridFunctions with the transient XDMF
  // grid. The XDMF writer dumps the moved mesh anew at every snapshot
  // (MeshPolicy::Transient) so ParaView sees the curved fit per frame.
  // -------------------------------------------------------------------------
  LocalMesh moved(mesh);
  ScalarP1 p1FesMoved(moved);
  ScalarP0 p0FesMoved(moved);

  GridFunction cellLabelHF(p0FesMoved); cellLabelHF.setName("cell_label");
  GridFunction jKgf(p0FesMoved);        jKgf.setName("j");
  GridFunction qShape(p0FesMoved);      qShape.setName("q_shape");
  GridFunction phiMoved(p1FesMoved);    phiMoved.setName("phi_moved");

  // -------------------------------------------------------------------------
  // XDMF writer in transient mode (two grids: LF and HF).
  // -------------------------------------------------------------------------
  IO::XDMF xdmf("LevelSetSDRWavyCircleSweep");
  auto lfGrid = xdmf.grid("LF");
  lfGrid.setMesh(mesh, IO::XDMF::MeshPolicy::Transient);
  lfGrid.add(cellLabel,   IO::XDMF::Center::Cell);
  lfGrid.add(phaseMoment, IO::XDMF::Center::Cell);
  lfGrid.add(sigmaKgf,    IO::XDMF::Center::Cell);
  lfGrid.add(phiGf,       IO::XDMF::Center::Node);
  lfGrid.add(sLF,         IO::XDMF::Center::Node);
  lfGrid.add(u,           IO::XDMF::Center::Node);

  auto hfGrid = xdmf.grid("HF");
  hfGrid.setMesh(moved, IO::XDMF::MeshPolicy::Transient);
  hfGrid.add(cellLabelHF, IO::XDMF::Center::Cell);
  hfGrid.add(jKgf,        IO::XDMF::Center::Cell);
  hfGrid.add(qShape,      IO::XDMF::Center::Cell);
  hfGrid.add(phiMoved,    IO::XDMF::Center::Node);

  // -------------------------------------------------------------------------
  // Frame loop.
  // -------------------------------------------------------------------------
  std::cout << "Wavy-circle SDR sweep on " << n << "x" << n
            << " unit-square mesh, " << nFrames << " frames\n";
  std::cout << "  R0=" << R0 << "  amp=" << amp << "  k=" << kLobes
            << "  orbit R=" << orbitR << '\n';

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
    levelSet.phase = angle;  // co-rotate the lobes for visual variety.

    std::cout << "\n--- Frame " << std::setw(2) << frame
              << " : c=(" << std::fixed << std::setprecision(4)
              << levelSet.cx << ", " << levelSet.cy << ")"
              << "  phase=" << std::setprecision(3) << angle << " rad\n";

    // Reset attributes so the new classification owns them.
    clearXDMFRegionAttributes(mesh);

    // Re-tag static boundary (clearXDMFRegionAttributes wipes it too).
    for (auto faceIt = mesh.getBoundary(); faceIt; ++faceIt)
      mesh.setAttribute({mesh.getDimension() - 1, faceIt->getIndex()},
                        boundaryAttribute);

    // ---- Stage 1: phase moments + classification ----
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
      const auto& incident = mesh.getConnectivity().getIncidence({1, 2}, facet);
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

    // ---- Stage 2: FMM signed distance s_h^LF ----
    sLF = Real(0);
    Distance::Eikonal<ScalarP1, Math::Vector<Real>>(sLF)
      .setInterface(interfaceAttribute)
      .setInterior(interiorAttribute)
      .solve()
      .sign();

    // ---- Stage 3: smooth-band measure for the SDR normalisation ----
    Real domainMeasure = 0;
    Real weightedBandMeasure = 0;
    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const auto& cell = *cellIt;
      for (const auto& bary : TriangleBarycentricQuadrature)
      {
        Math::SpatialPoint rc(2);
        rc(0) = bary[1];
        rc(1) = bary[2];
        const Geometry::Point pt(cell, rc);
        Math::SpatialMatrix<Real> J;
        cell.getTransformation().jacobian(J, rc);
        const Real triangleArea = Real(0.5) * std::abs(J.determinant());
        const Real w =
          triangleArea / static_cast<Real>(TriangleBarycentricQuadrature.size());
        domainMeasure += w;
        const Real s = sLF.getValue(pt);
        const Real W = std::exp(-s * s / (2 * deltaW * deltaW));
        weightedBandMeasure += W * w;
      }
    }
    if (weightedBandMeasure <= 0)
      throw std::runtime_error("Empty smooth band: M_w = 0 at frame "
                               + std::to_string(frame));

    // ---- Stage 4a: Initial guess for u (strategy-selected) ----

    SDRParameters params;
    params.rhoS = 1;
    params.deltaW = deltaW;
    params.hRef = h;
    params.normalizer = Real(1) / (weightedBandMeasure * h * h);

    RealFunction<Real> barrierGamma(Real(1e-1));
    RealFunction<Real> barrierBeta(Real(1e-2));
    BarrierParameters barrierParams;
    barrierParams.jMin  = 1e-8;
    barrierParams.jSafe = 1e-3;
    barrierParams.domainMeasure = domainMeasure;

    auto cellGeom = precomputeCellGeometry(mesh);
    auto& cellCache  = cellGeom.first;
    auto& geomToLocal = cellGeom.second;

    SignedDistanceRegistration sdr(du, v, u);
    RealFunction phiFn(
        [&](const Geometry::Point& p) -> Real
        {
          const auto& c = p.getPhysicalCoordinates();
          return levelSet.phi(vec2(c(0), c(1)));
        });
    AnalyticVectorFunction gradFn(
        [&](const Geometry::Point& p) -> Math::SpatialVector<Real>
        {
          const auto& c = p.getPhysicalCoordinates();
          return levelSet.grad(vec2(c(0), c(1)));
        },
        /*dimension=*/2);
    AnalyticMatrixFunction hessFn(
        [&](const Geometry::Point& p) -> Math::SpatialMatrix<Real>
        {
          const auto& c = p.getPhysicalCoordinates();
          return levelSet.hess(vec2(c(0), c(1)));
        },
        /*rows=*/2, /*cols=*/2);

    JacobianAdmissibilityBarrier barrier(
        du, v, u, cellCache, geomToLocal);

    // ---- Hilbert-lift initial guess (closed-form Riesz lift in V₀) ----
    //
    // Solve for u_H ∈ V₀ such that
    //
    //     a(u_H, v) = -⟨δE_SDR(0), v⟩      ∀ v ∈ V₀,
    //
    // with a(u, v) = ∫ ∇u : ∇v dx and Dirichlet zero on ∂Ω. The right-
    // hand side is the SDR variational derivative AT u = 0:
    //
    //     ⟨δE_SDR(0), v⟩
    //       = ρ_s · normalizer · ∫_Ω w_δ(s_h^LF) · (φ - s_h^LF) · ∇φ · v dx.
    //
    // By Lax–Milgram on (V₀, a), this defines u_H uniquely; it is the
    // gradient of E_SDR at u = 0 in the H¹₀ topology (i.e. a Sobolev
    // gradient / Hilbert extension), and is the textbook clean initial
    // guess. See `body.tex § sec:initial-guess`.
    auto hilbertLiftToInit = [&]() {
      TrialFunction duH(vectorFes);
      TestFunction  vH(vectorFes);

      // RHS scalar: ρ_s · w_δ(s_h^LF) · (φ - s_h^LF).
      //
      // The SDR-normaliser 1/(M_w·h²) is NOT included here. It enters
      // the variational Newton problem to balance the rank-1 SDR block
      // against the shape/floor blocks at quadrature points; but the
      // Hilbert lift's LHS is unscaled (harmonic or elasticity), so
      // applying the normaliser to the RHS alone would inflate u_H by
      // ~1/h². The un-normalised RHS gives u_H at the natural
      // displacement scale.
      RealFunction rhsScalar(
          [&](const Geometry::Point& p) -> Real {
            const auto& c = p.getPhysicalCoordinates();
            const Real phiVal = levelSet.phi(vec2(c(0), c(1)));
            const Real sVal   = sLF.getValue(p);
            const Real wt     =
              std::exp(-sVal * sVal / (2 * deltaW * deltaW));
            return params.rhoS * wt * (phiVal - sVal);
          });

      // For the ShapeHessian variant the barrier integrator evaluates
      // δ²(E_shape + E_floor) at the current state u. We want that
      // tangent at u = 0, so zero u before assembly. The other variants
      // are linear forms and don't read u; zeroing it doesn't affect
      // them. u is then overwritten with u_H below regardless.
      u.getData().setZero();

      Problem hilbert(duH, vH);
      switch (kHilbertBilinearForm)
      {
        case HilbertBilinearForm::Harmonic:
          hilbert =
              Integral(Jacobian(duH), Jacobian(vH))
            + Integral(rhsScalar * gradFn, vH)
            + DirichletBC(duH, zero).on(boundaryAttribute);
          break;
        case HilbertBilinearForm::Elasticity:
          hilbert =
              LinearElasticityIntegral(duH, vH)(
                  kHilbertElasticityLambda,
                  kHilbertElasticityMu)
            + Integral(rhsScalar * gradFn, vH)
            + DirichletBC(duH, zero).on(boundaryAttribute);
          break;
        case HilbertBilinearForm::ShapeHessian:
        {
          JacobianAdmissibilityBarrier shapeBarrier(
              duH, vH, u, cellCache, geomToLocal);
          hilbert =
              shapeBarrier.Tangent(
                  barrierGamma, barrierBeta, barrierParams)
            + Integral(rhsScalar * gradFn, vH)
            + DirichletBC(duH, zero).on(boundaryAttribute);
          break;
        }
      }
      Solver::SparseLU(hilbert).solve();
      u.getData() = duH.getSolution().getData();
    };

    // ---- Newton solve with the PSD-projected full-Newton SDR tangent ----
    //
    // Why PSD-projected and not plain GN: pure Gauss-Newton drops the
    // second-order term `r * hess(phi)` from the SDR tangent. Whenever
    // r* != 0 at the discrete minimum (always the case here, since the
    // wavy contour is not exactly fittable by P1 on a uniform-grid
    // triangulation), the GN iteration is non-contracting near the
    // minimum: spectral_radius(I - K_GN^{-1} K_full) > 0. Damped GN
    // (damping < 1) reduces that spectral radius but produces a
    // residual floor proportional to r* and oscillates there.
    //
    // The principled fix is to put `r * hess(phi)` back IN, but clamped
    // to its PSD part per quadrature point. The PSD projection:
    //
    //   - keeps the contribution wherever it is positive-semidefinite,
    //     restoring quadratic Newton convergence in those directions;
    //   - zeros out the indefinite eigenspace, where raw full Newton
    //     would otherwise produce indefinite steps that flip cells.
    //
    // Concretely we use `sdr.TangentPSDProjected(phi, grad, hess, sLF,
    // params)` with damping = 1.0 (no static damping needed because
    // the tangent is SPD by construction).
    //
    // Per-iteration diagnostics: dumped when kVerbose, kept quiet
    // otherwise so the summary stays readable.
    constexpr bool kVerbose = false;
    barrierMinJ().store(std::numeric_limits<Real>::infinity(),
                        std::memory_order_relaxed);
    barrierFloorActiveCount().exchange(0);
    barrierInadmissibleCount().exchange(0);

    RealFunction<Real> harmEpsFn(kHarmonicEps);
    RealFunction<Real> tikEpsFn(kTikhonovEps);

    Problem newton(du, v);
    newton =
        sdr.TangentPSDProjected(phiFn, gradFn, hessFn, sLF, params)
      + sdr.Residual(phiFn, gradFn, sLF, params)
      + barrier.Tangent(barrierGamma, barrierBeta, barrierParams)
      + barrier.Residual(barrierGamma, barrierBeta, barrierParams)
      + Integral(harmEpsFn * Jacobian(du), Jacobian(v))
      + Integral(tikEpsFn  * du,            v          )
      + DirichletBC(du, zero).on(boundaryAttribute);

    // Best-iterate tracker: PSD-projected Newton contracts super-linearly
    // far from the minimum but degenerates to GN-with-drift once
    // r * hess(phi) becomes small at the minimum. The cleanest exit is
    // therefore to remember the iterate with the smallest residual seen
    // so far, and roll back to it after the solver returns. We also stop
    // early as soon as the residual fails to decrease for `kStallPatience`
    // consecutive iterations — that signals we have hit the local
    // minimum and any further iteration is drift, not progress.
    Math::Vector<Real> uBest = u.getData();
    Real residualBest = std::numeric_limits<Real>::infinity();
    std::size_t stallCount = 0;
    bool earlyExit = false;
    constexpr std::size_t kStallPatience = 3;

    // ------------------------------------------------------------------
    // Outer Newton loop.
    //
    //   for k = 0, 1, ...
    //     NewtonSolver does one step with damping = 1.0.
    //     p = u_after - u_before          (the un-damped Newton direction)
    //     restore u and re-apply via admissibility-aware backtracking
    //     line search.
    //
    // The line search and admissibility evaluator live in the support
    // header (LevelSetSDRWavyCircleSweepDiagnostics.h); the example
    // only orchestrates.
    // ------------------------------------------------------------------
    constexpr Real        kDamping        = 0.5;   ///< static damping when line search disabled
    constexpr Real        kTau            = 0.5;
    constexpr std::size_t kMaxIt          = 40;
    constexpr Real        kAbsTol         = 1e-8;
    constexpr Real        kRelTol         = 1e-7;

    auto evaluator = [&](const Math::Vector<Real>& uTry)
    {
      return Diag::evaluateAdmissibilityP1(
          uTry, cellCache, vectorFes, barrierParams.jMin);
    };

    Solver::SparseLU linearSolver(newton);
    Solver::NewtonSolver newtonSolver(linearSolver);
    newtonSolver
      .setMaxIterations(1)
      .setDampingFactor(Real(1))
      .setAbsoluteTolerance(Real(0))   ///< outer loop owns convergence
      .setRelativeTolerance(Real(0));

    // -------------------------------------------------------------------
    // Newton solve with retry strategy: warm-start (already loaded above)
    // first, fall back to cold-start if the warm attempt fails to meet
    // the convergence criterion. Frame 0 only ever runs one attempt
    // because it is already cold-started.
    // -------------------------------------------------------------------
    Real        r0      = -1;
    std::size_t itCount = 0;
    Real        interfaceFit = 0;
    bool        convergedFlag    = false;
    const char* convergedReason  = "no";
    int         attemptUsed      = 0;

    for (int attempt = 0; attempt < 2; ++attempt)
    {
      // Attempt 0: strategy-selected initial guess.
      // Attempt 1: cold-start retry (u = 0).
      if (attempt == 0)
      {
        switch (kInitialGuessStrategy)
        {
          case InitialGuessStrategy::Cold:
            u.getData().setZero();
            break;
          case InitialGuessStrategy::Warm:
            if (frame == 0) u.getData().setZero();
            // else: leave u with previous frame's value.
            break;
          case InitialGuessStrategy::Hilbert:
            hilbertLiftToInit();
            break;
        }
      }
      else
      {
        u.getData().setZero();
      }
      uBest = u.getData();
      residualBest = std::numeric_limits<Real>::infinity();
      stallCount = 0;
      earlyExit  = false;
      attemptUsed = attempt;
      r0      = -1;
      itCount = 0;

    for (std::size_t it = 0; it < kMaxIt; ++it)
    {
      itCount = it + 1;
      const Math::Vector<Real> uOld = u.getData();

      barrierInadmissibleCount().exchange(0);
      barrierFloorActiveCount().exchange(0);
      barrierMinJ().store(std::numeric_limits<Real>::infinity(),
                          std::memory_order_relaxed);

      // One Newton step: assemble + solve + un-damped apply.
      newtonSolver.solve(u);
      const auto step = newtonSolver.getReport();
      const Real R = step.final_residual;

      if (R < residualBest)
      {
        residualBest = R;
        uBest = uOld;
        stallCount = 0;
      }
      else
      {
        ++stallCount;
        if (stallCount >= kStallPatience) earlyExit = true;
      }

      if (it == 0) r0 = R;
      if (R <= kAbsTol) break;
      if (r0 > 0 && R <= kRelTol * r0) break;

      // Reconstruct the un-damped Newton direction and back-substitute
      // through the admissibility line search.
      const Math::Vector<Real> p = u.getData() - uOld;
      u.getData() = uOld;

      Diag::LineSearchResult lsRes;
      Real alphaUsed = Real(0);
      if (kLineSearchEnabled)
      {
        // Tighten the line-search floor from j_min (catastrophic flip
        // protection only) to a safety multiple of j_safe (avoid the
        // floor barrier's active region). See `kLineSearchSafetyMargin`.
        const Real lsFloor =
          std::max(barrierParams.jMin,
                   kLineSearchSafetyMargin * barrierParams.jSafe);
        lsRes = Diag::admissibilityLineSearch(
            u.getData(), p, evaluator,
            lsFloor,
            kLineSearchAlphaInit, kLineSearchReduction, kLineSearchAlphaMin);
        alphaUsed = lsRes.alphaAccepted;
      }
      else
      {
        // Baseline static damping; admissibility checked for diagnostics
        // only (never used to reject).
        const auto adm1 = evaluator(u.getData() + Real(1) * p);
        lsRes.minJAtAlpha1              = adm1.minJ;
        lsRes.inadmissibleCountAtAlpha1 = adm1.inadmissibleCount;
        u.getData() += kDamping * p;
        const auto admUsed = evaluator(u.getData());
        lsRes.alphaAccepted             = kDamping;
        lsRes.minJAccepted              = admUsed.minJ;
        lsRes.inadmissibleCountAccepted = admUsed.inadmissibleCount;
        lsRes.succeeded                 = true;
        alphaUsed = kDamping;
      }
      const Real stepNorm = alphaUsed * p.norm();

      // Observation-only diagnostics.
      const Diag::AdmissibilityDiagnostics admDiag =
        Diag::computeAdmissibilityDiagnostics(
            uOld, p, cellCache, vectorFes, barrierParams.jMin, kTau);

      if (kVerbose)
      {
        const Real minJBar = barrierMinJ().load(std::memory_order_relaxed);
        const auto badN    = barrierInadmissibleCount().load(std::memory_order_relaxed);
        const auto floorN  = barrierFloorActiveCount().load(std::memory_order_relaxed);

        std::cout << "      it=" << std::setw(2) << it
                  << "  ||R||="  << std::scientific << std::setprecision(3) << R
                  << "  step="   << stepNorm << '\n';
        std::cout << "           min_j_before="    << std::setprecision(3) << admDiag.minJBefore
                  << "  min_j_after_full=" << admDiag.minJAfterFull
                  << "  inadm_after="      << admDiag.inadmAfter
                  << "  min_pred_margin="  << admDiag.minPredMargin
                  << "  pred_inadm="       << admDiag.predInadm;
        if (std::isfinite(minJBar)) std::cout << "  min_j(barrier)=" << minJBar;
        if (floorN > 0)             std::cout << "  floor_cells=" << floorN;
        if (badN > 0)               std::cout << "  singular=" << badN;
        std::cout << '\n';
        std::cout << "           line_search:"
                  << " alpha="                  << lsRes.alphaAccepted
                  << "  backtracks="            << lsRes.backtracks
                  << "  min_j@alpha=1="         << lsRes.minJAtAlpha1
                  << "  inadm@alpha=1="         << lsRes.inadmissibleCountAtAlpha1
                  << "  min_j_accepted="        << lsRes.minJAccepted
                  << "  inadm_accepted="        << lsRes.inadmissibleCountAccepted
                  << (lsRes.succeeded ? "" : "  [FAILED]")
                  << '\n';
      }

      if (kLineSearchEnabled && !lsRes.succeeded) break;
      if (stepNorm <= Real(0))                    break;
      if (earlyExit)                              break;
    }

    // Roll back to the best iterate seen.
    if (residualBest < std::numeric_limits<Real>::infinity())
      u.getData() = uBest;
    (void) earlyExit;

    // Compute the interface fit RMS at the rolled-back u, so it can
    // feed the convergence decision below.
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
        const Real val  = Real(0.5) * (phiA + phiB);
        const Real len  = (b - a).norm();
        interfacePhi += val * val * len;
        interfaceLen += len;
      }
      interfaceFit =
        std::sqrt(interfacePhi / std::max(interfaceLen, Real(1e-30)));
    }

    const bool residualOK =
         residualBest < kResidualAbsTol
      || (r0 > 0 && residualBest < kResidualRelTol * r0);
    const bool geometryOK = interfaceFit < kInterfaceFitTol;
    switch (kConvergenceMode)
    {
      case ConvergenceMode::Residual:
        convergedFlag   = residualOK;
        convergedReason = residualOK ? "yes (residual)" : "no";
        break;
      case ConvergenceMode::Geometry:
        convergedFlag   = geometryOK;
        convergedReason = geometryOK ? "yes (geometry)" : "no";
        break;
      case ConvergenceMode::Either:
        convergedFlag = geometryOK || residualOK;
        convergedReason =
            geometryOK && residualOK ? "yes (both)"
          : geometryOK               ? "yes (geometry)"
          : residualOK               ? "yes (residual)"
                                     : "no";
        break;
    }

      // Stop the retry loop on success or after the cold fallback.
      if (convergedFlag || attempt == 1) break;
    }
    finalFitPerFrame.push_back(interfaceFit);

    struct AggregateReport
    {
      std::size_t iterations;
      Real        final_residual;
      Real        interface_fit;
      bool        converged;
    } report{
      itCount,
      residualBest,
      interfaceFit,
      convergedFlag
    };

    // ---- Stage 5a: LF cell-P0 fields + nodal phi ----
    {
      const std::size_t d = mesh.getDimension();
      for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
      {
        const Index cellIdx = cellIt->getIndex();
        const std::size_t local = cellToLocal.at(cellIdx);
        const Index dof = p0Fes.getGlobalIndex({d, cellIdx}, 0);
        cellLabel.getData()(dof)   =
          static_cast<Real>(classified.labels[local]);
        phaseMoment.getData()(dof) = cellMoments[local].moment;
        sigmaKgf.getData()(dof)    =
          static_cast<Real>(cellCache[local].sigmaK);
      }
    }
    phiGf = [&](const Geometry::Point& p) -> Real
    {
      const auto& X = p.getCoordinates();
      return levelSet.phi(vec2(X(0), X(1)));
    };

    // ---- Stage 5b: build the moved (HF) mesh and its HF fields ----
    //
    // The HF mesh shares combinatorics with `mesh`; only vertex
    // coordinates are pushed by u. Region attributes (interior /
    // exterior / interface / boundary) carry over so the HF grid carries
    // the same classification ParaView sees on LF.
    {
      const auto& uFes = u.getFiniteElementSpace();
      const auto& uData = u.getData();
      const Index vn = mesh.getVertexCount();
      for (Index vertex = 0; vertex < vn; ++vertex)
      {
        const Vec2 x = mesh.getVertexCoordinates(vertex);
        const auto& dofs = uFes.getDOFs(0, vertex);
        const Real ux = uData(dofs[0]);
        const Real uy = uData(dofs[1]);
        moved.setVertexCoordinates(vertex, vec2(x(0) + ux, x(1) + uy));
      }
      const std::size_t D = mesh.getDimension();
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
    }

    // Per-cell j_K^u and Q_shape on the moved mesh — same closed-form
    // expressions as the parent example (the σ_K branch is fixed on the
    // original mesh; the moved Jacobian is read from `movedCache`).
    auto [movedCache, movedToLocal] = precomputeCellGeometry(moved);
    {
      const std::size_t d = moved.getDimension();
      for (auto cellIt = moved.getCell(); cellIt; ++cellIt)
      {
        const Index cellIdx = cellIt->getIndex();
        const std::size_t local = cellToLocal.at(cellIdx);
        const std::size_t movedLocal = movedToLocal.at(cellIdx);
        const Index dof = p0FesMoved.getGlobalIndex({d, cellIdx}, 0);
        const auto& src = cellCache[local];
        const auto& dst = movedCache[movedLocal];
        const Real sigDetAu = static_cast<Real>(src.sigmaK) * dst.detAK;
        jKgf.getData()(dof) = sigDetAu / src.Jscale;
        const Real dExp = Real(2) / Real(d);
        qShape.getData()(dof) =
          dst.A.squaredNorm()
            / (static_cast<Real>(d) * std::pow(sigDetAu, dExp));
        cellLabelHF.getData()(dof) =
          static_cast<Real>(classified.labels[local]);
      }
    }
    phiMoved = [&](const Geometry::Point& p) -> Real
    {
      const auto& X = p.getCoordinates();
      return levelSet.phi(vec2(X(0), X(1)));
    };

    const char* strategyName =
        kInitialGuessStrategy == InitialGuessStrategy::Cold    ? "cold"
      : kInitialGuessStrategy == InitialGuessStrategy::Warm    ? "warm"
      :                                                          "hilbert";
    std::cout << "    Newton it=" << report.iterations
              << "  ||R||=" << std::scientific << std::setprecision(3)
              << report.final_residual
              << "  fit=" << std::setprecision(3) << interfaceFit
              << "  init=" << (attemptUsed == 0 ? strategyName : "cold")
              << "  converged=" << convergedReason << '\n';
    if (report.converged)
      ++framesConverged;

    // ---- Stage 6: append XDMF snapshot ----
    xdmf.write(t).flush();
  }

  xdmf.close();

  // -------------------------------------------------------------------------
  // Summary.
  // -------------------------------------------------------------------------
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
    for (Real x : finalFitPerFrame) fitMean += x;
    fitMean /= static_cast<Real>(finalFitPerFrame.size());
    std::cout << "  ||phi(X+u)||_RMS  min=" << std::scientific
              << std::setprecision(3) << fitMin
              << "  mean=" << fitMean
              << "  max=" << fitMax << '\n';
  }

  return framesConverged == nFrames ? 0 : 1;
}
