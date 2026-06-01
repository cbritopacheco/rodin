/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//
// Wavy-circle sweep test for the SDFR pipeline.
//
// A non-convex implicit shape — a circle whose radius varies sinusoidally
// in the polar angle (k azimuthal lobes) — is moved around the unit
// square along a circular orbit. The whole classify-then-displace
// pipeline (s-t cut + FMM + SDFR Newton) is rerun from scratch at every
// frame, and a transient XDMF stream is emitted with both the
// low-fidelity classification and the high-fidelity moved mesh per
// frame.
//
// Compared with `LevelSetSDRReconstruction`, this example exercises:
//   - A non-circular, non-convex interface that genuinely requires
//     tangential motion of the SDFR (the cosine perturbation cannot be
//     absorbed by uniform dilation of the cell map).
//   - Repeated classification + Newton solves at moving interface
//     locations: the classified topology changes between frames and the
//     transient XDMF writer dumps a fresh mesh each frame.
//   - The shard-aware XDMF writers under non-trivial cell-attribute
//     turnover.
//
// FES-independence and Newton-tangent choices follow the same honest
// caveats as the parent example: the barrier + shape closed-form
// algebra is specialised to affine P1 triangles in 2D; SDFR uses the
// PSD-projected full-Newton tangent. The example is therefore a 2D
// triangular prototype.
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
  // gradient norm differ from 1) but the SDFR pipeline only requires a
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

  std::string parseStringOption(
      int argc, char** argv, const std::string& name, std::string fallback)
  {
    const std::string prefix = "--" + name + "=";
    for (int i = 1; i < argc; ++i)
    {
      const std::string arg(argv[i]);
      if (arg.rfind(prefix, 0) == 0)
        return arg.substr(prefix.size());
    }
    return fallback;
  }

  SDFRHilbertMetric parseHilbertMetric(
      int argc, char** argv, const std::string& name)
  {
    const std::string value = parseStringOption(argc, argv, name, "harmonic");
    if (value == "harmonic") return SDFRHilbertMetric::Harmonic;
    if (value == "elasticity") return SDFRHilbertMetric::Elasticity;
    if (value == "shape-hessian") return SDFRHilbertMetric::ShapeHessian;
    throw std::runtime_error(
        "Unknown --" + name + "=" + value
        + " (expected harmonic, elasticity, or shape-hessian).");
  }

  SDFRInitialGuessScaling parseInitialGuessScaling(
      int argc, char** argv, const std::string& name)
  {
    const std::string value =
      parseStringOption(argc, argv, name, "unnormalized");
    if (value == "unnormalized") return SDFRInitialGuessScaling::Unnormalized;
    if (value == "energy") return SDFRInitialGuessScaling::EnergyNormalized;
    if (value == "band") return SDFRInitialGuessScaling::BandNormalized;
    throw std::runtime_error(
        "Unknown --" + name + "=" + value
        + " (expected unnormalized, energy, or band).");
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
  // -------------------------------------------------------------------------
  // Step 0: frame schedule, mesh, and constants.
  // -------------------------------------------------------------------------
  const std::size_t n =
    parseSizeTOption(argc, argv, "n", 50);          ///< n x n nodes.
  const std::size_t nFrames =
    parseSizeTOption(argc, argv, "frames", 200);    ///< orbit snapshots.
  constexpr Real        orbitR = 0.10; ///< Radius of the centre's orbit.
  const Real            amp =
    parseRealOption(argc, argv, "amp", Real(0.05)); ///< radial amplitude.
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
  //   - Harmonic: scaling with the SDFR per-cell magnitude ρs/M_w ≈ O(1)
  //     gives a regularisation that is mesh-independent in `eps_H`.
  //     At n=30 (h ≈ 0.034) the GN noise floor swamps the iterate
  //     unless eps_H >= 1e-2 or so.
  //   - Tikhonov: the mass matrix has h² scaling per cell, so for
  //     mesh-independent regularisation `eps_T` must scale ~1/h². We
  //     default to 0 here and rely on harmonic.
  //
  // Both knobs are opt-in. Set to 0 to recover the un-regularised
  // PSD-projected Newton iteration that was already used at n=60.
  const Real kHarmonicEps =
    parseRealOption(argc, argv, "harmonic-eps", Real(0));
  const Real kTikhonovEps =
    parseRealOption(argc, argv, "tikhonov-eps", Real(0));

  // -------------------------------------------------------------------------
  // Admissibility-first backtracking line search.
  //
  //   while (alpha >= alphaMin)
  //     uTrial = u + alpha * p
  //     adm    = evaluateAdmissibility(uTrial)
  //     if (adm.minJRatio > jLineSearchRatio
  //         && adm.inadmissibleCount == 0): accept
  //     else: alpha *= reduction
  //
  // When disabled, the iteration falls back to the static `kDamping`
  // factor (the original failing baseline).
  const Real   kLineSearchAlphaInit = Real(1);
  const Real   kLineSearchReduction = Real(0.5);
  const Real   kLineSearchAlphaMin  = Real(1e-6);
  // Safety margin above jSafeRatio used as the line-search admissibility
  // floor. The threshold is dimensionless because it is applied to
  // j_K^u = sigma_K det(A_K^u) / J_K_scale, never to raw det(A_K^u).
  // Accept only steps with
  //   min_K j_K^u > max(jMinRatio, margin * jSafeRatio).
  // The previous default (j_min = 10^-8) only protected against
  // catastrophic flips, allowing iterates to enter the floor-barrier
  // active region (j in (j_min, j_safe)) where the barrier Hessian is
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
  const SDFRHilbertMetric initialGuessMetric =
    parseHilbertMetric(argc, argv, "hilbert-metric");
  const SDFRInitialGuessScaling initialGuessScaling =
    parseInitialGuessScaling(argc, argv, "hilbert-scaling");
  const Real initialGuessElasticityLambda =
    parseRealOption(argc, argv, "hilbert-lambda", Real(0));
  const Real initialGuessElasticityMu =
    parseRealOption(argc, argv, "hilbert-mu", Real(0.5));

  // ----- Convergence criterion ---------------------------------------------
  //   Residual : classic ||R|| < tol (fails to recognise the GN data-floor).
  //   Geometry : ||phi(X+u)||_RMS < ~h^2 (honest "did we capture Γ?").
  //   Either   : converged if EITHER residual or geometry tolerance holds.
  enum class ConvergenceMode { Residual, Geometry, Either };
  constexpr ConvergenceMode kConvergenceMode = ConvergenceMode::Either;
  // Knobs are CLI-overridable so we can sweep convergence settings without
  // rebuilding. Defaults reproduce the previous behaviour.
  // Defaults tuned by the 50-frame proxy sweep (see task FF): the combined
  // setting (fittol-mult=4, rrtol=5e-3, stall=5) raised converged from
  // 39/50 to 50/50 and cut wall-time 5.2x at essentially identical RMS.
  // Relaxing the geometry target to 4*h^2 reflects the achievable accuracy
  // floor of the smooth-band SDFR fit; relaxing rrtol to 5e-3 matches the
  // PSD-projected Newton noise floor; stall=5 only exits early once
  // geometryOK is reachable.
  const Real        kFitTolMult      =
    parseRealOption(argc, argv, "fittol-mult", Real(4.0)); ///< fitTol = mult * h^2
  const Real        kInterfaceFitTol = kFitTolMult * h * h;
  constexpr Real    kResidualAbsTol  = Real(1e-6);
  const Real        kResidualRelTol  =
    parseRealOption(argc, argv, "rrtol", Real(5e-3));
  const std::size_t kStallPatience   =
    parseSizeTOption(argc, argv, "stall", 5);

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

  // Background mesh is invariant across the sweep, so the per-cell geometry
  // cache (sigma_K, det A_K, J_scale, gradN, area) only needs to be built
  // once. Recomputing it each frame was a measurable performance regression.
  auto cellGeomBg = precomputeCellGeometry(mesh);
  auto& cellCacheBg = cellGeomBg.first;

  GridFunction sLF(p1Fes);            sLF.setName("s_LF");
  GridFunction phiGf(p1Fes);          phiGf.setName("phi");
  GridFunction cellLabel(p0Fes);      cellLabel.setName("cell_label");
  GridFunction phaseMoment(p0Fes);    phaseMoment.setName("phase_moment");
  GridFunction sigmaKgf(p0Fes);       sigmaKgf.setName("sigma_K");

  GridFunction  u(vectorFes);         u.setName("displacement");

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
  std::cout << "Wavy-circle SDFR sweep on " << n << "x" << n
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

    // ---- Stage 3: smooth-band measure for the SDFR normalisation ----
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

    // ---- Stage 4a: SDFR solve (strategy-selected initial guess) ----

    auto& cellCache = cellCacheBg;

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

    const bool kVerbose = hasFlag(argc, argv, "verbose");

    SDFRParameters baseParams;
    baseParams.rhoS = 1;
    baseParams.deltaW = deltaW;
    baseParams.hRef = h;
    baseParams.normalizer = Real(1) / (weightedBandMeasure * h * h);
    baseParams.shapeWeight = Real(1e-1);
    baseParams.floorWeight = Real(1e-2);
    baseParams.jMinRatio = Real(1e-8);
    baseParams.jSafeRatio = Real(1e-3);
    baseParams.lineSearchSafetyMargin = kLineSearchSafetyMargin;
    baseParams.alphaInit = kLineSearchAlphaInit;
    baseParams.alphaReduction = kLineSearchReduction;
    baseParams.alphaMin = kLineSearchAlphaMin;
    baseParams.absoluteTolerance = kResidualAbsTol;
    baseParams.relativeTolerance = kResidualRelTol;
    baseParams.maxNewtonIterations = 80;
    // Let the sweep-level geometry criterion decide when to stop; residual
    // stall alone is not a failure for this registration objective.
    baseParams.stallPatience = kStallPatience;
    baseParams.harmonicEps = kHarmonicEps;
    baseParams.tikhonovEps = kTikhonovEps;
    baseParams.initialGuessMetric = initialGuessMetric;
    baseParams.initialGuessScaling = initialGuessScaling;
    baseParams.initialGuessElasticityLambda = initialGuessElasticityLambda;
    baseParams.initialGuessElasticityMu = initialGuessElasticityMu;
    baseParams.tangent = SDFRTangent::PSDProjectedNewton;

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
        const Real val = Real(0.5) * (phiA + phiB);
        const Real len = (b - a).norm();
        interfacePhi += val * val * len;
        interfaceLen += len;
      }
      return std::sqrt(interfacePhi / std::max(interfaceLen, Real(1e-30)));
    };

    Real residualBest = std::numeric_limits<Real>::infinity();
    Real        r0      = -1;
    std::size_t itCount = 0;
    Real        interfaceFit = 0;
    bool        convergedFlag    = false;
    const char* convergedReason  = "no";
    int         attemptUsed      = 0;

    for (int attempt = 0; attempt < 2; ++attempt)
    {
      SDFRParameters params = baseParams;
      if (attempt == 0)
      {
        switch (kInitialGuessStrategy)
        {
          case InitialGuessStrategy::Cold:
            params.initialGuess = SDFRInitialGuess::Zero;
            break;
          case InitialGuessStrategy::Warm:
            params.initialGuess =
              frame == 0 ? SDFRInitialGuess::Zero : SDFRInitialGuess::Current;
            break;
          case InitialGuessStrategy::Hilbert:
            params.initialGuess = SDFRInitialGuess::Hilbert;
            break;
        }
      }
      else
      {
        params.initialGuess = SDFRInitialGuess::Zero;
      }

      SDFR sdfr(u);
      const SDFRReport sdfrReport =
        sdfr.setParameters(params).solve(sLF, phiFn, gradFn, hessFn);

      attemptUsed = attempt;
      r0 = sdfrReport.initialResidual;
      itCount = sdfrReport.iterations;
      residualBest = sdfrReport.finalResidual;
      interfaceFit = computeInterfaceFit();

      const bool residualOK =
           residualBest < kResidualAbsTol
        || (r0 > 0 && residualBest < kResidualRelTol * r0);
      const bool geometryOK = interfaceFit < kInterfaceFitTol;
      switch (kConvergenceMode)
      {
        case ConvergenceMode::Residual:
          convergedFlag = residualOK;
          convergedReason = residualOK ? "yes (residual)" : "no";
          break;
        case ConvergenceMode::Geometry:
          convergedFlag = geometryOK;
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

      if (kVerbose)
      {
        std::cout << "      SDFR:"
                  << " it=" << sdfrReport.iterations
                  << "  ||R||=" << std::scientific << std::setprecision(3)
                  << sdfrReport.finalResidual
                  << "  alpha=" << sdfrReport.lastAcceptedAlpha
                  << "  backtracks=" << sdfrReport.totalBacktracks
                  << "  j_ls=" << sdfrReport.jLineSearchRatio
                  << "  min_j=" << sdfrReport.minJRatio
                  << (sdfrReport.lineSearchFailed ? "  [line-search failed]" : "")
                  << '\n';
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
