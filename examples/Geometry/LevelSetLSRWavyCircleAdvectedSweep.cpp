/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//
// ADVECTED wavy-circle sweep test for the LSR push pipeline.
//
// Same orbit + lobed-circle geometry as `LevelSetLSRWavyCircleSweep`, but
// the data level set phi is NOT supplied analytically per frame. Instead:
//
//   Frame 0:  phi_h is built from the analytic wavy circle. Runtime options
//             select analytic sampling, screened-Poisson reconstruction, or
//             FMM reconstruction through --phi-init.
//
//   Frame k > 0:  phi_h is advected by ONE semi-Lagrangian step with a
//                 prescribed velocity v(x) using `Advection::Lagrangian`.
//
// The LSR push solve at every frame consumes phi_h (discrete P1) and
// Grad(phi_h) (discrete piecewise-constant) -- NO analytic phi enters the
// Newton tangent or residual. The shape barrier, the alpha-predictor,
// the line-search, and the H^1 Hilbert initial guess are unchanged from
// the analytic sweep.
//
// Diagnostic only: the analytic level set is still reconstructed each
// frame to compute a reference fit metric
//   ||phi_analytic(X + u_h(X))||_RMS  on the classified skeleton,
// alongside the discrete fit
//   ||phi_h(X + u_h(X))||_RMS         on the classified skeleton.
// The gap between the two reflects accumulated advection error in phi_h.
//
// Goal of this example: assess how the LSR fit degrades when phi_h drifts
// from its analytic reference under repeated advection -- i.e. whether
// the push + psi-reconstruction + alpha-predictor + line-search pipeline is
// stable when phi loses its SDF-like property and its zero level set stops
// matching the ground-truth geometry exactly. Runtime options include
//   --psi=screened|fmm|projected-phi-h1
//   --phi-redistance=off|screened-moved|fmm-moved|projected-phi-h1-moved
//
// Tangent. The discrete phi has no recoverable Hess; the LSR solve
// defaults to Gauss-Newton (--tangent=gn). A zero Hess function is
// passed as a placeholder so the existing `solve` overload still
// resolves; the Newton tangent path is not exercised here.
//
#include <Rodin/Adaptation.h>
#include <Rodin/Advection/Lagrangian.h>
#include <Rodin/Assembly.h>
#include <Rodin/Distance/Eikonal.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Math.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>
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
      const Real ux = uData(dofs[0]);
      const Real uy = uData(dofs[1]);
      moved.setVertexCoordinates(vertex, vec2(x(0) + ux, x(1) + uy));
    }

#ifdef RODIN_LSR_P2_DISPLACEMENT
    Variational::RealH1Element<2> geomFe(Polytope::Type::Triangle);
    Math::SpatialPoint X;
    const std::size_t D = mesh.getDimension();
    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const auto& cell = *cellIt;
      Geometry::PointCloud pm(2, geomFe.getCount());
      for (std::size_t a = 0; a < geomFe.getCount(); ++a)
      {
        const auto& rc = geomFe.getNode(a);
        cell.getTransformation().transform(X, rc);
        const Real ux = uData(uFes.getGlobalIndex({D, cell.getIndex()}, a * 2));
        const Real uy = uData(uFes.getGlobalIndex({D, cell.getIndex()}, a * 2 + 1));
        pm(0, a) = X(0) + ux;
        pm(1, a) = X(1) + uy;
      }
      moved.setPolytopeTransformation(
          {D, cell.getIndex()},
          new Geometry::ParametricTransformation<Variational::RealH1Element<2>>(
            std::move(pm), geomFe));
    }
#endif
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
  // gradient norm differ from 1) but the LSR pipeline only requires a
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
      const LocalMesh& mesh, PhiFn&& phi, Real epsilon)
  {
    std::vector<CellMomentInfo> cells;
    cells.reserve(mesh.getCellCount());

    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const auto& cellPolytope = *cellIt;
      const auto& vertices = cellPolytope.getVertices();
      if (vertices.size() != 3)
        throw std::runtime_error(
            "LevelSetLSRWavyCircleSweep expects triangular cells.");

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
        moment += applyPhaseMomentMap(phi(xq), epsilon);
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

  LSRHilbertMetric parseHilbertMetric(
      int argc, char** argv, const std::string& name)
  {
    const std::string value = parseStringOption(argc, argv, name, "harmonic");
    if (value == "harmonic") return LSRHilbertMetric::Harmonic;
    if (value == "elasticity") return LSRHilbertMetric::Elasticity;
    if (value == "shape-hessian") return LSRHilbertMetric::ShapeHessian;
    throw std::runtime_error(
        "Unknown --" + name + "=" + value
        + " (expected harmonic, elasticity, or shape-hessian).");
  }

  LSRInitialGuessScaling parseInitialGuessScaling(
      int argc, char** argv, const std::string& name)
  {
    const std::string value =
      parseStringOption(argc, argv, name, "unnormalized");
    if (value == "unnormalized") return LSRInitialGuessScaling::Unnormalized;
    if (value == "energy") return LSRInitialGuessScaling::EnergyNormalized;
    if (value == "band") return LSRInitialGuessScaling::BandNormalized;
    throw std::runtime_error(
        "Unknown --" + name + "=" + value
        + " (expected unnormalized, energy, or band).");
  }

  LSRTangent parseTangentMode(
      int argc, char** argv, const std::string& name)
  {
    const std::string value = parseStringOption(argc, argv, name, "psd");
    if (value == "psd") return LSRTangent::PSDProjectedNewton;
    if (value == "newton") return LSRTangent::Newton;
    if (value == "gn") return LSRTangent::GaussNewton;
    throw std::runtime_error(
        "Unknown --" + name + "=" + value
        + " (expected psd, newton, or gn).");
  }

  // Prescribed velocity field driving the per-frame semi-Lagrangian
  // advection of phi_h.
  //
  //   Rotation  : v(x, y) = omega * ( -(y - 0.5),  x - 0.5 ).
  //               With omega = 2*pi the orbit period 1 sweeps a full
  //               rotation around the unit square's centre, matching the
  //               analytic level set's centre orbit AND phase rotation
  //               exactly (modulo discretisation).
  //   None      : v = 0. phi_h stays at its frame-0 reconstruction; the
  //               sweep reduces to a non-advection consistency check.
  enum class VelocityField { Rotation, None };

  VelocityField parseVelocityField(
      int argc, char** argv, const std::string& name)
  {
    const std::string value = parseStringOption(argc, argv, name, "rotation");
    if (value == "rotation") return VelocityField::Rotation;
    if (value == "none")     return VelocityField::None;
    throw std::runtime_error(
        "Unknown --" + name + "=" + value
        + " (expected rotation or none).");
  }

  enum class LevelSetRepresentative { Screened, Fmm, ProjectedPhiH1 };

  LevelSetRepresentative parseLevelSetRepresentative(
      int argc, char** argv, const std::string& name,
      const std::string& fallback)
  {
    const std::string value = parseStringOption(argc, argv, name, fallback);
    if (value == "screened") return LevelSetRepresentative::Screened;
    if (value == "fmm")      return LevelSetRepresentative::Fmm;
    if (value == "projected-phi-h1" || value == "projection-h1")
      return LevelSetRepresentative::ProjectedPhiH1;
    throw std::runtime_error(
        "Unknown --" + name + "=" + value
        + " (expected screened, fmm, or projected-phi-h1).");
  }

  enum class PhiInitialRepresentative { Analytic, Screened, Fmm };

  PhiInitialRepresentative parsePhiInitialRepresentative(
      int argc, char** argv, const std::string& name)
  {
    const std::string value = parseStringOption(argc, argv, name, "analytic");
    if (value == "analytic") return PhiInitialRepresentative::Analytic;
    if (value == "screened") return PhiInitialRepresentative::Screened;
    if (value == "fmm")      return PhiInitialRepresentative::Fmm;
    throw std::runtime_error(
        "Unknown --" + name + "=" + value
        + " (expected analytic, screened, or fmm).");
  }

  enum class PhiRedistance { Off, ScreenedMoved, FmmMoved, ProjectedPhiH1Moved };

  PhiRedistance parsePhiRedistance(
      int argc, char** argv, const std::string& name)
  {
    const std::string value = parseStringOption(argc, argv, name, "off");
    if (value == "off") return PhiRedistance::Off;
    if (value == "screened-moved") return PhiRedistance::ScreenedMoved;
    if (value == "fmm-moved") return PhiRedistance::FmmMoved;
    if (value == "projected-phi-h1-moved" || value == "projection-h1-moved")
      return PhiRedistance::ProjectedPhiH1Moved;
    throw std::runtime_error(
        "Unknown --" + name + "=" + value
        + " (expected off, screened-moved, fmm-moved, or projected-phi-h1-moved).");
  }

  // Cold    : u0 = 0 every frame.
  // Warm    : u0 = previous frame's converged u (zero on frame 0).
  // Hilbert : u0 = Riesz lift of -delta E_LSR(0) in H^1_0.
  enum class InitialGuessStrategy { Cold, Warm, Hilbert };

  InitialGuessStrategy parseInitialGuessStrategy(
      int argc, char** argv, const std::string& name)
  {
    const std::string value = parseStringOption(argc, argv, name, "hilbert");
    if (value == "cold") return InitialGuessStrategy::Cold;
    if (value == "warm") return InitialGuessStrategy::Warm;
    if (value == "hilbert") return InitialGuessStrategy::Hilbert;
    throw std::runtime_error(
        "Unknown --" + name + "=" + value
        + " (expected cold, warm, or hilbert).");
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
  std::cout << std::unitbuf;
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

  const Real h = Real(1) / static_cast<Real>(n - 1);
  const Real epsilon = 1.25 * h;
  const Real lambdaC = 0.008;
  const Real delta   = 1.75 * h;
  const Real deltaW  = 1.5 * delta;

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
  // Relative-distortion cap for the line search (identity-neutral):
  //   Q_rel(F) = ||F||^2 / (d * det(F)^(2/d)),  F = I + grad u_h.
  // Default 10 leaves the historical effective behaviour untouched.
  const Real   kQRelMax =
    parseRealOption(argc, argv, "qrel-max", Real(10));
  const std::size_t kLSRQuadratureOrder =
    parseSizeTOption(argc, argv, "lsr-quad-order", 0);
  const Real kH1RegularizationWeight =
    parseRealOption(argc, argv, "h1-reg", Real(0));
  // BEST-QUALITY profile default: γ_shape = √h (dimensionless;
  // multiplied by normalizer · h_ref² inside LSR.h).
  const Real kShapeWeight =
    parseRealOption(argc, argv, "gamma", std::sqrt(h));

  // ----- Initial guess strategy --------------------------------------------
  //   Cold    : u₀ = 0
  //   Warm    : u₀ = previous frame's converged u
  //   Hilbert : u₀ = Riesz lift of −δE_LSR(0) in the H¹₀ inner product
  //             (a(u₀, v) = −⟨δE_LSR(0), v⟩  on V₀)
  // Selected at runtime via --init=cold|warm|hilbert (default: hilbert).
  const InitialGuessStrategy kInitialGuessStrategy =
      parseInitialGuessStrategy(argc, argv, "init");
  const LSRHilbertMetric initialGuessMetric =
    parseHilbertMetric(argc, argv, "hilbert-metric");
  const LSRInitialGuessScaling initialGuessScaling =
    parseInitialGuessScaling(argc, argv, "hilbert-scaling");
  const Real initialGuessElasticityLambda =
    parseRealOption(argc, argv, "hilbert-lambda", Real(0));
  const Real initialGuessElasticityMu =
    parseRealOption(argc, argv, "hilbert-mu", Real(0.5));
  const LSRTangent tangentMode =
    parseTangentMode(argc, argv, "tangent");
  const LevelSetRepresentative psiRepresentative =
    // ψ default = FMM for robustness; survives the discrete-φ
    // dissipation noise per frame without depending on φ-regularity.
    parseLevelSetRepresentative(argc, argv, "psi", "fmm");
  const PhiInitialRepresentative phiInitialRepresentative =
    parsePhiInitialRepresentative(argc, argv, "phi-init");

  // -----  phi_h initial reconstruction (screened-Poisson)  ---------------
  // Source magnitude and screening length identical to
  // LevelSetLSRReconstruction (psi path): psiEll = ellMult * h,
  // M = magMult * psiEll. Calibration rescales |grad phi_h| -> |grad
  // phi_analytic| in the band.
  const Real kPsiEllMult =
    parseRealOption(argc, argv, "phi-ell-mult", Real(2));
  const Real kPsiEll = kPsiEllMult * h;
  const Real kPsiSourceMagMult =
    parseRealOption(argc, argv, "phi-source-mag-mult", Real(5));
  const Real kPsiSourceMagnitude = kPsiSourceMagMult * kPsiEll;

  // -----  Advection velocity  ---------------------------------------------
  const VelocityField velocityField =
    parseVelocityField(argc, argv, "velocity");
  const PhiRedistance phiRedistance =
    parsePhiRedistance(argc, argv, "phi-redistance");
  // dt per frame is one orbit step (T = 1, n frames => dt = 1/n).
  const Real kAdvectionDt = Real(1) / static_cast<Real>(nFrames);

  // ----- Convergence criterion ---------------------------------------------
  //   Residual : classic ||R|| < tol (fails to recognise the GN data-floor).
  //   Geometry : ||phi(X+u)||_RMS < ~h^2 (honest "did we capture Γ?").
  //   Either   : converged if EITHER residual or geometry tolerance holds.
  enum class ConvergenceMode { Residual, Geometry, Either };
  constexpr ConvergenceMode kConvergenceMode = ConvergenceMode::Either;
  // Knobs are CLI-overridable so convergence settings can be swept without
  // rebuilding.
#ifdef RODIN_LSR_P2_DISPLACEMENT
  constexpr Real kDefaultFitTolMult = Real(5.0);
#else
  constexpr Real kDefaultFitTolMult = Real(4.0);
#endif
  const Real        kFitTolMult      =
    parseRealOption(argc, argv, "fittol-mult", kDefaultFitTolMult); ///< fitTol = mult * h^2
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
#ifdef RODIN_LSR_P2_DISPLACEMENT
  using VectorFES = H1<2, Math::SpatialVector<Real>, LocalMesh>;
#else
  using VectorFES = P1<Math::SpatialVector<Real>, LocalMesh>;
#endif
  using ScalarP0 = P0<Real, LocalMesh>;

  ScalarP1 p1Fes(mesh);
  ScalarP0 p0Fes(mesh);
  using VectorP1Phi = P1<Math::SpatialVector<Real>, LocalMesh>;
  VectorP1Phi gradPhiFes(mesh, /*vdim=*/2);
  GridFunction gradPhiSmoothed(gradPhiFes);
  gradPhiSmoothed.getData().setZero();
#ifdef RODIN_LSR_P2_DISPLACEMENT
  VectorFES vectorFes(std::integral_constant<std::size_t, 2>{}, mesh, 2);
#else
  VectorFES vectorFes(mesh, 2);
#endif

  // Background mesh is invariant across the sweep, so the per-cell geometry
  // cache (sigma_K, det A_K, J_scale, gradN, area) only needs to be built
  // once. Recomputing it each frame was a measurable performance regression.
  auto cellGeomBg = precomputeCellGeometry(mesh);
  auto& cellCacheBg = cellGeomBg.first;

  GridFunction psi(p1Fes);            psi.setName("psi");
  // phiGf: analytic phi sampled at original X (diagnostic only).
  GridFunction phiGf(p1Fes);          phiGf.setName("phi_analytic");
  // phiH: ADVECTED P1 grid function used by the LSR push solve. Built
  // from screened Poisson at frame 0, advected one step per subsequent
  // frame. NO analytic dependency once initialised.
  GridFunction phiH(p1Fes);           phiH.setName("phi");
  // Diagnostic: pointwise difference phi_h - phi_analytic.
  GridFunction phiHError(p1Fes);      phiHError.setName("phi_error");
  GridFunction cellLabel(p0Fes);      cellLabel.setName("cell_label");
  GridFunction phaseMoment(p0Fes);    phaseMoment.setName("phase_moment");
  GridFunction sigmaKgf(p0Fes);       sigmaKgf.setName("sigma_K");

  GridFunction  u(vectorFes);         u.setName("displacement");

  // -------------------------------------------------------------------------
  // phi (moved) side. Same combinatorics as `mesh`; coordinates and, for
  // P2 displacement, cell polytope transformations change per frame. The
  // XDMF writer dumps the moved mesh anew at every snapshot
  // (MeshPolicy::Transient) so ParaView sees the curved fit per frame.
  // -------------------------------------------------------------------------
  LocalMesh moved(mesh);
  ScalarP1 p1FesMoved(moved);
  ScalarP0 p0FesMoved(moved);

  GridFunction cellLabelPhi(p0FesMoved); cellLabelPhi.setName("cell_label");
  // phi_redist : FMM signed distance computed on the LSR-displaced
  // (moved) mesh using the classifier interior/skeleton attributes.
  // This IS the LSR-reconstructed-domain SDF; it becomes phi_h for the
  // next frame's advection.
  GridFunction phiRedist(p1FesMoved);    phiRedist.setName("phi_redist");
  // Background copy of phiRedist for side-by-side visualisation on
  // psiGrid (background mesh) next to phi_analytic and phi_drifted.
  GridFunction phiRedistBg(p1Fes);       phiRedistBg.setName("phi_redist");
  GridFunction jKgf(p0FesMoved);        jKgf.setName("j");
  // q_abs = Q_abs(A_K^u): absolute intrinsic shape quality of the moved
  //                       cell (similarity of REFERENCE cell at value 1).
  // q_rel = Q_rel(F),  F = A_K^u (A_K)^{-1} = I + grad u:
  //                       relative-distortion measure minimised by the
  //                       LSR shape barrier; equals 1 at u = 0 for every
  //                       cell regardless of background shape.
  GridFunction qAbs(p0FesMoved);        qAbs.setName("q_abs");
  GridFunction qRel(p0FesMoved);        qRel.setName("q_rel");
  GridFunction phiMoved(p1FesMoved);    phiMoved.setName("phi_moved");

  // -------------------------------------------------------------------------
  // XDMF writer in transient mode (two grids: psi and phi).
  // -------------------------------------------------------------------------
  IO::XDMF xdmf("LevelSetLSRWavyCircleAdvectedSweep");
  auto psiGrid = xdmf.grid("psi");
  psiGrid.setMesh(mesh, IO::XDMF::MeshPolicy::Transient);
  psiGrid.add(cellLabel,   IO::XDMF::Center::Cell);
  psiGrid.add(phaseMoment, IO::XDMF::Center::Cell);
  psiGrid.add(sigmaKgf,    IO::XDMF::Center::Cell);
  psiGrid.add(phiGf,        IO::XDMF::Center::Node);  // phi_analytic
  psiGrid.add(phiH,         IO::XDMF::Center::Node);  // phi (drifted advected, pre-redistance, what LSR consumed)
  psiGrid.add(phiHError,    IO::XDMF::Center::Node);  // phi - phi_analytic
  psiGrid.add(phiRedistBg,  IO::XDMF::Center::Node);  // phi_redist (FMM of LSR-moved domain), copied to background DOFs
  psiGrid.add(psi,          IO::XDMF::Center::Node);
  psiGrid.add(u,            IO::XDMF::Center::Node);

  auto phiGrid = xdmf.grid("phi");
  phiGrid.setMesh(moved, IO::XDMF::MeshPolicy::Transient);
  phiGrid.add(cellLabelPhi, IO::XDMF::Center::Cell);
  phiGrid.add(jKgf,        IO::XDMF::Center::Cell);
  phiGrid.add(qAbs,        IO::XDMF::Center::Cell);
  phiGrid.add(qRel,        IO::XDMF::Center::Cell);
  phiGrid.add(phiMoved,    IO::XDMF::Center::Node);
  phiGrid.add(phiRedist,   IO::XDMF::Center::Node);

  // -------------------------------------------------------------------------
  // Frame loop.
  // -------------------------------------------------------------------------
  std::cout << "Wavy-circle LSR sweep on " << n << "x" << n
            << " unit-square mesh, " << nFrames << " frames\n";
  std::cout << "  R0=" << R0 << "  amp=" << amp << "  k=" << kLobes
            << "  orbit R=" << orbitR << '\n';
  std::cout << "  tangent="
            << (tangentMode == LSRTangent::GaussNewton ? "gn"
                : tangentMode == LSRTangent::Newton ? "newton" : "psd")
            << "  psi="
            << (psiRepresentative == LevelSetRepresentative::Fmm
                  ? "fmm"
                  : psiRepresentative == LevelSetRepresentative::ProjectedPhiH1
                      ? "projected-phi-h1"
                      : "screened")
            << "  phi-init="
            << (phiInitialRepresentative == PhiInitialRepresentative::Analytic
                  ? "analytic"
                  : phiInitialRepresentative == PhiInitialRepresentative::Fmm
                      ? "fmm"
                      : "screened")
            << "  phi-redistance="
            << (phiRedistance == PhiRedistance::FmmMoved
                  ? "fmm-moved"
                  : phiRedistance == PhiRedistance::ProjectedPhiH1Moved
                      ? "projected-phi-h1-moved"
                  : phiRedistance == PhiRedistance::ScreenedMoved
                      ? "screened-moved"
                      : "off")
            << '\n';

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

    // ---- Stage 0: semi-Lagrangian advection of phi_h (frame > 0) ----
    if (frame > 0)
    {
      AnalyticVectorFunction velocity(
          [&](const Geometry::Point& p) -> Math::SpatialVector<Real>
          {
            const auto& X = p.getCoordinates();
            switch (velocityField)
            {
              case VelocityField::Rotation:
              {
                const Real omega = Real(2) * Real(M_PI);
                return vec2(-omega * (X(1) - Real(0.5)),
                             omega * (X(0) - Real(0.5)));
              }
              case VelocityField::None:
              default:
                return vec2(Real(0), Real(0));
            }
          },
          /*dimension=*/2);
      TrialFunction advect(p1Fes);
      TestFunction  advectTest(p1Fes);
      Advection::Lagrangian adv(advect, advectTest, phiH, velocity);
      adv.step(kAdvectionDt);
      phiH.getData() = advect.getSolution().getData();
    }

    // ---- Stage 1: phase moments + classification ----
    //
    // Source of phi values for the moments:
    //   Frame 0 : analytic level set (phi_h not yet built).
    //   Frame k : the ADVECTED phi_h grid function -- we test the working
    //             pipeline on the field we actually carry forward.
    //
    // The interface skeleton produced here is consumed by both the
    // screened-Poisson Dirichlet (frame 0 only) and the FMM psi.
    const bool useAnalyticForMoments = (frame == 0);
    const auto cellMoments =
      useAnalyticForMoments
        ? collectCellMomentInfo(
              mesh,
              [&](const Vec2& p) { return levelSet.phi(p); },
              epsilon)
        : [&]()
          {
            // Evaluate phi_h at quadrature points by walking cells and
            // calling GridFunction::getValue at the reference coord.
            std::vector<CellMomentInfo> cells;
            cells.reserve(mesh.getCellCount());
            for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
            {
              const auto& cellPolytope = *cellIt;
              const auto& vertices = cellPolytope.getVertices();
              if (vertices.size() != 3)
                throw std::runtime_error(
                    "expects triangular cells.");
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
                Math::SpatialPoint rc(2);
                rc(0) = bary[1];
                rc(1) = bary[2];
                const Geometry::Point pt(cellPolytope, rc);
                const Real phiq = phiH.getValue(pt);
                moment += applyPhaseMomentMap(phiq, epsilon);
              }
              info.moment = moment / TriangleBarycentricQuadrature.size();
              cells.push_back(std::move(info));
            }
            return cells;
          }();

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

#if 0
    // (Disabled) Re-screened-Poisson phi_h block.
    //
    // Conceptually attractive but collapses the LSR objective when both
    // phi_h and psi_h are screened-Poisson lifts of the same classifier
    // skeleton: phi_h ~= psi_h up to a calibration scalar, residual ~ 0
    // at u = 0, Newton trivially "converges" without doing any
    // registration. Replaced below by a no-op for the input phi (use
    // analytic at frame 0, advected discrete phi_h at frame > 0) and an
    // FMM redistance from the classifier AFTER the LSR solve.
    {
      TrialFunction phiTrial(p1Fes);
      TestFunction  phiTest(p1Fes);
      RealFunction phiSource(
          [&](const Geometry::Point& p) -> Real
          {
            const auto attr = p.getPolytope().getAttribute();
            if (attr && *attr == interiorAttribute)
              return -kPsiSourceMagnitude;
            return kPsiSourceMagnitude;
          });
      Problem phiProblem(phiTrial, phiTest);
      phiProblem =
          Integral((kPsiEll * kPsiEll) * Grad(phiTrial), Grad(phiTest))
        + Integral(phiTrial, phiTest)
        - Integral(phiSource, phiTest)
        + DirichletBC(phiTrial, RealFunction(Real(0)))
            .on(interfaceAttribute);
      Solver::SparseLU(phiProblem).solve();
      phiH.getData() = phiTrial.getSolution().getData();

      // Band-averaged |grad phi_h_raw| using the raw field itself as
      // the band-weight argument (W = exp(-phi_h^2 / 2 deltaW^2)).
      Real gradAcc = 0, weightAcc = 0;
      for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
      {
        const auto& cell = *cellIt;
        const auto& vertices = cell.getVertices();
        if (vertices.size() != 3) continue;
        std::array<Vec2, 3> x;
        std::array<Real, 3> phiNode;
        for (std::size_t a = 0; a < 3; ++a)
        {
          x[a] = mesh.getVertexCoordinates(vertices[a]);
          Math::SpatialPoint rc(2);
          rc.setZero();
          if (a == 1) rc(0) = 1;
          if (a == 2) rc(1) = 1;
          phiNode[a] = phiH.getValue(Geometry::Point(cell, rc));
        }
        Math::SpatialMatrix<Real> M(2, 2);
        M(0, 0) = x[1](0) - x[0](0);
        M(0, 1) = x[2](0) - x[0](0);
        M(1, 0) = x[1](1) - x[0](1);
        M(1, 1) = x[2](1) - x[0](1);
        const Real detM = M(0, 0) * M(1, 1) - M(0, 1) * M(1, 0);
        Vec2 rhs = vec2(phiNode[1] - phiNode[0], phiNode[2] - phiNode[0]);
        Vec2 gradPhiCell(2);
        gradPhiCell(0) = ( M(1, 1) * rhs(0) - M(0, 1) * rhs(1)) / detM;
        gradPhiCell(1) = (-M(1, 0) * rhs(0) + M(0, 0) * rhs(1)) / detM;
        const Real gradNorm = gradPhiCell.norm();
        const Vec2 e1 = x[1] - x[0];
        const Vec2 e2 = x[2] - x[0];
        const Real triangleArea =
          Real(0.5) * std::abs(e1(0) * e2(1) - e1(1) * e2(0));
        for (const auto& bary : TriangleBarycentricQuadrature)
        {
          const Real phiq =
            bary[0] * phiNode[0] + bary[1] * phiNode[1]
          + bary[2] * phiNode[2];
          const Real wBand =
            std::exp(-phiq * phiq / (Real(2) * deltaW * deltaW));
          const Real qw =
            triangleArea / static_cast<Real>(TriangleBarycentricQuadrature.size());
          gradAcc   += wBand * gradNorm * qw;
          weightAcc += wBand * qw;
        }
      }
      const Real meanGrad = weightAcc > 0 ? gradAcc / weightAcc : Real(1);
      const Real phiScale = meanGrad > Real(1e-30) ? Real(1) / meanGrad : Real(1);
      phiH.getData() *= phiScale;
      std::cout << "    phi_h SP-rebuild: ell=" << kPsiEll
                << ", scale=" << phiScale
                << ", mean|grad_raw|=" << meanGrad << '\n';
    }

    // ---- Stage 2: screened-Poisson psi_h ----
    //
    // Same construction as phi_h above, calibrated to unit gradient.
    // psi enters the LSR ONLY through the band weight W(psi(X)) and the
    // push residual r = phi(X+u) - psi(X); both want psi and phi to
    // share scale, which the unit-gradient calibration guarantees.
    {
      TrialFunction psiTrial(p1Fes);
      TestFunction  psiTest(p1Fes);
      RealFunction psiSource(
          [&](const Geometry::Point& p) -> Real
          {
            const auto attr = p.getPolytope().getAttribute();
            if (attr && *attr == interiorAttribute)
              return -kPsiSourceMagnitude;
            return kPsiSourceMagnitude;
          });
      Problem psiProblem(psiTrial, psiTest);
      psiProblem =
          Integral((kPsiEll * kPsiEll) * Grad(psiTrial), Grad(psiTest))
        + Integral(psiTrial, psiTest)
        - Integral(psiSource, psiTest)
        + DirichletBC(psiTrial, RealFunction(Real(0)))
            .on(interfaceAttribute);
      Solver::SparseLU(psiProblem).solve();
      psi.getData() = psiTrial.getSolution().getData();

      Real gradAcc = 0, weightAcc = 0;
      for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
      {
        const auto& cell = *cellIt;
        const auto& vertices = cell.getVertices();
        if (vertices.size() != 3) continue;
        std::array<Vec2, 3> x;
        std::array<Real, 3> psiNode;
        for (std::size_t a = 0; a < 3; ++a)
        {
          x[a] = mesh.getVertexCoordinates(vertices[a]);
          Math::SpatialPoint rc(2);
          rc.setZero();
          if (a == 1) rc(0) = 1;
          if (a == 2) rc(1) = 1;
          psiNode[a] = psi.getValue(Geometry::Point(cell, rc));
        }
        Math::SpatialMatrix<Real> M(2, 2);
        M(0, 0) = x[1](0) - x[0](0);
        M(0, 1) = x[2](0) - x[0](0);
        M(1, 0) = x[1](1) - x[0](1);
        M(1, 1) = x[2](1) - x[0](1);
        const Real detM = M(0, 0) * M(1, 1) - M(0, 1) * M(1, 0);
        Vec2 rhs = vec2(psiNode[1] - psiNode[0], psiNode[2] - psiNode[0]);
        Vec2 gradCell(2);
        gradCell(0) = ( M(1, 1) * rhs(0) - M(0, 1) * rhs(1)) / detM;
        gradCell(1) = (-M(1, 0) * rhs(0) + M(0, 0) * rhs(1)) / detM;
        const Real gradNorm = gradCell.norm();
        const Vec2 e1 = x[1] - x[0];
        const Vec2 e2 = x[2] - x[0];
        const Real triangleArea =
          Real(0.5) * std::abs(e1(0) * e2(1) - e1(1) * e2(0));
        for (const auto& bary : TriangleBarycentricQuadrature)
        {
          const Real psiq =
            bary[0] * psiNode[0] + bary[1] * psiNode[1]
          + bary[2] * psiNode[2];
          const Real wBand =
            std::exp(-psiq * psiq / (Real(2) * deltaW * deltaW));
          const Real qw =
            triangleArea / static_cast<Real>(TriangleBarycentricQuadrature.size());
          gradAcc   += wBand * gradNorm * qw;
          weightAcc += wBand * qw;
        }
      }
      const Real meanGrad = weightAcc > 0 ? gradAcc / weightAcc : Real(1);
      const Real psiScale = meanGrad > Real(1e-30) ? Real(1) / meanGrad : Real(1);
      psi.getData() *= psiScale;
    }
#endif

    // ---- Stage 2: screened-Poisson psi_h on background mesh ----
    //
    //   (I - ell^2 Delta) psi_raw = +- M_src,  psi_raw = 0 on Gamma_psi,h.
    // Calibrated to band-averaged |grad psi| = 1. Smoother than FMM
    // (no kinks at the medial axis), which interacts better with the
    // semi-Lagrangian + L2-projection advection step that consumes the
    // post-frame redistanced phi_h on the next iteration.
    {
      TrialFunction psiTrial(p1Fes);
      TestFunction  psiTest(p1Fes);
      RealFunction psiSource(
          [&](const Geometry::Point& p) -> Real
          {
            const auto attr = p.getPolytope().getAttribute();
            if (attr && *attr == interiorAttribute)
              return -kPsiSourceMagnitude;
            return kPsiSourceMagnitude;
          });
      Problem psiProblem(psiTrial, psiTest);
      psiProblem =
          Integral((kPsiEll * kPsiEll) * Grad(psiTrial), Grad(psiTest))
        + Integral(psiTrial, psiTest)
        - Integral(psiSource, psiTest)
        + DirichletBC(psiTrial, RealFunction(Real(0)))
            .on(interfaceAttribute);
      Solver::SparseLU(psiProblem).solve();
      psi.getData() = psiTrial.getSolution().getData();

      Real gradAcc = 0, weightAcc = 0;
      for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
      {
        const auto& cell = *cellIt;
        const auto& vertices = cell.getVertices();
        if (vertices.size() != 3) continue;
        std::array<Vec2, 3> x;
        std::array<Real, 3> psiNode;
        for (std::size_t a = 0; a < 3; ++a)
        {
          x[a] = mesh.getVertexCoordinates(vertices[a]);
          Math::SpatialPoint rc(2);
          rc.setZero();
          if (a == 1) rc(0) = 1;
          if (a == 2) rc(1) = 1;
          psiNode[a] = psi.getValue(Geometry::Point(cell, rc));
        }
        Math::SpatialMatrix<Real> M(2, 2);
        M(0, 0) = x[1](0) - x[0](0);
        M(0, 1) = x[2](0) - x[0](0);
        M(1, 0) = x[1](1) - x[0](1);
        M(1, 1) = x[2](1) - x[0](1);
        const Real detM = M(0, 0) * M(1, 1) - M(0, 1) * M(1, 0);
        Vec2 rhs = vec2(psiNode[1] - psiNode[0], psiNode[2] - psiNode[0]);
        Vec2 gradCell(2);
        gradCell(0) = ( M(1, 1) * rhs(0) - M(0, 1) * rhs(1)) / detM;
        gradCell(1) = (-M(1, 0) * rhs(0) + M(0, 0) * rhs(1)) / detM;
        const Real gradNorm = gradCell.norm();
        const Vec2 e1 = x[1] - x[0];
        const Vec2 e2 = x[2] - x[0];
        const Real triangleArea =
          Real(0.5) * std::abs(e1(0) * e2(1) - e1(1) * e2(0));
        for (const auto& bary : TriangleBarycentricQuadrature)
        {
          const Real psiq =
            bary[0] * psiNode[0] + bary[1] * psiNode[1]
          + bary[2] * psiNode[2];
          const Real wBand =
            std::exp(-psiq * psiq / (Real(2) * deltaW * deltaW));
          const Real qw =
            triangleArea / static_cast<Real>(TriangleBarycentricQuadrature.size());
          gradAcc   += wBand * gradNorm * qw;
          weightAcc += wBand * qw;
        }
      }
      const Real meanGrad = weightAcc > 0 ? gradAcc / weightAcc : Real(1);
      const Real psiScale = meanGrad > Real(1e-30) ? Real(1) / meanGrad : Real(1);
      psi.getData() *= psiScale;
    }
    if (frame == 0
        && phiInitialRepresentative == PhiInitialRepresentative::Screened)
    {
      phiH.getData() = psi.getData();
    }
    if (psiRepresentative == LevelSetRepresentative::ProjectedPhiH1)
    {
      TrialFunction psiTrial(p1Fes);
      TestFunction  psiTest(p1Fes);
      RealFunction phiData(
          [&](const Geometry::Point& p) -> Real
          {
            if (frame == 0)
            {
              const auto& X = p.getCoordinates();
              return levelSet.phi(vec2(X(0), X(1)));
            }
            return phiH.getValue(p);
          });
      auto gradPhiHData = Grad(phiH);
      AnalyticVectorFunction gradPhiData(
          [&](const Geometry::Point& p) -> Math::SpatialVector<Real>
          {
            if (frame == 0)
            {
              const auto& X = p.getCoordinates();
              return levelSet.grad(vec2(X(0), X(1)));
            }
            return gradPhiHData.getValue(p);
          },
          /*dimension=*/2);

      Problem psiProblem(psiTrial, psiTest);
      RealFunction<Real> ell2(kPsiEll * kPsiEll);
      psiProblem =
          Integral(psiTrial, psiTest)
        + Integral(ell2 * Grad(psiTrial), Grad(psiTest))
        - Integral(phiData, psiTest)
        - Integral(ell2 * gradPhiData, Grad(psiTest))
        + DirichletBC(psiTrial, RealFunction(Real(0)))
            .on(interfaceAttribute);
      Solver::SparseLU(psiProblem).solve();
      psi.getData() = psiTrial.getSolution().getData();
    }
    if (psiRepresentative == LevelSetRepresentative::Fmm
        || (frame == 0
            && phiInitialRepresentative == PhiInitialRepresentative::Fmm))
    {
      GridFunction psiFmm(p1Fes);
      psiFmm = Real(0);
      Distance::Eikonal<ScalarP1, Math::Vector<Real>>(psiFmm)
        .setInterface(interfaceAttribute)
        .setInterior(interiorAttribute)
        .solve()
        .sign();
      if (psiRepresentative == LevelSetRepresentative::Fmm)
        psi.getData() = psiFmm.getData();
      if (frame == 0
          && phiInitialRepresentative == PhiInitialRepresentative::Fmm)
        phiH.getData() = psiFmm.getData();
    }

    // ---- Diagnostic: classifier interior centroid vs analytic centre ---
    {
      Real wx = 0, wy = 0, wA = 0;
      std::size_t insideCells = 0;
      for (std::size_t local = 0; local < classified.labels.size(); ++local)
      {
        if (classified.labels[local] != MinSTCut::Inside) continue;
        const auto& info = cellMoments[local];
        const Vec2 ctr =
          (info.x[0] + info.x[1] + info.x[2]) / Real(3);
        wx += info.area * ctr(0);
        wy += info.area * ctr(1);
        wA += info.area;
        ++insideCells;
      }
      const Real cxClass = wA > 0 ? wx / wA : Real(0);
      const Real cyClass = wA > 0 ? wy / wA : Real(0);
      const Real dxClass = cxClass - levelSet.cx;
      const Real dyClass = cyClass - levelSet.cy;
      const Real lag = std::sqrt(dxClass * dxClass + dyClass * dyClass);
      std::cout << "    classifier interior: cells=" << insideCells
                << "  area=" << std::fixed << std::setprecision(4) << wA
                << "  centroid=(" << cxClass << ", " << cyClass << ")"
                << "  vs c_analytic=(" << levelSet.cx << ", " << levelSet.cy
                << ")"
                << "  lag=" << std::scientific << std::setprecision(3)
                << lag << '\n';
    }

    // ---- Stage 3: smooth-band measure for the LSR normalisation ----
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
        const Real s = psi.getValue(pt);
        const Real W = std::exp(-s * s / (2 * deltaW * deltaW));
        weightedBandMeasure += W * w;
      }
    }
    if (weightedBandMeasure <= 0)
      throw std::runtime_error("Empty smooth band: M_w = 0 at frame "
                               + std::to_string(frame));

    // ---- Stage 4a: LSR solve (strategy-selected initial guess) ----

    auto& cellCache = cellCacheBg;

    if (frame == 0
        && phiInitialRepresentative == PhiInitialRepresentative::Analytic)
    {
      phiH = [&](const Geometry::Point& p) -> Real
      {
        const auto& X = p.getCoordinates();
        return levelSet.phi(vec2(X(0), X(1)));
      };
    }

    auto gradPhiDiscrete = Grad(phiH);
    {
      TrialFunction gtrial(gradPhiFes);
      TestFunction  gtest(gradPhiFes);
      Problem gpb(gtrial, gtest);
      gpb = Integral(gtrial, gtest)
          - Integral(gradPhiDiscrete, gtest);
      Solver::SparseLU(gpb).solve();
      gradPhiSmoothed.getData() = gtrial.getSolution().getData();
    }
    auto hessFromSmoothedGrad = Jacobian(gradPhiSmoothed);

    RealFunction phi(
        [&](const Geometry::Point& p) -> Real
        {
          return phiH.getValue(p);
        });
    AnalyticVectorFunction gradPhi(
        [&](const Geometry::Point& p) -> Math::SpatialVector<Real>
        {
          return gradPhiSmoothed.getValue(p);
        },
        /*dimension=*/2);
    AnalyticMatrixFunction hessPhi(
        [&](const Geometry::Point& p) -> Math::SpatialMatrix<Real>
        {
          return hessFromSmoothedGrad.getValue(p);
        },
        /*rows=*/2, /*cols=*/2);

    const bool kVerbose = hasFlag(argc, argv, "verbose");

    LSRParameters baseParams;
    baseParams.rhoS = 1;
    baseParams.deltaW = deltaW;
    baseParams.hRef = h;
    baseParams.normalizer = Real(1) / (weightedBandMeasure * h * h);
    baseParams.quadratureOrder = kLSRQuadratureOrder;
    baseParams.h1RegularizationWeight = kH1RegularizationWeight;
    baseParams.shapeWeight = kShapeWeight;
    // BEST-QUALITY profile safety nets.
    baseParams.jBarrierWeight = parseRealOption(
        argc, argv, "j-barrier-weight", Real(1.0));
    baseParams.jBarrierSafeRatio = parseRealOption(
        argc, argv, "j-barrier-safe", Real(0.5));
    baseParams.jVolumeTetherWeight = parseRealOption(
        argc, argv, "volume-tether-weight", Real(0.01));
    baseParams.jMinRatio = Real(1e-8);
    baseParams.jSafeRatio = Real(1e-3);
    baseParams.lineSearchSafetyMargin = kLineSearchSafetyMargin;
    baseParams.qRelMax = kQRelMax;
    baseParams.alphaInit = kLineSearchAlphaInit;
    baseParams.alphaReduction = kLineSearchReduction;
    baseParams.alphaMin = kLineSearchAlphaMin;
    baseParams.absoluteTolerance = kResidualAbsTol;
    baseParams.relativeTolerance = kResidualRelTol;
    baseParams.maxNewtonIterations = 80;
    // Let the sweep-level geometry criterion decide when to stop; residual
    // stall alone is not a failure for this registration objective.
    baseParams.stallPatience = kStallPatience;
    baseParams.initialGuessMetric = initialGuessMetric;
    baseParams.initialGuessScaling = initialGuessScaling;
    baseParams.initialGuessElasticityLambda = initialGuessElasticityLambda;
    baseParams.initialGuessElasticityMu = initialGuessElasticityMu;
    baseParams.tangent = tangentMode;

    // Evaluate phi_h at a moved physical point (a + ua) by walking the
    // mesh and using its inverse polytope map. Simple O(N_cells) fallback
    // — fine here because we only sample ~|skeleton| points.
    auto phiHAtMovedVertex = [&](const Vec2& y) -> Real
    {
      for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
      {
        const auto& cell = *cellIt;
        Math::SpatialPoint rc(2);
        const auto& trf = cell.getTransformation();
        Math::SpatialPoint Y(2);
        Y(0) = y(0); Y(1) = y(1);
        trf.inverse(rc, Y);
        // Inside-triangle test in reference (xi, eta in [0,1], xi+eta<=1).
        const Real xi  = rc(0), eta = rc(1);
        const Real tol = Real(1e-9);
        if (xi >= -tol && eta >= -tol && xi + eta <= Real(1) + tol)
          return phiH.getValue(Geometry::Point(cell, rc, Y));
      }
      return Real(0);  // outside the mesh — should not happen for moved vertices
    };

    auto computeInterfaceFit = [&](bool discrete) -> Real
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
        Real phiA, phiB;
        if (discrete)
        {
          phiA = phiHAtMovedVertex(a + ua);
          phiB = phiHAtMovedVertex(b + ub);
        }
        else
        {
          phiA = levelSet.phi(a + ua);
          phiB = levelSet.phi(b + ub);
        }
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

    LSRParameters params = baseParams;
    switch (kInitialGuessStrategy)
    {
      case InitialGuessStrategy::Cold:
        params.initialGuess = LSRInitialGuess::Zero;
        break;
      case InitialGuessStrategy::Warm:
        params.initialGuess =
          frame == 0 ? LSRInitialGuess::Zero : LSRInitialGuess::Current;
        break;
      case InitialGuessStrategy::Hilbert:
        params.initialGuess = LSRInitialGuess::Hilbert;
        break;
    }
    // Convergence is judged on the DISCRETE fit (what the solver is
    // actually optimising against). The analytic fit is reported
    // alongside as ground-truth quality.
    params.acceptedStateConvergenceTest =
      [&](const LSRReport&) -> bool
      {
        interfaceFit = computeInterfaceFit(/*discrete=*/true);
        return interfaceFit < kInterfaceFitTol;
      };

    LSR lsr(u);
    const LSRReport lsrReport =
      lsr.setParameters(params).solve(psi, phi, gradPhi, hessPhi);

    r0 = lsrReport.initialResidual;
    itCount = lsrReport.iterations;
    residualBest = lsrReport.finalResidual;
    if (interfaceFit == Real(0))
      interfaceFit = computeInterfaceFit(/*discrete=*/true);
    const Real interfaceFitAnalytic = computeInterfaceFit(/*discrete=*/false);

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
      std::cout << "      LSR:"
                << " it=" << lsrReport.iterations
                << "  ||R||=" << std::scientific << std::setprecision(3)
                << lsrReport.finalResidual
                << "  alpha=" << lsrReport.lastAcceptedAlpha
                << "  backtracks=" << lsrReport.totalBacktracks
                << "  j_ls=" << lsrReport.jLineSearchRatio
                << "  min_j=" << lsrReport.minJRatio
                << (lsrReport.lineSearchFailed ? "  [line-search failed]" : "")
                << '\n';
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

    // ---- Stage 5a: psi-side cell-P0 fields + nodal phi ----
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
    // phiHError = phi_h - phi_analytic at vertices (diagnostic of
    // accumulated advection error).
    phiHError = [&](const Geometry::Point& p) -> Real
    {
      const auto& X = p.getCoordinates();
      const Real ana = levelSet.phi(vec2(X(0), X(1)));
      const Real dis = phiH.getValue(p);
      return dis - ana;
    };

    // ---- Stage 5b: build the moved (phi) mesh and its phi fields ----
    //
    // The phi mesh shares combinatorics with `mesh`; its geometry is pushed
    // by u. For P2 displacement this includes P2 cell transformations.
    // Region attributes (interior /
    // exterior / interface / boundary) carry over so the phi grid carries
    // the same classification ParaView sees on psi.
    {
      updateMovedMeshFromDisplacement(mesh, moved, u);
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

    // Per-cell centroid j_K^u and absolute Q diagnostics on the moved mesh.
    // The sigma_K branch is fixed on the original mesh; for curved P2 cells
    // the affine cache reports a centroid summary, while the line search uses
    // sampled admissibility on the actual displacement field.
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
        const Real jKval = sigDetAu / src.Jscale;
        jKgf.getData()(dof) = jKval;
        const Real dExp = Real(2) / Real(d);
        qAbs.getData()(dof) =
          dst.A.squaredNorm()
            / (static_cast<Real>(d) * std::pow(sigDetAu, dExp));
        const Math::SpatialMatrix<Real> F = dst.A * src.A.inverse();
        qRel.getData()(dof) =
          F.squaredNorm()
            / (static_cast<Real>(d) * std::pow(jKval, dExp));
        cellLabelPhi.getData()(dof) =
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
              << "  fit_h=" << std::setprecision(3) << interfaceFit
              << "  fit_an=" << std::setprecision(3) << interfaceFitAnalytic
              << "  init=" << strategyName
              << "  converged=" << convergedReason << '\n';
    if (report.converged)
      ++framesConverged;

#if 1
    // ---- Stage 5c: screened-Poisson redistance of the LSR-reconstructed
    //                domain (on the MOVED mesh).
    //
    // The classifier interior labels live on cells; their physical
    // positions are on the MOVED mesh (X + u(X)). Solving
    //   (I - ell^2 Delta) phi_redist = +- M_src,
    //   phi_redist = 0 on Gamma_psi,h
    // on the moved mesh therefore lifts the LSR-fit interface to a
    // smooth band-shaped scalar, calibrated to unit gradient. Compared
    // with FMM-redistancing the same domain, the screened-Poisson lift
    // has no medial-axis kinks and a wider, smoother band -- both
    // matter because phi_h is consumed by the L2-projection
    // semi-Lagrangian step at the next iteration, which dissipates
    // band-localised high-frequency content. The smoother input means
    // less per-step relative dissipation.
    if (phiRedistance == PhiRedistance::ScreenedMoved)
    {
      TrialFunction phiTrial(p1FesMoved);
      TestFunction  phiTest(p1FesMoved);
      RealFunction phiSource(
          [&](const Geometry::Point& p) -> Real
          {
            const auto attr = p.getPolytope().getAttribute();
            if (attr && *attr == interiorAttribute)
              return -kPsiSourceMagnitude;
            return kPsiSourceMagnitude;
          });
      Problem phiProblem(phiTrial, phiTest);
      phiProblem =
          Integral((kPsiEll * kPsiEll) * Grad(phiTrial), Grad(phiTest))
        + Integral(phiTrial, phiTest)
        - Integral(phiSource, phiTest)
        + DirichletBC(phiTrial, RealFunction(Real(0)))
            .on(interfaceAttribute);
      Solver::SparseLU(phiProblem).solve();
      phiRedist.getData() = phiTrial.getSolution().getData();

      Real gradAcc = 0, weightAcc = 0;
      for (auto cellIt = moved.getCell(); cellIt; ++cellIt)
      {
        const auto& cell = *cellIt;
        const auto& vertices = cell.getVertices();
        if (vertices.size() != 3) continue;
        std::array<Vec2, 3> x;
        std::array<Real, 3> phiNode;
        for (std::size_t a = 0; a < 3; ++a)
        {
          x[a] = moved.getVertexCoordinates(vertices[a]);
          Math::SpatialPoint rc(2);
          rc.setZero();
          if (a == 1) rc(0) = 1;
          if (a == 2) rc(1) = 1;
          phiNode[a] = phiRedist.getValue(Geometry::Point(cell, rc));
        }
        Math::SpatialMatrix<Real> M(2, 2);
        M(0, 0) = x[1](0) - x[0](0);
        M(0, 1) = x[2](0) - x[0](0);
        M(1, 0) = x[1](1) - x[0](1);
        M(1, 1) = x[2](1) - x[0](1);
        const Real detM = M(0, 0) * M(1, 1) - M(0, 1) * M(1, 0);
        Vec2 rhs = vec2(phiNode[1] - phiNode[0], phiNode[2] - phiNode[0]);
        Vec2 gradCell(2);
        gradCell(0) = ( M(1, 1) * rhs(0) - M(0, 1) * rhs(1)) / detM;
        gradCell(1) = (-M(1, 0) * rhs(0) + M(0, 0) * rhs(1)) / detM;
        const Real gradNorm = gradCell.norm();
        const Vec2 e1 = x[1] - x[0];
        const Vec2 e2 = x[2] - x[0];
        const Real triangleArea =
          Real(0.5) * std::abs(e1(0) * e2(1) - e1(1) * e2(0));
        for (const auto& bary : TriangleBarycentricQuadrature)
        {
          const Real phiq =
            bary[0] * phiNode[0] + bary[1] * phiNode[1]
          + bary[2] * phiNode[2];
          const Real wBand =
            std::exp(-phiq * phiq / (Real(2) * deltaW * deltaW));
          const Real qw =
            triangleArea / static_cast<Real>(TriangleBarycentricQuadrature.size());
          gradAcc   += wBand * gradNorm * qw;
          weightAcc += wBand * qw;
        }
      }
      const Real meanGrad = weightAcc > 0 ? gradAcc / weightAcc : Real(1);
      const Real phiScale = meanGrad > Real(1e-30) ? Real(1) / meanGrad : Real(1);
      phiRedist.getData() *= phiScale;
    }
    // Background-grid copy for psiGrid visualisation; same DOF layout.
    if (phiRedistance == PhiRedistance::ScreenedMoved)
    {
      phiRedistBg.getData() = phiRedist.getData();
      phiH.getData() = phiRedist.getData();
    }
#endif
    if (phiRedistance == PhiRedistance::ProjectedPhiH1Moved)
    {
      auto gradPhiHCurrent = Grad(phiH);
      auto localizedOnBackground =
        [&](const Geometry::Point& p, const auto& callback)
        {
          const auto& Yphys = p.getPhysicalCoordinates();
          for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
          {
            const auto& cell = *cellIt;
            Math::SpatialPoint rc(2);
            cell.getTransformation().inverse(rc, Yphys);
            const Real xi = rc(0), eta = rc(1);
            const Real tol = Real(1e-9);
            if (xi >= -tol && eta >= -tol
                && xi + eta <= Real(1) + tol)
              return callback(Geometry::Point(cell, rc, Yphys));
          }
          return callback(p);
        };
      RealFunction phiDataMoved(
          [&](const Geometry::Point& p) -> Real
          {
            return localizedOnBackground(
                p,
                [&](const Geometry::Point& q) -> Real
                {
                  return phiH.getValue(q);
                });
          });
      AnalyticVectorFunction gradPhiDataMoved(
          [&](const Geometry::Point& p) -> Math::SpatialVector<Real>
          {
            return localizedOnBackground(
                p,
                [&](const Geometry::Point& q) -> Math::SpatialVector<Real>
                {
                  return gradPhiHCurrent.getValue(q);
                });
          },
          /*dimension=*/2);

      TrialFunction phiTrial(p1FesMoved);
      TestFunction  phiTest(p1FesMoved);
      RealFunction<Real> ell2(kPsiEll * kPsiEll);
      Problem phiProblem(phiTrial, phiTest);
      phiProblem =
          Integral(phiTrial, phiTest)
        + Integral(ell2 * Grad(phiTrial), Grad(phiTest))
        - Integral(phiDataMoved, phiTest)
        - Integral(ell2 * gradPhiDataMoved, Grad(phiTest))
        + DirichletBC(phiTrial, RealFunction(Real(0)))
            .on(interfaceAttribute);
      Solver::SparseLU(phiProblem).solve();
      phiRedist.getData() = phiTrial.getSolution().getData();
      phiRedistBg.getData() = phiRedist.getData();
      phiH.getData() = phiRedist.getData();
    }
    if (phiRedistance == PhiRedistance::FmmMoved)
    {
      phiRedist = Real(0);
      Distance::Eikonal<ScalarP1, Math::Vector<Real>>(phiRedist)
        .setInterface(interfaceAttribute)
        .setInterior(interiorAttribute)
        .solve()
        .sign();
      phiRedistBg.getData() = phiRedist.getData();
      phiH.getData() = phiRedist.getData();
    }
    else if (phiRedistance == PhiRedistance::Off)
    {
      phiRedist = Real(0);
      phiRedistBg = Real(0);
    }

    // ---- Stage 6: append XDMF snapshot ----
    // XDMF write captures `phi` (the high-fidelity advected field). No
    // post-write reassignment: phi_h stays as the advected result so
    // the next frame's advection extends the high-fidelity transport
    // history rather than restarting from a reconstructed surrogate.
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
