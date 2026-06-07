/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//
// Two-wavy-circles merging test for the LSR pipeline.
//
// Two wavy circles travel on opposite (CCW and CW) circular orbits about
// the centre of the unit square. They meet head-on twice per cycle, at
// which point the implicit shape
//
//     phi(p) = min(phi_1(p), phi_2(p))
//
// changes topology: two disjoint blobs merge into one and split again. The
// classify-then-displace pipeline of LevelSetLSRWavyCircleSweep is rerun
// from scratch every frame, so the classifier sees a genuine topology
// change on the psi skeleton and the LSR is challenged with merging /
// splitting interfaces — a setting where warm-start is unsafe and the
// Hilbert lift is the principled initial guess.
//
// As in the parent sweep, the example exercises:
//   - Repeated classification + Newton solves at moving interface
//     locations, including frames straddling a topology change.
//   - The shard-aware XDMF writers under non-trivial cell-attribute
//     turnover.
//   - PSD-projected full-Newton tangent on a non-smooth phi (min of
//     two smooth phi's: continuous, piecewise C^2 with a ridge where
//     phi_1 = phi_2). The gradient and Hessian are taken from the
//     active branch (whichever phi_i is smaller); the ridge has
//     measure zero and is harmless at quadrature points.
//
// FES-independence and Newton-tangent choices follow the same honest
// caveats as the parent example: LSR uses sampled admissibility and a
// sampled relative-distortion energy, with PSD-projected Newton tangents.
//
#include <Rodin/Adaptation.h>
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
  // Single wavy-circle level set (same as the parent sweep example).
  // -------------------------------------------------------------------------
  struct WavyCircleLevelSet
  {
    Real cx = Real(0.5);
    Real cy = Real(0.5);
    Real R0 = Real(0.14);
    Real amp = Real(0.03);
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

      Math::SpatialMatrix<Real> Hr(2, 2);
      Hr(0, 0) = (Real(1) - dx * dx / r2) / r;
      Hr(1, 1) = (Real(1) - dy * dy / r2) / r;
      Hr(0, 1) = -dx * dy / r3;
      Hr(1, 0) = Hr(0, 1);

      Math::SpatialMatrix<Real> Ht(2, 2);
      Ht(0, 0) =  Real(2) * dx * dy / r4;
      Ht(1, 1) = -Real(2) * dx * dy / r4;
      Ht(0, 1) = (dy * dy - dx * dx) / r4;
      Ht(1, 0) = Ht(0, 1);

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
  // Composite level set: log-sum-exp soft-min of phi_a and phi_b.
  //
  //     phi_eps(p) = -eps * log( exp(-phi_a/eps) + exp(-phi_b/eps) ).
  //
  // As eps -> 0, phi_eps -> min(phi_a, phi_b) uniformly with error at
  // most eps*log(2). For eps > 0, phi_eps is C^infinity, and its
  // gradient / Hessian are smooth convex combinations of the individual
  // gradients / Hessians plus a (PSD) outer-product term in eps.
  //
  // The hard min was tried first and produced 0/50 convergence on the
  // n=50 sweep: cells straddling the ridge {phi_a = phi_b} saw
  // inconsistent gradients across quadrature points, and the
  // PSD-projected Newton tangent (which assumes phi is locally C^2 in
  // the LSR band) had no way to converge. Soft-min with eps = O(h)
  // removes the ridge singularity at the cost of an O(eps) bias in the
  // interface location, which is absorbed by the same O(h^2) FMM
  // discretisation error as the parent sweep.
  //
  // Numerically robust evaluation uses the standard "subtract the min"
  // trick to avoid overflow in exp() when |phi/eps| is large.
  // -------------------------------------------------------------------------
  struct TwoWavyCirclesLevelSet
  {
    WavyCircleLevelSet a;
    WavyCircleLevelSet b;
    Real eps = Real(1e-2); ///< softness; should be set ~ O(h) by caller.

    Real phi(const Vec2& p) const
    {
      const Real pa = a.phi(p);
      const Real pb = b.phi(p);
      const Real m = std::min(pa, pb);
      const Real ea = std::exp(-(pa - m) / eps);
      const Real eb = std::exp(-(pb - m) / eps);
      return m - eps * std::log(ea + eb);
    }

    // Branch weights w_a, w_b with w_a + w_b = 1.
    void weights(const Vec2& p, Real& wa, Real& wb) const
    {
      const Real pa = a.phi(p);
      const Real pb = b.phi(p);
      const Real m = std::min(pa, pb);
      const Real ea = std::exp(-(pa - m) / eps);
      const Real eb = std::exp(-(pb - m) / eps);
      const Real Z = ea + eb;
      wa = ea / Z;
      wb = eb / Z;
    }

    Vec2 grad(const Vec2& p) const
    {
      Real wa, wb;
      weights(p, wa, wb);
      return wa * a.grad(p) + wb * b.grad(p);
    }

    // Hessian of phi_eps:
    //   H = w_a H_a + w_b H_b
    //       - (1/eps) * [ w_a g_a g_a^T + w_b g_b g_b^T - gbar gbar^T ],
    //
    // with gbar = w_a g_a + w_b g_b. The bracketed term is the weighted
    // covariance of the branch gradients (symmetric PSD). The MINUS
    // sign is correct: at a convex corner of the union (where the two
    // branch normals disagree) the smoothed SDF has NEGATIVE curvature
    // along the ridge, reflecting the local concavity of the boundary
    // of the union. Sanity check: phi_a(x) = x, phi_b(x) = -x gives
    // phi_eps(x) = -eps log(2 cosh(x/eps)), phi_eps''(0) = -1/eps,
    // which agrees with the formula here.
    //
    // The PSD-projected Newton tangent will discard the negative
    // eigenvalues of r*Hess(phi_eps) in the LSR block, so the sign of
    // the Hessian here is allowed to be indefinite.
    Math::SpatialMatrix<Real> hess(const Vec2& p) const
    {
      Real wa, wb;
      weights(p, wa, wb);
      const Vec2 ga = a.grad(p);
      const Vec2 gb = b.grad(p);
      const Vec2 gbar = wa * ga + wb * gb;

      Math::SpatialMatrix<Real> H = wa * a.hess(p) + wb * b.hess(p);
      auto addOuter = [&](Real coef, const Vec2& u, const Vec2& v)
      {
        H(0, 0) += coef * u(0) * v(0);
        H(0, 1) += coef * u(0) * v(1);
        H(1, 0) += coef * u(1) * v(0);
        H(1, 1) += coef * u(1) * v(1);
      };
      const Real invEps = Real(1) / eps;
      addOuter(-wa * invEps, ga, ga);
      addOuter(-wb * invEps, gb, gb);
      addOuter( invEps,      gbar, gbar);
      return H;
    }
  };

  // -------------------------------------------------------------------------
  // Triangle-barycentric quadrature for the phase moments.
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
            "LevelSetLSRTwoWavyCirclesMerge expects triangular cells.");

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
      if (std::string(argv[i]) == flag) return true;
    return false;
  }

  LSRHilbertMetric parseHilbertMetric(
      int argc, char** argv, const std::string& name)
  {
    const std::string prefix = "--" + name + "=";
    for (int i = 1; i < argc; ++i)
    {
      const std::string arg(argv[i]);
      if (arg.rfind(prefix, 0) != 0) continue;
      const std::string v = arg.substr(prefix.size());
      if (v == "harmonic")   return LSRHilbertMetric::Harmonic;
      if (v == "elasticity") return LSRHilbertMetric::Elasticity;
      if (v == "shape")      return LSRHilbertMetric::ShapeHessian;
    }
    return LSRHilbertMetric::Harmonic;
  }

  LSRInitialGuessScaling parseInitialGuessScaling(
      int argc, char** argv, const std::string& name)
  {
    const std::string prefix = "--" + name + "=";
    for (int i = 1; i < argc; ++i)
    {
      const std::string arg(argv[i]);
      if (arg.rfind(prefix, 0) != 0) continue;
      const std::string v = arg.substr(prefix.size());
      if (v == "raw")        return LSRInitialGuessScaling::Unnormalized;
      if (v == "energy")     return LSRInitialGuessScaling::EnergyNormalized;
      if (v == "band")       return LSRInitialGuessScaling::BandNormalized;
    }
    return LSRInitialGuessScaling::Unnormalized;
  }

  LSRTangent parseTangentMode(
      int argc, char** argv, const std::string& name)
  {
    const std::string prefix = "--" + name + "=";
    for (int i = 1; i < argc; ++i)
    {
      const std::string arg(argv[i]);
      if (arg.rfind(prefix, 0) != 0) continue;
      const std::string v = arg.substr(prefix.size());
      if (v == "psd")    return LSRTangent::PSDProjectedNewton;
      if (v == "newton") return LSRTangent::Newton;
      if (v == "gn")     return LSRTangent::GaussNewton;
      throw std::runtime_error(
          "Unknown --" + name + "=" + v
          + " (expected psd, newton, or gn).");
    }
    return LSRTangent::PSDProjectedNewton;
  }

  // Cold    : u0 = 0 every frame.
  // Warm    : u0 = previous frame's converged u (zero on frame 0).
  // Hilbert : u0 = Riesz lift of -delta E_LSR(0) in H^1_0.
  enum class InitialGuessStrategy { Cold, Warm, Hilbert };

  InitialGuessStrategy parseInitialGuessStrategy(
      int argc, char** argv, const std::string& name)
  {
    const std::string prefix = "--" + name + "=";
    for (int i = 1; i < argc; ++i)
    {
      const std::string arg(argv[i]);
      if (arg.rfind(prefix, 0) != 0) continue;
      const std::string v = arg.substr(prefix.size());
      if (v == "cold")    return InitialGuessStrategy::Cold;
      if (v == "warm")    return InitialGuessStrategy::Warm;
      if (v == "hilbert") return InitialGuessStrategy::Hilbert;
      throw std::runtime_error(
          "Unknown --" + name + "=" + v
          + " (expected cold, warm, or hilbert).");
    }
    return InitialGuessStrategy::Hilbert;
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
  // Orbit and shape match LevelSetLSRWavyCircleSweep exactly: orbitR
  // = 0.10 about (0.5, 0.5), nominal radius R0 = 0.20, amplitude 0.05,
  // kLobes = 6. The two wavy circles follow the SAME orbit, only in
  // opposite directions (CCW for `a`, CW for `b`). Centre-to-centre
  // distance varies as 2*orbitR*|sin(angle)| from 0 (head-on at
  // angle=0, pi) to 2*orbitR (maximum separation at pi/2, 3pi/2).
  // Since 2*orbitR = 0.20 < 2*R0 = 0.40, the two regions OVERLAP for
  // the entire cycle, but the union geometry varies continuously from
  // a single near-circular blob (head-on) to a peanut-shaped union
  // (lateral separation). The implicit shape phi = min(phi_a, phi_b)
  // has a ridge along the locus phi_a = phi_b which sweeps across the
  // mesh as the centres orbit; the classifier and LSR see a smoothly
  // deforming connected region with non-trivial concave geometry near
  // the merge ridge.
  // Default geometry matches LevelSetLSRWavyCircleSweep: features are
  // well-resolved (R0 - amp = 0.15 = ~7h at n=50). With 2*orbitR = 0.20
  // < 2*(R0 - amp) = 0.30 the two regions OVERLAP throughout the cycle,
  // so this is a "merge ridge sweeps across the mesh" stress test, not
  // a topology-change one. For a genuine topology-change test, pass
  //   --orbitR=0.20 --R0=0.13 --amp=0.03 --n=100
  // (smaller, well-resolved circles on a wider orbit). At n=50 with the
  // default R0=0.20 you get the always-overlap union variant.
  const Real        orbitR =
    parseRealOption(argc, argv, "orbitR", Real(0.10));
  const Real        amp =
    parseRealOption(argc, argv, "amp", Real(0.05));
  const Real        R0 =
    parseRealOption(argc, argv, "R0", Real(0.20));
  constexpr Real    kLobes = 6;

  const Real h = Real(1) / static_cast<Real>(n - 1);
  const Real epsilon = 1.25 * h;
  const Real lambdaC = 0.008;
  const Real delta   = 1.75 * h;
  const Real deltaW  = 1.5 * delta;

  // Soft-min smoothing for phi = LSE(phi_a, phi_b). Default eps = 0.5*h
  // keeps the smoothing zone narrower than one cell, so the interface
  // bias is below the FMM consistency floor of Lemma~fit-floor.
  const Real kSoftMinEpsMult =
    parseRealOption(argc, argv, "softmin-eps-mult", Real(0.5));
  const Real kSoftMinEps = kSoftMinEpsMult * h;

  // Dimensionless γ_shape default = sqrt(h_ref).
  // The effective weight is γ_shape · normalizer · h_ref² inside LSR.h.
  // See LSR.h docstring for the dimensionless-energy formulation and
  // the (P_k, n, lobes) sweep that picked sqrt(h_ref) as the robust
  // default: at h ≈ 0.02 (n=50), γ_shape ≈ 0.14; at h ≈ 0.01 (n=100),
  // γ_shape ≈ 0.10. Same value works for P1 and P2.
  const Real kShapeWeight =
    parseRealOption(argc, argv, "gamma", std::sqrt(h));

  const Real   kLineSearchAlphaInit = Real(1);
  const Real   kLineSearchReduction = Real(0.5);
  const Real   kLineSearchAlphaMin  = Real(1e-6);
  const Real   kLineSearchSafetyMargin = Real(10);
  // Relative-distortion cap (identity-neutral): max Q_rel(I + grad u_h)
  // accepted by the line search. Default 10 reproduces the historical
  // effective behaviour.
  const Real   kQRelMax =
    parseRealOption(argc, argv, "qrel-max", Real(10));
  const std::size_t kLSRQuadratureOrder =
    parseSizeTOption(argc, argv, "lsr-quad-order", 0);
  const Real kH1RegularizationWeight =
    parseRealOption(argc, argv, "h1-reg", Real(0));

  // Smooth relative-Q barrier (interior-point penalty steering Newton
  // away from high relative distortion). With qbar-weight > 0 and a
  // finite qbar-max > qbar-act, the assembled tangent is bent so the
  // natural step length alpha=1 already respects the envelope. Default
  // disabled (qbar-weight=0).
  const Real kQBarrierWeight =
    parseRealOption(argc, argv, "qbar-weight", Real(0));
  const Real kQBarrierAct =
    parseRealOption(argc, argv, "qbar-act",
                    std::numeric_limits<Real>::infinity());
  const Real kQBarrierMax =
    parseRealOption(argc, argv, "qbar-max",
                    std::numeric_limits<Real>::infinity());

  // For a merging-topology sweep, warm-start across a topology change
  // is incorrect by construction: the previous frame's u lives on a
  // classifier with a different skeleton. Default to Hilbert and let
  // the cold-fallback rescue any frame that straddles a merge.
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

  enum class ConvergenceMode { Residual, Geometry, Either };
  constexpr ConvergenceMode kConvergenceMode = ConvergenceMode::Either;
#ifdef RODIN_LSR_P2_DISPLACEMENT
  constexpr Real kDefaultFitTolMult = Real(5.0);
#else
  constexpr Real kDefaultFitTolMult = Real(4.0);
#endif
  const Real        kFitTolMult      =
    parseRealOption(argc, argv, "fittol-mult", kDefaultFitTolMult);
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
#ifdef RODIN_LSR_P2_DISPLACEMENT
  using VectorFES = H1<2, Math::SpatialVector<Real>, LocalMesh>;
#else
  using VectorFES = P1<Math::SpatialVector<Real>, LocalMesh>;
#endif
  using ScalarP0 = P0<Real, LocalMesh>;

  ScalarP1 p1Fes(mesh);
  ScalarP0 p0Fes(mesh);
#ifdef RODIN_LSR_P2_DISPLACEMENT
  VectorFES vectorFes(std::integral_constant<std::size_t, 2>{}, mesh, 2);
#else
  VectorFES vectorFes(mesh, 2);
#endif

  auto cellGeomBg = precomputeCellGeometry(mesh);
  auto& cellCacheBg = cellGeomBg.first;

  GridFunction psi(p1Fes);            psi.setName("psi");
  GridFunction phiGf(p1Fes);          phiGf.setName("phi");
  GridFunction cellLabel(p0Fes);      cellLabel.setName("cell_label");
  GridFunction phaseMoment(p0Fes);    phaseMoment.setName("phase_moment");
  GridFunction sigmaKgf(p0Fes);       sigmaKgf.setName("sigma_K");

  GridFunction  u(vectorFes);         u.setName("displacement");

  LocalMesh moved(mesh);
  ScalarP1 p1FesMoved(moved);
  ScalarP0 p0FesMoved(moved);

  GridFunction cellLabelPhi(p0FesMoved); cellLabelPhi.setName("cell_label");
  GridFunction jKgf(p0FesMoved);        jKgf.setName("j");
  // q_abs = Q_abs(A_K^u) ; q_rel = Q_rel(I + grad u). See header.
  GridFunction qAbs(p0FesMoved);        qAbs.setName("q_abs");
  GridFunction qRel(p0FesMoved);        qRel.setName("q_rel");
  GridFunction phiMoved(p1FesMoved);    phiMoved.setName("phi_moved");

  IO::XDMF xdmf("LevelSetLSRTwoWavyCirclesMerge");
  auto psiGrid = xdmf.grid("psi");
  psiGrid.setMesh(mesh, IO::XDMF::MeshPolicy::Transient);
  psiGrid.add(cellLabel,   IO::XDMF::Center::Cell);
  psiGrid.add(phaseMoment, IO::XDMF::Center::Cell);
  psiGrid.add(sigmaKgf,    IO::XDMF::Center::Cell);
  psiGrid.add(phiGf,       IO::XDMF::Center::Node);
  psiGrid.add(psi,         IO::XDMF::Center::Node);
  psiGrid.add(u,           IO::XDMF::Center::Node);

  auto phiGrid = xdmf.grid("phi");
  phiGrid.setMesh(moved, IO::XDMF::MeshPolicy::Transient);
  phiGrid.add(cellLabelPhi, IO::XDMF::Center::Cell);
  phiGrid.add(jKgf,        IO::XDMF::Center::Cell);
  phiGrid.add(qAbs,        IO::XDMF::Center::Cell);
  phiGrid.add(qRel,        IO::XDMF::Center::Cell);
  phiGrid.add(phiMoved,    IO::XDMF::Center::Node);

  std::cout << "Two-wavy-circles merging LSR sweep on " << n << "x" << n
            << " unit-square mesh, " << nFrames << " frames\n";
  std::cout << "  R0=" << R0 << "  amp=" << amp << "  k=" << kLobes
            << "  orbit R=" << orbitR
            << "  merge-half-arc~" << std::asin(R0 / orbitR) << " rad\n";

  std::size_t framesConverged = 0;
  std::vector<Real> finalFitPerFrame;
  finalFitPerFrame.reserve(nFrames);

  for (std::size_t frame = 0; frame < nFrames; ++frame)
  {
    const Real t = static_cast<Real>(frame) / static_cast<Real>(nFrames);
    const Real angle = Real(2) * Real(M_PI) * t;

    // Circle A: CCW orbit. Circle B: CW orbit. They coincide at
    // angle=0 and angle=pi and are diametrically opposite at
    // angle=pi/2 and 3pi/2.
    TwoWavyCirclesLevelSet levelSet;
    levelSet.a.cx = Real(0.5) + orbitR * std::cos(angle);
    levelSet.a.cy = Real(0.5) + orbitR * std::sin(angle);
    levelSet.a.R0 = R0;  levelSet.a.amp = amp;
    levelSet.a.k = kLobes;
    levelSet.a.phase =  angle;

    levelSet.b.cx = Real(0.5) + orbitR * std::cos(-angle);
    levelSet.b.cy = Real(0.5) + orbitR * std::sin(-angle);
    levelSet.b.R0 = R0;  levelSet.b.amp = amp;
    levelSet.b.k = kLobes;
    levelSet.b.phase = -angle;
    levelSet.eps = kSoftMinEps;

    const Real centerDist =
      std::sqrt((levelSet.a.cx - levelSet.b.cx) *
                (levelSet.a.cx - levelSet.b.cx)
              + (levelSet.a.cy - levelSet.b.cy) *
                (levelSet.a.cy - levelSet.b.cy));
    // Head-on (centerDist near 0) vs maximally-separated lateral
    // configuration (centerDist near 2*orbitR). With orbitR=0.10 and
    // R0=0.20 the regions always overlap; the qualitative state names
    // refer to the centre-pair geometry, not to a topology change.
    const char* configName =
        centerDist < Real(0.5) * orbitR     ? "head-on"
      : centerDist < Real(1.5) * orbitR     ? "oblique"
      :                                       "lateral";

    std::cout << "\n--- Frame " << std::setw(3) << frame
              << " : angle=" << std::fixed << std::setprecision(3) << angle
              << " rad  d(a,b)=" << std::setprecision(3) << centerDist
              << "  config=" << configName << '\n';

    clearXDMFRegionAttributes(mesh);
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
      if (incident.size() != 2) continue;
      const auto itA = cellToLocal.find(incident[0]);
      const auto itB = cellToLocal.find(incident[1]);
      if (itA == cellToLocal.end() || itB == cellToLocal.end()) continue;
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

    // ---- Stage 2: FMM signed distance ----
    psi = Real(0);
    Distance::Eikonal<ScalarP1, Math::Vector<Real>>(psi)
      .setInterface(interfaceAttribute)
      .setInterior(interiorAttribute)
      .solve()
      .sign();

    // ---- Stage 3: smooth-band measure ----
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

    // ---- Stage 4a: LSR solve ----
    auto& cellCache = cellCacheBg;

    RealFunction phi(
        [&](const Geometry::Point& p) -> Real
        {
          const auto& c = p.getPhysicalCoordinates();
          return levelSet.phi(vec2(c(0), c(1)));
        });
    AnalyticVectorFunction gradPhi(
        [&](const Geometry::Point& p) -> Math::SpatialVector<Real>
        {
          const auto& c = p.getPhysicalCoordinates();
          return levelSet.grad(vec2(c(0), c(1)));
        },
        /*dimension=*/2);
    AnalyticMatrixFunction hessPhi(
        [&](const Geometry::Point& p) -> Math::SpatialMatrix<Real>
        {
          const auto& c = p.getPhysicalCoordinates();
          return levelSet.hess(vec2(c(0), c(1)));
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
    baseParams.jMinRatio = Real(1e-8);
    baseParams.jSafeRatio = Real(1e-3);
    baseParams.lineSearchSafetyMargin = kLineSearchSafetyMargin;
    baseParams.qRelMax = kQRelMax;
    baseParams.qBarrierWeight = kQBarrierWeight;
    baseParams.qBarrierAct    = kQBarrierAct;
    baseParams.qBarrierMax    = kQBarrierMax;
    // BEST-QUALITY profile (see LSR.h docstring + safety-net sweep):
    // γ_jBar=1, jBarSafe=0.5, γ_vol=0.01 minimises worst-Q_rel without
    // losing fit. All dimensionless via the normalizer · h_ref² factor
    // applied inside LSR.h.
    baseParams.jBarrierWeight = parseRealOption(
        argc, argv, "j-barrier-weight", Real(1.0));
    baseParams.jBarrierSafeRatio = parseRealOption(
        argc, argv, "j-barrier-safe", Real(0.5));
    baseParams.jVolumeTetherWeight = parseRealOption(
        argc, argv, "volume-tether-weight", Real(0.01));
    baseParams.alphaInit = kLineSearchAlphaInit;
    baseParams.alphaReduction = kLineSearchReduction;
    baseParams.alphaMin = kLineSearchAlphaMin;
    baseParams.absoluteTolerance = kResidualAbsTol;
    baseParams.relativeTolerance = kResidualRelTol;
    baseParams.maxNewtonIterations = 80;
    baseParams.stallPatience = kStallPatience;
    baseParams.initialGuessMetric = initialGuessMetric;
    baseParams.initialGuessScaling = initialGuessScaling;
    baseParams.initialGuessElasticityLambda = initialGuessElasticityLambda;
    baseParams.initialGuessElasticityMu = initialGuessElasticityMu;
    baseParams.tangent = tangentMode;

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

    Real        residualBest = std::numeric_limits<Real>::infinity();
    Real        r0           = -1;
    std::size_t itCount      = 0;
    Real        interfaceFit = 0;
    bool        convergedFlag = false;
    const char* convergedReason = "no";

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
    params.acceptedStateConvergenceTest =
      [&](const LSRReport&) -> bool
      {
        interfaceFit = computeInterfaceFit();
        return interfaceFit < kInterfaceFitTol;
      };

    LSR lsr(u);
    const LSRReport lsrReport =
      lsr.setParameters(params).solve(psi, phi, gradPhi, hessPhi);

    r0 = lsrReport.initialResidual;
    itCount = lsrReport.iterations;
    residualBest = lsrReport.finalResidual;
    if (interfaceFit == Real(0))
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

    // ---- Stage 5a: psi-side fields ----
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

    // ---- Stage 5b: phi (moved) fields ----
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

    auto [movedCache, movedToLocal] = precomputeCellGeometry(moved);
    Real qMin = std::numeric_limits<Real>::infinity();
    Real qMax = 0;
    Real qSum = 0;
    std::size_t qCount = 0;
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
        const Real q =                                      // = q_abs
          dst.A.squaredNorm()
            / (static_cast<Real>(d) * std::pow(sigDetAu, dExp));
        qAbs.getData()(dof) = q;
        const Math::SpatialMatrix<Real> F = dst.A * src.A.inverse();
        qRel.getData()(dof) =
          F.squaredNorm()
            / (static_cast<Real>(d) * std::pow(jKval, dExp));
        cellLabelPhi.getData()(dof) =
          static_cast<Real>(classified.labels[local]);
        // Printed Q[min/mean/max] diagnostic stays on q_abs (absolute
        // line-search cap which is an absolute quality bound).
        qMin = std::min(qMin, q);
        qMax = std::max(qMax, q);
        qSum += q;
        ++qCount;
      }
    }
    const Real qMean = qCount > 0 ? qSum / static_cast<Real>(qCount) : Real(0);
    phiMoved = [&](const Geometry::Point& p) -> Real
    {
      const auto& X = p.getCoordinates();
      return levelSet.phi(vec2(X(0), X(1)));
    };

    const char* strategyName =
        kInitialGuessStrategy == InitialGuessStrategy::Cold    ? "cold"
      : kInitialGuessStrategy == InitialGuessStrategy::Warm    ? "warm"
      :                                                          "hilbert";
    std::cout << "    Newton it=" << itCount
              << "  ||R||=" << std::scientific << std::setprecision(3)
              << residualBest
              << "  fit=" << std::setprecision(3) << interfaceFit
              << "  Q[min/mean/max]="
              << std::setprecision(2) << qMin << "/" << qMean << "/" << qMax
              << "  init=" << strategyName
              << "  converged=" << convergedReason << '\n';
    if (convergedFlag) ++framesConverged;

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
    for (Real x : finalFitPerFrame) fitMean += x;
    fitMean /= static_cast<Real>(finalFitPerFrame.size());
    std::cout << "  ||phi(X+u)||_RMS  min=" << std::scientific
              << std::setprecision(3) << fitMin
              << "  mean=" << fitMean
              << "  max=" << fitMax << '\n';
  }

  return framesConverged == nFrames ? 0 : 1;
}
