/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//
// LSR (level-set registration) reconstruction example — single-frame demo.
//
// Pipeline
//   1. Build an n x n triangular mesh of the unit square.
//   2. Classify cells using a min-s/t cut on phase moments of the
//      analytic level set (WavyCircleLevelSet).
//   3. Tag classifier-cut facets and solve the screened-Poisson
//      (I - psiEll^2 Delta) psi = +-M with psi = 0 on the cut skeleton,
//      then rescale psi by a narrow-band normal-gradient calibration.
//      psi is the smooth topological reference field.
//   4. Construct phi/gradPhi/hessPhi as plain analytic operands of
//      WavyCircleLevelSet. Adaptation::LSR owns the deformed evaluation
//      and tracing semantics. phi is the geometric reference field.
//   5. Solve the LSR problem via Adaptation::LSR (the facade) with the
//      "always converges" recipe described below.
//   6. Build the moved mesh, write the psi/phi XDMF grids with per-cell
//      diagnostics (j, q_abs, q_rel), print summary.
//
// "Always converges" recipe (hard-coded in main; no runtime switch)
//   - Tangent     : LSRTangent::PSDProjectedNewton.
//                   Full Newton diverges at the medial axis of phi; pure
//                   GaussNewton converges but linearly. PSD-projected
//                   gives the best of both.
//   - Initial guess: LSRInitialGuess::Hilbert with
//                    LSRHilbertMetric::Harmonic.
//                    The Hilbert lift gets the first step in close to
//                    one and skips most line-search backtracks.
//   - Quadrature  : lsrQuadOrderFor(feOrder) = 4 * feOrder (Rodin core).
//                   Bigger than the polynomial-exactness lower bound;
//                   makes the basis-at-qpt sampling map full-rank on the
//                   P2 internal modes.
//   - Sampled geometric energies: JacobianAdmissibilityBarrierSampled
//                   assembles E_shape, the optional Q_rel barrier, and the
//                   optional E_admissibility barrier on j = det(I + grad u).
//                   Built into the facade; no wiring needed here.
//   - Q_rel cap   : qRelMax = +infinity (off). Available but unnecessary.
//   - Shape weight: gamma = 1e-1.
//   - H1 reg      : disabled. The sampled barrier + quadrature bump are
//                   enough; explicit smoothing is not needed.
//   - Band width  : deltaW ~= 2.6 h (factor 1.5 * 1.75 baked in).
//   - Cut weight  : lambdaC = 0.008 (classifier, not LSR).
//
// With these defaults the example converges 3 Newton iterations on P1
// (fit ~1.3e-3 at n=50, h^2-rate) and ~10 iterations on P2.
//
// Scope and honest caveats
//   - Planar 2D, affine triangle cells. The P1 vector displacement is the
//     default; LevelSetLSRReconstructionP2 (compiled with
//     -DRODIN_LSR_P2_DISPLACEMENT) is the H1<2> vector variant.
//   - The classifier and phase-moment computation use a triangle
//     barycentric quadrature; this part of the example is NOT
//     FES-independent. The LSR integrators that consume the result are.
//   - phi here is analytic for clarity. Deformed evaluation is handled in
//     Adaptation::LSR, so replacing phi by a discrete field does not alter
//     the example-level solve call.
//   - Full Newton tangent is unsafe: phi's Hessian
//       hess(phi) = (I - n n^T) / ||p - c||
//     is singular on the medial axis (the centre for a disk). PSD
//     projection clips the indefinite piece and converges; raw Newton
//     drives a handful of cells past j_min and stalls. This is documented
//     in the file-header history rather than exposed as a runtime flag.
//
#include <Rodin/Adaptation.h>
#include <Rodin/Advection/Lagrangian.h>
#include <Rodin/Assembly.h>
#include <Rodin/Distance/Eikonal.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Math.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>
#include <Rodin/Solver/NewtonSolver.h>
#include <Rodin/Solver/SparseLU.h>
#include <Rodin/Variational.h>
#include <Rodin/Variational/IntegrationPoint.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
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
  using Clock = std::chrono::steady_clock;

  Vec2 vec2(Real x = 0, Real y = 0)
  {
    Vec2 out(2);
    out(0) = x;
    out(1) = y;
    return out;
  }

  Real elapsedSeconds(Clock::time_point start)
  {
    return std::chrono::duration<Real>(Clock::now() - start).count();
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
        const Real ux =
          uData(uFes.getGlobalIndex({D, cell.getIndex()}, a * 2));
        const Real uy =
          uData(uFes.getGlobalIndex({D, cell.getIndex()}, a * 2 + 1));
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
  // Single-frame wavy circle level set.
  // -------------------------------------------------------------------------
  struct WavyCircleLevelSet
  {
    Real cx = Real(0.5);
    Real cy = Real(0.5);
    Real R0 = Real(0.20);
    Real amp = Real(0.05);
    Real lobes = Real(6);
    Real phase = Real(0);

    Real phi(const Vec2& p) const
    {
      const Real dx = p(0) - cx;
      const Real dy = p(1) - cy;
      const Real r = std::sqrt(dx * dx + dy * dy);
      const Real theta = std::atan2(dy, dx);
      const Real R = R0 + amp * std::cos(lobes * theta + phase);
      return r - R;
    }

    Vec2 grad(const Vec2& p) const
    {
      const Real dx = p(0) - cx;
      const Real dy = p(1) - cy;
      const Real r2 = dx * dx + dy * dy;
      const Real r = std::max(std::sqrt(r2), Real(1e-14));
      const Real theta = std::atan2(dy, dx);
      const Real dRdtheta = -amp * lobes * std::sin(lobes * theta + phase);
      const Real r2safe = std::max(r2, Real(1e-28));
      return vec2(
          dx / r - dRdtheta * (-dy / r2safe),
          dy / r - dRdtheta * ( dx / r2safe));
    }

    Math::SpatialMatrix<Real> hess(const Vec2& p) const
    {
      const Real dx = p(0) - cx;
      const Real dy = p(1) - cy;
      const Real r2 = std::max(dx * dx + dy * dy, Real(1e-28));
      const Real r = std::max(std::sqrt(r2), Real(1e-14));
      const Real r3 = r * r * r;
      const Real r4 = r2 * r2;
      const Real theta = std::atan2(dy, dx);
      const Real angle = lobes * theta + phase;
      const Real dRdtheta = -amp * lobes * std::sin(angle);
      const Real d2Rdtheta2 = -amp * lobes * lobes * std::cos(angle);

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

      const Vec2 gradTheta = vec2(-dy / r2, dx / r2);
      Math::SpatialMatrix<Real> GG(2, 2);
      GG(0, 0) = gradTheta(0) * gradTheta(0);
      GG(0, 1) = gradTheta(0) * gradTheta(1);
      GG(1, 0) = GG(0, 1);
      GG(1, 1) = gradTheta(1) * gradTheta(1);

      Math::SpatialMatrix<Real> H(2, 2);
      H = Hr - dRdtheta * Ht - d2Rdtheta2 * GG;
      return H;
    }
  };

  // -------------------------------------------------------------------------
  // Triangle-barycentric quadrature used by the (triangle-only) phase
  // moment computation. NOT FES-independent.
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

  template <class LevelSet>
  std::vector<CellMomentInfo> collectCellMomentInfo(
      const Mesh<Context::Local>& mesh,
      const LevelSet& levelSet,
      Real epsilon)
  {
    std::vector<CellMomentInfo> cells;
    cells.reserve(mesh.getCellCount());

    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const auto& cellPolytope = *cellIt;
      const auto& vertices = cellPolytope.getVertices();
      if (vertices.size() != 3)
        throw std::runtime_error(
            "LevelSetLSRReconstruction expects triangular cells.");

      CellMomentInfo info;
      info.index = cellPolytope.getIndex();
      for (size_t i = 0; i < 3; ++i)
      {
        info.vertices[i] = vertices[i];
        info.x[i] = mesh.getVertexCoordinates(vertices[i]);
      }

      const Vec2 e1 = info.x[1] - info.x[0];
      const Vec2 e2 = info.x[2] - info.x[0];
      const Real signedArea =
        Real(0.5) * (e1(0) * e2(1) - e1(1) * e2(0));
      info.area = std::abs(signedArea);

      Real moment = 0;
      for (const auto& bary : TriangleBarycentricQuadrature)
      {
        const Vec2 xq = interpolateVec(info.x, bary);
        moment += applyPhaseMomentMap(levelSet.phi(xq), epsilon);
      }
      info.moment = moment / TriangleBarycentricQuadrature.size();
      cells.push_back(std::move(info));
    }
    return cells;
  }

  Vec2 cellP1Gradient(
      const std::array<Vec2, 3>& x,
      const std::array<Real, 3>& value)
  {
    const Vec2 e1 = x[1] - x[0];
    const Vec2 e2 = x[2] - x[0];
    const Real det = e1(0) * e2(1) - e1(1) * e2(0);
    if (std::abs(det) <= Real(1e-30))
      return vec2(0, 0);
    const Real dv1 = value[1] - value[0];
    const Real dv2 = value[2] - value[0];
    return vec2(
        (dv1 * e2(1) - e1(1) * dv2) / det,
        (e1(0) * dv2 - dv1 * e2(0)) / det);
  }

  Real facetLength(const Mesh<Context::Local>& mesh, Index facet)
  {
    const auto face = mesh.getFace(facet);
    const auto& vertices = face->getVertices();
    const Vec2 a = mesh.getVertexCoordinates(vertices[0]);
    const Vec2 b = mesh.getVertexCoordinates(vertices[1]);
    return (b - a).norm();
  }

  // -------------------------------------------------------------------------
  // Convergence-rate diagnostics. Given a non-increasing residual history,
  // estimate the local order
  //
  //   q_n = log(r_{n+1}/r_n) / log(r_n/r_{n-1})
  //
  // which approaches 1 for linear convergence, > 1 for super-linear, and
  // 2 for quadratic.
  // -------------------------------------------------------------------------
  void printConvergenceOrders(const std::vector<Real>& residuals)
  {
    std::cout << "\nObserved local convergence order q[n]:\n";
    for (std::size_t n = 1; n + 1 < residuals.size(); ++n)
    {
      const Real rPrev = residuals[n - 1];
      const Real rThis = residuals[n];
      const Real rNext = residuals[n + 1];
      if (rPrev <= 0 || rThis <= 0 || rNext <= 0)
        continue;
      const Real q =
        std::log(rNext / rThis) / std::log(rThis / rPrev);
      std::cout << "  q[" << std::setw(2) << (n + 1) << "] = "
                << std::fixed << std::setprecision(3) << q
                << "  (||R||_{n+1}/||R||_n = "
                << std::scientific << std::setprecision(3)
                << rNext / rThis << ")\n";
    }
  }

  bool hasFlag(int argc, char** argv, const std::string& name)
  {
    const std::string flag = "--" + name;
    for (int i = 1; i < argc; ++i)
      if (argv[i] == flag)
        return true;
    return false;
  }

  Real getOptionReal(
      int argc,
      char** argv,
      const std::string& name,
      Real defaultValue)
  {
    const std::string prefix = "--" + name + "=";
    for (int i = 1; i < argc; ++i)
    {
      const std::string arg(argv[i]);
      if (arg.rfind(prefix, 0) == 0)
        return std::stod(arg.substr(prefix.size()));
    }
    return defaultValue;
  }

  std::string getOptionString(
      int argc,
      char** argv,
      const std::string& name,
      const std::string& defaultValue)
  {
    const std::string prefix = "--" + name + "=";
    for (int i = 1; i < argc; ++i)
    {
      const std::string arg(argv[i]);
      if (arg.rfind(prefix, 0) == 0)
        return arg.substr(prefix.size());
    }
    return defaultValue;
  }

  std::size_t getOptionSizeT(
      int argc,
      char** argv,
      const std::string& name,
      std::size_t defaultValue)
  {
    const std::string prefix = "--" + name + "=";
    for (int i = 1; i < argc; ++i)
    {
      const std::string arg(argv[i]);
      if (arg.rfind(prefix, 0) == 0)
        return static_cast<std::size_t>(std::stoul(arg.substr(prefix.size())));
    }
    return defaultValue;
  }

}

int main(int argc, char** argv)
{
  // ---------------------------------------------------------------------------
  // Knobs.
  //
  // Compile-time:
  //   LevelSetLSRReconstruction     -> P1 vector displacement (default).
  //   LevelSetLSRReconstructionP2   -> H1<2> (P2) vector displacement;
  //                                    compile with -DRODIN_LSR_P2_DISPLACEMENT.
  //
  // Runtime:
  //   --n=<N>              mesh is N x N (default 50).
  //   --newton-steps=<N>   maximum Newton steps, default 40.
  //   --R0=<r>             base radius, default 0.20.
  //   --amp=<a>            radial wave amplitude, default 0.05.
  //   --lobes=<m>          radial wave lobe count, default 6.
  //   --phase=<theta>      radial wave phase, default 0.
  //   --flow-trace         evaluate phi through Flow tracing inside LSR.
  //   --energy=push|pull|push-swapped|push-swapped-forward
  //                         data term in the LSR sum:
  //                          push (default): (phi(X+u) - psi(X))^2,
  //                          pull          : (phi(X) - psi(X-u))^2,
  //                          push-swapped  : solve (psi(X+v)-phi(X))^2
  //                                          and report u=-v,
  //                          push-swapped-forward
  //                                        : use u=-v as warm start, then
  //                                          solve the default push problem.
  //   --tangent=gn|newton|psd
  //                         LSR data tangent, default psd.
  //   --initial-guess=hilbert|zero
  //                         zero selects a cold start, default hilbert.
  //   --pull-psi-grad=smoothed|raw
  //                         derivative source for E_pull, default smoothed.
  //   --hess-smoothing=none|l2|h1
  //                         optional componentwise Hessian smoothing.
  //   --hess-smooth-ell=<l> length for h1 Hessian smoothing, default h.
  //   --shape-weight=<w>   weight of E_shape, default 1e-1.
  //   --gamma-weight=<w>   weight of E_Gamma, default 1.
  //   --j-barrier-weight=<w>
  //                         weight of E_admissibility, default 0.
  //   --j-barrier-safe=<j> activation ratio for E_admissibility, default 0.
  //   --volume-tether-weight=<w>
  //                         weight of 0.5 (log j)^2, default 0.
  //   --no-alpha-predictor disable the sampled quadratic alpha predictor.
  //
  // Everything else is hard-coded below to the "always converges" recipe.
  // The file-header comment lists the choices and why they are the defaults.
  // ---------------------------------------------------------------------------
  const std::size_t n = getOptionSizeT(argc, argv, "n", 50);
  const Real h = Real(1) / static_cast<Real>(n - 1);
  const Real epsilon = 1.25 * h;
  const Real lambdaC = 0.008;
  constexpr Real cx = Real(0.5);
  constexpr Real cy = Real(0.5);
  const Real R0 = getOptionReal(argc, argv, "R0", Real(0.20));
  const Real amp = getOptionReal(argc, argv, "amp", Real(0.05));
  const Real lobes = getOptionReal(argc, argv, "lobes", Real(6));
  const Real phase = getOptionReal(argc, argv, "phase", Real(0));
  const Real psiEll = Real(2) * h;
  const Real psiSourceMagnitude = Real(5) * psiEll;
  const std::size_t maxNewtonSteps =
    getOptionSizeT(argc, argv, "newton-steps", 40);
  const Real shapeWeight = getOptionReal(argc, argv, "shape-weight", Real(1e-1));
  const Real gammaWeight = getOptionReal(argc, argv, "gamma-weight", Real(1));
  const Real jBarrierWeight =
    getOptionReal(argc, argv, "j-barrier-weight", Real(0));
  const Real jBarrierSafeRatio =
    getOptionReal(argc, argv, "j-barrier-safe", Real(0));
  const Real jVolumeTetherWeight =
    getOptionReal(argc, argv, "volume-tether-weight", Real(0));
  enum class DataEnergy { Push, Pull, PushSwapped, PushSwappedForward };
  const std::string energyStr =
    getOptionString(argc, argv, "energy", "push");
  const DataEnergy dataEnergy =
      energyStr == "pull"
        ? DataEnergy::Pull
    : energyStr == "push"
        ? DataEnergy::Push
    : energyStr == "push-swapped"
        ? DataEnergy::PushSwapped
    : energyStr == "push-swapped-forward"
        ? DataEnergy::PushSwappedForward
    : throw std::runtime_error(
        "Unknown --energy. Expected push, pull, push-swapped, "
        "or push-swapped-forward.");
  const std::string tangentStr =
    getOptionString(argc, argv, "tangent", "psd");
  const LSRTangent requestedTangent =
      tangentStr == "gn" || tangentStr == "gauss-newton"
        ? LSRTangent::GaussNewton
    : tangentStr == "newton" || tangentStr == "full-newton"
        ? LSRTangent::Newton
    : tangentStr == "psd" || tangentStr == "psd-projected-newton"
        ? LSRTangent::PSDProjectedNewton
    : throw std::runtime_error(
        "Unknown --tangent. Expected gn, newton, or psd.");
  const std::string initialGuessStr =
    getOptionString(argc, argv, "initial-guess", "hilbert");
  const LSRInitialGuess requestedInitialGuess =
      initialGuessStr == "zero" || initialGuessStr == "cold"
        ? LSRInitialGuess::Zero
    : initialGuessStr == "current" || initialGuessStr == "keep"
        ? LSRInitialGuess::Current
    : initialGuessStr == "hilbert"
        ? LSRInitialGuess::Hilbert
    : throw std::runtime_error(
        "Unknown --initial-guess. Expected hilbert, zero, or current.");
  enum class PullPsiDerivative { Smoothed, Raw };
  const std::string pullPsiGradStr =
    getOptionString(argc, argv, "pull-psi-grad", "smoothed");
  const PullPsiDerivative pullPsiDerivative =
      pullPsiGradStr == "smoothed"
        ? PullPsiDerivative::Smoothed
    : pullPsiGradStr == "raw"
        ? PullPsiDerivative::Raw
    : throw std::runtime_error(
        "Unknown --pull-psi-grad. Expected smoothed or raw.");
  enum class HessianSmoothing { None, L2, H1 };
  const std::string hessSmoothingStr =
    getOptionString(argc, argv, "hess-smoothing", "none");
  const HessianSmoothing hessSmoothing =
      hessSmoothingStr == "none"
        ? HessianSmoothing::None
    : hessSmoothingStr == "l2"
        ? HessianSmoothing::L2
    : hessSmoothingStr == "h1"
        ? HessianSmoothing::H1
    : throw std::runtime_error(
        "Unknown --hess-smoothing. Expected none, l2, or h1.");
  const Real hessSmoothEll =
    getOptionReal(argc, argv, "hess-smooth-ell", h);

  // psi reconstruction.
  //   poisson : screened-Poisson lift of a +/-M sign source. Smooth field
  //             with |grad psi| ~ M / psiEll, not a true SDF.
  //   fmm           : Fast Marching Method (Distance::Eikonal). True SDF
  //                   reconstruction with |grad psi| = 1 a.e.
  //   projected-phi    : L2 projection of the input phi with homogeneous
  //                      Dirichlet data on the classified skeleton.
  //   projected-phi-h1 : H1 projection of the input phi with the same
  //                      skeleton constraint.
  enum class PsiReconstruction { Poisson, Fmm, ProjectedPhi, ProjectedPhiH1 };
  const std::string psiStr =
    getOptionString(argc, argv, "psi", "poisson");
  const PsiReconstruction psiReconstruction =
      psiStr == "fmm"
        ? PsiReconstruction::Fmm
    : psiStr == "poisson"
        ? PsiReconstruction::Poisson
    : psiStr == "projected-phi" || psiStr == "projection"
        ? PsiReconstruction::ProjectedPhi
    : psiStr == "projected-phi-h1" || psiStr == "projection-h1"
        ? PsiReconstruction::ProjectedPhiH1
    : throw std::runtime_error(
        "Unknown --psi. Expected poisson, fmm, projected-phi, "
        "or projected-phi-h1.");
  const Real psiProjectEll =
    getOptionReal(argc, argv, "psi-project-ell", h);
  enum class PhiValueSource { Analytic, Screened };
  const std::string phiStr =
    getOptionString(argc, argv, "phi", "analytic");
  const PhiValueSource phiValueSource =
    (phiStr == "screened")
      ? PhiValueSource::Screened
      : (phiStr == "analytic")
          ? PhiValueSource::Analytic
          : throw std::runtime_error("Unknown --phi. Expected analytic or screened.");

  // -----  One-step advection test  ---------------------------------------
  //
  //   --advect=on|off                 (default off)
  //   --advect-vx=<value>             constant x velocity (default 0.05)
  //   --advect-vy=<value>             constant y velocity (default 0)
  //   --advect-dt=<value>             single semi-Lagrangian step length
  //                                    (default 1.0, so total displacement
  //                                     is (vx, vy))
  //
  // When `--advect=on`:
  //   1. The analytic level set is sampled at its INITIAL position
  //      (--R0, --amp, --lobes, --phase, centred at the domain centre)
  //      to a P1 grid function phi_h_0.
  //   2. phi_h_0 is advected by one semi-Lagrangian step under the
  //      constant velocity v = (vx, vy) for time dt, giving phi_h_adv.
  //   3. The classifier, psi reconstruction, and LSR data term ALL
  //      consume phi_h_adv instead of the analytic level set. Tangent
  //      is forced to Gauss-Newton (no discrete Hess available on P1).
  //   4. The ANALYTIC reference for the fit metric is the level set at
  //      the advected position (cx + vx*dt, cy + vy*dt), i.e. the
  //      ground-truth geometry that phi_h_adv approximates.
  //
  // The interface fit then reports
  //      ||phi_GT(X + u)||_RMS  on Gamma_psi,h
  // where phi_GT is the analytic level set at the advected position --
  // this is the honest measure of "does the LSR recover the true
  // geometry from a once-advected discrete phi?".
  const bool advectOn =
    getOptionString(argc, argv, "advect", "off") == "on";
  const Real advectVx = getOptionReal(argc, argv, "advect-vx", Real(0.05));
  const Real advectVy = getOptionReal(argc, argv, "advect-vy", Real(0));
  const Real advectDt = getOptionReal(argc, argv, "advect-dt", Real(1));
  const bool advectDiagnostics = hasFlag(argc, argv, "advect-diagnostics");

  std::cout << std::unitbuf;
  const auto totalStart = Clock::now();
  auto phaseStart = Clock::now();

  // Tangent: PSD-projected Newton (LSRTangent::PSDProjectedNewton) is the
  // canonical mode for this prototype. Full Newton is unsafe because phi's
  // Hessian is indefinite at the medial axis; GaussNewton converges but
  // slower. PSD is what the facade picks by default.

  // -------------------------------------------------------------------------
  // Step 1-3: mesh, classification, attribute tagging.
  // -------------------------------------------------------------------------
  LocalMesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { n, n });
  mesh.scale(h);
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 2);
  mesh.getConnectivity().compute(2, 2);
  mesh.getConnectivity().compute(0, 0);

  // Analytic level set at the INITIAL position (source for the optional
  // semi-Lagrangian advection step).
  WavyCircleLevelSet levelSetInit;
  levelSetInit.cx = cx;
  levelSetInit.cy = cy;
  levelSetInit.R0 = R0;
  levelSetInit.amp = amp;
  levelSetInit.lobes = lobes;
  levelSetInit.phase = phase;

  // Ground-truth analytic level set at the ADVECTED position. With
  // advect=off this equals levelSetInit; with advect=on the centre is
  // translated by (vx, vy) * dt under the prescribed constant
  // velocity. The lobe phase is unchanged because rigid translation
  // does not rotate the wavy-circle parameter space.
  WavyCircleLevelSet levelSet = levelSetInit;
  if (advectOn)
  {
    levelSet.cx = cx + advectVx * advectDt;
    levelSet.cy = cy + advectVy * advectDt;
  }

  // Build the discrete advected phi_h_adv on a scalar FES of
  // `mesh`:
  //   1. phi_h_0  = nodal interpolant of levelSetInit.
  //   2. phi_h_adv = one semi-Lagrangian step of phi_h_0 under
  //                  v(x) = (vx, vy), step length advectDt.
  // P1 is the default. Compile-time variants RODIN_LSR_PHI_P2 and
  // RODIN_LSR_PHI_P3 select H1<2> and H1<3> scalar interpolation.
#if defined(RODIN_LSR_PHI_P3)
  using ScalarPhiFes = H1<3, Real, LocalMesh>;
  const char* scalarPhiFesName = "advected P3 GridFunction";
#elif defined(RODIN_LSR_PHI_P2)
  using ScalarPhiFes = H1<2, Real, LocalMesh>;
  const char* scalarPhiFesName = "advected P2 GridFunction";
#else
  using ScalarPhiFes = P1<Real, LocalMesh>;
  const char* scalarPhiFesName = "advected P1 GridFunction";
#endif
#if defined(RODIN_LSR_PHI_P3)
  ScalarPhiFes phiHFes(std::integral_constant<std::size_t, 3>{}, mesh);
#elif defined(RODIN_LSR_PHI_P2)
  ScalarPhiFes phiHFes(std::integral_constant<std::size_t, 2>{}, mesh);
#else
  ScalarPhiFes phiHFes(mesh);
#endif
  GridFunction phiH0(phiHFes);
  GridFunction phiHAdv(phiHFes);

  // ---  Gradient/Hessian source for the LSR data integrand  ---------------
  //
  //   --phi-grad=smoothed  (default when advect=on)
  //       gradPhi = L2-projection of Grad(phi_h_adv) onto a P1 vector
  //                 finite-element space (Zienkiewicz-Zhu-style nodal
  //                 recovery on the same mesh).
  //       hessPhi = Jacobian(grad_smoothed)   (P0 from a P1 vector field)
  //       Tangent = PSDProjectedNewton.
  //
  //   --phi-grad=discrete
  //       gradPhi = Grad(phi_h_adv)
  //       hessPhi = 0
  //       Tangent = GN.
  //
  //   --phi-grad=analytic
  //       gradPhi = levelSet.grad(X)      (analytic at advected position)
  //       hessPhi = levelSet.hess(X)
  //       Tangent = PSDProjectedNewton.
  enum class PhiGradSource { Discrete, Analytic, Smoothed };
  const std::string phiGradStr =
    getOptionString(argc, argv, "phi-grad", advectOn ? "smoothed" : "analytic");
  const PhiGradSource phiGradSource =
      phiGradStr == "discrete" ? PhiGradSource::Discrete
    : phiGradStr == "analytic" ? PhiGradSource::Analytic
    : phiGradStr == "smoothed" ? PhiGradSource::Smoothed
    : throw std::runtime_error(
          "Unknown --phi-grad. Expected discrete, analytic, or smoothed.");

  // Vector P1 grid function for the smoothed gradient. Allocated
  // unconditionally for type stability; populated only when
  // phiGradSource == Smoothed.
  using VectorP1Phi = P1<Math::SpatialVector<Real>, LocalMesh>;
  VectorP1Phi gradPhiFes(mesh, /*vdim=*/2);
  GridFunction gradPhiSmoothed(gradPhiFes);
  gradPhiSmoothed.getData().setZero();
  if (advectOn)
  {
    // (1) Sample analytic levelSetInit onto P1.
    phiH0 = [&](const Geometry::Point& p) -> Real
    {
      const auto& X = p.getCoordinates();
      return levelSetInit.phi(vec2(X(0), X(1)));
    };
    // (2) Single semi-Lagrangian step under constant v.
    AnalyticVectorFunction velocity(
        [&](const Geometry::Point&) -> Math::SpatialVector<Real>
        {
          return vec2(advectVx, advectVy);
        },
        /*dimension=*/2);
    TrialFunction advectTrial(phiHFes);
    TestFunction  advectTest(phiHFes);
    Advection::Lagrangian adv(advectTrial, advectTest, phiH0, velocity);
    adv.step(advectDt);
    phiHAdv.getData() = advectTrial.getSolution().getData();

    if (advectDiagnostics)
    {
      auto gradPhiHAdvDiag = Grad(phiHAdv);
      const auto coeffDiff = phiHAdv.getData() - phiH0.getData();
      const Real coeffNorm = phiH0.getData().norm();
      Real measure = 0;
      Real err0L2 = 0;
      Real errAdvL2 = 0;
      Real advDiffL2 = 0;
      Real maxErr0 = 0;
      Real maxErrAdv = 0;
      Real maxAdvDiff = 0;
      Real bandMeasure = 0;
      Real bandGradSum = 0;
      Real bandGradExactSum = 0;
      Real bandGradMin = std::numeric_limits<Real>::infinity();
      Real bandGradMax = 0;
      const Real bandRadius = Real(2) * epsilon;
      for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
      {
        const auto& cell = *cellIt;
        const auto& qf =
          QF::PolytopeQuadratureFormula::get(8, cell.getGeometry());
        const auto& quadrature = cell.getQuadrature(qf);
        for (std::size_t q = 0; q < quadrature.getSize(); ++q)
        {
          const auto& pt = quadrature.getPoint(q);
          const Variational::IntegrationPoint ip(pt, &qf, q);
          const auto& X = pt.getPhysicalCoordinates();
          const Vec2 x = vec2(X(0), X(1));
          const Real w = qf.getWeight(q) * pt.getDistortion();
          const Real val0 = phiH0.getValue(ip);
          const Real valAdv = phiHAdv.getValue(ip);
          const Real exact0 = levelSetInit.phi(x);
          const Real exactAdv = levelSet.phi(x);
          const Real err0 = val0 - exact0;
          const Real errAdv = valAdv - exactAdv;
          const Real advDiff = valAdv - val0;
          measure += w;
          err0L2 += w * err0 * err0;
          errAdvL2 += w * errAdv * errAdv;
          advDiffL2 += w * advDiff * advDiff;
          maxErr0 = std::max(maxErr0, std::abs(err0));
          maxErrAdv = std::max(maxErrAdv, std::abs(errAdv));
          maxAdvDiff = std::max(maxAdvDiff, std::abs(advDiff));
          if (std::abs(valAdv) <= bandRadius)
          {
            const Real gn = gradPhiHAdvDiag.getValue(ip).norm();
            const Real gne = levelSet.grad(x).norm();
            bandMeasure += w;
            bandGradSum += w * gn;
            bandGradExactSum += w * gne;
            bandGradMin = std::min(bandGradMin, gn);
            bandGradMax = std::max(bandGradMax, gn);
          }
        }
      }
      std::cout << "  advect diagnostics:\n"
                << "    ||coeff(phi_h_adv-phi_h0)|| / ||coeff(phi_h0)|| = "
                << std::scientific << std::setprecision(6)
                << coeffDiff.norm() / std::max(coeffNorm, Real(1e-30))
                << '\n'
                << "    L2(phi_h0-phi_exact0) = "
                << std::sqrt(err0L2 / std::max(measure, Real(1e-30)))
                << "  max=" << maxErr0 << '\n'
                << "    L2(phi_h_adv-phi_exact_adv) = "
                << std::sqrt(errAdvL2 / std::max(measure, Real(1e-30)))
                << "  max=" << maxErrAdv << '\n'
                << "    L2(phi_h_adv-phi_h0) = "
                << std::sqrt(advDiffL2 / std::max(measure, Real(1e-30)))
                << "  max=" << maxAdvDiff << '\n';
      if (bandMeasure > Real(0))
      {
        std::cout << "    band |grad phi_h_adv| mean/min/max = "
                  << bandGradSum / bandMeasure
                  << " / " << bandGradMin
                  << " / " << bandGradMax
                  << "  exact mean="
                  << bandGradExactSum / bandMeasure
                  << '\n';
      }
    }

    // ---- L2-projection of Grad(phi_h_adv) onto a continuous P1 vector
    //      grid function (smoothed-gradient experiment).
    if (phiGradSource == PhiGradSource::Smoothed)
    {
      TrialFunction gtrial(gradPhiFes);
      TestFunction  gtest(gradPhiFes);
      Problem gpb(gtrial, gtest);
      gpb = Integral(gtrial, gtest)
          - Integral(Grad(phiHAdv), gtest);
      Solver::SparseLU(gpb).solve();
      gradPhiSmoothed.getData() = gtrial.getSolution().getData();
    }
  }

  // Phase moments. With advect=off we just call the existing helper on
  // the analytic level set. With advect=on we build CellMomentInfo
  // entries from phi_h_adv evaluated at cell+ref-coord, so the
  // classifier sees exactly the same field that LSR consumes.
  std::vector<CellMomentInfo> cellMoments;
  if (!advectOn)
  {
    cellMoments = collectCellMomentInfo(mesh, levelSet, epsilon);
  }
  else
  {
    cellMoments.reserve(mesh.getCellCount());
    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const auto& cellPolytope = *cellIt;
      const auto& vertices = cellPolytope.getVertices();
      if (vertices.size() != 3)
        throw std::runtime_error(
            "LevelSetLSRReconstruction (advect=on) expects triangular cells.");
      CellMomentInfo info;
      info.index = cellPolytope.getIndex();
      for (std::size_t i = 0; i < 3; ++i)
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
        moment += applyPhaseMomentMap(phiHAdv.getValue(pt), epsilon);
      }
      info.moment = moment / TriangleBarycentricQuadrature.size();
      cellMoments.push_back(std::move(info));
    }
  }

  // Explicit cell-index <-> local-index maps for everything downstream.
  // We never assume mesh cell index == index into cellMoments / cellCache /
  // classified.labels. Build the maps once from the iteration order in
  // collectCellMomentInfo, then use them everywhere.
  std::unordered_map<Index, std::size_t> cellToLocal;
  std::vector<Index> localToCell;
  cellToLocal.reserve(cellMoments.size());
  localToCell.reserve(cellMoments.size());
  for (std::size_t local = 0; local < cellMoments.size(); ++local)
  {
    cellToLocal[cellMoments[local].index] = local;
    localToCell.push_back(cellMoments[local].index);
  }

  // Build graph arrays in LOCAL index space.
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
    if (incident.size() == 2)
    {
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
  }

  const MinSTCut cut;
  const MinSTCut::Result classified =
    cut.classify(volumes, moments, graphEdges);
  std::cout << "Phase classification: "
            << elapsedSeconds(phaseStart) << " s\n";
  phaseStart = Clock::now();

  std::vector<Index> interfaceFacets;
  interfaceFacets.reserve(classified.cutEdges.size());
  for (const MinSTCut::Edge& edge : classified.cutEdges)
    if (edge.index != MinSTCut::InvalidIndex)
      interfaceFacets.push_back(edge.index);

  constexpr Attribute interiorAttribute  = 1;
  constexpr Attribute exteriorAttribute  = 2;
  constexpr Attribute interfaceAttribute = 10;
  constexpr Attribute boundaryAttribute  = 20;
  auto applyClassificationAttributes =
    [&](const auto& labels,
        const std::vector<Index>& cells,
        const std::vector<Index>& facets)
    {
      const std::size_t D = mesh.getDimension();
      for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
        mesh.setAttribute({D, cellIt->getIndex()}, Attribute{0});
      for (auto faceIt = mesh.getFace(); faceIt; ++faceIt)
        mesh.setAttribute({D - 1, faceIt->getIndex()}, Attribute{0});

      for (std::size_t local = 0; local < labels.size(); ++local)
      {
        const Index cellIdx = cells[local];
        mesh.setAttribute(
            {D, cellIdx},
            labels[local] == MinSTCut::Inside
              ? interiorAttribute
              : exteriorAttribute);
      }
      for (const Index facet : facets)
        mesh.setAttribute({D - 1, facet}, interfaceAttribute);
      for (auto faceIt = mesh.getBoundary(); faceIt; ++faceIt)
        mesh.setAttribute({D - 1, faceIt->getIndex()}, boundaryAttribute);
    };
  applyClassificationAttributes(
      classified.labels, localToCell, interfaceFacets);

  mesh.save("LevelSetLSRReconstruction_psi.mesh", IO::FileFormat::MFEM);

  const Real delta = 1.75 * h;
  const Real deltaW = 1.5 * delta;

  // -------------------------------------------------------------------------
  // Step 4: psi (topological reference level set) via screened Poisson.
  // The source is negative on classified interior cells and positive on
  // exterior cells, while psi = 0 is imposed on the classified interface
  // skeleton. This keeps psi tied to the graph-cut topology without using
  // signed-distance reconstruction as an assumption.
  // -------------------------------------------------------------------------
  using ScalarFES = P1<Real, LocalMesh>;
  ScalarFES psiFes(mesh);
  GridFunction psi(psiFes);
  Real psiScale = 1;
  Real psiScaleNumerator = 0;
  Real psiScaleDenominator = 0;
  Real psiScaleSignMoment = 0;
  if (psiReconstruction == PsiReconstruction::Fmm)
  {
    // True SDF reconstruction: solve |grad psi| = 1 with psi = 0 on
    // Gamma_psi, then sign by classifier label. Yields |grad psi| ~ 1
    // a.e., so the LSR data term (whether push or pull) inherits the
    // same scale as analytic SDF phi -- no calibration needed.
    psi = Real(0);
    Distance::Eikonal<ScalarFES, Math::Vector<Real>>(psi)
      .setInterface(interfaceAttribute)
      .setInterior(interiorAttribute)
      .solve()
      .sign();
  }
  else if (psiReconstruction == PsiReconstruction::ProjectedPhi
           || psiReconstruction == PsiReconstruction::ProjectedPhiH1)
  {
    TrialFunction psiTrial(psiFes);
    TestFunction  psiTest(psiFes);
    RealFunction phiData(
        [&](const Geometry::Point& p) -> Real
        {
          if (advectOn)
            return phiHAdv.getValue(p);
          const auto& X = p.getCoordinates();
          return levelSet.phi(vec2(X(0), X(1)));
        });
    auto gradPhiHAdvData = Grad(phiHAdv);
    AnalyticVectorFunction gradPhiData(
        [&](const Geometry::Point& p) -> Math::SpatialVector<Real>
        {
          if (advectOn)
            return gradPhiHAdvData.getValue(p);
          const auto& X = p.getCoordinates();
          return levelSet.grad(vec2(X(0), X(1)));
        },
        /*dimension=*/2);
    Problem psiProblem(psiTrial, psiTest);
    if (psiReconstruction == PsiReconstruction::ProjectedPhiH1)
    {
      RealFunction<Real> ell2(psiProjectEll * psiProjectEll);
      psiProblem =
          Integral(psiTrial, psiTest)
        + Integral(ell2 * Grad(psiTrial), Grad(psiTest))
        - Integral(phiData, psiTest)
        - Integral(ell2 * gradPhiData, Grad(psiTest))
        + DirichletBC(psiTrial, RealFunction(0.0))
            .on(interfaceAttribute);
    }
    else
    {
      psiProblem =
          Integral(psiTrial, psiTest)
        - Integral(phiData, psiTest)
        + DirichletBC(psiTrial, RealFunction(0.0))
            .on(interfaceAttribute);
    }
    Solver::SparseLU(psiProblem).solve();
    psi.getData() = psiTrial.getSolution().getData();
  }
  else
  {
    TrialFunction psiTrial(psiFes);
    TestFunction  psiTest(psiFes);
    // Sign source: psi is the screened-Poisson lift of an inside/outside
    // step function of magnitude `psiSourceMagnitude`, pinned to zero on
    // the classifier-cut skeleton (interfaceAttribute). This makes psi a
    // smooth topological reference that is independent of the geometric
    // reference phi.
    RealFunction psiSource(
        [&](const Geometry::Point& p) -> Real
        {
          const auto attr = p.getPolytope().getAttribute();
          if (attr && *attr == interiorAttribute)
            return -psiSourceMagnitude;
          return psiSourceMagnitude;
        });

    Problem psiProblem(psiTrial, psiTest);
    psiProblem =
        Integral((psiEll * psiEll) * Grad(psiTrial), Grad(psiTest))
      + Integral(psiTrial, psiTest)
      - Integral(psiSource, psiTest)
      + DirichletBC(psiTrial, RealFunction(0.0))
          .on(interfaceAttribute);
    Solver::SparseLU(psiProblem).solve();
    psi.getData() = psiTrial.getSolution().getData();

    // Normal-gradient calibration. The screened-Poisson reconstruction fixes
    // the interface and topology, but its amplitude is set by the chosen
    // source and screening length. A single scalar is therefore applied so
    // that |grad psi| matches |grad phi| in a narrow band around the target
    // interface. The sign is chosen by the weighted value correlation.
    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const auto& cell = *cellIt;
      const auto& vertices = cell.getVertices();
      if (vertices.size() != 3)
        throw std::runtime_error(
            "LevelSetLSRReconstruction expects triangular cells.");

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

      const Vec2 gradPsiCell = cellP1Gradient(x, psiNode);
      const Real gradPsiNorm = gradPsiCell.norm();
      const Vec2 e1 = x[1] - x[0];
      const Vec2 e2 = x[2] - x[0];
      const Real triangleArea =
        Real(0.5) * std::abs(e1(0) * e2(1) - e1(1) * e2(0));

      for (const auto& bary : TriangleBarycentricQuadrature)
      {
        const Vec2 xq = interpolateVec(x, bary);
        const Real psiq =
          bary[0] * psiNode[0] + bary[1] * psiNode[1]
        + bary[2] * psiNode[2];
        const Real phiq = levelSet.phi(xq);
        const Real weight =
          std::exp(-phiq * phiq / (Real(2) * deltaW * deltaW));
        const Real qWeight =
          triangleArea / static_cast<Real>(TriangleBarycentricQuadrature.size());
        psiScaleNumerator += weight * levelSet.grad(xq).norm() * qWeight;
        psiScaleDenominator += weight * gradPsiNorm * qWeight;
        psiScaleSignMoment += weight * psiq * phiq * qWeight;
      }
    }

    if (psiScaleDenominator > Real(1e-30))
    {
      psiScale = psiScaleNumerator / psiScaleDenominator;
      if (psiScaleSignMoment < 0)
        psiScale = -psiScale;
      psi.getData() *= psiScale;
    }
  }
  std::cout << "Phase psi reconstruction/calibration: "
            << elapsedSeconds(phaseStart) << " s\n";
  phaseStart = Clock::now();

  GridFunction phiScreened(psiFes);
  phiScreened.setName("phi_screened");
  Real phiScreenScale = 1;
  if (phiValueSource == PhiValueSource::Screened)
  {
    const auto phiMoments = collectCellMomentInfo(mesh, levelSet, epsilon);
    std::unordered_map<Index, std::size_t> phiCellToLocal;
    std::vector<Index> phiLocalToCell;
    phiCellToLocal.reserve(phiMoments.size());
    phiLocalToCell.reserve(phiMoments.size());
    for (std::size_t local = 0; local < phiMoments.size(); ++local)
    {
      phiCellToLocal[phiMoments[local].index] = local;
      phiLocalToCell.push_back(phiMoments[local].index);
    }

    std::vector<Real> phiVolumes(phiMoments.size());
    std::vector<Real> phiMomentValues(phiMoments.size());
    for (std::size_t local = 0; local < phiMoments.size(); ++local)
    {
      phiVolumes[local] = phiMoments[local].area;
      phiMomentValues[local] = phiMoments[local].moment;
    }

    std::vector<MinSTCut::Edge> phiGraphEdges;
    for (auto faceIt = mesh.getFace(); faceIt; ++faceIt)
    {
      const Index facet = faceIt->getIndex();
      const auto& incident = mesh.getConnectivity().getIncidence({1, 2}, facet);
      if (incident.size() == 2)
      {
        const auto itA = phiCellToLocal.find(incident[0]);
        const auto itB = phiCellToLocal.find(incident[1]);
        if (itA == phiCellToLocal.end() || itB == phiCellToLocal.end())
          continue;
        phiGraphEdges.push_back({
            static_cast<Index>(itA->second),
            static_cast<Index>(itB->second),
            lambdaC * facetLength(mesh, facet),
            facet});
      }
    }

    const MinSTCut::Result phiClassified =
      cut.classify(phiVolumes, phiMomentValues, phiGraphEdges);
    std::vector<Index> phiInterfaceFacets;
    phiInterfaceFacets.reserve(phiClassified.cutEdges.size());
    for (const MinSTCut::Edge& edge : phiClassified.cutEdges)
      if (edge.index != MinSTCut::InvalidIndex)
        phiInterfaceFacets.push_back(edge.index);

    applyClassificationAttributes(
        phiClassified.labels, phiLocalToCell, phiInterfaceFacets);

    TrialFunction phiTrial(psiFes);
    TestFunction  phiTest(psiFes);
    RealFunction phiSource(
        [&](const Geometry::Point& p) -> Real
        {
          const auto attr = p.getPolytope().getAttribute();
          if (attr && *attr == interiorAttribute)
            return -psiSourceMagnitude;
          return psiSourceMagnitude;
        });
    Problem phiProblem(phiTrial, phiTest);
    phiProblem =
        Integral((psiEll * psiEll) * Grad(phiTrial), Grad(phiTest))
      + Integral(phiTrial, phiTest)
      - Integral(phiSource, phiTest)
      + DirichletBC(phiTrial, RealFunction(0.0))
          .on(interfaceAttribute);
    Solver::SparseLU(phiProblem).solve();
    phiScreened.getData() = phiTrial.getSolution().getData();

    Real phiScaleNumerator = 0;
    Real phiScaleDenominator = 0;
    Real phiScaleSignMoment = 0;
    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const auto& cell = *cellIt;
      const auto& vertices = cell.getVertices();
      std::array<Vec2, 3> x;
      std::array<Real, 3> phiNode;
      for (std::size_t a = 0; a < 3; ++a)
      {
        x[a] = mesh.getVertexCoordinates(vertices[a]);
        Math::SpatialPoint rc(2);
        rc.setZero();
        if (a == 1) rc(0) = 1;
        if (a == 2) rc(1) = 1;
        phiNode[a] = phiScreened.getValue(Geometry::Point(cell, rc));
      }

      const Vec2 gradPhiCell = cellP1Gradient(x, phiNode);
      const Real gradPhiNorm = gradPhiCell.norm();
      const Vec2 e1 = x[1] - x[0];
      const Vec2 e2 = x[2] - x[0];
      const Real triangleArea =
        Real(0.5) * std::abs(e1(0) * e2(1) - e1(1) * e2(0));

      for (const auto& bary : TriangleBarycentricQuadrature)
      {
        const Vec2 xq = interpolateVec(x, bary);
        const Real phiScreenq =
          bary[0] * phiNode[0] + bary[1] * phiNode[1]
        + bary[2] * phiNode[2];
        const Real phiq = levelSet.phi(xq);
        const Real weight =
          std::exp(-phiq * phiq / (Real(2) * deltaW * deltaW));
        const Real qWeight =
          triangleArea / static_cast<Real>(TriangleBarycentricQuadrature.size());
        phiScaleNumerator += weight * levelSet.grad(xq).norm() * qWeight;
        phiScaleDenominator += weight * gradPhiNorm * qWeight;
        phiScaleSignMoment += weight * phiScreenq * phiq * qWeight;
      }
    }
    if (phiScaleDenominator > Real(1e-30))
    {
      phiScreenScale = phiScaleNumerator / phiScaleDenominator;
      if (phiScaleSignMoment < 0)
        phiScreenScale = -phiScreenScale;
      phiScreened.getData() *= phiScreenScale;
    }

    applyClassificationAttributes(
        classified.labels, localToCell, interfaceFacets);
  }

  // -------------------------------------------------------------------------
  // phi is analytic by default; with --phi=screened it is reconstructed by
  // the same screened-Poisson profile used for psi.
  // Adaptation::LSR evaluates these operands at deformed/traced points.
  // -------------------------------------------------------------------------


  // -------------------------------------------------------------------------
  // Step 5a: smooth band weight and weighted band measure M_w.
  // -------------------------------------------------------------------------
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
    throw std::runtime_error("Empty smooth band: M_w = 0.");
  std::cout << "Phase band measure: "
            << elapsedSeconds(phaseStart) << " s\n";
  phaseStart = Clock::now();

  // -------------------------------------------------------------------------
  // Step 5b: Vector FES (P1 by default, H1<2>/P2 when
  // RODIN_LSR_P2_DISPLACEMENT is defined), trial/test/state functions.
  // -------------------------------------------------------------------------
#ifdef RODIN_LSR_P2_DISPLACEMENT
  using VectorFES = H1<2, Math::SpatialVector<Real>, LocalMesh>;
  VectorFES vectorFes(std::integral_constant<std::size_t, 2>{}, mesh, 2);
#else
  using VectorFES = P1<Math::SpatialVector<Real>, LocalMesh>;
  VectorFES vectorFes(mesh, 2);
#endif
  TrialFunction du(vectorFes);
  TestFunction  v(vectorFes);
  GridFunction  u(vectorFes);
  u.getData().setZero();

  // LSRIntegratorParameters no longer carry a tangent-mode flag — that's a type
  // property of the LSRRegistration composite chosen below.
  LSRIntegratorParameters params;
  params.rhoS = 1;
  params.deltaW = deltaW;
  params.hRef = h;
  params.normalizer = Real(1) / (weightedBandMeasure * h * h);

  // gamma is now a form-language scalar field. Constant weights are
  // wrapped in `RealFunction<Real>(value)` exactly like any other
  // scalar constant in Rodin form language; spatially varying weights
  // can be any `Variational::RealFunctionBase<...>`. gamma weighs the
  // shape term E_shape.
  RealFunction<Real> barrierGamma(shapeWeight);
  BarrierParameters barrierParams;
  barrierParams.jMin  = 1e-8;
  barrierParams.domainMeasure = domainMeasure;

  // Per-cell geometry cache, indexed by local iteration order.
  auto [cellCache, geomToLocal] = precomputeCellGeometry(mesh);
  // The cellToLocal map built from cellMoments matches the iteration
  // order used by precomputeCellGeometry, but assert it for safety.
  for (const auto& [parent, local] : cellToLocal)
  {
    const auto it = geomToLocal.find(parent);
    if (it == geomToLocal.end() || it->second != local)
      throw std::runtime_error(
          "cellToLocal map disagrees with precomputeCellGeometry iteration "
          "order — code paths got out of sync.");
  }
  std::cout << "Phase FES/cache setup: "
            << elapsedSeconds(phaseStart) << " s\n";
  phaseStart = Clock::now();

  // -------------------------------------------------------------------------
  // The sampled relative-distortion barrier is verified separately by the
  // unit test
  //   tests/unit/Rodin/Adaptation/BarrierSampledResidualFDTest.cpp
  // We do not repeat that check at runtime here.
  //
  // Step 6: Problem assembly and Newton solve.
  // -------------------------------------------------------------------------
  auto zero = VectorFunction{ Zero(), Zero() };
  GridFunction gradPhiScreenedSmoothed(gradPhiFes);
  gradPhiScreenedSmoothed.getData().setZero();
  if (phiValueSource == PhiValueSource::Screened)
  {
    TrialFunction gtrial(gradPhiFes);
    TestFunction  gtest(gradPhiFes);
    Problem gpb(gtrial, gtest);
    gpb = Integral(gtrial, gtest)
        - Integral(Grad(phiScreened), gtest);
    Solver::SparseLU(gpb).solve();
    gradPhiScreenedSmoothed.getData() = gtrial.getSolution().getData();
  }

  // Data phi for the LSR push solve.
  //   advect=off : analytic levelSet everywhere (ground truth).
  //   advect=on  : phi value is the once-advected discrete P1 field
  //                phi_h_adv. Gradient / Hess come from one of three
  //                sources selected by --phi-grad:
  //                  discrete : Grad(phi_h_adv), hess = 0.
  //                  analytic : levelSet.grad/.hess at the advected
  //                             position (no discretisation of derivatives).
  //                  smoothed : L2-projection of Grad(phi_h_adv) onto P1
  //                             [continuous nodal recovery], hess via
  //                             Jacobian of the smoothed grad.
  RealFunction phi(
      [&](const Geometry::Point& p) -> Real
      {
        if (phiValueSource == PhiValueSource::Screened)
          return phiScreened.getValue(p);
        if (advectOn)
          return phiHAdv.getValue(p);
        const auto& c = p.getPhysicalCoordinates();
        return levelSet.phi(vec2(c(0), c(1)));
      });
  auto gradPhiDiscrete = Grad(phiHAdv);
  auto hessFromSmoothedGrad = Jacobian(gradPhiSmoothed);
  auto hessFromScreenedGrad = Jacobian(gradPhiScreenedSmoothed);
  AnalyticVectorFunction gradPhi(
      [&](const Geometry::Point& p) -> Math::SpatialVector<Real>
      {
        if (phiValueSource == PhiValueSource::Screened)
          return gradPhiScreenedSmoothed.getValue(p);
        if (advectOn)
        {
          switch (phiGradSource)
          {
            case PhiGradSource::Discrete:
              return gradPhiDiscrete.getValue(p);
            case PhiGradSource::Smoothed:
              return gradPhiSmoothed.getValue(p);
            case PhiGradSource::Analytic:
              break;
          }
        }
        const auto& c = p.getPhysicalCoordinates();
        return levelSet.grad(vec2(c(0), c(1)));
      },
      /*dimension=*/2);
  AnalyticMatrixFunction hessPhiRaw(
      [&](const Geometry::Point& p) -> Math::SpatialMatrix<Real>
      {
        if (phiValueSource == PhiValueSource::Screened)
          return hessFromScreenedGrad.getValue(p);
        if (advectOn)
        {
          switch (phiGradSource)
          {
            case PhiGradSource::Discrete:
            {
              Math::SpatialMatrix<Real> Z(2, 2);
              Z.setZero();
              return Z;
            }
            case PhiGradSource::Smoothed:
              return hessFromSmoothedGrad.getValue(p);
            case PhiGradSource::Analytic:
              break;
          }
        }
        const auto& c = p.getPhysicalCoordinates();
        return levelSet.hess(vec2(c(0), c(1)));
      },
      /*rows=*/2, /*cols=*/2);

  GridFunction gradPsiSmoothed(gradPhiFes);
  gradPsiSmoothed.getData().setZero();
  {
    TrialFunction gtrial(gradPhiFes);
    TestFunction  gtest(gradPhiFes);
    Problem gpb(gtrial, gtest);
    gpb = Integral(gtrial, gtest)
        - Integral(Grad(psi), gtest);
    Solver::SparseLU(gpb).solve();
    gradPsiSmoothed.getData() = gtrial.getSolution().getData();
  }
  auto gradPsiRaw = Grad(psi);
  auto hessPsiFromGrad = Jacobian(gradPsiSmoothed);
  AnalyticVectorFunction gradPsi(
      [&](const Geometry::Point& p) -> Math::SpatialVector<Real>
      {
        return pullPsiDerivative == PullPsiDerivative::Smoothed
          ? gradPsiSmoothed.getValue(p)
          : gradPsiRaw.getValue(p);
      },
      /*dimension=*/2);
  AnalyticMatrixFunction hessPsiRaw(
      [&](const Geometry::Point& p) -> Math::SpatialMatrix<Real>
      {
        if (pullPsiDerivative == PullPsiDerivative::Smoothed)
          return hessPsiFromGrad.getValue(p);
        Math::SpatialMatrix<Real> Z(2, 2);
        Z.setZero();
        return Z;
      },
      /*rows=*/2, /*cols=*/2);

  GridFunction hessPhi00(psiFes), hessPhi01(psiFes);
  GridFunction hessPhi10(psiFes), hessPhi11(psiFes);
  GridFunction hessPsi00(psiFes), hessPsi01(psiFes);
  GridFunction hessPsi10(psiFes), hessPsi11(psiFes);
  auto projectHessianComponent =
    [&](const auto& hess, std::size_t row, std::size_t col, auto& out)
    {
      TrialFunction htrial(psiFes);
      TestFunction  htest(psiFes);
      RealFunction component(
          [&](const Geometry::Point& p) -> Real
          {
            return hess.getValue(p)(row, col);
          });
      Problem hpb(htrial, htest);
      if (hessSmoothing == HessianSmoothing::H1)
      {
        RealFunction<Real> ell2(hessSmoothEll * hessSmoothEll);
        hpb = Integral(htrial, htest)
            + Integral(ell2 * Grad(htrial), Grad(htest))
            - Integral(component, htest);
      }
      else
      {
        hpb = Integral(htrial, htest)
            - Integral(component, htest);
      }
      Solver::SparseLU(hpb).solve();
      out.getData() = htrial.getSolution().getData();
    };
  if (hessSmoothing != HessianSmoothing::None)
  {
    projectHessianComponent(hessPhiRaw, 0, 0, hessPhi00);
    projectHessianComponent(hessPhiRaw, 0, 1, hessPhi01);
    projectHessianComponent(hessPhiRaw, 1, 0, hessPhi10);
    projectHessianComponent(hessPhiRaw, 1, 1, hessPhi11);
    projectHessianComponent(hessPsiRaw, 0, 0, hessPsi00);
    projectHessianComponent(hessPsiRaw, 0, 1, hessPsi01);
    projectHessianComponent(hessPsiRaw, 1, 0, hessPsi10);
    projectHessianComponent(hessPsiRaw, 1, 1, hessPsi11);
  }
  AnalyticMatrixFunction hessPhi(
      [&](const Geometry::Point& p) -> Math::SpatialMatrix<Real>
      {
        if (hessSmoothing == HessianSmoothing::None)
          return hessPhiRaw.getValue(p);
        Math::SpatialMatrix<Real> H(2, 2);
        H(0, 0) = hessPhi00.getValue(p);
        H(0, 1) = hessPhi01.getValue(p);
        H(1, 0) = hessPhi10.getValue(p);
        H(1, 1) = hessPhi11.getValue(p);
        return H;
      },
      /*rows=*/2, /*cols=*/2);
  AnalyticMatrixFunction hessPsi(
      [&](const Geometry::Point& p) -> Math::SpatialMatrix<Real>
      {
        if (hessSmoothing == HessianSmoothing::None)
          return hessPsiRaw.getValue(p);
        Math::SpatialMatrix<Real> H(2, 2);
        H(0, 0) = hessPsi00.getValue(p);
        H(0, 1) = hessPsi01.getValue(p);
        H(1, 0) = hessPsi10.getValue(p);
        H(1, 1) = hessPsi11.getValue(p);
        return H;
      },
      /*rows=*/2, /*cols=*/2);
  RealFunction psiFunction(
      [&](const Geometry::Point& p) -> Real
      {
        return psi.getValue(p);
      });

  bool solveCompleted = true;
  bool newtonConverged = false;
  std::size_t newtonIterations = 0;
  Real initialGuessAlpha = 0;
  Real initialGuessMinJ = 0;
  std::size_t initialGuessBacktracks = 0;
  std::string solveError;
  constexpr Real geometryFitTolMult = Real(4);
  const Real geometryTolerance = geometryFitTolMult * h * h;
  auto interfaceFacetSamples =
    [&](const Index facet, const auto& callback)
    {
      const auto face = mesh.getFace(facet);
      const auto& faceVertices = face->getVertices();
      const auto& incident =
        mesh.getConnectivity().getIncidence({1, 2}, facet);
      if (faceVertices.size() != 2 || incident.empty())
        return;
      const Index cellIdx = incident[0];
      const auto cellIt = mesh.getCell(cellIdx);
      const auto& cell = *cellIt;
      const auto& cellVertices = cell.getVertices();
      std::array<int, 2> local = {{-1, -1}};
      for (std::size_t a = 0; a < cellVertices.size(); ++a)
      {
        if (cellVertices[a] == faceVertices[0])
          local[0] = static_cast<int>(a);
        if (cellVertices[a] == faceVertices[1])
          local[1] = static_cast<int>(a);
      }
      if (local[0] < 0 || local[1] < 0)
        return;

      const Vec2 xa = mesh.getVertexCoordinates(faceVertices[0]);
      const Vec2 xb = mesh.getVertexCoordinates(faceVertices[1]);
      const Real len = (xb - xa).norm();
      constexpr std::array<Real, 2> sPts = {{
        Real(0.21132486540518713),
        Real(0.78867513459481287)
      }};
      for (const Real s : sPts)
      {
        std::array<Real, 3> lambda = {{Real(0), Real(0), Real(0)}};
        lambda[static_cast<std::size_t>(local[0])] = Real(1) - s;
        lambda[static_cast<std::size_t>(local[1])] = s;
        Math::SpatialPoint rc(2);
        rc(0) = lambda[1];
        rc(1) = lambda[2];
        Math::SpatialPoint X(2);
        cell.getTransformation().transform(X, rc);
        callback(Geometry::Point(cell, rc, X), len * Real(0.5));
      }
    };

  auto interfacePhiRMSFromDisplacement = [&]() -> Real
  {
    Real interfacePhi = 0;
    Real interfaceLength = 0;
    for (const Index facet : interfaceFacets)
    {
      interfaceFacetSamples(
          facet,
          [&](const Geometry::Point& pt, Real weight)
          {
            const auto& X = pt.getPhysicalCoordinates();
            const auto uX = u.getValue(pt);
            const Vec2 y = vec2(X(0) + uX(0), X(1) + uX(1));
            const Real val = levelSet.phi(y);
            interfacePhi += val * val * weight;
            interfaceLength += weight;
          });
    }
    return std::sqrt(interfacePhi / std::max(interfaceLength, Real(1e-30)));
  };

  // Discrete fit: phi_h evaluated at the LSR-displaced interface
  // midpoints. With advect=on the LSR data residual is
  //   (phi_h(X+u) - psi(X))^2,
  // and at convergence X+u sits on phi_h's zero contour. This metric
  // reports how well that was achieved, in contrast with
  // `interfacePhiRMSFromDisplacement` above which uses analytic phi as
  // ground truth (sensitive to phi_h vs analytic gap AND to classifier
  // chord error). For --advect=off the two metrics coincide.
  auto interfacePhiHRMSFromDisplacement = [&]() -> Real
  {
    if (!advectOn) return interfacePhiRMSFromDisplacement();
    Real interfacePhi = 0;
    Real interfaceLength = 0;
    for (const Index facet : interfaceFacets)
    {
      interfaceFacetSamples(
          facet,
          [&](const Geometry::Point& pt, Real weight)
          {
            const auto& X = pt.getPhysicalCoordinates();
            const auto uX = u.getValue(pt);
            const Vec2 y = vec2(X(0) + uX(0), X(1) + uX(1));
            Real val = 0;
            for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
            {
              const auto& cell = *cellIt;
              Math::SpatialPoint Y(2);
              Y(0) = y(0); Y(1) = y(1);
              Math::SpatialPoint rc(2);
              cell.getTransformation().inverse(rc, Y);
              const Real tol = Real(1e-9);
              if (rc(0) >= -tol && rc(1) >= -tol
                  && rc(0) + rc(1) <= Real(1) + tol)
              {
                val = phiHAdv.getValue(Geometry::Point(cell, rc, Y));
                break;
              }
            }
            interfacePhi += val * val * weight;
            interfaceLength += weight;
          });
    }
    return std::sqrt(interfacePhi / std::max(interfaceLength, Real(1e-30)));
  };

  LSRTangent effectiveTangent = requestedTangent;
  {
    LSRParameters lsrParams;
    lsrParams.rhoS = params.rhoS;
    lsrParams.deltaW = params.deltaW;
    lsrParams.hRef = params.hRef;
    lsrParams.normalizer = params.normalizer;
    lsrParams.shapeWeight = shapeWeight;
    lsrParams.interfaceWeight = gammaWeight;
    lsrParams.interfaceAttribute = interfaceAttribute;
    lsrParams.jBarrierWeight = jBarrierWeight;
    lsrParams.jBarrierSafeRatio = jBarrierSafeRatio;
    lsrParams.jVolumeTetherWeight = jVolumeTetherWeight;
    lsrParams.useSampledQuadraticAlphaPredictor =
      !hasFlag(argc, argv, "no-alpha-predictor");
    if (hasFlag(argc, argv, "flow-trace"))
      lsrParams.fieldEvaluation =
        LSRIntegratorParameters::FieldEvaluation::FlowTrace;
    lsrParams.jMinRatio = barrierParams.jMin;
    lsrParams.jSafeRatio = 1e-3;
    lsrParams.initialGuess = requestedInitialGuess;
    effectiveTangent =
        (dataEnergy != DataEnergy::Pull
         && advectOn
         && phiGradSource == PhiGradSource::Discrete
         && requestedTangent != LSRTangent::GaussNewton)
          ? LSRTangent::GaussNewton
          : requestedTangent;
    lsrParams.tangent = effectiveTangent;
    lsrParams.maxNewtonIterations = maxNewtonSteps;
    lsrParams.absoluteTolerance = 1e-8;
    lsrParams.relativeTolerance = 1e-7;
    lsrParams.alphaWarmStartGrowth = 4;
    std::size_t lastPrintedBacktracks = 0;
    lsrParams.acceptedStateConvergenceTest =
      [&](const LSRReport& report) -> bool
      {
        const Real fit = interfacePhiRMSFromDisplacement();
        const std::size_t stepBacktracks =
          report.totalBacktracks - lastPrintedBacktracks;
        lastPrintedBacktracks = report.totalBacktracks;
        std::cout << "Newton " << std::setw(2) << report.iterations
                  << ": ||R||=" << std::scientific << std::setprecision(6)
                  << report.finalResidual
                  << "  step=" << report.finalStepNorm
                  << "  alpha=" << report.lastAcceptedAlpha
                  << "  backtracks=" << stepBacktracks
                  << "  min_j=" << report.minJRatio
                  << "  fit=" << fit
                  << '\n';
        return fit <= geometryTolerance;
      };

    // Pull-back adapters: psi(X-u) and grad psi(X-u) via Flow with constant
    // velocity = -u(X) at each starting integration point. With the
    // constant velocity, Phi_1(X) = X + (-u(X)) = X - u(X) exactly, so the
    // pulled-back Lagrangian semantics is preserved. The cell-walking from
    // the original cell to the cell containing X - u is delegated to
    // Variational::Flow.
    auto recoverOriginalPoint = [](const Geometry::Point& p)
    {
      Math::SpatialPoint X;
      p.getPolytope().getTransformation().transform(
          X, p.getReferenceCoordinates());
      return Geometry::Point(
          p.getPolytope(), p.getReferenceCoordinates(), X);
    };
    RealFunction psiDisp(
        [&](const Geometry::Point& p) -> Real
        {
          const auto origPoint = recoverOriginalPoint(p);
          const auto uX = u.getValue(origPoint);
          Math::SpatialVector<Real> v(2);
          v(0) = -uX(0); v(1) = -uX(1);
          auto constV = [v](const Geometry::Point&) { return v; };
          Flow flow(Real(1), psi, constV);
          const auto tr = flow.trace(origPoint);
          if (tr.exited()) return Real(0);
          return psi.getValue(tr.getPoint()) + tr.getCorrection();
        });
    AnalyticVectorFunction gradPsiDisp(
        [&](const Geometry::Point& p) -> Math::SpatialVector<Real>
        {
          const auto origPoint = recoverOriginalPoint(p);
          const auto uX = u.getValue(origPoint);
          Math::SpatialVector<Real> v(2);
          v(0) = -uX(0); v(1) = -uX(1);
          auto constV = [v](const Geometry::Point&) { return v; };
          Flow flow(Real(1), psi, constV);
          const auto tr = flow.trace(origPoint);
          if (tr.exited())
          {
            Math::SpatialVector<Real> z(2);
            z.setZero();
            return z;
          }
          return gradPsi.getValue(tr.getPoint());
        },
        /*dimension=*/2);
    AnalyticMatrixFunction hessPsiDisp(
        [&](const Geometry::Point& p) -> Math::SpatialMatrix<Real>
        {
          const auto origPoint = recoverOriginalPoint(p);
          const auto uX = u.getValue(origPoint);
          Math::SpatialVector<Real> v(2);
          v(0) = -uX(0); v(1) = -uX(1);
          auto constV = [v](const Geometry::Point&) { return v; };
          Flow flow(Real(1), psi, constV);
          const auto tr = flow.trace(origPoint);
          if (tr.exited())
          {
            Math::SpatialMatrix<Real> Z(2, 2);
            Z.setZero();
            return Z;
          }
          return hessPsi.getValue(tr.getPoint());
        },
        /*rows=*/2, /*cols=*/2);

    try
    {
      std::cout << "Phase LSR solve: begin (energy=" << energyStr << ")\n";
      LSR lsr(u);
      LSRReport lsrReport;
      if (dataEnergy == DataEnergy::Pull)
      {
        lsrReport = lsr.setParameters(lsrParams).solvePull(
            psi, phi, gradPhi, psiDisp, gradPsiDisp, hessPsiDisp);
      }
      else if (dataEnergy == DataEnergy::Push)
      {
        lsrReport = lsr.setParameters(lsrParams).solve(
            psi, phi, gradPhi, hessPhi);
      }
      else
      {
        LSRParameters swappedParams = lsrParams;
        swappedParams.normalizer = 0;
        swappedParams.acceptedStateConvergenceTest = nullptr;
        lsrReport = lsr.setParameters(swappedParams).solve(
            phi, psiFunction, gradPsi, hessPsi);

        u.getData() *= Real(-1);
        std::cout << "Phase LSR solve: swapped push produced inverse "
                  << "candidate; using u=-v";
        if (dataEnergy == DataEnergy::PushSwappedForward)
        {
          std::cout << " and refining with forward push\n";
          LSRParameters forwardParams = lsrParams;
          forwardParams.initialGuess = LSRInitialGuess::Current;
          lsrReport = lsr.setParameters(forwardParams).solve(
              psi, phi, gradPhi, hessPhi);
        }
        else
        {
          std::cout << '\n';
        }
      }
      newtonConverged = lsrReport.converged;
      newtonIterations = lsrReport.iterations;
      initialGuessAlpha = lsrReport.initialGuessAlpha;
      initialGuessMinJ = lsrReport.initialGuessMinJRatio;
      initialGuessBacktracks = lsrReport.initialGuessBacktracks;
      solveCompleted = lsrReport.converged && !lsrReport.lineSearchFailed;
      if (!solveCompleted)
        solveError = "LSR facade did not converge";
      std::cout << "Phase LSR solve: "
                << elapsedSeconds(phaseStart) << " s\n";
      phaseStart = Clock::now();
    }
    catch (const std::exception& ex)
    {
      solveCompleted = false;
      solveError = ex.what();
      std::cout << "\nLSR facade aborted: " << ex.what() << '\n';
    }
  }

  // -------------------------------------------------------------------------
  // Differential-coordinate diagnostics.
  //
  // These quantities compare the two normal fields which define the
  // Gauss-Newton part of the push and pull linearizations:
  //   push : grad phi(X + u(X)),
  //   pull : grad psi(X - u(X)).
  // Equality of zero sets does not imply equality of these fields; the
  // reported angle and scale defects are direct diagnostics of the numerical
  // asymmetry between E_push and E_pull.
  // -------------------------------------------------------------------------
  struct DifferentialCoordinateStats
  {
    Real measure = 0;
    Real angleMean = 0;
    Real angleMax = 0;
    Real cosMean = 0;
    Real gradPhiMean = 0;
    Real gradPsiMean = 0;
    Real gradRatioMean = 0;
    Real gradRatioMin = std::numeric_limits<Real>::infinity();
    Real gradRatioMax = 0;
    Real hessPhiMean = 0;
    Real hessPsiMean = 0;
    Real hessRatioMean = 0;
    Real pushResidualRms = 0;
    Real pullResidualRms = 0;
    std::size_t samples = 0;
    std::size_t skipped = 0;
  } dc;

  auto withLocalizedPoint =
    [&](const Vec2& y, const auto& callback) -> bool
    {
      Math::SpatialPoint Y(2);
      Y(0) = y(0);
      Y(1) = y(1);
      for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
      {
        const auto& cell = *cellIt;
        Math::SpatialPoint rc(2);
        cell.getTransformation().inverse(rc, Y);
        const Real tol = Real(1e-9);
        if (rc(0) >= -tol && rc(1) >= -tol
            && rc(0) + rc(1) <= Real(1) + tol)
        {
          callback(Geometry::Point(cell, rc, Y));
          return true;
        }
      }
      return false;
    };

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
      const Real qWeight =
        triangleArea / static_cast<Real>(TriangleBarycentricQuadrature.size());
      const Real s = psi.getValue(pt);
      const Real W = std::exp(-s * s / (2 * deltaW * deltaW));
      const Real w = W * qWeight;
      if (w <= Real(0))
        continue;

      const auto& Xc = pt.getPhysicalCoordinates();
      const auto uX = u.getValue(pt);
      const Vec2 X = vec2(Xc(0), Xc(1));
      const Vec2 yPush = vec2(X(0) + uX(0), X(1) + uX(1));
      const Vec2 yPull = vec2(X(0) - uX(0), X(1) - uX(1));

      bool havePhi = false;
      bool havePsi = false;
      Real phiPush = 0;
      Real psiPull = 0;
      Vec2 gPhi(2), gPsi(2);
      Math::SpatialMatrix<Real> HPhi(2, 2), HPsi(2, 2);
      gPhi.setZero();
      gPsi.setZero();
      HPhi.setZero();
      HPsi.setZero();

      havePhi = withLocalizedPoint(
          yPush,
          [&](const Geometry::Point& yp)
          {
            phiPush = phi.getValue(yp);
            gPhi = gradPhi.getValue(yp);
            HPhi = hessPhi.getValue(yp);
          });
      havePsi = withLocalizedPoint(
          yPull,
          [&](const Geometry::Point& ym)
          {
            psiPull = psi.getValue(ym);
            gPsi = gradPsi.getValue(ym);
            HPsi = hessPsi.getValue(ym);
          });
      if (!havePhi || !havePsi)
      {
        dc.skipped++;
        continue;
      }

      const Real psiX = psi.getValue(pt);
      const Real phiX = phi.getValue(pt);
      const Real gPhiNorm = gPhi.norm();
      const Real gPsiNorm = gPsi.norm();
      const Real hPhiNorm = std::sqrt(HPhi.squaredNorm());
      const Real hPsiNorm = std::sqrt(HPsi.squaredNorm());
      if (gPhiNorm <= Real(1e-14) || gPsiNorm <= Real(1e-14))
      {
        dc.skipped++;
        continue;
      }

      const Real cosTheta =
        std::clamp(gPhi.dot(gPsi) / (gPhiNorm * gPsiNorm), Real(-1), Real(1));
      const Real angle = std::acos(cosTheta);
      const Real gradRatio = gPhiNorm / gPsiNorm;
      const Real hessRatio =
        hPsiNorm > Real(1e-14) ? hPhiNorm / hPsiNorm : Real(0);
      const Real rPush = phiPush - psiX;
      const Real rPull = phiX - psiPull;

      dc.measure += w;
      dc.angleMean += w * angle;
      dc.angleMax = std::max(dc.angleMax, angle);
      dc.cosMean += w * cosTheta;
      dc.gradPhiMean += w * gPhiNorm;
      dc.gradPsiMean += w * gPsiNorm;
      dc.gradRatioMean += w * gradRatio;
      dc.gradRatioMin = std::min(dc.gradRatioMin, gradRatio);
      dc.gradRatioMax = std::max(dc.gradRatioMax, gradRatio);
      dc.hessPhiMean += w * hPhiNorm;
      dc.hessPsiMean += w * hPsiNorm;
      dc.hessRatioMean += w * hessRatio;
      dc.pushResidualRms += w * rPush * rPush;
      dc.pullResidualRms += w * rPull * rPull;
      dc.samples++;
    }
  }
  if (dc.measure > Real(0))
  {
    dc.angleMean /= dc.measure;
    dc.cosMean /= dc.measure;
    dc.gradPhiMean /= dc.measure;
    dc.gradPsiMean /= dc.measure;
    dc.gradRatioMean /= dc.measure;
    dc.hessPhiMean /= dc.measure;
    dc.hessPsiMean /= dc.measure;
    dc.hessRatioMean /= dc.measure;
    dc.pushResidualRms = std::sqrt(dc.pushResidualRms / dc.measure);
    dc.pullResidualRms = std::sqrt(dc.pullResidualRms / dc.measure);
  }

  // -------------------------------------------------------------------------
  // Step 7: build moved mesh through FES dof mapping (no raw indexing).
  // -------------------------------------------------------------------------
  LocalMesh moved(mesh);
  updateMovedMeshFromDisplacement(mesh, moved, u);
  moved.save("LevelSetLSRReconstruction_phi.mesh", IO::FileFormat::MFEM);

  const Real interfacePhiRMS = interfacePhiRMSFromDisplacement();

  // -------------------------------------------------------------------------
  // Step 8: XDMF output. P0 cell writes go through FES.getGlobalIndex; the
  // cell-index <-> local-index translation goes through the explicit map.
  // -------------------------------------------------------------------------
  using ScalarP0 = P0<Real, LocalMesh>;
  ScalarP0 p0Fes(mesh);
  GridFunction cellLabel(p0Fes);
  GridFunction phaseMoment(p0Fes);
  GridFunction sigmaKgf(p0Fes);
  {
    const std::size_t d = mesh.getDimension();
    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const Index cellIdx = cellIt->getIndex();
      const std::size_t local = cellToLocal.at(cellIdx);
      const Index dof = p0Fes.getGlobalIndex({d, cellIdx}, 0);
      cellLabel.getData()(dof)   = static_cast<Real>(classified.labels[local]);
      phaseMoment.getData()(dof) = cellMoments[local].moment;
      sigmaKgf.getData()(dof)    =
        static_cast<Real>(cellCache[local].sigmaK);
    }
  }
  cellLabel.setName("cell_label");
  phaseMoment.setName("phase_moment");
  sigmaKgf.setName("sigma_K");

  P1<Real, LocalMesh> p1Fes(mesh);
  GridFunction phiGf(p1Fes);
  phiGf = [&](const Geometry::Point& p) -> Real
  {
    const auto& X = p.getCoordinates();
    return levelSet.phi(vec2(X(0), X(1)));
  };
  phiGf.setName("phi");
  psi.setName("psi");

  // ----- phi (moved) mesh -----
  //
  // We report two distinct shape quality diagnostics:
  //
  //   q_abs := Q_abs(A_K^u)     -- absolute intrinsic shape quality of
  //                                the moved cell. Equals 1 iff the moved
  //                                cell is a similarity of the REFERENCE
  //                                cell (right triangle in Rodin).
  //
  //   q_rel := Q_rel(F),  F = A_K^u (A_K)^{-1} = I + grad u_h
  //                             -- relative-distortion measure that the
  //                                shape barrier in LSR actually minimises.
  //                                Equals 1 at u = 0 for every cell
  //                                regardless of background shape; grows
  //                                only with deformation introduced by u.
  ScalarP0 p0FesMoved(moved);
  GridFunction jK(p0FesMoved);
  GridFunction qAbs(p0FesMoved);
  GridFunction qRel(p0FesMoved);
  GridFunction cellLabelPhi(p0FesMoved);
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
      jK.getData()(dof) = jKval;
      const Real dExp = Real(2) / Real(d);
      qAbs.getData()(dof) =
        dst.A.squaredNorm()
          / (static_cast<Real>(d) * std::pow(sigDetAu, dExp));
      // F = A_K^u (A_K)^{-1}, det F = j_K (>0 for admissible cells).
      const Math::SpatialMatrix<Real> F = dst.A * src.A.inverse();
      qRel.getData()(dof) =
        F.squaredNorm()
          / (static_cast<Real>(d) * std::pow(jKval, dExp));
      cellLabelPhi.getData()(dof) =
        static_cast<Real>(classified.labels[local]);
    }
  }
  jK.setName("j");
  qAbs.setName("q_abs");
  qRel.setName("q_rel");
  cellLabelPhi.setName("cell_label");

  P1<Real, LocalMesh> p1FesMoved(moved);
  GridFunction phiMoved(p1FesMoved);
  phiMoved = [&](const Geometry::Point& p) -> Real
  {
    const auto& X = p.getCoordinates();
    return levelSet.phi(vec2(X(0), X(1)));
  };
  phiMoved.setName("phi_moved");

  u.setName("displacement");

  IO::XDMF xdmf("LevelSetLSRReconstruction");
  auto psiGrid = xdmf.grid("psi");
  psiGrid.setMesh(mesh);
  psiGrid.add(cellLabel, IO::XDMF::Center::Cell);
  psiGrid.add(phaseMoment, IO::XDMF::Center::Cell);
  psiGrid.add(sigmaKgf, IO::XDMF::Center::Cell);
  psiGrid.add(phiGf, IO::XDMF::Center::Node);
  psiGrid.add(psi, IO::XDMF::Center::Node);
  psiGrid.add(u, IO::XDMF::Center::Node);

  auto phiGrid = xdmf.grid("phi");
  phiGrid.setMesh(moved);
  phiGrid.add(cellLabelPhi, IO::XDMF::Center::Cell);
  phiGrid.add(jK, IO::XDMF::Center::Cell);
  phiGrid.add(qAbs, IO::XDMF::Center::Cell);
  phiGrid.add(qRel, IO::XDMF::Center::Cell);
  phiGrid.add(phiMoved, IO::XDMF::Center::Node);

  xdmf.write().close();

  // -------------------------------------------------------------------------
  // Diagnostics.
  // -------------------------------------------------------------------------
  const char* dataEnergyName =
      dataEnergy == DataEnergy::Pull
        ? "E_pull"
    : dataEnergy == DataEnergy::Push
        ? "E_push"
    : dataEnergy == DataEnergy::PushSwapped
        ? "E_push_swapped_then_negate"
        : "E_push_swapped_then_forward_push";
  std::cout << "\nDiagnostics\n"
            << "  data energy: " << dataEnergyName << '\n'
            << "  shape energy: E_shape, weight=" << shapeWeight << '\n'
            << "  interface energy: E_Gamma, weight=" << gammaWeight
              << " (normalized phi/|grad phi|)\n"
            << "  admissibility energy: E_admissibility, weight="
              << jBarrierWeight
              << ", safe ratio=" << jBarrierSafeRatio << '\n'
            << "  volume tether: weight=" << jVolumeTetherWeight << '\n'
            << "  alpha predictor: "
              << (hasFlag(argc, argv, "no-alpha-predictor")
                    ? "sampled quadratic disabled"
                    : "sampled quadratic enabled")
              << '\n'
            << "  field evaluation: "
              << (hasFlag(argc, argv, "flow-trace")
                    ? "FlowTrace" : "PhysicalPoint")
              << '\n'
            << "  LSR tangent mode: "
              << (effectiveTangent == LSRTangent::GaussNewton
                    ? "GaussNewton"
                    : effectiveTangent == LSRTangent::Newton
                        ? "Newton"
                        : "PSDProjectedNewton")
              << '\n'
            << "  initial guess: "
              << (requestedInitialGuess == LSRInitialGuess::Zero
                    ? "zero"
                    : requestedInitialGuess == LSRInitialGuess::Current
                        ? "current"
                        : "Hilbert")
              << '\n'
            << "  pull psi derivatives: "
              << (pullPsiDerivative == PullPsiDerivative::Smoothed
                    ? "smoothed grad, recovered Hessian"
                    : "raw Grad(psi), zero Hessian")
              << '\n'
            << "  Hessian smoothing: "
              << (hessSmoothing == HessianSmoothing::None
                    ? "none"
                    : hessSmoothing == HessianSmoothing::L2
                        ? "componentwise L2 projection"
                        : "componentwise H1 projection")
              << ", ell=" << hessSmoothEll
              << '\n'
            << "  cells inside / outside: "
              << classified.insideCells.size() << " / "
              << classified.outsideCells.size() << '\n'
            << "  interface facets: " << interfaceFacets.size() << '\n'
            << "  vector dofs: " << u.getData().size() << '\n'
            << "  wavy circle: center=(" << levelSet.cx << ", "
              << levelSet.cy << "), R0=" << levelSet.R0
              << ", amp=" << levelSet.amp
              << ", lobes=" << levelSet.lobes
              << ", phase=" << levelSet.phase << '\n'
            << "  h: " << h << ", delta: " << delta
              << ", deltaW: " << deltaW << '\n'
            << "  psi reconstruction: "
              << (psiReconstruction == PsiReconstruction::Fmm
                    ? "FMM signed-distance"
                    : psiReconstruction == PsiReconstruction::ProjectedPhi
                        ? "L2 projection of phi with skeleton Dirichlet"
                        : psiReconstruction == PsiReconstruction::ProjectedPhiH1
                            ? "H1 projection of phi with skeleton Dirichlet"
                        : "screened Poisson sign source")
              << ", ell=" << psiEll
              << ", project ell=" << psiProjectEll
              << ", magnitude=" << psiSourceMagnitude
              << ", normal-gradient scale=" << psiScale << '\n'
            << "  psi scale diagnostics: numerator=" << psiScaleNumerator
              << ", denominator=" << psiScaleDenominator
              << ", sign moment=" << psiScaleSignMoment << '\n'
            << "  phi source: "
              << (phiValueSource == PhiValueSource::Screened
                    ? "screened Poisson"
                    : (advectOn ? scalarPhiFesName : "analytic WavyCircleLevelSet"))
              << ", grad source="
              << (phiValueSource == PhiValueSource::Screened
                    ? "screened-smoothed"
                    : phiGradStr)
              << ", phi scale=" << phiScreenScale << '\n'
            << "  |Omega_h|: " << domainMeasure
              << ", M_w (smooth band): " << weightedBandMeasure << '\n'
            << "  differential coordinate comparison:"
              << " samples=" << dc.samples
              << ", skipped=" << dc.skipped
              << ", band measure=" << dc.measure << '\n'
            << "    angle(phi_push, psi_pull) mean/max [deg]: "
              << (Real(180) / M_PI) * dc.angleMean
              << " / " << (Real(180) / M_PI) * dc.angleMax
              << ", mean cos=" << dc.cosMean << '\n'
            << "    |grad phi_push| mean: " << dc.gradPhiMean
              << ", |grad psi_pull| mean: " << dc.gradPsiMean
              << ", ratio mean/min/max: " << dc.gradRatioMean
              << " / " << dc.gradRatioMin
              << " / " << dc.gradRatioMax << '\n'
            << "    |H phi_push|_F mean: " << dc.hessPhiMean
              << ", |H psi_pull|_F mean: " << dc.hessPsiMean
              << ", ratio mean: " << dc.hessRatioMean << '\n'
            << "    residual RMS push/pull on same band: "
              << dc.pushResidualRms << " / " << dc.pullResidualRms << '\n'
            << "  final ||phi(X+u)||_RMS on Gamma_psi,h: "
              << interfacePhiRMS << '\n'
            << "  final ||phi_h(X+u)||_RMS on Gamma_psi,h (discrete): "
              << interfacePhiHRMSFromDisplacement() << '\n'
            << "  geometry tolerance: " << geometryTolerance << '\n'
            << "  initial guess line search: alpha=" << initialGuessAlpha
              << ", backtracks=" << initialGuessBacktracks
              << ", min_j=" << initialGuessMinJ << '\n'
            << "  Newton iterations: " << newtonIterations
              << " / max " << maxNewtonSteps
              << ", converged: " << (newtonConverged ? "yes" : "no") << '\n'
            << "  total runtime: " << elapsedSeconds(totalStart) << " s\n";

  return newtonConverged ? 0 : 1;
}
