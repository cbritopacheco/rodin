/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//
// 2D cantilever compliance minimization (level-set topology optimization)
// using WNGIR for the body-fitted mesh instead of MMG.
//
// This is the LevelSetCantilever2D loop with MMG remeshing replaced by the
// WNGIR pipeline: a fixed background grid carries the level set phi; each
// iteration the grid is classified (graph cut on phase moments) into
// interior/exterior with a cut-skeleton interface, WNGIR fits the mesh so the
// skeleton lands on phi=0, the exterior is trimmed away, and linear elasticity
// is solved on the fitted material submesh. The shape gradient (boundary strain
// energy density minus a fixed Lagrange multiplier) is Hilbert-regularized,
// the level set is reinitialized (FMM) and advected by the descent velocity.
//
//   J(Omega) = compliance(Omega) + ell * |Omega|.
//
#include <Rodin/Adaptation.h>
#include <Rodin/Advection/Lagrangian.h>
#include <Rodin/Assembly.h>
#include <Rodin/Distance/Eikonal.h>
#include <Rodin/Geometry.h>
#include <Rodin/Geometry/Region.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Math.h>
#include <Rodin/MMG.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>
#include <Rodin/Solver/CG.h>
#include <Rodin/Solver/SparseLU.h>
#include <Rodin/Solid.h>
#include <Rodin/Variational.h>

#include "../WNGIRExampleParameters.h"

#include <Eigen/IterativeLinearSolvers>

#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
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
  using WNGIRMesh = Geometry::Mesh<Context::Local>;

  // Cell / boundary attributes.
  // Numbering matches resources/.../LevelSetCantilever2D.mfem.mesh.
  constexpr Attribute interiorAttribute = 1;   // material Omega   (cell attr 1)
  constexpr Attribute exteriorAttribute = 2;   // void             (cell attr 2)
  constexpr Attribute Gamma0 = 1;               // free outer       (bdr attr 1)
  constexpr Attribute GammaD = 2;               // clamped left     (bdr attr 2)
  constexpr Attribute GammaN = 3;               // loaded right-mid (bdr attr 3)
  constexpr Attribute Gamma  = 4;               // free interface   (bdr attr 4)

  // Lame coefficients (plane, matching LevelSetCantilever2D).
  constexpr Real muLame = 0.3846;
  constexpr Real lambdaLame = 0.5769;

  Vec2 vec2(Real x = 0, Real y = 0)
  {
    Vec2 out(2);
    out(0) = x;
    out(1) = y;
    return out;
  }

  Real applyPhaseMomentMap(Real phi, Real epsilon)
  {
    return std::tanh(phi / epsilon);
  }

  constexpr std::array<std::array<Real, 3>, 3> TriangleBarycentricQuadrature = {{
    {{ Real(2) / 3, Real(1) / 6, Real(1) / 6 }},
    {{ Real(1) / 6, Real(2) / 3, Real(1) / 6 }},
    {{ Real(1) / 6, Real(1) / 6, Real(2) / 3 }}
  }};

  Real facetLength(const WNGIRMesh& mesh, Index facet)
  {
    const auto face = mesh.getFace(facet);
    const auto& v = face->getVertices();
    return (mesh.getVertexCoordinates(v[1]) - mesh.getVertexCoordinates(v[0]))
      .norm();
  }

  // Distance from point p to segment [a,b].
  Real pointSegmentDistance(const Vec2& p, const Vec2& a, const Vec2& b)
  {
    const Vec2 ab = b - a;
    const Real len2 = ab.squaredNorm();
    Real t = len2 > Real(1e-30) ? ((p - a).dot(ab) / len2) : Real(0);
    t = std::max(Real(0), std::min(Real(1), t));
    return (p - (a + t * ab)).norm();
  }

  // Interface-residual (IR) reinit on a uniform lattice. Starting from the
  // closest-point field `a` (exact zero set, medial-axis kinks), run K steps of
  //   d_tau:  a <- a - dtau [ S0 (|grad a|_Godunov - 1) ]          (eikonal)
  //                       + IR force from the fitted interface segments        (pin)
  // The eikonal term cleans |grad a|->1 (smooths the kinks); the IR force
  // -2*lambda*a(seg)*ds scattered to the 4 surrounding lattice nodes pins the
  // zero set to the fitted interface (the int_Gamma a^2 gradient).
  // `a` is row-major a[i*ny + j] on a gnx x gny lattice with spacing h.
  void reinitIR(std::vector<Real>& a, int gnx, int gny, Real h,
                const std::vector<std::array<Vec2, 2>>& seg,
                int K, Real lambdaIR)
  {
    auto AT = [&](int i, int j) -> Real& { return a[i * gny + j]; };
    auto bilinear = [&](Real x, Real y) -> Real {
      Real gx = x / h, gy = y / h;
      int i = std::max(0, std::min(gnx - 2, (int)gx));
      int j = std::max(0, std::min(gny - 2, (int)gy));
      Real u = gx - i, v = gy - j;
      return (1 - u) * (1 - v) * AT(i, j) + u * (1 - v) * AT(i + 1, j)
           + (1 - u) * v * AT(i, j + 1) + u * v * AT(i + 1, j + 1);
    };
    std::vector<Real> S0(a.size());
    for (std::size_t k = 0; k < a.size(); ++k)
      S0[k] = a[k] / std::sqrt(a[k] * a[k] + h * h);
    const Real dtau = Real(0.3) * h;
    std::vector<Real> np = a;
    for (int it = 0; it < K; ++it)
    {
      np = a;
      for (int i = 1; i < gnx - 1; ++i)
        for (int j = 1; j < gny - 1; ++j)
        {
          const Real S = S0[i * gny + j];
          const Real am = (AT(i, j) - AT(i - 1, j)) / h;
          const Real ap = (AT(i + 1, j) - AT(i, j)) / h;
          const Real cm = (AT(i, j) - AT(i, j - 1)) / h;
          const Real cp = (AT(i, j + 1) - AT(i, j)) / h;
          auto sq = [](Real x) { return x * x; };
          Real gx2, gy2;
          if (S > 0) {
            gx2 = std::max(sq(std::max(am, Real(0))), sq(std::min(ap, Real(0))));
            gy2 = std::max(sq(std::max(cm, Real(0))), sq(std::min(cp, Real(0))));
          } else {
            gx2 = std::max(sq(std::min(am, Real(0))), sq(std::max(ap, Real(0))));
            gy2 = std::max(sq(std::min(cm, Real(0))), sq(std::max(cp, Real(0))));
          }
          np[i * gny + j] -= dtau * S * (std::sqrt(gx2 + gy2) - Real(1));
        }
      // IR pin: scatter the int_Gamma a^2 gradient from each fitted segment.
      for (const auto& s : seg)
      {
        const Vec2 m = Real(0.5) * (s[0] + s[1]);
        const Real ds = (s[1] - s[0]).norm();
        Real gx = m(0) / h, gy = m(1) / h;
        int i = (int)gx, j = (int)gy;
        if (i < 1 || i > gnx - 2 || j < 1 || j > gny - 2) continue;
        Real u = gx - i, v = gy - j;
        const Real f = -dtau * lambdaIR * Real(2) * bilinear(m(0), m(1)) * ds / h;
        np[i * gny + j]         += f * (1 - u) * (1 - v);
        np[(i + 1) * gny + j]   += f * u * (1 - v);
        np[i * gny + j + 1]     += f * (1 - u) * v;
        np[(i + 1) * gny + j + 1] += f * u * v;
      }
      a.swap(np);
    }
  }

  struct CellMomentInfo
  {
    Index index = 0;
    Real area = 0;
    Real moment = 0;
    std::array<Index, 3> vertices{};
    std::array<Vec2, 3> x;
  };

  // Phase moments of the carried discrete level set phiH (P1) per cell.
  template <class PhiGf>
  std::vector<CellMomentInfo> collectCellMomentInfo(
      const WNGIRMesh& mesh, const PhiGf& phiH, Real epsilon)
  {
    std::vector<CellMomentInfo> cells;
    cells.reserve(mesh.getCellCount());
    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const auto& cell = *cellIt;
      const auto& vertices = cell.getVertices();
      CellMomentInfo info;
      info.index = cell.getIndex();
      for (std::size_t i = 0; i < 3; ++i)
      {
        info.vertices[i] = vertices[i];
        info.x[i] = mesh.getVertexCoordinates(vertices[i]);
      }
      const Vec2 e1 = info.x[1] - info.x[0];
      const Vec2 e2 = info.x[2] - info.x[0];
      info.area = std::abs(Real(0.5) * (e1(0) * e2(1) - e1(1) * e2(0)));
      Real moment = 0;
      for (const auto& bary : TriangleBarycentricQuadrature)
      {
        Math::SpatialPoint rc(2);
        rc(0) = bary[1];
        rc(1) = bary[2];
        const Geometry::Point pt(cell, rc);
        moment += applyPhaseMomentMap(phiH.getValue(pt), epsilon);
      }
      info.moment = moment / TriangleBarycentricQuadrature.size();
      cells.push_back(std::move(info));
    }
    return cells;
  }

  template <class Displacement>
  void updateMovedMeshFromDisplacement(
      const WNGIRMesh& mesh, WNGIRMesh& moved, const Displacement& u)
  {
    const auto& uFes = u.getFiniteElementSpace();
    const auto& uData = u.getData();
    for (Index v = 0; v < mesh.getVertexCount(); ++v)
    {
      const Vec2 x = mesh.getVertexCoordinates(v);
      const auto& dofs = uFes.getDOFs(0, v);
      moved.setVertexCoordinates(v, vec2(x(0) + uData(dofs[0]),
                                         x(1) + uData(dofs[1])));
    }
  }

  std::size_t parseSizeTOption(int argc, char** argv, const std::string& name,
                               std::size_t fallback)
  {
    const std::string prefix = "--" + name + "=";
    for (int i = 1; i < argc; ++i)
      if (std::string(argv[i]).rfind(prefix, 0) == 0)
        return static_cast<std::size_t>(std::stoul(std::string(argv[i]).substr(prefix.size())));
    return fallback;
  }

  Real parseRealOption(int argc, char** argv, const std::string& name,
                       Real fallback)
  {
    const std::string prefix = "--" + name + "=";
    for (int i = 1; i < argc; ++i)
      if (std::string(argv[i]).rfind(prefix, 0) == 0)
        return static_cast<Real>(std::stod(std::string(argv[i]).substr(prefix.size())));
    return fallback;
  }

  std::string parseStringOption(int argc, char** argv, const std::string& name,
                                const std::string& fallback)
  {
    const std::string prefix = "--" + name + "=";
    for (int i = 1; i < argc; ++i)
      if (std::string(argv[i]).rfind(prefix, 0) == 0)
        return std::string(argv[i]).substr(prefix.size());
    return fallback;
  }

  struct StageTimer
  {
    using Clock = std::chrono::steady_clock;

    Clock::time_point t = Clock::now();

    Real reset()
    {
      const auto now = Clock::now();
      const Real out = std::chrono::duration<Real>(now - t).count();
      t = now;
      return out;
    }
  };
}

int main(int argc, char** argv)
{
  const Real L = parseRealOption(argc, argv, "L", Real(2));       // length
  const Real H = parseRealOption(argc, argv, "H", Real(1));       // height
  const std::size_t maxIt = parseSizeTOption(argc, argv, "iters", 300);
  const std::size_t n = parseSizeTOption(argc, argv, "n", 80);  // cells along H
  const Real ell = parseRealOption(argc, argv, "ell", Real(0.4));
  const Real h = H / static_cast<Real>(n);   // uniform cell size = metric h
  const Real hmin = parseRealOption(argc, argv, "hmin", Real(0.1) * h);
  const bool initialMMG =
    parseSizeTOption(argc, argv, "initial-mmg", 0) != 0;
  const Real initialMMGHMin =
    parseRealOption(argc, argv, "initial-mmg-hmin", hmin);
  const Real initialMMGHMax =
    parseRealOption(argc, argv, "initial-mmg-hmax", h);
  const Real initialMMGHausd =
    parseRealOption(argc, argv, "initial-mmg-hausd", Real(0.5) * hmin);
  // Design-smoothness knobs (sensitive — raise gently):
  //   alpha  : Hilbertian shape-gradient length. alpha~2h smooths the boundary
  //            while keeping members; alpha>~4h washes features out (holes fill).
  //   lambdaC: classifier perimeter penalty. Too large collapses the min-cut to
  //            a single phase (no interface). Keep ~0.004; nudge up for smoother.
  const Real alphaReg = parseRealOption(argc, argv, "alpha", Real(2) * h);
  const Real epsilon = parseRealOption(argc, argv, "classifier-eps", Real(1.25) * h);
  const Real lambdaC = parseRealOption(argc, argv, "classifier-lambda", Real(0.004));
  // Adaptive step: start at dt, grow x1.3 on accepted steps up to dtMax, halve
  // on a rejected (blown-up) step.
  const Real dt    = parseRealOption(argc, argv, "dt", Real(2) * h);
  const Real dtMax = parseRealOption(argc, argv, "dt-max", Real(4) * h);
  const bool objectiveLineSearch =
    parseSizeTOption(argc, argv, "objective-linesearch", 1) != 0;
  const Real objectiveDecreaseTol =
    parseRealOption(argc, argv, "objective-decrease-tol", Real(1e-10));
  // Empty (default) => uniform grid; pass --mesh=<file> to load instead.
  const std::string meshFile = parseStringOption(argc, argv, "mesh", "");
  // Stage-6 redistance:
  //   "cp" (default) : SUB-CELL reinit — signed closest-point distance to the
  //                    WNGIR-fitted interface segments, signed by the outward
  //                    normal of the nearest fitted segment (NOT the stale FMM
  //                    sign, which produced the spurious crossings / singular
  //                    phi). Carries the sub-cell fit; stable.
  //   "fmm"          : fast-marching eikonal reinit pinned at the CLASSIFIED
  //                    (cell-staircase) interface. Robust but not sub-cell.
  //   "ir"           : Li-Xu distance-regularized level-set reinitialization.
  //                    It starts from the clean FMM signed distance and uses
  //                    the WNGIR-fitted interface only through an implicit
  //                    surface pin, so the zero set is pulled sub-cell without
  //                    feeding closest-point kinks to the distance regularizer.
  //   "ir-phi"       : same Li-Xu distance regularizer, but the zero set is
  //                    pinned by the carried level set phi_h through the
  //                    coarea weight delta_eps(phi_h)|grad phi_h|. This does
  //                    not use the fitted/classified mesh interface as the
  //                    geometric reference.
  //   "none"         : no reinitialization; advect the current carried phi_h.
  //                    This is the cheapest path and avoids DRLSE texture, but
  //                    the level-set scale is then controlled only by advection.
  const std::string reinitMode = parseStringOption(argc, argv, "reinit", "cp");
  const std::size_t irSteps = parseSizeTOption(argc, argv, "ir-steps", 30);
  const Real irLambda = parseRealOption(argc, argv, "ir-lambda", Real(20));
  const Real irDtauFactor =
    parseRealOption(argc, argv, "ir-dtau-factor", Real(0.1));
  const Real irDeltaEps =
    parseRealOption(argc, argv, "ir-delta-eps", Real(1.5) * h);
  const std::size_t reinitEvery =
    parseSizeTOption(argc, argv, "reinit-every", 1);
  const std::string redistanceMode =
    parseStringOption(argc, argv, "redistance-mode", "none");
  const std::size_t redistanceEvery =
    parseSizeTOption(argc, argv, "redistance-every", 0);
  const std::size_t outputEvery =
    parseSizeTOption(argc, argv, "output-every", 1);
  const std::size_t irDpUpdateEvery =
    std::max<std::size_t>(1, parseSizeTOption(argc, argv, "ir-dp-update-every", 1));
  const std::size_t repairEvery =
    parseSizeTOption(argc, argv, "repair-every", 5);
  const Real repairEta = parseRealOption(argc, argv, "repair-eta", Real(0.5) * h);
  const Real repairBand = parseRealOption(argc, argv, "repair-band", Real(4) * h);
  const bool repairCalibrate =
    parseSizeTOption(argc, argv, "repair-calibrate", 0) != 0;
  const bool printTiming = parseSizeTOption(argc, argv, "timing", 0) != 0;
  const bool trace = parseSizeTOption(argc, argv, "trace", 0) != 0;

  // ---- Background grid: uniform by default, finer = more detail -----------
  WNGIRMesh mesh;
  if (meshFile.empty())
  {
    const std::size_t nx = static_cast<std::size_t>(std::lround(L / h)) + 1;
    const std::size_t ny = n + 1;
    mesh = WNGIRMesh::UniformGrid(Polytope::Type::Triangle, { nx, ny });
    mesh.scale(h);
  }
  else
  {
    mesh.load(meshFile, IO::FileFormat::MFEM);
  }
  mesh.getConnectivity().compute(1, 2);
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(2, 2);
  mesh.getConnectivity().compute(0, 0);
  const std::size_t D = mesh.getDimension();

  // Tag outer boundary: left edge clamped (GammaD), small right-edge segment
  // loaded (GammaN), the rest free (Gamma0).
  auto tagBoundary = [&](WNGIRMesh& m)
  {
    const std::size_t D = m.getDimension();
    for (auto it = m.getBoundary(); it; ++it)
    {
      const Index f = it->getIndex();
      const auto& vv = it->getVertices();
      const Vec2 a = m.getVertexCoordinates(vv[0]);
      const Vec2 b = m.getVertexCoordinates(vv[1]);
      const Real xm = Real(0.5) * (a(0) + b(0));
      const Real ym = Real(0.5) * (a(1) + b(1));
      Attribute attr = Gamma0;
      if (xm < Real(1e-9))
        attr = GammaD;                                   // left edge
      else if (xm > L - Real(1e-9)
               && ym > Real(0.45) * H && ym < Real(0.55) * H)
        attr = GammaN;                                   // right-middle load
      m.setAttribute({D - 1, f}, attr);
    }
  };
  tagBoundary(mesh);

  if (initialMMG)
  {
    std::cout << "  initial MMG remesh:"
              << " hmin=" << initialMMGHMin
              << " hmax=" << initialMMGHMax
              << " hausd=" << initialMMGHausd << "\n";
    MMG::Mesh mmgMesh(std::move(mesh));
    MMG::Optimizer optimizer;
    optimizer.setHMin(initialMMGHMin)
             .setHMax(initialMMGHMax)
             .setAngleDetection(false);
    if (initialMMGHausd > 0)
      optimizer.setHausdorff(initialMMGHausd);
    optimizer.optimize(mmgMesh);
    mesh = std::move(static_cast<MMG::Mesh::Parent&>(mmgMesh));
    mesh.getConnectivity().compute(1, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(2, 2);
    mesh.getConnectivity().compute(0, 0);
    tagBoundary(mesh);
  }

  using ScalarP1 = P1<Real, WNGIRMesh>;
  using ScalarP0 = P0<Real, WNGIRMesh>;
  using VectorP1 = P1<Math::SpatialVector<Real>, WNGIRMesh>;

  ScalarP1 sh(mesh);
  // Uniform-grid lattice index (i,j)->vertex, for the FD IR reinit.
  const int gnx = static_cast<int>(std::lround(L / h)) + 1;
  const int gny = static_cast<int>(std::lround(H / h)) + 1;
  std::vector<Index> gridVtx(static_cast<std::size_t>(gnx) * gny, 0);
  if (meshFile.empty())
    for (Index v = 0; v < mesh.getVertexCount(); ++v)
    {
      const Vec2 Xv = mesh.getVertexCoordinates(v);
      const int i = static_cast<int>(std::lround(Xv(0) / h));
      const int j = static_cast<int>(std::lround(Xv(1) / h));
      if (i >= 0 && i < gnx && j >= 0 && j < gny)
        gridVtx[static_cast<std::size_t>(i) * gny + j] = v;
    }
  ScalarP0 p0(mesh);
  VectorP1 vh(mesh, 2);

  GridFunction phiH(sh);          phiH.setName("phi");
  TrialFunction wngirTrial(vh);
  TestFunction  wngirTest(vh);
  auto& u = wngirTrial.getSolution();
  u.setName("displacement");      // WNGIR fit
  GridFunction dJ(vh);            dJ.setName("dJ");            // shape velocity

  // Reusable background-space objects. The background mesh connectivity is
  // fixed; only coefficients, attributes, and vertex coordinates are updated.
  VectorP1 gradFes(mesh, 2);
  GridFunction gradPhiSmoothed(gradFes);
  auto gradPhiDiscrete = Grad(phiH);
  TrialFunction gradTrial(gradFes);
  TestFunction  gradTest(gradFes);
  Problem gradProjection(gradTrial, gradTest);
  gradProjection = Integral(gradTrial, gradTest)
                 - Integral(gradPhiDiscrete, gradTest);

  GridFunction dist(sh);
  GridFunction speed(sh);
  GridFunction phiR(sh);
  GridFunction dp(p0);
  TrialFunction irTrial(sh);
  TestFunction  irTest(sh);
  BilinearForm irMass(irTrial, irTest);
  irMass = Integral(irTrial, irTest);
  irMass.assemble();
  const Math::Vector<Real> irMassLump =
    irMass.getOperator() * Math::Vector<Real>::Ones(phiH.getData().size());
  BilinearForm irDiffusion(irTrial, irTest);
  irDiffusion = Integral(dp * Grad(irTrial), Grad(irTest));
  BilinearForm repairStiffness(irTrial, irTest);
  repairStiffness = Integral(Grad(irTrial), Grad(irTest));
  repairStiffness.assemble();
  Math::SparseMatrix<Real> repairOperator = irMass.getOperator();
  repairOperator += (repairEta * repairEta) * repairStiffness.getOperator();
  Eigen::ConjugateGradient<Math::SparseMatrix<Real>,
    Eigen::Lower | Eigen::Upper> repairCG;
  repairCG.setTolerance(Real(1e-10));
  repairCG.setMaxIterations(static_cast<int>(phiH.getData().size()));
  repairCG.compute(repairOperator);

  TrialFunction adv(sh);
  TestFunction  advTest(sh);

  Rodin::Examples::WNGIRExampleDefaults wngirDefaults;
  wngirDefaults.maxIterations = 60;
  wngirDefaults.activeRMSOverHTol = Real(0.2);
  WNGIRParameters wp =
    Rodin::Examples::makeWNGIRParameters(
        argc, argv, h, Gamma, wngirDefaults);
  wp.trace = trace;
  WNGIR wngir(wngirTrial, wngirTest);
  wngir.setParameters(wp);

  // ---- Initial level set ---------------------------------------------------
  if (meshFile.empty())
  {
    // Uniform grid: analytic array of holes (classic cantilever seed).
    phiH = [&](const Geometry::Point& p) -> Real
    {
      const auto& X = p.getCoordinates();
      const Real nxh = 8, nyh = 4;
      const Real s = std::cos(nxh * M_PI * X(0) / L)
                   * std::cos(nyh * M_PI * X(1) / H);
      return -(s + Real(0.4));   // phi<0 (material) where s > -0.4
    };
  }
  else
  {
    // Loaded mesh: FMM from the interior/exterior labels it carries.
    for (auto faceIt = mesh.getFace(); faceIt; ++faceIt)
    {
      const Index f = faceIt->getIndex();
      const auto& inc = mesh.getConnectivity().getIncidence({1, 2}, f);
      if (inc.size() != 2) continue;
      const auto a = mesh.getAttribute(D, inc[0]);
      const auto b = mesh.getAttribute(D, inc[1]);
      const bool ia = a && *a == interiorAttribute;
      const bool ib = b && *b == interiorAttribute;
      if (ia != ib) mesh.setAttribute({D - 1, f}, Gamma);
    }
    phiH = Real(0);
    Distance::Eikonal<ScalarP1, Math::Vector<Real>>(phiH)
      .setInterior(interiorAttribute)
      .setInterface(Gamma)
      .solve()
      .sign();
  }

  IO::XDMF xdmf("LevelSetWNGIRCantilever2D");
  auto domainGrid = xdmf.grid("domain");
  domainGrid.setMesh(mesh, IO::XDMF::MeshPolicy::Transient);
  domainGrid.add(phiH, IO::XDMF::Center::Node);
  domainGrid.add(dJ,   IO::XDMF::Center::Node);

  WNGIRMesh moved(mesh);
  VectorP1 vhMoved(moved, 2);
  TrialFunction g(vhMoved);
  TestFunction  w(vhMoved);
  auto stateGrid = xdmf.grid("state");

  const std::string meshDescription =
    meshFile.empty()
      ? (initialMMG ? "uniform grid + initial MMG" : "uniform grid")
      : (initialMMG ? meshFile + " + initial MMG" : meshFile);
  std::cout << "WNGIR cantilever on " << mesh.getCellCount() << " cells ("
            << meshDescription << ")"
            << "\n  domain [0," << L << "]x[0," << H << "]"
            << "  ell=" << ell << "  alpha=" << alphaReg
            << "  h=" << h << "  dt=" << dt
            << "  objectiveLineSearch=" << objectiveLineSearch
            << "\n  WNGIR: gammaM=" << wp.gammaM
            << " gammaDev=" << wp.gammaH
            << " gammaDiv=" << wp.gammaDiv
            << " ellM=" << wp.ellM
            << " rmsTol=" << wp.activeRMSTol
            << " supTol=" << wp.activeSupTol
            << " steps=" << wp.maxIterations
            << "  reinit=" << reinitMode << "/" << reinitEvery
            << "  redistance=" << redistanceMode << "/" << redistanceEvery
            << "  repairEvery=" << repairEvery
            << "  outputEvery=" << outputEvery << '\n';
  std::ofstream fObj("LevelSetWNGIRCantilever2D.obj.txt");

  // Backtracking guard: body-fitted topology optimization can sever the load
  // path, making the elasticity solve singular (compliance -> inf/garbage).
  // We checkpoint the last good level set and, on a blown-up step, revert and
  // halve the step length.
  Math::Vector<Real> phiGood = phiH.getData();
  Real dtCur = dt;
  Real complianceCap = std::numeric_limits<Real>::infinity();  // set after it 0
  Real lastAcceptedObjective = std::numeric_limits<Real>::infinity();
  auto backgroundCellGradientMagnitude = [&](const auto& gf,
                                             const Polytope& cell) -> Real
  {
    const auto& vv = cell.getVertices();
    const Vec2 x0 = mesh.getVertexCoordinates(vv[0]);
    const Vec2 x1 = mesh.getVertexCoordinates(vv[1]);
    const Vec2 x2 = mesh.getVertexCoordinates(vv[2]);
    const Real p0v = gf.getData()(sh.getGlobalIndex({0, vv[0]}, 0));
    const Real p1v = gf.getData()(sh.getGlobalIndex({0, vv[1]}, 0));
    const Real p2v = gf.getData()(sh.getGlobalIndex({0, vv[2]}, 0));
    const Real J00 = x1(0) - x0(0), J01 = x2(0) - x0(0);
    const Real J10 = x1(1) - x0(1), J11 = x2(1) - x0(1);
    const Real det = J00 * J11 - J01 * J10;
    if (std::abs(det) < Real(1e-30))
      return Real(1);
    const Real g0 = p1v - p0v;
    const Real g1 = p2v - p0v;
    const Real gx = ( J11 * g0 - J10 * g1) / det;
    const Real gy = (-J01 * g0 + J00 * g1) / det;
    return std::sqrt(gx * gx + gy * gy);
  };

  std::size_t shapeAttempts = 0;
  std::size_t acceptedShapeIterations = 0;
  std::string stopReason = "completed-requested-accepted-iterations";
  while (acceptedShapeIterations < maxIt)
  {
    const std::size_t it = acceptedShapeIterations;
    ++shapeAttempts;
    std::cout << "\n--- Iteration " << it << " ---\n";
    StageTimer iterTimer;
    Real tClassify = 0, tGrad = 0, tWNGIR = 0, tMoveTrim = 0;
    Real tElasticity = 0, tHilbert = 0, tRedistance = 0, tAdvect = 0;
    Real tRepair = 0, tWrite = 0;

    // ---- Stage 1: classify the carried level set -------------------------
    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
      mesh.setAttribute({D, c}, Attribute{0});
    // Clear ALL faces too: otherwise last iteration's Gamma interface faces
    // persist, and FaceIntegral.over(Gamma)+traceOf(interior) hits a stale
    // face with no interior neighbor -> crash. Re-tag the outer boundary after.
    for (auto it = mesh.getFace(); it; ++it)
      mesh.setAttribute({D - 1, it->getIndex()}, Attribute{0});
    tagBoundary(mesh);

    const auto cellMoments = collectCellMomentInfo(mesh, phiH, epsilon);
    std::unordered_map<Index, std::size_t> cellToLocal;
    std::vector<Index> localToCell;
    cellToLocal.reserve(cellMoments.size());
    for (std::size_t l = 0; l < cellMoments.size(); ++l)
    {
      cellToLocal[cellMoments[l].index] = l;
      localToCell.push_back(cellMoments[l].index);
    }
    std::vector<Real> volumes(cellMoments.size()), moments(cellMoments.size());
    for (std::size_t l = 0; l < cellMoments.size(); ++l)
    {
      volumes[l] = cellMoments[l].area;
      moments[l] = cellMoments[l].moment;
    }

    std::vector<MinSTCut::Edge> graphEdges;
    for (auto faceIt = mesh.getFace(); faceIt; ++faceIt)
    {
      const Index facet = faceIt->getIndex();
      const auto& inc = mesh.getConnectivity().getIncidence({1, 2}, facet);
      if (inc.size() != 2) continue;
      const auto a = cellToLocal.find(inc[0]);
      const auto b = cellToLocal.find(inc[1]);
      if (a == cellToLocal.end() || b == cellToLocal.end()) continue;
      graphEdges.push_back({ static_cast<Index>(a->second),
                             static_cast<Index>(b->second),
                             lambdaC * facetLength(mesh, facet), facet });
    }

    const MinSTCut cut;
    const MinSTCut::Result classified = cut.classify(volumes, moments, graphEdges);

    for (std::size_t l = 0; l < classified.labels.size(); ++l)
      mesh.setAttribute({D, localToCell[l]},
                        classified.labels[l] == MinSTCut::Inside
                          ? interiorAttribute : exteriorAttribute);

    // ---- Connected-component cleanup -------------------------------------
    // Remove interior material not connected, through interior cells, to the
    // clamped support GammaD. Floating islands and load-only fragments carry
    // rigid-body modes that make the elasticity solve singular; dropping them
    // keeps every iteration well-posed (the body-fitted analogue of the weak
    // ersatz material used on fixed grids).
    {
      std::vector<char> isSupportCell(mesh.getCellCount(), 0);
      for (auto it = mesh.getBoundary(); it; ++it)
      {
        if (mesh.getAttribute(D - 1, it->getIndex()) != Attribute{GammaD})
          continue;
        for (const Index c :
             mesh.getConnectivity().getIncidence({1, 2}, it->getIndex()))
          isSupportCell[c] = 1;
      }
      const auto comps = mesh.ccl(
          [](const Polytope& a, const Polytope& b)
          { return a.getAttribute() == b.getAttribute(); }).getComponents();
      std::size_t removed = 0;
      for (const auto& comp : comps)
      {
        if (comp.empty()) continue;
        const Index rep = *comp.begin();
        const auto ra = mesh.getAttribute(D, rep);
        if (!ra || *ra != interiorAttribute) continue;   // exterior component
        bool anchored = false;
        for (const Index c : comp)
          if (isSupportCell[c]) { anchored = true; break; }
        if (!anchored)
          for (const Index c : comp)
          { mesh.setAttribute({D, c}, exteriorAttribute); ++removed; }
      }
      if (removed)
        std::cout << "  CCL: removed " << removed
                  << " floating interior cells (unsupported components)\n";
    }

    // ---- Interface = internal faces between interior and exterior --------
    std::vector<Index> interfaceFacets;
    for (auto faceIt = mesh.getFace(); faceIt; ++faceIt)
    {
      const Index f = faceIt->getIndex();
      const auto& inc = mesh.getConnectivity().getIncidence({1, 2}, f);
      if (inc.size() != 2) continue;
      const auto a = mesh.getAttribute(D, inc[0]);
      const auto b = mesh.getAttribute(D, inc[1]);
      const bool ia = a && *a == interiorAttribute;
      const bool ib = b && *b == interiorAttribute;
      if (ia != ib)
      {
        interfaceFacets.push_back(f);
        mesh.setAttribute({D - 1, f}, Gamma);
      }
    }

    std::size_t nInside = 0;
    Real areaInterior = 0;
    for (std::size_t l = 0; l < classified.labels.size(); ++l)
    {
      const auto a = mesh.getAttribute(D, localToCell[l]);
      if (a && *a == interiorAttribute)
      { ++nInside; areaInterior += cellMoments[l].area; }
    }

    std::cout << "  classify: interior cells=" << nInside
              << "  interface facets=" << interfaceFacets.size()
              << "  area=" << std::fixed << std::setprecision(4)
              << areaInterior << '\n';
    tClassify = iterTimer.reset();

    // ---- Stage 2: WNGIR fit mesh so skeleton lands on phi=0 --------------
    gradProjection.assemble();
    Solver::SparseLU(gradProjection).solve();
    gradPhiSmoothed.getData() = gradTrial.getSolution().getData();
    tGrad = iterTimer.reset();
    RealFunction phiFn(
        [&](const Geometry::Point& p) -> Real { return phiH.getValue(p); });
    AnalyticVectorFunction gradPhiFn(
        [&](const Geometry::Point& p) -> Math::SpatialVector<Real>
        { return gradPhiSmoothed.getValue(p); }, /*dim=*/2);

    u.getData().setZero();
    {
      const auto rep = wngir.solve(mesh, interfaceFacets, phiFn, gradPhiFn);
      // Diagnostic: did WNGIR actually move the mesh, and did it converge?
      Real maxUoverH = 0;
      const auto& uFes = u.getFiniteElementSpace();
      for (Index v = 0; v < mesh.getVertexCount(); ++v)
      {
        const auto& d = uFes.getDOFs(0, v);
        const Real ux = u.getData()(d[0]), uy = u.getData()(d[1]);
        maxUoverH = std::max(maxUoverH, std::sqrt(ux * ux + uy * uy) / h);
      }
      const Real activeRMSOverH =
        h > Real(0) ? rep.activeRMS / h : Real(0);
      const Real activeSupOverH =
        h > Real(0) ? rep.activeSup / h : Real(0);
      std::cout << "  WNGIR: it=" << rep.iterations
                << "  exit=" << rep.exitReason
                << "  activeRMS=" << std::scientific << std::setprecision(2)
                << rep.activeRMS
                << "  activeRMS/h=" << activeRMSOverH
                << "  activeSup/h=" << activeSupOverH
                << "  max|u|/h=" << maxUoverH << '\n';
    }
    tWNGIR = iterTimer.reset();

    // Build the body-fitted moved mesh and carry attributes.
    updateMovedMeshFromDisplacement(mesh, moved, u);
    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
      if (const auto a = mesh.getAttribute(D, c))
        moved.setAttribute({D, c}, *a);
    for (auto it = mesh.getFace(); it; ++it)
    {
      const Index f = it->getIndex();
      if (const auto a = mesh.getAttribute(D - 1, f))
        moved.setAttribute({D - 1, f}, *a);
    }

    // ---- Stage 3: elasticity on the fitted material submesh --------------
    SubMesh<Context::Local> trimmed = moved.trim(exteriorAttribute);
    trimmed.getConnectivity().compute(1, 2);
    tMoveTrim = iterTimer.reset();

    VectorP1 vhInt(trimmed, 2);
    auto fLoad = VectorFunction{ Real(0), Real(-1) };
    TrialFunction us(vhInt);
    TestFunction  vs(vhInt);
    // Tikhonov soft-foundation: a tiny mass term k_reg*(u,v) removes the
    // near-rigid-body null space that a one-cell-wide (nearly disconnected)
    // member produces, so the solve stays well-posed (bounded compliance)
    // instead of blowing up and freezing the optimizer. Body-fitted analogue
    // of the weak ersatz material. k_reg is tiny relative to the stiffness.
    Problem elasticity(us, vs);
    elasticity = LinearElasticityIntegral(us, vs)(lambdaLame, muLame)
               - BoundaryIntegral(fLoad, vs).over(GammaN)
               + DirichletBC(us, VectorFunction{ Real(0), Real(0) }).on(GammaD);
    Solver::CG(elasticity).solve();

    // Compliance = a(u,u).
    Real compliance = 0;
    {
      TrialFunction cu(vhInt);
      TestFunction  cv(vhInt);
      BilinearForm bf(cu, cv);
      bf = LinearElasticityIntegral(cu, cv)(lambdaLame, muLame);
      bf.assemble();
      compliance = bf(us.getSolution(), us.getSolution());
    }

    const Real objective = compliance + ell * areaInterior;
    tElasticity = iterTimer.reset();

    // Reject a blown-up step or, when requested, a non-descending objective
    // value. The monotonicity check is delayed by one outer iteration: the
    // trial design created by the previous advection is evaluated here, and
    // rejected by restoring the last accepted carrier and halving dt.
    const bool blewUp = !std::isfinite(compliance) || compliance < Real(0)
                        || compliance > complianceCap;
    const bool objectiveIncreased =
      objectiveLineSearch && std::isfinite(lastAcceptedObjective)
      && objective > lastAcceptedObjective
                   + objectiveDecreaseTol * std::max(Real(1), std::abs(lastAcceptedObjective));
    if (blewUp)
    {
      dtCur *= Real(0.5);
      phiH.getData() = phiGood;
      std::cout << "  REJECT (blow-up, compliance=" << std::scientific
                << std::setprecision(2) << compliance
                << "): reverting, dt -> " << dtCur << '\n';
      if (dtCur < Real(1e-3) * h)
      {
        stopReason = "dt-floor-after-blow-up-rejections";
        break;
      }
      continue;
    }
    if (objectiveIncreased)
    {
      dtCur *= Real(0.5);
      phiH.getData() = phiGood;
      std::cout << "  REJECT (objective increase, J=" << std::scientific
                << std::setprecision(4) << objective
                << " > " << lastAcceptedObjective
                << "): reverting, dt -> " << dtCur << '\n';
      if (dtCur < Real(1e-3) * h)
      {
        stopReason = "dt-floor-after-objective-rejections";
        break;
      }
      continue;
    }
    phiGood = phiH.getData();   // accept current level set as the new best
    lastAcceptedObjective = objective;
    dtCur = std::min(dtMax, dtCur * Real(1.3));   // grow the step on success
    if (it == 0)
      complianceCap = Real(50) * compliance;   // blow-up threshold

    fObj << objective << "\n"; fObj.flush();
    std::cout << "  compliance=" << std::scientific << std::setprecision(4)
              << compliance << "  area=" << areaInterior
              << "  J=" << objective << "  dt=" << dtCur << '\n';

    // ---- Stage 4: shape gradient + Hilbert regularization ----------------
    // The state lives on `trimmed` (a SubMesh of `moved`), so the gradient and
    // the Hilbert extension must be posed on `moved` (its parent) — exactly as
    // the MMG reference poses them on `th`. The strain energy density Ae:e on
    // the interface is lifted to a velocity field on `moved`, pinned only on
    // the structural Dirichlet boundary and on the Neumann/load boundary.
    auto jac = Jacobian(us.getSolution());
    jac.traceOf(interiorAttribute);
    auto e = Real(0.5) * (jac + jac.T());
    auto Ae = Real(2) * muLame * e + lambdaLame * Trace(e) * IdentityMatrix(2);
    auto nrm = FaceNormal(moved);
    nrm.traceOf(interiorAttribute);
    Problem hilbert(g, w);
    hilbert = Integral(alphaReg * alphaReg * Jacobian(g), Jacobian(w))
            + Integral(g, w)
            - FaceIntegral(Dot(Ae, e) - ell, Dot(nrm, w)).over(Gamma)
            + DirichletBC(g, VectorFunction{ Real(0), Real(0) }).on(GammaD)
            + DirichletBC(g, VectorFunction{ Real(0), Real(0) }).on(GammaN);
    Solver::CG(hilbert).solve();
    // Velocity is computed on `moved`; carry it to the fixed background grid
    // (same connectivity) to advect phi there. The moved->background transfer
    // is a coefficient copy; since the fit displacement is zero on the outer
    // boundary and the velocity is H1-smooth, the discrepancy is O(u) near the
    // interface (acceptable for the descent direction).
    dJ.getData() = g.getSolution().getData();
    tHilbert = iterTimer.reset();

    // ---- XDMF snapshot of the accepted objective state -------------------
    // Write before redistance/advection mutates phiH into the next trial
    // carrier. Thus every frame corresponds to a passed shape objective step.
    const bool writeSnapshot =
      (outputEvery > 0 && it % outputEvery == 0) || it + 1 == maxIt;
    if (writeSnapshot)
    {
      domainGrid.clear();
      domainGrid.add("phi", phiH, IO::XDMF::Center::Node);
      domainGrid.add("dJ", dJ, IO::XDMF::Center::Node);
      stateGrid.clear();
      stateGrid.setMesh(trimmed, IO::XDMF::MeshPolicy::Transient);
      stateGrid.add("u", us.getSolution(), IO::XDMF::Center::Node);
      xdmf.write(static_cast<Real>(it)).flush();
    }
    tWrite = iterTimer.reset();

    // ---- Stage 5: redistance from the WNGIR-FITTED interface -------------
    // The FMM gives the correct *sign* field and far distances, but its zero
    // set is the cell-quantized classifier staircase. Inside a narrow band we
    // overwrite the magnitude with the exact distance to the WNGIR-fitted
    // interface facets (moved-mesh segments at X+u), so phi carries the
    // sub-cell interface position forward instead of re-quantizing every step.
    std::string activeReinitMode = reinitMode;
    bool doReinit =
      reinitMode != "none"
      && reinitEvery > 0
      && (it % reinitEvery == 0 || it + 1 == maxIt);
    if (!doReinit && reinitMode == "none" && redistanceMode != "none"
        && redistanceEvery > 0
        && ((it + 1) % redistanceEvery == 0 || it + 1 == maxIt))
    {
      activeReinitMode = redistanceMode;
      doReinit = true;
    }
    if (!doReinit)
    {
      dist.getData() = phiH.getData();
    }
    else
    {
      dist = Real(0);
      Distance::Eikonal<ScalarP1, Math::Vector<Real>>(dist)
        .setInterior(interiorAttribute)
        .setInterface(Gamma)
        .solve()
        .sign();
      const Math::Vector<Real> fmmDist = dist.getData();
      // reinit=fmm: use the plain FMM signed distance (the reference example's
      // reinit) with NO closest-point/IR overwrite — diagnostic baseline.
      if (activeReinitMode != "fmm")
      {
      // Each fitted segment gets an OUTWARD unit normal (pointing from the
      // interior neighbour to the exterior), so the closest-point sign is taken
      // from the fitted interface itself, NOT from the stale FMM/staircase sign.
      // Using the FMM sign caused a spurious zero crossing (phantom sliver ->
      // singular elasticity) wherever the fit moved the interface across a node.
      std::vector<std::array<Vec2, 2>> fittedSeg;
      std::vector<Vec2> fittedN;
      fittedSeg.reserve(interfaceFacets.size());
      fittedN.reserve(interfaceFacets.size());
      for (const Index f : interfaceFacets)
      {
        const auto& vv = moved.getFace(f)->getVertices();
        const Vec2 a = moved.getVertexCoordinates(vv[0]);
        const Vec2 b = moved.getVertexCoordinates(vv[1]);
        Vec2 nrm = vec2(b(1) - a(1), -(b(0) - a(0)));        // perp to segment
        const Real nn = nrm.norm();
        if (nn > Real(1e-30)) nrm = nrm / nn; else nrm = vec2(1, 0);
        // orient away from the interior incident cell's (moved) centroid
        const auto& inc = mesh.getConnectivity().getIncidence({1, 2}, f);
        Vec2 cInt = vec2(0, 0); bool haveInt = false;
        for (const Index c : inc)
        {
          const auto at = mesh.getAttribute(D, c);
          if (at && *at == interiorAttribute)
          {
            const auto& cvv = mesh.getCell(c)->getVertices();
            Vec2 ctr = vec2(0, 0);
            for (const Index cv : cvv) ctr += moved.getVertexCoordinates(cv);
            cInt = ctr / static_cast<Real>(cvv.size()); haveInt = true; break;
          }
        }
        const Vec2 mid = Real(0.5) * (a + b);
        if (haveInt && (mid - cInt).dot(nrm) < Real(0)) nrm = -nrm;
        fittedSeg.push_back({ a, b });
        fittedN.push_back(nrm);
      }
      auto overwriteClosestPointBand = [&]()
      {
        const Real band = Real(6) * h;
        for (Index v = 0; v < mesh.getVertexCount(); ++v)
        {
          const Index d = sh.getGlobalIndex({0, v}, 0);
          if (std::abs(fmmDist(d)) > band) continue;
          const Vec2 X = mesh.getVertexCoordinates(v);
          Real dmin = std::numeric_limits<Real>::infinity();
          std::size_t kmin = 0;
          for (std::size_t k = 0; k < fittedSeg.size(); ++k)
          {
            const Real dk =
              pointSegmentDistance(X, fittedSeg[k][0], fittedSeg[k][1]);
            if (dk < dmin) { dmin = dk; kmin = k; }
          }
          if (std::isfinite(dmin))
          {
            const Vec2 mid =
              Real(0.5) * (fittedSeg[kmin][0] + fittedSeg[kmin][1]);
            const Real s = (X - mid).dot(fittedN[kmin]) >= Real(0)
                             ? Real(1) : Real(-1);   // outward = exterior = +
            dist.getData()(d) = s * dmin;
          }
        }
      };

      auto cellGradientMagnitude = [&](const auto& gf,
                                       const Polytope& cell) -> Real
      {
        const auto& vv = cell.getVertices();
        const Vec2 x0 = mesh.getVertexCoordinates(vv[0]);
        const Vec2 x1 = mesh.getVertexCoordinates(vv[1]);
        const Vec2 x2 = mesh.getVertexCoordinates(vv[2]);
        const Real p0v = gf.getData()(sh.getGlobalIndex({0, vv[0]}, 0));
        const Real p1v = gf.getData()(sh.getGlobalIndex({0, vv[1]}, 0));
        const Real p2v = gf.getData()(sh.getGlobalIndex({0, vv[2]}, 0));
        const Real J00 = x1(0) - x0(0), J01 = x2(0) - x0(0);
        const Real J10 = x1(1) - x0(1), J11 = x2(1) - x0(1);
        const Real det = J00 * J11 - J01 * J10;
        if (std::abs(det) < Real(1e-30))
          return Real(1);
        const Real g0 = p1v - p0v;
        const Real g1 = p2v - p0v;
        const Real gx = ( J11 * g0 - J10 * g1) / det;
        const Real gy = (-J01 * g0 + J00 * g1) / det;
        return std::sqrt(gx * gx + gy * gy);
      };

      auto gradientStats = [&](const auto& gf, bool cutCellsOnly)
      {
        Real gmin = std::numeric_limits<Real>::infinity();
        Real gmax = Real(0);
        const Real band = Real(4) * h;
        std::size_t cnt = 0;
        for (auto cit = mesh.getCell(); cit; ++cit)
        {
          const auto& vv = cit->getVertices();
          Real av = Real(0);
          Real vmin = std::numeric_limits<Real>::infinity();
          Real vmax = -std::numeric_limits<Real>::infinity();
          for (const Index v : vv)
          {
            const Real pv = gf.getData()(sh.getGlobalIndex({0, v}, 0));
            av += std::abs(pv);
            vmin = std::min(vmin, pv);
            vmax = std::max(vmax, pv);
          }
          av /= static_cast<Real>(vv.size());
          if (cutCellsOnly)
          {
            if (!(vmin <= Real(0) && vmax >= Real(0)))
              continue;
          }
          else if (av > band)
          {
            continue;
          }
          const Real gm = cellGradientMagnitude(gf, *cit);
          gmin = std::min(gmin, gm);
          gmax = std::max(gmax, gm);
          ++cnt;
        }
        if (cnt == 0)
          return std::array<Real, 2>{ Real(0), Real(0) };
        return std::array<Real, 2>{ gmin, gmax };
      };

      if (activeReinitMode == "ir" || activeReinitMode == "ir-phi")
      {
        // Li-Xu DRLSE with the full signed bounded double-well diffusivity
        //
        //   d_p(s) = 1 - 1/s,                 s >= 1,
        //          = sin(2*pi*s)/(2*pi*s),    s <  1,
        //
        // and an implicit interface pin. For reinit=ir, the pin is
        // int_{Gamma_fit} phi^2 ds on the moved WNGIR segments. For
        // reinit=ir-phi, the pin is the coarea approximation
        // int delta_eps(phi_h) |grad phi_h| phi^2 dx, hence the geometric
        // reference is the level set itself and not the classified mesh.
        const Real dtauR = irDtauFactor * h * h;
        const int  Kr    = static_cast<int>(irSteps);
        if (activeReinitMode == "ir")
        {
          // Mesh-pinned variant: start from the clean FMM SDF and use the
          // fitted mesh only as a zero-trace constraint.
          dist.getData() = fmmDist;
          phiR.getData() = dist.getData();
        }
        else
        {
          // Phi-pinned variant: start from the carried level set, rescaled by
          // the average cut-cell gradient so the unknown has distance units.
          Real gsum = Real(0);
          std::size_t gcnt = 0;
          for (auto cit = mesh.getCell(); cit; ++cit)
          {
            const auto& vv = cit->getVertices();
            Real vmin = std::numeric_limits<Real>::infinity();
            Real vmax = -std::numeric_limits<Real>::infinity();
            for (const Index v : vv)
            {
              const Real pv = phiH.getData()(sh.getGlobalIndex({0, v}, 0));
              vmin = std::min(vmin, pv);
              vmax = std::max(vmax, pv);
            }
            if (vmin <= Real(0) && vmax >= Real(0))
            {
              gsum += cellGradientMagnitude(phiH, *cit);
              ++gcnt;
            }
          }
          const Real gscale =
            gcnt > 0 ? std::max(gsum / static_cast<Real>(gcnt), Real(1e-12))
                     : Real(1);
          phiR.getData() = phiH.getData() / gscale;
        }
        std::vector<Eigen::Triplet<Real>> pinTriplets;

        auto addPinPoint = [&](const Vec2& y, Real w)
        {
          Math::SpatialPoint pc(2);
          pc(0) = y(0);
          pc(1) = y(1);
          constexpr Real tol = Real(1e-10);
          for (auto cit = mesh.getCell(); cit; ++cit)
          {
            Math::SpatialPoint rc(2);
            cit->getTransformation().inverse(rc, pc);
            const Real l0 = Real(1) - rc(0) - rc(1);
            const Real l1 = rc(0);
            const Real l2 = rc(1);
            if (l0 < -tol || l1 < -tol || l2 < -tol)
              continue;
            const auto& vv = cit->getVertices();
            const std::array<Real, 3> N = {{ l0, l1, l2 }};
            for (std::size_t a = 0; a < 3; ++a)
            {
              const Index ia = sh.getGlobalIndex({0, vv[a]}, 0);
              for (std::size_t b = 0; b < 3; ++b)
              {
                const Index ib = sh.getGlobalIndex({0, vv[b]}, 0);
                pinTriplets.emplace_back(
                    static_cast<int>(ia),
                    static_cast<int>(ib),
                    w * N[a] * N[b]);
              }
            }
            return;
          }
        };

        if (activeReinitMode == "ir")
        {
          pinTriplets.reserve(18 * fittedSeg.size());
          // Two-point Gauss quadrature on every moved interface segment.
          constexpr Real q0 = Real(0.21132486540518711775);
          constexpr Real q1 = Real(0.78867513459481288225);
          for (const auto& seg : fittedSeg)
          {
            const Vec2 a = seg[0];
            const Vec2 b = seg[1];
            const Real ds = (b - a).norm();
            addPinPoint((Real(1) - q0) * a + q0 * b, Real(0.5) * ds);
            addPinPoint((Real(1) - q1) * a + q1 * b, Real(0.5) * ds);
          }
        }
        else
        {
          pinTriplets.reserve(27 * mesh.getCellCount());
          Real gsum = Real(0);
          std::size_t gcnt = 0;
          for (auto cit = mesh.getCell(); cit; ++cit)
          {
            const auto& vv = cit->getVertices();
            Real vmin = std::numeric_limits<Real>::infinity();
            Real vmax = -std::numeric_limits<Real>::infinity();
            for (const Index v : vv)
            {
              const Real pv = phiH.getData()(sh.getGlobalIndex({0, v}, 0));
              vmin = std::min(vmin, pv);
              vmax = std::max(vmax, pv);
            }
            if (vmin <= Real(0) && vmax >= Real(0))
            {
              gsum += cellGradientMagnitude(phiH, *cit);
              ++gcnt;
            }
          }
          const Real gavg =
            gcnt > 0 ? std::max(gsum / static_cast<Real>(gcnt), Real(1e-12))
                     : Real(1);
          const Real epsPhi = std::max(irDeltaEps * gavg, Real(1e-12));
          for (auto cit = mesh.getCell(); cit; ++cit)
          {
            const auto& vv = cit->getVertices();
            const Vec2 x0 = mesh.getVertexCoordinates(vv[0]);
            const Vec2 x1 = mesh.getVertexCoordinates(vv[1]);
            const Vec2 x2 = mesh.getVertexCoordinates(vv[2]);
            const Real area = std::abs(Real(0.5) *
              ((x1(0) - x0(0)) * (x2(1) - x0(1))
             - (x1(1) - x0(1)) * (x2(0) - x0(0))));
            const Real gmag = cellGradientMagnitude(phiH, *cit);
            std::array<Index, 3> gd;
            std::array<Real, 3> pv;
            for (std::size_t a = 0; a < 3; ++a)
            {
              gd[a] = sh.getGlobalIndex({0, vv[a]}, 0);
              pv[a] = phiH.getData()(gd[a]);
            }
            for (const auto& bary : TriangleBarycentricQuadrature)
            {
              const std::array<Real, 3> N = {{ bary[0], bary[1], bary[2] }};
              const Real phiq = N[0] * pv[0] + N[1] * pv[1] + N[2] * pv[2];
              if (std::abs(phiq) >= epsPhi)
                continue;
              const Real delta =
                (Real(1) / (Real(2) * epsPhi))
                * (Real(1) + std::cos(M_PI * phiq / epsPhi));
              const Real w = (area / Real(3)) * delta * gmag;
              for (std::size_t a = 0; a < 3; ++a)
                for (std::size_t b = 0; b < 3; ++b)
                  pinTriplets.emplace_back(
                      static_cast<int>(gd[a]),
                      static_cast<int>(gd[b]),
                      w * N[a] * N[b]);
            }
          }
        }

        Math::SparseMatrix<Real> A(phiR.getData().size(), phiR.getData().size());
        std::vector<Eigen::Triplet<Real>> aTriplets;
        aTriplets.reserve(
            static_cast<std::size_t>(phiR.getData().size()) + pinTriplets.size());
        for (Eigen::Index i = 0; i < phiR.getData().size(); ++i)
          aTriplets.emplace_back(
              static_cast<int>(i), static_cast<int>(i), irMassLump(i) / dtauR);
        for (const auto& t : pinTriplets)
          aTriplets.emplace_back(t.row(), t.col(), Real(2) * irLambda * t.value());
        A.setFromTriplets(aTriplets.begin(), aTriplets.end());
        Eigen::ConjugateGradient<Math::SparseMatrix<Real>,
          Eigen::Lower | Eigen::Upper> irCG;
        irCG.setTolerance(Real(1e-10));
        irCG.setMaxIterations(static_cast<int>(phiR.getData().size()));
        irCG.compute(A);
        const auto gInCut = gradientStats(phiR, true);
        const auto gInBand = gradientStats(phiR, false);
        Math::Vector<Real> rK(phiR.getData().size());
        for (int kk = 0; kk < Kr; ++kk)
        {
          if (kk == 0 || static_cast<std::size_t>(kk) % irDpUpdateEvery == 0)
          {
            for (auto cit = mesh.getCell(); cit; ++cit)
            {
              const Index dc = p0.getGlobalIndex({D, cit->getIndex()}, 0);
              const Real s = cellGradientMagnitude(phiR, *cit);
              Real dpv;
              if (s >= Real(1)) dpv = Real(1) - Real(1)/s;
              else { const Real x = Real(2)*M_PI*s;
                     dpv = s > Real(1e-6) ? std::sin(x)/x : Real(1); }
              dp.getData()(dc) = dpv;
            }
            irDiffusion.assemble();
          }
          rK = irDiffusion.getOperator() * phiR.getData();
          Math::Vector<Real> rhs(phiR.getData().size());
          for (Eigen::Index i = 0; i < phiR.getData().size(); ++i)
            rhs(i) = (irMassLump(i) / dtauR) * phiR.getData()(i) - rK(i);
          phiR.getData() = irCG.solve(rhs);
        }
        const auto gOutCut = gradientStats(phiR, true);
        const auto gOutBand = gradientStats(phiR, false);
        std::cout << "    IR-DRLSE"
                  << (activeReinitMode == "ir-phi" ? "-phi" : "")
                  << ": |grad phi| cut ["
                  << std::fixed << std::setprecision(3)
                  << gInCut[0] << "," << gInCut[1] << "] -> ["
                  << gOutCut[0] << "," << gOutCut[1] << "]"
                  << "  band [" << gInBand[0] << "," << gInBand[1]
                  << "] -> [" << gOutBand[0] << "," << gOutBand[1] << "]\n";
        dist.getData() = phiR.getData();
      }
      else
      {
        overwriteClosestPointBand();
      }
      }
    }
    tRedistance = iterTimer.reset();

    // Normalize the descent velocity, advect distance by one step.
    speed = Frobenius(dJ);
    const Real vmax = std::max(speed.max(), Real(1e-30));
    dJ /= vmax;
    Advection::Lagrangian(adv, advTest, dist, dJ).step(dtCur);
    phiH.getData() = adv.getSolution().getData();
    tAdvect = iterTimer.reset();

    const bool doRepair =
      repairEvery > 0 && ((it + 1) % repairEvery == 0 || it + 1 == maxIt);
    if (doRepair)
    {
      const Math::Vector<Real> rhs = irMass.getOperator() * phiH.getData();
      phiH.getData() = repairCG.solve(rhs);
      if (repairCalibrate)
      {
        Real gsum = Real(0);
        Real wsum = Real(0);
        for (auto cit = mesh.getCell(); cit; ++cit)
        {
          const auto& vv = cit->getVertices();
          Real av = Real(0);
          for (const Index v : vv)
            av += std::abs(phiH.getData()(sh.getGlobalIndex({0, v}, 0)));
          av /= static_cast<Real>(vv.size());
          if (av > repairBand)
            continue;
          const Vec2 x0 = mesh.getVertexCoordinates(vv[0]);
          const Vec2 x1 = mesh.getVertexCoordinates(vv[1]);
          const Vec2 x2 = mesh.getVertexCoordinates(vv[2]);
          const Real area = std::abs(Real(0.5) *
            ((x1(0) - x0(0)) * (x2(1) - x0(1))
           - (x1(1) - x0(1)) * (x2(0) - x0(0))));
          gsum += area * backgroundCellGradientMagnitude(phiH, *cit);
          wsum += area;
        }
        if (wsum > Real(0))
        {
          const Real c = std::max(gsum / wsum, Real(1e-12));
          phiH.getData() /= c;
          std::cout << "    repair: scale=" << std::scientific
                    << std::setprecision(3) << c << '\n';
        }
      }
    }
    tRepair = iterTimer.reset();

    if (printTiming)
    {
      const Real total = tClassify + tGrad + tWNGIR + tMoveTrim
                       + tElasticity + tHilbert + tRedistance + tAdvect
                       + tRepair + tWrite;
      std::cout << "  timing:"
                << " classify=" << std::scientific << std::setprecision(2) << tClassify
                << " grad=" << tGrad
                << " wngir=" << tWNGIR
                << " moveTrim=" << tMoveTrim
                << " elas=" << tElasticity
                << " hilbert=" << tHilbert
                << " redist=" << tRedistance
                << " advect=" << tAdvect
                << " repair=" << tRepair
                << " write=" << tWrite
                << " total=" << total << '\n';
    }
    ++acceptedShapeIterations;
  }

  xdmf.close();
  std::cout << "\nDone. Accepted shape iterations: "
            << acceptedShapeIterations << " / " << maxIt
            << "  attempts=" << shapeAttempts
            << "  stop=" << stopReason
            << ". Objective history in LevelSetWNGIRCantilever2D.obj.txt\n";
  return 0;
}
