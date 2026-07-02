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
// the level set is optionally redistanced and advected by the descent velocity.
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
#include <limits>
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

#ifdef RODIN_WNGIR_P2_DISPLACEMENT
  void installP2Transformations(WNGIRMesh& mesh)
  {
    Math::SpatialPoint X;
    const std::size_t D = mesh.getDimension();
    const std::size_t sdim = mesh.getSpaceDimension();
    for (std::size_t d = 1; d <= D; ++d)
    {
      for (auto pit = mesh.getPolytope(d); pit; ++pit)
      {
        const auto& polytope = *pit;
        Variational::RealH1Element<2> geomFe(polytope.getGeometry());
        Geometry::PointCloud pm(sdim, geomFe.getCount());
        for (std::size_t a = 0; a < geomFe.getCount(); ++a)
        {
          const auto& rc = geomFe.getNode(a);
          polytope.getTransformation().transform(X, rc);
          for (std::size_t c = 0; c < sdim; ++c)
            pm(c, a) = X(c);
        }
        mesh.setPolytopeTransformation(
            {d, polytope.getIndex()},
            new Geometry::ParametricTransformation<Variational::RealH1Element<2>>(
              std::move(pm), geomFe));
      }
    }
  }
#endif

  // Distance from point p to segment [a,b].
  Real pointSegmentDistance(const Vec2& p, const Vec2& a, const Vec2& b)
  {
    const Vec2 ab = b - a;
    const Real len2 = ab.squaredNorm();
    Real t = len2 > Real(1e-30) ? ((p - a).dot(ab) / len2) : Real(0);
    t = std::max(Real(0), std::min(Real(1), t));
    return (p - (a + t * ab)).norm();
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
#ifdef RODIN_WNGIR_P2_DISPLACEMENT
    Math::SpatialPoint X;
    const std::size_t D = mesh.getDimension();
    const std::size_t sdim = mesh.getSpaceDimension();
    for (std::size_t d = 1; d <= D; ++d)
    {
      for (auto pit = mesh.getPolytope(d); pit; ++pit)
      {
        const auto& polytope = *pit;
        Variational::RealH1Element<2> geomFe(polytope.getGeometry());
        Geometry::PointCloud pm(sdim, geomFe.getCount());
        for (std::size_t a = 0; a < geomFe.getCount(); ++a)
        {
          const auto& rc = geomFe.getNode(a);
          polytope.getTransformation().transform(X, rc);
          for (std::size_t c = 0; c < sdim; ++c)
          {
            const Real uc =
              uData(uFes.getGlobalIndex({d, polytope.getIndex()}, a * sdim + c));
            pm(c, a) = X(c) + uc;
          }
        }
        moved.setPolytopeTransformation(
            {d, polytope.getIndex()},
            new Geometry::ParametricTransformation<Variational::RealH1Element<2>>(
              std::move(pm), geomFe));
      }
    }
#endif
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
  const std::string shapeDerivativeMode =
    parseStringOption(argc, argv, "shape-derivative", "volume");
  // Adaptive redistance is enabled by default. Periodic redistance is only
  // requested through --redistance-every; use --redistance-mode=none to disable.
  const std::size_t classifyEvery =
    parseSizeTOption(argc, argv, "classify-every", 1);
  const std::string defaultRedistanceMode =
#ifdef RODIN_WNGIR_P2_DISPLACEMENT
    "cp";
#else
    "fmm";
#endif
  const std::string redistanceMode =
    parseStringOption(argc, argv, "redistance-mode", defaultRedistanceMode);
  const std::size_t redistanceEvery =
    parseSizeTOption(argc, argv, "redistance-every", 0);
  const bool adaptiveRedistance =
    parseSizeTOption(argc, argv, "redistance-adaptive", 1) != 0;
  const Real redistanceEikonalTol =
    parseRealOption(argc, argv, "redistance-eikonal-tol", Real(0.35));
  const std::size_t outputEvery =
    parseSizeTOption(argc, argv, "output-every", 1);
  const std::size_t repairEvery =
    parseSizeTOption(argc, argv, "repair-every", 0);
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
#ifdef RODIN_WNGIR_P2_DISPLACEMENT
  installP2Transformations(mesh);
#endif

#ifdef RODIN_WNGIR_P2_DISPLACEMENT
  using ScalarFES = H1<2, Real, WNGIRMesh>;
  using VectorFES = H1<2, Math::SpatialVector<Real>, WNGIRMesh>;
#else
  using ScalarFES = P1<Real, WNGIRMesh>;
  using VectorFES = P1<Math::SpatialVector<Real>, WNGIRMesh>;
#endif
  using ScalarP1 = P1<Real, WNGIRMesh>;
  using VectorP1 = P1<Math::SpatialVector<Real>, WNGIRMesh>;

  ScalarFES sh(
#ifdef RODIN_WNGIR_P2_DISPLACEMENT
      std::integral_constant<std::size_t, 2>{},
#endif
      mesh);
  VectorFES vh(
#ifdef RODIN_WNGIR_P2_DISPLACEMENT
      std::integral_constant<std::size_t, 2>{},
#endif
      mesh, 2);

  GridFunction phiH(sh);          phiH.setName("phi");
  TrialFunction wngirTrial(vh);
  TestFunction  wngirTest(vh);
  auto& u = wngirTrial.getSolution();
  u.setName("displacement");      // WNGIR fit
  GridFunction dJ(vh);            dJ.setName("dJ");            // shape velocity

  // Reusable background-space objects. The background mesh connectivity is
  // fixed; only coefficients, attributes, and vertex coordinates are updated.
  VectorFES gradFes(
#ifdef RODIN_WNGIR_P2_DISPLACEMENT
      std::integral_constant<std::size_t, 2>{},
#endif
      mesh, 2);
  GridFunction gradPhiSmoothed(gradFes);
  auto gradPhiDiscrete = Grad(phiH);
  TrialFunction gradTrial(gradFes);
  TestFunction  gradTest(gradFes);
  Problem gradProjection(gradTrial, gradTest);
  gradProjection = Integral(gradTrial, gradTest)
                 - Integral(gradPhiDiscrete, gradTest);

  GridFunction dist(sh);
  GridFunction speed(sh);
  TrialFunction repairTrial(sh);
  TestFunction  repairTest(sh);
  BilinearForm repairMass(repairTrial, repairTest);
  BilinearForm repairStiffness(repairTrial, repairTest);
  Eigen::ConjugateGradient<Math::SparseMatrix<Real>,
    Eigen::Lower | Eigen::Upper> repairCG;
  repairCG.setTolerance(Real(1e-10));
  repairCG.setMaxIterations(static_cast<int>(phiH.getData().size()));
  if (repairEvery > 0)
  {
    repairMass = Integral(repairTrial, repairTest);
    repairMass.assemble();
    repairStiffness = Integral(Grad(repairTrial), Grad(repairTest));
    repairStiffness.assemble();
    Math::SparseMatrix<Real> repairOperator = repairMass.getOperator();
    repairOperator += (repairEta * repairEta) * repairStiffness.getOperator();
    repairCG.compute(repairOperator);
  }

  TrialFunction adv(sh);
  TestFunction  advTest(sh);

  Rodin::Examples::WNGIRExampleDefaults wngirDefaults;
  wngirDefaults.maxIterations = 60;
  wngirDefaults.activeRMSOverHTol = Real(0.05);
  wngirDefaults.activeSupOverHTol = Real(0.25);
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
    ScalarP1 shP1(mesh);
    GridFunction phiP1(shP1);
    phiP1 = Real(0);
    Distance::Eikonal<ScalarP1, Math::Vector<Real>>(phiP1)
      .setInterior(interiorAttribute)
      .setInterface(Gamma)
      .solve()
      .sign();
    phiH = [&](const Geometry::Point& p) -> Real
    {
      return phiP1.getValue(p);
    };
  }

  IO::XDMF xdmf("LevelSetWNGIRCantilever2D");
  auto domainGrid = xdmf.grid("domain");
  domainGrid.setMesh(mesh, IO::XDMF::MeshPolicy::Transient);
  domainGrid.add(phiH, IO::XDMF::Center::Node);
  domainGrid.add(dJ,   IO::XDMF::Center::Node);

  WNGIRMesh moved(mesh);
  VectorFES vhMoved(
#ifdef RODIN_WNGIR_P2_DISPLACEMENT
      std::integral_constant<std::size_t, 2>{},
#endif
      moved, 2);
  TrialFunction g(vhMoved);
  TestFunction  w(vhMoved);
  auto fittedGrid = xdmf.grid("fitted");
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
            << "  classify=" << classifyEvery
            << "  redistance=" << redistanceMode << "/" << redistanceEvery
            << " adaptive=" << adaptiveRedistance
            << " eikTol=" << redistanceEikonalTol
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
  auto weightedEikonalDefect = [&](const auto& gf, Real eta) -> Real
  {
    const Real eta2 = std::max(eta * eta, Real(1e-30));
    Real weightedError = 0;
    Real weightSum = 0;
    for (auto cit = mesh.getCell(); cit; ++cit)
    {
      const auto& vv = cit->getVertices();
      const Vec2 x0 = mesh.getVertexCoordinates(vv[0]);
      const Vec2 x1 = mesh.getVertexCoordinates(vv[1]);
      const Vec2 x2 = mesh.getVertexCoordinates(vv[2]);
      const Real area = std::abs(Real(0.5)
        * ((x1(0) - x0(0)) * (x2(1) - x0(1))
         - (x1(1) - x0(1)) * (x2(0) - x0(0))));
      Real phiAvg = 0;
      for (const Index v : vv)
        phiAvg += gf.getData()(sh.getGlobalIndex({0, v}, 0));
      phiAvg /= static_cast<Real>(vv.size());
      const Real omega = std::exp(-(phiAvg * phiAvg) / eta2);
      const Real gradMag = backgroundCellGradientMagnitude(gf, *cit);
      weightedError += area * omega * (gradMag - Real(1)) * (gradMag - Real(1));
      weightSum += area * omega;
    }
    return weightSum > Real(0) ? std::sqrt(weightedError / weightSum) : Real(0);
  };

  std::size_t shapeAttempts = 0;
  std::size_t acceptedShapeIterations = 0;
  std::string stopReason = "completed-requested-accepted-iterations";
  bool haveClassification = false;
  std::vector<char> previousInterior(mesh.getCellCount(), 0);
  std::vector<Index> cachedInterfaceFacets;
  Real cachedAreaInterior = 0;
  std::size_t cachedInsideCount = 0;
  while (acceptedShapeIterations < maxIt)
  {
    const std::size_t it = acceptedShapeIterations;
    ++shapeAttempts;
    std::cout << "\n--- Iteration " << it << " ---\n";
    StageTimer iterTimer;
    Real tClassify = 0, tGrad = 0, tWNGIR = 0, tMoveTrim = 0;
    Real tElasticity = 0, tHilbert = 0, tRedistance = 0, tAdvect = 0;
    Real tRepair = 0, tWrite = 0;

    std::vector<Index> interfaceFacets;
    std::size_t nInside = cachedInsideCount;
    Real areaInterior = cachedAreaInterior;
    const bool doClassify =
      !haveClassification || classifyEvery == 0 || (it % classifyEvery == 0);
    if (doClassify)
    {
      for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
        mesh.setAttribute({D, c}, Attribute{0});
      for (auto fit = mesh.getFace(); fit; ++fit)
        mesh.setAttribute({D - 1, fit->getIndex()}, Attribute{0});
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
      const MinSTCut::Result classified =
        cut.classify(volumes, moments, graphEdges);

      for (std::size_t l = 0; l < classified.labels.size(); ++l)
        mesh.setAttribute({D, localToCell[l]},
                          classified.labels[l] == MinSTCut::Inside
                            ? interiorAttribute : exteriorAttribute);

      {
        std::vector<char> isSupportCell(mesh.getCellCount(), 0);
        for (auto bit = mesh.getBoundary(); bit; ++bit)
        {
          if (mesh.getAttribute(D - 1, bit->getIndex()) != Attribute{GammaD})
            continue;
          for (const Index c :
               mesh.getConnectivity().getIncidence({1, 2}, bit->getIndex()))
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
          if (!ra || *ra != interiorAttribute) continue;
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

      std::vector<char> currentInterior(mesh.getCellCount(), 0);
      nInside = 0;
      areaInterior = 0;
      for (std::size_t l = 0; l < classified.labels.size(); ++l)
      {
        const auto a = mesh.getAttribute(D, localToCell[l]);
        if (a && *a == interiorAttribute)
        {
          ++nInside;
          areaInterior += cellMoments[l].area;
          currentInterior[localToCell[l]] = 1;
        }
      }
      std::size_t changedCells = 0;
      if (haveClassification)
        for (std::size_t c = 0; c < currentInterior.size(); ++c)
          if (currentInterior[c] != previousInterior[c]) ++changedCells;
      std::vector<Index> oldFacets = cachedInterfaceFacets;
      std::vector<Index> newFacets = interfaceFacets;
      std::sort(oldFacets.begin(), oldFacets.end());
      std::sort(newFacets.begin(), newFacets.end());
      std::vector<Index> changedFacets;
      std::set_symmetric_difference(
          oldFacets.begin(), oldFacets.end(),
          newFacets.begin(), newFacets.end(),
          std::back_inserter(changedFacets));
      previousInterior = std::move(currentInterior);
      cachedInterfaceFacets = interfaceFacets;
      cachedAreaInterior = areaInterior;
      cachedInsideCount = nInside;
      haveClassification = true;
      std::cout << "  classify: interior cells=" << nInside
                << "  interface facets=" << interfaceFacets.size()
                << "  area=" << std::fixed << std::setprecision(4)
                << areaInterior
                << "  changedCells=" << changedCells
                << "  changedInterface=" << changedFacets.size() << '\n';
    }
    else
    {
      interfaceFacets = cachedInterfaceFacets;
      std::cout << "  classify: skipped"
                << "  interior cells=" << nInside
                << "  interface facets=" << interfaceFacets.size()
                << "  area=" << std::fixed << std::setprecision(4)
                << areaInterior << '\n';
    }
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
    WNGIRReport rep;
    {
      rep = wngir.solve(mesh, interfaceFacets, phiFn, gradPhiFn);
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
                << "  tolRMS/h=" << rep.effectiveRMSOverHTol
                << "  tolSup/h=" << rep.effectiveSupOverHTol
                << "  nJumpRMS=" << rep.normalJumpRMS
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
#ifdef RODIN_WNGIR_P2_DISPLACEMENT
    trimmed.getConnectivity().compute(2, 1);
    trimmed.getConnectivity().compute(1, 0);
    trimmed.getConnectivity().compute(2, 2);
#endif
    tMoveTrim = iterTimer.reset();

#ifdef RODIN_WNGIR_P2_DISPLACEMENT
    using StateVectorFES = H1<2, Math::SpatialVector<Real>, WNGIRMesh>;
    StateVectorFES vhInt(std::integral_constant<std::size_t, 2>{}, trimmed, 2);
#else
    using StateVectorFES = P1<Math::SpatialVector<Real>, WNGIRMesh>;
    StateVectorFES vhInt(trimmed, 2);
#endif
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
    auto shapeDensity = Dot(Ae, e) - ell;
    if (trace)
    {
      Real minFaceDist = std::numeric_limits<Real>::infinity();
      Real maxFaceDist = Real(0);
      Real maxAbsRho = Real(0);
      Real maxAbsWeightedRho = Real(0);
      Index minDistFacet = -1;
      Index maxRhoFacet = -1;
      Index maxWeightedFacet = -1;
      std::size_t faceCount = 0;
      for (auto fit = moved.getFace(); fit; ++fit)
      {
        const auto attr = fit->getAttribute();
        if (!attr || *attr != Gamma)
          continue;
        ++faceCount;
        const auto& qf =
          QF::PolytopeQuadratureFormula::get(8, fit->getGeometry());
        const auto& q = fit->getQuadrature(qf);
        for (std::size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          const IntegrationPoint ip(p, &qf, qp);
          const Real dist = p.getDistortion();
          const Real rho = shapeDensity(ip);
          const Real absRho = std::abs(rho);
          const Real absWeightedRho = absRho * dist;
          if (dist < minFaceDist)
          {
            minFaceDist = dist;
            minDistFacet = fit->getIndex();
          }
          if (dist > maxFaceDist)
            maxFaceDist = dist;
          if (absRho > maxAbsRho)
          {
            maxAbsRho = absRho;
            maxRhoFacet = fit->getIndex();
          }
          if (absWeightedRho > maxAbsWeightedRho)
          {
            maxAbsWeightedRho = absWeightedRho;
            maxWeightedFacet = fit->getIndex();
          }
        }
      }
      std::cout << "  hilbertFaceDiag: faces=" << faceCount
                << " minDist=" << std::scientific << minFaceDist
                << "@f" << minDistFacet
                << " maxDist=" << maxFaceDist
                << " max|rho|=" << maxAbsRho
                << "@f" << maxRhoFacet
                << " max|rho|ds=" << maxAbsWeightedRho
                << "@f" << maxWeightedFacet << '\n';
#ifdef RODIN_WNGIR_P2_DISPLACEMENT
      LinearForm faceRhs(w);
      faceRhs -= FaceIntegral(shapeDensity, Dot(nrm, w)).over(Gamma);
      faceRhs.assemble();

      const auto assembleManualFaceRhs =
        [&](std::size_t order)
      {
        Math::Vector<Real> manual(faceRhs.getVector().size());
        manual.setZero();
        const std::size_t vdim = vhMoved.getVectorDimension();
        for (auto fit = moved.getFace(); fit; ++fit)
        {
          const auto attr = fit->getAttribute();
          if (!attr || *attr != Gamma)
            continue;
          const auto& qf =
            QF::PolytopeQuadratureFormula::get(order, fit->getGeometry());
          const auto& q = fit->getQuadrature(qf);
          const Variational::RealH1Element<2> scalarFE(fit->getGeometry());
          const auto& tab = scalarFE.getTabulation(qf);
          const auto& dofs = vhMoved.getDOFs(fit->getDimension(), fit->getIndex());
          for (std::size_t qp = 0; qp < q.getSize(); ++qp)
          {
            const auto& p = q.getPoint(qp);
            const IntegrationPoint ip(p, &qf, qp);
            const Real wdet = qf.getWeight(qp) * p.getDistortion();
            const Real rho = shapeDensity(ip);
            const Vec2 normal = nrm(ip);
            for (std::size_t l = 0; l < static_cast<std::size_t>(dofs.size()); ++l)
            {
              const std::size_t scalarLocal = l / vdim;
              const std::size_t comp = l % vdim;
              manual(dofs(l)) -= wdet * rho * normal(comp)
                * tab.getBasis(qp, scalarLocal);
            }
          }
        }
        return manual;
      };
      Eigen::Index iForm = 0;
      const Real maxForm = faceRhs.getVector().cwiseAbs().maxCoeff(&iForm);
      std::cout << "  faceRhsCheck: maxForm=" << maxForm
                << "@dof" << iForm;
      for (const std::size_t order : {std::size_t(2), std::size_t(4), std::size_t(8)})
      {
        const Math::Vector<Real> manual = assembleManualFaceRhs(order);
        const Math::Vector<Real> diff = faceRhs.getVector() - manual;
        const Real manualNorm = std::max(manual.norm(), Real(1e-30));
        Eigen::Index iManual = 0, iDiff = 0;
        const Real maxManual = manual.cwiseAbs().maxCoeff(&iManual);
        const Real maxDiff = diff.cwiseAbs().maxCoeff(&iDiff);
        std::cout << " q" << order
                  << " relL2=" << diff.norm() / manualNorm
                  << " maxManual=" << maxManual << "@dof" << iManual
                  << " maxDiff=" << maxDiff << "@dof" << iDiff;
      }
      std::cout << '\n';
#endif
    }
    Problem hilbert(g, w);
    if (shapeDerivativeMode == "volume")
    {
      auto eshelby =
        shapeDensity * IdentityMatrix(2)
        - Real(2) * jac.T() * Ae;
      hilbert = Integral(alphaReg * alphaReg * Jacobian(g), Jacobian(w))
              + Integral(g, w)
              - Integral(eshelby, Jacobian(w)).over(interiorAttribute)
              + DirichletBC(g, VectorFunction{ Real(0), Real(0) }).on(GammaD)
              + DirichletBC(g, VectorFunction{ Real(0), Real(0) }).on(GammaN);
    }
    else
    {
      hilbert = Integral(alphaReg * alphaReg * Jacobian(g), Jacobian(w))
              + Integral(g, w)
              - FaceIntegral(shapeDensity, Dot(nrm, w)).over(Gamma)
              + DirichletBC(g, VectorFunction{ Real(0), Real(0) }).on(GammaD)
              + DirichletBC(g, VectorFunction{ Real(0), Real(0) }).on(GammaN);
    }
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
      fittedGrid.clear();
      fittedGrid.setMesh(moved, IO::XDMF::MeshPolicy::Transient);
      fittedGrid.add("dJ_fit", g.getSolution(), IO::XDMF::Center::Node);
      stateGrid.clear();
      stateGrid.setMesh(trimmed, IO::XDMF::MeshPolicy::Transient);
      stateGrid.add("u", us.getSolution(), IO::XDMF::Center::Node);
      xdmf.write(static_cast<Real>(it)).flush();
    }
    tWrite = iterTimer.reset();

    // ---- Stage 5: optional redistance --------------------------------------
    const Real redistanceEta = std::max(rep.sigma, Real(3) * h);
    const Real eikonalDefect =
      weightedEikonalDefect(phiH, redistanceEta);
    const bool periodicRedistance = redistanceEvery > 0
      && ((it + 1) % redistanceEvery == 0 || it + 1 == maxIt);
    const bool adaptiveRedistanceNow =
      adaptiveRedistance && eikonalDefect > redistanceEikonalTol;
    const bool doRedistance = redistanceMode != "none"
      && (periodicRedistance || adaptiveRedistanceNow);
    if (!doRedistance)
    {
      dist.getData() = phiH.getData();
      std::cout << "  redistance: skipped"
                << " eikWelsch=" << std::scientific << eikonalDefect
                << " eta/h=" << (h > Real(0) ? redistanceEta / h : Real(0))
                << '\n';
    }
    else
    {
      ScalarP1 shP1(mesh);
      GridFunction distP1(shP1);
      distP1 = Real(0);
      Distance::Eikonal<ScalarP1, Math::Vector<Real>>(distP1)
        .setInterior(interiorAttribute)
        .setInterface(Gamma)
        .solve()
        .sign();
      dist = [&](const Geometry::Point& p) -> Real
      {
        return distP1.getValue(p);
      };
      const Math::Vector<Real> fmmDist = dist.getData();
      if (redistanceMode == "cp")
      {
        std::vector<std::array<Vec2, 2>> fittedSeg;
        std::vector<Vec2> fittedN;
        fittedSeg.reserve(interfaceFacets.size());
        fittedN.reserve(interfaceFacets.size());
        for (const Index f : interfaceFacets)
        {
          const auto& vv = moved.getFace(f)->getVertices();
          const Vec2 a = moved.getVertexCoordinates(vv[0]);
          const Vec2 b = moved.getVertexCoordinates(vv[1]);
          Vec2 nrm = vec2(b(1) - a(1), -(b(0) - a(0)));
          const Real nn = nrm.norm();
          nrm = nn > Real(1e-30) ? nrm / nn : vec2(1, 0);
          const auto& inc = mesh.getConnectivity().getIncidence({1, 2}, f);
          Vec2 cInt = vec2(0, 0);
          bool haveInt = false;
          for (const Index c : inc)
          {
            const auto at = mesh.getAttribute(D, c);
            if (at && *at == interiorAttribute)
            {
              const auto& cvv = mesh.getCell(c)->getVertices();
              Vec2 ctr = vec2(0, 0);
              for (const Index cv : cvv) ctr += moved.getVertexCoordinates(cv);
              cInt = ctr / static_cast<Real>(cvv.size());
              haveInt = true;
              break;
            }
          }
          const Vec2 mid = Real(0.5) * (a + b);
          if (haveInt && (mid - cInt).dot(nrm) < Real(0)) nrm = -nrm;
          fittedSeg.push_back({ a, b });
          fittedN.push_back(nrm);
        }
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
              ? Real(1) : Real(-1);
            dist.getData()(d) = s * dmin;
          }
        }
      }
      else if (redistanceMode != "fmm")
      {
        std::cerr << "Unsupported redistance mode: " << redistanceMode << '\n';
        return 2;
      }
      std::cout << "  redistance: mode=" << redistanceMode
                << " reason="
                << (periodicRedistance && adaptiveRedistanceNow ? "periodic+adaptive"
                    : periodicRedistance ? "periodic" : "adaptive")
                << " eikWelsch=" << std::scientific << eikonalDefect
                << " eta/h=" << (h > Real(0) ? redistanceEta / h : Real(0))
                << '\n';
    }
    tRedistance = iterTimer.reset();

    // Normalize the descent velocity, advect distance by one step.
    speed = Frobenius(dJ);
    const Real vmax = std::max(speed.max(), Real(1e-30));
    const Math::Vector<Real> phiBeforeAdvection = phiH.getData();
    dJ /= vmax;
    Advection::Lagrangian(adv, advTest, dist, dJ).step(dtCur);
    phiH.getData() = adv.getSolution().getData();
    Real maxPhiChange = 0;
    for (Eigen::Index i = 0; i < phiH.getData().size(); ++i)
      maxPhiChange =
        std::max(maxPhiChange, std::abs(phiH.getData()(i) - phiBeforeAdvection(i)));
    std::cout << "  advect: dt/h=" << (h > Real(0) ? dtCur / h : Real(0))
              << " rawSpeedMax=" << vmax
              << " max|dphi|/h=" << (h > Real(0) ? maxPhiChange / h : Real(0))
              << '\n';
    tAdvect = iterTimer.reset();

    const bool doRepair =
      repairEvery > 0 && ((it + 1) % repairEvery == 0 || it + 1 == maxIt);
    if (doRepair)
    {
      const Math::Vector<Real> rhs = repairMass.getOperator() * phiH.getData();
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
