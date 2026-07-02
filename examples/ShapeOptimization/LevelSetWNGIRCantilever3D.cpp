/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//
// 3D cantilever compliance minimization on a fixed tetrahedral background
// mesh. A graph cut classifies the carried level set, WNGIR fits the cut
// skeleton to phi_h = 0, and elasticity is solved on the fitted material
// submesh. The shape derivative is extended to the background by an H1 metric
// and the level set is advected there.
#include <Rodin/Adaptation.h>
#include <Rodin/Advection/Lagrangian.h>
#include <Rodin/Assembly.h>
#include <Rodin/Distance/Eikonal.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Math.h>
#include <Rodin/MMG.h>
#include <Rodin/Solver/CG.h>
#include <Rodin/Solid.h>
#include <Rodin/Variational.h>

#include "../WNGIRExampleParameters.h"

#include <Eigen/IterativeLinearSolvers>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

using namespace Rodin;
using namespace Rodin::Adaptation;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  using Vec3 = Math::SpatialVector<Real>;
  using Mat3 = Math::SpatialMatrix<Real>;
  using WNGIRMesh = Geometry::Mesh<Context::Local>;

  constexpr Attribute Interior = 1;
  constexpr Attribute Exterior = 2;
  constexpr Attribute Gamma0 = 1;
  constexpr Attribute GammaD = 2;
  constexpr Attribute GammaN = 3;
  constexpr Attribute Gamma = 4;

  constexpr Real muLame = 0.3846;
  constexpr Real lambdaLame = 0.5769;

  Vec3 vec3(Real x = 0, Real y = 0, Real z = 0)
  {
    Vec3 out(3);
    out(0) = x;
    out(1) = y;
    out(2) = z;
    return out;
  }

  Vec3 cross3(const Vec3& a, const Vec3& b)
  {
    return vec3(a(1) * b(2) - a(2) * b(1),
                a(2) * b(0) - a(0) * b(2),
                a(0) * b(1) - a(1) * b(0));
  }

  Real det3(const Vec3& a, const Vec3& b, const Vec3& c)
  {
    return a.dot(cross3(b, c));
  }

  Mat3 edgeMatrix(const std::array<Vec3, 4>& x)
  {
    Mat3 A(3, 3);
    A.col(0) = x[1] - x[0];
    A.col(1) = x[2] - x[0];
    A.col(2) = x[3] - x[0];
    return A;
  }

  Real triangleArea(const WNGIRMesh& mesh, Index facet)
  {
    const auto& v = mesh.getFace(facet)->getVertices();
    const Vec3 a = mesh.getVertexCoordinates(v[0]);
    const Vec3 b = mesh.getVertexCoordinates(v[1]);
    const Vec3 c = mesh.getVertexCoordinates(v[2]);
    return Real(0.5) * cross3(b - a, c - a).norm();
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

  // Squared distance to a triangle, following the Voronoi-region test of
  // Christer Ericson, Real-Time Collision Detection, Section 5.1.5.
  Real pointTriangleDistance(
      const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c)
  {
    const Vec3 ab = b - a, ac = c - a, ap = p - a;
    const Real d1 = ab.dot(ap), d2 = ac.dot(ap);
    if (d1 <= 0 && d2 <= 0) return ap.norm();

    const Vec3 bp = p - b;
    const Real d3 = ab.dot(bp), d4 = ac.dot(bp);
    if (d3 >= 0 && d4 <= d3) return bp.norm();

    const Real vc = d1 * d4 - d3 * d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0)
      return (p - (a + (d1 / (d1 - d3)) * ab)).norm();

    const Vec3 cp = p - c;
    const Real d5 = ab.dot(cp), d6 = ac.dot(cp);
    if (d6 >= 0 && d5 <= d6) return cp.norm();

    const Real vb = d5 * d2 - d1 * d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0)
      return (p - (a + (d2 / (d2 - d6)) * ac)).norm();

    const Real va = d3 * d6 - d5 * d4;
    if (va <= 0 && d4 - d3 >= 0 && d5 - d6 >= 0)
      return (p - (b + ((d4 - d3) / ((d4 - d3) + (d5 - d6))) * (c - b))).norm();

    const Real denom = Real(1) / (va + vb + vc);
    return (p - (a + vb * denom * ab + vc * denom * ac)).norm();
  }

  constexpr std::array<std::array<Real, 4>, 4> TetraQuadrature = {{
    {{ 0.5854101966249685, 0.1381966011250105,
       0.1381966011250105, 0.1381966011250105 }},
    {{ 0.1381966011250105, 0.5854101966249685,
       0.1381966011250105, 0.1381966011250105 }},
    {{ 0.1381966011250105, 0.1381966011250105,
       0.5854101966249685, 0.1381966011250105 }},
    {{ 0.1381966011250105, 0.1381966011250105,
       0.1381966011250105, 0.5854101966249685 }}
  }};

  struct CellMomentInfo
  {
    Index index = 0;
    Real volume = 0;
    Real moment = 0;
    std::array<Vec3, 4> x;
  };

  template <class Phi>
  std::vector<CellMomentInfo> collectCellMoments(
      const WNGIRMesh& mesh, const Phi& phi, Real epsilon)
  {
    std::vector<CellMomentInfo> cells;
    cells.reserve(mesh.getCellCount());
    for (auto it = mesh.getCell(); it; ++it)
    {
      CellMomentInfo info;
      info.index = it->getIndex();
      const auto& v = it->getVertices();
      for (std::size_t i = 0; i < 4; ++i)
        info.x[i] = mesh.getVertexCoordinates(v[i]);
      info.volume = std::abs(det3(info.x[1] - info.x[0],
                                  info.x[2] - info.x[0],
                                  info.x[3] - info.x[0])) / Real(6);
      for (const auto& l : TetraQuadrature)
      {
        Math::SpatialPoint rc(3);
        rc(0) = l[1];
        rc(1) = l[2];
        rc(2) = l[3];
        info.moment += std::tanh(phi.getValue(Geometry::Point(*it, rc)) / epsilon);
      }
      info.moment /= TetraQuadrature.size();
      cells.push_back(std::move(info));
    }
    return cells;
  }

  template <class Displacement>
  void updateMovedMesh(
      const WNGIRMesh& mesh, WNGIRMesh& moved, const Displacement& u)
  {
    const auto& fes = u.getFiniteElementSpace();
    for (Index v = 0; v < mesh.getVertexCount(); ++v)
    {
      const Vec3 x = mesh.getVertexCoordinates(v);
      const auto& d = fes.getDOFs(0, v);
      moved.setVertexCoordinates(v, x + vec3(u[d[0]], u[d[1]], u[d[2]]));
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
              u.getData()(fes.getGlobalIndex({d, polytope.getIndex()}, a * sdim + c));
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

  std::size_t sizeOption(
      int argc, char** argv, const std::string& name, std::size_t fallback)
  {
    const std::string prefix = "--" + name + "=";
    for (int i = 1; i < argc; ++i)
      if (std::string(argv[i]).rfind(prefix, 0) == 0)
        return std::stoul(std::string(argv[i]).substr(prefix.size()));
    return fallback;
  }

  Real realOption(
      int argc, char** argv, const std::string& name, Real fallback)
  {
    const std::string prefix = "--" + name + "=";
    for (int i = 1; i < argc; ++i)
      if (std::string(argv[i]).rfind(prefix, 0) == 0)
        return std::stod(std::string(argv[i]).substr(prefix.size()));
    return fallback;
  }

  std::string stringOption(
      int argc, char** argv, const std::string& name, std::string fallback)
  {
    const std::string prefix = "--" + name + "=";
    for (int i = 1; i < argc; ++i)
      if (std::string(argv[i]).rfind(prefix, 0) == 0)
        return std::string(argv[i]).substr(prefix.size());
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

int run(int argc, char** argv)
{
  const Real L = realOption(argc, argv, "L", Real(2));
  const Real H = realOption(argc, argv, "H", Real(1));
  const Real W = realOption(argc, argv, "W", Real(1));
  const std::size_t n = sizeOption(argc, argv, "n", 20);
  const std::size_t maxIt = sizeOption(argc, argv, "iters", 200);
  const Real h = H / static_cast<Real>(n);
  const Real hmin = realOption(argc, argv, "hmin", Real(0.1) * h);
  const bool initialMMG = sizeOption(argc, argv, "initial-mmg", 0) != 0;
  const Real initialMMGHMin =
    realOption(argc, argv, "initial-mmg-hmin", hmin);
  const Real initialMMGHMax =
    realOption(argc, argv, "initial-mmg-hmax", h);
  const Real initialMMGHausd =
    realOption(argc, argv, "initial-mmg-hausd", Real(0));
  const Real ell = realOption(argc, argv, "ell", Real(0.5));
  const Real alphaReg = realOption(argc, argv, "alpha", Real(2) * h);
  const std::size_t holesX = sizeOption(argc, argv, "holes-x", 6);
  const std::size_t holesY = sizeOption(argc, argv, "holes-y", 2);
  const std::size_t holesZ = sizeOption(argc, argv, "holes-z", 2);
  const Real holeRadiusFactor =
    realOption(argc, argv, "hole-radius-factor", Real(0.12));
  const bool boundaryHoles = sizeOption(argc, argv, "boundary-holes", 1) != 0;
  const Real epsilon = realOption(argc, argv, "classifier-eps", Real(1.25) * h);
  const Real lambdaC = realOption(argc, argv, "classifier-lambda", Real(0.004));
  const Real dt0 = realOption(argc, argv, "dt", Real(2) * h);
  const Real dtMax = realOption(argc, argv, "dt-max", Real(4) * h);
  const bool objectiveGuard = sizeOption(argc, argv, "objective-linesearch", 1) != 0;
  const Real objectiveTol = realOption(argc, argv, "objective-decrease-tol", 1e-10);
  const std::string shapeDerivativeMode =
    stringOption(argc, argv, "shape-derivative", "volume");
  const std::size_t classifyEvery = sizeOption(argc, argv, "classify-every", 1);
  const std::string defaultRedistanceMode =
#ifdef RODIN_WNGIR_P2_DISPLACEMENT
    "cp";
#else
    "fmm";
#endif
  const std::string redistanceMode =
    stringOption(argc, argv, "redistance-mode", defaultRedistanceMode);
  const std::size_t redistanceEvery = sizeOption(argc, argv, "redistance-every", 0);
  const bool adaptiveRedistance =
    sizeOption(argc, argv, "redistance-adaptive", 1) != 0;
  const Real redistanceEikonalTol =
    realOption(argc, argv, "redistance-eikonal-tol", Real(0.35));
  const std::size_t repairEvery = sizeOption(argc, argv, "repair-every", 0);
  const Real repairEta = realOption(argc, argv, "repair-eta", Real(0.5) * h);
  const std::size_t outputEvery = sizeOption(argc, argv, "output-every", 1);
  const bool printTiming = sizeOption(argc, argv, "timing", 0) != 0;
  const bool trace = hasFlag(argc, argv, "trace");

  const std::size_t nx = static_cast<std::size_t>(std::lround(L / h)) + 1;
  const std::size_t ny = n + 1;
  const std::size_t nz = static_cast<std::size_t>(std::lround(W / h)) + 1;
  WNGIRMesh mesh = WNGIRMesh::UniformGrid(
      Polytope::Type::Tetrahedron, { nx, ny, nz });
  mesh.scale(h);
  mesh.getConnectivity().compute(2, 3);
  mesh.getConnectivity().compute(3, 2);
  mesh.getConnectivity().compute(3, 3);
  mesh.getConnectivity().compute(0, 0);
  const std::size_t D = mesh.getDimension();

  auto tagBoundary = [&](WNGIRMesh& m)
  {
    for (auto it = m.getBoundary(); it; ++it)
    {
      Vec3 c = vec3();
      for (const Index v : it->getVertices())
        c += m.getVertexCoordinates(v);
      c /= static_cast<Real>(it->getVertices().size());
      Attribute attr = Gamma0;
      if (c(0) < Real(1e-9))
        attr = GammaD;
      else if (c(0) > L - Real(1e-9)
               && std::abs(c(1) - H / 2) < Real(0.15) * H
               && std::abs(c(2) - W / 2) < Real(0.15) * W)
        attr = GammaN;
      m.setAttribute({D - 1, it->getIndex()}, attr);
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
    mesh.getConnectivity().compute(2, 3);
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(3, 3);
    mesh.getConnectivity().compute(0, 0);
    tagBoundary(mesh);
  }
#ifdef RODIN_WNGIR_P2_DISPLACEMENT
  // Affine copy of the background grid, taken before the parametric P2
  // transformations are installed. Lagrangian advection locates a
  // characteristic foot point per DOF; on parametric cells each location is
  // an inverse-map Newton solve, on this copy it is a single affine inverse.
  WNGIRMesh advectionMesh(mesh);
  mesh.getConnectivity().compute(1, 0);
  mesh.getConnectivity().compute(3, 1);
  mesh.getConnectivity().compute(2, 1);
  installP2Transformations(mesh);
#endif
  auto& conn = mesh.getConnectivity();

#ifdef RODIN_WNGIR_P2_DISPLACEMENT
  using ScalarFES = H1<2, Real, WNGIRMesh>;
  using VectorFES = H1<2, Math::SpatialVector<Real>, WNGIRMesh>;
#else
  using ScalarFES = P1<Real, WNGIRMesh>;
  using VectorFES = P1<Math::SpatialVector<Real>, WNGIRMesh>;
#endif
  using ScalarP1 = P1<Real, WNGIRMesh>;
  using ScalarP0 = P0<Real, WNGIRMesh>;
  using VectorP1 = P1<Math::SpatialVector<Real>, WNGIRMesh>;
  ScalarFES sh(
#ifdef RODIN_WNGIR_P2_DISPLACEMENT
      std::integral_constant<std::size_t, 2>{},
#endif
      mesh);
  ScalarP0 p0(mesh);
  VectorFES vh(
#ifdef RODIN_WNGIR_P2_DISPLACEMENT
      std::integral_constant<std::size_t, 2>{},
#endif
      mesh, 3);
  GridFunction phiH(sh); phiH.setName("phi");
  TrialFunction wngirTrial(vh);
  TestFunction wngirTest(vh);
  auto& u = wngirTrial.getSolution();
  u.setName("fit");
  GridFunction dJ(vh); dJ.setName("dJ");

  // A solid box perforated by a regular array of spherical voids.
  phiH = [&](const Geometry::Point& p) -> Real
  {
    const auto& x = p.getCoordinates();
    Real value = -std::numeric_limits<Real>::infinity();
    const Real radius = holeRadiusFactor * std::min(H, W);
    auto addVoid = [&](const Vec3& c)
    {
      value = std::max(value, radius - (x - c).norm());
    };
    for (std::size_t i = 1; i <= holesX; ++i)
    {
      const Real cx = L * static_cast<Real>(i) / static_cast<Real>(holesX + 1);
      for (std::size_t j = 1; j <= holesY; ++j)
      {
        const Real cy = H * static_cast<Real>(j) / static_cast<Real>(holesY + 1);
        for (std::size_t k = 1; k <= holesZ; ++k)
        {
          const Real cz = W * static_cast<Real>(k) / static_cast<Real>(holesZ + 1);
          addVoid(vec3(cx, cy, cz));
        }
      }
      if (boundaryHoles)
      {
        for (std::size_t k = 1; k <= holesZ; ++k)
        {
          const Real cz = W * static_cast<Real>(k) / static_cast<Real>(holesZ + 1);
          addVoid(vec3(cx, Real(0), cz));
          addVoid(vec3(cx, H, cz));
        }
        for (std::size_t j = 1; j <= holesY; ++j)
        {
          const Real cy = H * static_cast<Real>(j) / static_cast<Real>(holesY + 1);
          addVoid(vec3(cx, cy, Real(0)));
          addVoid(vec3(cx, cy, W));
        }
      }
    }
    return value; // negative material, positive void
  };

  VectorFES gradFes(
#ifdef RODIN_WNGIR_P2_DISPLACEMENT
      std::integral_constant<std::size_t, 2>{},
#endif
      mesh, 3);
  GridFunction gradPhi(gradFes);
  TrialFunction gp(gradFes);
  TestFunction gq(gradFes);
  BilinearForm gradMass(gp, gq);
  gradMass = Integral(gp, gq);
  gradMass.assemble();
  Problem gradProjection(gp, gq);
  gradProjection = gradMass - Integral(Grad(phiH), gq);
  Solver::CG gradSolver(gradProjection);

  TrialFunction sr(sh);
  TestFunction st(sh);
  BilinearForm mass(sr, st);
  BilinearForm stiffness(sr, st);
  Eigen::ConjugateGradient<Math::SparseMatrix<Real>, Eigen::Lower | Eigen::Upper>
    repairCG;
  repairCG.setTolerance(Real(1e-10));
  repairCG.setMaxIterations(static_cast<int>(phiH.getData().size()));
  if (repairEvery > 0)
  {
    mass = Integral(sr, st);
    mass.assemble();
    stiffness = Integral(Grad(sr), Grad(st));
    stiffness.assemble();
    Math::SparseMatrix<Real> repairOperator = mass.getOperator();
    repairOperator += repairEta * repairEta * stiffness.getOperator();
    repairCG.compute(repairOperator);
  }

  GridFunction dist(sh), speed(sh);
#ifdef RODIN_WNGIR_P2_DISPLACEMENT
  // Advection state lives on the affine background copy: transport the P1
  // vertex values there, then recover the mid-edge DOFs of phiH by
  // interpolating the advected field (same pattern as the redistance stage).
  ScalarP1 advFes(advectionMesh);
  VectorP1 advVelFes(advectionMesh, 3);
  ScalarP1 vertexFes(mesh);
  GridFunction advDist(advFes);
  GridFunction advVel(advVelFes);
  GridFunction phiAdvected(vertexFes);
  TrialFunction adv(advFes);
  TestFunction advTest(advFes);
#else
  TrialFunction adv(sh);
  TestFunction advTest(sh);
#endif

  Rodin::Examples::WNGIRExampleDefaults wngirDefaults;
  wngirDefaults.maxIterations = 60;
  wngirDefaults.activeRMSOverHTol = Real(0.12);
  wngirDefaults.activeSupOverHTol = Real(0.65);
  WNGIRParameters wp =
    Rodin::Examples::makeWNGIRParameters(
        argc, argv, h, Gamma, wngirDefaults);
  wp.trace = trace;
  Adaptation::WNGIR wngir(wngirTrial, wngirTest);
  wngir.setParameters(wp);

  WNGIRMesh moved(mesh);
  VectorFES vhMoved(
#ifdef RODIN_WNGIR_P2_DISPLACEMENT
      std::integral_constant<std::size_t, 2>{},
#endif
      moved, 3);
  TrialFunction g(vhMoved);
  TestFunction w(vhMoved);

  IO::XDMF xdmf("LevelSetWNGIRCantilever3D");
  auto domainGrid = xdmf.grid("domain");
  domainGrid.setMesh(mesh, IO::XDMF::MeshPolicy::Transient);
  auto fittedGrid = xdmf.grid("fitted");
  auto stateGrid = xdmf.grid("state");
  std::ofstream objectiveFile("LevelSetWNGIRCantilever3D.obj.txt");

  Math::Vector<Real> phiGood = phiH.getData();
  Real dt = std::min(dt0, dtMax);
  Real lastObjective = std::numeric_limits<Real>::infinity();
  Real complianceCap = std::numeric_limits<Real>::infinity();

  std::cout << "3D WNGIR cantilever: " << mesh.getCellCount() << " tetrahedra"
            << "  grid=" << nx << "x" << ny << "x" << nz
            << "\n  h=" << h << " ell=" << ell << " alpha=" << alphaReg
            << " holes=" << holesX << "x" << holesY << "x" << holesZ
            << " holeR=" << holeRadiusFactor * std::min(H, W)
            << " boundaryHoles=" << boundaryHoles
            << " dt=" << dt << " classify=" << classifyEvery
            << " redistance=" << redistanceMode << "/" << redistanceEvery
            << " adaptive=" << adaptiveRedistance
            << " eikTol=" << redistanceEikonalTol
            << " repair=" << repairEvery
            << "\n  WNGIR metric: gammaM=" << wp.gammaM
            << " gammaDev=" << wp.gammaH
            << " gammaDiv=" << wp.gammaDiv
            << " ellM=" << wp.ellM << '\n';

  auto cellGradientMagnitude = [&](const auto& gf, const Polytope& cell) -> Real
  {
    const auto& vv = cell.getVertices();
    const Vec3 x0 = mesh.getVertexCoordinates(vv[0]);
    const Vec3 x1 = mesh.getVertexCoordinates(vv[1]);
    const Vec3 x2 = mesh.getVertexCoordinates(vv[2]);
    const Vec3 x3 = mesh.getVertexCoordinates(vv[3]);
    const Real p0 = gf.getData()(sh.getGlobalIndex({0, vv[0]}, 0));
    const Real p1 = gf.getData()(sh.getGlobalIndex({0, vv[1]}, 0));
    const Real p2 = gf.getData()(sh.getGlobalIndex({0, vv[2]}, 0));
    const Real p3 = gf.getData()(sh.getGlobalIndex({0, vv[3]}, 0));
    Mat3 J(3, 3);
    J.col(0) = x1 - x0;
    J.col(1) = x2 - x0;
    J.col(2) = x3 - x0;
    const Real det = J.determinant();
    if (std::abs(det) < Real(1e-30))
      return Real(1);
    Vec3 dphi(3);
    dphi(0) = p1 - p0;
    dphi(1) = p2 - p0;
    dphi(2) = p3 - p0;
    return (J.inverse().transpose() * dphi).norm();
  };
  auto weightedEikonalDefect = [&](const auto& gf, Real eta) -> Real
  {
    const Real eta2 = std::max(eta * eta, Real(1e-30));
    Real weightedError = 0;
    Real weightSum = 0;
    for (auto cit = mesh.getCell(); cit; ++cit)
    {
      const auto& vv = cit->getVertices();
      std::array<Vec3, 4> x;
      Real phiAvg = 0;
      for (std::size_t i = 0; i < 4; ++i)
      {
        x[i] = mesh.getVertexCoordinates(vv[i]);
        phiAvg += gf.getData()(sh.getGlobalIndex({0, vv[i]}, 0));
      }
      phiAvg /= Real(4);
      const Real volume = std::abs(det3(x[1] - x[0], x[2] - x[0], x[3] - x[0]))
        / Real(6);
      const Real omega = std::exp(-(phiAvg * phiAvg) / eta2);
      const Real gradMag = cellGradientMagnitude(gf, *cit);
      weightedError += volume * omega * (gradMag - Real(1)) * (gradMag - Real(1));
      weightSum += volume * omega;
    }
    return weightSum > Real(0) ? std::sqrt(weightedError / weightSum) : Real(0);
  };

  std::size_t shapeAttempts = 0;
  std::size_t acceptedShapeIterations = 0;
  bool haveClassification = false;
  std::vector<char> previousInterior(mesh.getCellCount(), 0);
  std::vector<Index> cachedInterfaceFacets;
  Real cachedVolume = 0;
  std::size_t cachedInsideCount = 0;
  while (acceptedShapeIterations < maxIt)
  {
    StageTimer iterTimer;
    Real tClassify = 0;
    Real tGrad = 0;
    Real tWNGIR = 0;
    Real tMoveTrim = 0;
    Real tElasticity = 0;
    Real tHilbert = 0;
    Real tWrite = 0;
    Real tRedistance = 0;
    Real tAdvect = 0;
    Real tRepair = 0;

    const std::size_t itn = acceptedShapeIterations;
    ++shapeAttempts;
    std::cout << "\n--- Iteration " << itn << " ---\n";
    std::vector<Index> interfaceFacets;
    Real volume = cachedVolume;
    std::size_t insideCount = cachedInsideCount;
    const bool doClassify =
      !haveClassification || classifyEvery == 0 || (itn % classifyEvery == 0);
    if (doClassify)
    {
      for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
        mesh.setAttribute({D, c}, Attribute{0});
      for (auto it = mesh.getFace(); it; ++it)
        mesh.setAttribute({D - 1, it->getIndex()}, Attribute{0});
      tagBoundary(mesh);

      const auto cells = collectCellMoments(mesh, phiH, epsilon);
      std::unordered_map<Index, std::size_t> cellToLocal;
      std::vector<Index> localToCell;
      std::vector<Real> volumes(cells.size()), moments(cells.size());
      for (std::size_t i = 0; i < cells.size(); ++i)
      {
        cellToLocal[cells[i].index] = i;
        localToCell.push_back(cells[i].index);
        volumes[i] = cells[i].volume;
        moments[i] = cells[i].moment;
      }
      std::vector<MinSTCut::Edge> edges;
      for (auto fit = mesh.getFace(); fit; ++fit)
      {
        const Index f = fit->getIndex();
        const auto& inc = conn.getIncidence({2, 3}, f);
        if (inc.size() != 2) continue;
        edges.push_back({static_cast<Index>(cellToLocal.at(inc[0])),
                         static_cast<Index>(cellToLocal.at(inc[1])),
                         lambdaC * triangleArea(mesh, f), f});
      }
      const auto classified = MinSTCut().classify(volumes, moments, edges);
      for (std::size_t i = 0; i < classified.labels.size(); ++i)
        mesh.setAttribute({D, localToCell[i]},
            classified.labels[i] == MinSTCut::Inside ? Interior : Exterior);

      std::vector<char> support(mesh.getCellCount(), 0);
      for (auto bit = mesh.getBoundary(); bit; ++bit)
        if (mesh.getAttribute(D - 1, bit->getIndex()) == Attribute{GammaD})
          for (const Index c : conn.getIncidence({2, 3}, bit->getIndex()))
            support[c] = 1;
      const auto materialComponents = mesh.ccl(
            [](const Polytope& a, const Polytope& b)
            { return a.getAttribute() == b.getAttribute(); });
      for (const auto& component : materialComponents.getComponents())
      {
        if (component.empty()) continue;
        const auto attr = mesh.getAttribute(D, *component.begin());
        if (!attr || *attr != Interior) continue;
        bool anchored = false;
        for (const Index c : component) anchored = anchored || support[c];
        if (!anchored)
          for (const Index c : component) mesh.setAttribute({D, c}, Exterior);
      }

      for (auto fit = mesh.getFace(); fit; ++fit)
      {
        const Index f = fit->getIndex();
        const auto& inc = conn.getIncidence({2, 3}, f);
        if (inc.size() != 2) continue;
        const auto a = mesh.getAttribute(D, inc[0]);
        const auto b = mesh.getAttribute(D, inc[1]);
        if ((a && *a == Interior) != (b && *b == Interior))
        {
          interfaceFacets.push_back(f);
          mesh.setAttribute({D - 1, f}, Gamma);
        }
      }
      std::vector<char> currentInterior(mesh.getCellCount(), 0);
      volume = 0;
      insideCount = 0;
      for (std::size_t i = 0; i < cells.size(); ++i)
        if (mesh.getAttribute(D, localToCell[i]) == Attribute{Interior})
        {
          volume += cells[i].volume;
          ++insideCount;
          currentInterior[localToCell[i]] = 1;
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
      cachedVolume = volume;
      cachedInsideCount = insideCount;
      haveClassification = true;
      std::cout << "  classify: inside=" << insideCount
                << " interface=" << interfaceFacets.size()
                << " volume=" << volume
                << " changedCells=" << changedCells
                << " changedInterface=" << changedFacets.size() << '\n';
    }
    else
    {
      interfaceFacets = cachedInterfaceFacets;
      std::cout << "  classify: skipped"
                << " inside=" << insideCount
                << " interface=" << interfaceFacets.size()
                << " volume=" << volume << '\n';
    }
    if (interfaceFacets.empty())
    {
      std::cerr << "  empty classified interface; stopping\n";
      break;
    }
    tClassify = iterTimer.reset();

    gradProjection.assemble();
    gradProjection.solve(gradSolver);
    gradPhi.getData() = gp.getSolution().getData();
    tGrad = iterTimer.reset();

    RealFunction phiFn([&](const Geometry::Point& p) { return phiH.getValue(p); });
    AnalyticVectorFunction gradFn(
        [&](const Geometry::Point& p) { return gradPhi.getValue(p); }, 3);
    u.getData().setZero();
    const auto report = wngir.solve(mesh, interfaceFacets, phiFn, gradFn);
    std::cout << "  WNGIR: it=" << report.iterations
              << " exit=" << report.exitReason
              << " activeRMS=" << std::scientific << report.activeRMS
              << " activeRMS/h=" << (h > Real(0) ? report.activeRMS / h : Real(0))
              << " activeSup/h=" << (h > Real(0) ? report.activeSup / h : Real(0))
              << " tolRMS/h=" << report.effectiveRMSOverHTol
              << " tolSup/h=" << report.effectiveSupOverHTol
              << " nJumpRMS=" << report.normalJumpRMS
              << '\n';
    tWNGIR = iterTimer.reset();

    updateMovedMesh(mesh, moved, u);
    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
      if (const auto a = mesh.getAttribute(D, c)) moved.setAttribute({D, c}, *a);
    for (auto fit = mesh.getFace(); fit; ++fit)
      if (const auto a = mesh.getAttribute(D - 1, fit->getIndex()))
        moved.setAttribute({D - 1, fit->getIndex()}, *a);

    // Solve mechanics only on the fitted material submesh. Since the active
    // cell set may change, this FE graph is rebuilt for the current design.
    SubMesh<Context::Local> trimmed = moved.trim(Exterior);
    trimmed.getConnectivity().compute(2, 3);
#ifdef RODIN_WNGIR_P2_DISPLACEMENT
    trimmed.getConnectivity().compute(3, 2);
    trimmed.getConnectivity().compute(1, 0);
    trimmed.getConnectivity().compute(3, 1);
    trimmed.getConnectivity().compute(2, 1);
#endif
    tMoveTrim = iterTimer.reset();

#ifdef RODIN_WNGIR_P2_DISPLACEMENT
    using StateVectorFES = H1<2, Math::SpatialVector<Real>, WNGIRMesh>;
    StateVectorFES vhInterior(std::integral_constant<std::size_t, 2>{}, trimmed, 3);
#else
    using StateVectorFES = P1<Math::SpatialVector<Real>, WNGIRMesh>;
    StateVectorFES vhInterior(trimmed, 3);
#endif
    TrialFunction us(vhInterior);
    TestFunction vs(vhInterior);
    Problem elasticity(us, vs);
    elasticity = LinearElasticityIntegral(us, vs)(lambdaLame, muLame)
               - BoundaryIntegral(VectorFunction{0, 0, -1}, vs).over(GammaN)
               + DirichletBC(us, VectorFunction{0, 0, 0}).on(GammaD);
    Solver::CG elasticitySolver(elasticity);
    elasticity.assemble();
    elasticity.solve(elasticitySolver);
    tElasticity = iterTimer.reset();

    const Real compliance =
      elasticity.getLinearSystem().getVector().dot(
          elasticity.getLinearSystem().getSolution());
    const Real objective = compliance + ell * volume;
    const bool invalid = !std::isfinite(compliance) || compliance < 0
                         || compliance > complianceCap;
    const bool increased = objectiveGuard && std::isfinite(lastObjective)
      && objective > lastObjective
                   + objectiveTol * std::max(Real(1), std::abs(lastObjective));
    if (invalid || increased)
    {
      phiH.getData() = phiGood;
      dt *= Real(0.5);
      std::cout << "  REJECT: J=" << objective << " dt -> " << dt << '\n';
      if (dt < Real(1e-3) * h) break;
      continue;
    }
    phiGood = phiH.getData();
    lastObjective = objective;
    dt = std::min(dtMax, Real(1.3) * dt);
    if (itn == 0) complianceCap = Real(50) * compliance;
    objectiveFile << objective << '\n';
    std::cout << "  compliance=" << compliance << " volume=" << volume
              << " J=" << objective << " dt=" << dt << '\n';

    auto jac = Jacobian(us.getSolution());
    jac.traceOf(Interior);
    auto strain = Real(0.5) * (jac + jac.T());
    auto stress = Real(2) * muLame * strain
                + lambdaLame * Trace(strain) * IdentityMatrix(3);
    auto normal = FaceNormal(moved);
    normal.traceOf(Interior);
    auto shapeDensity = Dot(stress, strain) - ell;
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
                << " max|rho|dS=" << maxAbsWeightedRho
                << "@f" << maxWeightedFacet << '\n';
    }
    Problem hilbert(g, w);
    if (shapeDerivativeMode == "volume")
    {
      auto eshelby =
        shapeDensity * IdentityMatrix(3)
        - Real(2) * jac.T() * stress;
      hilbert = Integral(alphaReg * alphaReg * Jacobian(g), Jacobian(w))
              + Integral(g, w)
              - Integral(eshelby, Jacobian(w)).over(Interior)
              + DirichletBC(g, VectorFunction{0, 0, 0}).on(GammaD)
              + DirichletBC(g, VectorFunction{0, 0, 0}).on(GammaN);
    }
    else
    {
      hilbert = Integral(alphaReg * alphaReg * Jacobian(g), Jacobian(w))
              + Integral(g, w)
              - FaceIntegral(shapeDensity, Dot(normal, w)).over(Gamma)
              + DirichletBC(g, VectorFunction{0, 0, 0}).on(GammaD)
              + DirichletBC(g, VectorFunction{0, 0, 0}).on(GammaN);
    }
    Solver::CG hilbertSolver(hilbert);
    hilbert.assemble();
    hilbert.solve(hilbertSolver);
    dJ.getData() = g.getSolution().getData();
    tHilbert = iterTimer.reset();

    // Write before redistance/advection mutates phiH into the next trial
    // carrier. Thus every frame corresponds to a passed shape objective step.
    if ((outputEvery > 0 && itn % outputEvery == 0) || itn + 1 == maxIt)
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
      xdmf.write(static_cast<Real>(itn)).flush();
      domainGrid.clear();
      stateGrid.clear();
    }
    tWrite = iterTimer.reset();

    const Real redistanceEta = std::max(report.sigma, Real(3) * h);
    const Real eikonalDefect =
      weightedEikonalDefect(phiH, redistanceEta);
    const bool periodicRedistance = redistanceEvery > 0
      && ((itn + 1) % redistanceEvery == 0 || itn + 1 == maxIt);
    const bool adaptiveRedistanceNow =
      adaptiveRedistance && eikonalDefect > redistanceEikonalTol;
    const bool doRedistance = redistanceMode != "none"
      && (periodicRedistance || adaptiveRedistanceNow);
    if (doRedistance)
    {
      ScalarP1 shP1(mesh);
      GridFunction distP1(shP1);
      distP1 = Real(0);
      Distance::Eikonal<ScalarP1, Math::Vector<Real>>(distP1)
        .setInterior(Interior).setInterface(Gamma).solve().sign();
      dist = [&](const Geometry::Point& p) -> Real
      {
        return distP1.getValue(p);
      };
      if (redistanceMode == "cp")
      {
        struct FittedFace { std::array<Vec3, 3> x; Vec3 normal; };
        std::vector<FittedFace> fitted;
        fitted.reserve(interfaceFacets.size());
        for (const Index f : interfaceFacets)
        {
          const auto& fv = moved.getFace(f)->getVertices();
          FittedFace ff;
          for (int i = 0; i < 3; ++i) ff.x[i] = moved.getVertexCoordinates(fv[i]);
          ff.normal = cross3(ff.x[1] - ff.x[0], ff.x[2] - ff.x[0]);
          ff.normal /= std::max(ff.normal.norm(), Real(1e-30));
          Vec3 fc = (ff.x[0] + ff.x[1] + ff.x[2]) / Real(3);
          for (const Index c : conn.getIncidence({2, 3}, f))
            if (mesh.getAttribute(D, c) == Attribute{Interior})
            {
              Vec3 cc = vec3();
              for (const Index v : moved.getCell(c)->getVertices())
                cc += moved.getVertexCoordinates(v);
              cc /= Real(4);
              if ((fc - cc).dot(ff.normal) < 0) ff.normal = -ff.normal;
              break;
            }
          fitted.push_back(ff);
        }
        const Real band = Real(6) * h;
        for (Index v = 0; v < mesh.getVertexCount(); ++v)
        {
          const Index dof = sh.getGlobalIndex({0, v}, 0);
          if (std::abs(dist.getData()(dof)) > band) continue;
          const Vec3 x = mesh.getVertexCoordinates(v);
          Real dmin = std::numeric_limits<Real>::infinity();
          std::size_t nearest = 0;
          for (std::size_t i = 0; i < fitted.size(); ++i)
          {
            const Real d = pointTriangleDistance(
                x, fitted[i].x[0], fitted[i].x[1], fitted[i].x[2]);
            if (d < dmin) { dmin = d; nearest = i; }
          }
          const Vec3 fc = (fitted[nearest].x[0] + fitted[nearest].x[1]
                         + fitted[nearest].x[2]) / Real(3);
          dist.getData()(dof) =
            ((x - fc).dot(fitted[nearest].normal) >= 0 ? dmin : -dmin);
        }
      }
      else if (redistanceMode != "fmm")
      {
        std::cerr << "Unsupported 3D redistance mode: " << redistanceMode << '\n';
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
    else
    {
      dist.getData() = phiH.getData();
      std::cout << "  redistance: skipped"
                << " eikWelsch=" << std::scientific << eikonalDefect
                << " eta/h=" << (h > Real(0) ? redistanceEta / h : Real(0))
                << '\n';
    }
    tRedistance = iterTimer.reset();

    speed = Frobenius(dJ);
    const Real rawSpeedMax = std::max(speed.max(), Real(1e-30));
    const Math::Vector<Real> phiBeforeAdvection = phiH.getData();
    dJ /= rawSpeedMax;
#ifdef RODIN_WNGIR_P2_DISPLACEMENT
    for (Index v = 0; v < mesh.getVertexCount(); ++v)
    {
      advDist.getData()(advFes.getGlobalIndex({0, v}, 0)) =
        dist.getData()(sh.getGlobalIndex({0, v}, 0));
      const auto& src = vh.getDOFs(0, v);
      const auto& dst = advVelFes.getDOFs(0, v);
      for (std::size_t c = 0; c < 3; ++c)
        advVel.getData()(dst[c]) = dJ.getData()(src[c]);
    }
    Advection::Lagrangian(adv, advTest, advDist, advVel).step(dt);
    phiAdvected.getData() = adv.getSolution().getData();
    phiH = [&](const Geometry::Point& p) -> Real
    {
      return phiAdvected.getValue(p);
    };
#else
    Advection::Lagrangian(adv, advTest, dist, dJ).step(dt);
    phiH.getData() = adv.getSolution().getData();
#endif
    Real maxPhiChange = 0;
    for (Eigen::Index i = 0; i < phiH.getData().size(); ++i)
      maxPhiChange =
        std::max(maxPhiChange, std::abs(phiH.getData()(i) - phiBeforeAdvection(i)));
    std::cout << "  advect: dt/h=" << (h > Real(0) ? dt / h : Real(0))
              << " rawSpeedMax=" << rawSpeedMax
              << " max|dphi|/h=" << (h > Real(0) ? maxPhiChange / h : Real(0))
              << '\n';
    tAdvect = iterTimer.reset();

    if (repairEvery > 0
        && ((itn + 1) % repairEvery == 0 || itn + 1 == maxIt))
    {
      const Math::Vector<Real> rhs = mass.getOperator() * phiH.getData();
      phiH.getData() = repairCG.solve(rhs);
    }
    tRepair = iterTimer.reset();

    if (printTiming)
    {
      const Real total = tClassify + tGrad + tWNGIR + tMoveTrim
                       + tElasticity + tHilbert + tWrite + tRedistance
                       + tAdvect + tRepair;
      std::cout << "  timing:"
                << " classify=" << std::scientific << std::setprecision(2) << tClassify
                << " grad=" << tGrad
                << " wngir=" << tWNGIR
                << " wngirSetup=" << report.tSetup
                << " wngirBulk=" << report.tBulk
                << " wngirAsm=" << report.tAssembly
                << " wngirSolve=" << report.tSolve
                << " wngirLS=" << report.tLineSearch
                << " moveTrim=" << tMoveTrim
                << " elas=" << tElasticity
                << " hilbert=" << tHilbert
                << " write=" << tWrite
                << " redist=" << tRedistance
                << " advect=" << tAdvect
                << " repair=" << tRepair
                << " total=" << total << '\n';
    }

    ++acceptedShapeIterations;
  }

  std::cout << "\nDone. Accepted shape iterations: "
            << acceptedShapeIterations << " / " << maxIt
            << "  attempts=" << shapeAttempts
            << ". Objective history in LevelSetWNGIRCantilever3D.obj.txt\n";
  return 0;
}

int main(int argc, char** argv)
{
  return run(argc, argv);
}
