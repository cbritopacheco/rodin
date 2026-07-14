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
#include <Rodin/PETSc.h>
#include <Rodin/Solver/CG.h>
#include <Rodin/Solid.h>
#include <Rodin/Variational.h>

#include "../WNGIRExampleParameters.h"

#include <Eigen/IterativeLinearSolvers>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
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

  void copyPETScToEigen(Vec source, Math::Vector<Real>& destination)
  {
    PetscInt n = 0;
    VecGetSize(source, &n);
    destination.resize(n);
    const PetscScalar* values = nullptr;
    VecGetArrayRead(source, &values);
    for (PetscInt i = 0; i < n; ++i)
      destination(i) = PetscRealPart(values[i]);
    VecRestoreArrayRead(source, &values);
  }

  void enableFactorReuse(::KSP ksp)
  {
    PC pc = nullptr;
    KSPGetPC(ksp, &pc);
    const char* pcType = nullptr;
    PCGetType(pc, &pcType);
    if (pcType && (!std::strcmp(pcType, PCLU)
                || !std::strcmp(pcType, PCCHOLESKY)))
    {
      PCFactorSetReuseOrdering(pc, PETSC_TRUE);
      PCFactorSetReuseFill(pc, PETSC_TRUE);
    }
  }

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

  PetscInt kspIterations(::KSP ksp)
  {
    PetscInt iterations = 0;
    KSPGetIterationNumber(ksp, &iterations);
    return iterations;
  }
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
  const Real epsilon = realOption(argc, argv, "classifier-eps", Real(1.25) * h);
  const Real lambdaC = realOption(argc, argv, "classifier-lambda", Real(0.004));
  const Real dt0 = realOption(argc, argv, "dt", Real(2) * h);
  const Real dtMax = realOption(argc, argv, "dt-max", Real(4) * h);
  const bool objectiveGuard = sizeOption(argc, argv, "objective-linesearch", 1) != 0;
  const Real objectiveTol = realOption(argc, argv, "objective-decrease-tol", 1e-10);
  const std::string reinitMode = stringOption(argc, argv, "reinit", "none");
  const std::size_t reinitEvery = sizeOption(argc, argv, "reinit-every", 1);
  const std::string redistanceMode =
    stringOption(argc, argv, "redistance-mode", "none");
  const std::size_t redistanceEvery = sizeOption(argc, argv, "redistance-every", 0);
  const std::size_t repairEvery = sizeOption(argc, argv, "repair-every", 5);
  const Real repairEta = realOption(argc, argv, "repair-eta", Real(0.5) * h);
  const std::size_t outputEvery = sizeOption(argc, argv, "output-every", 5);
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
  auto& conn = mesh.getConnectivity();

  using ScalarP1 = P1<Real, WNGIRMesh>;
  using ScalarP0 = P0<Real, WNGIRMesh>;
  using VectorP1 = P1<Math::SpatialVector<Real>, WNGIRMesh>;
  ScalarP1 sh(mesh);
  ScalarP0 p0(mesh);
  VectorP1 vh(mesh, 3);
  GridFunction phiH(sh); phiH.setName("phi");
  PETSc::Variational::TrialFunction wngirTrial(vh);
  PETSc::Variational::TestFunction wngirTest(vh);
  auto& u = wngirTrial.getSolution();
  u.setName("fit");
  GridFunction dJ(vh); dJ.setName("dJ");

  // A solid box perforated by a regular array of spherical voids.
  phiH = [&](const Geometry::Point& p) -> Real
  {
    const auto& x = p.getCoordinates();
    Real value = -std::numeric_limits<Real>::infinity();
    const Real radius = Real(0.16) * std::min(H, W);
    for (int i = 1; i <= 6; ++i)
      for (int j = 1; j <= 2; ++j)
        for (int k = 1; k <= 2; ++k)
        {
          const Vec3 c = vec3(L * i / Real(7), H * j / Real(3), W * k / Real(3));
          value = std::max(value, radius - (x - c).norm());
        }
    return value; // negative material, positive void
  };

  VectorP1 gradFes(mesh, 3);
  GridFunction gradPhi(gradFes);
  PETSc::Variational::TrialFunction gp(gradFes);
  PETSc::Variational::TestFunction gq(gradFes);
  BilinearForm gradMass(gp, gq);
  gradMass = Integral(gp, gq);
  gradMass.assemble();
  Problem gradProjection(gp, gq);
  gradProjection = gradMass - Integral(Grad(phiH), gq);
  Solver::KSP gradSolver(gradProjection);
  gradSolver.setPrefix("grad_");

  TrialFunction sr(sh);
  TestFunction st(sh);
  BilinearForm mass(sr, st);
  mass = Integral(sr, st);
  mass.assemble();
  BilinearForm stiffness(sr, st);
  stiffness = Integral(Grad(sr), Grad(st));
  stiffness.assemble();
  Math::SparseMatrix<Real> repairOperator = mass.getOperator();
  repairOperator += repairEta * repairEta * stiffness.getOperator();
  Eigen::ConjugateGradient<Math::SparseMatrix<Real>, Eigen::Lower | Eigen::Upper>
    repairCG;
  repairCG.setTolerance(Real(1e-10));
  repairCG.setMaxIterations(static_cast<int>(phiH.getData().size()));
  repairCG.compute(repairOperator);

  GridFunction dist(sh), speed(sh);
  TrialFunction adv(sh);
  TestFunction advTest(sh);

  Rodin::Examples::WNGIRExampleDefaults wngirDefaults;
  wngirDefaults.maxIterations = 60;
  wngirDefaults.activeRMSOverHTol = Real(0.2);
  WNGIRParameters wp =
    Rodin::Examples::makeWNGIRParameters(
        argc, argv, h, Gamma, wngirDefaults);
  wp.trace = trace;
  Adaptation::WNGIR wngir(wngirTrial, wngirTest);
  wngir.setParameters(wp);

  WNGIRMesh moved(mesh);
  VectorP1 vhMoved(moved, 3);
  PETSc::Variational::TrialFunction g(vhMoved);
  PETSc::Variational::TestFunction w(vhMoved);

  IO::XDMF xdmf("LevelSetWNGIRCantilever3D");
  auto domainGrid = xdmf.grid("domain");
  domainGrid.setMesh(mesh, IO::XDMF::MeshPolicy::Transient);
  auto stateGrid = xdmf.grid("state");
  std::ofstream objectiveFile("LevelSetWNGIRCantilever3D.obj.txt");

  Math::Vector<Real> phiGood = phiH.getData();
  Real dt = std::min(dt0, dtMax);
  Real lastObjective = std::numeric_limits<Real>::infinity();
  Real complianceCap = std::numeric_limits<Real>::infinity();

  std::cout << "3D WNGIR cantilever: " << mesh.getCellCount() << " tetrahedra"
            << "  grid=" << nx << "x" << ny << "x" << nz
            << "\n  h=" << h << " ell=" << ell << " alpha=" << alphaReg
            << " dt=" << dt << " reinit=" << reinitMode << "/" << reinitEvery
            << " redistance=" << redistanceMode << "/" << redistanceEvery
            << " repair=" << repairEvery
            << "\n  WNGIR metric: gammaM=" << wp.gammaM
            << " gammaDev=" << wp.gammaH
            << " gammaDiv=" << wp.gammaDiv
            << " ellM=" << wp.ellM << '\n';

  std::size_t shapeAttempts = 0;
  std::size_t acceptedShapeIterations = 0;
  while (acceptedShapeIterations < maxIt)
  {
    const std::size_t itn = acceptedShapeIterations;
    ++shapeAttempts;
    std::cout << "\n--- Iteration " << itn << " ---\n";
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

    // Keep only material components connected to the clamped boundary.
    std::vector<char> support(mesh.getCellCount(), 0);
    for (auto bit = mesh.getBoundary(); bit; ++bit)
      if (mesh.getAttribute(D - 1, bit->getIndex()) == Attribute{GammaD})
        for (const Index c : conn.getIncidence({2, 3}, bit->getIndex()))
          support[c] = 1;
    for (const auto& component : mesh.ccl(
          [](const Polytope& a, const Polytope& b)
          { return a.getAttribute() == b.getAttribute(); }).getComponents())
    {
      if (component.empty()) continue;
      const auto attr = mesh.getAttribute(D, *component.begin());
      if (!attr || *attr != Interior) continue;
      bool anchored = false;
      for (const Index c : component) anchored = anchored || support[c];
      if (!anchored)
        for (const Index c : component) mesh.setAttribute({D, c}, Exterior);
    }

    std::vector<Index> interfaceFacets;
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
    Real volume = 0;
    std::size_t insideCount = 0;
    for (std::size_t i = 0; i < cells.size(); ++i)
      if (mesh.getAttribute(D, localToCell[i]) == Attribute{Interior})
      { volume += cells[i].volume; ++insideCount; }
    std::cout << "  classify: inside=" << insideCount
              << " interface=" << interfaceFacets.size()
              << " volume=" << volume << '\n';
    if (interfaceFacets.empty())
    {
      std::cerr << "  empty classified interface; stopping\n";
      break;
    }

    gradProjection.assemble();
    gradProjection.solve(gradSolver);
    enableFactorReuse(gradSolver.getHandle());
    if (trace)
      std::cout << "  PETSc grad: it="
                << kspIterations(gradSolver.getHandle()) << '\n';
    copyPETScToEigen(gp.getSolution().getData(), gradPhi.getData());
    RealFunction phiFn([&](const Geometry::Point& p) { return phiH.getValue(p); });
    AnalyticVectorFunction gradFn(
        [&](const Geometry::Point& p) { return gradPhi.getValue(p); }, 3);
    u = Real(0);
    const auto report = wngir.solve(mesh, interfaceFacets, phiFn, gradFn);
    std::cout << "  WNGIR: it=" << report.iterations
              << " exit=" << report.exitReason
              << " activeRMS=" << std::scientific << report.activeRMS
              << " activeRMS/h=" << (h > Real(0) ? report.activeRMS / h : Real(0))
              << " activeSup/h=" << (h > Real(0) ? report.activeSup / h : Real(0))
              << '\n';

    updateMovedMesh(mesh, moved, u);
    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
      if (const auto a = mesh.getAttribute(D, c)) moved.setAttribute({D, c}, *a);
    for (auto fit = mesh.getFace(); fit; ++fit)
      if (const auto a = mesh.getAttribute(D - 1, fit->getIndex()))
        moved.setAttribute({D - 1, fit->getIndex()}, *a);

    // Solve mechanics only on the fitted material submesh. Since the active
    // cell set may change, this FE graph and its PETSc symbolic structure are
    // necessarily rebuilt for the current design.
    SubMesh<Context::Local> trimmed = moved.trim(Exterior);
    trimmed.getConnectivity().compute(2, 3);
    VectorP1 vhInterior(trimmed, 3);
    PETSc::Variational::TrialFunction us(vhInterior);
    PETSc::Variational::TestFunction vs(vhInterior);
    Problem elasticity(us, vs);
    elasticity = LinearElasticityIntegral(us, vs)(lambdaLame, muLame)
               - BoundaryIntegral(VectorFunction{0, 0, -1}, vs).over(GammaN)
               + DirichletBC(us, VectorFunction{0, 0, 0}).on(GammaD);
    Solver::KSP elasticitySolver(elasticity);
    elasticitySolver.setPrefix("elasticity_");
    elasticity.assemble();
    elasticity.solve(elasticitySolver);
    enableFactorReuse(elasticitySolver.getHandle());
    if (trace)
      std::cout << "  PETSc elasticity: it="
                << kspIterations(elasticitySolver.getHandle()) << '\n';
    PetscScalar complianceValue = 0;
    VecDot(elasticity.getLinearSystem().getVector(),
           elasticity.getLinearSystem().getSolution(), &complianceValue);
    const Real compliance = PetscRealPart(complianceValue);
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
    Problem hilbert(g, w);
    hilbert = Integral(alphaReg * alphaReg * Jacobian(g), Jacobian(w))
            + Integral(g, w)
            - FaceIntegral(Dot(stress, strain) - ell, Dot(normal, w)).over(Gamma)
            + DirichletBC(g, VectorFunction{0, 0, 0}).on(GammaD)
            + DirichletBC(g, VectorFunction{0, 0, 0}).on(GammaN);
    Solver::KSP hilbertSolver(hilbert);
    hilbertSolver.setPrefix("hilbert_");
    hilbert.assemble();
    hilbert.solve(hilbertSolver);
    enableFactorReuse(hilbertSolver.getHandle());
    if (trace)
      std::cout << "  PETSc hilbert: it="
                << kspIterations(hilbertSolver.getHandle()) << '\n';
    copyPETScToEigen(g.getSolution().getData(), dJ.getData());

    // Write before redistance/advection mutates phiH into the next trial
    // carrier. Thus every frame corresponds to a passed shape objective step.
    if ((outputEvery > 0 && itn % outputEvery == 0) || itn + 1 == maxIt)
    {
      domainGrid.clear();
      domainGrid.add("phi", phiH, IO::XDMF::Center::Node);
      domainGrid.add("dJ", dJ, IO::XDMF::Center::Node);
      stateGrid.clear();
      stateGrid.setMesh(trimmed, IO::XDMF::MeshPolicy::Transient);
      stateGrid.add("u", us.getSolution(), IO::XDMF::Center::Node);
      xdmf.write(static_cast<Real>(itn)).flush();
    }

    std::string activeMode = reinitMode;
    bool doRedistance = reinitMode != "none" && reinitEvery > 0
      && (itn % reinitEvery == 0 || itn + 1 == maxIt);
    if (!doRedistance && reinitMode == "none" && redistanceMode != "none"
        && redistanceEvery > 0
        && ((itn + 1) % redistanceEvery == 0 || itn + 1 == maxIt))
    {
      activeMode = redistanceMode;
      doRedistance = true;
    }
    if (doRedistance)
    {
      Distance::Eikonal<ScalarP1, Math::Vector<Real>>(dist)
        .setInterior(Interior).setInterface(Gamma).solve().sign();
      if (activeMode == "cp")
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
      else if (activeMode != "fmm")
      {
        std::cerr << "Unsupported 3D redistance mode: " << activeMode << '\n';
        return 2;
      }
    }
    else
    {
      dist.getData() = phiH.getData();
    }

    speed = Frobenius(dJ);
    dJ /= std::max(speed.max(), Real(1e-30));
    Advection::Lagrangian(adv, advTest, dist, dJ).step(dt);
    phiH.getData() = adv.getSolution().getData();
    if (repairEvery > 0
        && ((itn + 1) % repairEvery == 0 || itn + 1 == maxIt))
    {
      const Math::Vector<Real> rhs = mass.getOperator() * phiH.getData();
      phiH.getData() = repairCG.solve(rhs);
    }

    ++acceptedShapeIterations;
  }

  xdmf.close();
  std::cout << "\nDone. Accepted shape iterations: "
            << acceptedShapeIterations << " / " << maxIt
            << "  attempts=" << shapeAttempts
            << ". Objective history in LevelSetWNGIRCantilever3D.obj.txt\n";
  return 0;
}

int main(int argc, char** argv)
{
  std::vector<char*> petscArguments;
  petscArguments.push_back(argv[0]);
  for (int i = 1; i < argc; ++i)
    if (std::string(argv[i]).rfind("--", 0) != 0)
      petscArguments.push_back(argv[i]);
  int petscArgc = static_cast<int>(petscArguments.size());
  char** petscArgv = petscArguments.data();
  PetscInitialize(&petscArgc, &petscArgv, PETSC_NULLPTR, PETSC_NULLPTR);
  const int status = run(argc, argv);
  PetscFinalize();
  return status;
}
