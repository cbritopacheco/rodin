/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//
// Lobed-sphere sweep test for a WNGIR-style level-set mesh motion in 3D.
//
// At every frame the background tetrahedral grid is classified from the
// analytic level set. WNGIR then registers the classified cut-surface
// skeleton directly to phi = 0 using the residual phi(X + u(X)) on the
// classified interior faces. Between frames the target orbits the cube
// center and its lobes rotate in phase, so each frame is a cold-start
// reconstruction of a moving, morphing implicit surface.
//
#include <Rodin/Adaptation.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Math.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>
#include <Rodin/Solver/CG.h>
#include <Rodin/Solver/SparseLU.h>
#include <Rodin/Variational.h>

#include "../WNGIRExampleParameters.h"

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
  using Vec3 = Math::SpatialVector<Real>;
  using Mat3 = Math::SpatialMatrix<Real>;
  using LocalMesh = Geometry::Mesh<Context::Local>;

  Vec3 vec3(Real x = 0, Real y = 0, Real z = 0)
  {
    Vec3 out(3);
    out(0) = x;
    out(1) = y;
    out(2) = z;
    return out;
  }

  Real det3(const Vec3& a, const Vec3& b, const Vec3& c)
  {
    return a(0) * (b(1) * c(2) - b(2) * c(1)) - a(1) * (b(0) * c(2) - b(2) * c(0)) +
      a(2) * (b(0) * c(1) - b(1) * c(0));
  }

  Vec3 cross3(const Vec3& a, const Vec3& b)
  {
    return vec3(
      a(1) * b(2) - a(2) * b(1), a(2) * b(0) - a(0) * b(2), a(0) * b(1) - a(1) * b(0));
  }

  Mat3 edgeMatrix(const std::array<Vec3, 4>& x)
  {
    Mat3 A(3, 3);
    const Vec3 e1 = x[1] - x[0];
    const Vec3 e2 = x[2] - x[0];
    const Vec3 e3 = x[3] - x[0];
    for (int r = 0; r < 3; ++r)
    {
      A(r, 0) = e1(r);
      A(r, 1) = e2(r);
      A(r, 2) = e3(r);
    }
    return A;
  }

  struct LobedSphereLevelSet
  {
      Vec3 c = vec3(Real(0.5), Real(0.5), Real(0.5));
      Real R0 = Real(0.25);
      Real amp = Real(0.05);
      Real lobes = Real(6);
      Real phase = Real(0);

      Real radius(const Vec3& p) const
      {
        const Vec3 d = p - c;
        const Real r = std::max(d.norm(), Real(1e-14));
        const Real theta = std::atan2(d(1), d(0));
        const Real mu = d(2) / r;
        return R0 + amp * std::cos(lobes * theta + phase) * (Real(1) - mu * mu);
      }

      Real phi(const Vec3& p) const
      {
        return (p - c).norm() - radius(p);
      }

      Vec3 grad(const Vec3& p) const
      {
        const Vec3 x = p - c;
        const Real x0 = x(0);
        const Real x1 = x(1);
        const Real x2 = x(2);
        const Real r2 = x.squaredNorm();
        const Real r = std::sqrt(r2);
        if (r <= Real(1e-14))
          return vec3(0, 0, 0);

        const Real rho2 = x0 * x0 + x1 * x1;
        const Real theta = std::atan2(x1, x0);
        const Real angle = lobes * theta + phase;
        const Real cosA = std::cos(angle);
        const Real sinA = std::sin(angle);

        const Real invR2 = Real(1) / r2;
        const Real invR4 = invR2 * invR2;
        const Vec3 gradR = amp *
          (-sinA * lobes * vec3(-x1 * invR2, x0 * invR2, 0) +
            cosA *
              vec3(Real(2) * x0 * x2 * x2 * invR4, Real(2) * x1 * x2 * x2 * invR4,
                -Real(2) * x2 * rho2 * invR4));

        return x / r - gradR;
      }
  };

  constexpr std::array<std::array<Real, 4>, 4> TetraBarycentricQuadrature = {
    {{{Real(0.5854101966249685), Real(0.1381966011250105), Real(0.1381966011250105),
       Real(0.1381966011250105)}},
      {{Real(0.1381966011250105), Real(0.5854101966249685), Real(0.1381966011250105),
        Real(0.1381966011250105)}},
      {{Real(0.1381966011250105), Real(0.1381966011250105), Real(0.5854101966249685),
        Real(0.1381966011250105)}},
      {{Real(0.1381966011250105), Real(0.1381966011250105), Real(0.1381966011250105),
        Real(0.5854101966249685)}}}};

  Real applyPhaseMomentMap(Real phi, Real epsilon)
  {
    return std::tanh(phi / epsilon);
  }

  Vec3 interpolateVec(const std::array<Vec3, 4>& values, const std::array<Real, 4>& bary)
  {
    Vec3 out = vec3(0, 0, 0);
    for (std::size_t i = 0; i < 4; ++i)
      out += bary[i] * values[i];
    return out;
  }

  struct CellMomentInfo
  {
      Index index = 0;
      Real volume = 0;
      Real moment = 0;
      std::array<Vec3, 4> x;
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
      if (vertices.size() != 4)
        throw std::runtime_error("LevelSetWNGIRSweep3D expects tetrahedral cells.");

      CellMomentInfo info;
      info.index = cell.getIndex();
      for (std::size_t i = 0; i < 4; ++i)
        info.x[i] = mesh.getVertexCoordinates(vertices[i]);

      info.volume = std::abs(det3(info.x[1] - info.x[0], info.x[2] - info.x[0],
                      info.x[3] - info.x[0])) /
        Real(6);

      Real moment = 0;
      for (const auto& bary : TetraBarycentricQuadrature)
        moment += applyPhaseMomentMap(phi(interpolateVec(info.x, bary)), epsilon);
      info.moment = moment / TetraBarycentricQuadrature.size();
      cells.push_back(std::move(info));
    }
    return cells;
  }

  Real faceArea(const LocalMesh& mesh, Index facet)
  {
    const auto face = mesh.getFace(facet);
    const auto& vertices = face->getVertices();
    const Vec3 a = mesh.getVertexCoordinates(vertices[0]);
    const Vec3 b = mesh.getVertexCoordinates(vertices[1]);
    const Vec3 c = mesh.getVertexCoordinates(vertices[2]);
    return Real(0.5) * cross3(b - a, c - a).norm();
  }

  void clearXDMFRegionAttributes(LocalMesh& mesh)
  {
    const std::size_t D = mesh.getDimension();
    for (Index i = 0; i < static_cast<Index>(mesh.getCellCount()); ++i)
      mesh.setAttribute({D, i}, Attribute{0});
    for (auto it = mesh.getFace(); it; ++it)
      mesh.setAttribute({D - 1, it->getIndex()}, Attribute{0});
  }

#ifdef RODIN_WNGIR_P2_DISPLACEMENT
  // Promote every tetrahedron to a quadratic (P2) isoparametric cell by
  // installing a degree-2 parametric transformation seeded from the affine
  // node positions. The WNGIR solve then integrates over curved geometry.
  void installP2CellTransformations(LocalMesh& mesh)
  {
    Variational::RealH1Element<2> geomFe(Polytope::Type::Tetrahedron);
    Math::SpatialPoint X;
    const std::size_t D = mesh.getDimension();
    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const auto& cell = *cellIt;
      Geometry::PointCloud pm(3, geomFe.getCount());
      for (std::size_t a = 0; a < geomFe.getCount(); ++a)
      {
        const auto& rc = geomFe.getNode(a);
        cell.getTransformation().transform(X, rc);
        pm(0, a) = X(0);
        pm(1, a) = X(1);
        pm(2, a) = X(2);
      }
      mesh.setPolytopeTransformation({D, cell.getIndex()},
        new Geometry::ParametricTransformation<Variational::RealH1Element<2>>(
          std::move(pm), geomFe));
    }
  }
#endif

  template <class Displacement>
  void updateMovedMeshFromDisplacement(
    const LocalMesh& mesh, LocalMesh& moved, const Displacement& u)
  {
    const auto& uFes = u.getFiniteElementSpace();
    const auto& uData = u.getData();
    const Index vn = mesh.getVertexCount();
    for (Index vertex = 0; vertex < vn; ++vertex)
    {
      const Vec3 x = mesh.getVertexCoordinates(vertex);
      const auto& dofs = uFes.getDOFs(0, vertex);
      moved.setVertexCoordinates(vertex,
        vec3(x(0) + uData(dofs[0]), x(1) + uData(dofs[1]), x(2) + uData(dofs[2])));
    }

#ifdef RODIN_WNGIR_P2_DISPLACEMENT
    // Carry the full P2 displacement (corner + mid-edge nodes) onto the moved
    // mesh as a curved parametric transformation.
    Variational::RealH1Element<2> geomFe(Polytope::Type::Tetrahedron);
    Math::SpatialPoint X;
    const std::size_t D = mesh.getDimension();
    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const auto& cell = *cellIt;
      Geometry::PointCloud pm(3, geomFe.getCount());
      for (std::size_t a = 0; a < geomFe.getCount(); ++a)
      {
        const auto& rc = geomFe.getNode(a);
        cell.getTransformation().transform(X, rc);
        const Real ux = uData(uFes.getGlobalIndex({D, cell.getIndex()}, a * 3));
        const Real uy = uData(uFes.getGlobalIndex({D, cell.getIndex()}, a * 3 + 1));
        const Real uz = uData(uFes.getGlobalIndex({D, cell.getIndex()}, a * 3 + 2));
        pm(0, a) = X(0) + ux;
        pm(1, a) = X(1) + uy;
        pm(2, a) = X(2) + uz;
      }
      moved.setPolytopeTransformation({D, cell.getIndex()},
        new Geometry::ParametricTransformation<Variational::RealH1Element<2>>(
          std::move(pm), geomFe));
    }
#endif
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

  Real parseRealOption(int argc, char** argv, const std::string& name, Real fallback)
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
  const std::size_t n = parseSizeTOption(argc, argv, "n", 24);
  const std::size_t nFrames = parseSizeTOption(argc, argv, "frames", 24);

  const Real orbitR = parseRealOption(argc, argv, "orbitR", Real(0.08));
  const Real amp = parseRealOption(argc, argv, "amp", Real(0.045));
  const Real R0 = parseRealOption(argc, argv, "R0", Real(0.24));
  const Real kLobes = parseRealOption(argc, argv, "lobes", Real(6));

  const Real h = Real(1) / static_cast<Real>(n - 1);
  const Real epsilon = parseRealOption(argc, argv, "classifier-eps", Real(1.25) * h);
  const Real lambdaC = parseRealOption(argc, argv, "classifier-lambda", Real(0.004));
  Rodin::Examples::WNGIRExampleDefaults wngirDefaults;
  wngirDefaults.maxIterations = 120;
  wngirDefaults.gammaMFactor = 0;
  wngirDefaults.gammaHFactor = Real(0.0125);
  wngirDefaults.gammaDivFactor = Real(0.0125);
  wngirDefaults.ellOverH = Real(0.75);
  wngirDefaults.activeRMSOverHTol = Real(0.03);
  wngirDefaults.activeSupOverHTol = Real(0.20);
#ifdef RODIN_WNGIR_P2_DISPLACEMENT
  wngirDefaults.betaMax = 50;
#else
  wngirDefaults.betaMax = 100;
#endif
  const bool verbose = hasFlag(argc, argv, "verbose");

  constexpr Attribute interiorAttribute = 1;
  constexpr Attribute exteriorAttribute = 2;
  constexpr Attribute interfaceAttribute = 10;
  constexpr Attribute boundaryAttribute = 20;

  const auto wngirParams = Rodin::Examples::makeWNGIRParameters(
    argc, argv, h, interfaceAttribute, wngirDefaults);
  const Real fitTol = parseRealOption(argc, argv, "fit-tol", wngirParams.activeRMSTol);
  const std::size_t qOrder = wngirParams.quadratureOrder;
  const bool trace = wngirParams.trace;

  LocalMesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {n, n, n});
  mesh.scale(h);
  mesh.getConnectivity().compute(3, 2);
  mesh.getConnectivity().compute(2, 3);
  mesh.getConnectivity().compute(0, 0);
#ifdef RODIN_WNGIR_P2_DISPLACEMENT
  // P2 needs edge entities (dim 1) for the mid-edge degrees of freedom.
  mesh.getConnectivity().compute(1, 0);
  mesh.getConnectivity().compute(3, 1);
  mesh.getConnectivity().compute(2, 1);
  installP2CellTransformations(mesh);
#endif

  for (auto faceIt = mesh.getBoundary(); faceIt; ++faceIt)
    mesh.setAttribute({mesh.getDimension() - 1, faceIt->getIndex()}, boundaryAttribute);

#ifdef RODIN_WNGIR_P2_DISPLACEMENT
  using ScalarFES = H1<2, Real, LocalMesh>;
  using VectorFES = H1<2, Math::SpatialVector<Real>, LocalMesh>;
#else
  using ScalarFES = P1<Real, LocalMesh>;
  using VectorFES = P1<Math::SpatialVector<Real>, LocalMesh>;
#endif
  using ScalarP0 = P0<Real, LocalMesh>;

  ScalarFES scalarFes(
#ifdef RODIN_WNGIR_P2_DISPLACEMENT
    std::integral_constant<std::size_t, 2>{},
#endif
    mesh);
  ScalarP0 p0Fes(mesh);
#ifdef RODIN_WNGIR_P2_DISPLACEMENT
  VectorFES vectorFes(std::integral_constant<std::size_t, 2>{}, mesh, 3);
#else
  VectorFES vectorFes(mesh, 3);
#endif

  GridFunction phiGf(scalarFes);
  phiGf.setName("phi");
  GridFunction cellLabel(p0Fes);
  cellLabel.setName("cell_label");
  GridFunction phaseMoment(p0Fes);
  phaseMoment.setName("phase_moment");
  TrialFunction wngirTrial(vectorFes);
  TestFunction wngirTest(vectorFes);
  auto& u = wngirTrial.getSolution();
  u.setName("displacement");
  auto wngirSolveParams = wngirParams;
  wngirSolveParams.activeRMSTol = fitTol;
  Rodin::Adaptation::WNGIR wngirSolver(wngirTrial, wngirTest);
  wngirSolver.setParameters(wngirSolveParams);

  LocalMesh moved(mesh);
  ScalarP0 p0FesMoved(moved);
  ScalarFES scalarFesMoved(
#ifdef RODIN_WNGIR_P2_DISPLACEMENT
    std::integral_constant<std::size_t, 2>{},
#endif
    moved);
  GridFunction movedLabel(p0FesMoved);
  movedLabel.setName("cell_label");
  GridFunction phiMoved(scalarFesMoved);
  phiMoved.setName("phi_moved");
  GridFunction jMoved(p0FesMoved);
  jMoved.setName("j");
  GridFunction qRelMoved(p0FesMoved);
  qRelMoved.setName("q_rel");

  IO::XDMF xdmf("LevelSetWNGIRSweep3D");
  auto backgroundGrid = xdmf.grid("background");
  backgroundGrid.setMesh(mesh, IO::XDMF::MeshPolicy::Transient);
  backgroundGrid.add(cellLabel, IO::XDMF::Center::Cell);
  backgroundGrid.add(phaseMoment, IO::XDMF::Center::Cell);
  backgroundGrid.add(phiGf, IO::XDMF::Center::Node);
  backgroundGrid.add(u, IO::XDMF::Center::Node);

  auto movedGrid = xdmf.grid("moved");
  movedGrid.setMesh(moved, IO::XDMF::MeshPolicy::Transient);
  movedGrid.add(movedLabel, IO::XDMF::Center::Cell);
  movedGrid.add(jMoved, IO::XDMF::Center::Cell);
  movedGrid.add(qRelMoved, IO::XDMF::Center::Cell);
  movedGrid.add(phiMoved, IO::XDMF::Center::Node);

  std::cout << "Lobed-sphere WNGIR sweep on " << n << "x" << n << "x" << n
            << " tetrahedral unit-cube mesh, " << nFrames << " frames\n";
  std::cout << "  R0=" << R0 << "  amp=" << amp << "  lobes=" << kLobes
            << "  orbit R=" << orbitR << "  wngirEll=" << wngirParams.ellM
            << "  betaMax=" << wngirParams.betaMax << '\n';

  std::size_t framesConverged = 0;
  std::vector<Real> finalFitPerFrame;
  finalFitPerFrame.reserve(nFrames);

  for (std::size_t frame = 0; frame < nFrames; ++frame)
  {
    const Real t = static_cast<Real>(frame) / static_cast<Real>(nFrames);
    const Real angle = Real(2) * Real(M_PI) * t;

    LobedSphereLevelSet levelSet;
    levelSet.c =
      vec3(Real(0.5) + orbitR * std::cos(angle), Real(0.5) + orbitR * std::sin(angle),
        Real(0.5) + Real(0.5) * orbitR * std::sin(Real(2) * angle));
    levelSet.R0 = R0;
    levelSet.amp = amp;
    levelSet.lobes = kLobes;
    levelSet.phase = angle;

    std::cout << "\n--- Frame " << std::setw(2) << frame << " : c=(" << std::fixed
              << std::setprecision(4) << levelSet.c(0) << ", " << levelSet.c(1) << ", "
              << levelSet.c(2) << ")"
              << "  phase=" << std::setprecision(3) << angle << " rad\n";

    clearXDMFRegionAttributes(mesh);
    for (auto faceIt = mesh.getBoundary(); faceIt; ++faceIt)
      mesh.setAttribute({mesh.getDimension() - 1, faceIt->getIndex()}, boundaryAttribute);

    const auto cellMoments = collectCellMomentInfo(
      mesh, [&](const Vec3& p) { return levelSet.phi(p); }, epsilon);

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
      volumes[local] = cellMoments[local].volume;
      moments[local] = cellMoments[local].moment;
    }

    std::vector<MinSTCut::Edge> graphEdges;
    for (auto faceIt = mesh.getFace(); faceIt; ++faceIt)
    {
      const Index facet = faceIt->getIndex();
      const auto& incident = mesh.getConnectivity().getIncidence({2, 3}, facet);
      if (incident.size() != 2)
        continue;
      const auto itA = cellToLocal.find(incident[0]);
      const auto itB = cellToLocal.find(incident[1]);
      if (itA == cellToLocal.end() || itB == cellToLocal.end())
        continue;
      graphEdges.push_back({static_cast<Index>(itA->second),
        static_cast<Index>(itB->second), lambdaC * faceArea(mesh, facet), facet});
    }

    const MinSTCut cut;
    const MinSTCut::Result classified = cut.classify(volumes, moments, graphEdges);

    std::vector<Index> interfaceFacets;
    interfaceFacets.reserve(classified.cutEdges.size());
    for (const MinSTCut::Edge& edge : classified.cutEdges)
      if (edge.index != MinSTCut::InvalidIndex)
        interfaceFacets.push_back(edge.index);

    for (std::size_t local = 0; local < classified.labels.size(); ++local)
    {
      const Index cellIdx = localToCell[local];
      mesh.setAttribute({mesh.getDimension(), cellIdx},
        classified.labels[local] == MinSTCut::Inside ? interiorAttribute
                                                     : exteriorAttribute);
    }
    for (const Index facet : interfaceFacets)
      mesh.setAttribute({mesh.getDimension() - 1, facet}, interfaceAttribute);

    phiGf = [&](const Geometry::Point& p) -> Real {
      const auto& X = p.getCoordinates();
      return levelSet.phi(vec3(X(0), X(1), X(2)));
    };

    RealFunction phi([&](const Geometry::Point& p) -> Real {
      const auto& X = p.getPhysicalCoordinates();
      return levelSet.phi(vec3(X(0), X(1), X(2)));
    });
    AnalyticVectorFunction gradPhi(
      [&](const Geometry::Point& p) -> Math::SpatialVector<Real> {
        const auto& X = p.getPhysicalCoordinates();
        return levelSet.grad(vec3(X(0), X(1), X(2)));
      },
      /*dimension=*/3);

    u.getData().setZero();

    auto computeInterfaceFit = [&]() -> Real {
      Real interfacePhi = 0;
      Real interfaceArea = 0;
      const auto& fes = u.getFiniteElementSpace();
      const std::size_t meshDim = mesh.getDimension();
      for (const Index facet : interfaceFacets)
      {
        const auto face = mesh.getFace(facet);
        const auto& fe = fes.getFiniteElement(meshDim - 1, facet);
        const std::size_t nLocal = fe.getCount();
        const std::size_t qFitOrder = std::max<std::size_t>(qOrder, 2 * fe.getOrder());
        const auto& qf =
          QF::PolytopeQuadratureFormula::get(qFitOrder, face->getGeometry());
        const auto& quad = face->getQuadrature(qf);
        std::vector<Index> dofs(nLocal);
        for (std::size_t l = 0; l < nLocal; ++l)
          dofs[l] = fes.getGlobalIndex({meshDim - 1, facet}, l);
        for (std::size_t q = 0; q < quad.getSize(); ++q)
        {
          const auto& pt = quad.getPoint(q);
          const auto& Xp = pt.getPhysicalCoordinates();
          const auto& rc = pt.getReferenceCoordinates();
          Vec3 ux = vec3(0, 0, 0);
          for (std::size_t l = 0; l < nLocal; ++l)
          {
            const auto bv = fe.getBasis(l)(rc);
            for (int a = 0; a < 3; ++a)
              ux(a) += bv(a) * u.getData()(dofs[l]);
          }
          const Vec3 y = vec3(Xp(0) + ux(0), Xp(1) + ux(1), Xp(2) + ux(2));
          const Real phiVal = levelSet.phi(y);
          const Real w = qf.getWeight(q) * pt.getDistortion();
          interfacePhi += w * phiVal * phiVal;
          interfaceArea += w;
        }
      }
      return std::sqrt(interfacePhi / std::max(interfaceArea, Real(1e-30)));
    };

    Real interfaceFit = computeInterfaceFit();
    if (verbose)
    {
      std::size_t insideCount = 0;
      for (int lbl : classified.labels)
        if (lbl == MinSTCut::Inside)
          ++insideCount;
      std::cout << "    debug: facets=" << interfaceFacets.size()
                << "  inside=" << insideCount
                << "  outside=" << (classified.labels.size() - insideCount)
                << "  fit0=" << interfaceFit << "\n";
    }

    Math::Vector<Real> bestU = u.getData();
    Real bestFit = interfaceFit;
    Real minJ = Real(1);
    Real maxQRel = Real(1);
    Real lastAlpha = Real(0);
    Real acceptedStep = Real(0);
    std::size_t iterations = 0;
    const char* exitReason = "iter-budget";
    {
      const auto wngirRep = wngirSolver.solve(mesh, interfaceFacets, phi, gradPhi);
      std::cout << "    wngir timing: it=" << wngirRep.iterations << std::scientific
                << std::setprecision(2) << "  assembly=" << wngirRep.tAssembly
                << "  setup=" << wngirRep.tFactor << "  solve=" << wngirRep.tSolve
                << "  cgIt=" << wngirRep.linearIterations
                << "  cgErr=" << wngirRep.linearError << "  ls=" << wngirRep.tLineSearch
                << "  exit=" << wngirRep.exitReason << '\n';
      iterations = wngirRep.iterations;
      lastAlpha = wngirRep.lastAlpha;
      acceptedStep = wngirRep.acceptedStep;
      minJ = wngirRep.minJ;
      maxQRel = wngirRep.maxQRel;
      exitReason = wngirRep.exitReason;
      interfaceFit = computeInterfaceFit();
      if (interfaceFit < bestFit)
      {
        bestFit = interfaceFit;
        bestU = u.getData();
      }
      if (trace)
        std::cout << "      wngir sigma=" << wngirRep.sigma << "  (3h=" << Real(3) * h
                  << ")\n";
    }

    u.getData() = bestU;
    interfaceFit = bestFit;

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
      cellLabel.getData()(dof) = static_cast<Real>(classified.labels[local]);
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

    std::unordered_map<Index, std::size_t> movedCellToLocal;
    std::vector<std::array<Vec3, 4>> movedCells;
    movedCells.reserve(moved.getCellCount());
    movedCellToLocal.reserve(moved.getCellCount());
    for (auto cellIt = moved.getCell(); cellIt; ++cellIt)
    {
      const auto& vertices = cellIt->getVertices();
      std::array<Vec3, 4> x;
      for (std::size_t a = 0; a < 4; ++a)
        x[a] = moved.getVertexCoordinates(vertices[a]);
      movedCellToLocal[cellIt->getIndex()] = movedCells.size();
      movedCells.push_back(x);
    }

    for (auto cellIt = moved.getCell(); cellIt; ++cellIt)
    {
      const Index cellIdx = cellIt->getIndex();
      const std::size_t local = cellToLocal.at(cellIdx);
      const std::size_t movedLocal = movedCellToLocal.at(cellIdx);
      const Index dof = p0FesMoved.getGlobalIndex({D, cellIdx}, 0);
      const Mat3 A0 = edgeMatrix(cellMoments[local].x);
      const Mat3 A1 = edgeMatrix(movedCells[movedLocal]);
      const Real j = A1.determinant() / A0.determinant();
      const Mat3 F = A1 * A0.inverse();
      jMoved.getData()(dof) = j;
      qRelMoved.getData()(dof) = j > Real(0)
        ? F.squaredNorm() / (Real(3) * std::pow(j, Real(2) / Real(3)))
        : std::numeric_limits<Real>::infinity();
      movedLabel.getData()(dof) =
        static_cast<Real>(classified.labels[cellToLocal.at(cellIdx)]);
    }

    phiMoved = [&](const Geometry::Point& p) -> Real {
      const auto& X = p.getCoordinates();
      return levelSet.phi(vec3(X(0), X(1), X(2)));
    };

    std::cout << "    WNGIR it=" << iterations << "  fit=" << std::scientific
              << std::setprecision(3) << interfaceFit << "  alpha=" << lastAlpha
              << "  step=" << acceptedStep << "  min_j=" << minJ
              << "  max_qrel=" << maxQRel
              << "  converged=" << (converged ? "yes" : "best-effort")
              << "  exit=" << exitReason << '\n';

    xdmf.write(t).flush();
  }

  xdmf.close();

  std::cout << "\nSummary\n";
  std::cout << "  frames converged: " << framesConverged << " / " << nFrames << '\n';
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
    std::cout << "  ||phi(X+u)||_RMS  min=" << std::scientific << std::setprecision(3)
              << fitMin << "  mean=" << fitMean << "  max=" << fitMax << '\n';
  }

  return 0;
}
