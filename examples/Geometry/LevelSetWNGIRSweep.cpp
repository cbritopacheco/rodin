/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//
// Wavy-circle sweep test for a WNGIR-style level-set mesh motion.
//
// At every frame the background grid is classified from the analytic level
// set. WNGIR then registers the classified skeleton directly to phi = 0 using
// the residual phi(X + u(X)) on the classified interior faces.
//
#include <Rodin/Adaptation.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Math.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>
#include <Rodin/Solver/CG.h>
#include <Rodin/Solver/SparseLU.h>
#include <Rodin/Solid.h>
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
  using Vec2 = Math::SpatialVector<Real>;
  using LocalMesh = Geometry::Mesh<Context::Local>;

  Vec2 vec2(Real x = 0, Real y = 0)
  {
    Vec2 out(2);
    out(0) = x;
    out(1) = y;
    return out;
  }

#ifdef RODIN_WNGIR_P2_DISPLACEMENT
  // Promote every triangle to a quadratic (P2) isoparametric cell by installing
  // a degree-2 parametric transformation seeded from the affine node positions.
  // The WNGIR solve then integrates over genuinely curved geometry.
  void installP2CellTransformations(LocalMesh& mesh)
  {
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
        pm(0, a) = X(0);
        pm(1, a) = X(1);
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
      const Vec2 x = mesh.getVertexCoordinates(vertex);
      const auto& dofs = uFes.getDOFs(0, vertex);
      moved.setVertexCoordinates(
        vertex, vec2(x(0) + uData(dofs[0]), x(1) + uData(dofs[1])));
    }

#ifdef RODIN_WNGIR_P2_DISPLACEMENT
    // Carry the full P2 displacement (corner + mid-edge nodes) onto the moved
    // mesh as a curved parametric transformation, so the displaced geometry is
    // genuinely quadratic rather than its affine corner approximation.
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
      moved.setPolytopeTransformation({D, cell.getIndex()},
        new Geometry::ParametricTransformation<Variational::RealH1Element<2>>(
          std::move(pm), geomFe));
    }
#endif
  }

  struct AdmissibilityDistribution
  {
      Real minJ = std::numeric_limits<Real>::infinity();
      Real p10J = 0, p50J = 0, p90J = 0;
      Real maxQ = 0;
      Real p10Q = 0, p50Q = 0, p90Q = 0;
      std::size_t numQuadPoints = 0;
      std::size_t numNearJSafe = 0;   // j_K < 1.5·j_safe
      std::size_t numNearQMax = 0;    // Q_rel > 0.7·Q_max
      std::size_t numInadmissible = 0;
  };

  /// Sampled per-cell admissibility distribution (j, Q_rel) at all
  /// quadrature points of the displacement field. Returns percentiles
  /// (10/50/90) of both margins plus counts of cells "near" each
  /// constraint boundary. Lighter than sampled admissibility diagnostic
  /// only conceptually — same dominant cost (cell-loop gradU eval),
  /// but accumulates a distribution rather than min/max.
  template <class Displacement>
  AdmissibilityDistribution evaluateAdmissibilityDistribution(
    Displacement& u, Real jSafe, Real qMax, std::size_t qOrder = 2)
  {
    using Variational::IntegrationPoint;
    using Variational::Jacobian;

    AdmissibilityDistribution dist;
    const auto& fes = u.getFiniteElementSpace();
    const auto& mesh = fes.getMesh();
    const std::size_t dim = mesh.getDimension();
    auto gradU = Jacobian(u);

    std::vector<Real> jSamples;
    std::vector<Real> qSamples;
    jSamples.reserve(mesh.getCellCount() * 4);
    qSamples.reserve(mesh.getCellCount() * 4);

    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const auto& cell = *cellIt;
      const auto& fe = fes.getFiniteElement(cell.getDimension(), cell.getIndex());
      const auto& qf = QF::PolytopeQuadratureFormula::get(qOrder, cell.getGeometry());
      const auto& quadrature = cell.getQuadrature(qf);
      for (std::size_t q = 0; q < quadrature.getSize(); ++q)
      {
        const auto& pt = quadrature.getPoint(q);
        const IntegrationPoint ip(pt, &qf, q);
        Math::SpatialMatrix<Real> F =
          Math::SpatialMatrix<Real>::Identity(dim, dim) + gradU.getValue(ip);
        const Real j = F.determinant();
        Real qRel = std::numeric_limits<Real>::infinity();
        if (j > Real(0))
          qRel = F.squaredNorm() /
            (static_cast<Real>(dim) * std::pow(j, Real(2) / static_cast<Real>(dim)));
        jSamples.push_back(j);
        qSamples.push_back(qRel);
        if (j < dist.minJ)
          dist.minJ = j;
        if (qRel > dist.maxQ)
          dist.maxQ = qRel;
        if (j <= jSafe)
          ++dist.numInadmissible;
        if (j < Real(1.5) * jSafe)
          ++dist.numNearJSafe;
        if (qRel > Real(0.7) * qMax)
          ++dist.numNearQMax;
      }
    }
    dist.numQuadPoints = jSamples.size();
    if (!jSamples.empty())
    {
      auto pct = [](std::vector<Real>& v, double p) -> Real {
        const std::size_t k = static_cast<std::size_t>(p * static_cast<double>(v.size()));
        std::nth_element(v.begin(), v.begin() + k, v.end());
        return v[k];
      };
      dist.p10J = pct(jSamples, 0.10);
      dist.p50J = pct(jSamples, 0.50);
      dist.p90J = pct(jSamples, 0.90);
      dist.p10Q = pct(qSamples, 0.10);
      dist.p50Q = pct(qSamples, 0.50);
      dist.p90Q = pct(qSamples, 0.90);
    }
    return dist;
  }

  struct WavyCircleLevelSet
  {
      Real cx = Real(0.5);
      Real cy = Real(0.5);
      Real R0 = Real(0.20);
      Real amp = Real(0.05);
      Real k = Real(6);
      Real phase = Real(0);

      Real phi(const Vec2& p) const
      {
        const Real dx = p(0) - cx;
        const Real dy = p(1) - cy;
        const Real r = std::sqrt(dx * dx + dy * dy);
        const Real theta = std::atan2(dy, dx);
        return r - (R0 + amp * std::cos(k * theta + phase));
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
          dx / r - dRdtheta * (-dy / r2safe), dy / r - dRdtheta * (dx / r2safe));
      }
  };

  constexpr std::array<std::array<Real, 3>, 3> TriangleBarycentricQuadrature = {
    {{{Real(2) / 3, Real(1) / 6, Real(1) / 6}}, {{Real(1) / 6, Real(2) / 3, Real(1) / 6}},
      {{Real(1) / 6, Real(1) / 6, Real(2) / 3}}}};

  Real applyPhaseMomentMap(Real phi, Real epsilon)
  {
    return std::tanh(phi / epsilon);
  }

  Vec2 interpolateVec(const std::array<Vec2, 3>& values, const std::array<Real, 3>& bary)
  {
    return bary[0] * values[0] + bary[1] * values[1] + bary[2] * values[2];
  }

  struct CellMomentInfo
  {
      Index index = 0;
      Real area = 0;
      Real moment = 0;
      std::array<Vec2, 3> x;
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
      if (vertices.size() != 3)
        throw std::runtime_error("LevelSetWNGIRSweep expects triangular cells.");

      CellMomentInfo info;
      info.index = cell.getIndex();
      for (std::size_t i = 0; i < 3; ++i)
        info.x[i] = mesh.getVertexCoordinates(vertices[i]);

      const Vec2 e1 = info.x[1] - info.x[0];
      const Vec2 e2 = info.x[2] - info.x[0];
      info.area = std::abs(Real(0.5) * (e1(0) * e2(1) - e1(1) * e2(0)));

      Real moment = 0;
      for (const auto& bary : TriangleBarycentricQuadrature)
        moment += applyPhaseMomentMap(phi(interpolateVec(info.x, bary)), epsilon);
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
  const std::size_t n = parseSizeTOption(argc, argv, "n", 50);
  const std::size_t nFrames = parseSizeTOption(argc, argv, "frames", 40);

  const Real orbitR = parseRealOption(argc, argv, "orbitR", Real(0.10));
  const Real amp = parseRealOption(argc, argv, "amp", Real(0.05));
  const Real R0 = parseRealOption(argc, argv, "R0", Real(0.20));
  const Real kLobes = parseRealOption(argc, argv, "lobes", Real(6));

  const Real h = Real(1) / static_cast<Real>(n - 1);
  const Real epsilon = parseRealOption(argc, argv, "classifier-eps", Real(1.25) * h);
  const Real lambdaC = parseRealOption(argc, argv, "classifier-lambda", Real(0.008));
  const bool verbose = hasFlag(argc, argv, "verbose");

  constexpr Attribute interiorAttribute = 1;
  constexpr Attribute exteriorAttribute = 2;
  constexpr Attribute interfaceAttribute = 10;
  constexpr Attribute boundaryAttribute = 20;

  Rodin::Examples::WNGIRExampleDefaults wngirDefaults;
  wngirDefaults.maxIterations = 120;
#ifdef RODIN_WNGIR_P2_DISPLACEMENT
  wngirDefaults.betaMax = 10;
#endif
  const auto wngirParams = Rodin::Examples::makeWNGIRParameters(
    argc, argv, h, interfaceAttribute, wngirDefaults);
  const Real fitTol = parseRealOption(argc, argv, "fit-tol", wngirParams.activeRMSTol);
  const std::size_t qOrder = wngirParams.quadratureOrder;
  const bool trace = wngirParams.trace;

  LocalMesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {n, n});
  mesh.scale(h);
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 2);
  mesh.getConnectivity().compute(2, 2);
  mesh.getConnectivity().compute(0, 0);
#ifdef RODIN_WNGIR_P2_DISPLACEMENT
  installP2CellTransformations(mesh);
#endif

  for (auto faceIt = mesh.getBoundary(); faceIt; ++faceIt)
    mesh.setAttribute({mesh.getDimension() - 1, faceIt->getIndex()}, boundaryAttribute);

  using ScalarP1 = P1<Real, LocalMesh>;
  using ScalarP0 = P0<Real, LocalMesh>;
#ifdef RODIN_WNGIR_P2_DISPLACEMENT
  using VectorFES = H1<2, Math::SpatialVector<Real>, LocalMesh>;
#else
  using VectorFES = P1<Math::SpatialVector<Real>, LocalMesh>;
#endif

  ScalarP1 p1Fes(mesh);
  ScalarP0 p0Fes(mesh);
#ifdef RODIN_WNGIR_P2_DISPLACEMENT
  VectorFES vectorFes(std::integral_constant<std::size_t, 2>{}, mesh, 2);
#else
  VectorFES vectorFes(mesh, 2);
#endif

  GridFunction phiGf(p1Fes);
  phiGf.setName("phi");
  GridFunction cellLabel(p0Fes);
  cellLabel.setName("cell_label");
  GridFunction phaseMoment(p0Fes);
  phaseMoment.setName("phase_moment");
  TrialFunction wngirTrial(vectorFes);
  TestFunction wngirTest(vectorFes);
  auto& u = wngirTrial.getSolution();
  u.setName("displacement");
  GridFunction du(vectorFes);
  du.setName("wngir_step");
  auto wngirSolveParams = wngirParams;
  wngirSolveParams.activeRMSTol = fitTol;
  Rodin::Adaptation::WNGIR wngirSolver(wngirTrial, wngirTest);
  wngirSolver.setParameters(wngirSolveParams);

  LocalMesh moved(mesh);
  ScalarP0 p0FesMoved(moved);
  ScalarP1 p1FesMoved(moved);
  GridFunction movedLabel(p0FesMoved);
  movedLabel.setName("cell_label");
  GridFunction phiMoved(p1FesMoved);
  phiMoved.setName("phi_moved");
  GridFunction jMoved(p0FesMoved);
  jMoved.setName("j");
  GridFunction qRelMoved(p0FesMoved);
  qRelMoved.setName("q_rel");

  IO::XDMF xdmf("LevelSetWNGIRSweep");
  auto backgroundGrid = xdmf.grid("background");
  backgroundGrid.setMesh(mesh, IO::XDMF::MeshPolicy::Transient);
  backgroundGrid.add(cellLabel, IO::XDMF::Center::Cell);
  backgroundGrid.add(phaseMoment, IO::XDMF::Center::Cell);
  backgroundGrid.add(phiGf, IO::XDMF::Center::Node);
  backgroundGrid.add(u, IO::XDMF::Center::Node);
  backgroundGrid.add(du, IO::XDMF::Center::Node);

  auto movedGrid = xdmf.grid("moved");
  movedGrid.setMesh(moved, IO::XDMF::MeshPolicy::Transient);
  movedGrid.add(movedLabel, IO::XDMF::Center::Cell);
  movedGrid.add(jMoved, IO::XDMF::Center::Cell);
  movedGrid.add(qRelMoved, IO::XDMF::Center::Cell);
  movedGrid.add(phiMoved, IO::XDMF::Center::Node);

  std::cout << "Wavy-circle WNGIR sweep on " << n << "x" << n << " unit-square mesh, "
            << nFrames << " frames\n";
  std::cout << "  R0=" << R0 << "  amp=" << amp << "  k=" << kLobes
            << "  orbit R=" << orbitR << "  wngirEll=" << wngirParams.ellM << '\n';

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
    levelSet.phase = angle;

    std::cout << "\n--- Frame " << std::setw(2) << frame << " : c=(" << std::fixed
              << std::setprecision(4) << levelSet.cx << ", " << levelSet.cy << ")"
              << "  phase=" << std::setprecision(3) << angle << " rad\n";

    clearXDMFRegionAttributes(mesh);
    for (auto faceIt = mesh.getBoundary(); faceIt; ++faceIt)
      mesh.setAttribute({mesh.getDimension() - 1, faceIt->getIndex()}, boundaryAttribute);

    const auto cellMoments = collectCellMomentInfo(
      mesh, [&](const Vec2& p) { return levelSet.phi(p); }, epsilon);

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
      graphEdges.push_back({static_cast<Index>(itA->second),
        static_cast<Index>(itB->second), lambdaC * facetLength(mesh, facet), facet});
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
      return levelSet.phi(vec2(X(0), X(1)));
    };

    RealFunction phi([&](const Geometry::Point& p) -> Real {
      const auto& X = p.getPhysicalCoordinates();
      return levelSet.phi(vec2(X(0), X(1)));
    });
    AnalyticVectorFunction gradPhi(
      [&](const Geometry::Point& p) -> Math::SpatialVector<Real> {
        const auto& X = p.getPhysicalCoordinates();
        return levelSet.grad(vec2(X(0), X(1)));
      },
      /*dimension=*/2);

    u.getData().setZero();
    du.getData().setZero();

    auto computeInterfaceFit = [&]() -> Real {
      Real interfacePhi = 0;
      Real interfaceLen = 0;
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
          Vec2 ux = vec2(0, 0);
          for (std::size_t l = 0; l < nLocal; ++l)
          {
            const auto bv = fe.getBasis(l)(rc);
            ux(0) += bv(0) * u.getData()(dofs[l]);
            ux(1) += bv(1) * u.getData()(dofs[l]);
          }
          const Vec2 y = vec2(Xp(0) + ux(0), Xp(1) + ux(1));
          const Real phiVal = levelSet.phi(y);
          const Real w = qf.getWeight(q) * pt.getDistortion();
          interfacePhi += w * phiVal * phiVal;
          interfaceLen += w;
        }
      }
      return std::sqrt(interfacePhi / std::max(interfaceLen, Real(1e-30)));
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
    Real bestFit = interfaceFit;
    Math::Vector<Real> bestU = u.getData();
    Real minJ = Real(1);
    Real maxQRel = Real(1);
    Real lastAlpha = Real(0);
    Real maxStep = Real(0);
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
    const auto bestAdm =
      evaluateWNGIRAdmissibilitySampled(u, u.getData(), wngirParams.jMinRatio, qOrder);
    minJ = bestAdm.minJ;
    maxQRel = bestAdm.maxQRel;

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

    auto [srcCache, srcToLocal] = precomputeCellGeometry(mesh);
    auto [dstCache, dstToLocal] = precomputeCellGeometry(moved);
    for (auto cellIt = moved.getCell(); cellIt; ++cellIt)
    {
      const Index cellIdx = cellIt->getIndex();
      const std::size_t srcLocal = srcToLocal.at(cellIdx);
      const std::size_t dstLocal = dstToLocal.at(cellIdx);
      const Index dof = p0FesMoved.getGlobalIndex({D, cellIdx}, 0);
      const auto& src = srcCache[srcLocal];
      const auto& dst = dstCache[dstLocal];
      const Real sigDetAu = static_cast<Real>(src.sigmaK) * dst.detAK;
      const Real jK = sigDetAu / src.Jscale;
      jMoved.getData()(dof) = jK;
      const Math::SpatialMatrix<Real> F = dst.A * src.A.inverse();
      qRelMoved.getData()(dof) = F.squaredNorm() / (Real(2) * std::max(jK, Real(1e-30)));
      movedLabel.getData()(dof) =
        static_cast<Real>(classified.labels[cellToLocal.at(cellIdx)]);
    }

    phiMoved = [&](const Geometry::Point& p) -> Real {
      const auto& X = p.getCoordinates();
      return levelSet.phi(vec2(X(0), X(1)));
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
