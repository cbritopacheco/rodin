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
#include <Rodin/Assembly.h>
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
  //   --shape-weight=<w>   weight of E_shape, default 1e-1.
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
  const Real jBarrierWeight =
    getOptionReal(argc, argv, "j-barrier-weight", Real(0));
  const Real jBarrierSafeRatio =
    getOptionReal(argc, argv, "j-barrier-safe", Real(0));
  const Real jVolumeTetherWeight =
    getOptionReal(argc, argv, "volume-tether-weight", Real(0));
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

  WavyCircleLevelSet levelSet;
  levelSet.cx = cx;
  levelSet.cy = cy;
  levelSet.R0 = R0;
  levelSet.amp = amp;
  levelSet.lobes = lobes;
  levelSet.phase = phase;
  const auto cellMoments = collectCellMomentInfo(mesh, levelSet, epsilon);

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
  for (auto faceIt = mesh.getBoundary(); faceIt; ++faceIt)
    mesh.setAttribute({mesh.getDimension() - 1, faceIt->getIndex()},
                      boundaryAttribute);

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

  // -------------------------------------------------------------------------
  // phi is analytic in this prototype; no discrete phi-field is built.
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

  bool solveCompleted = true;
  bool newtonConverged = false;
  std::size_t newtonIterations = 0;
  Real initialGuessAlpha = 0;
  Real initialGuessMinJ = 0;
  std::size_t initialGuessBacktracks = 0;
  std::string solveError;
  constexpr Real geometryFitTolMult = Real(4);
  const Real geometryTolerance = geometryFitTolMult * h * h;
  auto interfacePhiRMSFromDisplacement = [&]() -> Real
  {
    const auto& uFes = u.getFiniteElementSpace();
    const auto& uData = u.getData();
    Real interfacePhi = 0;
    Real interfaceLength = 0;
    for (const Index facet : interfaceFacets)
    {
      const auto face = mesh.getFace(facet);
      const auto& vertices = face->getVertices();
      const Vec2 xa = mesh.getVertexCoordinates(vertices[0]);
      const Vec2 xb = mesh.getVertexCoordinates(vertices[1]);
      const auto& dofsA = uFes.getDOFs(0, vertices[0]);
      const auto& dofsB = uFes.getDOFs(0, vertices[1]);
      const Vec2 a = vec2(xa(0) + uData(dofsA[0]), xa(1) + uData(dofsA[1]));
      const Vec2 b = vec2(xb(0) + uData(dofsB[0]), xb(1) + uData(dofsB[1]));
      const Vec2 mid = Real(0.5) * (a + b);
      const Real val = levelSet.phi(mid);
      const Real len = (b - a).norm();
      interfacePhi += val * val * len;
      interfaceLength += len;
    }
    return std::sqrt(interfacePhi / std::max(interfaceLength, Real(1e-30)));
  };

  {
    LSRParameters lsrParams;
    lsrParams.rhoS = params.rhoS;
    lsrParams.deltaW = params.deltaW;
    lsrParams.hRef = params.hRef;
    lsrParams.normalizer = params.normalizer;
    lsrParams.shapeWeight = shapeWeight;
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
    lsrParams.initialGuess = LSRInitialGuess::Hilbert;
    lsrParams.tangent = LSRTangent::PSDProjectedNewton;
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

    try
    {
      std::cout << "Phase LSR solve: begin\n";
      LSR lsr(u);
      const LSRReport lsrReport =
        lsr.setParameters(lsrParams).solve(
            psi, phi, gradPhi, hessPhi);
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
  // Step 7: build moved mesh through FES dof mapping (no raw indexing).
  // -------------------------------------------------------------------------
  LocalMesh moved(mesh);
  updateMovedMeshFromDisplacement(mesh, moved, u);
  moved.save("LevelSetLSRReconstruction_phi.mesh", IO::FileFormat::MFEM);

  // Interface ||phi(X+u)|| RMS on Gamma_psi,h.
  Real interfacePhi = 0;
  Real interfaceLength = 0;
  for (const Index facet : interfaceFacets)
  {
    const auto face = mesh.getFace(facet);
    const auto& vertices = face->getVertices();
    const Vec2 a = moved.getVertexCoordinates(vertices[0]);
    const Vec2 b = moved.getVertexCoordinates(vertices[1]);
    const Vec2 mid = Real(0.5) * (a + b);
    const Real val = levelSet.phi(mid);
    const Real len = (b - a).norm();
    interfacePhi += val * val * len;
    interfaceLength += len;
  }
  const Real interfacePhiRMS =
    std::sqrt(interfacePhi / std::max(interfaceLength, Real(1e-30)));

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
  std::cout << "\nDiagnostics\n"
            << "  data energy: E_LSR\n"
            << "  shape energy: E_shape, weight=" << shapeWeight << '\n'
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
            << "  LSR tangent mode: PSDProjectedNewton\n"
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
            << "  psi reconstruction: screened Poisson sign source,"
              << " ell=" << psiEll
              << ", magnitude=" << psiSourceMagnitude
              << ", normal-gradient scale=" << psiScale << '\n'
            << "  psi scale diagnostics: numerator=" << psiScaleNumerator
              << ", denominator=" << psiScaleDenominator
              << ", sign moment=" << psiScaleSignMoment << '\n'
            << "  phi source: analytic (WavyCircleLevelSet)\n"
            << "  |Omega_h|: " << domainMeasure
              << ", M_w (smooth band): " << weightedBandMeasure << '\n'
            << "  final ||phi(X+u)||_RMS on Gamma_psi,h: "
              << interfacePhiRMS << '\n'
            << "  geometry tolerance: " << geometryTolerance << '\n'
            << "  Hilbert initial guess: alpha=" << initialGuessAlpha
              << ", backtracks=" << initialGuessBacktracks
              << ", min_j=" << initialGuessMinJ << '\n'
            << "  Newton iterations: " << newtonIterations
              << " / max " << maxNewtonSteps
              << ", converged: " << (newtonConverged ? "yes" : "no") << '\n'
            << "  total runtime: " << elapsedSeconds(totalStart) << " s\n";

  return newtonConverged ? 0 : 1;
}
