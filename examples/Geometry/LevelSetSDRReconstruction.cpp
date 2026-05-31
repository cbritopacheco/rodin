/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//
// 2D triangular static prototype of SDFR (signed distance function
// registration) displacement, on a P1 vector finite-element space,
// expressed through Rodin::Variational::Problem.
//
// Scope.
//   - Planar 2D meshes; affine triangle cells.
//   - A single static SDFR solve.
//   - An analytic level set with closed-form phi, grad phi, hess phi.
//
// FES-independence — honest version.
//   - The SDFR integrators (Rodin::Adaptation::SDFR{Residual,Tangent}Integrator)
//     use Rodin's FES interface for per-cell finite-element lookup and
//     vector-valued basis evaluation; they are FES-generic within the
//     planar 2D vector setting and have not been exercised on higher-order
//     or non-affine elements.
//   - The shape + singularity-floor (E_floor) integrators
//     (Rodin::Adaptation::Barrier{Residual,Tangent}Integrator) are NOT
//     FES-independent: they are specialised to affine P1 triangles in 2D
//     (six local dofs, constant grad_X u_h per cell). Extending them
//     requires replacing the 6-dof closed-form algebra by a
//     quadrature-driven evaluation of A_K^u = D(F_K + u_h o F_K).
//   - The phase moment computation in this example is NOT FES-independent
//     either: it uses a hard-coded triangle barycentric quadrature and
//     vertex interpolation.
//   - The output side is FES-correct (P0 cell writes go through
//     `fes.getGlobalIndex({d, cellIdx}, 0)`; vector displacement reads
//     go through `fes.getDOFs(0, vertex)`); the cell-index <-> local-index
//     translation goes through an explicit map (no compact-index
//     assumption).
//
// Energy and well-posedness, in words.
//   - The SDFR data term controls (in the sense of dominating the Hessian)
//     the component of u along grad_y phi at every quadrature point with
//     non-negligible weight. It does NOT control motion tangent to the
//     level set, nor mesh motion in the far field. The Gaussian smoothing
//     widens the support of the SDFR Hessian but does not make the SDFR
//     functional coercive on its own.
//   - The intrinsic shape quality energy (Q_shape - 1, cell-area-weighted)
//     supplies tangential control and far-field regularisation.
//   - The singularity floor E_floor is the only piece protecting
//       j_K^u(xhat) = sigma_K det(A_K^u) / J_K_scale
//     from collapsing through the cell's initial sigma_K branch
//     (j_K^u -> j_min). It is identically zero (and so are its
//     derivatives) whenever j_K^u >= j_safe, so it does NOT bias the
//     element toward any particular volume; it is NOT an orientation
//     barrier; it does NOT prefer j_K^u = 1. It is only a one-sided
//     safeguard against the singular set, sufficient to keep
//     Q_shape well-defined on the same sigma_K branch.
//   - Zero Dirichlet on the mesh boundary pins the perimeter dofs and
//     kills rigid-body modes.
// All four ingredients together make the discrete Newton system
// positive-definite; none of them alone does so.
//
// FullNewton vs GaussNewton — what we actually observe here.
//
//   Mode selection follows the Rodin::Variational::LinearElasticityIntegral
//   pattern. The SDFRRegistration helper captures the
//   VARIATIONAL SKELETON (trial du, test v, state u) at construction,
//   and the DATA (level set, optional Hessian callable, signed-distance
//   field, parameters) is passed at `operator()` / `Tangent` /
//   `Residual` time. Overload resolution on the call signature picks
//   the right SDFRTangentIntegrator specialisation:
//
//     SDFRRegistration sdfrTerm(du, v, u);
//     newton = sdfrTerm(phi, grad,        sLF, params);   // -> GaussNewton
//     newton = sdfrTerm(phi, grad, hess,  sLF, params);   // -> Newton
//
//   The three derivatives of the level set (phi, grad, hess) are
//   first-class equal arguments at the call site, exactly like
//   `LinearElasticityIntegral` passes (lambda, mu) at its operator().
//   No `LevelSet` bundle object hides any of them; the symmetry follows
//   from the Rodin form-language pattern of separating the variational
//   skeleton (held by the helper) from the data (passed per call).
//   There is no runtime `tangentMode` flag and no setter on the helper. The barrier and shape tangents
//   are always the analytical full Hessian (verified against central FD
//   by the unit tests in
//   tests/unit/Rodin/Adaptation/BarrierLocalHessianFDTest.cpp;
//   relative error ~1e-10 at the optimum eps with V-shape signature).
//   Only the SDFR tangent flavour changes between modes.
//
//   GaussNewton (default).
//     - K_GN = (grad phi)(grad phi)^T per quadrature point.
//       Rank-1 PSD; robust; monotone.
//     - The observed local order q[n] -> 1.000 with contraction ratio
//       ~0.43 in the tail. That is plain LINEAR convergence, not
//       super-linear: the dropped second-order term
//
//           r * hess(phi)
//
//       does not vanish at the GN minimum because the SDFR L2 fit is
//       approximate (r* != 0). The linear contraction floor is set by
//       the spectral radius of the GN-vs-Newton tangent mismatch on
//       that nonzero residual.
//
//   Newton (full Hessian).
//     - K_N = (grad phi)(grad phi)^T + r * hess(phi).
//     - phi is the signed distance to a circle. Its Hessian
//
//           hess(phi) = (I - n n^T) / ||p - c||,  n = (p - c)/||p - c||
//
//       is positive semi-definite EVERYWHERE except at the medial axis
//       (the centre point, the only singularity). The interface itself
//       is perfectly smooth.
//     - The Newton tangent is therefore PSD whenever r > 0 and goes
//       INDEFINITE in the tangential direction whenever r < 0: the
//       rank-1 (grad phi)(grad phi)^T fills only the radial direction,
//       leaving the tangential eigenvalue at r / ||p - c|| with the
//       sign of r. This is a textbook outcome for nonlinear
//       least-squares with a non-vanishing residual; nothing about
//       phi or s_h^LF is irregular.
//     - On this prototype Newton matches GaussNewton for the first ~5
//       iterations, then the indefinite Newton step drives a handful
//       of cells past j_K^u = j_min (the singular floor). Those cells
//       are flagged via `barrierInadmissibleCount()` and contribute a
//       small identity
//       block (instead of throwing, which under OpenMP assembly is
//       std::terminate); Newton continues but no longer converges.
//
//   Takeaway: the textbook quadratic Newton regime is not reachable
//   on this objective without one of (a) line search / trust region,
//   (b) damping of the r * hess(phi) term whenever r changes sign,
//   or (c) exact-fit data so r* = 0. None of those are in scope here.
//
#include <Rodin/Adaptation.h>
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

  Vec2 vec2(Real x = 0, Real y = 0)
  {
    Vec2 out(2);
    out(0) = x;
    out(1) = y;
    return out;
  }

  // -------------------------------------------------------------------------
  // Analytic circle level set.
  // -------------------------------------------------------------------------
  struct CircleLevelSet
  {
    Real cx = 0.51;
    Real cy = 0.48;
    Real radius = 0.31;

    Real phi(const Vec2& p) const
    {
      const Real dx = p(0) - cx;
      const Real dy = p(1) - cy;
      return std::sqrt(dx * dx + dy * dy) - radius;
    }

    Vec2 grad(const Vec2& p) const
    {
      const Real dx = p(0) - cx;
      const Real dy = p(1) - cy;
      const Real r = std::max(std::sqrt(dx * dx + dy * dy), Real(1e-14));
      return vec2(dx / r, dy / r);
    }

    Math::SpatialMatrix<Real> hess(const Vec2& p) const
    {
      const Real dx = p(0) - cx;
      const Real dy = p(1) - cy;
      const Real r = std::max(std::sqrt(dx * dx + dy * dy), Real(1e-14));
      Math::SpatialMatrix<Real> H(2, 2);
      H(0, 0) = 1 / r - dx * dx / (r * r * r);
      H(0, 1) = -dx * dy / (r * r * r);
      H(1, 0) = H(0, 1);
      H(1, 1) = 1 / r - dy * dy / (r * r * r);
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

  std::vector<CellMomentInfo> collectCellMomentInfo(
      const Mesh<Context::Local>& mesh,
      const CircleLevelSet& levelSet,
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
            "LevelSetSDRReconstruction expects triangular cells.");

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

}

int main(int argc, char** argv)
{
  constexpr std::size_t n = 50;
  const Real h = Real(1) / static_cast<Real>(n - 1);
  const Real epsilon = 1.25 * h;
  const Real lambdaC = 0.008;

  // Pick the SDFR tangent flavour for this run. Flip to
  // SDFRIntegratorTangentMode::Newton to audit the full-Newton tangent.
  //
  // For this prototype, empirically:
  //   - GaussNewton converges monotonically with factor ~3 per iter early
  //     (super-linear) and continues all the way down to ||R|| ~ 1e-9.
  //   - Newton matches GaussNewton for the first few iterations, then the
  //     indefinite r * hess(phi) term takes over and produces residual
  //     spikes. On this circle level set inside the sphere phi has a
  //     non-PSD Hessian, so the SDFR full Newton tangent is genuinely
  //     indefinite away from the solution.
  // The SDFR tangent mode is selected by the constructor call site of
  // `SDFRRegistration` below, not by a runtime flag. Flip
  // `kUseFullNewton` to instantiate the Newton specialisation (which
  // also requires passing the level-set Hessian as a callable).
  constexpr bool kUseFullNewton = false;

  // -------------------------------------------------------------------------
  // Step 1-3: mesh, classification, attribute tagging.
  // -------------------------------------------------------------------------
  LocalMesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { n, n });
  mesh.scale(h);
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 2);
  mesh.getConnectivity().compute(2, 2);
  mesh.getConnectivity().compute(0, 0);

  const CircleLevelSet levelSet;
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

  std::vector<Index> interfaceFacets;
  interfaceFacets.reserve(classified.cutEdges.size());
  for (const MinSTCut::Edge& edge : classified.cutEdges)
    if (edge.index != MinSTCut::InvalidIndex)
      interfaceFacets.push_back(edge.index);

  constexpr Attribute interiorAttribute = 1;
  constexpr Attribute exteriorAttribute = 2;
  constexpr Attribute interfaceAttribute = 10;
  constexpr Attribute boundaryAttribute = 20;
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

  mesh.save("LevelSetSDRReconstruction_LF.mesh", IO::FileFormat::MFEM);

  // -------------------------------------------------------------------------
  // Step 4: signed distance s_h^LF via Distance::Eikonal (FMM).
  // -------------------------------------------------------------------------
  using ScalarFES = P1<Real, LocalMesh>;
  ScalarFES sLFFes(mesh);
  GridFunction sLF(sLFFes);
  sLF = Real(0);
  Distance::Eikonal<ScalarFES, Math::Vector<Real>>(sLF)
    .setInterface(interfaceAttribute)
    .setInterior(interiorAttribute)
    .solve()
    .sign();

  // -------------------------------------------------------------------------
  // Step 5a: smooth band weight and weighted band measure M_w.
  // -------------------------------------------------------------------------
  const Real delta = 1.75 * h;
  const Real deltaW = 1.5 * delta;
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
      const Real s = sLF.getValue(pt);
      const Real W = std::exp(-s * s / (2 * deltaW * deltaW));
      weightedBandMeasure += W * w;
    }
  }
  if (weightedBandMeasure <= 0)
    throw std::runtime_error("Empty smooth band: M_w = 0.");

  // -------------------------------------------------------------------------
  // Step 5b: P1 vector FES, trial/test/state functions.
  // -------------------------------------------------------------------------
  using VectorFES = P1<Math::SpatialVector<Real>, LocalMesh>;
  VectorFES vectorFes(mesh, 2);
  TrialFunction du(vectorFes);
  TestFunction  v(vectorFes);
  GridFunction  u(vectorFes);
  u.getData().setZero();

  // SDFRIntegratorParameters no longer carry a tangent-mode flag — that's a type
  // property of the SDFRRegistration composite chosen below.
  SDFRIntegratorParameters params;
  params.rhoS = 1;
  params.deltaW = deltaW;
  params.hRef = h;
  params.normalizer = Real(1) / (weightedBandMeasure * h * h);

  // gamma and beta are now form-language scalar fields. Constant
  // weights are wrapped in `RealFunction<Real>(value)` exactly like
  // any other scalar constant in Rodin form language; spatially
  // varying weights can be any `Variational::RealFunctionBase<...>`.
  // gamma weighs the shape term E_shape; beta weighs the singularity
  // floor term E_floor. The floor barrier is active ONLY when
  // j_K^u < jSafe, so beta has no effect away from the singular set.
  RealFunction<Real> barrierGamma(Real(1e-1));
  RealFunction<Real> barrierBeta(Real(1e-2));
  BarrierParameters barrierParams;
  barrierParams.jMin  = 1e-8;
  barrierParams.jSafe = 1e-3;
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

  // -------------------------------------------------------------------------
  // The analytical local barrier Hessian is verified separately by the
  // unit test
  //   tests/unit/Rodin/Adaptation/BarrierLocalHessianFDTest.cpp
  // which exercises the V-shape signature of a correct central-FD check
  // across decreasing eps. We do not repeat that check at runtime here.
  //
  // Step 6: Problem assembly and Newton solve.
  // -------------------------------------------------------------------------
  auto zero = VectorFunction{ Zero(), Zero() };

  // Wrap CircleLevelSet's three derivatives as Rodin form-language
  // function objects. Each takes a `Geometry::Point`; the SDFR
  // integrators construct that Point at the deformed location
  // y = X + u_h(X) per quadrature point and the lambdas read
  // p.getPhysicalCoordinates() to recover y.
  RealFunction phiFn(
      [&](const Geometry::Point& p) -> Real
      {
        const auto& c = p.getPhysicalCoordinates();
        return levelSet.phi(vec2(c(0), c(1)));
      });
  AnalyticVectorFunction gradFn(
      [&](const Geometry::Point& p) -> Math::SpatialVector<Real>
      {
        const auto& c = p.getPhysicalCoordinates();
        return levelSet.grad(vec2(c(0), c(1)));
      },
      /*dimension=*/2);
  AnalyticMatrixFunction hessFn(
      [&](const Geometry::Point& p) -> Math::SpatialMatrix<Real>
      {
        const auto& c = p.getPhysicalCoordinates();
        return levelSet.hess(vec2(c(0), c(1)));
      },
      /*rows=*/2, /*cols=*/2);

  bool solveCompleted = true;
  bool newtonConverged = false;
  std::size_t newtonIterations = 0;
  std::string solveError;
  const bool useSDFRFacade = hasFlag(argc, argv, "use-sdfr-facade");

  if (useSDFRFacade)
  {
    SDFRParameters sdfrParams;
    sdfrParams.rhoS = params.rhoS;
    sdfrParams.deltaW = params.deltaW;
    sdfrParams.hRef = params.hRef;
    sdfrParams.normalizer = params.normalizer;
    sdfrParams.shapeWeight = 1e-1;
    sdfrParams.floorWeight = 1e-2;
    sdfrParams.jMinRatio = barrierParams.jMin;
    sdfrParams.jSafeRatio = barrierParams.jSafe;
    sdfrParams.initialGuess = SDFRInitialGuess::Zero;
    sdfrParams.tangent =
      kUseFullNewton ? SDFRTangent::Newton : SDFRTangent::GaussNewton;
    sdfrParams.maxNewtonIterations = 20;
    sdfrParams.absoluteTolerance = 1e-10;
    sdfrParams.relativeTolerance = 1e-8;

    try
    {
      SDFR sdfr(u);
      const SDFRReport sdfrReport =
        sdfr.setParameters(sdfrParams).solve(sLF, phiFn, gradFn, hessFn);
      newtonConverged = sdfrReport.converged;
      newtonIterations = sdfrReport.iterations;
      solveCompleted = sdfrReport.converged && !sdfrReport.lineSearchFailed;
      if (!solveCompleted)
        solveError = "SDFR facade did not converge";
    }
    catch (const std::exception& ex)
    {
      solveCompleted = false;
      solveError = ex.what();
      std::cout << "\nSDFR facade aborted: " << ex.what() << '\n';
    }
  }
  else
  {
    // The composite captures only the variational skeleton (du, v, u).
    // The data — phi, grad, [hess], sLF, params — is supplied at the
    // operator() call below.
    SDFRRegistration sdfrTerm(du, v, u);

    // Barrier helper, same shape as `SDFRRegistration`:
    // captures (du, v, u) and the per-cell cache + index map at
    // construction; gamma, beta and params arrive at `operator()`.
    JacobianAdmissibilityBarrier barrier(
        du, v, u, cellCache, geomToLocal);

    // We use the DECOMPOSED form (one Tangent + one Residual per helper)
    // because Rodin's form language composes individual integrators, not
    // problem-body fragments. This is the same idiom the Solid examples
    // use: explicit `+ tangent + residual + DirichletBC(...)`. The
    // composite `sdfrTerm()` / `barrier()` operator() forms are convenient
    // when they are the SOLE contribution to a Problem; combining
    // multiple helpers into one Problem goes through Tangent/Residual.
    Problem newton(du, v);
    if constexpr (kUseFullNewton)
    {
      // 4 data args (phi, grad, hess, sLF) + params
      // -> SDFRTangentIntegrator<Newton, ...>
      newton =
            sdfrTerm.Tangent(phiFn, gradFn, hessFn, sLF, params)
          + sdfrTerm.Residual(phiFn, gradFn, sLF, params)
          + barrier.Tangent(barrierGamma, barrierBeta, barrierParams)
          + barrier.Residual(barrierGamma, barrierBeta, barrierParams)
          + DirichletBC(du, zero).on(boundaryAttribute);
    }
    else
    {
      // 3 data args (phi, grad, sLF) + params
      // -> SDFRTangentIntegrator<GaussNewton, ...>
      newton =
            sdfrTerm.Tangent(phiFn, gradFn, sLF, params)
          + sdfrTerm.Residual(phiFn, gradFn, sLF, params)
          + barrier.Tangent(barrierGamma, barrierBeta, barrierParams)
          + barrier.Residual(barrierGamma, barrierBeta, barrierParams)
          + DirichletBC(du, zero).on(boundaryAttribute);
    }

    Solver::SparseLU linearSolver(newton);
    Solver::NewtonSolver solver(linearSolver);
    std::vector<Real> residualHistory;
    solver
      .setMaxIterations(20)
      .setDampingFactor(1.0)
      .setAbsoluteTolerance(1e-10)
      .setRelativeTolerance(1e-8)
      .setMonitor([&](const auto& report)
      {
        const auto bad = barrierInadmissibleCount().exchange(0);
        const auto floorN = barrierFloorActiveCount().exchange(0);
        const Real minJ = barrierMinJ().exchange(
            std::numeric_limits<Real>::infinity(),
            std::memory_order_relaxed);
        std::cout << "Newton " << std::setw(2) << report.iterations
                  << ": ||R||=" << std::scientific << std::setprecision(6)
                  << report.final_residual
                  << "  step=" << report.final_step_norm
                  << "  damping=" << report.damping_factor;
        if (std::isfinite(minJ))
          std::cout << "  min_j=" << std::setprecision(3) << minJ;
        if (floorN > 0)
          std::cout << "  active_floor_cells=" << floorN;
        if (bad > 0)
          std::cout << "  [singular cells=" << bad << "]";
        std::cout << '\n';
        residualHistory.push_back(report.final_residual);
      });

    try
    {
      solver.solve(u);
    }
    catch (const std::exception& ex)
    {
      solveCompleted = false;
      solveError = ex.what();
      std::cout << "\nNewton aborted: " << ex.what() << '\n';
    }
    const auto report = solver.getReport();
    newtonConverged = report.converged;
    newtonIterations = report.iterations;
    printConvergenceOrders(residualHistory);
  }

  // -------------------------------------------------------------------------
  // Step 7: build moved mesh through FES dof mapping (no raw indexing).
  // -------------------------------------------------------------------------
  LocalMesh moved(mesh);
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
  }
  moved.save("LevelSetSDRReconstruction_HF.mesh", IO::FileFormat::MFEM);

  // Interface ||phi(X+u)|| RMS on Gamma_h^LF.
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
  GridFunction phi(p1Fes);
  phi = [&](const Geometry::Point& p) -> Real
  {
    const auto& X = p.getCoordinates();
    return levelSet.phi(vec2(X(0), X(1)));
  };
  phi.setName("phi");
  sLF.setName("s_LF");

  // ----- HF (moved) mesh -----
  ScalarP0 p0FesMoved(moved);
  GridFunction jK(p0FesMoved);
  GridFunction qShape(p0FesMoved);
  GridFunction cellLabelHF(p0FesMoved);
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
      jK.getData()(dof) = sigDetAu / src.Jscale;
      const Real dExp = Real(2) / Real(d);
      qShape.getData()(dof) =
        dst.A.squaredNorm()
          / (static_cast<Real>(d) * std::pow(sigDetAu, dExp));
      cellLabelHF.getData()(dof) =
        static_cast<Real>(classified.labels[local]);
    }
  }
  jK.setName("j");
  qShape.setName("q_shape");
  cellLabelHF.setName("cell_label");

  P1<Real, LocalMesh> p1FesMoved(moved);
  GridFunction phiMoved(p1FesMoved);
  phiMoved = [&](const Geometry::Point& p) -> Real
  {
    const auto& X = p.getCoordinates();
    return levelSet.phi(vec2(X(0), X(1)));
  };
  phiMoved.setName("phi_moved");

  u.setName("displacement");

  IO::XDMF xdmf("LevelSetSDRReconstruction");
  auto lfGrid = xdmf.grid("LF");
  lfGrid.setMesh(mesh);
  lfGrid.add(cellLabel, IO::XDMF::Center::Cell);
  lfGrid.add(phaseMoment, IO::XDMF::Center::Cell);
  lfGrid.add(sigmaKgf, IO::XDMF::Center::Cell);
  lfGrid.add(phi, IO::XDMF::Center::Node);
  lfGrid.add(sLF, IO::XDMF::Center::Node);
  lfGrid.add(u, IO::XDMF::Center::Node);

  auto hfGrid = xdmf.grid("HF");
  hfGrid.setMesh(moved);
  hfGrid.add(cellLabelHF, IO::XDMF::Center::Cell);
  hfGrid.add(jK, IO::XDMF::Center::Cell);
  hfGrid.add(qShape, IO::XDMF::Center::Cell);
  hfGrid.add(phiMoved, IO::XDMF::Center::Node);

  xdmf.write().close();

  // -------------------------------------------------------------------------
  // Diagnostics.
  // -------------------------------------------------------------------------
  std::cout << "\nDiagnostics\n";
  std::cout << "  SDFR tangent mode: "
            << (kUseFullNewton ? sdfrIntegratorTangentModeName(SDFRIntegratorTangentMode::Newton)
                               : sdfrIntegratorTangentModeName(SDFRIntegratorTangentMode::GaussNewton))
            << '\n';
  std::cout << "  cells inside / outside: "
            << classified.insideCells.size() << " / "
            << classified.outsideCells.size() << '\n';
  std::cout << "  interface facets: " << interfaceFacets.size() << '\n';
  std::cout << "  vector dofs: " << u.getData().size() << '\n';
  std::cout << "  h: " << h << ", delta: " << delta
            << ", deltaW: " << deltaW << '\n';
  std::cout << "  |Omega_h|: " << domainMeasure
            << ", M_w (smooth band): " << weightedBandMeasure << '\n';
  std::cout << "  final ||phi(X+u)||_RMS on Gamma_h^LF: "
            << interfacePhiRMS << '\n';
  std::cout << "  solve path: "
            << (useSDFRFacade ? "Adaptation::SDFR" : "manual forms")
            << '\n';
  std::cout << "  Newton iterations: " << newtonIterations
            << ", converged: " << (newtonConverged ? "yes" : "no") << '\n';

  return newtonConverged ? 0 : 1;
}
