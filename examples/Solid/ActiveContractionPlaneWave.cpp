/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ActiveContractionPlaneWave.cpp
 * @brief Dynamic Hill-Maxwell active fiber contraction driven by repeated
 *        electrical activation plane waves with relaxation gaps.
 *
 * A unit square is clamped on its top edge (homogeneous Dirichlet) and is
 * otherwise free.  A passive NeoHookean response is augmented with the active
 * fiber stress of @ref Solid::ActiveContraction with a horizontal fiber
 * direction @f$ n_f = (1, 0) @f$.
 *
 * The driver is the *electrical activation* @f$ u_1(y, t) @f$, modelled as
 * repeated Gaussian plane waves traveling upward from the bottom of the domain:
 *
 *   u_1(y, t) = U * exp( -((y - (y_0 + c*(t - t_k))) / w)^2 ) ,    U > 0.
 *
 * The launches @f$t_k@f$ are separated by a relaxation gap so the active
 * extension @f$ e_c @f$ has time to decay before the second wave arrives.
 *
 * The displacement solve is advanced with an implicit dynamic relaxation term:
 *
 *   rho * BDF2(u) + div(P(u, e_c)) = 0.
 *
 * The first time step uses backward Euler, then the scheme switches to BDF2
 * once two previous displacement states are available.
 *
 * Output fields written to XDMF:
 *   - Displacement (transient vector) - warp by vector in ParaView
 *   - Activation   (transient scalar) - Gaussian travelling plane wave
 *   - FiberDirection (static vector)  - constant (1,0) fiber field
 *
 * Optional command-line arguments:
 *   ActiveContractionPlaneWave [nc=64] [nSteps=full]
 *                                [targetMaxDisplacement=0]
 *                                [activeScale=1]
 */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <vector>

#include <Rodin/Geometry.h>
#include <Rodin/Assembly.h>
#include <Rodin/Variational.h>
#include <Rodin/Solid.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Solver/NewtonSolver.h>
#include <Rodin/Solver/CG.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

namespace
{
  struct LocalState
  {
    Real ec    = 0.0;
    Real gamma = 0.0;
    Real beta  = 0.0;
  };

  using StateBuffer = std::vector<std::vector<LocalState>>;

  struct StepDiagnostics
  {
    Real minActiveExtension = std::numeric_limits<Real>::infinity();
    Real maxActiveExtension = -std::numeric_limits<Real>::infinity();
    Real maxGamma = 0.0;
    Real maxBeta = 0.0;
    size_t maxLocalIterations = 0;
  };

  size_t parseSize(const char* value, const char* name)
  {
    char* end = nullptr;
    const auto parsed = std::strtoull(value, &end, 10);
    if (end == value || *end != '\0')
    {
      std::cerr << "Invalid " << name << ": " << value << std::endl;
      std::exit(2);
    }
    return static_cast<size_t>(parsed);
  }

  Real parseReal(const char* value, const char* name)
  {
    char* end = nullptr;
    const Real parsed = std::strtod(value, &end);
    if (end == value || *end != '\0')
    {
      std::cerr << "Invalid " << name << ": " << value << std::endl;
      std::exit(2);
    }
    return parsed;
  }
}

int main(int argc, char** argv)
{
  size_t nc = 64;
  size_t requestedSteps = 0;
  Real targetMaxDisplacement = 0.0;
  Real activeScale = 1.0;

  if (argc > 1)
    nc = parseSize(argv[1], "mesh resolution");
  if (argc > 2)
    requestedSteps = parseSize(argv[2], "step count");
  if (argc > 3)
    targetMaxDisplacement = parseReal(argv[3], "target max displacement");
  if (argc > 4)
    activeScale = parseReal(argv[4], "active scale");
  if (argc > 5)
  {
    std::cerr << "Usage: " << argv[0] << " [nc=64] [nSteps=full]"
              << " [targetMaxDisplacement=0] [activeScale=1]" << std::endl;
    return 2;
  }
  if (nc < 2)
  {
    std::cerr << "Mesh resolution must be at least 2." << std::endl;
    return 2;
  }
  if (targetMaxDisplacement < 0.0)
  {
    std::cerr << "Target max displacement must be nonnegative." << std::endl;
    return 2;
  }
  if (activeScale <= 0.0)
  {
    std::cerr << "Active scale must be positive." << std::endl;
    return 2;
  }

  // ---- geometry -----------------------------------------------------------
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, { nc, nc, nc });
  mesh.scale(1.0 / static_cast<Real>(nc - 1));
  mesh.getConnectivity().compute(2, 3);
  mesh.save("miaow.mesh", IO::FileFormat::MEDIT);

  std::cout << "Mesh: " << mesh.getVertexCount() << " vertices, "
            << mesh.getCellCount() << " cells." << std::endl;

  // ---- label top boundary -------------------------------------------------
  constexpr Attribute topBC = 1;
  constexpr Real eps = 1e-10;

  for (auto it = mesh.getBoundary(); !it.end(); ++it)
  {
    const auto& verts = it->getVertices();
    const size_t nv = verts.size();
    Real ySum = 0;
    for (size_t i = 0; i < nv; ++i)
      ySum += mesh.getVertexCoordinates(verts[i])(1);
    const Real yMid = ySum / static_cast<Real>(nv);
    if (yMid > 1.0 - eps)
      mesh.setAttribute({ 2, it->getIndex() }, topBC);
  }

  // ---- finite-element spaces ----------------------------------------------
  const size_t dim = mesh.getSpaceDimension();
  P1 Vh(mesh, dim);   // vector (displacement, fiber direction)
  P1 Sh(mesh);        // scalar (electrical activation)

  // ---- materials ----------------------------------------------------------
  const Real E      = 50.0;
  const Real nu     = 0.3;
  const Real lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  const Real mu     = E / (2.0 * (1.0 + nu));
  const Real rho    = 1.0;
  Solid::NeoHookean passive(lambda, mu);

  Solid::ActiveFiberLaw::Parameters activeParams;
  activeParams.stiffness = activeScale * 20.0;
  activeParams.damping              = 300;
  activeParams.destructionRate      = 0.4;
  activeParams.crossBridgeStiffness = activeScale * 100.0;
  activeParams.contractility = activeScale * 20.0;
  activeParams.initial.extension    = 0.0;
  activeParams.initial.stiffness    = 0.0;
  activeParams.initial.stress       = 0.0;
  Solid::ActiveFiberLaw activeLaw(activeParams);

  Solid::ActiveContraction<Solid::NeoHookean, Solid::ActiveFiberLaw>
    law(passive, activeLaw);
  law.setLocalTolerance(1e-14).setLocalMaxIterations(80);

  // ---- repeated electrical activation plane waves -------------------------
  // u_1(y, t) = U * exp(-((y - (y0 + c*(t - tk))) / w)^2).
  // The launch interval is intentionally longer than the domain traversal
  // time so there is a quiet relaxation gap between waves.
  struct PlaneWaveParams
  {
    Real amplitude;
    Real speed;
    Real width;
    Real start;
    Real gap;
  };

  const size_t nWaves = 2;

  PlaneWaveParams wave;
  wave.amplitude = activeScale * 20.0;
  wave.speed          = 0.5;
  wave.width          = 0.05;
  wave.start          = 0;
  wave.gap            = 1.5;

  auto activationAt = [&](Real y, Real t) -> Real
  {
    Real activation = 0.0;
    for (size_t k = 0; k < nWaves; ++k)
    {
      const Real localTime = t - static_cast<Real>(k) * wave.gap;
      if (localTime < 0.0)
        continue;

      const Real center = wave.start + wave.speed * localTime;
      const Real arg = (y - center) / wave.width;
      activation = std::max(activation, wave.amplitude * std::exp(-arg * arg));
    }
    return activation - 0.5 * wave.amplitude;   // offset so activation is
                                                // negative at rest
  };

  const Real   dtStep = 0.02;
  const Real traversalTime = (1.0 - wave.start) / wave.speed;
  const Real tEnd =
    static_cast<Real>(nWaves - 1) * wave.gap
    + traversalTime
    + 0.4 / wave.speed;
  const size_t fullSteps = static_cast<size_t>(tEnd / dtStep) + 1;
  const size_t nSteps =
    requestedSteps == 0 ? fullSteps : std::min(requestedSteps, fullSteps);

  Real currentTime = 0.0;

  // ---- per-quadrature-point state buffer ----------------------------------
  StateBuffer state(mesh.getCellCount());

  auto activeInput = [&](Solid::ConstitutivePoint& cp)
  {
    const Index cellIdx   = cp.get<Solid::Tags::CellIndex>();
    const std::size_t qpIdx = cp.get<Solid::Tags::QuadraturePointIndex>();

    if (cellIdx >= state.size())
      state.resize(cellIdx + 1);
    if (qpIdx >= state[cellIdx].size())
      state[cellIdx].resize(qpIdx + 1, LocalState{});
    const LocalState& s = state[cellIdx][qpIdx];

    Math::SpatialVector<Real> n_f(static_cast<std::uint8_t>(dim));
    n_f.setZero();
    n_f[0] = 1.0;
    cp.set<Solid::Tags::FiberDirection>(n_f);

    Real y = 0.0;
    if (cp.getPoint())
      y = cp.getPoint()->get().getPhysicalCoordinates()(1);

    cp.set<Solid::Tags::TimeStep>(dtStep);
    cp.set<Solid::Tags::PreviousActiveExtension>(s.ec);
    cp.set<Solid::Tags::PreviousActiveGamma>(s.gamma);
    cp.set<Solid::Tags::PreviousActiveBeta>(s.beta);
    cp.set<Solid::Tags::ElectricalActivation>(activationAt(y, currentTime));
  };

  // ---- solution fields ----------------------------------------------------
  GridFunction u(Vh);
  u.setName("Displacement");
  u = VectorFunction{ Zero(), Zero() };

  GridFunction uPrevious(Vh);
  uPrevious = u;

  GridFunction uPreviousPrevious(Vh);
  uPreviousPrevious = u;

  GridFunction activationField(Sh);
  activationField.setName("Activation");
  activationField = [&](const Geometry::Point& p)
  {
    return activationAt(p.getPhysicalCoordinates()(1), currentTime);
  };

  GridFunction fiberField(Vh);
  fiberField.setName("FiberDirection");
  fiberField = VectorFunction{ RealFunction(1.0), Zero() };

  auto computeMaxNodalDisplacement = [&]() -> Real {
    Real maxDisplacement = 0.0;
    const auto& uData = u.getData();
    for (auto it = mesh.getVertex(); it; ++it)
    {
      const auto& dofs = Vh.getDOFs(0, it->getIndex());
      Real squaredNorm = 0.0;
      for (size_t c = 0; c < dim && c < static_cast<size_t>(dofs.size()); ++c)
      {
        const Real component = uData(dofs(static_cast<Index>(c)));
        squaredNorm += component * component;
      }
      maxDisplacement = std::max(maxDisplacement, std::sqrt(squaredNorm));
    }
    return maxDisplacement;
  };

  // ---- XDMF output --------------------------------------------------------
  IO::XDMF xdmf("ActiveContractionPlaneWave");
  auto grid = xdmf.grid();
  grid.setMesh(mesh);
  grid.add(u);
  grid.add(activationField);
  grid.add(fiberField,
           IO::XDMF::Center::Node,
           IO::XDMF::AttributePolicy::Static);   // written once, reused
  xdmf.write(0.0);

  // ---- commit sweep -------------------------------------------------------
  auto commitState = [&]()
  {
    StepDiagnostics diagnostics;

    const auto& uData = u.getData();
    for (auto it = mesh.getCell(); it; ++it)
    {
      const auto& polytope = *it;
      const Index idx = polytope.getIndex();
      const auto& fe = Vh.getFiniteElement(dim, idx);
      const size_t ndof = fe.getCount();
      const size_t order = 2 * fe.getOrder();
      const auto& qf = QF::PolytopeQuadratureFormula::get(order, polytope.getGeometry());
      const auto& quadrature = polytope.getQuadrature(qf);
      const size_t nqp = quadrature.getSize();

      if (idx >= state.size())
        state.resize(idx + 1);
      if (state[idx].size() != nqp)
        state[idx].assign(nqp, LocalState{});

      for (size_t q = 0; q < nqp; ++q)
      {
        const auto& pt = quadrature.getPoint(q);
        const auto& rc = qf.getPoint(q);
        const auto& JacInv = pt.getJacobianInverse();

        Math::SpatialMatrix<Real> H(static_cast<std::uint8_t>(dim),
                                    static_cast<std::uint8_t>(dim));
        H.setZero();
        for (size_t dof = 0; dof < ndof; ++dof)
        {
          Math::SpatialMatrix<Real> refJac =
            fe.getBasis(dof).getJacobian()(rc);
          Math::SpatialMatrix<Real> physJac = refJac * JacInv;
          const Real u_dof = uData(Vh.getGlobalIndex({dim, idx}, dof));
          for (size_t c = 0; c < dim; ++c)
            for (size_t k = 0; k < dim; ++k)
              H(c, k) += u_dof * physJac(c, k);
        }

        Solid::KinematicState ks(dim);
        ks.setDisplacementGradient(H);
        Solid::ConstitutivePoint cp(pt, ks);
        cp.set<Solid::Tags::CellIndex>(idx);
        cp.set<Solid::Tags::QuadraturePointIndex>(q);
        activeInput(cp);

        typename decltype(law)::Cache cache;
        law.setCache(cache, cp);

        state[idx][q].ec    = cache.activeExtension;
        state[idx][q].gamma = cache.newState.gamma;
        state[idx][q].beta  = cache.newState.beta;

        diagnostics.minActiveExtension =
          std::min(diagnostics.minActiveExtension, cache.activeExtension);
        diagnostics.maxActiveExtension =
          std::max(diagnostics.maxActiveExtension, cache.activeExtension);
        diagnostics.maxGamma =
          std::max(diagnostics.maxGamma, std::abs(cache.newState.gamma));
        diagnostics.maxBeta =
          std::max(diagnostics.maxBeta, std::abs(cache.newState.beta));
        diagnostics.maxLocalIterations =
          std::max(diagnostics.maxLocalIterations, cache.localIterations);
      }
    }

    return diagnostics;
  };

  // ---- time loop ----------------------------------------------------------
  TrialFunction du(Vh);
  TestFunction  v(Vh);
  auto zero = VectorFunction{ Zero(), Zero() };

  std::cout << std::scientific << std::setprecision(6);
  std::cout << "Time steps: " << nSteps;
  if (nSteps != fullSteps)
    std::cout << " (truncated from " << fullSteps << ")";
  std::cout << std::endl;
  std::cout << "Active scale: " << activeScale;
  if (targetMaxDisplacement > 0.0)
    std::cout << ", target max displacement: " << targetMaxDisplacement;
  std::cout << std::endl;

  Real peakMaxNodalDisplacement = 0.0;

  for (size_t step = 1; step <= nSteps; ++step)
  {
    currentTime = step * dtStep;

    const bool useBDF2 = step > 1;
    const Real bdfA0 = useBDF2 ? 1.5 : 1.0;
    const Real bdfA1 = useBDF2 ? -2.0 : -1.0;
    const Real bdfA2 = useBDF2 ? 0.5 : 0.0;

    if (useBDF2)
      u.getData() = 2.0 * uPrevious.getData() - uPreviousPrevious.getData();
    else
      u.getData() = uPrevious.getData();


    Real activationMin = std::numeric_limits<Real>::infinity();
    Real activationMax = -std::numeric_limits<Real>::infinity();
    for (auto it = mesh.getVertex(); it; ++it)
    {
      const Real a = activationAt(it->getCoordinates()(1), currentTime);
      activationMin = std::min(activationMin, a);
      activationMax = std::max(activationMax, a);
    }

    std::cout
      << "\n[step " << step << " / " << nSteps << "]"
      << " t = " << currentTime
      << " scheme = " << (useBDF2 ? "BDF2" : "BE")
      << " activation = [" << activationMin << ", " << activationMax << "]"
      << " |u_guess| = " << u.getData().norm()
      << std::endl;

    auto ivw = Solid::InternalVirtualWork(law, u).setInput(activeInput);

    Problem newton(du, v);
    newton =
             ivw(du, v)
           + (rho * bdfA0 / dtStep) * Integral(du, v)
           + (rho / dtStep)
             * Integral(bdfA0 * u + bdfA1 * uPrevious
                        + bdfA2 * uPreviousPrevious, v)
           + DirichletBC(du, zero).on(topBC);

    CG linearSolver(newton);
    linearSolver.setTolerance(1e-8).setMaxIterations(2000);
    NewtonSolver solver(linearSolver);
    solver.setMaxIterations(50)
          .setAbsoluteTolerance(1e-9)
          .setRelativeTolerance(1e-7);
    Real previousNewtonResidual = std::numeric_limits<Real>::quiet_NaN();
    solver.setMonitor([&](const auto& report) {
      std::cout << "  Newton " << std::setw(2) << report.iterations
                << " residual = " << report.finalResidual
                << " step = " << report.finalStepNorm
                << " damping = " << report.dampingFactor;
      if (std::isfinite(previousNewtonResidual)
          && previousNewtonResidual > 0.0)
      {
        std::cout << " ratio = " << report.finalResidual / previousNewtonResidual;
      }
      if (report.converged)
        std::cout << " converged";
      std::cout << std::endl;
      previousNewtonResidual = report.finalResidual;
    });

    solver.solve(u);

    if (!solver.getReport().converged)
    {
      std::cerr << "Newton solver failed to converge at step " << step
                << " (t = " << currentTime << ")" << std::endl;
      return 1;
    }

    const auto diagnostics = commitState();
    const auto& report = solver.getReport();
    const Real maxNodalDisplacement = computeMaxNodalDisplacement();
    peakMaxNodalDisplacement = std::max(peakMaxNodalDisplacement, maxNodalDisplacement);

    const auto bdfData =
      bdfA0 * u.getData()
      + bdfA1 * uPrevious.getData()
      + bdfA2 * uPreviousPrevious.getData();

    std::cout << "  summary:"
              << " converged = " << (report.converged ? "yes" : "no")
              << " iterations = " << report.iterations
              << " finalResidual = " << report.finalResidual
              << " |u| = " << u.getData().norm()
              << " max_node|u| = " << maxNodalDisplacement
              << " peak_max_node|u| = " << peakMaxNodalDisplacement
              << " |BDF(u)| = " << bdfData.norm() << " ec = ["
              << diagnostics.minActiveExtension << ", " << diagnostics.maxActiveExtension
              << "]"
              << " max|gamma| = " << diagnostics.maxGamma
              << " max|beta| = " << diagnostics.maxBeta
              << " max_local_newton = " << diagnostics.maxLocalIterations;
    if (targetMaxDisplacement > 0.0 && peakMaxNodalDisplacement > 0.0)
    {
      const Real suggestedScale =
        activeScale * targetMaxDisplacement / peakMaxNodalDisplacement;
      std::cout << " suggested_active_scale = " << suggestedScale;
    }
    std::cout << std::endl;

    uPreviousPrevious = uPrevious;
    uPrevious = u;

    // Update activation field from current time and mesh node positions
    activationField = [&](const Geometry::Point& p)
    {
      return activationAt(p.getPhysicalCoordinates()(1), currentTime);
    };

    xdmf.write(currentTime).flush();
  }

  if (targetMaxDisplacement > 0.0 && peakMaxNodalDisplacement > 0.0)
  {
    const Real suggestedScale =
      activeScale * targetMaxDisplacement / peakMaxNodalDisplacement;
    std::cout << "\nCalibration:"
              << " target_max_node|u| = " << targetMaxDisplacement
              << " observed_peak_max_node|u| = " << peakMaxNodalDisplacement
              << " active_scale = " << activeScale
              << " suggested_active_scale = " << suggestedScale << std::endl;
  }

  xdmf.close();
  return 0;
}
