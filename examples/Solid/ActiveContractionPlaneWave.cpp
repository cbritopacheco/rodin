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
 */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <vector>

#include <Rodin/Geometry.h>
#include <Rodin/Assembly.h>
#include <Rodin/Variational.h>
#include <Rodin/Solid.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Solver/NewtonSolver.h>
#include <Rodin/Solver/SparseLU.h>
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
}

int main(int, char**)
{
  // ---- geometry -----------------------------------------------------------
  constexpr size_t nc = 33;
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, { nc, nc });
  mesh.scale(1.0 / static_cast<Real>(nc - 1));
  mesh.getConnectivity().compute(1, 2);

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
      mesh.setAttribute({ 1, it->getIndex() }, topBC);
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
  activeParams.stiffness            = 200.0;
  activeParams.damping              = 0.5;
  activeParams.destructionRate      = 0.4;
  activeParams.crossBridgeStiffness = 100.0;
  activeParams.contractility        = 80.0;
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
  wave.amplitude      = 20;
  wave.speed          = 1.0;
  wave.width          = 0.05;
  wave.start          = -0.3;
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
  const size_t nSteps = static_cast<size_t>(tEnd / dtStep) + 1;

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

    Solid::MaterialTangent tangent(law, du, v);
    tangent.setDisplacement(u);
    tangent.setInput(activeInput);

    Solid::InternalForce residual(law, v);
    residual.setDisplacement(u);
    residual.setInput(activeInput);

    Problem newton(du, v);
    newton =
             (rho * bdfA0 / dtStep) * Integral(du, v)
           + tangent
           + (rho / dtStep)
             * Integral(bdfA0 * u + bdfA1 * uPrevious
                        + bdfA2 * uPreviousPrevious, v)
           + residual
           + DirichletBC(du, zero).on(topBC);

    SparseLU linearSolver(newton);
    NewtonSolver solver(linearSolver);
    solver.setMaxIterations(50)
          .setAbsoluteTolerance(1e-9)
          .setRelativeTolerance(1e-7);
    Real previousNewtonResidual = std::numeric_limits<Real>::quiet_NaN();
    solver.setMonitor([&](const auto& report)
    {
      std::cout
        << "  Newton " << std::setw(2) << report.iterations
        << " residual = " << report.final_residual
        << " step = " << report.final_step_norm
        << " damping = " << report.damping_factor;
      if (std::isfinite(previousNewtonResidual)
          && previousNewtonResidual > 0.0)
      {
        std::cout
          << " ratio = "
          << report.final_residual / previousNewtonResidual;
      }
      if (report.converged)
        std::cout << " converged";
      std::cout << std::endl;
      previousNewtonResidual = report.final_residual;
    });

    solver.solve(u);

    const auto diagnostics = commitState();
    const auto& report = solver.getReport();

    const auto bdfData =
      bdfA0 * u.getData()
      + bdfA1 * uPrevious.getData()
      + bdfA2 * uPreviousPrevious.getData();

    std::cout
      << "  summary:"
      << " converged = " << (report.converged ? "yes" : "no")
      << " iterations = " << report.iterations
      << " final_residual = " << report.final_residual
      << " |u| = " << u.getData().norm()
      << " |BDF(u)| = " << bdfData.norm()
      << " ec = [" << diagnostics.minActiveExtension
      << ", " << diagnostics.maxActiveExtension << "]"
      << " max|gamma| = " << diagnostics.maxGamma
      << " max|beta| = " << diagnostics.maxBeta
      << " max_local_newton = " << diagnostics.maxLocalIterations
      << std::endl;

    uPreviousPrevious = uPrevious;
    uPrevious = u;

    // Update activation field from current time and mesh node positions
    activationField = [&](const Geometry::Point& p)
    {
      return activationAt(p.getPhysicalCoordinates()(1), currentTime);
    };

    xdmf.write(currentTime).flush();
  }

  xdmf.close();
  return 0;
}
