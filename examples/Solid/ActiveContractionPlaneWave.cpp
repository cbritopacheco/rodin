/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ActiveContractionPlaneWave.cpp
 * @brief Dynamic Hill-Maxwell active fiber contraction driven by an electrical
 *        activation plane wave.
 *
 * A unit square is clamped on its top edge (homogeneous Dirichlet) and is
 * otherwise free.  A passive NeoHookean response is augmented with the active
 * fiber stress of @ref Solid::ActiveContraction with a horizontal fiber
 * direction @f$ n_f = (1, 0) @f$.
 *
 * The driver is the *electrical activation* @f$ u_1(y, t) @f$, modelled as a
 * smooth rectangular pulse traveling upward from the bottom of the domain:
 *
 *   u_1(y, t) = U * exp( -((y - c t) / w)^2 ) ,    U > 0.
 *
 * The active extension @f$ e_c @f$ and the Hill-Maxwell internal variables
 * @f$ \gamma, \beta @f$ live at every quadrature point and are advanced by the
 * per-quadrature-point local Newton inside @ref Solid::ActiveContraction
 * during the global Newton iteration.  After global convergence at each
 * time step we sweep the quadrature points and commit the converged state to
 * the per-cell, per-QP buffer used as the previous state for the next step.
 *
 * Output is written to XDMF for visualization in ParaView (apply "Warp by
 * Vector" on Displacement to see the contraction).
 *
 * The same code works in 2D or 3D: Rodin's Solid module determines the spatial
 * dimension at runtime from the mesh.
 */
#include <cmath>
#include <cstddef>
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
  // Per-cell, per-quadrature-point Hill-Maxwell state at time n.
  struct LocalState
  {
    Real ec    = 0.0;   // active extension e_c^n
    Real gamma = 0.0;   // sqrt(k_c^n)
    Real beta  = 0.0;   // tau_c^n / gamma^n
  };

  using StateBuffer = std::vector<std::vector<LocalState>>;
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

  // ---- finite-element space -----------------------------------------------
  const size_t dim = mesh.getSpaceDimension();
  P1 Vh(mesh, dim);

  // ---- materials ----------------------------------------------------------
  const Real E      = 50.0;
  const Real nu     = 0.3;
  const Real lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  const Real mu     = E / (2.0 * (1.0 + nu));
  Solid::NeoHookean passive(lambda, mu);

  Solid::ActiveFiberLaw::Parameters activeParams;
  activeParams.stiffness            = 200.0;  // E_s
  activeParams.damping              = 0.5;    // mu (Hill-Maxwell viscous damping)
  activeParams.destructionRate      = 0.4;    // alpha
  activeParams.crossBridgeStiffness = 100.0;  // k_0
  activeParams.contractility        = 80.0;   // sigma_0
  activeParams.initial.extension    = 0.0;
  activeParams.initial.stiffness    = 0.0;
  activeParams.initial.stress       = 0.0;
  Solid::ActiveFiberLaw activeLaw(activeParams);

  Solid::ActiveContraction<Solid::NeoHookean, Solid::ActiveFiberLaw>
    law(passive, activeLaw);
  law.setLocalTolerance(1e-12).setLocalMaxIterations(50);

  // ---- electrical activation pulse ----------------------------------------
  // u_1(y, t) = U * exp( -((y - c t - y0) / w)^2 ).  Pulse travels upward.
  const Real activationAmplitude = 1.5;
  const Real waveSpeed           = 1.0;
  const Real waveWidth           = 0.05;
  const Real waveStart           = -0.3;

  Real currentTime = 0.0;
  const Real dtStep = 0.02;
  const size_t nSteps = static_cast<size_t>((1.0 - waveStart + 0.4) / (waveSpeed * dtStep));

  auto activationAt = [&](Real y, Real t) -> Real
  {
    const Real arg = (y - (waveStart + waveSpeed * t)) / waveWidth;
    return activationAmplitude * std::exp(-arg * arg);
  };

  // ---- per-quadrature-point state buffer ----------------------------------
  StateBuffer state(mesh.getCellCount());

  auto ensureCell = [&](Index cellIdx, size_t nqp)
  {
    if (cellIdx >= state.size())
      state.resize(cellIdx + 1);
    if (state[cellIdx].size() != nqp)
      state[cellIdx].assign(nqp, LocalState{});
  };

  // ---- input lambda used by both integrators and commit sweep -------------
  auto activeInput = [&](Solid::ConstitutivePoint& cp)
  {
    const Index cellIdx = cp.get<Solid::Tags::CellIndex>();
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

  // ---- solution -----------------------------------------------------------
  GridFunction u(Vh);
  u.setName("Displacement");
  u = VectorFunction{ Zero(), Zero() };

  // ---- XDMF output --------------------------------------------------------
  IO::XDMF xdmf("ActiveContractionPlaneWave");
  auto grid = xdmf.grid();
  grid.setMesh(mesh);
  grid.add(u);
  xdmf.write(0.0);

  // ---- commit sweep: re-evaluate cache at converged u and store new state -
  auto commitState = [&]()
  {
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
      ensureCell(idx, nqp);

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
      }
    }
  };

  // ---- time loop ----------------------------------------------------------
  TrialFunction du(Vh);
  TestFunction  v(Vh);
  auto zero = VectorFunction{ Zero(), Zero() };

  for (size_t step = 1; step <= nSteps; ++step)
  {
    currentTime = step * dtStep;

    Solid::MaterialTangent tangent(law, du, v);
    tangent.setDisplacement(u);
    tangent.setInput(activeInput);

    Solid::InternalForce residual(law, v);
    residual.setDisplacement(u);
    residual.setInput(activeInput);

    Problem newton(du, v);
    newton = tangent
           + residual
           + DirichletBC(du, zero).on(topBC);

    SparseLU linearSolver(newton);
    NewtonSolver solver(linearSolver);
    solver.setMaxIterations(50)
          .setAbsoluteTolerance(1e-9)
          .setRelativeTolerance(1e-7);

    solver.solve(u);

    commitState();
    xdmf.write(currentTime).flush();
  }

  xdmf.close();
  return 0;
}
