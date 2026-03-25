/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file DampedTransientCantileverBeam.cpp
 * @brief Transient hyperelastic cantilever beam with inertia, damping, and release.
 *
 * A rectangular 2D beam (4 × 1) is clamped on the left edge and subject to
 * a downward traction on the right edge. The traction ramps up linearly over
 * the first half of the simulation and is then suddenly removed.
 *
 * The beam is modeled as a compressible NeoHookean solid with inertia and
 * simple mass-proportional viscous damping:
 *
 *   rho * u_tt + c * u_t - div(P(F(u))) = 0
 *
 * Time integration is performed with the Newmark-beta method using the
 * average-acceleration choice beta = 1/4, gamma = 1/2.
 *
 * Output is written to XDMF for visualization in ParaView.
 *
 * This version also instruments a discrete derivative check for the assembled
 * nonlinear residual R(u) and tangent J(u):
 *
 *   FD(eps) = (R(u + eps * eta) - R(u)) / eps
 *
 * and compares FD(eps) against J(u) * eta.
 *
 * If the tangent is consistent, ||FD(eps) - J(u) eta|| should scale like O(eps)
 * for sufficiently small eps, until roundoff dominates.
 */
#include "Rodin/IO/ForwardDecls.h"
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>

#include <Rodin/Geometry.h>
#include <Rodin/Assembly.h>
#include <Rodin/Variational.h>
#include <Rodin/Solid.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Solver/NewtonSolver.h>
#include <Rodin/Solver/SparseLU.h>
#include <Rodin/Solver/CG.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

namespace
{
  template <class VectorType>
  VectorType copyVector(const VectorType& x)
  {
    VectorType y;
    y.resize(x.size());
    for (Index i = 0; i < x.size(); ++i)
      y(i) = x(i);
    return y;
  }

  template <class VectorType>
  VectorType makeZeroLike(const VectorType& x)
  {
    VectorType y;
    y.resize(x.size());
    y.setZero();
    return y;
  }

  template <class VectorType>
  VectorType addScaled(const VectorType& x, Real alpha, const VectorType& y)
  {
    VectorType z;
    z.resize(x.size());
    for (Index i = 0; i < x.size(); ++i)
      z(i) = x(i) + alpha * y(i);
    return z;
  }

  template <class VectorType>
  VectorType scaledDifference(const VectorType& x, const VectorType& y, Real alpha)
  {
    VectorType z;
    z.resize(x.size());
    for (Index i = 0; i < x.size(); ++i)
      z(i) = alpha * (x(i) - y(i));
    return z;
  }

  template <class VectorType>
  Real infNorm(const VectorType& x)
  {
    Real n = 0.0;
    for (Index i = 0; i < x.size(); ++i)
      n = std::max(n, std::abs(x(i)));
    return n;
  }

  template <class VectorType>
  VectorType makeDeterministicPerturbation(const VectorType& x)
  {
    VectorType eta;
    eta.resize(x.size());
    for (Index i = 0; i < x.size(); ++i)
    {
      const Real ii = static_cast<Real>(i + 1);
      eta(i) = std::sin(0.173 * ii) + 0.5 * std::cos(0.071 * ii);
    }
    return eta;
  }

  template <class VectorType>
  void applyHomogeneousDirichletToVector(VectorType& x, const IndexMap<Real>& dbc)
  {
    for (const auto& [local, value] : dbc)
    {
      (void) value;
      x(local) = 0.0;
    }
  }

  template <class GridFunctionType>
  void applyHomogeneousDirichletToGridFunction(GridFunctionType& gf, const IndexMap<Real>& dbc)
  {
    auto& data = gf.getData();
    applyHomogeneousDirichletToVector(data, dbc);
  }

  template <class ProblemType, class GridFunctionType>
  void assembleResidualAndJacobian(ProblemType& pb, GridFunctionType& u)
  {
    (void) u;
    pb.assemble();
  }

  template <class ProblemType, class GridFunctionType>
  void checkDiscreteDerivative(
      ProblemType& pb,
      GridFunctionType& u,
      const IndexMap<Real>& dbc)
  {
    using VectorType = std::decay_t<decltype(u.getData())>;

    auto& x = u.getData();

    // Base assembly at current state
    assembleResidualAndJacobian(pb, u);
    const auto& A = pb.getLinearSystem().getOperator();
    const auto& RbaseRef = pb.getLinearSystem().getVector();

    VectorType Rbase = copyVector(RbaseRef);

    VectorType eta = makeDeterministicPerturbation(x);
    applyHomogeneousDirichletToVector(eta, dbc);

    const VectorType Jeta = A * eta;

    std::cout << "\n============================================================\n";
    std::cout << "Discrete derivative check\n";
    std::cout << "  |R(u)|_2      = " << Rbase.norm() << '\n';
    std::cout << "  |R(u)|_inf    = " << infNorm(Rbase) << '\n';
    std::cout << "  |eta|_2       = " << eta.norm() << '\n';
    std::cout << "  |J(u) eta|_2  = " << Jeta.norm() << '\n';
    std::cout << "  |J(u) eta|_inf= " << infNorm(Jeta) << '\n';

    const std::array<Real, 8> epsilons = {
      1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8
    };

    const VectorType xBase = copyVector(x);

    for (const Real eps : epsilons)
    {
      x = addScaled(xBase, eps, eta);
      applyHomogeneousDirichletToVector(x, dbc);

      assembleResidualAndJacobian(pb, u);
      const auto& RepsRef = pb.getLinearSystem().getVector();
      VectorType Reps = copyVector(RepsRef);

      VectorType FD = scaledDifference(Reps, Rbase, 1.0 / eps);
      VectorType err = scaledDifference(FD, Jeta, 1.0);
      for (Index i = 0; i < err.size(); ++i)
        err(i) = FD(i) + Jeta(i);

      const Real nFD2   = FD.norm();
      const Real nErr2  = err.norm();
      const Real nFDInf = infNorm(FD);
      const Real nErrInf= infNorm(err);
      const Real rel2   = (nFD2   > 0.0) ? (nErr2   / nFD2)   : nErr2;
      const Real relInf = (nFDInf > 0.0) ? (nErrInf / nFDInf) : nErrInf;

      std::cout
        << "  eps = " << eps
        << "  |FD|_2 = " << nFD2
        << "  |FD + Jeta|_2 = " << nErr2
        << "  rel2 = " << rel2
        << "  |FD|_inf = " << nFDInf
        << "  |FD + Jeta|_inf = " << nErrInf
        << "  relInf = " << relInf
        << '\n';
    }

    x = xBase;
    assembleResidualAndJacobian(pb, u);
    std::cout << "============================================================\n\n";
  }
}

int main(int, char**)
{
  // ---- geometry -----------------------------------------------------------
  constexpr size_t nx = 66;
  constexpr size_t ny = 18;
  constexpr Real Lx = static_cast<Real>(nx - 1) / static_cast<Real>(ny - 1);

  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, { nx, ny });
  mesh.scale(1.0 / static_cast<Real>(ny - 1));
  mesh.getConnectivity().compute(1, 2);

  // ---- label boundary edges -----------------------------------------------
  constexpr Attribute leftBC  = 1;
  constexpr Attribute rightBC = 2;
  constexpr Real eps = 1e-10;

  for (auto it = mesh.getBoundary(); !it.end(); ++it)
  {
    const auto& verts = it->getVertices();
    const size_t nv = verts.size();

    Real xSum = 0;
    for (size_t i = 0; i < nv; ++i)
      xSum += mesh.getVertexCoordinates(verts[i])(0);
    const Real xMid = xSum / static_cast<Real>(nv);

    if (xMid < eps)
      mesh.setAttribute({ 1, it->getIndex() }, leftBC);
    else if (xMid > Lx - eps)
      mesh.setAttribute({ 1, it->getIndex() }, rightBC);
  }

  mesh.save("CantileverBeam.mesh", IO::FileFormat::MEDIT);

  // ---- Finite-element space -----------------------------------------------
  const size_t dim = mesh.getSpaceDimension();
  P1 Vh(mesh, dim);

  // ---- material -----------------------------------------------------------
  const Real E  = 200.0;
  const Real nu = 0.3;
  const Real lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  const Real mu     = E / (2.0 * (1.0 + nu));
  Solid::NeoHookean law(lambda, mu);

  // ---- dynamics parameters ------------------------------------------------
  const Real rho = 1.0;
  const Real c   = 0.5;   // mass-proportional viscous damping coefficient
  const Real dt  = 1e-2;
  const size_t nSteps = 2000;

  // Newmark average acceleration
  const Real beta  = 0.25;
  const Real gamma = 0.5;

  // Effective coefficients
  const Real aMass = rho / (beta * dt * dt);
  const Real aDamp = c * gamma / (beta * dt);

  // ---- fields -------------------------------------------------------------
  GridFunction u(Vh);      // displacement
  GridFunction vel(Vh);    // velocity
  GridFunction acc(Vh);    // acceleration

  u.setName("Displacement");
  vel.setName("Velocity");
  acc.setName("Acceleration");

  auto zero = VectorFunction{ Zero(), Zero() };
  u   = zero;
  vel = zero;
  acc = zero;

  GridFunction uPred(Vh);
  GridFunction vPred(Vh);
  GridFunction aNew(Vh);
  GridFunction rhsDamp(Vh);

  uPred.setName("PredictedDisplacement");
  vPred.setName("PredictedVelocity");
  aNew.setName("AccelerationNew");
  rhsDamp.setName("DampingRHS");

  // Build homogeneous Dirichlet map once, for perturbation filtering
  TrialFunction uBC(Vh);
  DirichletBC dbc(uBC, zero);
  dbc.on(leftBC);
  dbc.assemble();
  const IndexMap<Real> dbcMap = dbc.getDOFs();

  // ---- output -------------------------------------------------------------
  IO::XDMF xdmf("CantileverBeam");
  auto grid = xdmf.grid();
  grid.setMesh(mesh);
  grid.add(u);
  grid.add(vel);
  grid.add(acc);

  xdmf.write(0.0);

  // ---- time loop ----------------------------------------------------------
  constexpr Real maxTraction = -1.0; // downward

  TrialFunction du(Vh);
  TestFunction  w(Vh);

  for (size_t step = 1; step <= nSteps; ++step)
  {
    const Real t = step * dt;

    Real ty = 0.0;
    if (t > 0.5 && t < 1)
      ty = maxTraction;

    auto traction = VectorFunction{ Zero(), RealFunction(ty) };

    // ---- Newmark predictors -----------------------------------------------
    uPred = u + dt * vel + (dt * dt * (0.5 - beta)) * acc;
    vPred = vel + (dt * (1.0 - gamma)) * acc;

    // Damping RHS term:
    // rhsDamp = vPred - gamma/(beta dt) * uPred
    rhsDamp = vPred - (gamma / (beta * dt)) * uPred;

    // Use predictor as initial guess for u^{n+1}
    u = uPred;

    // ---- nonlinear solid operators ----------------------------------------
    Solid::MaterialTangent tangent(law, du, w);
    tangent.setLinearizationPoint(u);

    Solid::InternalForce internal(law, w);
    internal.setLinearizationPoint(u);

    // Effective nonlinear problem:
    //
    //   R(u; w)
    //   = F_int(u; w)
    //   + aMass * (u, w)
    //   + aDamp * (u, w)
    //   - F_ext(w)
    //   - aMass * (uPred, w)
    //   + c * (rhsDamp, w)
    //
    // Newton increment equation:
    //
    //   J(u^k)[du] + R(u^k) = 0
    //
    Problem newton(du, w);
    newton =
          tangent
           + aMass * Integral(du, w)
           + aDamp * Integral(du, w)
           + internal
           + aMass * Integral(u, w)
           + aDamp * Integral(u, w)
           - BoundaryIntegral(traction, w).over(rightBC)
           - aMass * Integral(uPred, w)
           + c * Integral(rhsDamp, w)
           + DirichletBC(du, zero).on(leftBC);

    // Discrete derivative check once, at the first time step and first Newton base point
    if (step == 1)
      checkDiscreteDerivative(newton, u, dbcMap);

    SparseLU linearSolver(newton);
    NewtonSolver solver(linearSolver);
    solver.setMaxIterations(50)
          .setDampingFactor(1.0)
          .setAbsoluteTolerance(1e-10)
          .setRelativeTolerance(1e-8);
    solver.solve(u);

    std::cout << "Step " << step << ", time " << t << std::endl;

    // ---- Newmark correctors -----------------------------------------------
    aNew = (1.0 / (beta * dt * dt)) * (u - uPred);
    vel  = vPred + (gamma * dt) * aNew;
    acc  = aNew;

    xdmf.write(t).flush();
  }

  xdmf.close();
  return 0;
}
