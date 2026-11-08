/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file CantileverBeam.cpp
 * @brief Hyperelastic cantilever beam with Newmark time stepping, assembled on
 * the device.
 *
 * Device port of examples/Solid/CantileverBeam.cpp. The Newmark predictors
 * and correctors act on the coefficient vectors, and the Neo-Hookean tangent
 * and residual are assembled by a device kernel at each Newton iteration,
 * and solved on the device by cuDSS or RocSolverRF when Rodin has them.
 * The discrete derivative check of the host example is kept.
 */
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>

#include <Rodin/Geometry.h>
#include <Rodin/Assembly.h>
#include <Rodin/Variational.h>
#include <Rodin/Solid.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Solver/NewtonSolver.h>
#include <Rodin/Solver.h>
#include <Rodin/Kokkos.h>

#if defined(KOKKOS_ENABLE_CUDA)
using Backend = Rodin::Context::CUDA;

#elif defined(KOKKOS_ENABLE_HIP)
using Backend = Rodin::Context::HIP;
#elif defined(KOKKOS_ENABLE_SYCL)
using Backend = Rodin::Context::SYCL;
#else
using Backend = ::Kokkos::DefaultHostExecutionSpace;
#endif

// A global using directive on Rodin would make Kokkos ambiguous in the Kokkos headers.
using DeviceContext = Rodin::Context::Kokkos<Rodin::Context::Local, Backend>;

namespace
{
  using Rodin::Index;
  using Rodin::IndexMap;
  using Rodin::Real;

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
      (void)value;
      x(local) = 0.0;
    }
  }

  template <class GridFunctionType>
  void applyHomogeneousDirichletToGridFunction(
    GridFunctionType& gf, const IndexMap<Real>& dbc)
  {
    auto& data = gf.getData();
    applyHomogeneousDirichletToVector(data, dbc);
  }

  template <class ProblemType, class GridFunctionType>
  void assembleResidualAndJacobian(ProblemType& pb, GridFunctionType& u)
  {
    (void)u;
    pb.assemble();
  }

  template <class ProblemType, class GridFunctionType>
  void checkDiscreteDerivative(
    ProblemType& pb, GridFunctionType& u, const IndexMap<Real>& dbc)
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

    const std::array<Real, 8> epsilons = {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8};

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

      const Real nFD2 = FD.norm();
      const Real nErr2 = err.norm();
      const Real nFDInf = infNorm(FD);
      const Real nErrInf = infNorm(err);
      const Real rel2 = (nFD2 > 0.0) ? (nErr2 / nFD2) : nErr2;
      const Real relInf = (nFDInf > 0.0) ? (nErrInf / nFDInf) : nErrInf;

      std::cout << "  eps = " << eps << "  |FD|_2 = " << nFD2
                << "  |FD + Jeta|_2 = " << nErr2 << "  rel2 = " << rel2
                << "  |FD|_inf = " << nFDInf << "  |FD + Jeta|_inf = " << nErrInf
                << "  relInf = " << relInf << '\n';
    }

    x = xBase;
    assembleResidualAndJacobian(pb, u);
    std::cout << "============================================================\n\n";
  }
}

int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Geometry;
  using namespace Rodin::Variational;

  ::Kokkos::ScopeGuard runtime(argc, argv);

  constexpr size_t nx = 66;
  constexpr size_t ny = 18;
  constexpr Real Lx = static_cast<Real>(nx - 1) / static_cast<Real>(ny - 1);

  Mesh host = Mesh().UniformGrid(Polytope::Type::Triangle, {nx, ny});
  host.scale(1.0 / static_cast<Real>(ny - 1));
  host.getConnectivity().compute(1, 2);

  constexpr Attribute leftBC = 1;
  constexpr Attribute rightBC = 2;
  constexpr Real eps = 1e-10;
  for (auto it = host.getBoundary(); !it.end(); ++it)
  {
    const auto& verts = it->getVertices();
    Real xSum = 0;
    for (size_t i = 0; i < verts.size(); ++i)
      xSum += host.getVertexCoordinates(verts[i])(0);
    const Real xMid = xSum / static_cast<Real>(verts.size());
    if (xMid < eps)
      host.setAttribute({1, it->getIndex()}, leftBC);
    else if (xMid > Lx - eps)
      host.setAttribute({1, it->getIndex()}, rightBC);
  }
  const Mesh<DeviceContext> mesh(DeviceContext{}, host);

  P1 vh(mesh, mesh.getSpaceDimension());

  const Real E = 200.0;
  const Real nu = 0.3;
  const Real lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  const Real mu = E / (2.0 * (1.0 + nu));
  Solid::NeoHookean law(lambda, mu);

  const Real rho = 1.0;
  const Real c = 0.5; // mass-proportional viscous damping coefficient
  const Real dt = 1e-2;
  const size_t nSteps = 2000;

  // Newmark average acceleration
  const Real beta = 0.25;
  const Real gamma = 0.5;
  const Real aMass = rho / (beta * dt * dt);
  const Real aDamp = c * gamma / (beta * dt);

  GridFunction u(vh);
  GridFunction vel(vh);
  GridFunction acc(vh);
  GridFunction uPred(vh);
  GridFunction vPred(vh);
  GridFunction rhsDamp(vh);
  u.getData().setZero();
  vel.getData().setZero();
  acc.getData().setZero();

  // Homogeneous Dirichlet map, for perturbation filtering
  TrialFunction uBC(vh);
  DirichletBC dbc(uBC, Zero(2));
  dbc.on(leftBC);
  dbc.assemble();
  const IndexMap<Real> dbcMap = std::get<IndexMap<Real>>(dbc.getDOFs());

  IO::XDMF xdmf("CantileverBeam");
  xdmf.grid()
    .setMesh(mesh)
    .add("Displacement", u)
    .add("Velocity", vel)
    .add("Acceleration", acc);
  xdmf.write(0.0).flush();

  constexpr Real maxTraction = -1.0;

  TrialFunction du(vh);
  TestFunction w(vh);
  auto ivw = Solid::InternalVirtualWork(law, u);

  // The problem and the solvers outlive the time loop: the sparsity pattern
  // is the same at every Newton iteration of every step, so cuDSS and
  // RocSolverRF keep their symbolic analysis and SparseLU its ordering, and
  // each iteration only refactorizes the tangent. With a device solver the
  // tangent stays in the memory of the backend, where the solver reads it
  // in place.
#if defined(RODIN_USE_CUDSS) || defined(RODIN_USE_ROCSOLVER)
  using System =
    Math::LinearSystem<Rodin::Kokkos::SparseMatrix<Real, Backend>, Math::Vector<Real>>;
  Problem<System, decltype(du), decltype(w)> newton(du, w);
#ifdef RODIN_USE_CUDSS
  Solver::CuDSS linearSolver(newton);
#else
  Solver::RocSolverRF linearSolver(newton);
#endif
#else
  Problem newton(du, w);
  Solver::SparseLU linearSolver(newton);
#endif
  Solver::NewtonSolver solver(linearSolver);
  solver.setMaxIterations(50)
    .setDampingFactor(1.0)
    .setAbsoluteTolerance(1e-10)
    .setRelativeTolerance(1e-8);

  const auto define = [&](auto& pb, Real ty) {
    pb = ivw(du, w) + aMass * Integral(du, w) + aDamp * Integral(du, w) +
      aMass * Integral(u, w) + aDamp * Integral(u, w) -
      BoundaryIntegral(VectorFunction{0, ty}, w).over(rightBC) -
      aMass * Integral(uPred, w) + c * Integral(rhsDamp, w) +
      DirichletBC(du, Zero(2)).on(leftBC);
  };

  for (size_t step = 1; step <= nSteps; ++step)
  {
    const Real t = step * dt;
    const Real ty = (t > 0.5 && t < 1) ? maxTraction : 0.0;

    // Newmark predictors
    uPred.getData() =
      u.getData() + dt * vel.getData() + (dt * dt * (0.5 - beta)) * acc.getData();
    vPred.getData() = vel.getData() + (dt * (1.0 - gamma)) * acc.getData();
    rhsDamp.getData() = vPred.getData() - (gamma / (beta * dt)) * uPred.getData();
    u.getData() = uPred.getData();

    if (step == 1)
    {
      // The check multiplies the tangent on the host.
      Problem check(du, w);
      define(check, ty);
      checkDiscreteDerivative(check, u, dbcMap);
    }

    define(newton, ty);
    solver.solve(u);
    if (!linearSolver.success())
      return EXIT_FAILURE;

    std::cout << "Step " << step << ", time " << t << std::endl;

    // Newmark correctors
    acc.getData() = (1.0 / (beta * dt * dt)) * (u.getData() - uPred.getData());
    vel.getData() = vPred.getData() + (gamma * dt) * acc.getData();

    xdmf.write(t).flush();
  }

  xdmf.close();
  return 0;
}
