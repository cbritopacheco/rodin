/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file Seq_NavierStokes.cpp
 * @brief Lid driven cavity, assembled and solved on the device through PETSc.
 *
 * Each time step solves an Oseen problem for a Taylor-Hood pair. The system
 * is assembled into a PETSc matrix of the Kokkos type and solved with a
 * Krylov method, so the iterations run where the matrix lives. The pressure
 * is defined up to a constant, which is given to PETSc as a null space.
 *
 * The first argument is the number of time steps and the second the number
 * of grid nodes per side. Usage, with a Schur complement field split:
 * @code
 *   KokkosSeq_NavierStokes 20 33 -ksp_type fgmres -pc_type fieldsplit
 *     -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type full
 *     -pc_fieldsplit_schur_precondition selfp
 *     -fieldsplit_u_ksp_type preonly -fieldsplit_u_pc_type gamg
 *     -fieldsplit_p_pc_type jacobi -fieldsplit_p_ksp_rtol 1e-2
 *     -ksp_converged_reason
 * @endcode
 * At such sizes the device time is the launch latency of small kernels, so
 * the split keeps few nested iterations: one V-cycle for each velocity solve
 * and a loose pressure solve, the outer FGMRES absorbing both. Solving the
 * velocity to 1e-5 inside each Schur complement product makes about 15000
 * V-cycles and 2 million launches for 5 steps at 149k unknowns, 3.6 times
 * slower on an RTX 5090.
 */
#include <chrono>
#include <cmath>
#include <cstdlib>

#include <Rodin/PETSc.h>
#include <Rodin/Types.h>
#include <Rodin/Solver.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Alert.h>
#include <Rodin/Kokkos/PETSc.h>

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

/**
 * @brief Elapsed milliseconds since a time point.
 * @param t0 Start time.
 * @return Milliseconds since @p t0.
 */
static auto since(std::chrono::high_resolution_clock::time_point t0)
{
  return std::chrono::duration_cast<std::chrono::milliseconds>(
    std::chrono::high_resolution_clock::now() - t0)
    .count();
}

/**
 * @brief Copies a PETSc vector into a host vector.
 * @param[in] x Sequential PETSc vector.
 * @return Host copy of @p x.
 */
static Rodin::Math::Vector<Rodin::Real> toHost(::Vec x)
{
  PetscInt n;
  PetscCallAbort(PETSC_COMM_SELF, VecGetLocalSize(x, &n));
  const PetscScalar* data;
  PetscCallAbort(PETSC_COMM_SELF, VecGetArrayRead(x, &data));
  Rodin::Math::Vector<Rodin::Real> res =
    Eigen::Map<const Rodin::Math::Vector<Rodin::Real>>(data, n);
  PetscCallAbort(PETSC_COMM_SELF, VecRestoreArrayRead(x, &data));
  return res;
}

int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Geometry;
  using namespace Rodin::Variational;

  ::Kokkos::ScopeGuard runtime(argc, argv);
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);
  {
    const size_t steps = argc > 1 && argv[1][0] != '-' ? std::stoul(argv[1]) : 20;
    const size_t nodes = argc > 2 && argv[2][0] != '-' ? std::stoul(argv[2]) : 33;
    constexpr Attribute lid = 2;
    constexpr Attribute walls = 3;
    constexpr Real boundaryTolerance = 1e-12; // Coordinate tolerance on this unit grid.
    Mesh host = Mesh().UniformGrid(Polytope::Type::Triangle, {nodes, nodes});
    host.scale(1.0 / (nodes - 1));
    host.getConnectivity().compute(1, 2);
    for (auto edge = host.getFace(); !edge.end(); ++edge)
    {
      if (!edge->isBoundary())
        continue;
      bool onTop = true;
      for (auto vertex = edge->getVertex(); !vertex.end(); ++vertex)
        onTop = onTop && std::abs(vertex->getCoordinates().y() - 1) < boundaryTolerance;
      host.setAttribute({1, edge->getIndex()}, onTop ? lid : walls);
    }
    const Mesh<DeviceContext> mesh(DeviceContext{}, host);

    H1 velocitySpace(std::integral_constant<size_t, 2>{}, mesh, 2);
    H1 pressureSpace(std::integral_constant<size_t, 1>{}, mesh);
    PETSc::Variational::TrialFunction u(velocitySpace);
    PETSc::Variational::TrialFunction p(pressureSpace);
    PETSc::Variational::TestFunction v(velocitySpace);
    PETSc::Variational::TestFunction q(pressureSpace);
    u.setName("u");
    p.setName("p");
    Alert::Info() << "Velocity and pressure degrees of freedom: "
                  << velocitySpace.getSize() << " + " << pressureSpace.getSize()
                  << Alert::Raise;

    constexpr Real peakLidSpeed = 2;
    constexpr Real viscosity = 0.01; // Reynolds number 200 for unit cavity width.
    constexpr Real dt = 0.01;
    constexpr Real rampTime = 1;
    GridFunction velocity(velocitySpace);
    velocity = Zero(2);

    size_t iterations = 0;
    long assembly = 0, solve = 0;
    for (size_t step = 1; step <= steps; ++step)
    {
      const Real time = step * dt;
      const Real speed = peakLidSpeed * std::min<Real>(1, time / rampTime);
      VectorFunction lidVelocity(
        KOKKOS_LAMBDA(const Rodin::Kokkos::Geometry::PointView& point)->Real {
          const Real offset = (point.x() - 0.5) / 0.25;
          const Real profile = offset * offset < 1 ? 1 - offset * offset : Real(0);
          return speed * profile * profile;
        },
        Real(0));
      Problem flow(u, p, v, q);
      flow = (1 / dt) * Integral(u, v) + Integral(Dot(Mult(Jacobian(u), velocity), v)) +
        viscosity * Integral(Jacobian(u), Jacobian(v)) - Integral(p, Div(v)) +
        Integral(Div(u), q) - (1 / dt) * Integral(velocity, v) +
        DirichletBC(u, lidVelocity).on(lid) + DirichletBC(u, Zero(2)).on(walls);

      auto t0 = std::chrono::high_resolution_clock::now();
      flow.assemble();
      flow.setFieldSplits();
      assembly += since(t0);

      // The constant pressure spans the null space of the operator.
      ::Mat& A = flow.getLinearSystem().getOperator();
      ::Vec constant;
      PetscCallAbort(PETSC_COMM_SELF, MatCreateVecs(A, &constant, PETSC_NULLPTR));
      PetscCallAbort(PETSC_COMM_SELF, VecSet(constant, 0));
      for (PetscInt i = velocitySpace.getSize();
        i < static_cast<PetscInt>(velocitySpace.getSize() + pressureSpace.getSize()); ++i)
        PetscCallAbort(PETSC_COMM_SELF, VecSetValue(constant, i, 1, INSERT_VALUES));
      PetscCallAbort(PETSC_COMM_SELF, VecAssemblyBegin(constant));
      PetscCallAbort(PETSC_COMM_SELF, VecAssemblyEnd(constant));
      PetscCallAbort(PETSC_COMM_SELF, VecNormalize(constant, PETSC_NULLPTR));
      MatNullSpace nullSpace;
      PetscCallAbort(PETSC_COMM_SELF,
        MatNullSpaceCreate(PETSC_COMM_SELF, PETSC_FALSE, 1, &constant, &nullSpace));
      PetscCallAbort(PETSC_COMM_SELF, MatSetNullSpace(A, nullSpace));
      PetscCallAbort(PETSC_COMM_SELF, MatNullSpaceDestroy(&nullSpace));
      PetscCallAbort(PETSC_COMM_SELF, VecDestroy(&constant));

      t0 = std::chrono::high_resolution_clock::now();
      Solver::GMRES solver(flow);
      solver.setTolerances(1e-8, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT).solve();
      solve += since(t0);
      iterations += solver.getIterationNumber();
      velocity.setData(toHost(u.getSolution().getData()));
      Alert::Info() << "Step " << step << " / " << steps << ", "
                    << solver.getIterationNumber() << " iterations." << Alert::Raise;
    }
    Alert::Success() << "Assembly " << assembly << " ms, solve " << solve << " ms, "
                     << iterations << " iterations in " << steps << " steps."
                     << Alert::Raise;
  }
  PetscFinalize();
  return 0;
}
