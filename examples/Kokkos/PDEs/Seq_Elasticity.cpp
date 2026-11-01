/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file Seq_Elasticity.cpp
 * @brief Cantilever in linear elasticity, assembled and solved on the device.
 *
 * The system is assembled into a PETSc matrix of the Kokkos type, so the
 * Krylov iterations and the preconditioner run where the matrix lives. The
 * first argument refines the 4 x 1 beam: the grid has 40 n x 10 n cells.
 *
 * Usage, with an algebraic multigrid that runs on the device:
 * @code
 *   KokkosSeq_Elasticity 8 -pc_type gamg -ksp_converged_reason
 *   KokkosSeq_Elasticity 8 -pc_type hypre -ksp_converged_reason
 * @endcode
 */
#include <array>
#include <chrono>
#include <cmath>
#include <numeric>
#include <vector>

#include <Rodin/PETSc.h>
#include <Rodin/Types.h>
#include <Rodin/Solver.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Alert.h>
#include <Rodin/Kokkos/PETSc.h>

#ifdef KOKKOS_ENABLE_CUDA
using Backend = Rodin::Context::CUDA;
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

int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Geometry;
  using namespace Rodin::Variational;

  ::Kokkos::ScopeGuard runtime(argc, argv);
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);
  {
    const size_t n = argc > 1 && argv[1][0] != '-' ? std::stoul(argv[1]) : 1;
    constexpr Real length = 4;
    constexpr Real boundaryTolerance = 1e-12; // Coordinate tolerance for this grid.
    constexpr Attribute fixed = 2;
    constexpr Attribute loaded = 3;
    constexpr Attribute free = 4;

    Mesh host = Mesh().UniformGrid(Polytope::Type::Triangle, {40 * n + 1, 10 * n + 1});
    host.scale(0.1 / n);
    host.getConnectivity().compute(1, 2);
    for (auto face = host.getFace(); !face.end(); ++face)
    {
      if (!face->isBoundary())
        continue;
      bool onLeft = true;
      bool onRight = true;
      for (auto vertex = face->getVertex(); !vertex.end(); ++vertex)
      {
        const Real x = vertex->getCoordinates().x();
        onLeft = onLeft && std::abs(x) < boundaryTolerance;
        onRight = onRight && std::abs(x - length) < boundaryTolerance;
      }
      host.setAttribute({1, face->getIndex()}, onLeft ? fixed : onRight ? loaded : free);
    }
    Alert::Info() << "Number of vertices in the mesh: " << host.getVertexCount()
                  << Alert::Raise;

    const Mesh<DeviceContext> mesh(DeviceContext{}, host);
    P1 vh(mesh, mesh.getSpaceDimension());
    PETSc::Variational::TrialFunction u(vh);
    PETSc::Variational::TestFunction v(vh);

    constexpr Real lambda = 0.5769;
    constexpr Real mu = 0.3846;
    VectorFunction traction{0, -1e-4};

    Problem elasticity(u, v);
    elasticity = LinearElasticityIntegral(u, v)(lambda, mu) -
      BoundaryIntegral(traction, v).over(loaded) + DirichletBC(u, Zero(2)).on(fixed);

    auto t0 = std::chrono::high_resolution_clock::now();
    elasticity.assemble();
    Alert::Success() << "Assembled the linear system on the device in " << since(t0)
                     << " ms." << Alert::Raise;

    // An algebraic multigrid needs the rigid body modes to coarsen elasticity.
    // The components are stored one after the other, so the modes are built
    // from the vertex coordinates instead of MatNullSpaceCreateRigidBody.
    ::Mat& A = elasticity.getLinearSystem().getOperator();
    const size_t vertices = host.getVertexCount();
    Math::Vector<Real> xs(vertices), ys(vertices);
    for (size_t i = 0; i < vertices; ++i)
    {
      xs[i] = host.getVertexCoordinates(i).x();
      ys[i] = host.getVertexCoordinates(i).y();
    }
    xs.array() -= xs.mean();
    ys.array() -= ys.mean();
    std::array<::Vec, 3> modes;
    std::vector<PetscInt> rows(2 * vertices);
    std::iota(rows.begin(), rows.end(), PetscInt(0));
    for (size_t k = 0; k < modes.size(); ++k)
    {
      Math::Vector<Real> mode = Math::Vector<Real>::Zero(2 * vertices);
      if (k < 2)
        mode.segment(k * vertices, vertices).setOnes();
      else
      {
        mode.head(vertices) = -ys;
        mode.tail(vertices) = xs;
      }
      mode.normalize();
      PetscCallAbort(PETSC_COMM_SELF, MatCreateVecs(A, &modes[k], PETSC_NULLPTR));
      PetscCallAbort(PETSC_COMM_SELF,
        VecSetValues(modes[k], rows.size(), rows.data(), mode.data(), INSERT_VALUES));
      PetscCallAbort(PETSC_COMM_SELF, VecAssemblyBegin(modes[k]));
      PetscCallAbort(PETSC_COMM_SELF, VecAssemblyEnd(modes[k]));
    }
    MatNullSpace nearNullSpace;
    PetscCallAbort(PETSC_COMM_SELF,
      MatNullSpaceCreate(
        PETSC_COMM_SELF, PETSC_FALSE, modes.size(), modes.data(), &nearNullSpace));
    PetscCallAbort(PETSC_COMM_SELF, MatSetNearNullSpace(A, nearNullSpace));
    PetscCallAbort(PETSC_COMM_SELF, MatNullSpaceDestroy(&nearNullSpace));
    for (auto& mode : modes)
      PetscCallAbort(PETSC_COMM_SELF, VecDestroy(&mode));

    t0 = std::chrono::high_resolution_clock::now();
    Solver::CG solver(elasticity);
    solver.setTolerances(1e-10, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT).solve();
    Alert::Success() << "Solved the linear system in " << since(t0) << " ms and "
                     << solver.getIterationNumber() << " iterations." << Alert::Raise;

    u.getSolution().save("Elasticity.gf", IO::FileFormat::MFEM);
    host.save("Elasticity.mesh", IO::FileFormat::MFEM);
  }
  PetscFinalize();
  return 0;
}
