/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Geometry.h>

#if defined(RODIN_CONSUMER_PETSC)
#include <Rodin/PETSc/Kokkos.h>
#elif defined(RODIN_CONSUMER_MPI)
#include <Rodin/Kokkos/MPI.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#elif defined(RODIN_CONSUMER_KOKKOS)
#include <Rodin/Kokkos.h>
#endif

#if !defined(RODIN_CONSUMER_PETSC) && defined(PETSC_VERSION_MAJOR)
#error "The mesh consumer must not include PETSc."
#endif

/// @brief Exercises installed headers, transitive linking and mesh/assembly construction.
int main(int argc, char** argv)
{
  using namespace Rodin;
#if defined(RODIN_CONSUMER_MPI)
  boost::mpi::environment environment(argc, argv);
  boost::mpi::communicator world;
  Context::MPI distribution(environment, world);
#endif
#if defined(RODIN_CONSUMER_KOKKOS) || defined(RODIN_CONSUMER_MPI) ||                     \
  defined(RODIN_CONSUMER_PETSC)
  ::Kokkos::ScopeGuard runtime(argc, argv);
#if defined(RODIN_CONSUMER_PETSC)
  if (PetscInitialize(&argc, &argv, nullptr, nullptr) != PETSC_SUCCESS)
    return 1;
#endif
  {
    using Backend = ::Kokkos::DefaultExecutionSpace;
#if defined(RODIN_CONSUMER_MPI)
    Geometry::Sharder<Context::MPI> sharder(distribution);
    if (world.rank() == 0)
    {
      auto source = Geometry::Mesh<Context::Local>::UniformGrid(
        Geometry::Polytope::Type::Triangle, {3, 3});
      source.getConnectivity().compute(2, 2);
      source.getConnectivity().compute(1, 2);
      Geometry::BalancedCompactPartitioner partitioner(source);
      partitioner.partition(static_cast<size_t>(world.size()));
      sharder.shard(partitioner);
      sharder.scatter(0);
    }
    const auto host = sharder.gather(0);
    using ContextType = Context::Kokkos<Context::MPI, Backend>;
    Geometry::Mesh<ContextType> mesh(ContextType(distribution), host);
    Variational::P1 space(mesh);
    if (space.getSize() != 9)
      return 1;
#else
    using ContextType = Context::Kokkos<Context::Local, Backend>;
    const auto source = Geometry::Mesh<Context::Local>::UniformGrid(
      Geometry::Polytope::Type::Triangle, {3, 3});
    Geometry::Mesh<ContextType> mesh(ContextType{}, source);
    Variational::P1 space(mesh);
#if defined(RODIN_CONSUMER_PETSC)
    PETSc::Variational::TrialFunction u(space);
    PETSc::Variational::TestFunction v(space);
    PETSc::Variational::Problem problem(u, v);
#else
    Variational::TrialFunction u(space);
    Variational::TestFunction v(space);
    Variational::Problem problem(u, v);
#endif
    problem = Variational::Integral(u, v);
    problem.assemble();
    if (space.getSize() != 9)
      return 1;
#endif
  }
#else
  (void)argc;
  (void)argv;
  const auto mesh = Geometry::Mesh<Context::Local>::UniformGrid(
    Geometry::Polytope::Type::Triangle, {3, 3});
  if (mesh.getCellCount() != 8)
    return 1;
#endif
#if defined(RODIN_CONSUMER_PETSC)
  return PetscFinalize() == PETSC_SUCCESS ? 0 : 1;
#else
  return 0;
#endif
}
