/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
/// @file @brief Two-sided interface assembly into local and distributed PETSc systems.
#include <gtest/gtest.h>
#include "Rodin/PETSc.h"
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"
#ifdef RODIN_USE_MPI
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include "Rodin/MPI/Context/MPI.h"
#include "Rodin/MPI/Geometry/Sharder.h"
#include "Rodin/MPI/Variational.h"
#include "Rodin/Geometry/BalancedCompactPartitioner.h"
#endif

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
#ifdef RODIN_USE_MPI
  boost::mpi::environment* environment;
  boost::mpi::communicator* world;
#endif

  /// @returns A triangular mesh with all face and neighboring-cell incidence available.
  LocalMesh makeMesh()
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 2);
    mesh.getConnectivity().compute(2, 2);
    mesh.getConnectivity().compute(2, 0);
    return mesh;
  }

  /**
   * @brief Compares interface scatter with explicit two-cell constant traces.
   * @tparam MeshType Local or distributed mesh type.
   * @param mesh Mesh with owned and ghost cells where applicable.
   */
  template <class MeshType>
  void checkInterfaceAssembly(const MeshType& mesh)
  {
    P0<PetscScalar, MeshType> fes(mesh);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction v(fes);
    BilinearForm form(u, v);
    LinearForm load(v);
    form = InterfaceIntegral(Jump(u), Jump(v));
    load = InterfaceIntegral(Average(v));
    Math::Matrix<PetscScalar> expected(fes.getSize(), fes.getSize());
    Math::Vector<PetscScalar> expectedLoad(fes.getSize());
    expected.setZero();
    expectedLoad.setZero();
    for (auto face = mesh.getInterface(); face; ++face)
    {
      if constexpr (requires { mesh.getShard(); })
      {
        if (!mesh.getShard().isOwned(1, face->getIndex()))
          continue;
      }
      const auto& incident =
        face->getMesh().getConnectivity().getIncidence({1, 2}, face->getIndex());
      ASSERT_EQ(incident.size(), 2);
      const Index first = fes.getDOFs(2, *incident.begin())[0];
      const Index second = fes.getDOFs(2, *std::next(incident.begin()))[0];
      const Real measure = face->getMeasure();
      expected(first, first) += measure;
      expected(second, second) += measure;
      expected(first, second) -= measure;
      expected(second, first) -= measure;
      expectedLoad[first] += measure / 2;
      expectedLoad[second] += measure / 2;
    }
#ifdef RODIN_USE_MPI
    if constexpr (requires { mesh.getShard(); })
    {
      MPI_Allreduce(MPI_IN_PLACE, expected.data(), static_cast<int>(expected.size()),
        MPIU_SCALAR, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, expectedLoad.data(),
        static_cast<int>(expectedLoad.size()), MPIU_SCALAR, MPI_SUM, MPI_COMM_WORLD);
    }
#endif
    for (size_t repeat = 0; repeat < 2; ++repeat)
    {
      form.assemble();
      load.assemble();
      PetscInt begin, end;
      ASSERT_EQ(MatGetOwnershipRange(form.getOperator(), &begin, &end), PETSC_SUCCESS);
      for (PetscInt row = begin; row < end; ++row)
      {
        PetscScalar value;
        ASSERT_EQ(VecGetValues(load.getVector(), 1, &row, &value), PETSC_SUCCESS);
        EXPECT_LT(PetscAbsScalar(value - expectedLoad[row]), 1e-12);
        for (PetscInt column = 0; column < static_cast<PetscInt>(fes.getSize()); ++column)
        {
          ASSERT_EQ(
            MatGetValues(form.getOperator(), 1, &row, 1, &column, &value), PETSC_SUCCESS);
          EXPECT_LT(PetscAbsScalar(value - expected(row, column)), 1e-12);
        }
      }
    }
  }
}

/// @brief Local PETSc assembly uses both cell supports and safely reassembles the graph.
TEST(PETScInterfaceShape, Local)
{
  const auto mesh = makeMesh();
  checkInterfaceAssembly(mesh);
}

#ifdef RODIN_USE_MPI
/// @brief Owned faces scatter to both incident cells, including off-process matrix entries.
TEST(PETScInterfaceShape, Distributed)
{
  Context::MPI context(*environment, *world);
  Sharder<Context::MPI> sharder(context);
  if (world->rank() == 0)
  {
    const auto mesh = makeMesh();
    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(static_cast<size_t>(world->size()));
    sharder.shard(partitioner);
    sharder.scatter(0);
  }
  const auto mesh = sharder.gather(0);
  checkInterfaceAssembly(mesh);
}
#endif

int main(int argc, char** argv)
{
#ifdef RODIN_USE_MPI
  boost::mpi::environment mpiEnvironment(argc, argv);
  boost::mpi::communicator communicator;
  environment = &mpiEnvironment;
  world = &communicator;
#endif
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  ::testing::InitGoogleTest(&argc, argv);
  const int local = RUN_ALL_TESTS();
  int result = 0;
  MPI_Allreduce(&local, &result, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  PetscFinalize();
  return result;
}
