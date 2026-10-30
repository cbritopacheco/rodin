/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file NonlocalTest.cpp
 * @brief Manufactured nonlocal solves on isolated MPI communicators.
 */
#include <gtest/gtest.h>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <Rodin/Geometry.h>
#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#include <Rodin/Variational.h>
#include <Rodin/MPI.h>
#include <Rodin/PETSc.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  /// @brief MPI environment owned by main and borrowed by the tests.
  boost::mpi::environment* environment = nullptr;

  /// @brief Preassembled mass-plus-nonlocal operators recover exact constant fields.
  TEST(Rodin_Manufactured_PETSc_MPI_Nonlocal, IsolatedConstantSolution)
  {
    boost::mpi::communicator self(MPI_COMM_SELF, boost::mpi::comm_attach);
    Context::MPI context(*environment, self);
#if defined(PETSC_USE_COMPLEX)
    const PetscScalar exact = PetscCMPLX(2, -3);
    const PetscScalar coefficient = PetscCMPLX(0.125, 0.25);
    const auto load = ComplexFunction(exact * (PetscScalar(1) + coefficient));
#else
    const PetscScalar exact = 2;
    const PetscScalar coefficient = 0.25;
    const auto load = RealFunction(exact * (PetscScalar(1) + coefficient));
#endif
    const auto kernel = [coefficient](const Point&, const Point&) { return coefficient; };
    for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      const size_t dimension = Polytope::Traits(geometry).getDimension();
      for (const size_t nodes : {size_t(2), size_t(3)})
      {
        SCOPED_TRACE(nodes);
        auto local =
          LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, nodes));
        local.scale(Real(1) / (nodes - 1));
        local.getConnectivity().compute(dimension, dimension);
        BalancedCompactPartitioner partitioner(local);
        partitioner.partition(1);
        Sharder<Context::MPI> sharder(context);
        sharder.shard(partitioner);
        sharder.scatter(0);
        auto mesh = sharder.gather(0);
        P1<PetscScalar, Mesh<Context::MPI>> space(mesh);
        PETSc::Variational::TrialFunction u(space);
        PETSc::Variational::TestFunction v(space);
        BilinearForm form(u, v);
        form += Integral(u, v);
        form += Integral(Potential(kernel, u), v);
        form.assemble();
        auto copied = form;
        auto moved = std::move(copied);
        const auto check = [&](auto& candidate) {
          for (size_t repeat = 0; repeat < 2; ++repeat)
          {
            SCOPED_TRACE(repeat);
            candidate.assemble();
            Problem problem(u, v);
            problem = candidate - Integral(load, v);
            problem.assemble();
            PETSc::Solver::KSP solver(problem);
            solver.setType(KSPPREONLY);
            ::PC preconditioner = nullptr;
            ASSERT_EQ(KSPGetPC(solver.getHandle(), &preconditioner), PETSC_SUCCESS);
            ASSERT_EQ(PCSetType(preconditioner, PCLU), PETSC_SUCCESS);
            solver.solve();
            const PetscScalar* values = nullptr;
            ASSERT_EQ(VecGetArrayRead(u.getSolution().getData(), &values), PETSC_SUCCESS);
            for (size_t i = 0; i < space.getSize(); ++i)
            {
              EXPECT_LT(PetscAbsScalar(values[i] - exact), 1e-10) << "DOF=" << i;
            }
            ASSERT_EQ(
              VecRestoreArrayRead(u.getSolution().getData(), &values), PETSC_SUCCESS);
          }
        };
        check(form);
        check(moved);
      }
    }
  }
}

int main(int argc, char** argv)
{
  boost::mpi::environment env(argc, argv);
  environment = &env;
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
