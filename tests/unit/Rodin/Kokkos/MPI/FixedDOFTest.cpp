/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Common.h"

namespace Rodin::Tests::Unit
{
  /**
   * @brief Changing fixed values on a reused distribution preserves a constant solution.
   * @param type Cell geometry of the distributed grid.
   *
   * The reaction-diffusion problem with constant source and matching Dirichlet
   * data has that same constant solution. Reassembling with different constants
   * exercises the cached fixed-DOF exchange, including clearing previous values.
   */
  void checkChangedFixedValues(Polytope::Type type)
  {
    using namespace Rodin::Variational;
    Context::MPI ctx(*g_env, *g_world);
    const HostMesh host = distributeFromRoot(ctx, type, 8);
    const DeviceMesh mesh(DeviceContext(ctx), host);
    P1<Real, DeviceMesh> fes(mesh);
    ::Kokkos::View<Real*, TestMemorySpace> prescribed("changed fixed values", 1);
    RealFunction g(
      KOKKOS_LAMBDA(const DeviceGeometry::PointView&)->Real { return prescribed(0); });
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction v(fes);
    Problem problem(u, v);
    problem =
      Integral(Grad(u), Grad(v)) + Integral(u, v) - Integral(g, v) + DirichletBC(u, g);
    for (const Real value : {Real(1), Real(3), Real(0), Real(-2), Real(1)})
    {
      SCOPED_TRACE(::testing::Message() << "prescribed value " << value);
      ::Kokkos::deep_copy(prescribed, value);
      problem.assemble();
      solve(problem);
      const auto x = u.getSolution().getData();
      PetscInt count = 0;
      ASSERT_EQ(VecGetLocalSize(x, &count), PETSC_SUCCESS);
      const PetscScalar* values = nullptr;
      ASSERT_EQ(VecGetArrayRead(x, &values), PETSC_SUCCESS);
      for (PetscInt i = 0; i < count; ++i)
        EXPECT_NEAR(PetscRealPart(values[i]), value, Real(1e-9));
      ASSERT_EQ(VecRestoreArrayRead(x, &values), PETSC_SUCCESS);
    }
  }

  TEST(Kokkos_MPI_FixedValues, Triangles)
  {
    checkChangedFixedValues(Polytope::Type::Triangle);
  }

  TEST(Kokkos_MPI_FixedValues, Quadrilaterals)
  {
    checkChangedFixedValues(Polytope::Type::Quadrilateral);
  }

}
