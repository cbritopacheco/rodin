/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Single-field target updates for all supported nodal value ranges. */
#include "SingleFieldTargetAssembly.h"

using namespace Rodin;
using namespace Rodin::Geometry;

namespace
{
  /// @brief Creates the unit interval with distinct endpoint boundary attributes.
  LocalMesh intervalMesh()
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    mesh.scale(Real(0.5));
    mesh.getConnectivity().compute(0, 1);
    for (auto vertex = mesh.getVertex(); vertex; ++vertex)
    {
      const Real x = mesh.getVertexCoordinates(vertex->getIndex()).x();
      if (x == 0)
        mesh.setAttribute(vertex.key(), 11);
      if (x == 1)
        mesh.setAttribute(vertex.key(), 12);
    }
    return mesh;
  }
}

/// @brief Sequential targets preserve exact nonzero boundary data in every value range.
TEST(PETSc_SingleFieldTargets, Sequential)
{
  auto mesh = intervalMesh();
  Rodin::Tests::Unit::checkAllSingleFieldTargets<PETSc::Assembly::Sequential>(
    mesh, PETSC_COMM_SELF);
}

#ifdef RODIN_USE_OPENMP
/// @brief OpenMP targets preserve exact nonzero boundary data in every value range.
TEST(PETSc_SingleFieldTargets, OpenMP)
{
  auto mesh = intervalMesh();
  Rodin::Tests::Unit::checkAllSingleFieldTargets<PETSc::Assembly::OpenMP>(
    mesh, PETSC_COMM_SELF);
}
#endif

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
