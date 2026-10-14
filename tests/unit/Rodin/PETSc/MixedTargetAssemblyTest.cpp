/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Exact block and boundary contracts for targeted PETSc assembly. */
#include <gtest/gtest.h>
#include <petsc.h>
#include <algorithm>
#include <array>
#include <functional>
#include <memory>
#include <numeric>
#include <type_traits>
#include <utility>
#include <vector>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/PETSc.h>
#include "MixedTargetAssembly.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  /// @brief Builds the exact interval fixture for local assembly backends.
  template <template <class, class> class Assembler>
  void checkLocalMixedTargets()
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    mesh.scale(Real(0.5));
    mesh.getConnectivity().compute(0, 1);
    P1<PetscScalar> space(mesh);
    ASSERT_EQ(space.getSize(), 3u);
    std::vector<std::pair<Real, Index>> positions;
    for (auto vertex = mesh.getVertex(); vertex; ++vertex)
    {
      const Real x = mesh.getVertexCoordinates(vertex->getIndex()).x();
      positions.emplace_back(x, space.getDOFs(0, vertex->getIndex())[0]);
      if (x == 0)
        mesh.setAttribute(vertex.key(), 11);
      if (x == 1)
        mesh.setAttribute(vertex.key(), 12);
    }
    std::sort(positions.begin(), positions.end());
    Rodin::Tests::Unit::checkMixedTargetAssembly<Assembler>(
      space, positions, PETSC_COMM_SELF);
  }
}

/// @brief Sequential assembly preserves the other target and exact constrained block values.
TEST(PETSc_MixedTargetAssembly, Sequential)
{
  checkLocalMixedTargets<PETSc::Assembly::Sequential>();
}

#ifdef RODIN_USE_OPENMP
/// @brief OpenMP assembly preserves the other target and exact constrained block values.
TEST(PETSc_MixedTargetAssembly, OpenMP)
{
  checkLocalMixedTargets<PETSc::Assembly::OpenMP>();
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
