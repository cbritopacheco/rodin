/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Analytic and coefficient oracles for PETSc MFEM persistence. */
#include <gtest/gtest.h>
#include <petsc.h>
#include <memory>
#include <sstream>
#include <type_traits>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/PETSc.h>
#include <Rodin/PETSc/IO.h>
#include "MFEMTest.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::IO;

namespace
{
  template <size_t Family, size_t Rank>
  void checkExport(Polytope::Type geometry)
  {
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    auto mesh = LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 2));
    for (size_t from = 0; from <= dimension; ++from)
    {
      for (size_t to = 0; to <= dimension; ++to)
      {
        if (from != to)
          mesh.getConnectivity().compute(from, to);
      }
    }
    MFEMTest::checkExport<Family, Rank>(mesh);
  }

  class PETScMFEM : public ::testing::TestWithParam<Polytope::Type>
  {};

  TEST_P(PETScMFEM, CubicScalar)
  {
    checkExport<3, 0>(GetParam());
  }

  TEST_P(PETScMFEM, CubicVector)
  {
    checkExport<3, 1>(GetParam());
  }

  TEST_P(PETScMFEM, CubicMatrix)
  {
    checkExport<3, 2>(GetParam());
  }
  TEST_P(PETScMFEM, LinearScalar)
  {
    checkExport<1, 0>(GetParam());
  }
  TEST_P(PETScMFEM, LinearVector)
  {
    checkExport<1, 1>(GetParam());
  }
  TEST_P(PETScMFEM, LinearMatrix)
  {
    checkExport<1, 2>(GetParam());
  }
  TEST_P(PETScMFEM, ConstantScalar)
  {
    checkExport<0, 0>(GetParam());
  }
  TEST_P(PETScMFEM, ConstantVector)
  {
    checkExport<0, 1>(GetParam());
  }
  TEST_P(PETScMFEM, ConstantMatrix)
  {
    checkExport<0, 2>(GetParam());
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, PETScMFEM,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge));
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
