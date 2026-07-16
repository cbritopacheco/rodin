/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <petsc.h>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/PETSc/Variational/GridFunction.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  /// @brief Verifies sequential operator bracket read write for PET sc grid function by checking tolerance-based numerical results.
  TEST(PETSc_GridFunction, SequentialOperatorBracketReadWrite)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    Rodin::PETSc::Variational::GridFunction gf(fes);

    for (Index i = 0; i < gf.getSize(); ++i)
      gf[i] = static_cast<PetscScalar>(1.0 + 0.25 * static_cast<Real>(i));
    gf.flush();

    const auto& cgf = gf;
    for (Index i = 0; i < gf.getSize(); ++i)
    {
      EXPECT_DOUBLE_EQ(
        static_cast<double>(PetscRealPart(cgf[i])),
        1.0 + 0.25 * static_cast<double>(i));
    }
    cgf.flush();
  }

  /// @brief Verifies sequential vec set keeps operator bracket readable for PET sc grid function by checking tolerance-based numerical results.
  TEST(PETSc_GridFunction, SequentialVecSetKeepsOperatorBracketReadable)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Quadrilateral, { 3, 3 });
    P1 fes(mesh);
    Rodin::PETSc::Variational::GridFunction gf(fes);

    gf = static_cast<PetscScalar>(7.0);

    const auto& cgf = gf;
    for (Index i = 0; i < gf.getSize(); ++i)
      EXPECT_DOUBLE_EQ(static_cast<double>(PetscRealPart(cgf[i])), 7.0);
    cgf.flush();
  }

  TEST(PETSc_GridFunction, SequentialPointEvaluationMatchesMappedBasisExpansion)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Tetrahedron, {2, 2, 2});
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    mesh.getConnectivity().compute(2, 0);
    mesh.getConnectivity().compute(3, 0);
    P1 fes(mesh, 3);
    Rodin::PETSc::Variational::GridFunction gf(fes);

    for (Index i = 0; i < gf.getSize(); ++i)
      gf[i] = static_cast<PetscScalar>(Real(-0.4) + Real(0.07) * i);
    gf.flush();

    const auto cell = mesh.getCell(0);
    const Point p(*cell, Math::SpatialPoint{0.17, 0.23, 0.19});
    const auto& fe = fes.getFiniteElement(3, cell->getIndex());
    const auto& dofs = fes.getDOFs(3, cell->getIndex());
    const auto& cgf = gf;
    Math::SpatialVector<PetscScalar> expected;
    for (size_t local = 0; local < fe.getCount(); ++local)
    {
      const auto basis = fes.getPushforward({3, cell->getIndex()}, fe.getBasis(local));
      const auto term = cgf[dofs[local]] * basis(p);
      if (local == 0)
        expected = term;
      else
        expected += term;
    }

    const auto value = cgf(p);
    EXPECT_NEAR((value - expected).norm(), 0, 1e-14);
    cgf.flush();
  }

  /// @brief Verifies sequential min returns value and index for PET sc grid function by checking tolerance-based numerical results, exact expected values.
  TEST(PETSc_GridFunction, SequentialMinReturnsValueAndIndex)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {4, 4});
    P1 fes(mesh);
    Rodin::PETSc::Variational::GridFunction gf(fes);

    for (Index i = 0; i < gf.getSize(); ++i)
      gf[i] = static_cast<PetscScalar>(10.0 + static_cast<Real>(i));

    constexpr Index expectedIdx = 3;
    gf[expectedIdx] = static_cast<PetscScalar>(-2.5);
    gf.flush();

    Index idx = gf.getSize();
    const auto value = gf.min(idx);
    EXPECT_EQ(idx, expectedIdx);
    EXPECT_DOUBLE_EQ(static_cast<double>(PetscRealPart(value)), -2.5);
  }

  /// @brief Verifies sequential max returns value and index for PET sc grid function by checking tolerance-based numerical results, exact expected values.
  TEST(PETSc_GridFunction, SequentialMaxReturnsValueAndIndex)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {4, 4});
    P1 fes(mesh);
    Rodin::PETSc::Variational::GridFunction gf(fes);

    for (Index i = 0; i < gf.getSize(); ++i)
      gf[i] = static_cast<PetscScalar>(-10.0 - static_cast<Real>(i));

    constexpr Index expectedIdx = 5;
    gf[expectedIdx] = static_cast<PetscScalar>(4.25);
    gf.flush();

    Index idx = gf.getSize();
    const auto value = gf.max(idx);
    EXPECT_EQ(idx, expectedIdx);
    EXPECT_DOUBLE_EQ(static_cast<double>(PetscRealPart(value)), 4.25);
  }
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
