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
