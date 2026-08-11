/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <cassert>
#include <cmath>

#include <petsc.h>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/PETSc/Variational/GridFunction.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  /// @brief Returns the PETSc object id of a vector.
  ///
  /// Unlike the raw handle, an id is never reused by a later object, which
  /// makes it a reliable witness of whether a vector was reallocated.
  PetscObjectId objectId(const ::Vec& vec)
  {
    PetscObjectId id = 0;
    const PetscErrorCode ierr = PetscObjectGetId((PetscObject)vec, &id);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
    return id;
  }

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

  /// @brief Verifies that sync releases pending write access before raw PETSc
  ///        calls and releases pending read access afterwards.
  TEST(PETSc_GridFunction, SequentialSyncReleasesReadAndWriteAccess)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {3, 3});
    P1 fes(mesh);
    Rodin::PETSc::Variational::GridFunction gf(fes);

    gf[0] = static_cast<PetscScalar>(2.0);
    ASSERT_TRUE(gf.getArrayWrite().acquired);
    gf.sync();
    EXPECT_FALSE(gf.getArrayWrite().acquired);

    PetscReal norm = 0.0;
    PetscErrorCode ierr = VecNorm(gf.getData(), NORM_2, &norm);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
    EXPECT_DOUBLE_EQ(norm, 2.0);

    const auto& cgf = gf;
    EXPECT_DOUBLE_EQ(static_cast<double>(PetscRealPart(cgf[0])), 2.0);
    ASSERT_TRUE(cgf.getArrayRead().acquired);
    gf.sync();
    EXPECT_FALSE(cgf.getArrayRead().acquired);

    ierr = VecScale(gf.getData(), static_cast<PetscScalar>(3.0));
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;

    EXPECT_DOUBLE_EQ(gf.norm(), 6.0);
  }

  /// @brief Verifies that reductions synchronize DOFs written through
  ///        operator[] without requiring an explicit flush.
  TEST(PETSc_GridFunction, SequentialReductionsSynchronizePendingWrites)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {4, 4});
    P1 fes(mesh);
    Rodin::PETSc::Variational::GridFunction gf(fes);

    gf = static_cast<PetscScalar>(0.0);
    gf[0] = static_cast<PetscScalar>(3.0);
    gf[1] = static_cast<PetscScalar>(-4.0);
    gf[2] = static_cast<PetscScalar>(9.0);

    Index idx = gf.getSize();
    EXPECT_DOUBLE_EQ(gf.norm(), std::sqrt(106.0));
    EXPECT_DOUBLE_EQ(gf.norm(NORM_1), 16.0);
    EXPECT_DOUBLE_EQ(static_cast<double>(PetscRealPart(gf.min(idx))), -4.0);
    EXPECT_EQ(idx, 1);
    EXPECT_DOUBLE_EQ(static_cast<double>(PetscRealPart(gf.max(idx))), 9.0);
    EXPECT_EQ(idx, 2);
  }

  /// @brief Verifies that axpy accumulates a scaled grid function without
  ///        disturbing the operand, by checking exact expected values.
  TEST(PETSc_GridFunction, SequentialAxpyAccumulatesScaledGridFunction)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {4, 4});
    P1 fes(mesh);
    Rodin::PETSc::Variational::GridFunction y(fes);
    Rodin::PETSc::Variational::GridFunction x(fes);

    for (Index i = 0; i < y.getSize(); ++i)
    {
      y[i] = static_cast<PetscScalar>(1.0 + static_cast<Real>(i));
      x[i] = static_cast<PetscScalar>(2.0 * static_cast<Real>(i));
    }
    y.flush();
    x.flush();

    y.axpy(static_cast<PetscScalar>(0.5), x);

    const auto& cy = y;
    const auto& cx = x;
    for (Index i = 0; i < y.getSize(); ++i)
    {
      EXPECT_DOUBLE_EQ(static_cast<double>(PetscRealPart(cy[i])),
        1.0 + static_cast<double>(i) + 0.5 * 2.0 * static_cast<double>(i));
      // The operand is left untouched.
      EXPECT_DOUBLE_EQ(
        static_cast<double>(PetscRealPart(cx[i])), 2.0 * static_cast<double>(i));
    }
    cy.flush();
    cx.flush();
  }

  /// @brief Verifies that axpy flushes DOFs written through operator[] on the
  ///        accumulator before calling PETSc, by checking exact expected values.
  TEST(PETSc_GridFunction, SequentialAxpyFlushesPendingWrites)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {3, 3});
    P1 fes(mesh);
    Rodin::PETSc::Variational::GridFunction y(fes);
    Rodin::PETSc::Variational::GridFunction x(fes);

    for (Index i = 0; i < x.getSize(); ++i)
      x[i] = static_cast<PetscScalar>(4.0);
    x.flush();

    // Deliberately no explicit flush() on the accumulator: axpy must restore
    // its array before handing the vector to PETSc.
    for (Index i = 0; i < y.getSize(); ++i)
      y[i] = static_cast<PetscScalar>(3.0);

    y.axpy(static_cast<PetscScalar>(2.0), x);

    const auto& cy = y;
    for (Index i = 0; i < y.getSize(); ++i)
      EXPECT_DOUBLE_EQ(static_cast<double>(PetscRealPart(cy[i])), 11.0);
    cy.flush();
  }

  /// @brief Verifies that operator+= and operator-= remain the unit-coefficient
  ///        cases of axpy, by checking exact expected values.
  TEST(PETSc_GridFunction, SequentialPlusEqualsAndMinusEqualsAgreeWithAxpy)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Quadrilateral, {3, 3});
    P1 fes(mesh);
    Rodin::PETSc::Variational::GridFunction viaOperators(fes);
    Rodin::PETSc::Variational::GridFunction viaAxpy(fes);
    Rodin::PETSc::Variational::GridFunction x(fes);

    for (Index i = 0; i < x.getSize(); ++i)
      x[i] = static_cast<PetscScalar>(0.5 + static_cast<Real>(i));
    x.flush();

    viaOperators = static_cast<PetscScalar>(10.0);
    viaAxpy = static_cast<PetscScalar>(10.0);

    viaOperators += x;
    viaOperators -= x;
    viaOperators += x;

    viaAxpy.axpy(static_cast<PetscScalar>(1.0), x);
    viaAxpy.axpy(static_cast<PetscScalar>(-1.0), x);
    viaAxpy.axpy(static_cast<PetscScalar>(1.0), x);

    const auto& cOperators = viaOperators;
    const auto& cAxpy = viaAxpy;
    for (Index i = 0; i < x.getSize(); ++i)
    {
      const double expected = 10.0 + 0.5 + static_cast<double>(i);
      EXPECT_DOUBLE_EQ(static_cast<double>(PetscRealPart(cOperators[i])), expected);
      EXPECT_DOUBLE_EQ(static_cast<double>(PetscRealPart(cAxpy[i])), expected);
    }
    cOperators.flush();
    cAxpy.flush();
  }

  /// @brief Verifies the 2-, 1- and infinity-norms of the coefficient vector by
  ///        checking exact expected values.
  TEST(PETSc_GridFunction, SequentialNormMatchesHandComputedValues)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {4, 4});
    P1 fes(mesh);
    Rodin::PETSc::Variational::GridFunction gf(fes);

    gf = static_cast<PetscScalar>(0.0);
    gf[0] = static_cast<PetscScalar>(3.0);
    gf[1] = static_cast<PetscScalar>(-4.0);
    gf.flush();

    EXPECT_DOUBLE_EQ(gf.norm(), 5.0);
    EXPECT_DOUBLE_EQ(gf.norm(NORM_2), 5.0);
    EXPECT_DOUBLE_EQ(gf.norm(NORM_1), 7.0);
    EXPECT_DOUBLE_EQ(gf.norm(NORM_INFINITY), 4.0);
  }

  /// @brief Verifies that the norm of a zero grid function vanishes.
  TEST(PETSc_GridFunction, SequentialNormOfZeroVanishes)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {2, 2});
    P1 fes(mesh);
    Rodin::PETSc::Variational::GridFunction gf(fes);
    gf = static_cast<PetscScalar>(0.0);
    EXPECT_DOUBLE_EQ(gf.norm(), 0.0);
  }

  /// @brief Verifies that copy assignment reuses the destination PETSc vector
  ///        instead of reallocating it, keeping the Vec handle stable.
  TEST(PETSc_GridFunction, SequentialCopyAssignmentReusesVectorHandle)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {3, 3});
    P1 fes(mesh);
    Rodin::PETSc::Variational::GridFunction source(fes);
    Rodin::PETSc::Variational::GridFunction destination(fes);

    source = static_cast<PetscScalar>(2.5);
    destination = static_cast<PetscScalar>(-1.0);

    // A cached handle -- as a KSP, a LinearForm, or user code holds -- must
    // survive the assignment.  The object id is compared as well as the
    // pointer: PETSc readily hands a freshly created Vec the address of the
    // one just destroyed, but never its id.
    const ::Vec before = destination.getData();
    const PetscObjectId beforeId = objectId(before);
    destination = source;
    EXPECT_EQ(before, destination.getData());
    EXPECT_EQ(beforeId, objectId(destination.getData()));

    const auto& cDestination = destination;
    for (Index i = 0; i < destination.getSize(); ++i)
      EXPECT_DOUBLE_EQ(static_cast<double>(PetscRealPart(cDestination[i])), 2.5);
    cDestination.flush();
  }

  /// @brief Verifies that copy assignment is a deep copy: the operands do not
  ///        alias after the assignment.
  TEST(PETSc_GridFunction, SequentialCopyAssignmentIsDeep)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {3, 3});
    P1 fes(mesh);
    Rodin::PETSc::Variational::GridFunction source(fes);
    Rodin::PETSc::Variational::GridFunction destination(fes);

    source = static_cast<PetscScalar>(1.0);
    destination = source;

    EXPECT_NE(source.getData(), destination.getData());

    source = static_cast<PetscScalar>(9.0);

    const auto& cDestination = destination;
    for (Index i = 0; i < destination.getSize(); ++i)
      EXPECT_DOUBLE_EQ(static_cast<double>(PetscRealPart(cDestination[i])), 1.0);
    cDestination.flush();
  }

  /// @brief Verifies that copy assignment flushes DOFs written through
  ///        operator[] on the source before copying.
  TEST(PETSc_GridFunction, SequentialCopyAssignmentFlushesPendingWrites)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {3, 3});
    P1 fes(mesh);
    Rodin::PETSc::Variational::GridFunction source(fes);
    Rodin::PETSc::Variational::GridFunction destination(fes);

    for (Index i = 0; i < source.getSize(); ++i)
      source[i] = static_cast<PetscScalar>(6.0);
    source.flush();

    // Deliberately no explicit flush() on the destination: the assignment
    // must restore its array before overwriting the vector.
    for (Index i = 0; i < destination.getSize(); ++i)
      destination[i] = static_cast<PetscScalar>(-6.0);

    destination = source;

    const auto& cDestination = destination;
    for (Index i = 0; i < destination.getSize(); ++i)
      EXPECT_DOUBLE_EQ(static_cast<double>(PetscRealPart(cDestination[i])), 6.0);
    cDestination.flush();
  }

  /// @brief Verifies that copy assignment into a moved-from grid function
  ///        allocates a vector rather than copying into a null handle.
  TEST(PETSc_GridFunction, SequentialCopyAssignmentIntoMovedFromAllocates)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {3, 3});
    P1 fes(mesh);
    Rodin::PETSc::Variational::GridFunction source(fes);
    source = static_cast<PetscScalar>(4.5);

    Rodin::PETSc::Variational::GridFunction movedFrom(fes);
    Rodin::PETSc::Variational::GridFunction sink(std::move(movedFrom));
    ASSERT_EQ(movedFrom.getData(), PETSC_NULLPTR);

    movedFrom = source;

    ASSERT_NE(movedFrom.getData(), PETSC_NULLPTR);
    EXPECT_NE(movedFrom.getData(), source.getData());

    const auto& cMovedFrom = movedFrom;
    for (Index i = 0; i < movedFrom.getSize(); ++i)
      EXPECT_DOUBLE_EQ(static_cast<double>(PetscRealPart(cMovedFrom[i])), 4.5);
    cMovedFrom.flush();
  }

  /// @brief Verifies that self-assignment leaves the grid function unchanged.
  TEST(PETSc_GridFunction, SequentialSelfAssignmentIsANoOp)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {3, 3});
    P1 fes(mesh);
    Rodin::PETSc::Variational::GridFunction gf(fes);
    gf = static_cast<PetscScalar>(8.0);

    const ::Vec before = gf.getData();
    const PetscObjectId beforeId = objectId(before);
    const auto& self = gf;
    gf = self;

    EXPECT_EQ(before, gf.getData());
    EXPECT_EQ(beforeId, objectId(gf.getData()));
    const auto& cgf = gf;
    for (Index i = 0; i < gf.getSize(); ++i)
      EXPECT_DOUBLE_EQ(static_cast<double>(PetscRealPart(cgf[i])), 8.0);
    cgf.flush();
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
