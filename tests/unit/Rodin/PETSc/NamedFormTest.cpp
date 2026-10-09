/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <type_traits>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Solid/Linear/LinearElasticityForm.h>
#include <Rodin/Solid/Linear/LinearElasticityIntegral.h>
#include <Rodin/PETSc.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  void expectMatrixNear(::Mat actual, ::Mat expected)
  {
    ::Mat difference = nullptr;
    ASSERT_EQ(MatDuplicate(actual, MAT_COPY_VALUES, &difference), PETSC_SUCCESS);
    ASSERT_EQ(
      MatAXPY(difference, -1, expected, DIFFERENT_NONZERO_PATTERN), PETSC_SUCCESS);
    PetscReal error, norm;
    ASSERT_EQ(MatNorm(difference, NORM_FROBENIUS, &error), PETSC_SUCCESS);
    ASSERT_EQ(MatNorm(expected, NORM_FROBENIUS, &norm), PETSC_SUCCESS);
    EXPECT_LE(error, 1e-11 * std::max(PetscReal(1), norm));
    ASSERT_EQ(MatDestroy(&difference), PETSC_SUCCESS);
  }

  LocalMesh geometryMesh(Polytope::Type geometry)
  {
    LocalMesh mesh;
    if (geometry == Polytope::Type::Segment)
      mesh = LocalMesh::UniformGrid(geometry, {3});
    else if (geometry == Polytope::Type::Triangle ||
      geometry == Polytope::Type::Quadrilateral)
      mesh = LocalMesh::UniformGrid(geometry, {3, 3});
    else
      mesh = LocalMesh::UniformGrid(geometry, {2, 2, 2});
    for (size_t d = mesh.getDimension(); d > 0; --d)
    {
      mesh.getConnectivity().compute(d, d - 1);
      mesh.getConnectivity().compute(d - 1, d);
    }
    return mesh;
  }

  template <class Form>
  void expectBackends(Form& form, ::Mat expected)
  {
    static_assert(std::is_same_v<typename Form::OperatorType, ::Mat>);
    expectMatrixNear(form.getOperator(), expected);
    ::Mat sequential = nullptr;
    Assembly::Sequential<::Mat, Form> seq;
    seq.execute(sequential, form);
    expectMatrixNear(sequential, expected);
    seq.execute(sequential, form);
    expectMatrixNear(sequential, expected);
    ASSERT_EQ(MatDestroy(&sequential), PETSC_SUCCESS);
#ifdef RODIN_USE_OPENMP
    ::Mat parallel = nullptr;
    Assembly::OpenMP<::Mat, Form> omp;
    omp.setThreadCount(2).execute(parallel, form);
    expectMatrixNear(parallel, expected);
    omp.execute(parallel, form);
    expectMatrixNear(parallel, expected);
    ASSERT_EQ(MatDestroy(&parallel), PETSC_SUCCESS);
#endif
  }

  template <class FES>
  void checkScalarForms(FES& fes)
  {
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction v(fes);
    Real scale = 2;
    RealFunction coefficient([&scale](const Point& p) { return scale + p.x(); });
    MassForm mass(coefficient, u, v);
    DiffusionForm diffusion(coefficient, u, v);
    HelmholtzForm helmholtz(coefficient, Real(3), u, v);
    BilinearForm expected(u, v);
    for (const Real value : {2., 5.})
    {
      scale = value;
      mass.assemble();
      expected = Integral(coefficient * u, v);
      expected.assemble();
      expectBackends(mass, expected.getOperator());
      diffusion.assemble();
      expected = Integral(coefficient * Grad(u), Grad(v));
      expected.assemble();
      expectBackends(diffusion, expected.getOperator());
      helmholtz.assemble();
      expected = Integral(coefficient * Grad(u), Grad(v)) + Integral(Real(3) * u, v);
      expected.assemble();
      expectBackends(helmholtz, expected.getOperator());
    }
    MassForm unweighted(u, v);
    expected = Integral(u, v);
    expected.assemble();
    expectBackends(unweighted, expected.getOperator());
  }

  template <class FES>
  void checkVectorForms(FES& fes)
  {
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction v(fes);
    MassForm mass(Real(2), u, v);
    BilinearForm expected(u, v);
    expected = Integral(Real(2) * u, v);
    expected.assemble();
    expectBackends(mass, expected.getOperator());
    if (fes.getMesh().getDimension() > 1)
    {
      LinearElasticityForm elasticity(Real(2), Real(3), u, v);
      expected = LinearElasticityIntegral(u, v)(Real(2), Real(3));
      expected.assemble();
      expectBackends(elasticity, expected.getOperator());
    }
  }

  class PETScNamedFormGeometry : public ::testing::TestWithParam<Polytope::Type>
  {};

  TEST_P(PETScNamedFormGeometry, ScalarP1AndP2)
  {
    auto mesh = geometryMesh(GetParam());
    P1 linear(mesh);
    checkScalarForms(linear);
    H1 quadratic(std::integral_constant<size_t, 2>{}, mesh);
    checkScalarForms(quadratic);
  }

  TEST_P(PETScNamedFormGeometry, VectorP1AndP2)
  {
    auto mesh = geometryMesh(GetParam());
    P1 linear(mesh, mesh.getSpaceDimension());
    checkVectorForms(linear);
    H1 quadratic(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    checkVectorForms(quadratic);
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, PETScNamedFormGeometry,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge));

  TEST(PETSc_NamedForm, ReinterpolatedP2Coefficient)
  {
    auto mesh = geometryMesh(Polytope::Type::Triangle);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction v(fes);
    PETSc::Variational::GridFunction coefficient(fes);
    coefficient = RealFunction(1.);
    MassForm mass(coefficient, u, v);
    DiffusionForm diffusion(coefficient, u, v);
    HelmholtzForm helmholtz(coefficient, coefficient, u, v);
    BilinearForm expected(u, v);
    for (const Real scale : {2., 5., 0.})
    {
      coefficient = RealFunction([scale](const Point& p) { return scale * (1 + p.x()); });
      mass.assemble();
      expected = Integral(coefficient * u, v);
      expected.assemble();
      expectBackends(mass, expected.getOperator());
      diffusion.assemble();
      expected = Integral(coefficient * Grad(u), Grad(v));
      expected.assemble();
      expectBackends(diffusion, expected.getOperator());
      helmholtz.assemble();
      expected = Integral(coefficient * Grad(u), Grad(v)) + Integral(coefficient * u, v);
      expected.assemble();
      expectBackends(helmholtz, expected.getOperator());
    }
  }

  TEST(PETSc_NamedForm, CopyMoveAndPolymorphicLifetime)
  {
    auto mesh = geometryMesh(Polytope::Type::Triangle);
    P1 fes(mesh);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction v(fes);
    MassForm original(u, v);
    MassForm copied(original);
    EXPECT_NE(copied.getOperator(), original.getOperator());
    expectMatrixNear(copied.getOperator(), original.getOperator());
    const auto handle = copied.getOperator();
    auto moved = std::move(copied);
    EXPECT_EQ(moved.getOperator(), handle);
    EXPECT_EQ(copied.getOperator(), nullptr);
    copied.assemble();
    expectMatrixNear(copied.getOperator(), original.getOperator());
    std::unique_ptr<BilinearFormBase<::Mat>> clone(original.copy());
    EXPECT_NE(clone->getOperator(), original.getOperator());
    ASSERT_EQ(MatZeroEntries(original.getOperator()), PETSC_SUCCESS);
    expectMatrixNear(clone->getOperator(), moved.getOperator());
  }

  TEST(PETSc_NamedForm, AttributeChangesKeepPatternAndProblemComposition)
  {
    auto mesh = geometryMesh(Polytope::Type::Triangle);
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      mesh.setAttribute(
        {cell->getDimension(), cell->getIndex()}, cell->getIndex() % 2 + 1);
    }
    P1 fes(mesh);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction v(fes);
    MassForm mass(u, v);
    PetscObjectState pattern;
    ASSERT_EQ(MatGetNonzeroState(mass.getOperator(), &pattern), PETSC_SUCCESS);
    ASSERT_EQ(
      MatSetOption(mass.getOperator(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE),
      PETSC_SUCCESS);
    ::Mat initiallyRestricted = nullptr;
    Assembly::Sequential<::Mat, decltype(mass)> sequential;
    BilinearForm expected(u, v);
    for (const Attribute attribute : {1, 2, 1})
    {
      mass.over(attribute);
      expected = Integral(u, v).over(attribute);
      expected.assemble();
      expectMatrixNear(mass.getOperator(), expected.getOperator());
      sequential.execute(initiallyRestricted, mass);
      expectMatrixNear(initiallyRestricted, expected.getOperator());
      ASSERT_EQ(
        MatSetOption(initiallyRestricted, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE),
        PETSC_SUCCESS);
      PetscObjectState current;
      ASSERT_EQ(MatGetNonzeroState(mass.getOperator(), &current), PETSC_SUCCESS);
      EXPECT_EQ(current, pattern);
    }
    ASSERT_EQ(MatDestroy(&initiallyRestricted), PETSC_SUCCESS);
    Problem actual(u, v);
    actual = HelmholtzForm(Real(1), Real(2), u, v) - Integral(RealFunction(1.), v) +
      DirichletBC(u, RealFunction(0));
    actual.assemble();
    Problem reference(u, v);
    reference = Integral(Grad(u), Grad(v)) + Integral(Real(2) * u, v) -
      Integral(RealFunction(1.), v) + DirichletBC(u, RealFunction(0));
    reference.assemble();
    expectMatrixNear(
      actual.getLinearSystem().getOperator(), reference.getLinearSystem().getOperator());
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
