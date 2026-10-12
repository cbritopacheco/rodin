/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <type_traits>
#include <vector>
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

/// @brief Nonlocal assembly must retain test DOFs while visiting other trial cells.
TEST(Rodin_PETSc_NamedForm, NonlocalConstantKernelMatchesExactRankOneMatrix)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
  P1<PetscScalar> space(mesh);
  PETSc::Variational::TrialFunction u(space);
  PETSc::Variational::TestFunction v(space);
  auto kernel = [](const Point&, const Point&) { return PetscScalar(2); };
  Problem problem(u, v);
  problem = Integral(Potential(kernel, u), v);
  problem.assemble();
  const PetscInt indices[] = {0, 1, 2};
  const PetscReal moments[] = {0.5, 1, 0.5};
  PetscScalar entries[9];
  ASSERT_EQ(MatGetValues(
              problem.getLinearSystem().getOperator(), 3, indices, 3, indices, entries),
    PETSC_SUCCESS);
  for (size_t row = 0; row < 3; ++row)
  {
    for (size_t col = 0; col < 3; ++col)
    {
      EXPECT_LT(
        PetscAbsScalar(entries[3 * row + col] - 2 * moments[row] * moments[col]), 1e-12);
    }
  }
}

/// @brief Every PETSc nonlocal block matches exact, independently restricted cell moments.
TEST(Rodin_PETSc_NamedForm, NonlocalBlocksMatchExactRestrictedMoments)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
  mesh.setAttribute({1, 0}, 7);
  mesh.setAttribute({1, 1}, 8);
  P1<PetscScalar> space(mesh);
  PETSc::Variational::TrialFunction u(space), p(space);
  PETSc::Variational::TestFunction v(space), q(space);
  Problem problem(u, v, p, q);
  const auto k00 = [](const Point&, const Point&) { return PetscScalar(2); };
  const auto k01 = [](const Point&, const Point&) { return PetscScalar(-3); };
  const auto k10 = [](const Point&, const Point&) { return PetscScalar(0.5); };
  const auto k11 = [](const Point&, const Point&) { return PetscScalar(4); };
  const PetscScalar coefficients[2][2] = {{2, -3}, {0.5, 4}};
  const PetscInt indices[] = {0, 1, 2, 3, 4, 5};
  for (size_t restriction = 0; restriction < 5; ++restriction)
  {
    const bool restrictTrial = restriction == 1 || restriction >= 3;
    const bool restrictTest = restriction >= 2;
    typename decltype(problem)::ProblemBodyType body;
    body.getGlobalBFIs().add(Integral(Potential(k00, u), v));
    body.getGlobalBFIs().add(Integral(Potential(k01, p), v));
    body.getGlobalBFIs().add(Integral(Potential(k10, u), q));
    body.getGlobalBFIs().add(Integral(Potential(k11, p), q));
    for (auto& integral : body.getGlobalBFIs())
    {
      if (restrictTrial)
        integral.setTrialAttributes({restriction == 4 ? 999u : 7u});
      if (restrictTest)
        integral.setTestAttributes({restriction == 4 ? 999u : 8u});
    }
    problem = body;
    const PetscReal trialMoments[] = {
      0.5, restrictTrial ? 0.5 : 1, restrictTrial ? 0 : 0.5};
    const PetscReal testMoments[] = {restrictTest ? 0 : 0.5, restrictTest ? 0.5 : 1, 0.5};
    for (size_t repeat = 0; repeat < 2; ++repeat)
    {
      SCOPED_TRACE(restriction);
      SCOPED_TRACE(repeat);
      problem.assemble();
      PetscScalar entries[36];
      ASSERT_EQ(MatGetValues(problem.getLinearSystem().getOperator(), 6, indices, 6,
                  indices, entries),
        PETSC_SUCCESS);
      for (size_t row = 0; row < 6; ++row)
      {
        for (size_t col = 0; col < 6; ++col)
        {
          const PetscScalar expected = restriction == 4 ? PetscScalar(0)
                                                        : coefficients[row / 3][col / 3] *
              testMoments[row % 3] * trialMoments[col % 3];
          EXPECT_LT(PetscAbsScalar(entries[6 * row + col] - expected), 1e-12)
            << "row=" << row << ", col=" << col;
        }
      }
    }
  }
}

/// @brief Standalone PETSc nonlocal forms preserve independent trial/test cell moments.
TEST(Rodin_PETSc_NamedForm, StandaloneNonlocalFormsMatchRestrictedMoments)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
  mesh.setAttribute({1, 0}, 7);
  mesh.setAttribute({1, 1}, 8);
  P1<PetscScalar> testSpace(mesh);
  const auto check = [&](auto& trialSpace, const std::vector<PetscReal>& fullMoments,
                       const std::vector<PetscReal>& leftMoments) {
    PETSc::Variational::TrialFunction u(trialSpace);
    PETSc::Variational::TestFunction v(testSpace);
    const auto kernel = [](const Point&, const Point&) { return PetscScalar(2); };
    const PetscInt indices[] = {0, 1, 2};
    for (size_t restriction = 0; restriction < 5; ++restriction)
    {
      SCOPED_TRACE(restriction);
      const bool restrictTrial = restriction == 1 || restriction >= 3;
      const bool restrictTest = restriction >= 2;
      auto integral = Integral(Potential(kernel, u), v);
      if (restrictTrial)
      {
        integral.setTrialAttributes({restriction == 4 ? 999u : 7u});
      }
      if (restrictTest)
      {
        integral.setTestAttributes({restriction == 4 ? 999u : 8u});
      }
      BilinearForm form(u, v);
      form += integral;
      const auto& trialMoments = restrictTrial ? leftMoments : fullMoments;
      const PetscReal testMoments[] = {
        restrictTest ? 0 : 0.5, restrictTest ? 0.5 : 1, 0.5};
      const PetscInt columns = static_cast<PetscInt>(trialMoments.size());
      for (size_t repeat = 0; repeat < 2; ++repeat)
      {
        SCOPED_TRACE(repeat);
        form.assemble();
        PetscScalar entries[9];
        ASSERT_EQ(MatGetValues(form.getOperator(), 3, indices, columns, indices, entries),
          PETSC_SUCCESS);
        for (size_t row = 0; row < 3; ++row)
        {
          for (size_t col = 0; col < trialMoments.size(); ++col)
          {
            const PetscScalar expected = restriction == 4
              ? PetscScalar(0)
              : 2 * testMoments[row] * trialMoments[col];
            EXPECT_LT(PetscAbsScalar(entries[row * columns + col] - expected), 1e-12)
              << "row=" << row << ", col=" << col;
          }
        }
      }
    }
  };
  P1<PetscScalar> linearTrial(mesh);
  check(linearTrial, {0.5, 1, 0.5}, {0.5, 0.5, 0});
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
