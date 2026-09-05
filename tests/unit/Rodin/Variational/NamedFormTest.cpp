/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <type_traits>

#include "Rodin/Variational.h"
#include "Rodin/Assembly/Default.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  namespace
  {
    LocalMesh unitSquare()
    {
      LocalMesh mesh =
        LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      return mesh;
    }

    LocalMesh geometryMesh(Polytope::Type geometry)
    {
      LocalMesh mesh;
      switch (geometry)
      {
        case Polytope::Type::Segment:
          mesh = LocalMesh::UniformGrid(geometry, {2});
          break;
        case Polytope::Type::Triangle:
        case Polytope::Type::Quadrilateral:
          mesh = LocalMesh::UniformGrid(geometry, {2, 2});
          break;
        case Polytope::Type::Tetrahedron:
        case Polytope::Type::Hexahedron:
        case Polytope::Type::Pyramid:
        case Polytope::Type::Wedge:
          mesh = LocalMesh::UniformGrid(geometry, {2, 2, 2});
          break;
        default:
          assert(false);
          return {};
      }
      for (size_t d = mesh.getDimension(); d > 0; --d)
        mesh.getConnectivity().compute(d, d - 1);
      return mesh;
    }

    template <class Actual, class Expected>
    void expectNear(const Actual& actual, const Expected& expected)
    {
      const Math::Matrix<Real> a = actual;
      const Math::Matrix<Real> e = expected;
      ASSERT_EQ(a.rows(), e.rows());
      ASSERT_EQ(a.cols(), e.cols());
      for (Eigen::Index i = 0; i < a.rows(); ++i)
        for (Eigen::Index j = 0; j < a.cols(); ++j)
          EXPECT_NEAR(a(i, j), e(i, j), 1e-12)
            << "entry (" << i << ", " << j << ")";
    }

    template <class FES>
    void expectMassFormMatches(FES& fes)
    {
      TrialFunction u(fes);
      TestFunction v(fes);

      MassForm mass(u, v);
      BilinearForm expectedMass(u, v);
      expectedMass = Integral(u, v);
      expectedMass.assemble();
      expectNear(mass.getOperator(), expectedMass.getOperator());
    }

    template <class FES>
    void expectDiffusionFormMatches(FES& fes)
    {
      TrialFunction u(fes);
      TestFunction v(fes);

      DiffusionForm diffusion(u, v);
      BilinearForm expectedDiffusion(u, v);
      expectedDiffusion = Integral(Grad(u), Grad(v));
      expectedDiffusion.assemble();
      expectNear(diffusion.getOperator(), expectedDiffusion.getOperator());
    }
  }

  class Rodin_Variational_NamedFormGeometry
    : public ::testing::TestWithParam<Polytope::Type>
  {};

  TEST_P(Rodin_Variational_NamedFormGeometry, P0MatchesIntegrals)
  {
    auto mesh = geometryMesh(GetParam());
    P0 fes(mesh);
    expectMassFormMatches(fes);
  }

  TEST_P(Rodin_Variational_NamedFormGeometry, P1MatchesIntegrals)
  {
    auto mesh = geometryMesh(GetParam());
    P1 fes(mesh);
    expectMassFormMatches(fes);
    expectDiffusionFormMatches(fes);
  }

  TEST_P(Rodin_Variational_NamedFormGeometry, H1P2MatchesIntegrals)
  {
    auto mesh = geometryMesh(GetParam());
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    expectMassFormMatches(fes);
    expectDiffusionFormMatches(fes);
  }

  INSTANTIATE_TEST_SUITE_P(
    AllCellGeometries,
    Rodin_Variational_NamedFormGeometry,
    ::testing::Values(
      Polytope::Type::Segment,
      Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron,
      Polytope::Type::Pyramid,
      Polytope::Type::Wedge));

  TEST(Rodin_Variational_NamedForm, SelectsContextAssembly)
  {
    auto mesh = unitSquare();
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    MassForm mass(u, v);

    using Form = decltype(mass);
    using Expected =
      typename Assembly::Default<Context::Local, Context::Local>
        ::template Type<typename Form::OperatorType, Form>;
    static_assert(std::is_same_v<typename Form::AssemblyType, Expected>);

    typename Form::OperatorType sequentialOperator;
    Assembly::Sequential<typename Form::OperatorType, Form> sequential;
    sequential.execute(sequentialOperator, mass);
    expectNear(sequentialOperator, mass.getOperator());

#ifdef RODIN_USE_OPENMP
    typename Form::OperatorType openMPOperator;
    Assembly::OpenMP<typename Form::OperatorType, Form> openMP;
    openMP.execute(openMPOperator, mass);
    expectNear(openMPOperator, mass.getOperator());
#endif
  }

  TEST(Rodin_Variational_NamedForm, MassMatchesIntegral)
  {
    auto mesh = unitSquare();
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    TrialFunction u(fes);
    TestFunction v(fes);

    MassForm actual(u, v);
    BilinearForm expected(u, v);
    expected = Integral(u, v);
    expected.assemble();

    expectNear(actual.getOperator(), expected.getOperator());
  }

  TEST(Rodin_Variational_NamedForm, DiffusionMatchesIntegral)
  {
    auto mesh = unitSquare();
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    TrialFunction u(fes);
    TestFunction v(fes);

    DiffusionForm actual(u, v);
    BilinearForm expected(u, v);
    expected = Integral(Grad(u), Grad(v));
    expected.assemble();

    expectNear(actual.getOperator(), expected.getOperator());
  }

  TEST(Rodin_Variational_NamedForm, P0gMatchesIntegrals)
  {
    auto mesh = unitSquare();
    P0g fes(mesh);
    expectMassFormMatches(fes);
  }

  TEST(Rodin_Variational_NamedForm, VectorP1MatchesIntegrals)
  {
    auto mesh = unitSquare();
    P1 fes(mesh, 2);
    expectMassFormMatches(fes);
  }

  TEST(Rodin_Variational_NamedForm, MixedH1OrdersMatchIntegrals)
  {
    auto mesh = unitSquare();
    H1 trialFES(std::integral_constant<size_t, 2>{}, mesh);
    H1 testFES(std::integral_constant<size_t, 1>{}, mesh);
    TrialFunction u(trialFES);
    TestFunction v(testFES);

    MassForm mass(u, v);
    BilinearForm expectedMass(u, v);
    expectedMass = Integral(u, v);
    expectedMass.assemble();
    expectNear(mass.getOperator(), expectedMass.getOperator());

    DiffusionForm diffusion(u, v);
    BilinearForm expectedDiffusion(u, v);
    expectedDiffusion = Integral(Grad(u), Grad(v));
    expectedDiffusion.assemble();
    expectNear(diffusion.getOperator(), expectedDiffusion.getOperator());
  }

  TEST(Rodin_Variational_NamedForm, ReassemblyMatchesIntegral)
  {
    auto mesh = unitSquare();
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    TrialFunction u(fes);
    TestFunction v(fes);

    MassForm mass(u, v);
    DiffusionForm diffusion(u, v);

    BilinearForm expectedMass(u, v);
    expectedMass = Integral(u, v);
    expectedMass.assemble();

    BilinearForm expectedDiffusion(u, v);
    expectedDiffusion = Integral(Grad(u), Grad(v));
    expectedDiffusion.assemble();

    // The first assembly builds the sparsity pattern; every later one
    // scatters into the value array it cached, so the reuse path is only
    // exercised from the second call onwards.
    const auto massNonZeros = mass.getOperator().nonZeros();
    const auto diffusionNonZeros = diffusion.getOperator().nonZeros();
    for (size_t pass = 0; pass < 3; ++pass)
    {
      mass.assemble();
      diffusion.assemble();
      expectNear(mass.getOperator(), expectedMass.getOperator());
      expectNear(diffusion.getOperator(), expectedDiffusion.getOperator());
      EXPECT_EQ(mass.getOperator().nonZeros(), massNonZeros);
      EXPECT_EQ(diffusion.getOperator().nonZeros(), diffusionNonZeros);
    }
  }

  TEST(Rodin_Variational_NamedForm, ReassemblyAgreesAcrossBackends)
  {
    auto mesh = unitSquare();
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    MassForm mass(u, v);

    using Form = decltype(mass);
    typename Form::OperatorType sequentialOperator;
    Assembly::Sequential<typename Form::OperatorType, Form> sequential;
    sequential.execute(sequentialOperator, mass);
    sequential.execute(sequentialOperator, mass);
    expectNear(sequentialOperator, mass.getOperator());

#ifdef RODIN_USE_OPENMP
    typename Form::OperatorType openMPOperator;
    Assembly::OpenMP<typename Form::OperatorType, Form> openMP;
    openMP.execute(openMPOperator, mass);
    openMP.execute(openMPOperator, mass);
    expectNear(openMPOperator, mass.getOperator());
#endif
  }

  TEST(Rodin_Variational_NamedForm, ComposesWithIntegralInProblem)
  {
    auto mesh = unitSquare();
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);

    Problem problem(u, v);
    problem = MassForm(u, v) + Integral(Grad(u), Grad(v));
    problem.assemble();

    BilinearForm expected(u, v);
    expected = Integral(u, v) + Integral(Grad(u), Grad(v));
    expected.assemble();

    expectNear(
      problem.getLinearSystem().getOperator(), expected.getOperator());
  }
}
