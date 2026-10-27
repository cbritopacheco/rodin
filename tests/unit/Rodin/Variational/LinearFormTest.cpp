/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <sstream>
#include <string>
#include <type_traits>
#include "Rodin/Test/Random.h"

#include "Rodin/Variational.h"
#include "Rodin/Assembly/Default.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Test::Random;

namespace Rodin::Tests::Unit
{
  /// @brief Verifies sanity test build for variational real P1 linear form by checking exact expected values.
  TEST(Rodin_Variational_Real_P1_LinearForm, SanityTest_Build)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    P1 fes(mesh);
    TestFunction v(fes);
    LinearForm lf(v);

    EXPECT_EQ(&lf.getTestFunction(), &v);
  }

  /// @brief Verifies fuzzy test uniform grid 4 x 4 for variational real P1 linear form by checking form assembly.
  TEST(Rodin_Variational_Real_P1_LinearForm, FuzzyTest_UniformGrid_4x4)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    P1 fes(mesh);
    TestFunction v(fes);
    LinearForm lf(v);
    lf = Integral(v);
    lf.assemble();

    const auto& vector = lf.getVector();
    EXPECT_GT(vector.size(), 0);
  }

  /// @brief Verifies copy constructor for variational real P1 linear form by checking exact expected values, copy semantics.
  TEST(Rodin_Variational_Real_P1_LinearForm, CopyConstructor)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    P1 fes(mesh);
    TestFunction v(fes);
    LinearForm lf(v);

    LinearForm lf_copy(lf);
    EXPECT_EQ(&lf_copy.getTestFunction().getUUID(), &lf.getTestFunction().getUUID());
  }

  /// @brief Verifies move constructor for variational real P1 linear form by checking exact expected values, move semantics.
  TEST(Rodin_Variational_Real_P1_LinearForm, MoveConstructor)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    P1 fes(mesh);
    TestFunction v(fes);
    LinearForm lf(v);
    LinearForm lf_moved(std::move(lf));
    EXPECT_EQ(&lf_moved.getTestFunction(), &v);
  }

  /// @brief Verifies assignment for variational real P1 linear form by checking false predicates.
  TEST(Rodin_Variational_Real_P1_LinearForm, Assignment)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    P1 fes(mesh);
    TestFunction v(fes);
    LinearForm lf(v);
    lf = Integral(v);
    EXPECT_FALSE(lf.getIntegrators().empty());
  }

  /// @brief Verifies addition assignment for variational real P1 linear form by checking false predicates.
  TEST(Rodin_Variational_Real_P1_LinearForm, AdditionAssignment)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    P1 fes(mesh);
    TestFunction v(fes);
    LinearForm lf(v);
    lf += Integral(v);
    EXPECT_FALSE(lf.getIntegrators().empty());
  }

  /// @brief Verifies subtraction assignment for variational real P1 linear form by checking false predicates.
  TEST(Rodin_Variational_Real_P1_LinearForm, SubtractionAssignment)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    P1 fes(mesh);
    TestFunction v(fes);
    LinearForm lf(v);
    lf -= Integral(v);
    EXPECT_FALSE(lf.getIntegrators().empty());
  }

  /// @brief Verifies assemble and get vector for variational real P1 linear form by checking exact expected values, form assembly.
  TEST(Rodin_Variational_Real_P1_LinearForm, AssembleAndGetVector)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    P1 fes(mesh);
    TestFunction v(fes);
    LinearForm lf(v);
    lf = Integral(v);
    lf.assemble();
    const auto& vector = lf.getVector();
    auto& mutable_vector = lf.getVector();
    EXPECT_GT(vector.size(), 0);
    EXPECT_EQ(vector.size(), mutable_vector.size());
  }

  /// @brief Verifies copy for variational real P1 linear form by checking copy semantics.
  TEST(Rodin_Variational_Real_P1_LinearForm, Copy)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    P1 fes(mesh);
    TestFunction v(fes);
    LinearForm lf(v);
    lf = Integral(v);
    auto copied = lf.copy();
    EXPECT_NE(copied, nullptr);
    delete copied;
  }

  /// @brief Verifies sanity test build for variational vector P1 linear form by checking exact expected values.
  TEST(Rodin_Variational_Vector_P1_LinearForm, SanityTest_Build)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    P1 fes(mesh, vdim);
    TestFunction v(fes);
    LinearForm lf(v);
    EXPECT_EQ(&lf.getTestFunction(), &v);
  }

  /// @brief Verifies assemble vector function for variational vector P1 linear form by checking form assembly.
  TEST(Rodin_Variational_Vector_P1_LinearForm, AssembleVectorFunction)
  {
    constexpr size_t vdim = 2;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    P1 fes(mesh, vdim);
    TestFunction v(fes);
    LinearForm lf(v);
    lf = Integral(Dot(v, VectorFunction{1, 1}));
    lf.assemble();
    const auto& vector = lf.getVector();
    EXPECT_GT(vector.size(), 0);
  }

  /// @brief Verifies complex linear-form action conjugates the test coefficients.
  TEST(Rodin_Variational_Complex_P1_LinearForm, ActionConjugatesTest)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    P1<Complex> fes(mesh);
    TestFunction v(fes);
    LinearForm form(v);
    form.getVector().resize(fes.getSize());
    form.getVector().setZero();
    form.getVector()(0) = Complex(1, 2);
    form.getVector()(1) = Complex(-3, 1);

    GridFunction test(fes);
    test.getData().setZero();
    test.getData()(0) = Complex(2, -1);
    test.getData()(1) = Complex(1, 4);

    const Complex expected =
      Complex(1, 2) * Complex(2, 1) + Complex(-3, 1) * Complex(1, -4);
    EXPECT_NEAR(std::abs(form(test) - expected), 0.0, 1e-14);
  }

  /// @brief Verifies clear integrators for variational real P1 linear form by checking false predicates.
  TEST(Rodin_Variational_Real_P1_LinearForm, ClearIntegrators)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    P1 fes(mesh);
    TestFunction v(fes);
    LinearForm lf(v);
    lf += Integral(v);
    EXPECT_FALSE(lf.getIntegrators().empty());
    lf = Integral(v); // This should clear and add new integrator
    EXPECT_FALSE(lf.getIntegrators().empty());
  }
}

namespace Rodin::Tests::Unit
{
  /// @brief Rejected assignment, addition and subtraction preserve the assembled load.
  TEST(FormOwnership, LinearContributionsPreserveLoad)
  {
    auto checkScalar = []<class Scalar>() {
      auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
      mesh.scale(Real(0.5));
      auto checkSpace = [&](const auto& space) {
        TestFunction v(space), foreignV(space);
        LinearForm form(v);
        form = Integral(v);
        form.assemble();
        const Math::Vector<Scalar> original = form.getVector();
        GridFunction constant(space);
        constant = Scalar(1);
        EXPECT_NEAR(std::abs(constant.getData().dot(original) - Scalar(1)), 0, 1e-12);
        const auto foreign = Integral(foreignV);
        std::ostringstream uuid;
        uuid << foreignV.getUUID();
        auto reject = [&](auto&& operation) {
          try
          {
            operation();
            FAIL() << "A foreign test-function contribution must be rejected.";
          }
          catch (const Alert::Exception& error)
          {
            const std::string message = error.what();
            EXPECT_NE(message.find("TestFunction"), std::string::npos);
            EXPECT_NE(message.find(uuid.str()), std::string::npos);
            EXPECT_NE(
              message.find("does not appear in the declaration"), std::string::npos);
          }
          EXPECT_EQ(form.getIntegrators().size(), 1u);
          form.assemble();
          EXPECT_NEAR((form.getVector() - original).norm(), 0, 1e-12);
        };
        reject([&] { form = foreign; });
        reject([&] { form += foreign; });
        reject([&] { form -= foreign; });
      };
      const P1<Scalar> p1(mesh);
      const H1<2, Scalar> h1(std::integral_constant<size_t, 2>{}, mesh);
      checkSpace(p1);
      checkSpace(h1);
    };
    checkScalar.operator()<Real>();
    checkScalar.operator()<Complex>();
  }
}

/// @brief Copy/move assignment transfers the source integrators for subsequent assembly.
TEST(LinearFormAssignment, IntegratorsSurviveReassembly)
{
  const auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
  P1 fes(mesh);
  TestFunction v(fes);
  LinearForm source(v), destination(v), moved(v);
  source = Integral(Real(3) * v);
  destination = Integral(Real(7) * v);
  source.assemble();
  const auto expected = source.getVector();
  destination = source;
  destination.assemble();
  EXPECT_LE((destination.getVector() - expected).norm(), 1e-13);
  moved = std::move(destination);
  moved.assemble();
  EXPECT_LE((moved.getVector() - expected).norm(), 1e-13);
}
