/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Independent moment oracles for PETSc region-restricted nonlocal assembly. */
#include <gtest/gtest.h>
#include <petsc.h>
#include <array>
#include <string>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/PETSc.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  /// @brief Trial/test restrictions give the exact outer product of their moment vectors.
  template <template <class, class> class Assembler, class Scalar = PetscScalar>
  void checkRestrictedNonlocalAssembly()
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    mesh.setAttribute({1, 0}, 3);
    mesh.setAttribute({1, 1}, 5);
    const Scalar factor = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(2, -3);
      else
        return Scalar(2);
    }();
    auto check = [&](auto& space) {
      PETSc::Variational::TrialFunction u(space);
      PETSc::Variational::TestFunction v(space);
      using System = PETSc::Math::LinearSystem;
      using ProblemType = Problem<System, decltype(u), decltype(v)>;
      using Body = typename ProblemType::ProblemBodyType;
      const auto kernel = [factor](const Point&, const Point&) { return factor; };
      const auto one = [] {
        if constexpr (std::is_same_v<Scalar, Complex>)
          return ComplexFunction(Scalar(1));
        else
          return RealFunction(Scalar(1));
      }();
      const auto first = [] {
        if constexpr (std::is_same_v<Scalar, Complex>)
        {
          return ComplexFunction(
            [](const Point& point) { return Scalar(1 + point.x()); });
        }
        else
          return RealFunction([](const Point& point) { return Scalar(1 + point.x()); });
      }();
      const auto second = [] {
        if constexpr (std::is_same_v<Scalar, Complex>)
        {
          return ComplexFunction(
            [](const Point& point) { return Scalar(3 + 2 * point.x()); });
        }
        else
        {
          return RealFunction(
            [](const Point& point) { return Scalar(3 + 2 * point.x()); });
        }
      }();
      GridFunction trialField(space), testField(space);
      trialField = first;
      testField = second;
      TestFunction momentTest(space);
      auto moments = [&](Optional<Attribute> attribute) {
        LinearForm form(momentTest);
        auto integral = Integral(one, momentTest);
        if (attribute)
          integral.over(*attribute);
        form = integral;
        form.assemble();
        return Math::Vector<Scalar>(form.getVector());
      };
      auto exactMoment = [](Optional<Attribute> attribute, bool trial) -> Real {
        if (!attribute)
          return trial ? 4 : 10;
        if (attribute == 3)
          return trial ? 1.5 : 4;
        if (attribute == 5)
          return trial ? 2.5 : 6;
        return 0;
      };
      System system(PETSC_COMM_SELF);
      Assembler<System, ProblemType> assembler;
      const std::array<Optional<Attribute>, 4> attributes = {std::nullopt, 3, 5, 999};
      for (const auto trialAttribute : attributes)
      {
        for (const auto testAttribute : attributes)
        {
          SCOPED_TRACE(trialAttribute ? std::to_string(*trialAttribute) : "all");
          SCOPED_TRACE(testAttribute ? std::to_string(*testAttribute) : "all");
          auto integral = Integral(Potential(kernel, u), v);
          if (trialAttribute)
            integral.setTrialAttributes({*trialAttribute});
          if (testAttribute)
            integral.setTestAttributes({*testAttribute});
          Body body(integral);
          Assembly::ProblemAssemblyInput<Body, decltype(u), decltype(v)> input(
            body, u, v);
          const auto trialMoments = moments(trialAttribute);
          const auto testMoments = moments(testAttribute);
          const Math::Matrix<Scalar> expected =
            factor * testMoments * trialMoments.transpose();
          for (size_t repeat = 0; repeat < 2; ++repeat)
          {
            assembler.execute(system, input);
            Math::Matrix<Scalar> actual(space.getSize(), space.getSize());
            for (PetscInt row = 0; row < static_cast<PetscInt>(space.getSize()); ++row)
            {
              for (PetscInt column = 0; column < static_cast<PetscInt>(space.getSize());
                   ++column)
              {
                Scalar entry;
                ASSERT_EQ(MatGetValues(system.getOperator(), 1, &row, 1, &column, &entry),
                  PETSC_SUCCESS);
                actual(row, column) = entry;
              }
            }
            EXPECT_LT((actual - expected).norm(), 1e-11);
            const Scalar action =
              Math::dot(actual * trialField.getData(), testField.getData());
            const Scalar exact = factor *
              Scalar(
                exactMoment(trialAttribute, true) * exactMoment(testAttribute, false));
            EXPECT_LT(std::abs(action - exact), 1e-10);
          }
        }
      }
    };
    P1<Scalar> linear(mesh);
    check(linear);
  }
}

TEST(Rodin_PETSc_RestrictedNonlocal, Sequential)
{
  checkRestrictedNonlocalAssembly<PETSc::Assembly::Sequential>();
}
#ifdef RODIN_USE_OPENMP
TEST(Rodin_PETSc_RestrictedNonlocal, OpenMP)
{
  checkRestrictedNonlocalAssembly<PETSc::Assembly::OpenMP>();
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
