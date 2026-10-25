/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <array>
#include <memory>
#include <string>
#include <vector>

#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  template <class Scalar = Real, class Operator = Math::SparseMatrix<Scalar>,
    bool ExplicitSequential = false>
  void checkNonlocalBlocks(
    bool restrictTrial = false, bool restrictTest = false, bool missingAttribute = false)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    mesh.setAttribute({1, 0}, 7);
    mesh.setAttribute({1, 1}, 8);
    P1<Scalar> fes(mesh);
    TrialFunction u(fes);
    TrialFunction p(fes);
    TestFunction v(fes);
    TestFunction q(fes);
    using System = Math::LinearSystem<Operator, Math::Vector<Scalar>>;
    Problem<System, decltype(u), decltype(v), decltype(p), decltype(q)> problem(
      u, v, p, q);
    const Scalar factor = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(1, -2);
      else
        return Scalar(1);
    }();
    const Scalar coefficients[2][2] = {
      {Real(2) * factor, Real(-3) * factor}, {Real(0.5) * factor, Real(4) * factor}};
    auto kernel = [](Scalar value) {
      return [value](const Point&, const Point&) { return value; };
    };
    const auto k00 = kernel(coefficients[0][0]);
    const auto k01 = kernel(coefficients[0][1]);
    const auto k10 = kernel(coefficients[1][0]);
    const auto k11 = kernel(coefficients[1][1]);
    typename decltype(problem)::ProblemBodyType body;
    body.getGlobalBFIs().add(Integral(Potential(k00, u), v));
    body.getGlobalBFIs().add(Integral(Potential(k01, p), v));
    body.getGlobalBFIs().add(Integral(Potential(k10, u), q));
    body.getGlobalBFIs().add(Integral(Potential(k11, p), q));
    for (auto& integral : body.getGlobalBFIs())
    {
      if (restrictTrial)
        integral.setTrialAttributes({missingAttribute ? 999u : 7u});
      if (restrictTest)
        integral.setTestAttributes({missingAttribute ? 999u : 8u});
    }
    problem = body;

    // The two unit-length cells give exact nodal integrals (1/2, 1, 1/2).
    // A constant kernel separates the double integral into their outer product.
    const Real trialMoments[] = {
      0.5, restrictTrial ? 0.5 : 1.0, restrictTrial ? 0.0 : 0.5};
    const Real testMoments[] = {restrictTest ? 0.0 : 0.5, restrictTest ? 0.5 : 1.0, 0.5};
    auto check = [&](auto& instance, size_t stage) {
      instance.assemble();
      System sequentialSystem;
      if constexpr (ExplicitSequential)
      {
        using ProblemType = std::decay_t<decltype(problem)>;
        Tuple<std::reference_wrapper<decltype(u)>, std::reference_wrapper<decltype(p)>>
          trials(std::ref(u), std::ref(p));
        Tuple<std::reference_wrapper<decltype(v)>, std::reference_wrapper<decltype(q)>>
          tests(std::ref(v), std::ref(q));
        std::array<size_t, 2> offsets{0, fes.getSize()};
        boost::bimap<FormLanguage::Base::UUID, size_t> trialMap, testMap;
        trialMap.insert({u.getUUID(), 0});
        trialMap.insert({p.getUUID(), 1});
        testMap.insert({v.getUUID(), 0});
        testMap.insert({q.getUUID(), 1});
        auto sequentialBody = instance.getBody();
        typename ProblemType::AssemblyInput input(sequentialBody, trials, tests, offsets,
          offsets, trialMap, testMap, 2 * fes.getSize(), 2 * fes.getSize());
        Assembly::Sequential<System, ProblemType> assembler;
        assembler.execute(sequentialSystem, input);
      }
      const auto& matrix = [&]() -> const Operator& {
        if constexpr (ExplicitSequential)
          return sequentialSystem.getOperator();
        else
          return instance.getLinearSystem().getOperator();
      }();
      ASSERT_EQ(matrix.rows(), 6);
      ASSERT_EQ(matrix.cols(), 6);
      for (size_t row = 0; row < 6; ++row)
      {
        for (size_t col = 0; col < 6; ++col)
        {
          const Scalar expected = missingAttribute ? Scalar(0)
                                                   : coefficients[row / 3][col / 3] *
              testMoments[row % 3] * trialMoments[col % 3];
          EXPECT_NEAR(std::abs(matrix.coeff(row, col) - expected), 0.0, 1e-12)
            << "row=" << row << ", col=" << col << ", assembly=" << stage;
        }
      }
    };
    for (size_t stage = 0; stage < 2; ++stage)
      check(problem, stage);
    auto copied = problem;
    check(copied, 2);
    auto moved = std::move(copied);
    check(moved, 3);
    auto assigned = problem;
    assigned = problem;
    check(assigned, 4);
    auto moveSource = problem;
    auto moveAssigned = problem;
    moveAssigned = std::move(moveSource);
    check(moveAssigned, 5);
    std::unique_ptr<decltype(problem)> clone(problem.copy());
    check(*clone, 6);
  }

  /// @brief Checks every sparse block against exact double integrals.
  TEST(NonlocalBlock, SparseReal)
  {
    checkNonlocalBlocks();
  }

  /// @brief Trial-region restriction changes the column moments of every block.
  TEST(NonlocalBlock, RestrictedTrial)
  {
    checkNonlocalBlocks(true, false);
  }

  /// @brief Test-region restriction changes the row moments of every block.
  TEST(NonlocalBlock, RestrictedTest)
  {
    checkNonlocalBlocks(false, true);
  }

  /// @brief Independent region restrictions produce the selected cross-cell coupling.
  TEST(NonlocalBlock, RestrictedBoth)
  {
    checkNonlocalBlocks(true, true);
  }

  /// @brief Nonmatching attributes produce zero entries without retaining old assembly data.
  TEST(NonlocalBlock, MissingAttributes)
  {
    checkNonlocalBlocks(true, true, true);
  }

  /// @brief Dense real mixed assembly uses every nonlocal block and both independent regions.
  TEST(NonlocalBlock, DenseReal)
  {
    checkNonlocalBlocks<Real, Math::Matrix<Real>>();
    checkNonlocalBlocks<Real, Math::Matrix<Real>>(true, true);
    checkNonlocalBlocks<Real, Math::Matrix<Real>>(true, true, true);
  }

  /// @brief Dense complex mixed assembly preserves the kernel's nonreal coefficient.
  TEST(NonlocalBlock, DenseComplex)
  {
    checkNonlocalBlocks<Complex, Math::Matrix<Complex>>();
    checkNonlocalBlocks<Complex, Math::Matrix<Complex>>(true, true);
    checkNonlocalBlocks<Complex, Math::Matrix<Complex>>(true, true, true);
  }

  /// @brief Sparse complex mixed assembly agrees with exact separated double integrals.
  TEST(NonlocalBlock, SparseComplex)
  {
    checkNonlocalBlocks<Complex>();
    checkNonlocalBlocks<Complex>(true, true);
    checkNonlocalBlocks<Complex>(true, true, true);
  }

  /// @brief Explicit Sequential mixed assembly checks each independently restricted block.
  TEST(NonlocalBlock, SequentialSparseReal)
  {
    checkNonlocalBlocks<Real, Math::SparseMatrix<Real>, true>();
    checkNonlocalBlocks<Real, Math::SparseMatrix<Real>, true>(true, false);
    checkNonlocalBlocks<Real, Math::SparseMatrix<Real>, true>(false, true);
    checkNonlocalBlocks<Real, Math::SparseMatrix<Real>, true>(true, true);
    checkNonlocalBlocks<Real, Math::SparseMatrix<Real>, true>(true, true, true);
  }

  /// @brief Dense real Sequential blocks obey the exact separated moment formula.
  TEST(NonlocalBlock, SequentialDenseReal)
  {
    checkNonlocalBlocks<Real, Math::Matrix<Real>, true>(true, true);
    checkNonlocalBlocks<Real, Math::Matrix<Real>, true>(true, true, true);
  }

  /// @brief Sparse complex Sequential blocks retain the nonreal kernel coefficient.
  TEST(NonlocalBlock, SequentialSparseComplex)
  {
    checkNonlocalBlocks<Complex, Math::SparseMatrix<Complex>, true>(true, true);
    checkNonlocalBlocks<Complex, Math::SparseMatrix<Complex>, true>(true, true, true);
  }

  /// @brief Dense complex Sequential blocks obey the same independent moment oracle.
  TEST(NonlocalBlock, SequentialDenseComplex)
  {
    checkNonlocalBlocks<Complex, Math::Matrix<Complex>, true>(true, true);
    checkNonlocalBlocks<Complex, Math::Matrix<Complex>, true>(true, true, true);
  }
}

namespace
{
  /// @brief Closed-form interval moments check independent region filters and targeted assembly.
  template <class Scalar, class Operator, template <class, class> class Assembler>
  void checkSingleFieldRestrictions()
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    mesh.setAttribute({1, 0}, 3);
    mesh.setAttribute({1, 1}, 5);
    P1<Scalar> space(mesh);
    TrialFunction u(space);
    TestFunction v(space);
    using System = Math::LinearSystem<Operator, Math::Vector<Scalar>>;
    using ProblemType = Problem<System, decltype(u), decltype(v)>;
    using Body = typename ProblemType::ProblemBodyType;
    using Base = Assembly::AssemblyBase<System, ProblemType>;
    const Scalar factor = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(2, -3);
      else
        return Scalar(2);
    }();
    const auto one = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return ComplexFunction(Scalar(1));
      else
        return RealFunction(Scalar(1));
    }();
    GridFunction trialField(space), testField(space);
    trialField = [](const Point& p) { return Scalar(1 + p.x()); };
    testField = [](const Point& p) { return Scalar(3 + 2 * p.x()); };
    const auto kernel = [factor](const Point&, const Point&) { return factor; };
    bool missingSecondAttribute = false;
    auto selected = [&](Optional<Attribute> attribute, size_t cell) {
      return !attribute ||
        ((cell == 0 || !missingSecondAttribute) && *attribute == (cell == 0 ? 3 : 5));
    };
    auto moments = [&](Optional<Attribute> attribute) {
      Math::Vector<Scalar> result = Math::Vector<Scalar>::Zero(space.getSize());
      for (size_t cell = 0; cell < 2; ++cell)
      {
        if (!selected(attribute, cell))
          continue;
        const auto& dofs = space.getDOFs(1, cell);
        for (const auto dof : dofs)
          result[dof] += Scalar(0.5);
      }
      return result;
    };
    auto mass = [&](Optional<Attribute> attribute) {
      Math::Matrix<Scalar> result =
        Math::Matrix<Scalar>::Zero(space.getSize(), space.getSize());
      for (size_t cell = 0; cell < 2; ++cell)
      {
        if (!selected(attribute, cell))
          continue;
        const auto& dofs = space.getDOFs(1, cell);
        for (size_t row = 0; row < 2; ++row)
        {
          for (size_t column = 0; column < 2; ++column)
          {
            result(dofs[row], dofs[column]) +=
              Scalar(row == column ? Real(1) / 3 : Real(1) / 6);
          }
        }
      }
      return result;
    };
    auto exactMoment = [&](Optional<Attribute> attribute, bool trial) -> Real {
      if (!attribute)
        return trial ? 4 : 10;
      if (*attribute == 3)
        return trial ? 1.5 : 4;
      if (*attribute == 5 && !missingSecondAttribute)
        return trial ? 2.5 : 6;
      return 0;
    };
    std::vector<std::unique_ptr<Base>> assemblers;
    {
      Assembler<System, ProblemType> original;
      auto copied(original);
      auto moved(std::move(copied));
      assemblers.emplace_back(moved.copy());
    }
    assemblers.emplace_back(assemblers.front()->copy());
    Assembler<System, ProblemType> targeted;
    System system;
    const std::array<Optional<Attribute>, 4> attributes = {std::nullopt, 3, 5, 999};
    for (const bool missing : {false, true})
    {
      SCOPED_TRACE(missing);
      missingSecondAttribute = missing;
      mesh.setAttribute({1, 1}, missing ? Optional<Attribute>{} : Optional<Attribute>{5});
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
          auto local = Integral(u, v);
          if (trialAttribute)
            local.over(*trialAttribute);
          auto load = Integral(one, v);
          if (testAttribute)
            load.over(*testAttribute);
          Body body(integral);
          body.getLocalBFIs().add(local);
          body.getLFIs().add(UnaryMinus(load));
          Assembly::ProblemAssemblyInput<Body, decltype(u), decltype(v)> input(
            body, u, v);
          const auto trialMoments = moments(trialAttribute);
          const auto testMoments = moments(testAttribute);
          const Math::Matrix<Scalar> expected =
            factor * testMoments * trialMoments.transpose() + mass(trialAttribute);
          auto verify = [&] {
            const Math::Matrix<Scalar> actual = system.getOperator();
            EXPECT_LT((actual - expected).norm(), 1e-11);
            EXPECT_LT((system.getVector() - testMoments).norm(), 1e-11);
            Real massAction = 0;
            for (size_t cell = 0; cell < 2; ++cell)
            {
              if (selected(trialAttribute, cell))
                massAction += cell == 0 ? Real(37) / 6 : Real(91) / 6;
            }
            const Scalar exact = factor *
                Scalar(
                  exactMoment(trialAttribute, true) * exactMoment(testAttribute, false)) +
              Scalar(massAction);
            EXPECT_LT(
              std::abs(
                Math::dot(actual * trialField.getData(), testField.getData()) - exact),
              1e-10);
            EXPECT_LT(std::abs(Math::dot(system.getVector(), testField.getData()) -
                        Scalar(exactMoment(testAttribute, false))),
              1e-10);
          };
          for (auto& assembler : assemblers)
          {
            assembler->execute(system, input);
            verify();
            assembler->execute(system, input);
            verify();
          }
          system.getVector().setConstant(Scalar(17));
          const Math::Vector<Scalar> savedVector = system.getVector();
          targeted.execute(system, input, AssemblyTarget::LHS);
          EXPECT_LT(
            (Math::Matrix<Scalar>(system.getOperator()) - expected).norm(), 1e-11);
          EXPECT_EQ((system.getVector() - savedVector).norm(), 0);
          system.getOperator() *= Scalar(7);
          const Operator savedOperator = system.getOperator();
          targeted.execute(system, input, AssemblyTarget::RHS);
          EXPECT_EQ((system.getOperator() - savedOperator).norm(), 0);
          EXPECT_LT((system.getVector() - testMoments).norm(), 1e-11);
        }
      }
    }
  }
}

TEST(NonlocalSingleField, SequentialSparseReal)
{
  checkSingleFieldRestrictions<Real, Math::SparseMatrix<Real>, Assembly::Sequential>();
}
TEST(NonlocalSingleField, SequentialDenseReal)
{
  checkSingleFieldRestrictions<Real, Math::Matrix<Real>, Assembly::Sequential>();
}
TEST(NonlocalSingleField, SequentialSparseComplex)
{
  checkSingleFieldRestrictions<Complex, Math::SparseMatrix<Complex>,
    Assembly::Sequential>();
}
TEST(NonlocalSingleField, SequentialDenseComplex)
{
  checkSingleFieldRestrictions<Complex, Math::Matrix<Complex>, Assembly::Sequential>();
}
#ifdef RODIN_USE_OPENMP
TEST(NonlocalSingleField, OpenMPSparseReal)
{
  checkSingleFieldRestrictions<Real, Math::SparseMatrix<Real>, Assembly::OpenMP>();
}
TEST(NonlocalSingleField, OpenMPDenseReal)
{
  checkSingleFieldRestrictions<Real, Math::Matrix<Real>, Assembly::OpenMP>();
}
TEST(NonlocalSingleField, OpenMPSparseComplex)
{
  checkSingleFieldRestrictions<Complex, Math::SparseMatrix<Complex>, Assembly::OpenMP>();
}
TEST(NonlocalSingleField, OpenMPDenseComplex)
{
  checkSingleFieldRestrictions<Complex, Math::Matrix<Complex>, Assembly::OpenMP>();
}
#endif
