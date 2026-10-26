/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <array>
#include <memory>
#include <boost/bimap.hpp>

#include "Rodin/Variational.h"
#include "Rodin/Assembly/Sequential.h"
#ifdef RODIN_USE_OPENMP
#include "Rodin/Assembly/OpenMP.h"
#endif

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  LocalMesh interval()
  {
    auto mesh = LocalMesh::Builder()
                  .initialize(1)
                  .nodes(4)
                  .vertex({0})
                  .vertex({1})
                  .vertex({2})
                  .vertex({3})
                  .polytope(Polytope::Type::Segment, {0, 1})
                  .polytope(Polytope::Type::Segment, {1, 2})
                  .polytope(Polytope::Type::Segment, {2, 3})
                  .finalize();
    mesh.getConnectivity().compute(0, 1);
    mesh.setAttribute({1, 0}, {});
    mesh.setAttribute({1, 1}, 17);
    mesh.setAttribute({1, 2}, 23);
    for (Index vertex = 0; vertex < 4; ++vertex)
      mesh.setAttribute({0, vertex}, {});
    mesh.setAttribute({0, 0}, 31);
    mesh.setAttribute({0, 1}, 41);
    return mesh;
  }
}

/// @brief Dense standalone forms filter both topology and optional attributes on repeated assembly.
TEST(SequentialRegionSelection, DenseStandaloneForm)
{
  auto mesh = interval();
  P1 space(mesh);
  TrialFunction u(space);
  TestFunction v(space);
  using Space = decltype(space);
  using Form =
    BilinearForm<typename decltype(u)::SolutionType, Space, Space, Math::Matrix<Real>>;
  Form form(u, v);
  Assembly::BilinearFormAssemblyInput input(
    space, space, form.getLocalIntegrators(), form.getGlobalIntegrators());
  Assembly::Sequential<Math::Matrix<Real>, Form> assembler;
  Math::Matrix<Real> result;
  for (size_t selection = 0; selection < 4; ++selection)
  {
    SCOPED_TRACE(selection);
    Math::Matrix<Real> expected =
      Math::Matrix<Real>::Zero(space.getSize(), space.getSize());
    if (selection == 0)
    {
      form = Integral(u, v).over(17);
      const auto dofs = space.getDOFs(1, 1);
      for (Index i = 0; i < 2; ++i)
      {
        for (Index j = 0; j < 2; ++j)
          expected(dofs[i], dofs[j]) = i == j ? Real(1) / 3 : Real(1) / 6;
      }
    }
    else if (selection == 1)
    {
      form = BoundaryIntegral(u, v).over(31);
      const Index dof = space.getDOFs(0, 0)[0];
      expected(dof, dof) = 1;
    }
    else if (selection == 2)
    {
      form = InterfaceIntegral(u, v).over(41);
      const Index dof = space.getDOFs(0, 1)[0];
      expected(dof, dof) = 1;
    }
    else
      form = Integral(u, v).over(999);
    assembler.execute(result, input);
    EXPECT_LT((result - expected).norm(), 1e-12);
  }
}

template <class Config>
class MixedRegionSelection : public ::testing::Test
{};
using SelectionOperators =
  ::testing::Types<std::pair<Math::Matrix<Real>, std::false_type>,
    std::pair<Math::SparseMatrix<Real>, std::false_type>
#ifdef RODIN_USE_OPENMP
    ,
    std::pair<Math::Matrix<Real>, std::true_type>,
    std::pair<Math::SparseMatrix<Real>, std::true_type>
#endif
    >;
TYPED_TEST_SUITE(MixedRegionSelection, SelectionOperators);

/// @brief Mixed forms preserve exact block offsets, selected loads and the unrequested assembly side.
TYPED_TEST(MixedRegionSelection, FiltersAndTargetedReplacement)
{
  auto mesh = interval();
  P1 space(mesh);
  TrialFunction u(space), p(space);
  TestFunction v(space), q(space);
  using Operator = typename TypeParam::first_type;
  using System = Math::LinearSystem<Operator, Math::Vector<Real>>;
  using ProblemType = Problem<System, decltype(u), decltype(v), decltype(p), decltype(q)>;
  using Space = decltype(space);
  using Form = BilinearForm<typename decltype(p)::SolutionType, Space, Space, Operator>;
  Form preassembled(p, q);
  const size_t n = space.getSize();
  const Index first = space.getDOFs(0, 0)[0];
  const Index interior = space.getDOFs(0, 1)[0];
  const Index last = space.getDOFs(0, 3)[0];
  ProblemType problem(u, v, p, q);
#ifdef RODIN_USE_OPENMP
  using Parallel = typename TypeParam::second_type;
  using Assembler = std::conditional_t<Parallel::value,
    Assembly::OpenMP<System, ProblemType>, Assembly::Sequential<System, ProblemType>>;
#else
  using Assembler = Assembly::Sequential<System, ProblemType>;
#endif
  Assembler assembler;
  std::unique_ptr<decltype(assembler)> clone(assembler.copy());
  ASSERT_NE(clone, nullptr);
  auto trials = Tuple{std::ref(u), std::ref(p)};
  auto tests = Tuple{std::ref(v), std::ref(q)};
  std::array<size_t, 2> offsets{0, n};
  boost::bimap<FormLanguage::Base::UUID, size_t> trialMap, testMap;
  trialMap.right.insert({0, u.getUUID()});
  trialMap.right.insert({1, p.getUUID()});
  testMap.right.insert({0, v.getUUID()});
  testMap.right.insert({1, q.getUUID()});
  System result;
  for (size_t selection = 0; selection < 3; ++selection)
  {
    SCOPED_TRACE(selection);
    preassembled.getOperator().resize(n, n);
    preassembled.getOperator().setZero();
    Math::Matrix<Real> expected = Math::Matrix<Real>::Zero(2 * n, 2 * n);
    Math::Vector<Real> rhs = Math::Vector<Real>::Zero(2 * n);
    if (selection != 2)
    {
      preassembled.getOperator().coeffRef(last, last) = 5;
      preassembled.getOperator().coeffRef(first, last) = -2;
      expected(n + last, n + last) = 5;
      expected(n + first, n + last) = -2;
    }
    if (selection == 1)
    {
      problem = BoundaryIntegral(u, v).over(31) + InterfaceIntegral(p, q).over(41) +
        preassembled - BoundaryIntegral(RealFunction(2), v).over(31) -
        InterfaceIntegral(RealFunction(3), q).over(41);
      expected(first, first) = 1;
      expected(n + interior, n + interior) = 1;
      rhs[first] = 2;
      rhs[n + interior] = 3;
    }
    else
    {
      const Attribute attribute = selection == 0 ? 17 : 999;
      problem = Integral(u, v).over(attribute) + preassembled -
        Integral(RealFunction(2), v).over(attribute) -
        Integral(RealFunction(3), q).over(selection == 0 ? 23 : 999);
      if (selection == 0)
      {
        const auto dofs = space.getDOFs(1, 1);
        for (Index i = 0; i < 2; ++i)
        {
          rhs[dofs[i]] = 1;
          for (Index j = 0; j < 2; ++j)
            expected(dofs[i], dofs[j]) = i == j ? Real(1) / 3 : Real(1) / 6;
        }
        const auto loadDofs = space.getDOFs(1, 2);
        for (Index i = 0; i < 2; ++i)
          rhs[n + loadDofs[i]] = Real(1.5);
      }
    }
    auto body = problem.getBody();
    Assembly::ProblemAssemblyInput<std::decay_t<decltype(body)>, decltype(u), decltype(v),
      decltype(p), decltype(q)>
      input(body, trials, tests, offsets, offsets, trialMap, testMap, 2 * n, 2 * n);
    assembler.execute(result, input);
    EXPECT_LT((Math::Matrix<Real>(result.getOperator()) - expected).norm(), 1e-12);
    EXPECT_LT((result.getVector() - rhs).norm(), 1e-12);
    System copied;
    clone->execute(copied, input);
    EXPECT_LT((Math::Matrix<Real>(copied.getOperator()) - expected).norm(), 1e-12);
    EXPECT_LT((copied.getVector() - rhs).norm(), 1e-12);
    result.getVector().setConstant(91);
    assembler.execute(result, input, AssemblyTarget::LHS);
    EXPECT_LT((Math::Matrix<Real>(result.getOperator()) - expected).norm(), 1e-12);
    EXPECT_EQ(result.getVector(), Math::Vector<Real>::Constant(2 * n, 91));
    result.getOperator().coeffRef(0, 0) = 73;
    const Math::Matrix<Real> preserved = result.getOperator();
    assembler.execute(result, input, AssemblyTarget::RHS);
    EXPECT_EQ(Math::Matrix<Real>(result.getOperator()), preserved);
    EXPECT_LT((result.getVector() - rhs).norm(), 1e-12);
  }
}
