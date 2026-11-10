/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
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
  /// @brief Compare periodic assembly with an independently condensed unconstrained system.
  template <class Operator, template <class, class> class Backend>
  void checkPeriodicAssembly()
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {5});
    mesh.scale(0.25);
    mesh.getConnectivity().compute(0, 1);
    P1 space(mesh);
    TrialFunction u(space);
    TestFunction v(space);
    using System = Math::LinearSystem<Operator, Math::Vector<Real>>;
    using ProblemType = Problem<System, decltype(u), decltype(v)>;
    using Body = typename ProblemType::ProblemBodyType;
    Backend<System, ProblemType> assembler;
    const Index first = space.getDOFs(0, 0)[0];
    const Index last = space.getDOFs(0, 4)[0];
    const IndexMap<IndexSet> pairs{{last, {first}}};
    Body freeBody = Integral(Grad(u), Grad(v)) + Integral(u, v) - Integral(F::x, v);
    Body periodicBody = freeBody + PeriodicBC(u, pairs);
    Assembly::ProblemAssemblyInput<Body, decltype(u), decltype(v)> freeInput(
      freeBody, u, v);
    Assembly::ProblemAssemblyInput<Body, decltype(u), decltype(v)> periodicInput(
      periodicBody, u, v);
    System original;
    System constrained;
    assembler.execute(original, freeInput);
    assembler.execute(constrained, periodicInput);
    Math::Matrix<Real> prolongation = Math::Matrix<Real>::Zero(5, 4);
    prolongation.topRows(4).setIdentity();
    prolongation(last, first) = 1;
    const Math::Matrix<Real> matrix = original.getOperator();
    const Math::Matrix<Real> reduced = prolongation.transpose() * matrix * prolongation;
    const Math::Vector<Real> load = prolongation.transpose() * original.getVector();
    const Math::Vector<Real> expected = prolongation * reduced.fullPivLu().solve(load);
    const Math::Matrix<Real> periodicMatrix = constrained.getOperator();
    const Math::Vector<Real> actual =
      periodicMatrix.fullPivLu().solve(constrained.getVector());
    EXPECT_LT((actual - expected).norm(), 1e-12);
    EXPECT_NEAR(actual[first], actual[last], 1e-12);
  }

  /// @brief Periodic constraints on the second trial block must use its global offset.
  template <template <class, class> class Backend>
  void checkMixedPeriodicAssembly()
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {5});
    mesh.scale(0.25);
    mesh.getConnectivity().compute(0, 1);
    P1 space(mesh);
    P0g meanSpace(mesh);
    TrialFunction mean(meanSpace);
    TrialFunction u(space);
    TestFunction meanTest(meanSpace);
    TestFunction v(space);
    using System = Math::LinearSystem<Math::Matrix<Real>, Math::Vector<Real>>;
    using ProblemType =
      Problem<System, decltype(mean), decltype(u), decltype(meanTest), decltype(v)>;
    using Body = typename ProblemType::ProblemBodyType;
    const Index first = space.getDOFs(0, 0)[0];
    const Index last = space.getDOFs(0, 4)[0];
    const IndexMap<IndexSet> pairs{{last, {first}}};
    Body body = Integral(Grad(u), Grad(v)) + Integral(mean, v) + Integral(u, meanTest) -
      Integral(F::x - 0.5, v) + PeriodicBC(u, pairs);
    Tuple<std::reference_wrapper<decltype(mean)>, std::reference_wrapper<decltype(u)>>
      trials(std::ref(mean), std::ref(u));
    Tuple<std::reference_wrapper<decltype(meanTest)>, std::reference_wrapper<decltype(v)>>
      tests(std::ref(meanTest), std::ref(v));
    std::array<size_t, 2> trialOffsets{0, 1};
    std::array<size_t, 2> testOffsets{0, 1};
    boost::bimap<FormLanguage::Base::UUID, size_t> trialMap;
    boost::bimap<FormLanguage::Base::UUID, size_t> testMap;
    trialMap.insert({mean.getUUID(), 0});
    trialMap.insert({u.getUUID(), 1});
    testMap.insert({meanTest.getUUID(), 0});
    testMap.insert({v.getUUID(), 1});
    Assembly::ProblemAssemblyInput<Body, decltype(mean), decltype(u), decltype(meanTest),
      decltype(v)>
      input(body, trials, tests, trialOffsets, testOffsets, trialMap, testMap, 6, 6);
    Backend<System, ProblemType> assembler;
    System system;
    assembler.execute(system, input);
    const Math::Vector<Real> solution =
      system.getOperator().fullPivLu().solve(system.getVector());
    EXPECT_LT((system.getOperator() * solution - system.getVector()).norm(), 1e-12);
    EXPECT_NEAR(solution[1], solution[5], 1e-12);
    // The exact integral of a linear field on each segment is its endpoint average.
    const Real integral = 0.125 *
      (solution[1] + 2 * solution[2] + 2 * solution[3] + 2 * solution[4] + solution[5]);
    EXPECT_NEAR(integral, 0, 1e-12);
  }
}

/// @brief Dense periodic assembly agrees with algebraic condensation.
TEST(PeriodicConstraint, SequentialDense)
{
  checkPeriodicAssembly<Math::Matrix<Real>, Assembly::Sequential>();
}

/// @brief Sparse periodic assembly agrees with algebraic condensation.
TEST(PeriodicConstraint, SequentialSparse)
{
  checkPeriodicAssembly<Math::SparseMatrix<Real>, Assembly::Sequential>();
}

/// @brief Mixed periodic assembly offsets the constrained field correctly.
TEST(PeriodicConstraint, SequentialMixed)
{
  checkMixedPeriodicAssembly<Assembly::Sequential>();
}

#ifdef RODIN_USE_OPENMP
/// @brief OpenMP dense periodic assembly agrees with algebraic condensation.
TEST(PeriodicConstraint, OpenMPDense)
{
  checkPeriodicAssembly<Math::Matrix<Real>, Assembly::OpenMP>();
}

/// @brief OpenMP sparse periodic assembly agrees with algebraic condensation.
TEST(PeriodicConstraint, OpenMPSparse)
{
  checkPeriodicAssembly<Math::SparseMatrix<Real>, Assembly::OpenMP>();
}

/// @brief OpenMP mixed periodic assembly offsets the constrained field correctly.
TEST(PeriodicConstraint, OpenMPMixed)
{
  checkMixedPeriodicAssembly<Assembly::OpenMP>();
}
#endif
