/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/// @file @brief Partial linear-system updates preserve nonzero boundary lifting.
#include <gtest/gtest.h>
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"
#include "Rodin/Variational/H1.h"
#include "Rodin/Assembly/Sequential.h"
#ifdef RODIN_USE_OPENMP
#include "Rodin/Assembly/OpenMP.h"
#endif

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  /// @brief Solves the exact mass problem after refreshing only one assembled side.
  template <bool HigherOrder, template <class, class> class Backend>
  void checkTargetedBoundary(AssemblyTarget target)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {5});
    mesh.scale(Real(0.25));
    mesh.getConnectivity().compute(0, 1);
    auto space = [&] {
      if constexpr (HigherOrder)
        return H1<2, Real>(std::integral_constant<size_t, 2>{}, mesh);
      else
        return P1<Real>(mesh);
    }();
    TrialFunction u(space);
    TestFunction v(space);
    RealFunction prescribed([](const Point& point) { return 1 + 2 * point.x(); });
    GridFunction exact(space);
    exact = prescribed;
    using System = Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;
    using ProblemType = Problem<System, decltype(u), decltype(v)>;
    using Body = typename ProblemType::ProblemBodyType;
    Body body = Integral(u, v) - Integral(prescribed, v) + DirichletBC(u, prescribed);
    Assembly::ProblemAssemblyInput<Body, decltype(u), decltype(v)> input(body, u, v);
    Backend<System, ProblemType> assembler;
    if constexpr (requires { assembler.setThreadCount(size_t(2)); })
      assembler.setThreadCount(2);
    System system;
    assembler.execute(system, input);
    const Math::Matrix<Real> fullMatrix = system.getOperator();
    const Math::Vector<Real> fullLoad = system.getVector();
    ASSERT_LT((fullMatrix * exact.getData() - fullLoad).norm(), 1e-12);
    assembler.execute(system, input, target);
    const Math::Matrix<Real> actual = system.getOperator();
    EXPECT_LT((actual - fullMatrix).norm(), 1e-12);
    EXPECT_LT((system.getVector() - fullLoad).norm(), 1e-12);
    const Math::Vector<Real> solution = actual.fullPivLu().solve(system.getVector());
    EXPECT_LT((solution - exact.getData()).norm(), 1e-11);
    EXPECT_LT((actual * exact.getData() - system.getVector()).norm(), 1e-12);
  }
}

/// @brief A sequential P1 RHS refresh retains boundary lifting from full assembly.
TEST(TargetedBoundary, SequentialP1RHS)
{
  checkTargetedBoundary<false, Assembly::Sequential>(AssemblyTarget::RHS);
}

/// @brief A sequential P1 LHS refresh remains compatible with the retained full load.
TEST(TargetedBoundary, SequentialP1LHS)
{
  checkTargetedBoundary<false, Assembly::Sequential>(AssemblyTarget::LHS);
}

/// @brief A sequential H1 RHS refresh retains boundary lifting from full assembly.
TEST(TargetedBoundary, SequentialH1RHS)
{
  checkTargetedBoundary<true, Assembly::Sequential>(AssemblyTarget::RHS);
}

/// @brief A sequential H1 LHS refresh remains compatible with the retained full load.
TEST(TargetedBoundary, SequentialH1LHS)
{
  checkTargetedBoundary<true, Assembly::Sequential>(AssemblyTarget::LHS);
}

#ifdef RODIN_USE_OPENMP
/// @brief OpenMP P1 RHS refresh is a compatible-boundary control.
TEST(TargetedBoundary, OpenMPP1RHS)
{
  checkTargetedBoundary<false, Assembly::OpenMP>(AssemblyTarget::RHS);
}

/// @brief OpenMP P1 LHS refresh is a compatible-boundary control.
TEST(TargetedBoundary, OpenMPP1LHS)
{
  checkTargetedBoundary<false, Assembly::OpenMP>(AssemblyTarget::LHS);
}

/// @brief OpenMP H1 RHS refresh is a compatible-boundary control.
TEST(TargetedBoundary, OpenMPH1RHS)
{
  checkTargetedBoundary<true, Assembly::OpenMP>(AssemblyTarget::RHS);
}

/// @brief OpenMP H1 LHS refresh is a compatible-boundary control.
TEST(TargetedBoundary, OpenMPH1LHS)
{
  checkTargetedBoundary<true, Assembly::OpenMP>(AssemblyTarget::LHS);
}
#endif
