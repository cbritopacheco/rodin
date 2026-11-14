/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Configure.h"

// Exercise core-header isolation even in CPU-only CI: these feature macros
// must not make the core solver umbrella include a vendor GPU implementation.
#ifndef RODIN_USE_CUDSS
#define RODIN_USE_CUDSS
#endif
#ifndef RODIN_USE_ROCSOLVER
#define RODIN_USE_ROCSOLVER
#endif

#include <gtest/gtest.h>

#include "Rodin/Solver.h"
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Solver_BackendHeaders, NativeSolveWithGpuFeaturesEnabled)
  {
    using namespace Rodin::Geometry;
    using namespace Rodin::Variational;
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {3, 3});
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    Problem pb(u, v);
    pb = Integral(u, v) - Integral(RealFunction(1), v);
    Solver::SparseLU solver(pb);
    solver.solve();
    ASSERT_TRUE(solver.success());
    EXPECT_LE((u.getSolution().getData().array() - Real(1)).matrix().norm(), Real(1e-12));
  }
}
