/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Checks that optimized QuadratureRule handlers are the ones selected.
 *
 * An integrand no optimized handler matches is still integrated, by the generic
 * handler, and still gives the right answer. So a specialization that stops
 * matching --- because an expression type moved, or a template parameter
 * changed, or an operator now returns something slightly different --- breaks
 * nothing that any other test can see. Results stay correct, every suite stays
 * green, and only the running time moves.
 *
 * These assertions are what makes that visible. They are compile-time, so a
 * pattern falling out of use is a build failure rather than an unexplained
 * regression in the benchmarks.
 *
 * Both directions are asserted. An integrand that no handler specializes must
 * report *unspecialized*, because optimized and generic handlers share base
 * classes: were a @c Specialized member ever declared on
 * LocalBilinearFormIntegratorBase or LinearFormIntegratorBase, every handler
 * would inherit it, every assertion below would pass, and the suite would
 * detect nothing ever again. The negative cases are what catch that.
 */
#include <gtest/gtest.h>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Variational/H1.h>

#include "Rodin/FormLanguage/IsSpecialized.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  using LocalMesh = Mesh<Context::Local>;

  /// @brief The handler a given integrand selects.
  template <class Integral>
  using HandlerOf = QuadratureRule<typename Integral::IntegrandType>;

  /// @brief Whether that handler reports itself optimized.
  template <class Integral>
  constexpr bool selectsOptimized()
  {
    return FormLanguage::IsSpecialized<HandlerOf<Integral>>::Value;
  }
}

/**
 * @brief Each P1 expression with a specialization selects it.
 *
 * The expressions mirror those benchmarked under Integrator/Kernel/P1, so a
 * pattern that stops matching shows up here as a build failure rather than as
 * a silent slowdown there.
 */
TEST(QuadratureRuleDispatchTest, P1ExpressionsSelectTheirSpecializations)
{
  LocalMesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);

  P1<Real, LocalMesh> fes(mesh);
  P1<Math::SpatialVector<Real>, LocalMesh> vfes(mesh, 2);

  TrialFunction u(fes);
  TestFunction v(fes);
  TrialFunction w(vfes);
  TestFunction z(vfes);

  RealFunction c([](const Point&) { return 1.0; });
  VectorFunction beta{1.0, 1.0};

  static_assert(
    selectsOptimized<decltype(Integral(v))>(), "P1 basis load is no longer specialized");
  static_assert(selectsOptimized<decltype(Integral(c, v))>(),
    "P1 source load is no longer specialized");
  static_assert(
    selectsOptimized<decltype(Integral(u, v))>(), "P1 mass is no longer specialized");
  static_assert(selectsOptimized<decltype(Integral(c * u, v))>(),
    "P1 weighted mass is no longer specialized");
  static_assert(selectsOptimized<decltype(Integral(c * Dot(u, v)))>(),
    "P1 outer-weighted mass is no longer specialized");
  static_assert(selectsOptimized<decltype(Integral(Grad(u), Grad(v)))>(),
    "P1 gradient-gradient is no longer specialized");
  static_assert(selectsOptimized<decltype(Integral(c * Grad(u), Grad(v)))>(),
    "P1 weighted gradient-gradient is no longer specialized");
  static_assert(selectsOptimized<decltype(Integral(Jacobian(w), Jacobian(z)))>(),
    "P1 Jacobian-Jacobian is no longer specialized");
  static_assert(selectsOptimized<decltype(Integral(c * Jacobian(w), Jacobian(z)))>(),
    "P1 weighted Jacobian-Jacobian is no longer specialized");
  static_assert(selectsOptimized<decltype(Integral(Dot(Jacobian(w) * beta, z)))>(),
    "P1 Jacobian-advection is no longer specialized");
  static_assert(selectsOptimized<decltype(Integral(Div(w), v))>(),
    "P1 divergence-pressure is no longer specialized");
  static_assert(selectsOptimized<decltype(Integral(u, Div(z)))>(),
    "P1 pressure-divergence is no longer specialized");

  SUCCEED();
}

/**
 * @brief Each H1 expression with a specialization selects it, at every order.
 *
 * The handlers are templated on the polynomial order, so one order passing
 * does not establish the others: a specialization that accidentally required a
 * particular degree would still match at that degree and silently fall through
 * elsewhere. First, second and third are checked.
 *
 * The element is not a variable here. No handler signature mentions a
 * geometry, so which handler is selected cannot depend on one; the elements
 * matter to what the handlers compute, which is checked separately.
 */
template <size_t K>
static void checkH1Dispatch()
{
  LocalMesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);
  mesh.getConnectivity().compute(1, 2);

  H1<K, Real, LocalMesh> fes(std::integral_constant<size_t, K>{}, mesh);
  H1<K, Math::SpatialVector<Real>, LocalMesh> vfes(
    std::integral_constant<size_t, K>{}, mesh, 2);

  TrialFunction u(fes);
  TestFunction v(fes);
  TrialFunction w(vfes);
  TestFunction z(vfes);

  RealFunction c([](const Point&) { return 1.0; });

  static_assert(
    selectsOptimized<decltype(Integral(v))>(), "H1 basis load is no longer specialized");
  static_assert(selectsOptimized<decltype(Integral(c, v))>(),
    "H1 source load is no longer specialized");
  static_assert(
    selectsOptimized<decltype(Integral(u, v))>(), "H1 mass is no longer specialized");
  static_assert(selectsOptimized<decltype(Integral(c * u, v))>(),
    "H1 weighted mass is no longer specialized");
  static_assert(selectsOptimized<decltype(Integral(c * Dot(u, v)))>(),
    "H1 outer-weighted mass is no longer specialized");
  static_assert(selectsOptimized<decltype(Integral(Grad(u), Grad(v)))>(),
    "H1 gradient-gradient is no longer specialized");
  static_assert(selectsOptimized<decltype(Integral(c * Grad(u), Grad(v)))>(),
    "H1 weighted gradient-gradient is no longer specialized");
  static_assert(selectsOptimized<decltype(Integral(Jacobian(w), Jacobian(z)))>(),
    "H1 Jacobian-Jacobian is no longer specialized");
  static_assert(selectsOptimized<decltype(Integral(c * Jacobian(w), Jacobian(z)))>(),
    "H1 weighted Jacobian-Jacobian is no longer specialized");
  static_assert(selectsOptimized<decltype(Integral(Div(w), v))>(),
    "H1 divergence-pressure is no longer specialized");
  static_assert(selectsOptimized<decltype(Integral(u, Div(z)))>(),
    "H1 pressure-divergence is no longer specialized");
}

TEST(QuadratureRuleDispatchTest, H1ExpressionsSelectTheirSpecializations)
{
  checkH1Dispatch<1>();
  checkH1Dispatch<2>();
  checkH1Dispatch<3>();
  SUCCEED();
}

/**
 * @brief Expressions no handler specializes report themselves generic.
 *
 * Without these the suite cannot tell "the specialization fires" from "the
 * trait says yes to everything", which is the state it would be left in by a
 * @c Specialized member appearing on a shared base class.
 */
TEST(QuadratureRuleDispatchTest, UnspecializedExpressionsReportGeneric)
{
  LocalMesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);

  P1<Math::SpatialVector<Real>, LocalMesh> vfes(mesh, 2);
  TrialFunction w(vfes);
  TestFunction z(vfes);

  static_assert(!selectsOptimized<decltype(Integral(Div(w), Div(z)))>(),
    "divergence-divergence is expected to have no specialization; if one was "
    "added, assert it above instead of here");

  SUCCEED();
}
