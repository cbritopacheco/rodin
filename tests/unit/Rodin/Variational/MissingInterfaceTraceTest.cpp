/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Missing interface trace diagnostics and subsequent valid rebinding. */
#include <gtest/gtest.h>
#include <string>
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
#ifdef NDEBUG
  constexpr bool AssertionsEnabled = false;
#else
  constexpr bool AssertionsEnabled = true;
#endif

  /**
   * @brief Checks the missing-trace diagnostic, then reuses the expression with a valid trace.
   * @tparam AssertBeforeException Whether Debug asserts before the release diagnostic branch.
   * @param operation Derivative expression to evaluate and rebind.
   * @param point Interior face point shared by the two attributed cells.
   */
  template <bool AssertBeforeException, class Operator>
  void checkMissingTrace(Operator& operation, const Point& point)
  {
    EXPECT_TRUE(operation.getTraceDomain().empty());
    // P0/P1 keep a Debug assertion before their explicit release exception.
    // H1 diagnoses missing input in both configurations. Accuracy of the
    // recovered derivatives is covered by DerivativeTraceTest and GradTest.
    if constexpr (!AssertBeforeException || !AssertionsEnabled)
    {
      try
      {
        operation.getValue(point);
        FAIL() << "An interface derivative requires a trace domain.";
      }
      catch (const Alert::Exception& error)
      {
        const std::string message = error.what();
        EXPECT_NE(
          message.find("In member function interpolate of class"), std::string::npos);
        EXPECT_NE(message.find("No trace domain provided"), std::string::npos);
        EXPECT_NE(message.find("getTraceDomain().size() == 0"), std::string::npos);
      }
    }
    EXPECT_TRUE(operation.getTraceDomain().empty());
    operation.traceOf(Attribute(11));
    EXPECT_TRUE(operation.getTraceDomain().contains(11));
    EXPECT_NO_THROW(operation.getValue(point));
  }

  /// @brief Supplies a shared interface and explicit attributes on both neighboring cells.
  LocalMesh makeTraceMesh()
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    mesh.getConnectivity().compute(1, 2);
    for (auto cell = mesh.getCell(); cell; ++cell)
      mesh.setAttribute({2, cell->getIndex()}, 11 + cell->getIndex());
    return mesh;
  }
}

/// @brief P0/P1 release diagnostics preserve field values and allow subsequent trace selection.
TEST(MissingInterfaceTrace, P0AndP1ReleaseDiagnosticsAndRecovery)
{
  auto mesh = makeTraceMesh();
  const auto face = *mesh.getInterface();
  const Point point(face, Polytope::Traits(face.getGeometry()).getCentroid());
  P0<Real> constants(mesh);
  P1<Real> linears(mesh);
  GridFunction constantField(constants);
  GridFunction linearField(linears);
  constantField = Real(7);
  linearField = [](const Point& p) { return 2 * p(0) - 3 * p(1); };
  const Math::Vector<Real> savedConstants = constantField.getData();
  const Math::Vector<Real> savedLinears = linearField.getData();
  auto constantGradient = Grad(constantField);
  auto linearGradient = Grad(linearField);
  Derivative derivative(0, linearField);
  checkMissingTrace<true>(constantGradient, point);
  checkMissingTrace<true>(linearGradient, point);
  checkMissingTrace<true>(derivative, point);
  EXPECT_EQ(constantField.getData(), savedConstants);
  EXPECT_EQ(linearField.getData(), savedLinears);
}

/// @brief Quadratic H1 rejects absent traces for gradient, partial derivative and Jacobian.
TEST(MissingInterfaceTrace, QuadraticH1DiagnosticsAndRecovery)
{
  auto mesh = makeTraceMesh();
  const auto face = *mesh.getInterface();
  const Point point(face, Polytope::Traits(face.getGeometry()).getCentroid());
  H1 scalars(std::integral_constant<size_t, 2>{}, mesh);
  H1 vectors(std::integral_constant<size_t, 2>{}, mesh, 2);
  GridFunction scalarField(scalars);
  GridFunction vectorField(vectors);
  scalarField = [](const Point& p) { return 2 * p(0) - 3 * p(1); };
  vectorField = [](const Point& p) { return Math::SpatialVector<Real>{p(0), 2 * p(1)}; };
  const Math::Vector<Real> savedScalars = scalarField.getData();
  const Math::Vector<Real> savedVectors = vectorField.getData();
  auto gradient = Grad(scalarField);
  Derivative derivative(0, scalarField);
  auto jacobian = Jacobian(vectorField);
  checkMissingTrace<false>(gradient, point);
  checkMissingTrace<false>(derivative, point);
  checkMissingTrace<false>(jacobian, point);
  EXPECT_EQ(scalarField.getData(), savedScalars);
  EXPECT_EQ(vectorField.getData(), savedVectors);
}
