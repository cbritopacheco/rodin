/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/// @file @brief Physical-coordinate oracles for weighted scalar diffusion.
#include <gtest/gtest.h>
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"
#include "Rodin/Variational/H1.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  /// @brief Compares every entry with analytic physical gradients of linear bases.
  template <bool HigherOrder, bool Tensor>
  void checkPhysicalDiffusion(bool embedded, bool stretched, bool symmetric)
  {
    const Real width = stretched ? Real(2) : Real(1);
    auto mesh = [&] {
      if (embedded)
      {
        return LocalMesh::Builder()
          .initialize(3)
          .nodes(3)
          .vertex({0, 0, 0})
          .vertex({width, 0, 0})
          .vertex({0, 1, 1})
          .polytope(Polytope::Type::Triangle, {0, 1, 2})
          .finalize();
      }
      return LocalMesh::Builder()
        .initialize(2)
        .nodes(3)
        .vertex({0, 0})
        .vertex({width, 0})
        .vertex({0, 1})
        .polytope(Polytope::Type::Triangle, {0, 1, 2})
        .finalize();
    }();
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    auto space = [&] {
      if constexpr (HigherOrder)
        return H1<1, Real>(std::integral_constant<size_t, 1>{}, mesh);
      else
        return P1<Real>(mesh);
    }();
    TrialFunction u(space);
    TestFunction v(space);
    Math::SpatialMatrix<Real> tensor(embedded ? 3 : 2, embedded ? 3 : 2);
    tensor.setIdentity();
    tensor(0, 0) = 2;
    tensor(0, 1) = symmetric ? Real(0.5) : Real(3);
    tensor(1, 0) = symmetric ? Real(0.5) : Real(0);
    auto coefficient = [&] {
      if constexpr (Tensor)
        return MatrixFunction(tensor);
      else
        return RealFunction(Real(1));
    }();
    auto integral = Integral(coefficient * Grad(u), Grad(v));
    const auto cell = mesh.getCell();
    integral.setPolytope(*cell);
    const auto& element = space.getFiniteElement(2, cell->getIndex());
    const Math::SpatialPoint reference{Real(1) / 3, Real(1) / 3};
    const Real area = width * (embedded ? std::sqrt(Real(2)) : Real(1)) / 2;
    std::vector<Math::SpatialVector<Real>> gradients;
    for (size_t local = 0; local < element.getCount(); ++local)
    {
      auto& gradient = gradients.emplace_back(embedded ? 3 : 2);
      gradient[0] =
        element.getBasis(local).template getDerivative<1>(0)(reference) / width;
      gradient[1] = element.getBasis(local).template getDerivative<1>(1)(reference);
      if (embedded)
      {
        // F(r,s) = (width*r,s,s): the tangential s-gradient is (0,1/2,1/2).
        gradient[1] /= 2;
        gradient[2] = gradient[1];
      }
    }
    for (size_t test = 0; test < element.getCount(); ++test)
    {
      for (size_t trial = 0; trial < element.getCount(); ++trial)
      {
        SCOPED_TRACE(::testing::Message() << "test=" << test << ", trial=" << trial);
        const Real expected = [&] {
          if constexpr (Tensor)
            return area * Math::dot(tensor * gradients[trial], gradients[test]);
          else
            return area * Math::dot(gradients[trial], gradients[test]);
        }();
        EXPECT_NEAR(integral.integrate(trial, test), expected, 1e-12);
      }
    }
  }
}

/// @brief A nonsymmetric tensor must retain its independent test/trial entries.
TEST(PhysicalDiffusion, P1NonsymmetricTensorOnReferenceTriangle)
{
  checkPhysicalDiffusion<false, true>(false, false, false);
}

/// @brief Physical tensors must be pulled back on both sides of the gradient contraction.
TEST(PhysicalDiffusion, P1NonsymmetricTensorOnStretchedTriangle)
{
  checkPhysicalDiffusion<false, true>(false, true, false);
}

/// @brief Symmetric positive-definite anisotropy also requires the tensor pullback.
TEST(PhysicalDiffusion, P1SymmetricTensorOnStretchedTriangle)
{
  checkPhysicalDiffusion<false, true>(false, true, true);
}

/// @brief The H1 physical-gradient path is an independent full-dimensional control.
TEST(PhysicalDiffusion, H1NonsymmetricTensorOnStretchedTriangle)
{
  checkPhysicalDiffusion<true, true>(false, true, false);
}

/// @brief Isotropic P1 diffusion retains all ambient directions through its metric.
TEST(PhysicalDiffusion, P1ScalarOnTiltedTriangle)
{
  checkPhysicalDiffusion<false, false>(true, true, false);
}

/// @brief Weighted H1 diffusion must retain the third tangential-gradient component.
TEST(PhysicalDiffusion, H1ScalarOnTiltedTriangle)
{
  checkPhysicalDiffusion<true, false>(true, true, false);
}

/// @brief Ambient tensors act on three physical gradient components of embedded P1 cells.
TEST(PhysicalDiffusion, P1TensorOnTiltedTriangle)
{
  checkPhysicalDiffusion<false, true>(true, true, false);
}

/// @brief Ambient tensors act on three physical gradient components of embedded H1 cells.
TEST(PhysicalDiffusion, H1TensorOnTiltedTriangle)
{
  checkPhysicalDiffusion<true, true>(true, true, false);
}
