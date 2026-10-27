/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Energy derivatives and mixed displacement-pressure element tangents. */
#include <gtest/gtest.h>
#include <memory>
#include <cmath>
#include <type_traits>
#include "Rodin/Variational.h"
#include "Rodin/Solid/Constitutive/Yeoh.h"
#include "Rodin/Solid/Integrators/InternalVirtualWorkResidual.h"
#include "Rodin/Solid/Integrators/InternalVirtualWorkTangent.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  /**
   * @brief Checks all mixed blocks against the independent isochoric energy.
   * @param space Displacement space, with its vector dimension matching the mesh dimension.
   */
  template <class Space>
  void checkMixedBlocks(const Space& space)
  {
    const auto& mesh = space.getMesh();
    const size_t dimension = mesh.getDimension();
    const auto cell = *mesh.getCell();
    const auto cellId = std::make_pair(dimension, cell.getIndex());
    P1 pressureSpace(mesh);
    TrialFunction u(space);
    TrialFunction p(pressureSpace);
    TestFunction v(space);
    TestFunction q(pressureSpace);
    GridFunction displacement(space);
    GridFunction pressure(pressureSpace);
    displacement = [dimension](const Point& x) {
      Math::SpatialVector<Real> value(dimension);
      for (size_t i = 0; i < dimension; ++i)
      {
        value[i] = 0;
        for (size_t j = 0; j < dimension; ++j)
          value[i] += Real(i == j ? 0.12 : 0.03) * x(j);
      }
      return value;
    };
    pressure = [](const Point& x) { return 0.7 + 0.2 * x(0); };
    // Zero bulk modulus leaves only the isochoric energy in the constitutive law.
    Solid::Yeoh law(1, 0, 0, 0);
    Solid::InternalVirtualWorkResidual ru(law, v, displacement, pressure);
    Solid::InternalVirtualWorkResidualP rp(q, displacement);
    Solid::InternalVirtualWorkTangent kuu(law, u, v, displacement, pressure);
    Solid::InternalVirtualWorkTangentUP kup(p, v, displacement);
    Solid::InternalVirtualWorkTangentPU kpu(u, q, displacement);
    const size_t nu = space.getFiniteElement(dimension, cell.getIndex()).getCount();
    const size_t np =
      pressureSpace.getFiniteElement(dimension, cell.getIndex()).getCount();
    const auto bind = [&] {
      ru.setPolytope(cell);
      rp.setPolytope(cell);
      kuu.setPolytope(cell);
      kup.setPolytope(cell);
      kpu.setPolytope(cell);
    };
    bind();
    Math::Vector<Real> residualU(nu), residualP(np);
    Math::Matrix<Real> uu(nu, nu), up(nu, np), pu(np, nu);
    for (size_t i = 0; i < nu; ++i)
    {
      residualU(i) = ru.integrate(i);
      for (size_t j = 0; j < nu; ++j)
        uu(i, j) = kuu.integrate(j, i);
      for (size_t j = 0; j < np; ++j)
        up(i, j) = kup.integrate(j, i);
    }
    for (size_t i = 0; i < np; ++i)
    {
      residualP(i) = rp.integrate(i);
      for (size_t j = 0; j < nu; ++j)
      {
        pu(i, j) = kpu.integrate(j, i);
        EXPECT_NEAR(pu(i, j), up(j, i), 1e-12);
      }
    }

    /** Independent density @f$ W=J^{-2/d}\operatorname{tr}(F^TF)-d+p(J-1) @f$. */
    const auto energy = [&] {
      const auto& formula = QF::QuadratureFormula::get(6, cell.getGeometry());
      const auto& quadrature = cell.getQuadrature(formula);
      auto gradient = Jacobian(displacement);
      Real result = 0;
      for (size_t i = 0; i < quadrature.getSize(); ++i)
      {
        const auto& point = quadrature.getPoint(i);
        const IntegrationPoint ip(point, &formula, i);
        const auto h = gradient.getValue(ip);
        Eigen::MatrixXd f = Eigen::MatrixXd::Identity(dimension, dimension);
        for (size_t row = 0; row < dimension; ++row)
        {
          for (size_t column = 0; column < dimension; ++column)
            f(row, column) += h(row, column);
        }
        const Real jacobian = f.determinant();
        const Real density = std::pow(jacobian, -Real(2) / dimension) * f.squaredNorm() -
          dimension + pressure.getValue(ip) * (jacobian - 1);
        result += formula.getWeight(i) * point.getDistortion() * density;
      }
      return result;
    };
    constexpr Real Step = 1e-6;
    constexpr Real DerivativeTolerance = 2e-7;
    for (size_t j = 0; j < nu; ++j)
    {
      const Index global = space.getGlobalIndex(cellId, j);
      const Real saved = displacement[global];
      displacement[global] = saved + Step;
      const Real plusEnergy = energy();
      ru.setPolytope(cell);
      rp.setPolytope(cell);
      Math::Vector<Real> plusU(nu), plusP(np);
      for (size_t i = 0; i < nu; ++i)
        plusU(i) = ru.integrate(i);
      for (size_t i = 0; i < np; ++i)
        plusP(i) = rp.integrate(i);
      displacement[global] = saved - Step;
      const Real minusEnergy = energy();
      EXPECT_NEAR(
        residualU(j), (plusEnergy - minusEnergy) / (2 * Step), DerivativeTolerance);
      ru.setPolytope(cell);
      rp.setPolytope(cell);
      for (size_t i = 0; i < nu; ++i)
      {
        EXPECT_NEAR(
          uu(i, j), (plusU(i) - ru.integrate(i)) / (2 * Step), DerivativeTolerance);
      }
      for (size_t i = 0; i < np; ++i)
      {
        EXPECT_NEAR(
          pu(i, j), (plusP(i) - rp.integrate(i)) / (2 * Step), DerivativeTolerance);
      }
      displacement[global] = saved;
    }
    for (size_t j = 0; j < np; ++j)
    {
      const Index global = pressureSpace.getGlobalIndex(cellId, j);
      const Real saved = pressure[global];
      pressure[global] = saved + Step;
      const Real plusEnergy = energy();
      ru.setPolytope(cell);
      Math::Vector<Real> plus(nu);
      for (size_t i = 0; i < nu; ++i)
        plus(i) = ru.integrate(i);
      pressure[global] = saved - Step;
      const Real minusEnergy = energy();
      EXPECT_NEAR(
        residualP(j), (plusEnergy - minusEnergy) / (2 * Step), DerivativeTolerance);
      ru.setPolytope(cell);
      for (size_t i = 0; i < nu; ++i)
      {
        EXPECT_NEAR(
          up(i, j), (plus(i) - ru.integrate(i)) / (2 * Step), DerivativeTolerance);
      }
      pressure[global] = saved;
    }
    // Clones retain borrowed states and can be configured and rebound independently.
    GridFunction replacement(displacement);
    GridFunction replacementPressure(pressure);
    size_t callbackCalls = 0;
    const auto input = [&](Solid::ConstitutivePoint& point) {
      ++callbackCalls;
      EXPECT_EQ(point.get<Solid::Tags::CellIndex>(), cell.getIndex());
      EXPECT_TRUE(point.has<Solid::Tags::QuadraturePointIndex>());
    };
    const auto clone = [&](auto& original) {
      std::unique_ptr<std::remove_reference_t<decltype(original)>> copied(
        original.copy());
      EXPECT_EQ(copied->getRegion(), Region::Cells);
      EXPECT_EQ(&copied->setQuadratureOrder(6), copied.get());
      EXPECT_EQ(&copied->setDisplacement(replacement), copied.get());
      if constexpr (requires { copied->setPressure(replacementPressure); })
      {
        EXPECT_EQ(&copied->setPressure(replacementPressure), copied.get());
      }
      if constexpr (requires { copied->setInput(input); })
      {
        EXPECT_EQ(&copied->setInput(input), copied.get());
      }
      copied->setPolytope(cell);
      EXPECT_EQ(&copied->getPolytope(), &cell);
      return copied;
    };
    auto ruCopy = clone(ru);
    auto rpCopy = clone(rp);
    auto uuCopy = clone(kuu);
    auto upCopy = clone(kup);
    auto puCopy = clone(kpu);
    EXPECT_GT(callbackCalls, 0u);
    for (size_t i = 0; i < nu; ++i)
    {
      EXPECT_NEAR(ruCopy->integrate(i), residualU(i), 1e-12);
      for (size_t j = 0; j < nu; ++j)
        EXPECT_NEAR(uuCopy->integrate(j, i), uu(i, j), 1e-12);
      for (size_t j = 0; j < np; ++j)
        EXPECT_NEAR(upCopy->integrate(j, i), up(i, j), 1e-12);
    }
    for (size_t i = 0; i < np; ++i)
    {
      EXPECT_NEAR(rpCopy->integrate(i), residualP(i), 1e-12);
      for (size_t j = 0; j < nu; ++j)
        EXPECT_NEAR(puCopy->integrate(j, i), pu(i, j), 1e-12);
    }
  }
}

/// @brief Linear mixed blocks obey energy derivatives in two and three dimensions.
TEST(MixedInternalVirtualWork, LinearDisplacement)
{
  for (const auto geometry : {Polytope::Type::Triangle, Polytope::Type::Tetrahedron})
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    auto mesh = LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 2));
    P1 space(mesh, dimension);
    checkMixedBlocks(space);
  }
}

/// @brief Quadratic displacement modes obey the same full mixed linearization.
TEST(MixedInternalVirtualWork, QuadraticDisplacement)
{
  for (const auto geometry : {Polytope::Type::Triangle, Polytope::Type::Tetrahedron})
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    auto mesh = LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 2));
    for (size_t d = dimension; d > 0; --d)
      mesh.getConnectivity().compute(d, d - 1);
    H1 space(std::integral_constant<size_t, 2>{}, mesh, dimension);
    checkMixedBlocks(space);
  }
}
