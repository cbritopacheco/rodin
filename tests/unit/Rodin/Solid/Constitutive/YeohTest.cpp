/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/// @file @brief Independent energy and finite-difference oracles for Yeoh material.
#include <gtest/gtest.h>
#include "Rodin/Solid/Constitutive/Yeoh.h"

using namespace Rodin;
using namespace Rodin::Solid;

namespace
{
  struct Response
  {
    Real energy;
    Math::SpatialMatrix<Real> stress, tangent;
  };

  Response evaluate(const HyperElasticLaw<Yeoh>& law, const Math::SpatialMatrix<Real>& h,
    const Math::SpatialMatrix<Real>& increment)
  {
    KinematicState state(h.rows());
    state.setDisplacementGradient(h);
    ConstitutivePoint point(state);
    Yeoh::Cache cache;
    law.setCache(cache, point);
    Response result;
    result.energy = law.getStrainEnergyDensity(cache, point);
    law.getFirstPiolaKirchhoffStress(result.stress, cache, point);
    law.getMaterialTangent(result.tangent, cache, point, increment);
    return result;
  }
}

TEST(Rodin_Solid_Yeoh, EnergyStressAndTangentInTwoAndThreeDimensions)
{
  const Yeoh law(2, 0.7, 0.3, 11);
  EXPECT_DOUBLE_EQ(law.getMaterialConstantC1(), 2);
  EXPECT_DOUBLE_EQ(law.getMaterialConstantC2(), 0.7);
  EXPECT_DOUBLE_EQ(law.getMaterialConstantC3(), 0.3);
  EXPECT_DOUBLE_EQ(law.getBulkModulus(), 11);
  for (size_t dimension : {2, 3})
  {
    SCOPED_TRACE(dimension);
    Math::SpatialMatrix<Real> h(dimension, dimension), increment(dimension, dimension);
    for (size_t row = 0; row < dimension; ++row)
    {
      for (size_t column = 0; column < dimension; ++column)
      {
        h(row, column) = row == column ? 0.1 * (row + 1) : 0.04 * (Real(row) - column);
        increment(row, column) = 0.13 * (row + 1) - 0.07 * column;
      }
    }
    const auto actual = evaluate(law, h, increment);
    Eigen::MatrixXd f(dimension, dimension);
    for (size_t row = 0; row < dimension; ++row)
    {
      for (size_t column = 0; column < dimension; ++column)
        f(row, column) = h(row, column) + (row == column ? 1 : 0);
    }
    const Real jacobian = f.determinant();
    const Real invariant =
      std::pow(jacobian, -2.0 / dimension) * f.squaredNorm() - dimension;
    const Real expected = 2 * invariant + 0.7 * invariant * invariant +
      0.3 * invariant * invariant * invariant + 5.5 * (jacobian - 1) * (jacobian - 1);
    EXPECT_NEAR(actual.energy, expected, 1e-12);
    const Real step = 1e-6;
    for (size_t row = 0; row < dimension; ++row)
    {
      for (size_t column = 0; column < dimension; ++column)
      {
        auto plus = h, minus = h;
        plus(row, column) += step;
        minus(row, column) -= step;
        const Real derivative = (evaluate(law, plus, increment).energy -
                                  evaluate(law, minus, increment).energy) /
          (2 * step);
        EXPECT_NEAR(actual.stress(row, column), derivative, 1e-7);
      }
    }
    const auto plus = evaluate(law, h + step * increment, increment);
    const auto minus = evaluate(law, h - step * increment, increment);
    EXPECT_LT((actual.tangent - (plus.stress - minus.stress) / (2 * step)).norm(), 1e-7);
    Yeoh copied(law);
    Yeoh moved(std::move(copied));
    EXPECT_DOUBLE_EQ(evaluate(moved, h, increment).energy, actual.energy);
    h.setZero();
    const auto rest = evaluate(law, h, increment);
    EXPECT_NEAR(rest.energy, 0, 1e-14);
    EXPECT_LT(rest.stress.norm(), 1e-13);
    h.setIdentity();
    h *= 0.2;
    const auto dilation = evaluate(law, h, increment);
    const Real j = std::pow(1.2, dimension);
    EXPECT_NEAR(dilation.energy, 5.5 * (j - 1) * (j - 1), 1e-12);
    for (size_t row = 0; row < dimension; ++row)
    {
      for (size_t column = 0; column < dimension; ++column)
      {
        EXPECT_NEAR(dilation.stress(row, column),
          row == column ? 11 * (j - 1) * j / 1.2 : 0, 1e-12);
      }
    }
  }
}
