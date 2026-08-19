/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include <gtest/gtest.h>

#include "Rodin/Adaptation/WNGIRValidationWeights.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"

using namespace Rodin;

namespace Rodin::Tests::Unit
{
  class WNGIRValidationWeightsTest
    : public testing::TestWithParam<Geometry::Polytope::Type>
  {};

  /// @brief Test case.
  TEST_P(WNGIRValidationWeightsTest, DefinesPositiveMeasureWithReferenceMass)
  {
    const auto& qf = QF::PolytopeQuadratureFormula::get(4, GetParam());
    const Adaptation::Detail::WNGIRValidationWeights weights(qf);
    Real signedWeight = 0;
    Real validationWeight = 0;
    bool hasNegativeWeight = false;
    for (std::size_t q = 0; q < qf.getSize(); ++q)
    {
      signedWeight += qf.getWeight(q);
      validationWeight += weights.getWeight(q);
      hasNegativeWeight |= qf.getWeight(q) < Real(0);
      EXPECT_GE(weights.getWeight(q), Real(0));
      EXPECT_NEAR(qf.getWeight(q) *
          Adaptation::Detail::WNGIRValidationWeights::getCorrection(qf, q),
        weights.getWeight(q), 1e-14);
    }
    EXPECT_TRUE(hasNegativeWeight);
    EXPECT_NEAR(validationWeight, signedWeight, 1e-14);
  }

  /// @brief Test case.
  INSTANTIATE_TEST_SUITE_P(Simplex, WNGIRValidationWeightsTest,
    testing::Values(
      Geometry::Polytope::Type::Triangle, Geometry::Polytope::Type::Tetrahedron));
}
