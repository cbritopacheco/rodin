/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/// @file @brief CRTP callback adapters preserve owned data and constitutive references.
#include <gtest/gtest.h>
#include <memory>
#include <type_traits>

#include "Rodin/Solid/Local/Input.h"
#include "Rodin/Solid/Local/Output.h"
#include "Rodin/Solid/Constitutive/NeoHookean.h"
#include "Rodin/Solid/Constitutive/Yeoh.h"

using namespace Rodin;
using namespace Rodin::Solid;

namespace
{
  struct Injection : Input<Injection>
  {
    Injection(Real value, const Math::SpatialVector<Real>& direction)
      : activation(value),
        fiber(direction)
    {}

    void populate(ConstitutivePoint& point) const
    {
      point.set<Tags::Activation>(activation);
      point.set<Tags::FiberDirection>(fiber);
    }

    Real activation;
    Math::SpatialVector<Real> fiber;
  };

  struct Sample
  {
    Index cell, quadrature;
    Real jacobian, invariant;
    const ConstitutivePoint* point;
    const void* prepared;
  };

  struct Extraction : Output<Extraction>
  {
    explicit Extraction(const std::shared_ptr<List<Sample>>& destination)
      : samples(destination)
    {}

    template <class Law>
    void extract(const ConstitutivePoint& point, const Law& prepared) const
    {
      samples->push_back(
        {point.get<Tags::CellIndex>(), point.get<Tags::QuadraturePointIndex>(),
          prepared.getKinematicState().getJacobian(),
          prepared.getKinematicState().getRightCauchyGreenTensor().trace(), &point,
          &prepared});
    }

    std::shared_ptr<List<Sample>> samples;
  };
}

/// @brief Copied input callbacks retain their injection payload after the original is destroyed.
TEST(Solid_LocalCallbacks, InputOwnsPayloadAndMutatesOriginalPoint)
{
  for (size_t dimension : {2u, 3u})
  {
    Math::SpatialVector<Real> fiber(dimension);
    for (size_t i = 0; i < dimension; ++i)
      fiber(i) = Real(i + 1);
    InputFunction callback;
    {
      const Injection injection(Real(4.25), fiber);
      callback = [injection](ConstitutivePoint& point) {
        const Input<Injection>& base = injection;
        base.populate(point);
      };
    }
    auto copied = callback;
    callback = {};
    KinematicState state(dimension);
    ConstitutivePoint point(state);
    point.set<Tags::CellIndex>(17);
    copied(point);
    ASSERT_TRUE(point.has<Tags::Activation>());
    ASSERT_TRUE(point.has<Tags::FiberDirection>());
    EXPECT_DOUBLE_EQ(point.get<Tags::Activation>(), Real(4.25));
    EXPECT_EQ((point.get<Tags::FiberDirection>() - fiber).norm(), Real(0));
    EXPECT_EQ(point.get<Tags::CellIndex>(), 17u);
    EXPECT_EQ(&point.getKinematicState(), &state);
  }
}

template <class Law>
class Solid_LocalOutput : public testing::Test
{};

using Laws = testing::Types<NeoHookean, Yeoh>;
TYPED_TEST_SUITE(Solid_LocalOutput, Laws);

/// @brief Output adapters dispatch both law prepared types without copying point or prepared.
TYPED_TEST(Solid_LocalOutput, PreservesLawAndPointReferencesAfterCallbackCopy)
{
  for (size_t dimension : {2u, 3u})
  {
    Math::SpatialMatrix<Real> gradient(dimension, dimension);
    gradient.setZero();
    gradient(0, 0) = Real(0.2);
    gradient(1, 1) = Real(-0.1);
    if (dimension == 3)
      gradient(2, 2) = Real(0.1);
    KinematicState state(dimension);
    state.setDisplacementGradient(gradient);
    ConstitutivePoint point(state);
    auto prepared = [] {
      if constexpr (std::is_same_v<TypeParam, NeoHookean>)
      {
        return NeoHookean(Real(2), Real(3));
      }
      else
      {
        return Yeoh(Real(1), Real(0.2), Real(0.1), Real(10));
      }
    }();
    prepared.setConstitutivePoint(point);

    auto samples = std::make_shared<List<Sample>>();
    OutputFunction<TypeParam> callback;
    {
      const Extraction extraction(samples);
      callback = [extraction](const ConstitutivePoint& point, const TypeParam& prepared) {
        const Output<Extraction>& base = extraction;
        base.extract(point, prepared);
      };
    }
    auto copied = callback;
    callback = {};
    for (Index i = 0; i < 2; ++i)
    {
      point.set<Tags::CellIndex>(17 + i);
      point.set<Tags::QuadraturePointIndex>(3 + i);
      copied(point, prepared);
    }
    ASSERT_EQ(samples->size(), 2u);
    Index i = 0;
    for (const auto& sample : *samples)
    {
      EXPECT_EQ(sample.cell, 17 + i);
      EXPECT_EQ(sample.quadrature, 3 + i);
      EXPECT_EQ(sample.point, &point);
      EXPECT_EQ(sample.prepared, &prepared);
      EXPECT_NEAR(
        sample.jacobian, dimension == 2 ? Real(1.08) : Real(1.188), Real(1e-12));
      EXPECT_NEAR(
        sample.invariant, dimension == 2 ? Real(2.25) : Real(3.46), Real(1e-12));
      ++i;
    }
    EXPECT_EQ(point.get<Tags::CellIndex>(), 18u);
    EXPECT_EQ(point.get<Tags::QuadraturePointIndex>(), 4u);
  }
}
