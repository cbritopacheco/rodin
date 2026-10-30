/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/// @file @brief Constitutive binding owns inputs and invalidates changed responses.
#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <thread>
#include <type_traits>
#include <Rodin/Solid.h>
#include "../../../../../examples/Heart/CoronaryArtery/GradedYeoh.h"

namespace Rodin::Tests::Unit
{
  namespace
  {
    template <class Law>
    concept PublicCache = requires { typename Law::Cache; };

    template <class Law>
    concept ConstBinding = requires(const Law& law,
      const Solid::ConstitutivePoint& point) { law.setConstitutivePoint(point); };

    template <class Law>
    Law material()
    {
      if constexpr (std::is_same_v<Law, Solid::NeoHookean> ||
        std::is_same_v<Law, Solid::SaintVenantKirchhoff>)
      {
        return Law(2, 3);
      }
      else if constexpr (std::is_same_v<Law, Solid::MooneyRivlin>)
      {
        return Law(1, 0.4, 12);
      }
      else if constexpr (std::is_same_v<Law, Solid::Yeoh>)
      {
        return Law(1, 0.2, 0.1, 12);
      }
      else
      {
        return Law(1, 0.3, 0.4, 0.2, 0.7, 0.3, 12);
      }
    }

    Math::SpatialMatrix<Real> gradient(size_t dimension, Real value)
    {
      Math::SpatialMatrix<Real> result(dimension, dimension);
      for (size_t i = 0; i < dimension; ++i)
      {
        for (size_t j = 0; j < dimension; ++j)
        {
          result(i, j) = i == j ? value * (i + 1) : value * (Real(i) - Real(j)) / 4;
        }
      }
      return result;
    }

    TEST(ConstitutiveBinding, SpatialMultiplierIsCapturedAgainOnEveryBinding)
    {
      auto mesh =
        Geometry::LocalMesh::UniformGrid(Geometry::Polytope::Type::Triangle, {2, 2});
      const auto cell = *mesh.getCell();
      Geometry::Point location(
        cell, Geometry::Polytope::Traits(cell.getGeometry()).getCentroid());
      Real scale = 2;
      size_t calls = 0;
      Examples::Heart::GradedYeoh law(
        material<Solid::Yeoh>(), [&](const Geometry::Point&) {
          ++calls;
          return scale;
        });
      Solid::KinematicState state(2);
      state.setDisplacementGradient(gradient(2, 0.08));
      Solid::ConstitutivePoint point(location, state);
      law.setConstitutivePoint(point);
      const Real original = law.getStrainEnergyDensity();
      auto retained = law;
      Math::SpatialMatrix<Real> stress, tangent;
      law.getFirstPiolaKirchhoffStress(stress);
      law.getMaterialTangent(tangent, gradient(2, 0.01));
      EXPECT_EQ(calls, 1u);
      scale = 3;
      EXPECT_DOUBLE_EQ(law.getStrainEnergyDensity(), original);
      law.setConstitutivePoint(point);
      EXPECT_EQ(calls, 2u);
      EXPECT_DOUBLE_EQ(law.getStrainEnergyDensity(), original * 1.5);
      EXPECT_DOUBLE_EQ(retained.getStrainEnergyDensity(), original);
      law.setConstitutivePoint(Solid::ConstitutivePoint(state));
      EXPECT_EQ(calls, 2u);
      EXPECT_DOUBLE_EQ(law.getStrainEnergyDensity(), original / 2);
      Solid::ActiveContraction active(law, Solid::ActiveFiberLaw{});
      active.setConstitutivePoint(point);
      const Real activeEnergy = active.getStrainEnergyDensity();
      scale = 4;
      active.setConstitutivePoint(point);
      EXPECT_EQ(calls, 4u);
      EXPECT_NEAR(active.getStrainEnergyDensity() - activeEnergy, original / 2, 1e-14);
    }

    template <class Law>
    class ConstitutiveBinding : public ::testing::Test
    {};

    using PassiveLaws = ::testing::Types<Solid::NeoHookean, Solid::MooneyRivlin,
      Solid::Yeoh, Solid::SaintVenantKirchhoff, Solid::HolzapfelOgden>;
    TYPED_TEST_SUITE(ConstitutiveBinding, PassiveLaws);

    TYPED_TEST(ConstitutiveBinding, OwnsInputsAcrossPointLifetimeCopiesAndRebinding)
    {
      static_assert(!PublicCache<TypeParam>);
      static_assert(!ConstBinding<TypeParam>);
      auto law = material<TypeParam>();
      Math::SpatialMatrix<Real> stress;
      EXPECT_THROW(law.getFirstPiolaKirchhoffStress(stress), Alert::Exception);
      for (size_t dimension : {2u, 3u})
      {
        {
          Solid::KinematicState state(dimension);
          state.setDisplacementGradient(gradient(dimension, 0.08));
          Solid::ConstitutivePoint point(state);
          Solid::HyperElasticLaw<TypeParam>& base = law;
          base.setConstitutivePoint(point);
          const Solid::HyperElasticLaw<TypeParam>& prepared = base;
          prepared.getFirstPiolaKirchhoffStress(stress);
          EXPECT_TRUE(std::isfinite(prepared.getStrainEnergyDensity()));
        }
        Math::SpatialMatrix<Real> retained;
        law.getFirstPiolaKirchhoffStress(retained);
        EXPECT_LT((stress.getData() - retained.getData()).norm(), 1e-14);
        const auto snapshot = law;
        Solid::KinematicState state(dimension);
        state.setDisplacementGradient(gradient(dimension, 0.08));
        Solid::ConstitutivePoint point(state);
        law.setConstitutivePoint(point);
        law.getFirstPiolaKirchhoffStress(retained);
        EXPECT_LT((stress.getData() - retained.getData()).norm(), 1e-14);
        // Mutating the borrowed input affects the law only at the next binding.
        state.setDisplacementGradient(gradient(dimension, 0.13));
        law.getFirstPiolaKirchhoffStress(retained);
        EXPECT_LT((stress.getData() - retained.getData()).norm(), 1e-14);
        law.setConstitutivePoint(point);
        law.getFirstPiolaKirchhoffStress(retained);
        EXPECT_GT((stress.getData() - retained.getData()).norm(), 1e-5);
        auto fresh = material<TypeParam>();
        fresh.setConstitutivePoint(point);
        Math::SpatialMatrix<Real> expected, increment = gradient(dimension, 0.03);
        fresh.getFirstPiolaKirchhoffStress(expected);
        EXPECT_LT((expected.getData() - retained.getData()).norm(), 1e-14);
        law.getMaterialTangent(retained, increment);
        fresh.getMaterialTangent(expected, increment);
        EXPECT_LT((expected.getData() - retained.getData()).norm(), 1e-14);
        snapshot.getFirstPiolaKirchhoffStress(expected);
        EXPECT_LT((stress.getData() - expected.getData()).norm(), 1e-14);
      }
    }

    TYPED_TEST(ConstitutiveBinding, IndependentWorkerCopiesDoNotSharePreparedState)
    {
      const auto prototype = material<TypeParam>();
      std::array<Real, 2> energies{};
      // Two workers exercise independent laws; numerical test cases remain serial.
      std::thread first([&] {
        auto local = prototype;
        Solid::KinematicState state(3);
        state.setDisplacementGradient(gradient(3, 0.04));
        local.setConstitutivePoint(Solid::ConstitutivePoint(state));
        energies[0] = local.getStrainEnergyDensity();
      });
      std::thread second([&] {
        auto local = prototype;
        Solid::KinematicState state(3);
        state.setDisplacementGradient(gradient(3, 0.12));
        local.setConstitutivePoint(Solid::ConstitutivePoint(state));
        energies[1] = local.getStrainEnergyDensity();
      });
      first.join();
      second.join();
      EXPECT_NE(energies[0], energies[1]);
      EXPECT_THROW(prototype.getStrainEnergyDensity(), Alert::Exception);
    }

    TEST(ConstitutiveBinding, FiberDirectionChangesAtIdenticalKinematics)
    {
      auto law = material<Solid::HolzapfelOgden>();
      Solid::KinematicState state(3);
      state.setDisplacementGradient(gradient(3, 0.12));
      Solid::ConstitutivePoint point(state);
      Math::SpatialVector<Real> fiber(3);
      fiber.setZero();
      fiber(0) = 1;
      point.set<Solid::Tags::FiberDirection>(fiber);
      law.setConstitutivePoint(point);
      const Real energy = law.getStrainEnergyDensity();
      fiber(0) = 0;
      fiber(1) = 1;
      point.set<Solid::Tags::FiberDirection>(fiber);
      law.setConstitutivePoint(point);
      EXPECT_NE(energy, law.getStrainEnergyDensity());
      auto fresh = material<Solid::HolzapfelOgden>();
      fresh.setConstitutivePoint(point);
      EXPECT_DOUBLE_EQ(fresh.getStrainEnergyDensity(), law.getStrainEnergyDensity());
    }

    TEST(ConstitutiveBinding, ActiveInputsPresenceAndFailedBindings)
    {
      using Law = Solid::ActiveContraction<Solid::NeoHookean>;
      static_assert(!PublicCache<Law>);
      static_assert(!ConstBinding<Law>);
      Solid::ActiveFiberLaw::Parameters parameters;
      parameters.stiffness = 2;
      parameters.damping = 0.3;
      parameters.contractility = 0.8;
      parameters.crossBridgeStiffness = 0.5;
      Law law(material<Solid::NeoHookean>(), Solid::ActiveFiberLaw(parameters));
      const auto prototype = law;
      Solid::KinematicState state(3);
      state.setDisplacementGradient(gradient(3, 0.1));
      Solid::ConstitutivePoint point(state);
      law.setConstitutivePoint(point);
      EXPECT_FALSE(law.isDynamic());
      point.set<Solid::Tags::TimeStep>(0.01);
      EXPECT_THROW(law.setConstitutivePoint(point), Alert::Exception);
      EXPECT_THROW(law.getStrainEnergyDensity(), Alert::Exception);
      EXPECT_THROW(law.getActiveResidual(), Alert::Exception);
      point.set<Solid::Tags::PreviousActiveExtension>(-0.02);
      point.set<Solid::Tags::PreviousActiveGamma>(1.2);
      point.set<Solid::Tags::PreviousActiveBeta>(0.7);
      point.set<Solid::Tags::ElectricalActivation>(0.9);
      point.set<Solid::Tags::ActiveExtension>(-0.01);
      point.set<Solid::Tags::ActiveGamma>(1.3);
      point.set<Solid::Tags::ActiveBeta>(0.8);
      auto compare = [&] {
        law.setConstitutivePoint(point);
        auto fresh = prototype;
        fresh.setConstitutivePoint(point);
        EXPECT_TRUE(law.isDynamic());
        EXPECT_DOUBLE_EQ(law.getStrainEnergyDensity(), fresh.getStrainEnergyDensity());
        Math::SpatialMatrix<Real> actual, expected;
        law.getFirstPiolaKirchhoffStress(actual);
        fresh.getFirstPiolaKirchhoffStress(expected);
        EXPECT_LT((actual.getData() - expected.getData()).norm(), 1e-14);
        const auto& a = law.getActiveResidual();
        const auto& b = fresh.getActiveResidual();
        for (size_t i = 0; i < 3; ++i)
        {
          EXPECT_DOUBLE_EQ(a.value[i], b.value[i]);
          for (size_t j = 0; j < 3; ++j)
          {
            EXPECT_DOUBLE_EQ(a.derivative(i, j), b.derivative(i, j));
          }
        }
      };
      compare();
      compare();
      point.set<Solid::Tags::TimeStep>(0.02);
      compare();
      point.set<Solid::Tags::PreviousActiveExtension>(-0.03);
      compare();
      point.set<Solid::Tags::PreviousActiveGamma>(1.1);
      compare();
      point.set<Solid::Tags::PreviousActiveBeta>(0.6);
      compare();
      point.set<Solid::Tags::ElectricalActivation>(1.0);
      compare();
      point.set<Solid::Tags::ActiveExtension>(-0.015);
      compare();
      point.set<Solid::Tags::ActiveGamma>(1.4);
      compare();
      point.set<Solid::Tags::ActiveBeta>(0.85);
      compare();
      point.set<Solid::Tags::PreviousFiberStrain>(0.03);
      compare();
      point.set<Solid::Tags::PreviousFiberStrain>(0.05);
      compare();
      state.setDisplacementGradient(gradient(3, 0.12));
      compare();
      Math::SpatialVector<Real> direction(3);
      direction(0) = 0.6;
      direction(1) = 0.8;
      direction(2) = 0;
      point.set<Solid::Tags::FiberDirection>(direction);
      compare();
      auto retained = law;
      law.setConstitutivePoint(Solid::ConstitutivePoint(state));
      EXPECT_FALSE(law.isDynamic());
      EXPECT_TRUE(retained.isDynamic());
    }
  }
}
