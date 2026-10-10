/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Tests for convergence-study infrastructure shared by all strategies.
 */

#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>

#include <gtest/gtest.h>
#include <gtest/gtest-spi.h>

#include "Convergence.h"
#include "LiftedConvergence.h"
#include "FieldConvergence.h"
#include "Conductivity.h"
#include "LinearElasticity.h"
#include "Stokes.h"
#include "NonlinearPoisson.h"

namespace Rodin::Tests::Convergence
{
  /** @brief Independent physical quadratic semilinear data on all cell families. */
  TEST(ConvergenceUtilities, QuadraticSemilinearDataOnAllGeometries)
  {
    for (const auto geometry : {Geometry::Polytope::Type::Segment,
           Geometry::Polytope::Type::Triangle, Geometry::Polytope::Type::Quadrilateral,
           Geometry::Polytope::Type::Tetrahedron, Geometry::Polytope::Type::Pyramid,
           Geometry::Polytope::Type::Hexahedron, Geometry::Polytope::Type::Wedge})
    {
      SCOPED_TRACE(UniformGrid::getGeometryName(geometry));
      const auto mesh = UniformGrid(geometry).makeMesh(2);
      const Real amplitude = Real(0.25);
      const NonlinearPoissonData data(
        mesh.getDimension(), amplitude, NonlinearPoissonData::Field::Quadratic);
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        const Geometry::Point point(
          *cell, Geometry::Polytope::Traits(geometry).getCentroid());
        const auto& x = point.getPhysicalCoordinates();
        Real polynomial = 1;
        for (size_t j = 0; j < mesh.getDimension(); ++j)
          polynomial += x(j) * x(j);
        const Real exact = amplitude * polynomial;
        EXPECT_EQ(data.getSolution(x), exact);
        EXPECT_EQ(data.getSolution()(point), exact);
        const auto derivative = data.getGradient(x);
        const auto evaluated = data.getGradient()(point);
        for (size_t j = 0; j < mesh.getDimension(); ++j)
        {
          EXPECT_EQ(derivative(j), 2 * amplitude * x(j));
          EXPECT_EQ(evaluated(j), derivative(j));
        }
        EXPECT_EQ(data.getSource()(point),
          -2 * amplitude * Real(mesh.getDimension()) + exact + exact * exact * exact);
      }
    }
  }

  /** @brief A pressure offset changes traction data, not the Stokes source.
   * @f$p_c=p_0+c@f$ implies @f$\nabla p_c=\nabla p_0@f$.
   * Exact component comparisons require no entity-matching tolerance.
   */
  TEST(ConvergenceUtilities, StokesPressureOffsetPreservesVelocityAndSource)
  {
    constexpr Real PressureOffset = 2;
    using Field = StokesData::Field;
    using Type = Geometry::Polytope::Type;
    for (const auto geometry : {Type::Triangle, Type::Quadrilateral, Type::Tetrahedron,
           Type::Pyramid, Type::Hexahedron, Type::Wedge})
    {
      SCOPED_TRACE(UniformGrid::getGeometryName(geometry));
      const auto mesh = UniformGrid(geometry).makeMesh(2);
      const size_t dim = mesh.getDimension();
      for (const auto field :
        {Field::Affine, Field::Quadratic, Field::Cubic, Field::Quartic, Field::Smooth})
      {
        for (size_t axis = 1; axis < dim; ++axis)
        {
          const StokesData original(dim, field, axis);
          const StokesData shifted(dim, field, axis, PressureOffset);
          for (auto cell = mesh.getCell(); cell; ++cell)
          {
            const auto rc = Geometry::Polytope::Traits(geometry).getCentroid();
            const Geometry::Point point(*cell, rc);
            const auto& x = point.getPhysicalCoordinates();
            EXPECT_EQ(shifted.getPressure(x), original.getPressure(x) + PressureOffset);
            EXPECT_EQ(shifted.getPressure()(point), shifted.getPressure(x));
            const auto velocity = original.getVelocity(x);
            const auto shiftedVelocity = shifted.getVelocity(x);
            const auto force = original.getForcing()(point);
            const auto shiftedForce = shifted.getForcing()(point);
            const auto gradient = original.getPressureGradient(x);
            const auto shiftedGradient = shifted.getPressureGradient(x);
            const auto jacobian = original.getVelocityJacobian(x);
            const auto shiftedJacobian = shifted.getVelocityJacobian(x);
            const auto stress = original.getStress(x);
            const auto shiftedStress = shifted.getStress(x);
            for (size_t i = 0; i < dim; ++i)
            {
              EXPECT_EQ(shiftedVelocity(i), velocity(i));
              EXPECT_EQ(shiftedForce(i), force(i));
              EXPECT_EQ(shiftedGradient(i), gradient(i));
              for (size_t j = 0; j < dim; ++j)
              {
                EXPECT_EQ(shiftedJacobian(i, j), jacobian(i, j));
                EXPECT_EQ(stress(i, j),
                  jacobian(i, j) + jacobian(j, i) -
                    (i == j ? original.getPressure(x) : Real(0)));
                EXPECT_EQ(
                  shiftedStress(i, j), i == j ? -shifted.getPressure(x) : stress(i, j));
              }
            }
          }
        }
      }
    }
  }

  /** The analytic divergence-free shear field does not depend on lambda. */
  TEST(ConvergenceUtilities, DivergenceFreeElasticityHasKnownStressAndSource)
  {
    using Data = LinearElasticity::ManufacturedSolution;
    for (const auto geometry :
      {Geometry::Polytope::Type::Triangle, Geometry::Polytope::Type::Quadrilateral,
        Geometry::Polytope::Type::Tetrahedron, Geometry::Polytope::Type::Pyramid,
        Geometry::Polytope::Type::Hexahedron, Geometry::Polytope::Type::Wedge})
    {
      const auto mesh = UniformGrid(geometry).makeMesh(2);
      const size_t dim = mesh.getDimension();
      const Data soft(dim, 1.5, 1, Data::Field::DivergenceFree);
      const Data stiff(dim, 1e4, 1, Data::Field::DivergenceFree);
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        const auto rc = Geometry::Polytope::Traits(geometry).getCentroid();
        const Geometry::Point p(*cell, rc);
        const auto u = soft.exact(p), f = soft.forcing(p);
        const auto gradient = soft.jacobian(p), strain = soft.strain(p),
                   stress = soft.stress(p);
        const Real pi = Math::Constants::pi();
        EXPECT_EQ(u(0), std::sin(pi * p(1)));
        EXPECT_EQ(f(0), pi * pi * u(0));
        for (size_t i = 0; i < dim; ++i)
        {
          EXPECT_EQ(gradient(i, i), 0);
          EXPECT_EQ(u(i), stiff.exact(p)(i));
          EXPECT_EQ(f(i), stiff.forcing(p)(i));
          if (i > 0)
          {
            EXPECT_EQ(u(i), 0);
            EXPECT_EQ(f(i), 0);
          }
          for (size_t j = 0; j < dim; ++j)
          {
            const Real shear = (i == 0 && j == 1) || (i == 1 && j == 0)
              ? pi * std::cos(pi * p(1))
              : Real(0);
            EXPECT_EQ(gradient(i, j), i == 0 && j == 1 ? shear : Real(0));
            EXPECT_EQ(strain(i, j), shear / 2);
            EXPECT_EQ(stress(i, j), shear);
            EXPECT_EQ(stress(i, j), stiff.stress(p)(i, j));
          }
        }
      }
    }
  }

  TEST(ConductivityDataTest, ExponentialPhysicalDataAndSourcesOnAllGeometries)
  {
    // Floating algebra budget, unrelated to entity identification or topology.
    constexpr Real RoundoffBudget = 64 * std::numeric_limits<Real>::epsilon();
    using Type = Geometry::Polytope::Type;
    for (const auto geometry : {Type::Segment, Type::Triangle, Type::Quadrilateral,
           Type::Tetrahedron, Type::Pyramid, Type::Hexahedron, Type::Wedge})
    {
      SCOPED_TRACE(UniformGrid::getGeometryName(geometry));
      const auto mesh = UniformGrid(geometry).makeMesh(3);
      const size_t dim = mesh.getDimension();
      const auto cell = mesh.getPolytope(dim, 0);
      const Math::Vector<Real> rc = Math::Vector<Real>::Constant(dim, 0.25);
      const Geometry::Point point(*cell, rc);
      const auto& x = point.getPhysicalCoordinates();
      Real sum = 0;
      for (size_t j = 0; j < dim; ++j)
        sum += x(j);
      const Real exact = std::exp(sum);
      const ConductivityData data(dim, ConductivityData::Field::Exponential);
      EXPECT_DOUBLE_EQ(data.getSolution(x), exact);
      EXPECT_DOUBLE_EQ(data.getSolution()(point), exact);
      for (size_t j = 0; j < dim; ++j)
      {
        EXPECT_DOUBLE_EQ(data.getGradient(x)(j), exact);
        EXPECT_DOUBLE_EQ(data.getGradient()(point)(j), exact);
      }
      for (bool poisson : {false, true})
      {
        const Real gamma = poisson ? 1 : 1 + sum;
        const Real source = -Real(dim) * (gamma + (poisson ? 0 : 1)) * exact;
        EXPECT_DOUBLE_EQ(data.getCoefficient(poisson)(point), gamma);
        EXPECT_NEAR(data.getSource(poisson)(point), source,
          RoundoffBudget * std::max(Real(1), std::abs(source)));
      }
    }
  }

  TEST(FieldConvergenceTest, AcceptsEveryFieldOnNonuniformRefinementPaths)
  {
    FieldConvergence<2> algebraic, exponential;
    FieldConvergence<1> scalar;
    for (Real h : {0.5, 0.125, 0.0625})
    {
      algebraic.append(h, {ErrorNorms(h * h, h), ErrorNorms(h * h * h, h * h)});
      scalar.append(h, {ErrorNorms(h * h, h)});
    }
    for (Real p : {1, 3, 4})
    {
      exponential.append(p,
        {ErrorNorms(std::exp(-2 * p), std::exp(-p)),
          ErrorNorms(std::exp(-3 * p), std::exp(-2 * p))});
    }
    algebraic.expectAlgebraicFloor({1.9, 0.9});
    exponential.expectExponentialFloor({1.9, 0.9});
    scalar.expectAlgebraicFloor({1.9, 0.9});
    EXPECT_EQ(algebraic.getSize(), 3u);
    EXPECT_EQ(exponential.getSize(), 3u);
    EXPECT_EQ(scalar.getSize(), 3u);
  }

  /**
   * @brief Extreme finite errors retain their finite mathematical rates.
   *
   * The sequence @f$(10^{300},10^{-200},10^{-300})@f$ has finite,
   * positive, decreasing entries, but its first floating-point quotient
   * overflows. Its logarithmic reductions are nevertheless
   * @f$500\log(10)@f$ and @f$100\log(10)@f$. These independently known
   * reductions certify both paired norms and the single-norm history.
   * Reversing the errors tests quotient underflow and finite negative rates;
   * such a history measures error growth, not convergence.
   */
  TEST(ErrorHistoryTest, ExtremeFiniteErrorsRetainKnownRates)
  {
    ErrorHistory algebraic, exponential, increasing;
    const std::array<Real, 3> errors = {1e300, 1e-200, 1e-300};
    const std::array<Real, 3> scales = {0.5, 0.25, 0.125};
    for (size_t i = 0; i < errors.size(); ++i)
    {
      algebraic.append(scales[i], ErrorNorms(errors[i], errors[i]));
      exponential.append(Real(i + 1), ErrorNorms(errors[i], errors[i]));
      increasing.append(scales[i],
        ErrorNorms(errors[errors.size() - 1 - i], errors[errors.size() - 1 - i]));
    }
    for (size_t i = 1; i < errors.size(); ++i)
    {
      const Real reduction = Real(i == 1 ? 500 : 100) * std::log(Real(10));
      const auto hRate = algebraic.getAlgebraicRates(i);
      const auto pRate = exponential.getExponentialRates(i);
      EXPECT_DOUBLE_EQ(hRate.getL2(), reduction / std::log(Real(2)));
      EXPECT_DOUBLE_EQ(hRate.getH1Seminorm(), reduction / std::log(Real(2)));
      EXPECT_DOUBLE_EQ(pRate.getL2(), reduction);
      EXPECT_DOUBLE_EQ(pRate.getH1Seminorm(), reduction);
      const Real growth = Real(i == 1 ? 100 : 500) * std::log(Real(10));
      EXPECT_DOUBLE_EQ(increasing.getAlgebraicRates(i).getL2(),
        -growth / std::log(Real(2)));
      EXPECT_DOUBLE_EQ(increasing.getAlgebraicRates(i).getH1Seminorm(),
        -growth / std::log(Real(2)));
    }
    NormHistory scalar;
    ErrorHistory paired;
    for (Real error : errors)
    {
      scalar.append(error, error);
      paired.append(error, ErrorNorms(error, error));
    }
    for (size_t i = 1; i < errors.size(); ++i)
    {
      // Equal logarithmic reductions in scale and error give rate one.
      EXPECT_EQ(scalar.getAlgebraicRate(i), 1);
      EXPECT_EQ(paired.getAlgebraicRates(i).getL2(), 1);
      EXPECT_EQ(paired.getAlgebraicRates(i).getH1Seminorm(), 1);
    }
  }

  /**
   * @brief The actual binary floating-point endpoints retain finite rates.
   *
   * For a binary scalar with precision @f$t@f$ and exponent limits
   * @f$e_{\min},e_{\max}@f$, its largest finite value and smallest
   * subnormal value are @f$(2-\epsilon)2^{e_{\max}-1}@f$ and
   * @f$2^{e_{\min}-t}@f$. Their logarithmic reduction is evaluated
   * independently from these integer exponents, not by their quotient.
   */
  TEST(ErrorHistoryTest, FloatingPointEndpointsRetainKnownRates)
  {
    using Limits = std::numeric_limits<Real>;
    static_assert(Limits::radix == 2 && Limits::has_denorm == std::denorm_present);
    constexpr Real Largest = Limits::max(), Smallest = Limits::denorm_min();
    const Real reduction =
      Real(Limits::max_exponent - 1 - Limits::min_exponent + Limits::digits)
        * std::log(Real(2)) + std::log(Real(2) - Limits::epsilon());
    ASSERT_TRUE(std::isfinite(reduction));
    ErrorHistory degree, paired;
    NormHistory scalar;
    degree.append(1, ErrorNorms(Largest, Smallest))
      .append(2, ErrorNorms(Smallest, Largest));
    paired.append(Largest, ErrorNorms(Largest, Largest))
      .append(Smallest, ErrorNorms(Smallest, Smallest));
    scalar.append(Largest, Largest).append(Smallest, Smallest);
    EXPECT_DOUBLE_EQ(degree.getExponentialRates(1).getL2(), reduction);
    EXPECT_DOUBLE_EQ(degree.getExponentialRates(1).getH1Seminorm(), -reduction);
    EXPECT_EQ(paired.getAlgebraicRates(1).getL2(), 1);
    EXPECT_EQ(paired.getAlgebraicRates(1).getH1Seminorm(), 1);
    EXPECT_EQ(scalar.getAlgebraicRate(1), 1);
  }

  /**
   * @brief Nearby large samples retain their resolved logarithmic reduction.
   *
   * Subtracting logarithms unconditionally loses this reduction because
   * @f$\log(10^{300})@f$ and the logarithm of its adjacent representable
   * neighbour round to the same value. A representable quotient must retain
   * its ordinary logarithmic evaluation.
   */
  TEST(ErrorHistoryTest, NearbyLargeSamplesRetainResolvedRates)
  {
    constexpr Real Coarse = 1e300;
    const Real fine = std::nextafter(Coarse, Real(0));
    const Real reduction = std::log(Coarse / fine);
    ASSERT_GT(reduction, 0);
    ErrorHistory scales, errors;
    NormHistory scalar;
    scales.append(Coarse, ErrorNorms(2, 2)).append(fine, ErrorNorms(1, 1));
    scalar.append(Coarse, 2).append(fine, 1);
    errors.append(0.5, ErrorNorms(Coarse, Coarse))
      .append(0.25, ErrorNorms(fine, fine));
    EXPECT_DOUBLE_EQ(scales.getAlgebraicRates(1).getL2(), std::log(Real(2)) / reduction);
    EXPECT_DOUBLE_EQ(scales.getAlgebraicRates(1).getH1Seminorm(),
      std::log(Real(2)) / reduction);
    EXPECT_DOUBLE_EQ(scalar.getAlgebraicRate(1), std::log(Real(2)) / reduction);
    EXPECT_DOUBLE_EQ(errors.getAlgebraicRates(1).getL2(), reduction / std::log(Real(2)));
    EXPECT_DOUBLE_EQ(errors.getAlgebraicRates(1).getH1Seminorm(),
      reduction / std::log(Real(2)));
  }

  /** @brief Resolved quotient evaluation retains the existing rate bits. */
  TEST(ErrorHistoryTest, ResolvedQuotientsPreserveExistingEvaluation)
  {
    using Bits = std::array<std::byte, sizeof(Real)>;
    const std::array<Real, 7> errors = {1e-200, 1e-20, 0.125, 1, 8, 1e20, 1e200};
    for (Real coarse : errors)
      for (Real fine : errors)
      {
        const Real ratio = coarse / fine, inverse = fine / coarse;
        // Extreme unrepresentable quotients have their own known-rate oracle.
        if (!std::isfinite(ratio) || ratio <= 0 || !std::isfinite(inverse) || inverse <= 0)
          continue;
        for (Real scale : {Real(0.5), Real(0.125), Real(1) / 6, Real(1) / 32})
        {
          SCOPED_TRACE(::testing::Message() << "coarse=" << coarse
            << " fine=" << fine << " scale=" << scale);
          ErrorHistory paired, exponential;
          NormHistory scalar;
          paired.append(1, ErrorNorms(coarse, fine))
            .append(scale, ErrorNorms(fine, coarse));
          scalar.append(1, coarse).append(scale, fine);
          exponential.append(1, ErrorNorms(coarse, fine))
            .append(4, ErrorNorms(fine, coarse));
          const auto hRate = paired.getAlgebraicRates(1);
          const auto pRate = exponential.getExponentialRates(1);
          const Real denominator = std::log(Real(1) / scale);
          for (const auto& [actual, expected] :
            {std::pair{hRate.getL2(), std::log(ratio) / denominator},
              std::pair{hRate.getH1Seminorm(), std::log(inverse) / denominator},
              std::pair{scalar.getAlgebraicRate(1), std::log(ratio) / denominator},
              std::pair{pRate.getL2(), std::log(ratio) / 3},
              std::pair{pRate.getH1Seminorm(), std::log(inverse) / 3}})
            EXPECT_EQ(std::bit_cast<Bits>(actual), std::bit_cast<Bits>(expected));
        }
      }
  }

  /**
   * @brief A genuinely unrepresentable rate cannot certify convergence.
   *
   * With @f$\Delta p=\mathrm{min}_{\mathrm{normal}}@f$, the rate
   * @f$500\log(10)/\Delta p@f$ exceeds the floating-point range, whereas
   * @f$\log(2)/\Delta p@f$ remains finite. Each norm is isolated in turn.
   */
  TEST(FieldConvergenceTest, RejectsNonfiniteExponentialRates)
  {
    for (bool h1 : {false, true})
    {
      SCOPED_TRACE(::testing::Message() << "h1=" << h1);
      FieldConvergence<1> study;
      const std::array<Real, 3> extreme = {1e300, 1e-200, 1e-300};
      const std::array<Real, 3> ordinary = {1, 0.5, 0.25};
      constexpr Real ParameterStep = std::numeric_limits<Real>::min();
      const std::array<Real, 3> parameters = {0, ParameterStep, 2 * ParameterStep};
      for (size_t i = 0; i < parameters.size(); ++i)
        study.append(parameters[i], {ErrorNorms(
          h1 ? ordinary[i] : extreme[i], h1 ? extreme[i] : ordinary[i])});
      ::testing::TestPartResultArray failures;
      {
        ::testing::ScopedFakeTestPartResultReporter intercept(
          ::testing::ScopedFakeTestPartResultReporter::INTERCEPT_ONLY_CURRENT_THREAD,
          &failures);
        study.expectExponentialFloor({0.5, 0.5});
      }
      EXPECT_EQ(failures.size(), 1);
      if (failures.size() != 1)
        continue;
      EXPECT_TRUE(failures.GetTestPartResult(0).fatally_failed());
      EXPECT_NE(std::string(failures.GetTestPartResult(0).message()).find(
        h1 ? "std::isfinite(rate.getH1Seminorm())" : "std::isfinite(rate.getL2())"),
        std::string::npos);
    }
  }

  TEST(FieldConvergenceTest, RejectsBadFinalIntervalInOtherField)
  {
    FieldConvergence<2> study;
    for (Real h : {0.5, 0.25, 0.125})
    {
      const Real other = h == 0.125 ? Real(0.25) : h;
      study.append(h, {ErrorNorms(h * h, h), ErrorNorms(2 * other * other, 2 * other)});
    }
    ::testing::TestPartResultArray failures;
    {
      ::testing::ScopedFakeTestPartResultReporter intercept(
        ::testing::ScopedFakeTestPartResultReporter::INTERCEPT_ONLY_CURRENT_THREAD,
        &failures);
      study.expectAlgebraicFloor({1.9, 0.9});
    }
    ASSERT_EQ(failures.size(), 4);
    for (int i = 0; i < failures.size(); ++i)
    {
      EXPECT_NE(
        std::string(failures.GetTestPartResult(i).message()).find("field=1 interval=2"),
        std::string::npos);
    }
  }

  TEST(FieldConvergenceTest, RejectsTwoLevelStudy)
  {
    FieldConvergence<1> study;
    study.append(0.5, {ErrorNorms(0.25, 0.5)}).append(0.25, {ErrorNorms(0.0625, 0.25)});
    ::testing::TestPartResultArray failures;
    {
      ::testing::ScopedFakeTestPartResultReporter intercept(
        ::testing::ScopedFakeTestPartResultReporter::INTERCEPT_ONLY_CURRENT_THREAD,
        &failures);
      study.expectAlgebraicFloor({1.9, 0.9});
    }
    ASSERT_EQ(failures.size(), 1);
    EXPECT_TRUE(failures.GetTestPartResult(0).fatally_failed());
  }

  TEST(FieldConvergenceTest, RejectsUnchangedFinalRefinementParameter)
  {
    FieldConvergence<1> algebraic, exponential;
    algebraic.append(0.5, {ErrorNorms(1, 1)})
      .append(0.25, {ErrorNorms(0.5, 0.5)})
      .append(0.25, {ErrorNorms(0.25, 0.25)});
    exponential.append(1, {ErrorNorms(1, 1)})
      .append(2, {ErrorNorms(0.5, 0.5)})
      .append(2, {ErrorNorms(0.25, 0.25)});
    ::testing::TestPartResultArray failures;
    {
      ::testing::ScopedFakeTestPartResultReporter intercept(
        ::testing::ScopedFakeTestPartResultReporter::INTERCEPT_ONLY_CURRENT_THREAD,
        &failures);
      algebraic.expectAlgebraicFloor({0.1, 0.1});
      exponential.expectExponentialFloor({0.1, 0.1});
    }
    ASSERT_EQ(failures.size(), 2);
    for (int i = 0; i < failures.size(); ++i)
    {
      EXPECT_TRUE(failures.GetTestPartResult(i).fatally_failed());
      EXPECT_NE(
        std::string(failures.GetTestPartResult(i).message()).find("field=0 interval=2"),
        std::string::npos);
    }
  }

  /** @brief Independent algebraic histories exercise every lifted rate interval. */
  TEST(LiftedConvergenceTest, AcceptsSeparateFieldAndGeometryOrders)
  {
    LiftedConvergence study;
    for (Real h : {0.25, 0.125, 0.0625})
    {
      const ErrorNorms field(h * h, h), geometry(h * h * h, h * h);
      const ErrorNorms total(field.getL2() + geometry.getL2(),
        field.getH1Seminorm() + geometry.getH1Seminorm());
      study.append(h, field, {field, geometry, total});
    }
    study.expectRates(1, 2);
  }

  TEST(LiftedConvergenceTest, MixedOrdersDoNotRequireGeometryDominance)
  {
    LiftedConvergence study;
    // A large field-error constant keeps the faster field term dominant on
    // these meshes, although the geometry has the slower asymptotic order.
    constexpr Real FieldConstant = 256;
    for (Real h : {0.25, 0.125, 0.0625})
    {
      const ErrorNorms field(
        FieldConstant * std::pow(h, 4), FieldConstant * std::pow(h, 3)),
        geometry(std::pow(h, 3), h * h);
      const ErrorNorms total(field.getL2() + geometry.getL2(),
        field.getH1Seminorm() + geometry.getH1Seminorm());
      study.append(h, field, {field, geometry, total});
    }
    ::testing::TestPartResultArray failures;
    {
      ::testing::ScopedFakeTestPartResultReporter intercept(
        ::testing::ScopedFakeTestPartResultReporter::INTERCEPT_ONLY_CURRENT_THREAD,
        &failures);
      study.expectRates(3, 2);
    }
    ASSERT_GT(failures.size(), 0);
    for (int i = 0; i < failures.size(); ++i)
    {
      EXPECT_NE(std::string(failures.GetTestPartResult(i).message()).find("component=3"),
        std::string::npos);
    }
    study.expectMixedRates(3, 2);
  }

  TEST(LiftedConvergenceTest, MixedOrdersRejectBadFinalFieldInterval)
  {
    LiftedConvergence study;
    for (Real h : {0.25, 0.125, 0.0625})
    {
      const Real scale = h == 0.0625 ? Real(0.125) : h;
      const ErrorNorms field(std::pow(scale, 4), std::pow(scale, 3)),
        geometry(std::pow(h, 3), h * h);
      const ErrorNorms total(field.getL2() + geometry.getL2(),
        field.getH1Seminorm() + geometry.getH1Seminorm());
      study.append(h, field, {field, geometry, total});
    }
    ::testing::TestPartResultArray failures;
    {
      ::testing::ScopedFakeTestPartResultReporter intercept(
        ::testing::ScopedFakeTestPartResultReporter::INTERCEPT_ONLY_CURRENT_THREAD,
        &failures);
      study.expectMixedRates(3, 2);
    }
    ASSERT_GT(failures.size(), 0);
    for (int i = 0; i < failures.size(); ++i)
    {
      EXPECT_NE(std::string(failures.GetTestPartResult(i).message()).find("interval=2"),
        std::string::npos);
    }
  }

  TEST(LiftedConvergenceTest, MixedOrdersRejectIncreasingTotalAfterCancellation)
  {
    LiftedConvergence study;
    // Nearly opposite defects on the first mesh do not certify monotone total
    // convergence, even when both individual defects have their exact orders.
    constexpr Real FieldConstant = 3.9;
    for (Real h : {0.25, 0.125, 0.0625})
    {
      const ErrorNorms field(
        FieldConstant * std::pow(h, 4), FieldConstant * std::pow(h, 3)),
        geometry(std::pow(h, 3), h * h);
      const Real sign = h == 0.25 ? Real(-1) : Real(1);
      const ErrorNorms total(std::abs(field.getL2() + sign * geometry.getL2()),
        std::abs(field.getH1Seminorm() + sign * geometry.getH1Seminorm()));
      study.append(h, field, {field, geometry, total});
    }
    ::testing::TestPartResultArray failures;
    {
      ::testing::ScopedFakeTestPartResultReporter intercept(
        ::testing::ScopedFakeTestPartResultReporter::INTERCEPT_ONLY_CURRENT_THREAD,
        &failures);
      study.expectMixedRates(3, 2);
    }
    ASSERT_EQ(failures.size(), 2);
    for (int i = 0; i < failures.size(); ++i)
    {
      EXPECT_NE(std::string(failures.GetTestPartResult(i).message())
                  .find("mixed total interval=1"),
        std::string::npos);
    }
  }

  TEST(LiftedConvergenceTest, RejectsBadFinalInterval)
  {
    LiftedConvergence study;
    for (Real h : {0.25, 0.125, 0.0625})
    {
      const Real scale = h == 0.0625 ? Real(0.25) : h;
      const ErrorNorms field(scale * scale, scale), geometry(h * h * h, h * h);
      const ErrorNorms total(field.getL2() + geometry.getL2(),
        field.getH1Seminorm() + geometry.getH1Seminorm());
      study.append(h, field, {field, geometry, total});
    }
    ::testing::TestPartResultArray failures;
    {
      ::testing::ScopedFakeTestPartResultReporter intercept(
        ::testing::ScopedFakeTestPartResultReporter::INTERCEPT_ONLY_CURRENT_THREAD,
        &failures);
      study.expectRates(1, 2);
    }
    ASSERT_GT(failures.size(), 0);
    for (int i = 0; i < failures.size(); ++i)
    {
      EXPECT_NE(std::string(failures.GetTestPartResult(i).message()).find("interval=2"),
        std::string::npos);
    }
  }

  TEST(LiftedConvergenceTest, AcceptsRepresentableFieldGeometryOrders)
  {
    constexpr Real PatchTolerance = 1e-9;
    for (size_t q : {1u, 2u, 3u})
    {
      SCOPED_TRACE(::testing::Message() << "geometry degree=" << q);
      LiftedConvergence study;
      for (Real h : {0.25, 0.125, 0.0625})
      {
        const ErrorNorms geometry(std::pow(h, q + 1), std::pow(h, q));
        study.appendRepresentable(
          h, {0, 0}, {{0, 0}, geometry, geometry}, PatchTolerance);
      }
      study.expectGeometryRates(q);
    }
  }

  TEST(LiftedConvergenceTest, RejectsRepresentedFieldAbovePatchBudget)
  {
    constexpr Real PatchTolerance = 1e-9;
    EXPECT_NONFATAL_FAILURE(
      LiftedConvergence::expectRepresentable(
        {2 * PatchTolerance, 0}, {{0, 0}, {0.01, 0.1}, {0.01, 0.1}}, PatchTolerance),
      "error.getL2()");
  }

  TEST(LiftedConvergenceTest, RejectsLiftedFieldAbovePatchBudget)
  {
    constexpr Real PatchTolerance = 1e-9;
    EXPECT_NONFATAL_FAILURE(
      LiftedConvergence::expectRepresentable(
        {0, 0}, {{2 * PatchTolerance, 0}, {0.01, 0.1}, {0.01, 0.1}}, PatchTolerance),
      "error.getL2()");
  }

  TEST(LiftedConvergenceTest, RejectsRepresentableBadFinalInterval)
  {
    constexpr Real PatchTolerance = 1e-9;
    LiftedConvergence study;
    for (Real h : {0.25, 0.125, 0.0625})
    {
      const Real scale = h == 0.0625 ? Real(0.125) : h;
      const ErrorNorms geometry(std::pow(scale, 4), std::pow(scale, 3));
      study.appendRepresentable(h, {0, 0}, {{0, 0}, geometry, geometry}, PatchTolerance);
    }
    ::testing::TestPartResultArray failures;
    {
      ::testing::ScopedFakeTestPartResultReporter intercept(
        ::testing::ScopedFakeTestPartResultReporter::INTERCEPT_ONLY_CURRENT_THREAD,
        &failures);
      study.expectGeometryRates(3);
    }
    ASSERT_GT(failures.size(), 0);
    for (int i = 0; i < failures.size(); ++i)
    {
      EXPECT_NE(std::string(failures.GetTestPartResult(i).message()).find("interval=2"),
        std::string::npos);
    }
  }

  TEST(LiftedConvergenceTest, RejectsTwoLevelRepresentableStudy)
  {
    constexpr Real PatchTolerance = 1e-9;
    LiftedConvergence study;
    for (Real h : {0.25, 0.125})
    {
      const ErrorNorms geometry(std::pow(h, 4), std::pow(h, 3));
      study.appendRepresentable(h, {0, 0}, {{0, 0}, geometry, geometry}, PatchTolerance);
    }
    ::testing::TestPartResultArray failures;
    {
      ::testing::ScopedFakeTestPartResultReporter intercept(
        ::testing::ScopedFakeTestPartResultReporter::INTERCEPT_ONLY_CURRENT_THREAD,
        &failures);
      study.expectGeometryRates(3);
    }
    ASSERT_EQ(failures.size(), 1);
    EXPECT_TRUE(failures.GetTestPartResult(0).fatally_failed());
    EXPECT_NE(
      std::string(failures.GetTestPartResult(0).message()).find("history.getSize()"),
      std::string::npos);
  }

  TEST(LiftedConvergenceTest, AcceptsGeometrySensitivityWithZeroField)
  {
    LiftedConvergence::expectGeometrySensitivity(
      {{0, 0}, {0.01, 0.1}, {0.01, 0.1}}, {{0, 0}, {0.01, 0.1}, {0.01, 0.1}});
  }

  TEST(LiftedConvergenceTest, RejectsGeometrySensitivityContamination)
  {
    constexpr Real Perturbation = 2 * LiftedConvergence::SensitivityTolerance;
    EXPECT_NONFATAL_FAILURE(
      LiftedConvergence::expectGeometrySensitivity(
        {{0, 0}, {1, 1}, {1, 1}}, {{0, 0}, {1 + Perturbation, 1}, {1, 1}}),
      "std::abs");
  }

  TEST(LiftedConvergenceTest, RejectsInvalidNormDecomposition)
  {
    EXPECT_NONFATAL_FAILURE(
      LiftedConvergence::expectDecomposition({{1, 1}, {2, 2}, {4, 3}}), "values[2]");
    EXPECT_NONFATAL_FAILURE(
      LiftedConvergence::expectDecomposition({{1, 1}, {2, 2}, {0.5, 1}}), "values[2]");
  }

  /**
   * @brief Verifies that algebraic rates use the actual scale ratio.
   *
   * Here @f$h@f$ decreases by four rather than two. Errors proportional to
   * @f$h^2@f$ and @f$h@f$ must still report rates two and one respectively.
   */
  TEST(ErrorHistoryTest, ComputesAlgebraicRatesForNonDyadicSpacing)
  {
    ErrorHistory history;
    history.append(0.5, ErrorNorms(0.25, 0.5)).append(0.125, ErrorNorms(0.015625, 0.125));

    const auto rates = history.getAlgebraicRates(1);
    EXPECT_DOUBLE_EQ(rates.getL2(), 2);
    EXPECT_DOUBLE_EQ(rates.getH1Seminorm(), 1);
  }

  /**
   * @brief Verifies exponential rates for non-unit degree spacing.
   *
   * Errors @f$e^{-2p}@f$ and @f$e^{-p}@f$ must report decay constants two
   * and one even when two polynomial degrees separate the samples.
   */
  TEST(ErrorHistoryTest, ComputesExponentialRatesForNonUnitSpacing)
  {
    ErrorHistory history;
    history.append(1, ErrorNorms(std::exp(-2), std::exp(-1)))
      .append(3, ErrorNorms(std::exp(-6), std::exp(-3)));

    const auto rates = history.getExponentialRates(1);
    EXPECT_DOUBLE_EQ(rates.getL2(), 2);
    EXPECT_DOUBLE_EQ(rates.getH1Seminorm(), 1);
  }

  /** @brief Verifies geometry-independent unit-box boundary partitioning. */
  TEST(UnitBoxBoundaryTest, LabelsCoordinateSidesAndRemainder)
  {
    UniformGrid grid(Geometry::Polytope::Type::Triangle);
    auto mesh = grid.makeMesh(4);
    constexpr Geometry::Attribute lower = 11;
    constexpr Geometry::Attribute upper = 12;
    constexpr Geometry::Attribute remainder = 13;
    UnitBoxBoundary::labelCoordinatePartition(mesh, 0, lower, upper, remainder);

    size_t lowerCount = 0;
    size_t upperCount = 0;
    size_t remainderCount = 0;
    for (auto boundary = mesh.getBoundary(); boundary; ++boundary)
    {
      const auto attribute = boundary->getAttribute();
      ASSERT_TRUE(attribute.has_value());
      if (*attribute == lower)
        ++lowerCount;
      else if (*attribute == upper)
        ++upperCount;
      else if (*attribute == remainder)
        ++remainderCount;
      else
        FAIL() << "Unexpected boundary attribute " << *attribute;
    }

    EXPECT_GT(lowerCount, 0);
    EXPECT_GT(upperCount, 0);
    EXPECT_GT(remainderCount, 0);
  }
}
