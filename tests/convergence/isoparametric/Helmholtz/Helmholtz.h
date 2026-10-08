/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef RODIN_TESTS_CONVERGENCE_ISOPARAMETRIC_HELMHOLTZ_H
#define RODIN_TESTS_CONVERGENCE_ISOPARAMETRIC_HELMHOLTZ_H

#include <array>
#include <type_traits>
#include <utility>
#include <gtest/gtest.h>

#include "../../CurvedGeometry.h"
#include "../../Helmholtz.h"
#include "../../LiftedErrorNorm.h"
#include "../../LiftedConvergence.h"
#include "../../SineMap.h"

namespace Rodin::Tests::Convergence::Isoparametric::Helmholtz
{
  inline constexpr size_t AssemblyOrder = 11;

  /**
   * @brief Complex Helmholtz acceptance on exact and approximated maps.
   * @par Architecture
   * Each workload owns a freshly mapped mesh and creates fresh field spaces
   * per solve. The continuous data and geometry remain fixed while nominal
   * grid spacing decreases. Physical field errors, map/metric checks,
   * representable patches, an omitted-mass control, and separate quadrature
   * and solver sensitivity checks have independent acceptance conditions.
   * Workload types retain their native or PETSc solver policy.
   */
  template <class Workload>
  class HelmholtzTest : public ::testing::TestWithParam<Geometry::Polytope::Type>
  {
    protected:
      using Map = typename std::remove_cvref_t<
        decltype(std::declval<Workload>().getGeometry())>::Map;
      static constexpr Real L2Margin = 0.55, H1Margin = 0.45;
      // Independent quadrature and algebraic contamination policies.
      static constexpr size_t NormOrder = 13, RefinedOrder = 16, RefinedNormOrder = 18;
      static constexpr Real SolverTolerance = 1e-13, RefinedTolerance = 1e-14;
      static constexpr Real SensitivityTolerance = 1e-6, RoundoffTolerance = 1e-11;
      static constexpr Real PatchTolerance = 1e-9;
      // Nondimensional omitted-mass error floors; not mathematical lower bounds.
      static constexpr Real ControlL2 = 1e-3, ControlH1 = 1e-2, ControlRatio = 2;

      void checkDecomposition(const LiftedErrorNorm::Result& e) const
      {
        for (const auto& norm : {e.field, e.geometry, e.total})
          EXPECT_TRUE(norm.isFinite());
        for (const auto& values :
          {std::array{e.field.getL2(), e.geometry.getL2(), e.total.getL2()},
            std::array{e.field.getH1Seminorm(), e.geometry.getH1Seminorm(),
              e.total.getH1Seminorm()}})
        {
          EXPECT_LE(values[2], values[0] + values[1] + RoundoffTolerance);
          EXPECT_GE(values[2], std::abs(values[0] - values[1]) - RoundoffTolerance);
        }
      }

      // An affine physical field pulls back into the geometry family. Raising
      // field degree to max(2,Q) isolates the nonpolynomial map interpolation.
      void checkMatchedGeometryRates() const
      {
        constexpr size_t Q = Workload::GeometryDegree;
        constexpr size_t K = std::max(size_t(2), Q);
        const auto levels = this->GetParam() == Geometry::Polytope::Type::Segment
          ? std::array<size_t, 3>{5, 9, 17}
          : std::array<size_t, 3>{3, 5, 9};
        std::array<ErrorHistory, 2> histories;
        for (size_t n : levels)
        {
          SCOPED_TRACE(::testing::Message() << "geometry degree=" << Q << " n=" << n);
          Workload problem(this->GetParam(), n, Map::Sine, true);
          LiftedErrorNorm::Result e;
          const auto represented = problem.template solve<K>(HelmholtzData::Field::Affine,
            false, AssemblyOrder, SolverTolerance, NormOrder, &e);
          checkDecomposition(e);
          for (const auto& field : {represented, e.field})
          {
            EXPECT_LT(field.getL2(), PatchTolerance);
            EXPECT_LT(field.getH1Seminorm(), PatchTolerance);
          }
          EXPECT_NEAR(e.total.getL2(), e.geometry.getL2(), PatchTolerance);
          EXPECT_NEAR(
            e.total.getH1Seminorm(), e.geometry.getH1Seminorm(), PatchTolerance);
          const std::array norms{e.geometry, e.total};
          for (size_t component = 0; component < norms.size(); ++component)
          {
            ASSERT_TRUE(norms[component].isFinite());
            ASSERT_GT(norms[component].getL2(), 0);
            ASSERT_GT(norms[component].getH1Seminorm(), 0);
            histories[component].append(Real(1) / Real(n - 1), norms[component]);
          }
        }
        for (const auto& history : histories)
          for (size_t i = 1; i < history.getSize(); ++i)
          {
            const auto& coarse = history.getSample(i - 1).error;
            const auto& fine = history.getSample(i).error;
            const auto rate = history.getAlgebraicRates(i);
            SCOPED_TRACE(::testing::Message()
              << "interval=" << i << " rates=" << rate.getL2() << ","
              << rate.getH1Seminorm());
            EXPECT_GT(coarse.getL2(), fine.getL2());
            EXPECT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
            EXPECT_GT(rate.getL2(), Q + 1 - L2Margin);
            EXPECT_LT(rate.getL2(), Q + 1 + L2Margin);
            EXPECT_GT(rate.getH1Seminorm(), Q - H1Margin);
            EXPECT_LT(rate.getH1Seminorm(), Q + H1Margin);
          }
      }

      void checkMatchedGeometrySensitivity() const
      {
        constexpr size_t K = std::max(size_t(2), Workload::GeometryDegree);
        Workload problem(this->GetParam(), 5, Map::Sine, true);
        std::array<LiftedErrorNorm::Result, 4> errors;
        for (size_t i = 0; i < errors.size(); ++i)
        {
          const auto represented = problem.template solve<K>(HelmholtzData::Field::Affine,
            false, i == 1 ? RefinedOrder : AssemblyOrder,
            i == 2 ? RefinedTolerance : SolverTolerance,
            i == 3 ? RefinedNormOrder : NormOrder, &errors[i]);
          checkDecomposition(errors[i]);
          for (const auto& field : {represented, errors[i].field})
          {
            EXPECT_LT(field.getL2(), PatchTolerance);
            EXPECT_LT(field.getH1Seminorm(), PatchTolerance);
          }
        }
        for (size_t i = 1; i < errors.size(); ++i)
        {
          const std::array base{errors[0].geometry, errors[0].total};
          const std::array refined{errors[i].geometry, errors[i].total};
          for (size_t component = 0; component < base.size(); ++component)
            for (const auto& pair :
              {std::pair{base[component].getL2(), refined[component].getL2()},
                std::pair{
                  base[component].getH1Seminorm(), refined[component].getH1Seminorm()}})
            {
              ASSERT_GT(pair.first, 0);
              ASSERT_TRUE(std::isfinite(pair.second));
              EXPECT_LT(std::abs(pair.second / pair.first - 1), SensitivityTolerance);
            }
        }
      }

      template <size_t K>
      void checkApproximatedRates() const
      {
        std::array<ErrorHistory, 4> histories;
        const auto levels = K == 1 ? std::initializer_list<size_t>{5, 9, 17}
          : this->GetParam() == Geometry::Polytope::Type::Segment
          ? std::initializer_list<size_t>{5, 9, 17, 33}
          : std::initializer_list<size_t>{3, 5, 9};
        for (size_t n : levels)
        {
          SCOPED_TRACE(::testing::Message() << "degree=" << K << " n=" << n);
          Workload problem(this->GetParam(), n, Map::Sine, true);
          LiftedErrorNorm::Result e;
          const auto represented = problem.template solve<K>(HelmholtzData::Field::Smooth,
            false, AssemblyOrder, SolverTolerance, NormOrder, &e);
          checkDecomposition(e);
          const std::array norms{represented, e.field, e.geometry, e.total};
          for (size_t component = 0; component < norms.size(); ++component)
          {
            ASSERT_TRUE(norms[component].isFinite());
            ASSERT_GT(norms[component].getL2(), 0);
            ASSERT_GT(norms[component].getH1Seminorm(), 0);
            histories[component].append(Real(1) / Real(n - 1), norms[component]);
          }
        }
        for (size_t component = 0; component < histories.size(); ++component)
          for (size_t i = 1; i < histories[component].getSize(); ++i)
          {
            const auto& history = histories[component];
            const auto& coarse = history.getSample(i - 1).error;
            const auto& fine = history.getSample(i).error;
            const auto rate = history.getAlgebraicRates(i);
            const size_t degree = component < 2 ? K
              : component == 2                  ? 2
                                                : std::min(K, size_t(2));
            SCOPED_TRACE(::testing::Message()
              << "component=" << component << " interval=" << i << " L2 "
              << coarse.getL2() << " -> " << fine.getL2() << " H1 "
              << coarse.getH1Seminorm() << " -> " << fine.getH1Seminorm()
              << " rates=" << rate.getL2() << "," << rate.getH1Seminorm());
            EXPECT_GT(coarse.getL2(), fine.getL2());
            EXPECT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
            EXPECT_GT(rate.getL2(), degree + 1 - L2Margin);
            EXPECT_LT(rate.getL2(), degree + 1 + L2Margin);
            EXPECT_GT(rate.getH1Seminorm(), degree - H1Margin);
            EXPECT_LT(rate.getH1Seminorm(), degree + H1Margin);
          }
      }

      template <size_t K>
      void checkApproximatedSensitivity() const
      {
        Workload problem(this->GetParam(), 5, Map::Sine, true);
        std::array<LiftedErrorNorm::Result, 4> errors;
        std::array<ErrorNorms, 4> represented{{{0, 0}, {0, 0}, {0, 0}, {0, 0}}};
        for (size_t i = 0; i < errors.size(); ++i)
        {
          represented[i] = problem.template solve<K>(HelmholtzData::Field::Smooth, false,
            i == 1 ? RefinedOrder : AssemblyOrder,
            i == 2 ? RefinedTolerance : SolverTolerance,
            i == 3 ? RefinedNormOrder : NormOrder, &errors[i]);
          checkDecomposition(errors[i]);
        }
        const std::array base{
          represented[0], errors[0].field, errors[0].geometry, errors[0].total};
        for (size_t i = 1; i < errors.size(); ++i)
        {
          const std::array refined{
            represented[i], errors[i].field, errors[i].geometry, errors[i].total};
          for (size_t component = 0; component < base.size(); ++component)
            for (const auto& pair :
              {std::pair{base[component].getL2(), refined[component].getL2()},
                std::pair{
                  base[component].getH1Seminorm(), refined[component].getH1Seminorm()}})
            {
              ASSERT_GT(pair.first, 0);
              ASSERT_TRUE(std::isfinite(pair.second));
              EXPECT_LT(std::abs(pair.second / pair.first - 1), SensitivityTolerance);
            }
        }
      }

      void checkHigherOrderApproximatedRates() const
      {
        static_assert(Workload::GeometryDegree == 2);
        const auto levels = this->GetParam() == Geometry::Polytope::Type::Segment
          ? std::array<size_t, 3>{5, 9, 17}
          : std::array<size_t, 3>{3, 5, 9};
        LiftedConvergence history;
        for (size_t n : levels)
        {
          SCOPED_TRACE(::testing::Message() << "field degree=3 geometry degree=2 n=" << n);
          Workload problem(this->GetParam(), n, Map::Sine, true);
          LiftedErrorNorm::Result lifted;
          const auto represented = problem.template solve<3>(HelmholtzData::Field::Smooth,
            false, AssemblyOrder, SolverTolerance, NormOrder, &lifted);
          ASSERT_FALSE(::testing::Test::HasFatalFailure());
          history.append(Real(1) / Real(n - 1), represented, lifted);
        }
        history.expectMixedRates(3, 2);
      }

      template <size_t K = 2>
      void checkApproximatedPatchAndControl() const
      {
        Workload problem(this->GetParam(), 5, Map::Sine, true);
        LiftedErrorNorm::Result base, wrong;
        const auto patch = problem.template solve<K>(HelmholtzData::Field::Affine, false,
          AssemblyOrder, SolverTolerance, NormOrder, &base);
        const auto control = problem.template solve<K>(HelmholtzData::Field::Affine, true,
          AssemblyOrder, SolverTolerance, NormOrder, &wrong);
        checkDecomposition(base);
        checkDecomposition(wrong);
        for (const auto& e : {patch, base.field})
        {
          EXPECT_LT(e.getL2(), PatchTolerance);
          EXPECT_LT(e.getH1Seminorm(), PatchTolerance);
        }
        EXPECT_EQ(base.geometry.getL2(), wrong.geometry.getL2());
        EXPECT_EQ(base.geometry.getH1Seminorm(), wrong.geometry.getH1Seminorm());
        EXPECT_GT(control.getL2(), ControlL2);
        EXPECT_GT(control.getH1Seminorm(), ControlH1);
        EXPECT_GT(wrong.field.getL2(), ControlL2);
        EXPECT_GT(wrong.field.getH1Seminorm(), ControlH1);
        EXPECT_GT(wrong.total.getL2(), ControlRatio * base.total.getL2());
        EXPECT_GT(wrong.total.getH1Seminorm(), ControlRatio * base.total.getH1Seminorm());
      }

      void checkComplexLiftedMetric() const
      {
        // Identity represented domain, independently prescribed exact sine map.
        Workload problem(this->GetParam(), 3, Map::Sine, true, 0);
        LiftedErrorNorm::Result error;
        problem.template solve<2>(HelmholtzData::Field::Affine, false, AssemblyOrder,
          SolverTolerance, RefinedNormOrder, &error);
        constexpr Real Amplitude = 0.1;
        const Real coefficient = std::abs(Complex(1, 0.5));
        const Real slope = Amplitude * Math::Constants::pi();
        const Real l2 = coefficient * Amplitude / std::sqrt(Real(2));
        const Real h1 = coefficient *
          (UniformGrid(this->GetParam()).getDimension() == 1
              ? std::sqrt(Real(1) / std::sqrt(Real(1) - slope * slope) - 1)
              : slope / std::sqrt(Real(2)));
        const auto& reference = problem.getReference();
        for (auto cell = problem.getMesh().getCell(); cell; ++cell)
        {
          const auto original = reference.getCell(cell->getIndex());
          EXPECT_EQ(cell->getGeometry(), original->getGeometry());
          const auto vertices = cell->getVertices(),
                     originalVertices = original->getVertices();
          ASSERT_EQ(vertices.size(), originalVertices.size());
          for (size_t i = 0; i < vertices.size(); ++i)
            EXPECT_EQ(vertices[i], originalVertices[i]);
          if constexpr (requires { reference.getShard(); })
          {
            EXPECT_EQ(reference.getShard().isOwned(
                        reference.getSpaceDimension(), cell->getIndex()),
              problem.getMesh().getShard().isOwned(
                reference.getSpaceDimension(), cell->getIndex()));
          }
        }
        EXPECT_TRUE(error.field.isFinite());
        EXPECT_LT(error.field.getL2(), PatchTolerance);
        EXPECT_LT(error.field.getH1Seminorm(), PatchTolerance);
        for (const auto& e : {error.geometry, error.total})
        {
          EXPECT_TRUE(e.isFinite());
          EXPECT_NEAR(e.getL2(), l2, PatchTolerance);
          EXPECT_NEAR(e.getH1Seminorm(), h1, PatchTolerance);
        }
      }

      template <size_t K>
      void checkRates() const
      {
        const auto levels = K == 1 ? std::initializer_list<size_t>{5, 9, 17}
                                   : std::initializer_list<size_t>{3, 5, 9};
        ErrorHistory history;
        for (size_t n : levels)
        {
          SCOPED_TRACE(::testing::Message() << "degree=" << K << " n=" << n);
          Workload problem(this->GetParam(), n);
          history.append(Real(1) / Real(n - 1),
            problem.template solve<K>(HelmholtzData::Field::Smooth));
        }
        ASSERT_EQ(history.getSize(), 3u);
        for (size_t i = 1; i < history.getSize(); ++i)
        {
          const auto& coarse = history.getSample(i - 1).error;
          const auto& fine = history.getSample(i).error;
          SCOPED_TRACE(::testing::Message()
            << "interval=" << i << " L2 " << coarse.getL2() << " -> " << fine.getL2()
            << " H1 " << coarse.getH1Seminorm() << " -> " << fine.getH1Seminorm());
          ASSERT_TRUE(coarse.isFinite());
          ASSERT_TRUE(fine.isFinite());
          ASSERT_GT(fine.getL2(), 0);
          ASSERT_GT(fine.getH1Seminorm(), 0);
          ASSERT_GT(coarse.getL2(), fine.getL2());
          ASSERT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
          const auto rate = history.getAlgebraicRates(i);
          SCOPED_TRACE(::testing::Message()
            << "rates " << rate.getL2() << ", " << rate.getH1Seminorm());
          EXPECT_GT(rate.getL2(), K == 1 ? 1.65 : 2.45);
          EXPECT_LT(rate.getL2(), K == 1 ? 2.35 : 3.55);
          EXPECT_GT(rate.getH1Seminorm(), K == 1 ? 0.75 : 1.55);
          EXPECT_LT(rate.getH1Seminorm(), K == 1 ? 1.25 : 2.45);
        }
      }

      template <size_t K>
      void checkPatch() const
      {
        Workload problem(this->GetParam(), 3);
        const auto error = problem.template solve<K>(
          K == 1 ? HelmholtzData::Field::Constant : HelmholtzData::Field::Affine);
        EXPECT_LT(error.getL2(), 1e-9);
        EXPECT_LT(error.getH1Seminorm(), 1e-9);
      }

      void checkNegativeControl() const
      {
        Workload problem(this->GetParam(), 5);
        const auto error = problem.template solve<2>(HelmholtzData::Field::Affine, true);
        ASSERT_TRUE(error.isFinite());
        EXPECT_GT(error.getL2(), 1e-3);
        EXPECT_GT(error.getH1Seminorm(), 1e-2);
      }

      template <size_t K>
      void checkSensitivity() const
      {
        Workload problem(this->GetParam(), 5);
        const auto baseline = problem.template solve<K>(HelmholtzData::Field::Smooth);
        const auto quadrature =
          problem.template solve<K>(HelmholtzData::Field::Smooth, false, 16);
        const auto solver = problem.template solve<K>(
          HelmholtzData::Field::Smooth, false, AssemblyOrder, 1e-14);
        ASSERT_TRUE(baseline.isFinite());
        ASSERT_GT(baseline.getL2(), 0);
        ASSERT_GT(baseline.getH1Seminorm(), 0);
        for (const auto& refined : {quadrature, solver})
        {
          ASSERT_TRUE(refined.isFinite());
          EXPECT_LT(std::abs(refined.getL2() / baseline.getL2() - 1), 1e-6);
          EXPECT_LT(
            std::abs(refined.getH1Seminorm() / baseline.getH1Seminorm() - 1), 1e-6);
        }
      }

      void checkGeometry() const
      {
        Workload problem(this->GetParam(), 2);
        const auto& mesh = problem.getMesh();
        const auto& mapping = problem.getGeometry();
        const size_t dimension = UniformGrid(this->GetParam()).getDimension();
        for (size_t d = 1; d <= dimension; ++d)
        {
          for (auto polytope = mesh.getPolytope(d); polytope; ++polytope)
          {
            const Variational::RealH1Element<2> element(polytope->getGeometry());
            EXPECT_EQ(polytope->getTransformation().getOrder(), element.getOrder());
            const auto& qf =
              QF::PolytopeQuadratureFormula::get(4, polytope->getGeometry());
            const Geometry::PolytopeQuadrature quadrature(*polytope, qf);
            for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
            {
              const auto& point = quadrature.getPoint(qp);
              const auto expected = mapping.mapToPhysical(
                mapping.referencePosition(*polytope, point.getReferenceCoordinates()));
              EXPECT_LT((point.getPhysicalCoordinates() - expected).norm(), 1e-11);
              EXPECT_TRUE(std::isfinite(point.getDistortion()));
              EXPECT_GT(point.getDistortion(), 0);
            }
          }
        }
        EXPECT_NEAR(mesh.getMeasure(dimension), dimension == 1 ? 1.1 : 1, 1e-12);
      }
  };
}

#endif
