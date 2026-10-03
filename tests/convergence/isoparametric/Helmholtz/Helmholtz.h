/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef RODIN_TESTS_CONVERGENCE_ISOPARAMETRIC_HELMHOLTZ_H
#define RODIN_TESTS_CONVERGENCE_ISOPARAMETRIC_HELMHOLTZ_H

#include <gtest/gtest.h>

#include "../../CurvedGeometry.h"
#include "../../Helmholtz.h"

namespace Rodin::Tests::Convergence::Isoparametric::Helmholtz
{
  inline constexpr size_t AssemblyOrder = 11;

  /**
   * @brief Backend-independent acceptance for complex Helmholtz on exact P2 maps.
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
            const auto& quadrature = polytope->getQuadrature(qf);
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
