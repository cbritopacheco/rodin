/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Native finite-mesh Stokes pressure stability independent of a solve. */
#include "Rodin/Assembly.h"
#include <limits>
#include "../../Convergence.h"
#include "../../CurvedGeometry.h"
#include "../../MixedStability.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Convergence::H::Stokes
{
  TEST(MixedStabilityTest, SubnormalQRRotationPreservesSpectrum)
  {
    // The pressure constant removes the third coordinate. The remaining
    // coupling is [[tiny,tiny],[1,0]], with a unit nonzero squared singular
    // value at working precision. Tiny leading entries must not turn a
    // Givens rotation into a scaling of the order-one trailing column.
    for (Real tiny : {std::numeric_limits<Real>::denorm_min(),
           Real(29) * std::numeric_limits<Real>::denorm_min(),
           std::numeric_limits<Real>::min()})
    {
      SCOPED_TRACE(::testing::Message() << "tiny=" << tiny);
      Math::SparseMatrix<Real> A(2, 2), B(3, 2), M(3, 3);
      A.setIdentity();
      M.setIdentity();
      B.insert(0, 0) = tiny;
      B.insert(0, 1) = tiny;
      B.insert(1, 0) = 1;
      MixedStability::Result result;
      MixedStability::compute(result, A, B, M, {}, Math::Vector<Real>::Unit(3, 2));
      ASSERT_FALSE(::testing::Test::HasFatalFailure());
      MixedStability::expectConsistent(result);
      ASSERT_EQ(result.freeVelocity, 2);
      ASSERT_EQ(result.zeroMeanPressure, 2);
      EXPECT_NEAR(result.eigenvalues(0), 0, MixedStability::ConsistencyTolerance);
      EXPECT_NEAR(result.eigenvalues(1), 1, MixedStability::ConsistencyTolerance);
    }
  }

  TEST(MixedStabilityTest, LogicalConstraintAndIndependentPressureMetric)
  {
    Math::SparseMatrix<Real> A(2, 2), B(2, 2), M(2, 2);
    A.insert(0, 0) = 4;
    A.insert(1, 1) = 9;
    B.insert(0, 0) = 2;
    B.insert(1, 0) = -2;
    B.insert(0, 1) = 7;
    B.insert(1, 1) = -7;
    M.insert(0, 0) = 1;
    M.insert(1, 1) = 3;
    const IndexMap<Real> constrained{{1, 0}};
    const Math::Vector<Real> constant = Math::Vector<Real>::Ones(2);
    MixedStability::Result result;
    MixedStability::compute(result, A, B, M, constrained, constant);
    ASSERT_FALSE(::testing::Test::HasFatalFailure());
    MixedStability::expectConsistent(result);
    ASSERT_EQ(result.freeVelocity, 1);
    ASSERT_EQ(result.zeroMeanPressure, 1);
    // mean=(1,3), T=(-3,1)^T, B_free=(2,-2)^T and A_free=4:
    // Schur=16, reduced mass=12, so beta^2=4/3 independently.
    EXPECT_NEAR(result.eigenvalues(0), Real(4) / 3, MixedStability::ConsistencyTolerance);
    EXPECT_TRUE(MixedStability::hasResolvedPositiveSpectrum(result));
  }

  TEST(MixedStabilityTest, NonNodalConstantKnownSpectrum)
  {
    Math::SparseMatrix<Real> A(2, 2), B(2, 2), M(2, 2);
    A.insert(0, 0) = 4;
    A.insert(1, 1) = 9;
    B.insert(0, 0) = 1;
    B.insert(1, 0) = 2;
    M.insert(0, 0) = 1;
    M.insert(1, 1) = 3;
    Math::Vector<Real> constant(2);
    constant << 2, -1;
    MixedStability::Result result;
    MixedStability::compute(result, A, B, M, {{1, 0}}, constant);
    ASSERT_FALSE(::testing::Test::HasFatalFailure());
    MixedStability::expectConsistent(result);
    // Mean=(2,-3), T=(3,2)^T: Schur=49/4 and mass=21.
    EXPECT_NEAR(result.eigenvalues(0), Real(7) / 12, MixedStability::ConsistencyTolerance);
  }

  TEST(MixedStabilityTest, RectangularObstructionAndZeroDivergence)
  {
    Math::SparseMatrix<Real> A(2, 2), B(3, 2), M(3, 3);
    A.setIdentity();
    M.setIdentity();
    B.insert(0, 0) = 1;
    B.insert(1, 0) = -1;
    for (bool omitted : {false, true})
    {
      if (omitted)
        B.setZero();
      MixedStability::Result result;
      MixedStability::compute(result, A, B, M, {{1, 0}}, Math::Vector<Real>::Ones(3));
      ASSERT_FALSE(::testing::Test::HasFatalFailure());
      MixedStability::expectConsistent(result);
      ASSERT_EQ(result.freeVelocity, 1);
      ASSERT_EQ(result.zeroMeanPressure, 2);
      EXPECT_TRUE(result.isDimensionObstructed());
      EXPECT_FALSE(MixedStability::hasResolvedPositiveSpectrum(result));
      EXPECT_NEAR(result.eigenvalues(0), 0, MixedStability::ConsistencyTolerance);
      EXPECT_NEAR(result.eigenvalues(1), omitted ? 0 : 2,
        MixedStability::ConsistencyTolerance);
      if (omitted)
      {
        EXPECT_TRUE(result.eigenvalues.isZero(0));
      }
    }
  }

  TEST(MixedStabilityWorkspaceTest, SparseCouplingExceedsDenseWorkspace)
  {
    // The old dense 19999-by-650 coupling exceeds 96 MiB. Sparse identity
    // blocks give Schur=T^T T=M0, so every eigenvalue is independently one.
    constexpr Eigen::Index VelocitySize = 20000, PressureSize = 650;
    ASSERT_GT(size_t(VelocitySize - 1) * size_t(PressureSize) * sizeof(Real),
      MixedStability::WorkspaceBytes);
    Math::SparseMatrix<Real> A(VelocitySize, VelocitySize);
    Math::SparseMatrix<Real> B(PressureSize, VelocitySize), M(PressureSize, PressureSize);
    A.setIdentity();
    M.setIdentity();
    for (Eigen::Index i = 0; i < PressureSize; ++i)
      B.insert(i, i) = 1;
    MixedStability::Result result;
    MixedStability::compute(result, A, B, M, {{Index(VelocitySize - 1), 0}},
      Math::Vector<Real>::Ones(PressureSize));
    ASSERT_FALSE(::testing::Test::HasFatalFailure());
    MixedStability::expectConsistent(result);
    ASSERT_EQ(result.freeVelocity, VelocitySize - 1);
    ASSERT_EQ(result.zeroMeanPressure, PressureSize - 1);
    EXPECT_TRUE(MixedStability::hasResolvedPositiveSpectrum(result));
    EXPECT_LT((result.eigenvalues.array() - 1).abs().maxCoeff(),
      MixedStability::ConsistencyTolerance);
  }

  class StokesStabilityTest : public ::testing::TestWithParam<Polytope::Type>
  {
    protected:
      static constexpr size_t AssemblyOrder = 16, RefinedOrder = 20;
      // Dimensionless relative quadrature-contamination policy for a positive
      // smallest eigenvalue, not an inferred uniform stability constant.
      static constexpr Real QuadratureTolerance = 1e-8;

      template <size_t K = 2>
      void measure(MixedStability::Result& result, size_t n, bool curved,
        size_t order, bool omitDivergence = false) const
      {
        auto mesh = UniformGrid(GetParam()).makeMesh(n);
        if (curved)
        {
          CurvedGeometry mapping(mesh);
          mapping.template install<2>();
        }
        H1 velocitySpace(std::integral_constant<size_t, K>{}, mesh, mesh.getDimension());
        H1 pressureSpace(std::integral_constant<size_t, K - 1>{}, mesh);
        TrialFunction u(velocitySpace);
        TrialFunction p(pressureSpace);
        TestFunction v(velocitySpace);
        TestFunction q(pressureSpace);
        auto diffusion = Integral(Jacobian(u), Jacobian(v));
        auto divergence = Integral(Real(omitDivergence ? 0 : 1) * Div(u), q);
        auto pressureMass = Integral(p, q);
        diffusion.setOrder(order);
        divergence.setOrder(order);
        pressureMass.setOrder(order);
        BilinearForm a(u, v);
        BilinearForm b(u, q);
        BilinearForm mass(p, q);
        a = diffusion;
        b = divergence;
        mass = pressureMass;
        a.assemble();
        b.assemble();
        mass.assemble();
        auto boundary = DirichletBC(u, Zero());
        boundary.assemble();
        GridFunction constant(pressureSpace);
        constant = RealFunction(Real(1));
        MixedStability::compute(result, a.getOperator(), b.getOperator(),
          mass.getOperator(), std::get<IndexMap<Real>>(boundary.getDOFs()),
          constant.getData());
      }

      template <size_t K>
      void checkHigherOrderHierarchy() const
      {
        // Three resolved levels are independent of the degree-two coarse-grid
        // rank obstruction. The pressure constant is interpolated in its actual
        // basis; its higher-order coefficients are not assumed to be ones.
        for (bool curved : {false, true})
        {
          for (size_t n : {3u, 4u, 5u})
          {
            SCOPED_TRACE(::testing::Message() << "velocity degree=" << K
              << " curved=" << curved << " n=" << n);
            std::array<MixedStability::Result, 2> measurements;
            for (size_t i = 0; i < measurements.size(); ++i)
            {
              measure<K>(measurements[i], n, curved, i == 0 ? AssemblyOrder : RefinedOrder);
              ASSERT_FALSE(::testing::Test::HasFatalFailure());
              SCOPED_TRACE(::testing::Message() << "freeVelocity=" << measurements[i].freeVelocity
                << " zeroMeanPressure=" << measurements[i].zeroMeanPressure
                << " minimumEigenvalue=" << measurements[i].eigenvalues.minCoeff());
              MixedStability::expectConsistent(measurements[i]);
              EXPECT_FALSE(measurements[i].isDimensionObstructed());
              EXPECT_TRUE(MixedStability::hasResolvedPositiveSpectrum(measurements[i]));
            }
            ASSERT_EQ(measurements[0].freeVelocity, measurements[1].freeVelocity);
            ASSERT_EQ(measurements[0].zeroMeanPressure, measurements[1].zeroMeanPressure);
            ASSERT_GT(measurements[0].eigenvalues.minCoeff(), 0);
            EXPECT_LT(std::abs(measurements[1].eigenvalues.minCoeff() /
              measurements[0].eigenvalues.minCoeff() - 1), QuadratureTolerance);
          }
        }
      }

      template <size_t K>
      void checkHigherOrderMissingDivergence() const
      {
        MixedStability::Result correct, wrong;
        measure<K>(correct, 3, true, AssemblyOrder);
        ASSERT_FALSE(::testing::Test::HasFatalFailure());
        measure<K>(wrong, 3, true, AssemblyOrder, true);
        ASSERT_FALSE(::testing::Test::HasFatalFailure());
        MixedStability::expectConsistent(correct);
        MixedStability::expectConsistent(wrong);
        EXPECT_TRUE(MixedStability::hasResolvedPositiveSpectrum(correct));
        EXPECT_FALSE(MixedStability::hasResolvedPositiveSpectrum(wrong));
        EXPECT_EQ(correct.freeVelocity, wrong.freeVelocity);
        EXPECT_EQ(correct.zeroMeanPressure, wrong.zeroMeanPressure);
        EXPECT_TRUE(wrong.eigenvalues.isZero(0));
      }
  };

  TEST_P(StokesStabilityTest, PressureSpectrumAcrossRefinementLevels)
  {
    for (bool curved : {false, true})
    {
      for (size_t n : {2u, 3u, 5u})
      {
        SCOPED_TRACE(::testing::Message() << "curved=" << curved << " n=" << n);
        std::array<MixedStability::Result, 2> measurements;
        for (size_t i = 0; i < measurements.size(); ++i)
        {
          measure(measurements[i], n, curved, i == 0 ? AssemblyOrder : RefinedOrder);
          ASSERT_FALSE(::testing::Test::HasFatalFailure());
          MixedStability::expectConsistent(measurements[i]);
          ASSERT_FALSE(::testing::Test::HasFatalFailure());
          SCOPED_TRACE(::testing::Message()
            << "freeVelocity=" << measurements[i].freeVelocity
            << "zeroMeanPressure=" << measurements[i].zeroMeanPressure
            << "minimumEigenvalue=" << measurements[i].eigenvalues.minCoeff());
          const bool obstructed = n == 2 && GetParam() != Polytope::Type::Pyramid;
          EXPECT_EQ(measurements[i].isDimensionObstructed(), obstructed);
          EXPECT_EQ(MixedStability::hasResolvedPositiveSpectrum(measurements[i]),
            !obstructed);
        }
        ASSERT_EQ(measurements[0].freeVelocity, measurements[1].freeVelocity);
        ASSERT_EQ(measurements[0].zeroMeanPressure, measurements[1].zeroMeanPressure);
        if (!measurements[0].isDimensionObstructed())
        {
          ASSERT_GT(measurements[0].eigenvalues.minCoeff(), 0);
          EXPECT_LT(std::abs(measurements[1].eigenvalues.minCoeff() /
              measurements[0].eigenvalues.minCoeff() - 1), QuadratureTolerance);
        }
      }
    }
  }

  TEST_P(StokesStabilityTest, MissingDivergenceRejected)
  {
    MixedStability::Result correct, wrong;
    measure(correct, 3, true, AssemblyOrder);
    ASSERT_FALSE(::testing::Test::HasFatalFailure());
    measure(wrong, 3, true, AssemblyOrder, true);
    ASSERT_FALSE(::testing::Test::HasFatalFailure());
    MixedStability::expectConsistent(correct);
    MixedStability::expectConsistent(wrong);
    ASSERT_FALSE(::testing::Test::HasFatalFailure());
    EXPECT_TRUE(MixedStability::hasResolvedPositiveSpectrum(correct));
    EXPECT_FALSE(MixedStability::hasResolvedPositiveSpectrum(wrong));
    EXPECT_EQ(correct.freeVelocity, wrong.freeVelocity);
    EXPECT_EQ(correct.zeroMeanPressure, wrong.zeroMeanPressure);
    // An explicitly zero divergence matrix has exactly zero pressure spectrum;
    // no numerical rank threshold is used to establish this negative control.
    EXPECT_TRUE(wrong.eigenvalues.isZero(0));
  }

  TEST_P(StokesStabilityTest, P3P2PressureSpectrumAcrossRefinementLevels)
  {
    checkHigherOrderHierarchy<3>();
  }

  TEST_P(StokesStabilityTest, P3P2MissingDivergenceRejected)
  {
    checkHigherOrderMissingDivergence<3>();
  }

  TEST_P(StokesStabilityTest, P4P3PressureSpectrumAcrossRefinementLevels)
  {
    checkHigherOrderHierarchy<4>();
  }

  TEST_P(StokesStabilityTest, P4P3MissingDivergenceRejected)
  {
    checkHigherOrderMissingDivergence<4>();
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, StokesStabilityTest,
    ::testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
}
