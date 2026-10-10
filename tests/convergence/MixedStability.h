/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TESTS_CONVERGENCE_MIXED_STABILITY_H
#define RODIN_TESTS_CONVERGENCE_MIXED_STABILITY_H

#include <algorithm>
#include <cmath>
#include <vector>
#include <Eigen/Eigenvalues>
#include <Eigen/Jacobi>
#include <Eigen/SparseCholesky>
#include <Eigen/SVD>
#include <gtest/gtest.h>
#include "Rodin/Math.h"
#include "Rodin/Types.h"

namespace Rodin::Tests::Convergence
{
  /** @brief Independent finite-dimensional mixed-operator stability measurements.
   * @par Contract
   * The caller supplies velocity energy, divergence and pressure mass matrices,
   * logical essential-DOF indices, and coefficients of the pressure constant.
   * Boundary elimination and the zero-mean pressure basis are constructed
   * independently of a saddle-point solve. Generalized Schur eigenvalues are
   * compared with squared singular values from independent LLT whitening.
   * @par Architecture
   * The reduced divergence remains sparse. LDLT solves construct the Schur
   * complement in pressure-column blocks. Independent LLT dual solves supply
   * rows of the whitened divergence to incremental Givens QR; only its
   * pressure-sized triangular factor is retained for the singular-value
   * calculation. No full dense velocity-pressure coupling is allocated.
   * These finite spectra do not prove mesh-uniform inf-sup stability.
   * @par Distribution
   * Inputs are complete matrices on one process. This measurement introduces
   * no collectives; distributed callers must explicitly provide global inputs.
   */
  class MixedStability
  {
    public:
      // Bounds one blocked velocity workspace before its allocation.
      static constexpr size_t WorkspaceBytes = 96 * 1024 * 1024;
      // The n=5 P4/P3 pyramid has 3925 pressure DOFs; its square matrix
      // requires 123245000 bytes. The next binary allowance is 128 MiB.
      // This is an individual pressure-matrix policy, not total process RSS.
      static constexpr size_t PressureWorkspaceBytes = 128 * 1024 * 1024;
      // Scheduling policy for triangular solves, not a rank or error threshold.
      static constexpr size_t BlockSize = 64;
      // Dimensionless algebraic consistency budget, not a uniform beta bound.
      static constexpr Real ConsistencyTolerance = 1e-10;

      struct Result
      {
          Result()
            : freeVelocity(0),
              zeroMeanPressure(0),
              eigenvalues(),
              meanBasisDefect(0),
              eigenResidual(0),
              spectralDifference(0),
              velocityResidual(0)
          {}

          Index freeVelocity, zeroMeanPressure;
          Math::Vector<Real> eigenvalues;
          Real meanBasisDefect, eigenResidual, spectralDifference, velocityResidual;

          bool isDimensionObstructed() const
          {
            return freeVelocity < zeroMeanPressure;
          }
      };

      static void compute(Result& result, const Math::SparseMatrix<Real>& fullA,
        const Math::SparseMatrix<Real>& fullB, const Math::SparseMatrix<Real>& fullM,
        const IndexMap<Real>& constrained, const Math::Vector<Real>& constant)
      {
        ASSERT_EQ(fullA.rows(), fullA.cols());
        ASSERT_EQ(fullB.cols(), fullA.cols());
        ASSERT_EQ(fullM.rows(), fullM.cols());
        ASSERT_EQ(fullB.rows(), fullM.rows());
        ASSERT_EQ(constant.size(), fullM.rows());
        ASSERT_GT(fullM.rows(), 1);
        std::vector<Eigen::Index> reduced(fullA.rows(), -1);
        for (const auto& [index, value] : constrained)
          ASSERT_LT(index, static_cast<Index>(fullA.rows()));
        Eigen::Index nv = 0;
        for (Index i = 0; i < static_cast<Index>(fullA.rows()); ++i)
          if (constrained.find(i) == constrained.end())
            reduced[i] = nv++;
        const Eigen::Index np = fullM.rows();
        ASSERT_GT(nv, 0);
        ASSERT_LE(size_t(np), PressureWorkspaceBytes / sizeof(Real) / size_t(np));
        const Eigen::Index blockSize = static_cast<Eigen::Index>(
          std::min(BlockSize, WorkspaceBytes / sizeof(Real) / size_t(nv)));
        ASSERT_GT(blockSize, 0);
        std::vector<Eigen::Triplet<Real>> entries;
        for (Eigen::Index col = 0; col < fullA.outerSize(); ++col)
          for (Math::SparseMatrix<Real>::InnerIterator it(fullA, col); it; ++it)
            if (reduced[it.row()] >= 0 && reduced[it.col()] >= 0)
              entries.emplace_back(reduced[it.row()], reduced[it.col()], it.value());
        Math::SparseMatrix<Real> A(nv, nv);
        A.setFromTriplets(entries.begin(), entries.end());
        entries.clear();
        for (Eigen::Index col = 0; col < fullB.outerSize(); ++col)
          for (Math::SparseMatrix<Real>::InnerIterator it(fullB, col); it; ++it)
            if (reduced[it.col()] >= 0)
              entries.emplace_back(it.row(), reduced[it.col()], it.value());
        Math::SparseMatrix<Real> B(np, nv);
        B.setFromTriplets(entries.begin(), entries.end());
        const Math::Matrix<Real> M = fullM;
        const Math::Vector<Real> mean = M * constant;
        Eigen::Index pivot = 0;
        mean.cwiseAbs().maxCoeff(&pivot);
        ASSERT_GT(std::abs(mean(pivot)), 0);
        Math::Matrix<Real> T = Math::Matrix<Real>::Zero(np, np - 1);
        Eigen::Index column = 0;
        for (Eigen::Index row = 0; row < np; ++row)
          if (row != pivot)
          {
            T(row, column) = 1;
            T(pivot, column++) = -mean(row) / mean(pivot);
          }
        Eigen::SimplicialLDLT<Math::SparseMatrix<Real>> factor(A);
        ASSERT_EQ(factor.info(), Eigen::Success);
        ASSERT_GT(factor.vectorD().minCoeff(), 0);
        Math::Matrix<Real> S0(np - 1, np - 1);
        Real residualSquared = 0, rhsSquared = 0;
        for (Eigen::Index first = 0; first < np - 1; first += blockSize)
        {
          const Eigen::Index count = std::min(blockSize, np - 1 - first);
          const Math::Matrix<Real> rhs = B.transpose() * T.middleCols(first, count);
          const Math::Matrix<Real> X = factor.solve(rhs);
          ASSERT_EQ(factor.info(), Eigen::Success);
          S0.middleCols(first, count) = T.transpose() * (B * X);
          residualSquared += (A * X - rhs).squaredNorm();
          rhsSquared += rhs.squaredNorm();
        }
        const Math::Matrix<Real> M0 = T.transpose() * M * T;
        Eigen::GeneralizedSelfAdjointEigenSolver<Math::Matrix<Real>> spectrum(S0, M0);
        ASSERT_EQ(spectrum.info(), Eigen::Success);
        // P A P^T = L L^T and M0 = C C^T give the independently whitened
        // divergence C^-1 T^T B P^T L^-T, without forming a Schur complement.
        Eigen::SimplicialLLT<Math::SparseMatrix<Real>> cholesky(A);
        Eigen::LLT<Math::Matrix<Real>> pressureCholesky(M0);
        ASSERT_EQ(cholesky.info(), Eigen::Success);
        ASSERT_EQ(pressureCholesky.info(), Eigen::Success);
        // Z=L^-1 P B^T T. For a block of coordinate columns E,
        // E^T Z=(T^T B P^T L^-T E)^T. Z=Q R then preserves singular
        // values after pressure whitening: Z C^-T=Q (R C^-T).
        Math::Matrix<Real> R = Math::Matrix<Real>::Zero(np - 1, np - 1);
        for (Eigen::Index first = 0; first < nv; first += blockSize)
        {
          const Eigen::Index count = std::min(blockSize, nv - first);
          Math::Matrix<Real> coordinates = Math::Matrix<Real>::Zero(nv, count);
          for (Eigen::Index col = 0; col < count; ++col)
            coordinates(first + col, col) = 1;
          const Math::Matrix<Real> dual = cholesky.matrixU().solve(coordinates);
          Math::Matrix<Real> rows =
            (T.transpose() * (B * (cholesky.permutationPinv() * dual))).transpose();
          for (Eigen::Index row = 0; row < count; ++row)
          {
            for (Eigen::Index column = 0; column < np - 1; ++column)
            {
              // Exact zero skips an identity rotation, not a rank decision.
              if (rows(row, column) == 0)
                continue;
              // Form the rotation from scaled ratios. Dividing by a rounded
              // subnormal hypot can violate c*c+s*s=1 and corrupt trailing
              // entries whose magnitude is unrelated to the tiny pivot.
              Eigen::JacobiRotation<Real> rotation;
              Real radius;
              rotation.makeGivens(R(column, column), rows(row, column), &radius);
              const Real cosine = rotation.c();
              const Real sine = -rotation.s();
              for (Eigen::Index next = column + 1; next < np - 1; ++next)
              {
                const Real upper = R(column, next), lower = rows(row, next);
                R(column, next) = cosine * upper + sine * lower;
                rows(row, next) = -sine * upper + cosine * lower;
              }
              R(column, column) = radius;
              rows(row, column) = 0;
            }
          }
        }
        const Math::Matrix<Real> whitened = pressureCholesky.matrixL().solve(R.transpose());
        // Divide-and-conquer retains the independent singular-value oracle;
        // forming another normal-equation spectrum would square conditioning.
        Eigen::BDCSVD<Math::Matrix<Real>> singularSpectrum(whitened);
        ASSERT_EQ(singularSpectrum.info(), Eigen::Success);
        Math::Vector<Real> squaredSingular = Math::Vector<Real>::Zero(np - 1);
        for (Eigen::Index i = 0; i < std::min(nv, np - 1); ++i)
          squaredSingular(np - 2 - i) =
            singularSpectrum.singularValues()(i) * singularSpectrum.singularValues()(i);
        // Padding follows the rectangular dimensions exactly, not a numerical
        // rank threshold. A coarse pressure dimension may exceed velocity's.
        const Real scale = std::max(spectrum.eigenvalues().cwiseAbs().maxCoeff(),
          squaredSingular.cwiseAbs().maxCoeff());
        Real maximumResidual = 0;
        for (Eigen::Index i = 0; i < spectrum.eigenvalues().size(); ++i)
        {
          const Real lambda = spectrum.eigenvalues()(i);
          const auto vector = spectrum.eigenvectors().col(i);
          const Real denominator = (S0.norm() + std::abs(lambda) * M0.norm()) * vector.norm();
          const Real residual = (S0 * vector - lambda * M0 * vector).norm();
          maximumResidual = std::max(maximumResidual,
            denominator > 0 ? residual / denominator : residual);
        }
        result.freeVelocity = nv;
        result.zeroMeanPressure = np - 1;
        result.eigenvalues = spectrum.eigenvalues();
        result.meanBasisDefect = (mean.transpose() * T).norm() /
          std::max(Real(1), mean.norm() * T.norm());
        result.eigenResidual = maximumResidual;
        const Real difference =
          (result.eigenvalues - squaredSingular).cwiseAbs().maxCoeff();
        result.spectralDifference = scale > 0 ? difference / scale : difference;
        result.velocityResidual =
          std::sqrt(residualSquared) / std::max(Real(1), std::sqrt(rhsSquared));
      }

      static void expectConsistent(const Result& result)
      {
        SCOPED_TRACE(::testing::Message()
          << "mean-basis defect=" << result.meanBasisDefect
          << " eigen-residual=" << result.eigenResidual
          << " spectral difference=" << result.spectralDifference
          << " velocity residual=" << result.velocityResidual);
        ASSERT_EQ(result.eigenvalues.size(), result.zeroMeanPressure);
        EXPECT_TRUE(result.eigenvalues.allFinite());
        for (Real defect : {result.meanBasisDefect, result.eigenResidual,
               result.spectralDifference, result.velocityResidual})
        {
          EXPECT_TRUE(std::isfinite(defect));
          EXPECT_GE(defect, 0);
          EXPECT_LT(defect, ConsistencyTolerance);
        }
      }

      /** @brief Resolve positivity above the algebraic consistency scale.
       * This finite-spectrum predicate is not a numerical ownership/rank test
       * and supplies no mesh-independent lower bound for the inf-sup constant.
       */
      static bool hasResolvedPositiveSpectrum(const Result& result)
      {
        return result.eigenvalues.size() > 0 && result.eigenvalues.allFinite() &&
          result.eigenvalues.minCoeff() >
            ConsistencyTolerance * result.eigenvalues.cwiseAbs().maxCoeff();
      }
  };
}
#endif
