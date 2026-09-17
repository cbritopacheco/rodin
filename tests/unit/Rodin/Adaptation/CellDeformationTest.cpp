#include <gtest/gtest.h>

#include <cmath>

#include "Rodin/Adaptation/CellDeformation.h"

using namespace Rodin;
using namespace Rodin::Adaptation;

namespace Rodin::Tests::Unit
{
  namespace
  {
    Math::SpatialMatrix<Real> matrix2(Real a, Real b, Real c, Real d)
    {
      Math::SpatialMatrix<Real> m = Math::SpatialMatrix<Real>::Identity(2, 2);
      m(0, 0) = a;
      m(0, 1) = b;
      m(1, 0) = c;
      m(1, 1) = d;
      return m;
    }

    /// @brief Rotation by @p theta.
    Math::SpatialMatrix<Real> rotation2(Real theta)
    {
      return matrix2(std::cos(theta), -std::sin(theta), std::sin(theta), std::cos(theta));
    }
  }

  /// @brief The Jacobian is the determinant of the deformation gradient.
  TEST(Rodin_Adaptation_CellDeformation, Jacobian_IsDeterminant)
  {
    CellDeformation def(2);
    def.setDeformationGradient(matrix2(2, 0, 0, 3));
    EXPECT_NEAR(def.getJacobian(), 6.0, 1e-14);
    EXPECT_TRUE(def.isAdmissible());
    EXPECT_TRUE(def.isInvertible());
  }

  /// @brief A displacement gradient contributes @f$F=I+H@f$.
  TEST(Rodin_Adaptation_CellDeformation, DisplacementGradient_AddsIdentity)
  {
    CellDeformation def(2);
    Math::SpatialMatrix<Real> h = Math::SpatialMatrix<Real>::Identity(2, 2);
    h(0, 0) = 1;
    h(0, 1) = 0;
    h(1, 0) = 0;
    h(1, 1) = 1;
      // H = I gives F = 2I, so j = 4.
    def.setDisplacementGradient(h);
    EXPECT_NEAR(def.getJacobian(), 4.0, 1e-14);
  }

  /**
   * @brief An inverted cell is invertible but not admissible.
   *
   * The distinction is load-bearing: the size-control terms act on inverted
   * cells precisely to pull them back, and need @f$F^{-\top}@f$ there.
   */
  TEST(Rodin_Adaptation_CellDeformation, Inverted_IsInvertibleButNotAdmissible)
  {
    CellDeformation def(2);
    def.setDeformationGradient(matrix2(-1, 0, 0, 1));

    EXPECT_LT(def.getJacobian(), 0.0);
    EXPECT_FALSE(def.isAdmissible());
    EXPECT_TRUE(def.isInvertible());

      // Defined on an inverted cell; must not trip the admissibility assert.
    const auto& FinvT = def.getInverseTranspose();
    EXPECT_NEAR(FinvT(0, 0), -1.0, 1e-14);
  }

  /// @brief A singular deformation gradient is neither invertible nor admissible.
  TEST(Rodin_Adaptation_CellDeformation, Singular_IsNotInvertible)
  {
    CellDeformation def(2);
    def.setDeformationGradient(matrix2(1, 0, 0, 0));
    EXPECT_FALSE(def.isInvertible());
    EXPECT_FALSE(def.isAdmissible());
  }

  /// @brief The relative distortion is minimal, and equal to one, at a similarity.
  TEST(Rodin_Adaptation_CellDeformation, RelativeDistortion_MinimalAtIdentity)
  {
    CellDeformation def(2);
    def.setDeformationGradient(Math::SpatialMatrix<Real>::Identity(2, 2));
    EXPECT_NEAR(def.getRelativeDistortion(), 1.0, 1e-14);

    def.setDeformationGradient(matrix2(4, 0, 0, 1));
    EXPECT_GT(def.getRelativeDistortion(), 1.0);
  }

  /**
   * @brief The relative distortion is invariant under @f$F\mapsto sRF@f$.
   *
   * This is what makes it a measure of *shape* alone: neither a uniform
   * rescaling nor a rotation changes it.
   */
  TEST(Rodin_Adaptation_CellDeformation, RelativeDistortion_InvariantUnderSimilarity)
  {
    const auto F = matrix2(3, 1, 0.5, 2);

    CellDeformation plain(2);
    plain.setDeformationGradient(F);
    const Real reference = plain.getRelativeDistortion();

    for (const Real s : {0.25, 1.0, 7.0})
    {
      for (const Real theta : {0.0, 0.7, 2.9})
      {
        CellDeformation transformed(2);
        transformed.setDeformationGradient(
          Math::SpatialMatrix<Real>(s * (rotation2(theta) * F)));
        EXPECT_NEAR(transformed.getRelativeDistortion(), reference, 1e-12)
          << "scaled by " << s << ", rotated by " << theta;
      }
    }
  }

  /**
   * @brief The Jacobian action is the directional derivative of @f$\det F@f$.
   *
   * Checked against a central difference, which is the definition the
   * admissibility linearisation relies on.
   */
  TEST(Rodin_Adaptation_CellDeformation, JacobianAction_MatchesFiniteDifference)
  {
    const auto F = matrix2(1.3, 0.2, -0.4, 0.9);
    const auto G = matrix2(0.7, -0.3, 0.15, 0.45);
    constexpr Real eps = 1e-6;

    CellDeformation def(2);
    def.setDeformationGradient(F);
    const Real analytic = def.getJacobianAction(G);

    CellDeformation plus(2), minus(2);
    plus.setDeformationGradient(Math::SpatialMatrix<Real>(F + eps * G));
    minus.setDeformationGradient(Math::SpatialMatrix<Real>(F - eps * G));
    const Real numeric = (plus.getJacobian() - minus.getJacobian()) / (2 * eps);

    EXPECT_NEAR(analytic, numeric, 1e-6);
  }

  /// @brief The distortion action is the directional derivative of @f$Q@f$.
  TEST(Rodin_Adaptation_CellDeformation, DistortionAction_MatchesFiniteDifference)
  {
    const auto F = matrix2(1.3, 0.2, -0.4, 0.9);
    const auto G = matrix2(0.7, -0.3, 0.15, 0.45);
    constexpr Real eps = 1e-6;

    CellDeformation def(2);
    def.setDeformationGradient(F);
    const Real analytic = def.getRelativeDistortionAction(G);

    CellDeformation plus(2), minus(2);
    plus.setDeformationGradient(Math::SpatialMatrix<Real>(F + eps * G));
    minus.setDeformationGradient(Math::SpatialMatrix<Real>(F - eps * G));
    const Real numeric =
      (plus.getRelativeDistortion() - minus.getRelativeDistortion()) / (2 * eps);

    EXPECT_NEAR(analytic, numeric, 1e-6);
  }

  /**
   * @brief Setting a new gradient invalidates every cached quantity.
   *
   * The class caches lazily, so a stale cache would silently report the
   * previous cell's state.
   */
  TEST(Rodin_Adaptation_CellDeformation, Cache_InvalidatedOnNewGradient)
  {
    CellDeformation def(2);

    def.setDeformationGradient(matrix2(2, 0, 0, 2));
      // Populate every cache before changing state.
    const Real j0 = def.getJacobian();
    const Real q0 = def.getRelativeDistortion();
    const auto finvT0 = def.getInverseTranspose();
    const auto dQdF0 = def.getRelativeDistortionGradient();
    EXPECT_NEAR(j0, 4.0, 1e-14);
    EXPECT_NEAR(q0, 1.0, 1e-14);

    def.setDeformationGradient(matrix2(5, 0, 0, 1));
    EXPECT_NEAR(def.getJacobian(), 5.0, 1e-14);
    EXPECT_GT(def.getRelativeDistortion(), q0);
    EXPECT_NEAR(def.getInverseTranspose()(0, 0), 0.2, 1e-14);
    EXPECT_GT(std::abs(def.getRelativeDistortionGradient()(0, 0) - dQdF0(0, 0)), 1e-9);
    EXPECT_NEAR(finvT0(0, 0), 0.5, 1e-14);
  }

  /// @brief Repeated reads return the cached value unchanged.
  TEST(Rodin_Adaptation_CellDeformation, Cache_StableAcrossReads)
  {
    CellDeformation def(3);
    Math::SpatialMatrix<Real> F = Math::SpatialMatrix<Real>::Identity(3, 3);
    F(0, 0) = 2;
    F(1, 1) = 3;
    F(2, 2) = 4;
    def.setDeformationGradient(F);

    EXPECT_NEAR(def.getJacobian(), 24.0, 1e-13);
    EXPECT_NEAR(def.getJacobian(), 24.0, 1e-13);
    const Real q1 = def.getRelativeDistortion();
    const Real q2 = def.getRelativeDistortion();
    EXPECT_EQ(q1, q2);
  }
}
