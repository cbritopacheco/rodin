#include <gtest/gtest.h>

#include <cmath>

#include "Rodin/Alert/Exception.h"
#include "Rodin/Adaptation/CellDeformation.h"
#include "Rodin/Adaptation/CellGeomCache.h"

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

  /// @brief One-dimensional deformation has scalar stretch and no shape distortion.
  TEST(Rodin_Adaptation_CellDeformation, SegmentStretchAndDerivatives)
  {
    constexpr Real Tolerance = Real(1e-12);
    CellDeformation deformation(1);
    Math::SpatialMatrix<Real> gradient(1, 1), direction(1, 1);
    gradient(0, 0) = Real(0.5);
    direction(0, 0) = Real(0.3);
    deformation.setDisplacementGradient(gradient);
    EXPECT_EQ(deformation.getDimension(), 1u);
    EXPECT_NEAR(deformation.getJacobian(), Real(1.5), Tolerance);
    EXPECT_NEAR(deformation.getInverseTranspose()(0, 0), Real(2) / Real(3), Tolerance);
    EXPECT_NEAR(deformation.getJacobianAction(direction), direction(0, 0), Tolerance);
    EXPECT_NEAR(deformation.getJacobianSecondAction(direction, direction), 0, Tolerance);
    EXPECT_NEAR(deformation.getRelativeDistortion(), 1, Tolerance);
    EXPECT_NEAR(deformation.getRelativeDistortionGradient()(0, 0), 0, Tolerance);
    EXPECT_NEAR(deformation.getRelativeDistortionAction(direction), 0, Tolerance);
    EXPECT_NEAR(
      deformation.getRelativeDistortionSecondAction(direction, direction), 0, Tolerance);

    gradient(0, 0) = Real(-2);
    deformation.setDeformationGradient(gradient);
    EXPECT_FALSE(deformation.isAdmissible());
    EXPECT_TRUE(deformation.isInvertible());
    EXPECT_NEAR(deformation.getInverseTranspose()(0, 0), Real(-0.5), Tolerance);

    gradient(0, 0) = 0;
    deformation.setDeformationGradient(gradient);
    EXPECT_FALSE(deformation.isAdmissible());
    EXPECT_FALSE(deformation.isInvertible());
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

  TEST(Rodin_Adaptation_CellDeformation, NormalizedInversionEnergyAndCurvature)
  {
    for (const size_t d : {2u, 3u})
    {
      Math::SpatialMatrix<Real> F = Math::SpatialMatrix<Real>::Identity(d, d);
      F(0, 0) = Real(1.3);
      F(0, 1) = Real(0.2);
      F(1, 1) = Real(0.8);
      Math::SpatialMatrix<Real> G = Math::SpatialMatrix<Real>::Identity(d, d);
      G(0, 1) = Real(0.3);
      const auto energy = [d](const Math::SpatialMatrix<Real>& matrix) {
        return -std::log(matrix.determinant()) +
          Real(d) / Real(2) * std::log(matrix.squaredNorm() / Real(d));
      };
      CellDeformation deformation(d);
      deformation.setDeformationGradient(F);
      const Real Q = deformation.getRelativeDistortion();
      EXPECT_NEAR(energy(F), Real(d) / Real(2) * std::log(Q), Real(1e-14));
      for (const Real scale : {Real(0.1), Real(1), Real(5)})
        EXPECT_NEAR(energy(Math::SpatialMatrix<Real>(scale * F)), energy(F), Real(1e-14));
      const auto hessian = [&](const Math::SpatialMatrix<Real>& direction) {
        const Real action = deformation.getRelativeDistortionAction(direction);
        return Real(d) / Real(2) *
          (deformation.getRelativeDistortionSecondAction(direction, direction) / Q -
            action * action / (Q * Q));
      };
      EXPECT_NEAR(hessian(F), Real(0), Real(1e-13));
      constexpr Real eps = Real(1e-4);
      const Real difference =
        (energy(Math::SpatialMatrix<Real>(F + eps * G)) - Real(2) * energy(F) +
          energy(Math::SpatialMatrix<Real>(F - eps * G))) /
        (eps * eps);
      EXPECT_NEAR(hessian(G), difference, Real(1e-6));
      F(1, 1) = Real(1e-6);
      EXPECT_GT(energy(F), Real(5));
    }
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
  /// @brief The affine prototype cache agrees with exact skew-triangle geometry in both orientations.
  TEST(Rodin_Adaptation_CellGeomCache, ExactAffineGeometry)
  {
    for (bool reversed : {false, true})
    {
      Geometry::LocalMesh::Builder builder;
      builder.initialize(2).nodes(3).vertex({1, 2}).vertex({3, 3}).vertex({2, 5});
      if (reversed)
        builder.polytope(Geometry::Polytope::Type::Triangle, {0, 2, 1});
      else
        builder.polytope(Geometry::Polytope::Type::Triangle, {0, 1, 2});
      auto mesh = builder.finalize();
      const auto [cache, indices] = precomputeCellGeometry(mesh);
      ASSERT_EQ(cache.size(), 1);
      ASSERT_EQ(indices.size(), 1);
      EXPECT_EQ(indices.at(0), 0);
      const auto& cell = cache[0];
      EXPECT_EQ(cell.index, 0);
      EXPECT_EQ(cell.vertices,
        (reversed ? std::array<Index, 3>{0, 2, 1} : std::array<Index, 3>{0, 1, 2}));
      EXPECT_NEAR(cell.detAK, reversed ? -5.0 : 5.0, 1e-13);
      EXPECT_NEAR(cell.area, 2.5, 1e-13);
      EXPECT_NEAR(cell.Jscale, 5, 1e-13);
      EXPECT_EQ(cell.sigmaK, reversed ? -1 : 1);
      const Real jacobian[2][2] = {{2, 1}, {1, 3}};
      const Real inverse[2][2] = {{0.6, -0.2}, {-0.2, 0.4}};
      const Real metric[2][2] = {{5, 5}, {5, 10}};
      const Real gradients[3][2] = {{-0.4, -0.2}, {0.6, -0.2}, {-0.2, 0.4}};
      for (size_t r = 0; r < 2; ++r)
      {
        for (size_t c = 0; c < 2; ++c)
        {
          EXPECT_NEAR(cell.A(r, c), jacobian[r][reversed ? 1 - c : c], 1e-13);
          EXPECT_NEAR(cell.Ainv(r, c), inverse[reversed ? 1 - r : r][c], 1e-13);
          EXPECT_NEAR(cell.AinvT(r, c), inverse[reversed ? 1 - c : c][r], 1e-13);
          EXPECT_NEAR(cell.C(r, c), metric[r][c], 1e-13);
        }
      }
      for (size_t r = 0; r < 3; ++r)
      {
        const size_t basis = reversed && r > 0 ? 3 - r : r;
        for (size_t c = 0; c < 2; ++c)
          EXPECT_NEAR(cell.gradN(r, c), gradients[basis][c], 1e-13);
      }
    }
  }

  /// @brief Distinct cells retain their mesh indices, vertex keys and independent measures.
  TEST(Rodin_Adaptation_CellGeomCache, MultipleCells)
  {
    auto mesh = Geometry::LocalMesh::Builder()
                  .initialize(2)
                  .nodes(4)
                  .vertex({0, 0})
                  .vertex({2, 0})
                  .vertex({0, 1})
                  .vertex({2, 3})
                  .polytope(Geometry::Polytope::Type::Triangle, {0, 1, 2})
                  .polytope(Geometry::Polytope::Type::Triangle, {1, 3, 2})
                  .finalize();
    const auto [cache, indices] = precomputeCellGeometry(mesh);
    ASSERT_EQ(cache.size(), 2);
    ASSERT_EQ(indices.size(), 2);
    const auto& first = cache.at(indices.at(0));
    const auto& second = cache.at(indices.at(1));
    EXPECT_EQ(first.index, 0);
    EXPECT_EQ(second.index, 1);
    EXPECT_EQ(first.vertices, (std::array<Index, 3>{0, 1, 2}));
    EXPECT_EQ(second.vertices, (std::array<Index, 3>{1, 3, 2}));
    EXPECT_NEAR(first.area, 1, 1e-13);
    EXPECT_NEAR(second.area, 3, 1e-13);
    EXPECT_NEAR(first.Jscale, 2, 1e-13);
    EXPECT_NEAR(second.Jscale, 6, 1e-13);
  }

  /// @brief Singular affine maps fail rather than yielding an inverse with nonfinite entries.
  TEST(Rodin_Adaptation_CellGeomCache, DegenerateTriangle)
  {
    auto mesh = Geometry::LocalMesh::Builder()
                  .initialize(2)
                  .nodes(3)
                  .vertex({0, 0})
                  .vertex({1, 0})
                  .vertex({2, 0})
                  .polytope(Geometry::Polytope::Type::Triangle, {{0, 1, 2}})
                  .finalize();
    EXPECT_THROW(precomputeCellGeometry(mesh), Alert::Exception);
  }

  /// @brief The triangular prototype rejects a quadrilateral rather than caching wrong geometry.
  TEST(Rodin_Adaptation_CellGeomCache, RejectsQuadrilateral)
  {
    auto mesh = Geometry::LocalMesh::Builder()
                  .initialize(2)
                  .nodes(4)
                  .vertex({0, 0})
                  .vertex({1, 0})
                  .vertex({1, 1})
                  .vertex({0, 1})
                  .polytope(Geometry::Polytope::Type::Quadrilateral, {0, 1, 2, 3})
                  .finalize();
    EXPECT_THROW(precomputeCellGeometry(mesh), Alert::Exception);
  }

  /// @brief Topological and embedding dimensions must both match the planar prototype.
  TEST(Rodin_Adaptation_CellGeomCache, RejectsOtherDimensions)
  {
    auto segment = Geometry::LocalMesh::Builder()
                     .initialize(1)
                     .nodes(2)
                     .vertex({0})
                     .vertex({1})
                     .polytope(Geometry::Polytope::Type::Segment, {0, 1})
                     .finalize();
    EXPECT_THROW(precomputeCellGeometry(segment), Alert::Exception);
    auto surface = Geometry::LocalMesh::Builder()
                     .initialize(3)
                     .nodes(3)
                     .vertex({0, 0, 0})
                     .vertex({1, 0, 0})
                     .vertex({0, 1, 0})
                     .polytope(Geometry::Polytope::Type::Triangle, {0, 1, 2})
                     .finalize();
    EXPECT_THROW(precomputeCellGeometry(surface), Alert::Exception);
  }

}
