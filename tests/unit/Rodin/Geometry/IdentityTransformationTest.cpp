/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <Rodin/Geometry.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  /// @brief Mutable affine map used to check geometric value semantics.
  class StateTransformation final : public PolytopeTransformation
  {
    public:
      /// @brief Constructs a unit-scale segment map.
      constexpr StateTransformation()
        : PolytopeTransformation(1, 1),
          m_scale(1)
      {}

      /// @brief Preserves constant-evaluated destruction.
      constexpr ~StateTransformation() override = default;

      /**
       * @brief Copies geometry into a fresh state.
       * @param other Source map.
       */
      constexpr StateTransformation(const StateTransformation& other) = default;

      /**
       * @brief Moves geometry into a fresh state.
       * @param other Source map.
       */
      constexpr StateTransformation(StateTransformation&& other) noexcept = default;

      /**
       * @brief Replaces the geometry state by move assignment.
       * @param other Source map.
       * @returns This map.
       */
      StateTransformation& operator=(StateTransformation&& other) noexcept = default;

      /**
       * @brief Changes scale and invalidates previous physical results.
       * @param scale New scale factor.
       */
      void setScale(Real scale)
      {
        m_scale = scale;
      }

      /**
       * @brief Returns the affine map order.
       * @returns One.
       */
      size_t getOrder() const override
      {
        return 1;
      }

      /**
       * @brief Applies the scalar dilation.
       * @param pc Output physical coordinates.
       * @param rc Reference coordinates.
       */
      void transform(Math::SpatialPoint& pc, const Math::SpatialPoint& rc) const override
      {
        pc = rc * m_scale;
      }

      /**
       * @brief Evaluates the constant map derivative.
       * @param jacobian Output derivative matrix.
       * @param rc Unused reference coordinates.
       */
      void jacobian(Math::SpatialMatrix<Real>& jacobian,
        [[maybe_unused]] const Math::SpatialPoint& rc) const override
      {
        jacobian.resize(1, 1);
        jacobian(0, 0) = m_scale;
      }

      /**
       * @brief Copies the map polymorphically.
       * @returns A new map owned by the caller.
       */
      StateTransformation* copy() const noexcept override
      {
        return new StateTransformation(*this);
      }

    private:
      Real m_scale;
  };

  /// @brief Transformations retain dimensions in constant-evaluated construction.
  constexpr bool constantMapDimensions()
  {
    StateTransformation source;
    StateTransformation moved(std::move(source));
    return moved.getReferenceDimension() == 1 && moved.getPhysicalDimension() == 1;
  }
  static_assert(constantMapDimensions());
  static_assert(std::is_nothrow_move_constructible_v<StateTransformation>);
  static_assert(std::is_nothrow_move_assignable_v<StateTransformation>);

  /// @brief Copying and moving preserve independently owned geometry parameters.
  TEST(Geometry_IdentityTransformation, StateGeometryLifecycle)
  {
    StateTransformation source;
    source.setScale(2);
    StateTransformation copied(source);
    source.setScale(4);
    StateTransformation moved(std::move(copied));
    Math::SpatialPoint result;
    moved.transform(result, Math::SpatialPoint{Real(0.5)});
    EXPECT_EQ(result(0), 1);
    moved = std::move(source);
    moved.transform(result, Math::SpatialPoint{Real(0.5)});
    EXPECT_EQ(result(0), 2);
  }

  /// @brief Archive loading visits both geometric dimensions without runtime state tokens.
  TEST(Geometry_IdentityTransformation, LoadingDimensions)
  {
    struct LoadingArchive
    {
      size_t visits = 0;
      LoadingArchive& operator&(size_t& dimension)
      {
        dimension = 2;
        ++visits;
        return *this;
      }
    };
    StateTransformation map;
    LoadingArchive archive;
    map.serialize(archive, 0);
    EXPECT_EQ(archive.visits, 2);
    EXPECT_EQ(map.getReferenceDimension(), 2);
    EXPECT_EQ(map.getPhysicalDimension(), 2);
  }

  // ---- IdentityTransformation basic properties ----

  /// @brief Verifies construction for geometry identity transformation by checking exact expected values.
  TEST(Geometry_IdentityTransformation, Construction)
  {
    IdentityTransformation t(2);
    EXPECT_EQ(t.getReferenceDimension(), 2);
    EXPECT_EQ(t.getPhysicalDimension(), 2);
    EXPECT_EQ(t.getOrder(), 1);
  }

  /// @brief Verifies construction 3 D for geometry identity transformation by checking exact expected values.
  TEST(Geometry_IdentityTransformation, Construction3D)
  {
    IdentityTransformation t(3);
    EXPECT_EQ(t.getReferenceDimension(), 3);
    EXPECT_EQ(t.getPhysicalDimension(), 3);
  }

  /// @brief Verifies transform 2 D for geometry identity transformation by checking tolerance-based numerical results.
  TEST(Geometry_IdentityTransformation, Transform2D)
  {
    IdentityTransformation t(2);
    Math::SpatialPoint rc{0.3, 0.7};
    Math::SpatialPoint pc(2);
    t.transform(pc, rc);
    EXPECT_NEAR(pc(0), 0.3, 1e-14);
    EXPECT_NEAR(pc(1), 0.7, 1e-14);
  }

  /// @brief Verifies transform 3 D for geometry identity transformation by checking tolerance-based numerical results.
  TEST(Geometry_IdentityTransformation, Transform3D)
  {
    IdentityTransformation t(3);
    Math::SpatialPoint rc{0.1, 0.2, 0.3};
    Math::SpatialPoint pc(3);
    t.transform(pc, rc);
    EXPECT_NEAR(pc(0), 0.1, 1e-14);
    EXPECT_NEAR(pc(1), 0.2, 1e-14);
    EXPECT_NEAR(pc(2), 0.3, 1e-14);
  }

  /// @brief Verifies jacobian 2 D for geometry identity transformation by checking tolerance-based numerical results.
  TEST(Geometry_IdentityTransformation, Jacobian2D)
  {
    IdentityTransformation t(2);
    Math::SpatialPoint rc{0.5, 0.5};
    Math::SpatialMatrix<Real> jac(2, 2);
    t.jacobian(jac, rc);
    // Should be identity matrix
    EXPECT_NEAR(jac(0, 0), 1.0, 1e-14);
    EXPECT_NEAR(jac(0, 1), 0.0, 1e-14);
    EXPECT_NEAR(jac(1, 0), 0.0, 1e-14);
    EXPECT_NEAR(jac(1, 1), 1.0, 1e-14);
  }

  /// @brief Verifies jacobian 3 D for geometry identity transformation by checking tolerance-based numerical results.
  TEST(Geometry_IdentityTransformation, Jacobian3D)
  {
    IdentityTransformation t(3);
    Math::SpatialPoint rc{0.1, 0.2, 0.3};
    Math::SpatialMatrix<Real> jac(3, 3);
    t.jacobian(jac, rc);
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
        EXPECT_NEAR(jac(i, j), (i == j) ? 1.0 : 0.0, 1e-14);
    }
  }

  /**
   * @brief Sizes identity Jacobians independently of caller storage.
   *
   * Both empty and incorrectly sized output matrices must become the identity
   * on the transformation's actual reference domain, including dimension zero.
   * Exact dimensions and entries are checked without a numerical tolerance.
   */
  TEST(Geometry_IdentityTransformation, JacobianResizesCallerStorage)
  {
    for (size_t dimension = 0; dimension <= 3; ++dimension)
    {
      SCOPED_TRACE(dimension);
      IdentityTransformation transformation(dimension);
      Math::SpatialPoint reference(dimension);
      reference.setZero();
      for (const size_t initialDimension : {size_t(0), size_t(1), size_t(3)})
      {
        Math::SpatialMatrix<Real> matrix(initialDimension, initialDimension);
        transformation.jacobian(matrix, reference);
        ASSERT_EQ(static_cast<size_t>(matrix.rows()), dimension);
        ASSERT_EQ(static_cast<size_t>(matrix.cols()), dimension);
        for (size_t j = 0; j < dimension; ++j)
        {
          for (size_t i = 0; i < dimension; ++i)
            EXPECT_EQ(matrix(j, i), Real(j == i));
        }
      }
    }
  }

  /// @brief Verifies copy for geometry identity transformation by checking tolerance-based numerical results, exact expected values, copy semantics.
  TEST(Geometry_IdentityTransformation, Copy)
  {
    IdentityTransformation t(2);
    auto* c = t.copy();
    EXPECT_EQ(c->getReferenceDimension(), 2);
    EXPECT_EQ(c->getPhysicalDimension(), 2);
    EXPECT_EQ(c->getOrder(), 1);

    Math::SpatialPoint rc{0.4, 0.6};
    Math::SpatialPoint pc(2);
    c->transform(pc, rc);
    EXPECT_NEAR(pc(0), 0.4, 1e-14);
    EXPECT_NEAR(pc(1), 0.6, 1e-14);
    delete c;
  }

  /// @brief Verifies move construction for geometry identity transformation by checking exact expected values, move semantics.
  TEST(Geometry_IdentityTransformation, MoveConstruction)
  {
    IdentityTransformation t(3);
    IdentityTransformation moved(std::move(t));
    EXPECT_EQ(moved.getReferenceDimension(), 3);
    EXPECT_EQ(moved.getPhysicalDimension(), 3);
  }

  /// @brief Verifies copy construction for geometry identity transformation by checking exact expected values, copy semantics.
  TEST(Geometry_IdentityTransformation, CopyConstruction)
  {
    IdentityTransformation t(2);
    IdentityTransformation copy(t);
    EXPECT_EQ(copy.getReferenceDimension(), 2);
    EXPECT_EQ(copy.getPhysicalDimension(), 2);
  }
}
