/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <array>
#include <thread>

#include <Rodin/Geometry.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  /// @brief Mutable test map obeying the transformation-state invalidation contract.
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
      invalidateCacheIdentity();
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

  /**
   * @brief Checks the constant-evaluated construction contract.
   * @returns Whether dimensions are retained without inventing a runtime token.
   */
  constexpr bool constantMapHasNoRuntimeIdentity()
  {
    StateTransformation map;
    return map.getReferenceDimension() == 1 && !map.getCacheIdentity();
  }
  static_assert(constantMapHasNoRuntimeIdentity());
  /**
   * @brief Checks constant-evaluated move construction without a runtime token.
   * @returns Whether both map states remain non-reusable in constant evaluation.
   */
  constexpr bool constantMapMoveHasNoRuntimeIdentity()
  {
    StateTransformation source;
    StateTransformation moved(std::move(source));
    return !source.getCacheIdentity() && !moved.getCacheIdentity() &&
      moved.getReferenceDimension() == 1;
  }
  static_assert(constantMapMoveHasNoRuntimeIdentity());
  static_assert(std::is_nothrow_move_constructible_v<StateTransformation>);
  static_assert(std::is_nothrow_move_assignable_v<StateTransformation>);

  /// @brief Construction, value transfer and mutation never reuse stale map tokens.
  TEST(Geometry_IdentityTransformation, StateIdentityLifecycle)
  {
    StateTransformation source;
    ASSERT_TRUE(source.getCacheIdentity());
    const auto initial = source.getCacheIdentity();
    StateTransformation copied(source);
    EXPECT_NE(copied.getCacheIdentity(), initial);
    StateTransformation moved(std::move(copied));
    ASSERT_TRUE(moved.getCacheIdentity());
    EXPECT_FALSE(copied.getCacheIdentity());
    EXPECT_NE(moved.getCacheIdentity(), initial);
    const auto previous = moved.getCacheIdentity();
    moved.setScale(2);
    EXPECT_NE(moved.getCacheIdentity(), previous);
    Math::SpatialMatrix<Real> jacobian;
    moved.jacobian(jacobian, Math::SpatialPoint{0});
    EXPECT_EQ(jacobian(0, 0), 2);
    const auto beforeAssignment = moved.getCacheIdentity();
    moved = std::move(source);
    EXPECT_NE(moved.getCacheIdentity(), beforeAssignment);
    EXPECT_FALSE(source.getCacheIdentity());
    EXPECT_NE(moved.getCacheIdentity(), initial);
  }

  /// @brief Archive loading refreshes local state without serializing its token.
  TEST(Geometry_IdentityTransformation, LoadingRefreshesStateIdentity)
  {
    struct LoadingArchive
    {
      using is_loading = std::true_type;
      LoadingArchive()
        : m_visits(0)
      {}
      LoadingArchive& operator&([[maybe_unused]] size_t& dimension)
      {
        ++m_visits;
        return *this;
      }
      size_t m_visits;
    };
    static_assert(LoadingArchive::is_loading::value);
    StateTransformation map;
    const auto before = map.getCacheIdentity();
    LoadingArchive archive;
    map.serialize(archive, 0);
    EXPECT_EQ(archive.m_visits, 2);
    EXPECT_NE(map.getCacheIdentity(), before);
    EXPECT_EQ(map.getReferenceDimension(), 1);
    EXPECT_EQ(map.getPhysicalDimension(), 1);
  }

  /// @brief Runtime state tokens remain distinct across construction threads.
  TEST(Geometry_IdentityTransformation, StateIdentityConcurrentConstruction)
  {
    constexpr size_t Threads = 4;
    constexpr size_t MapsPerThread = 64;
    std::array<std::array<size_t, MapsPerThread>, Threads> values{};
    std::array<std::thread, Threads> threads;
    for (size_t r = 0; r < Threads; ++r)
    {
      threads[r] = std::thread([&, r]() {
        for (size_t i = 0; i < MapsPerThread; ++i)
        {
          StateTransformation map;
          ASSERT_TRUE(map.getCacheIdentity());
          values[r][i] = *map.getCacheIdentity();
        }
      });
    }
    for (auto& thread : threads)
      thread.join();
    FlatSet<size_t> unique;
    for (const auto& row : values)
    {
      for (const auto identity : row)
        EXPECT_TRUE(unique.insert(identity).second);
    }
    EXPECT_EQ(unique.size(), Threads * MapsPerThread);
  }

  /// @brief Reconstructing a map in the same storage receives a different token.
  TEST(Geometry_IdentityTransformation, StateIdentityReusedStorage)
  {
    Optional<IdentityTransformation> map;
    map.emplace(2);
    const auto* address = &*map;
    const auto before = map->getCacheIdentity();
    map.emplace(2);
    EXPECT_EQ(&*map, address);
    ASSERT_TRUE(map->getCacheIdentity());
    EXPECT_NE(map->getCacheIdentity(), before);
  }

  /// @brief A failed archive load cannot leave the previous cache token reusable.
  TEST(Geometry_IdentityTransformation, FailedLoadingRefreshesStateIdentity)
  {
    struct FailingArchive
    {
      using is_loading = std::true_type;
      FailingArchive& operator&(size_t& dimension)
      {
        dimension = 2;
        Alert::Exception() << "Deliberately rejected test archive" << Alert::Raise;
        return *this;
      }
    };
    static_assert(FailingArchive::is_loading::value);
    StateTransformation map;
    const auto before = map.getCacheIdentity();
    FailingArchive archive;
    EXPECT_THROW(map.serialize(archive, 0), Alert::Exception);
    EXPECT_NE(map.getCacheIdentity(), before);
    EXPECT_EQ(map.getReferenceDimension(), 2);
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
