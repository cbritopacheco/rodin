/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <algorithm>
#include <random>
#include <Rodin/Kokkos/Geometry/Connectivity.h>
#include "../TestContext.h"

namespace Rodin::Tests::Unit
{
  using HostConnectivity = Geometry::Connectivity<Context::Local>;
  using ResidentConnectivity = Geometry::Connectivity<TestContext>;
  using GeometryType = Geometry::Polytope::Type;

  /**
   * @param host Independent CPU reference.
   * @param resident Resident implementation under test.
   */
  void checkConnectivity(
    const HostConnectivity& host, const ResidentConnectivity& resident)
  {
    ASSERT_EQ(host.getDimension(), resident.getDimension());
    for (size_t d = 0; d <= host.getDimension(); ++d)
    {
      ASSERT_EQ(host.getCount(d), resident.getCount(d));
      for (Index i = 0; i < host.getCount(d); ++i)
      {
        EXPECT_EQ(host.getGeometry(d, i), resident.getGeometry(d, i));
        const auto& a = host.getPolytope(d, i);
        const auto& b = resident.getPolytope(d, i);
        ASSERT_EQ(a.size(), b.size());
        for (Index k = 0; k < a.size(); ++k)
          EXPECT_EQ(a[k], b[k]);
      }
      for (size_t dp = 0; dp <= host.getDimension(); ++dp)
      {
        SCOPED_TRACE(::testing::Message() << "relation " << d << ", " << dp);
        EXPECT_EQ(host.isDirty(d, dp), resident.isDirty(d, dp));
        EXPECT_EQ(host.getIncidence(d, dp), resident.getIncidence(d, dp));
      }
    }
    for (size_t g = 0; g < 8; ++g)
    {
      EXPECT_EQ(host.getCount(static_cast<GeometryType>(g)),
        resident.getCount(static_cast<GeometryType>(g)));
    }
  }

  /**
   * @param host Independent CPU relation oracle.
   * @param resident Resident tables read through their device descriptor interface.
   */
  void checkDeviceIncidence(
    const HostConnectivity& host, const ResidentConnectivity& resident)
  {
    const auto view = resident.getDeviceView();
    for (size_t d = 0; d <= host.getDimension(); ++d)
    {
      for (size_t dp = 0; dp <= host.getDimension(); ++dp)
      {
        const auto& expected = host.getIncidence(d, dp);
        Index width = 1;
        for (const auto& row : expected)
          width = std::max(width, Index(row.size() + 1));
        ::Kokkos::View<Index**, TestMemorySpace> output(
          "computed device incidence", expected.size(), width);
        ::Kokkos::parallel_for(
          "read computed device CSR",
          ::Kokkos::RangePolicy<TestExecutionSpace>(
            resident.getContext().getExecutionSpace(), 0, expected.size()),
          KOKKOS_LAMBDA(Index i) {
            const auto row = view.getIncidence(d, dp, i);
            output(i, 0) = row.size();
            for (Index j = 0; j < row.size(); ++j)
              output(i, j + 1) = row[j];
          });
        const auto values =
          ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, output);
        for (Index i = 0; i < expected.size(); ++i)
        {
          ASSERT_EQ(values(i, 0), expected[i].size());
          for (Index j = 0; j < expected[i].size(); ++j)
            EXPECT_EQ(values(i, j + 1), expected[i][j]);
        }
      }
    }
  }

  /// @param geometry Reference family. @returns Two cells sharing vertices in different local orientations.
  HostConnectivity sharedCells(GeometryType geometry)
  {
    const Geometry::Polytope::Traits traits(geometry);
    HostConnectivity result;
    const Index n = traits.getVertexCount();
    result.initialize(traits.getDimension()).nodes(n + 1);
    if (geometry != GeometryType::Point)
    {
      Geometry::Polytope::Key a(n), b(n);
      for (Index v = 0; v < n; ++v)
        a[v] = b[v] = v;
      std::swap(b[0], b[1]);
      b[n - 1] = n;
      result.polytope(geometry, a).polytope(geometry, b);
    }
    return result;
  }

  /// @brief Every dimension pair, query order and reference family agrees exactly with Local connectivity.
  TEST(Kokkos_Connectivity, AllReferenceFamiliesAndQueryOrders)
  {
    for (size_t g = 0; g < 8; ++g)
    {
      for (size_t order = 0; order < 3; ++order)
      {
        SCOPED_TRACE(
          ::testing::Message() << "geometry " << g << ", query order " << order);
        auto host = sharedCells(static_cast<GeometryType>(g));
        ResidentConnectivity resident(TestContext{}, host);
        checkConnectivity(host, resident);
        std::vector<std::pair<size_t, size_t>> queries;
        for (size_t d = 0; d <= host.getDimension(); ++d)
        {
          for (size_t dp = 0; dp <= host.getDimension(); ++dp)
            queries.emplace_back(d, dp);
        }
        if (order == 1)
          std::reverse(queries.begin(), queries.end());
        if (order == 2)
        {
          std::mt19937 random(17);
          std::shuffle(queries.begin(), queries.end(), random);
        }
        for (const auto& [d, dp] : queries)
        {
          SCOPED_TRACE(::testing::Message() << "compute " << d << ", " << dp);
          host.compute(d, dp);
          EXPECT_EQ(&resident.compute(d, dp), &resident);
          checkConnectivity(host, resident);
          // Clean repeated queries must retain IDs and ordering.
          resident.compute(d, dp);
          checkConnectivity(host, resident);
        }
        checkDeviceIncidence(host, resident);
      }
    }
  }

  /// @brief Mixed cell families retain first-seen numbering and ordered vertices.
  TEST(Kokkos_Connectivity, MixedThreeDimensionalComplex)
  {
    HostConnectivity host;
    host.initialize(3)
      .nodes(12)
      .polytope(GeometryType::Hexahedron, {0, 1, 2, 3, 4, 5, 6, 7})
      .polytope(GeometryType::Wedge, {0, 1, 8, 4, 5, 9})
      .polytope(GeometryType::Pyramid, {4, 5, 6, 7, 10})
      .polytope(GeometryType::Tetrahedron, {4, 5, 10, 11});
    ResidentConnectivity resident(TestContext{}, host);
    for (const auto& pair :
      {std::pair<size_t, size_t>{2, 3}, {1, 2}, {0, 1}, {2, 2}, {1, 1}, {0, 0}, {3, 3}})
    {
      host.compute(pair.first, pair.second);
      resident.compute(pair.first, pair.second);
      checkConnectivity(host, resident);
    }
  }

  /// @brief Restriction uses only pre-existing entities, preserves their orientation, and obeys cached-query semantics.
  TEST(Kokkos_Connectivity, RestrictionAndExplicitRecomputation)
  {
    HostConnectivity host;
    host.initialize(2)
      .nodes(4)
      .polytope(GeometryType::Segment, {1, 0})
      .polytope(GeometryType::Segment, {2, 1})
      .polytope(GeometryType::Triangle, {0, 1, 2})
      .polytope(GeometryType::Triangle, {1, 3, 2});
    ResidentConnectivity resident(TestContext{}, host);
    host.restrict(2, 1);
    resident.restrict(2, 1);
    checkConnectivity(host, resident);
    EXPECT_EQ(resident.getCount(1), 2u);
    host.discover(2, 1);
    resident.discover(2, 1);
    checkConnectivity(host, resident);
    EXPECT_EQ(resident.getCount(1), 2u);
    host.clear(2, 1).discover(2, 1);
    resident.clear(2, 1).discover(2, 1);
    checkConnectivity(host, resident);
    EXPECT_EQ(resident.getCount(1), 5u);
    host.clear(1, 2).compute(1, 2);
    resident.clear(1, 2).compute(1, 2);
    checkConnectivity(host, resident);
    host.clear(2, 1).restrict(2, 1);
    resident.clear(2, 1).restrict(2, 1);
    checkConnectivity(host, resident);
  }

  /// @brief Transposition removes duplicate contributions even with nonmanifold high-valence incidence.
  TEST(Kokkos_Connectivity, NonmanifoldAndDuplicateIncidences)
  {
    constexpr Index count = 70;
    HostConnectivity host;
    host.initialize(2).nodes(count + 2);
    for (Index i = 0; i < count; ++i)
      host.polytope(GeometryType::Triangle, {0, 1, i + 2});
    ResidentConnectivity resident(TestContext{}, host);
    for (const auto& pair : {std::pair<size_t, size_t>{1, 2}, {2, 2}, {0, 0}, {1, 1}})
    {
      host.compute(pair.first, pair.second);
      resident.compute(pair.first, pair.second);
      checkConnectivity(host, resident);
    }
    auto rows = host.getIncidence(2, 0);
    for (auto& row : rows)
      row.push_back(row.front());
    host.setIncidence({2, 0}, std::move(rows));
    ResidentConnectivity duplicates(TestContext{}, host);
    host.clear(0, 2).transpose(0, 2);
    duplicates.clear(0, 2).transpose(0, 2);
    checkConnectivity(host, duplicates);
    checkDeviceIncidence(host, duplicates);
  }

  /// @brief Copies independently compute new relations and dirty dimensions reject invalid input.
  TEST(Kokkos_Connectivity, CopyAndInvalidDimensions)
  {
    auto host = sharedCells(GeometryType::Tetrahedron);
    ResidentConnectivity resident(TestContext{}, host);
    ResidentConnectivity copy(resident);
    copy.compute(2, 3);
    EXPECT_EQ(resident.getCount(2), 0u);
    host.compute(2, 3);
    checkConnectivity(host, copy);
    ResidentConnectivity moved(std::move(copy));
    moved.compute(1, 2);
    host.compute(1, 2);
    checkConnectivity(host, moved);
    EXPECT_THROW(resident.compute(4, 0), Alert::Exception);
    EXPECT_THROW(resident.isDirty(0, 4), Alert::Exception);
    EXPECT_THROW(resident.setDirty(4, 0), Alert::Exception);
    EXPECT_THROW(resident.clear(4, 0), Alert::Exception);
    resident.setDirty(2, 1, false);
    EXPECT_FALSE(resident.isDirty(2, 1));
  }

  /// @brief Empty resident connectivity needs no host topology or special discovery path.
  TEST(Kokkos_Connectivity, EmptyComplex)
  {
    ResidentConnectivity resident;
    EXPECT_EQ(&resident.compute(0, 0), &resident);
    EXPECT_EQ(resident.getCount(0), 0u);
    EXPECT_TRUE(resident.getIncidence(0, 0).empty());
  }
}
