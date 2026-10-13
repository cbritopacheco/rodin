/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#include <array>
#include <memory>

#include <gtest/gtest.h>

#include <Rodin/Variational.h>

namespace Rodin::Tests::Unit
{
  /// @brief Selects the integration region through the public extension interface.
  template <class Field>
  class GridFunctionRegionQuadrature final : public Variational::QuadratureRule<Field>
  {
    public:
      using Parent = Variational::QuadratureRule<Field>;

      /**
       * @brief Constructs a rule borrowing a field over the selected region.
       * @param field Field integrated by the rule.
       * @param region Region of the mesh to integrate.
       */
      GridFunctionRegionQuadrature(const Field& field, Geometry::Region region)
        : Parent(field),
          m_region(region)
      {}

      /**
       * @brief Copies the rule and its borrowed field reference.
       * @param other Rule to copy.
       */
      GridFunctionRegionQuadrature(const GridFunctionRegionQuadrature& other)
        : Parent(other),
          m_region(other.m_region)
      {}

      /**
       * @brief Moves the rule and its borrowed field reference.
       * @param other Rule to move.
       */
      GridFunctionRegionQuadrature(GridFunctionRegionQuadrature&& other)
        : Parent(std::move(other)),
          m_region(other.m_region)
      {}

      /**
       * @brief Returns the selected mesh region.
       * @returns Integration region.
       */
      Geometry::Region getRegion() const override
      {
        return m_region;
      }

      /**
       * @brief Clones the rule through its polymorphic interface.
       * @returns Owned copy of this rule.
       */
      GridFunctionRegionQuadrature* copy() const noexcept override
      {
        return new GridFunctionRegionQuadrature(*this);
      }

    private:
      const Geometry::Region m_region;  ///< Mesh region selected for integration.
  };

  template <class Scalar>
  class GridFunctionRegionQuadratureTest : public ::testing::Test
  {};

  using RegionQuadratureScalars = ::testing::Types<Real, Complex>;
  TYPED_TEST_SUITE(GridFunctionRegionQuadratureTest, RegionQuadratureScalars);

  /// @brief Exact affine moments distinguish all regions and attribute subsets.
  TYPED_TEST(GridFunctionRegionQuadratureTest, AffineMomentsAndLifecycle)
  {
    using Scalar = TypeParam;
    using Region = Geometry::Region;
    constexpr Real Tolerance = 1e-12;
    constexpr Geometry::Attribute Selected = 7, Other = 8;
    auto mesh = Geometry::LocalMesh::Builder()
                  .initialize(2)
                  .nodes(6)
                  .vertex({0, 0})
                  .vertex({0.5, 0})
                  .vertex({1, 0})
                  .vertex({0, 1})
                  .vertex({0.5, 1})
                  .vertex({1, 1})
                  .polytope(Geometry::Polytope::Type::Quadrilateral, {0, 1, 4, 3})
                  .polytope(Geometry::Polytope::Type::Quadrilateral, {1, 2, 5, 4})
                  .finalize();
    mesh.getConnectivity().compute(1, 2);
    ASSERT_EQ(mesh.getCellCount(), 2);
    ASSERT_EQ(mesh.getFaceCount(), 7);
    for (auto cell = mesh.getPolytope(2); cell; ++cell)
    {
      const Geometry::Point center(
        *cell, Geometry::Polytope::Traits(cell->getGeometry()).getCentroid());
      mesh.setAttribute({2, cell->getIndex()}, center.x() < 0.5 ? Selected : Other);
    }
    for (auto face = mesh.getPolytope(1); face; ++face)
    {
      const Geometry::Point center(
        *face, Geometry::Polytope::Traits(face->getGeometry()).getCentroid());
      mesh.setAttribute(
        {1, face->getIndex()}, center.x() == 0 || center.y() == 0 ? Selected : Other);
    }
    Variational::P1<Scalar> space(mesh);
    Variational::GridFunction field(space);
    const Scalar imaginary = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(0, 1);
      else
        return Scalar(0);
    }();
    const Scalar a = Scalar(2) + imaginary;
    const Scalar b = Scalar(3) - Scalar(2) * imaginary;
    const Scalar c = Scalar(-1) + Scalar(4) * imaginary;
    field = [=](
              const Geometry::Point& point) { return a + b * point.x() + c * point.y(); };

    // Moments of (1,x,y) on the unit square with its vertical midline.
    // Selected cells form the left half; selected faces are left and bottom.
    const std::array<Region, 4> regions{
      Region::Cells, Region::Boundary, Region::Faces, Region::Interface};
    const std::array<std::array<Real, 3>, 4> full{
      {{1, 0.5, 0.5}, {4, 2, 2}, {5, 2.5, 2.5}, {1, 0.5, 0.5}}};
    const std::array<std::array<Real, 3>, 4> selected{
      {{0.5, 0.125, 0.25}, {2, 0.5, 0.5}, {2, 0.5, 0.5}, {0, 0, 0}}};
    using Rule = GridFunctionRegionQuadrature<decltype(field)>;
    for (size_t index = 0; index < regions.size(); ++index)
    {
      SCOPED_TRACE(index);
      const auto moment = [&](const std::array<Real, 3>& values) {
        return a * values[0] + b * values[1] + c * values[2];
      };
      for (const bool restricted : {false, true})
      {
        SCOPED_TRACE(restricted);
        Rule rule(field, regions[index]);
        rule.setOrder(1);
        if (restricted)
          rule.over(Selected);
        EXPECT_FALSE(rule.getValue().has_value());
        EXPECT_EQ(rule.getType(), Variational::Integrator::Type::Linear);
        const Scalar expected = moment(restricted ? selected[index] : full[index]);
        EXPECT_LE(std::abs(rule.compute() - expected), Tolerance);
        ASSERT_TRUE(rule.getValue().has_value());
        EXPECT_LE(std::abs(*rule.getValue() - expected), Tolerance);
        EXPECT_LE(std::abs(static_cast<Scalar>(rule) - expected), Tolerance);

        Rule copied(rule);
        EXPECT_FALSE(copied.getValue().has_value());
        EXPECT_LE(std::abs(static_cast<Scalar>(copied) - expected), Tolerance);
        std::unique_ptr<Rule> clone(rule.copy());
        Rule moved(std::move(*clone));
        clone.reset();
        EXPECT_LE(std::abs(moved.compute() - expected), Tolerance);
      }
      Rule all(field, regions[index]);
      all.over(FlatSet<Geometry::Attribute>{Selected, Other});
      EXPECT_LE(std::abs(all.compute() - moment(full[index])), Tolerance);
      Rule empty(field, regions[index]);
      empty.over(Geometry::Attribute(999));
      EXPECT_LE(std::abs(empty.compute()), Tolerance);
    }
  }
}
