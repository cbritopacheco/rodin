/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <Rodin/Geometry.h>
#include <Rodin/Alert/Exception.h>

#include <algorithm>
#include <functional>
#include <limits>
#include <type_traits>
#include <utility>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  class Rodin_Geometry_MinSTCut : public testing::Test
  {
    protected:
      static std::function<Real(const Polytope&)> average(std::vector<Real> values)
      {
        return [values = std::move(values)](const Polytope& cell)
        {
          return values.at(cell.getIndex());
        };
      }
  };

  class MinSTCutMeshTest : public Rodin_Geometry_MinSTCut,
    public testing::WithParamInterface<Polytope::Type>
  {
    protected:
      MinSTCutMeshTest()
        : m_mesh(GetParam() == Polytope::Type::Segment
            ? LocalMesh::UniformGrid(GetParam(), {3})
            : GetParam() == Polytope::Type::Triangle ||
              GetParam() == Polytope::Type::Quadrilateral
              ? LocalMesh::UniformGrid(GetParam(), {3, 2})
              : GetParam() == Polytope::Type::Hexahedron
                ? LocalMesh::UniformGrid(GetParam(), {3, 2, 2})
                : LocalMesh::UniformGrid(GetParam(), {2, 2, 2}))
      {
        m_mesh.getConnectivity().compute(m_mesh.getDimension() - 1, m_mesh.getDimension());
      }

      LocalMesh m_mesh;
  };

  TEST_P(MinSTCutMeshTest, UniformSignsAndBorrowedMesh)
  {
    MinSTCut classifier(m_mesh);
    EXPECT_EQ(&classifier.getMesh(), &m_mesh);
    for (int sign : {-1, 1})
    {
      const auto result = classifier.classify(
        average(std::vector<Real>(m_mesh.getCellCount(), sign)));
      EXPECT_TRUE(result.cut.empty());
      EXPECT_DOUBLE_EQ(result.energy, 0);
      EXPECT_EQ(result.inside.size(), sign < 0 ? m_mesh.getCellCount() : 0);
      EXPECT_EQ(result.outside.size(), sign > 0 ? m_mesh.getCellCount() : 0);
    }
  }

  TEST_P(MinSTCutMeshTest, CapturingLambdaReceivesInteriorMeshFacetsOnce)
  {
    IndexVector visited;
    const Real weight = 0.01;
    MinSTCut<LocalMesh>::Parameters parameters;
    parameters.fidelity = 2;
    parameters.smoothing = [&visited, weight, this](const Polytope& facet)
    {
      EXPECT_EQ(&facet.getMesh(), &m_mesh);
      EXPECT_EQ(facet.getDimension(), m_mesh.getDimension() - 1);
      const auto& cells = m_mesh.getConnectivity().getIncidence(
        {m_mesh.getDimension() - 1, m_mesh.getDimension()}, facet.getIndex());
      EXPECT_EQ(cells.size(), 2u);
      visited.push_back(facet.getIndex());
      return weight;
    };
    std::vector<Real> averages(m_mesh.getCellCount());
    for (Index i = 0; i < averages.size(); ++i)
      averages[i] = i % 2 ? 1 : -1;
    MinSTCut classifier(m_mesh);
    classifier.setParameters(parameters);
    IndexVector cellsVisited;
    const auto result = classifier.classify([&](const Polytope& cell)
    {
      EXPECT_EQ(&cell.getMesh(), &m_mesh);
      EXPECT_EQ(cell.getDimension(), m_mesh.getDimension());
      cellsVisited.push_back(cell.getIndex());
      return averages[cell.getIndex()];
    });
    ASSERT_EQ(cellsVisited.size(), m_mesh.getCellCount());
    for (Index i = 0; i < cellsVisited.size(); ++i)
      EXPECT_EQ(cellsVisited[i], i);
    IndexVector expected;
    Real energy = 0;
    const auto& incidence = m_mesh.getConnectivity().getIncidence(
      m_mesh.getDimension() - 1, m_mesh.getDimension());
    for (auto face = m_mesh.getFace(); face; ++face)
    {
      const auto& cells = incidence.at(face->getIndex());
      if (cells.size() != 2)
        continue;
      expected.push_back(face->getIndex());
      if ((std::find(result.inside.begin(), result.inside.end(), cells[0]) !=
          result.inside.end()) !=
          (std::find(result.inside.begin(), result.inside.end(), cells[1]) !=
            result.inside.end()))
      {
        EXPECT_NE(std::find(result.cut.begin(), result.cut.end(), face->getIndex()),
          result.cut.end());
        energy += weight * (m_mesh.getDimension() == 1 ? Real(1) : face->getMeasure());
      }
    }
    EXPECT_EQ(visited, expected);
    for (auto cell = m_mesh.getCell(); cell; ++cell)
    {
      const Index i = cell->getIndex();
      const bool inside = std::find(result.inside.begin(), result.inside.end(), i) !=
        result.inside.end();
      energy += parameters.fidelity * (inside
        ? cell->getMeasure() * std::max(Real(0), averages[i])
        : cell->getMeasure() * std::max(Real(0), -averages[i]));
    }
    EXPECT_NEAR(result.energy, energy, 1e-14);
  }

  TEST_P(MinSTCutMeshTest, CutMatchesExhaustiveEnumeration)
  {
    MinSTCut<LocalMesh>::Parameters parameters;
    parameters.fidelity = 1.7;
    parameters.smoothing = [](const Polytope& facet)
    {
      return Real(0.03) * (1 + facet.getIndex() % 3);
    };
    std::vector<Real> averages(m_mesh.getCellCount());
    for (Index i = 0; i < averages.size(); ++i)
      averages[i] = (i % 2 ? 1 : -1) * Real(i + 1) / averages.size();
    Real minimum = std::numeric_limits<Real>::infinity();
    const auto& incidence = m_mesh.getConnectivity().getIncidence(
      m_mesh.getDimension() - 1, m_mesh.getDimension());
    for (size_t mask = 0; mask < (size_t(1) << averages.size()); ++mask)
    {
      Real energy = 0;
      for (auto cell = m_mesh.getCell(); cell; ++cell)
      {
        const Index i = cell->getIndex();
        energy += parameters.fidelity * ((mask & (size_t(1) << i))
          ? cell->getMeasure() * std::max(Real(0), averages[i])
          : cell->getMeasure() * std::max(Real(0), -averages[i]));
      }
      for (auto face = m_mesh.getFace(); face; ++face)
      {
        const auto& cells = incidence.at(face->getIndex());
        if (cells.size() == 2 &&
            bool(mask & (size_t(1) << cells[0])) != bool(mask & (size_t(1) << cells[1])))
        {
          energy += parameters.smoothing(*face) *
            (m_mesh.getDimension() == 1 ? Real(1) : face->getMeasure());
        }
      }
      minimum = std::min(minimum, energy);
    }
    const auto result = MinSTCut(m_mesh).setParameters(parameters).classify(average(averages));
    EXPECT_NEAR(result.energy, minimum, 1e-14);
    IndexVector cells = result.inside;
    cells.insert(cells.end(), result.outside.begin(), result.outside.end());
    std::sort(cells.begin(), cells.end());
    ASSERT_EQ(cells.size(), m_mesh.getCellCount());
    for (Index i = 0; i < cells.size(); ++i)
      EXPECT_EQ(cells[i], i);
  }

  TEST_P(MinSTCutMeshTest, ClassificationDoesNotModifyAttributes)
  {
    for (auto cell = m_mesh.getCell(); cell; ++cell)
      m_mesh.setAttribute({m_mesh.getDimension(), cell->getIndex()}, 17);
    for (auto face = m_mesh.getFace(); face; ++face)
      m_mesh.setAttribute({m_mesh.getDimension() - 1, face->getIndex()}, 23);
    MinSTCut(m_mesh).classify(average(std::vector<Real>(m_mesh.getCellCount(), -1)));
    for (auto cell = m_mesh.getCell(); cell; ++cell)
      EXPECT_EQ(cell->getAttribute(), 17);
    for (auto face = m_mesh.getFace(); face; ++face)
      EXPECT_EQ(face->getAttribute(), 23);
  }

  INSTANTIATE_TEST_SUITE_P(Rodin_Geometry_MinSTCut, MinSTCutMeshTest,
    testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron));

  TEST_F(Rodin_Geometry_MinSTCut, OwnedParametersAreCopiedAndChainable)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    mesh.getConnectivity().compute(1, 2);
    MinSTCut<LocalMesh>::Parameters parameters;
    const Real weight = 0.01;
    parameters.smoothing = [weight](const Polytope&) { return weight; };
    MinSTCut classifier(mesh);
    classifier.setParameters(parameters);
    static_assert(std::is_same_v<decltype(classifier), MinSTCut<LocalMesh>>);
    parameters.smoothing = [](const Polytope&) { return 1; };
    EXPECT_DOUBLE_EQ(classifier.getParameters().smoothing(*mesh.getFace()), weight);
    const auto owned = classifier.classify(average({-0.4, 0.6}));
    EXPECT_EQ(owned.inside, IndexVector{0});
    EXPECT_EQ(owned.outside, IndexVector{1});
    const auto updated = classifier.setParameters(parameters).classify(average({-0.4, 0.6}));
    EXPECT_TRUE(updated.inside.empty() || updated.outside.empty());
    parameters.fidelity = 2;
    EXPECT_EQ(&classifier.setParameters(parameters), &classifier);
    parameters.fidelity = 10;
    EXPECT_DOUBLE_EQ(classifier.getParameters().fidelity, 2);
    EXPECT_DOUBLE_EQ(MinSTCut<LocalMesh>::Parameters{}.smoothing(*mesh.getFace()), 1);
  }

  TEST_F(Rodin_Geometry_MinSTCut, FidelityControlsSmoothingBalanceWithoutPinning)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    mesh.getConnectivity().compute(1, 2);
    MinSTCut classifier(mesh);
    const auto smooth = classifier.classify(average({-0.9, 0.9}));
    EXPECT_TRUE(smooth.inside.empty() || smooth.outside.empty());
    EXPECT_TRUE(smooth.cut.empty());
    EXPECT_NEAR(smooth.energy, 0.45, 1e-14);
    MinSTCut<LocalMesh>::Parameters parameters;
    parameters.fidelity = 100;
    const auto split = classifier.setParameters(parameters).classify(average({-0.9, 0.9}));
    EXPECT_EQ(split.inside, IndexVector{0});
    EXPECT_EQ(split.outside, IndexVector{1});
    ASSERT_EQ(split.cut.size(), 1u);
    EXPECT_NEAR(split.energy, mesh.getFace(split.cut[0])->getMeasure(), 1e-14);
  }

  TEST_F(Rodin_Geometry_MinSTCut, CellVolumesAreAppliedOnce)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    mesh.getConnectivity().compute(1, 2);
    MinSTCut<LocalMesh>::Parameters parameters;
    parameters.smoothing = [](const Polytope&) { return 10; };
    const auto result = MinSTCut(mesh).setParameters(parameters).classify(average({-0.25, 0.5}));
    EXPECT_TRUE(result.inside.empty());
    EXPECT_EQ(result.outside, (IndexVector{0, 1}));
    EXPECT_NEAR(result.energy, 0.125, 1e-14);
  }

  TEST_F(Rodin_Geometry_MinSTCut, FacetAttributesCanControlTheCut)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, {5, 2});
    mesh.getConnectivity().compute(1, 2);
    Index weak = 0;
    for (auto face = mesh.getFace(); face; ++face)
    {
      const auto& cells = mesh.getConnectivity().getIncidence({1, 2}, face->getIndex());
      if (cells.size() == 2 && std::min(cells[0], cells[1]) == 1)
      {
        weak = face->getIndex();
        mesh.setAttribute({1, weak}, 42);
      }
    }
    MinSTCut<LocalMesh>::Parameters parameters;
    parameters.smoothing = [](const Polytope& facet)
    {
      return facet.getAttribute().value_or(0) == 42 ? 0.001 : 5.0;
    };
    const auto result = MinSTCut(mesh).setParameters(parameters).classify(
      average({-1, -0.4, 0.4, 1}));
    EXPECT_EQ(result.inside, (IndexVector{0, 1}));
    EXPECT_EQ(result.outside, (IndexVector{2, 3}));
    EXPECT_EQ(result.cut, IndexVector{weak});
  }

  TEST_F(Rodin_Geometry_MinSTCut, ZeroSmoothingFollowsData)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    mesh.getConnectivity().compute(1, 2);
    MinSTCut<LocalMesh>::Parameters parameters;
    parameters.smoothing = [](const Polytope&) { return 0; };
    std::vector<Real> averages(mesh.getCellCount());
    IndexVector inside;
    IndexVector outside;
    for (Index i = 0; i < averages.size(); ++i)
    {
      averages[i] = i % 2 ? 1 : -1;
      (i % 2 ? outside : inside).push_back(i);
    }
    const auto result = MinSTCut(mesh).setParameters(parameters).classify(average(averages));
    EXPECT_EQ(result.inside, inside);
    EXPECT_EQ(result.outside, outside);
    EXPECT_DOUBLE_EQ(result.energy, 0);
  }

  TEST_F(Rodin_Geometry_MinSTCut, InvalidInputsUseRodinAlerts)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    mesh.getConnectivity().compute(1, 2);
    MinSTCut classifier(mesh);
    EXPECT_THROW(classifier.classify({}), Alert::Exception);
    EXPECT_THROW(classifier.classify(average({std::numeric_limits<Real>::quiet_NaN(), 1})),
      Alert::Exception);
    for (Real value : {-1.0, std::numeric_limits<Real>::infinity(),
        std::numeric_limits<Real>::quiet_NaN()})
    {
      MinSTCut<LocalMesh>::Parameters parameters;
      parameters.fidelity = value;
      EXPECT_THROW(classifier.setParameters(parameters).classify(average({-1, 1})),
        Alert::Exception);
      parameters.fidelity = 1;
      parameters.smoothing = [value](const Polytope&) { return value; };
      EXPECT_THROW(classifier.setParameters(parameters).classify(average({-1, 1})),
        Alert::Exception);
    }
    MinSTCut<LocalMesh>::Parameters parameters;
    parameters.smoothing = {};
    EXPECT_THROW(classifier.setParameters(parameters).classify(average({-1, 1})), Alert::Exception);
  }

  TEST_F(Rodin_Geometry_MinSTCut, MissingConnectivityIsReported)
  {
    auto mesh = LocalMesh::Builder().initialize(2).nodes(3)
      .vertex({0, 0}).vertex({1, 0}).vertex({0, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2}).finalize();
    EXPECT_THROW(MinSTCut(mesh).classify(average({-1})), Alert::Exception);
    mesh.getConnectivity().compute(1, 2);
    const auto result = MinSTCut(mesh).classify(average({-1}));
    EXPECT_EQ(result.inside, IndexVector{0});
    EXPECT_TRUE(result.outside.empty());
  }

  TEST_F(Rodin_Geometry_MinSTCut, EmptyMeshHasEmptyClassification)
  {
    auto mesh = LocalMesh::Builder().initialize(2).finalize();
    const auto result = MinSTCut(mesh).classify(average({}));
    EXPECT_TRUE(result.inside.empty());
    EXPECT_TRUE(result.outside.empty());
    EXPECT_TRUE(result.cut.empty());
    EXPECT_DOUBLE_EQ(result.energy, 0);
  }

  TEST_F(Rodin_Geometry_MinSTCut, ZeroDimensionalMeshIsNoOp)
  {
    auto mesh = LocalMesh::Builder().initialize(1).nodes(1).vertex({0}).finalize();
    ASSERT_EQ(mesh.getDimension(), 0u);
    MinSTCut classifier(mesh);
    Index calls = 0;
    const auto result = classifier.classify([&](const Polytope&)
    {
      ++calls;
      return Real(0);
    });
    EXPECT_EQ(calls, 0u);
    EXPECT_TRUE(result.inside.empty());
    EXPECT_TRUE(result.outside.empty());
    EXPECT_TRUE(result.cut.empty());
    EXPECT_DOUBLE_EQ(result.energy, 0);
    EXPECT_NO_THROW(classifier.classify({}));
  }

  TEST_F(Rodin_Geometry_MinSTCut, UnequalCellMeasuresDetermineTheCommonLabel)
  {
    auto mesh = LocalMesh::Builder().initialize(2).nodes(6)
      .vertex({0, 0}).vertex({1, 0}).vertex({3, 0})
      .vertex({0, 1}).vertex({1, 1}).vertex({3, 1})
      .polytope(Polytope::Type::Quadrilateral, {0, 1, 4, 3})
      .polytope(Polytope::Type::Quadrilateral, {1, 2, 5, 4}).finalize();
    mesh.getConnectivity().compute(1, 2);
    MinSTCut<LocalMesh>::Parameters parameters;
    parameters.smoothing = [](const Polytope&) { return 10; };
    const auto result = MinSTCut(mesh).setParameters(parameters).classify(average({-0.8, 0.5}));
    EXPECT_TRUE(result.inside.empty());
    EXPECT_EQ(result.outside, (IndexVector{0, 1}));
    EXPECT_NEAR(result.energy, 0.8, 1e-14);
  }

  TEST_F(Rodin_Geometry_MinSTCut, NonmanifoldFacetsAreRejected)
  {
    auto mesh = LocalMesh::Builder().initialize(2).nodes(5)
      .vertex({0, 0}).vertex({1, 0}).vertex({0, 1})
      .vertex({0, -1}).vertex({1, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .polytope(Polytope::Type::Triangle, {1, 0, 3})
      .polytope(Polytope::Type::Triangle, {0, 1, 4}).finalize();
    mesh.getConnectivity().compute(1, 2);
    EXPECT_THROW(MinSTCut(mesh).classify(average({-1, 1, 1})), Alert::Exception);
  }

  TEST_F(Rodin_Geometry_MinSTCut, GeometryChangesAreNotCachedInTheClassifier)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    mesh.getConnectivity().compute(1, 2);
    MinSTCut<LocalMesh>::Parameters parameters;
    parameters.fidelity = 100;
    MinSTCut classifier(mesh);
    classifier.setParameters(parameters);
    const auto before = classifier.classify(average({-1, 1}));
    for (Index i = 0; i < mesh.getVertexCount(); ++i)
    {
      const auto coordinates = mesh.getVertexCoordinates(i);
      mesh.setVertexCoordinates(i, Real(2) * coordinates(0), 0);
      mesh.setVertexCoordinates(i, Real(2) * coordinates(1), 1);
    }
    const auto after = classifier.classify(average({-1, 1}));
    EXPECT_EQ(before.inside, after.inside);
    EXPECT_EQ(before.outside, after.outside);
    EXPECT_EQ(before.cut, after.cut);
    EXPECT_NEAR(after.energy, 2 * before.energy, 1e-14);
  }
}
