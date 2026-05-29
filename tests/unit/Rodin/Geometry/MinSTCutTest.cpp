/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <Rodin/Geometry.h>

#include <algorithm>
#include <vector>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  namespace
  {
    // Build the edges of a (rows x cols) 4-connected grid of cells with
    // uniform capacity.
    std::vector<MinSTCut::Edge> gridEdges(
        std::size_t rows, std::size_t cols, Real capacity)
    {
      std::vector<MinSTCut::Edge> edges;
      auto idx = [cols](std::size_t r, std::size_t c) -> Index {
        return r * cols + c;
      };
      for (std::size_t r = 0; r < rows; ++r)
      {
        for (std::size_t c = 0; c < cols; ++c)
        {
          if (c + 1 < cols)
            edges.push_back({ idx(r, c), idx(r, c + 1), capacity });
          if (r + 1 < rows)
            edges.push_back({ idx(r, c), idx(r + 1, c), capacity });
        }
      }
      return edges;
    }
  }

  TEST(Rodin_Geometry_MinSTCut, SignConventions)
  {
    EXPECT_EQ(MinSTCut::Inside, -1);
    EXPECT_EQ(MinSTCut::Outside, +1);
  }

  TEST(Rodin_Geometry_MinSTCut, AllNegativeMomentsAreInside)
  {
    const std::vector<Real> volumes(4, 1);
    const std::vector<Real> moments(4, -1);
    const std::vector<MinSTCut::Edge> edges;

    const auto result = MinSTCut().classify(volumes, moments, edges);

    ASSERT_EQ(result.labels.size(), 4u);
    for (const int label : result.labels)
      EXPECT_EQ(label, MinSTCut::Inside);
    EXPECT_EQ(result.insideCells.size(), 4u);
    EXPECT_EQ(result.outsideCells.size(), 0u);
    EXPECT_TRUE(result.cutEdges.empty());
  }

  TEST(Rodin_Geometry_MinSTCut, AllPositiveMomentsAreOutside)
  {
    const std::vector<Real> volumes(4, 1);
    const std::vector<Real> moments(4, +1);
    const std::vector<MinSTCut::Edge> edges;

    const auto result = MinSTCut().classify(volumes, moments, edges);

    for (const int label : result.labels)
      EXPECT_EQ(label, MinSTCut::Outside);
    EXPECT_EQ(result.outsideCells.size(), 4u);
    EXPECT_EQ(result.insideCells.size(), 0u);
  }

  TEST(Rodin_Geometry_MinSTCut, NegativeMomentIsInside)
  {
    // Single isolated cell with a negative moment is labelled Inside.
    const std::vector<Real> volumes{ 1.0 };
    const std::vector<Real> moments{ -0.7 };
    const auto result =
      MinSTCut().classify(volumes, moments, {});
    ASSERT_EQ(result.labels.size(), 1u);
    EXPECT_EQ(result.labels[0], MinSTCut::Inside);
  }

  TEST(Rodin_Geometry_MinSTCut, PositiveMomentIsOutside)
  {
    // Single isolated cell with a positive moment is labelled Outside.
    const std::vector<Real> volumes{ 1.0 };
    const std::vector<Real> moments{ +0.7 };
    const auto result =
      MinSTCut().classify(volumes, moments, {});
    ASSERT_EQ(result.labels.size(), 1u);
    EXPECT_EQ(result.labels[0], MinSTCut::Outside);
  }

  TEST(Rodin_Geometry_MinSTCut, TwoCellStrongEdgeForcesCommonLabel)
  {
    // Two cells with conflicting moments (one negative, one positive)
    // connected by an edge whose capacity dominates the unary terms.
    // The optimal cut keeps both on the same side of the source/sink
    // partition, i.e. they share a label and no cut edge appears.
    const std::vector<Real> volumes{ 1.0, 1.0 };
    const std::vector<Real> moments{ -0.4, +0.4 };
    const std::vector<MinSTCut::Edge> edges{ { 0, 1, 1e3 } };
    const auto result =
      MinSTCut().classify(volumes, moments, edges);
    ASSERT_EQ(result.labels.size(), 2u);
    EXPECT_EQ(result.labels[0], result.labels[1]);
    EXPECT_TRUE(result.cutEdges.empty());
  }

  TEST(Rodin_Geometry_MinSTCut, TwoCellWeakEdgeFollowsData)
  {
    // Same conflicting-moment setup as above but with a weak smoothing
    // edge: the data term wins and the labels split along the moment
    // sign, producing one cut edge.
    const std::vector<Real> volumes{ 1.0, 1.0 };
    const std::vector<Real> moments{ -0.4, +0.4 };
    const std::vector<MinSTCut::Edge> edges{ { 0, 1, 1e-3 } };
    const auto result =
      MinSTCut().classify(volumes, moments, edges);
    ASSERT_EQ(result.labels.size(), 2u);
    EXPECT_EQ(result.labels[0], MinSTCut::Inside);
    EXPECT_EQ(result.labels[1], MinSTCut::Outside);
    ASSERT_EQ(result.cutEdges.size(), 1u);
  }

  TEST(Rodin_Geometry_MinSTCut, Checkerboard2x2SmoothsForLargeLambda)
  {
    // 2x2 grid with checkerboard sign pattern:
    //
    //   (0,0) = -  (0,1) = +
    //   (1,0) = +  (1,1) = -
    //
    // With weak smoothing the cut produces a checkerboard label
    // pattern; with strong smoothing the optimal labelling is uniform
    // (all Inside or all Outside, here Inside because the negative
    // unary costs slightly dominate the matched positives in our toy
    // capacities).
    constexpr std::size_t rows = 2;
    constexpr std::size_t cols = 2;
    const std::vector<Real> volumes(rows * cols, 1.0);
    const std::vector<Real> moments{
      -1.0,  +1.0,
      +1.0,  -1.0
    };
    const auto edges = gridEdges(rows, cols, 1.0);

    {
      MinSTCut::Options weak;
      weak.lambdaScale = 1e-6;
      const auto result =
        MinSTCut().classify(volumes, moments, edges, weak);
      // Checkerboard preserved (each cell matches its own sign).
      EXPECT_EQ(result.labels[0], MinSTCut::Inside);
      EXPECT_EQ(result.labels[1], MinSTCut::Outside);
      EXPECT_EQ(result.labels[2], MinSTCut::Outside);
      EXPECT_EQ(result.labels[3], MinSTCut::Inside);
    }
    {
      MinSTCut::Options strong;
      strong.lambdaScale = 1e6;
      const auto result =
        MinSTCut().classify(volumes, moments, edges, strong);
      // Strong smoothing collapses the 2x2 to a single label.
      EXPECT_EQ(result.labels[0], result.labels[1]);
      EXPECT_EQ(result.labels[1], result.labels[2]);
      EXPECT_EQ(result.labels[2], result.labels[3]);
      EXPECT_TRUE(result.cutEdges.empty());
    }
  }

  TEST(Rodin_Geometry_MinSTCut, UnaryCostMatchesAdvertisedConvention)
  {
    // Inside cost charges only positive moments, outside only negative.
    EXPECT_DOUBLE_EQ(MinSTCut::getInsideCost(2.0, +0.5), 1.0);
    EXPECT_DOUBLE_EQ(MinSTCut::getInsideCost(2.0, -0.5), 0.0);
    EXPECT_DOUBLE_EQ(MinSTCut::getOutsideCost(2.0, -0.5), 1.0);
    EXPECT_DOUBLE_EQ(MinSTCut::getOutsideCost(2.0, +0.5), 0.0);
  }

  TEST(Rodin_Geometry_MinSTCut, SingleInteriorFlipIsHealedByStrongLambda)
  {
    // 3x3 grid: center cell has near-zero negative moment, surrounded by
    // strongly positive cells. Without smoothing the cut yields a
    // checkerboard with a single Inside cell. Sufficient lambda heals it.
    constexpr std::size_t rows = 3;
    constexpr std::size_t cols = 3;
    std::vector<Real> volumes(rows * cols, 1);
    std::vector<Real> moments(rows * cols, +1);
    const Index center = 1 * cols + 1;
    moments[center] = -1e-3;

    const auto baseEdges = gridEdges(rows, cols, 1.0);

    MinSTCut::Options weak;
    weak.lambdaScale = 1e-6;
    const auto weakResult = MinSTCut().classify(volumes, moments, baseEdges, weak);
    EXPECT_EQ(weakResult.labels[center], MinSTCut::Inside);

    MinSTCut::Options strong;
    strong.lambdaScale = 10;
    const auto strongResult =
      MinSTCut().classify(volumes, moments, baseEdges, strong);
    EXPECT_EQ(strongResult.labels[center], MinSTCut::Outside);
  }

  TEST(Rodin_Geometry_MinSTCut, FarFieldThresholdPinsLabels)
  {
    // Two cells: moment +0.9 (definitely Outside) and moment -0.9
    // (definitely Inside), connected by an enormous edge that would
    // otherwise collapse them onto the same label.
    const std::vector<Real> volumes{ 1, 1 };
    const std::vector<Real> moments{ -0.9, +0.9 };
    const std::vector<MinSTCut::Edge> edges{ { 0, 1, 1e6 } };

    // Without pinning, the huge edge forces both labels to whichever
    // side wins on total unary; certainly the cut does not split them.
    const auto unpinned = MinSTCut().classify(volumes, moments, edges);
    EXPECT_EQ(unpinned.labels[0], unpinned.labels[1]);

    MinSTCut::Options pinned;
    pinned.farFieldThreshold = 0.5;
    const auto pinnedResult =
      MinSTCut().classify(volumes, moments, edges, pinned);
    EXPECT_EQ(pinnedResult.labels[0], MinSTCut::Inside);
    EXPECT_EQ(pinnedResult.labels[1], MinSTCut::Outside);
    EXPECT_EQ(pinnedResult.cutEdges.size(), 1u);
  }

  TEST(Rodin_Geometry_MinSTCut, CellInBandFixesFarFieldCells)
  {
    // Three cells in a row. Outer cells are out-of-band (fixed to their
    // moment sign), the middle cell is free with a tiny opposing moment.
    // The fixed cells force the middle to flip along with the smoothing
    // term.
    const std::vector<Real> volumes(3, 1);
    const std::vector<Real> moments{ -1, +1e-3, -1 };
    const auto edges = gridEdges(1, 3, 10.0);

    MinSTCut::Options options;
    options.cellInBand = { false, true, false };
    options.lambdaScale = 1.0;

    const auto result = MinSTCut().classify(volumes, moments, edges, options);
    EXPECT_EQ(result.labels[0], MinSTCut::Inside);
    EXPECT_EQ(result.labels[2], MinSTCut::Inside);
    EXPECT_EQ(result.labels[1], MinSTCut::Inside);
  }

  TEST(Rodin_Geometry_MinSTCut, PerEdgeLambdaSteersCutToWeakestFacets)
  {
    // 1x4 row with a clean sign change between cells 1 and 2.
    // Default pairwise capacities are uniform and small; an additional
    // tiny per-edge lambda on edge (1,2) makes it the cheapest cut, while
    // boosting the other two edges. The cut must land on (1,2).
    const std::vector<Real> volumes(4, 1);
    const std::vector<Real> moments{ -1, -0.4, +0.4, +1 };
    const auto edges = gridEdges(1, 4, 1.0);
    ASSERT_EQ(edges.size(), 3u);

    MinSTCut::Options options;
    options.perEdgeLambda = { 5.0, 1e-3, 5.0 };

    const auto result = MinSTCut().classify(volumes, moments, edges, options);
    EXPECT_EQ(result.labels[0], MinSTCut::Inside);
    EXPECT_EQ(result.labels[1], MinSTCut::Inside);
    EXPECT_EQ(result.labels[2], MinSTCut::Outside);
    EXPECT_EQ(result.labels[3], MinSTCut::Outside);
    ASSERT_EQ(result.cutEdges.size(), 1u);
    EXPECT_EQ(result.cutEdges[0].first, 1u);
    EXPECT_EQ(result.cutEdges[0].second, 2u);
  }

  TEST(Rodin_Geometry_MinSTCut, UnaryScaleShiftsDataVsSmoothingBalance)
  {
    // Two cells with weak moments straddling zero, joined by a small
    // smoothing edge. Tiny unary lets smoothing dominate (same label);
    // large unary lets the data dominate (split).
    const std::vector<Real> volumes(2, 1);
    const std::vector<Real> moments{ -0.1, +0.1 };
    const std::vector<MinSTCut::Edge> edges{ { 0, 1, 1.0 } };

    MinSTCut::Options dataLight;
    dataLight.unaryScale = 1e-3;
    const auto smoothed =
      MinSTCut().classify(volumes, moments, edges, dataLight);
    EXPECT_EQ(smoothed.labels[0], smoothed.labels[1]);

    MinSTCut::Options dataHeavy;
    dataHeavy.unaryScale = 1e3;
    const auto split =
      MinSTCut().classify(volumes, moments, edges, dataHeavy);
    EXPECT_NE(split.labels[0], split.labels[1]);
    EXPECT_EQ(split.labels[0], MinSTCut::Inside);
    EXPECT_EQ(split.labels[1], MinSTCut::Outside);
  }

  TEST(Rodin_Geometry_MinSTCut, MismatchedInputSizesThrow)
  {
    EXPECT_THROW(
        MinSTCut().classify({ 1.0 }, { 1.0, 1.0 }, {}),
        std::invalid_argument);

    MinSTCut::Options options;
    options.perEdgeLambda = { 1.0 };
    EXPECT_THROW(
        MinSTCut().classify({ 1.0, 1.0 }, { 0.0, 0.0 }, {}, options),
        std::invalid_argument);

    MinSTCut::Options bandOptions;
    bandOptions.cellInBand = { true };
    EXPECT_THROW(
        MinSTCut().classify({ 1.0, 1.0 }, { 0.0, 0.0 }, {}, bandOptions),
        std::invalid_argument);
  }

  TEST(Rodin_Geometry_MinSTCut, NegativeCapacitiesThrow)
  {
    EXPECT_THROW(
        MinSTCut().classify(
            { 1.0, 1.0 },
            { -1.0, +1.0 },
            { { 0, 1, -1.0 } }),
        std::invalid_argument);

    MinSTCut::Options options;
    options.perEdgeLambda = { -1.0 };
    EXPECT_THROW(
        MinSTCut().classify(
            { 1.0, 1.0 },
            { -1.0, +1.0 },
            { { 0, 1, 1.0 } },
            options),
        std::invalid_argument);
  }

  TEST(Rodin_Geometry_MinSTCut, OutOfRangeEdgeThrows)
  {
    EXPECT_THROW(
        MinSTCut().classify(
            { 1.0, 1.0 },
            { -1.0, +1.0 },
            { { 0, 5, 1.0 } }),
        std::out_of_range);
  }

  TEST(Rodin_Geometry_MinSTCut, EnergyMatchesLabelChoiceForUnconnectedCells)
  {
    // Without edges the unary contribution equals the sum of the
    // selected-label costs.
    const std::vector<Real> volumes{ 2.0, 3.0 };
    const std::vector<Real> moments{ -0.25, +0.5 };
    const auto result = MinSTCut().classify(volumes, moments, {});
    const Real expected =
      MinSTCut::getInsideCost(volumes[0], moments[0])
      + MinSTCut::getOutsideCost(volumes[1], moments[1]);
    EXPECT_DOUBLE_EQ(result.energy, expected);
  }
}
