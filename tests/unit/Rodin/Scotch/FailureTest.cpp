/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <cstdio>
#include <algorithm>
#include <array>
#include <memory>
#include <string>
#include <vector>
#include "Rodin/Scotch/MeshPartitioner.h"

using namespace Rodin;
using namespace Rodin::Geometry;

namespace
{
  enum class Operation
  {
    GraphInit,
    GraphBuild,
    StrategyInit,
    StrategyMap,
    Partition,
    TemporaryFile,
    StrategySave,
    Rewind,
    Read
  };

  struct Failure
  {
    Failure()
      : failedCalls(0),
        graphsClosed(0),
        strategiesClosed(0),
        filesClosed(0),
        mappingCalls(0),
        readCalls(0)
    {}
    Optional<Operation> operation;
    Optional<FILE*> stream;
    size_t failedCalls, graphsClosed, strategiesClosed, filesClosed;
    size_t mappingCalls, readCalls;
    std::string mapping;

    void arm(Operation selected)
    {
      operation = selected;
      failedCalls = graphsClosed = strategiesClosed = filesClosed = readCalls = 0;
    }

    bool reject(Operation selected)
    {
      if (operation != selected)
        return false;
      ++failedCalls;
      return true;
    }
  };
  Failure failure;
}

extern "C"
{
  int __real_SCOTCH_graphInit(SCOTCH_Graph* graph);
  int __wrap_SCOTCH_graphInit(SCOTCH_Graph* graph)
  {
    return failure.reject(Operation::GraphInit) ? 1 : __real_SCOTCH_graphInit(graph);
  }

  int __real_SCOTCH_graphBuild(SCOTCH_Graph* graph, SCOTCH_Num base, SCOTCH_Num vertices,
    const SCOTCH_Num* starts, const SCOTCH_Num* ends, const SCOTCH_Num* weights,
    const SCOTCH_Num* labels, SCOTCH_Num edges, const SCOTCH_Num* indices,
    const SCOTCH_Num* loads);
  int __wrap_SCOTCH_graphBuild(SCOTCH_Graph* graph, SCOTCH_Num base, SCOTCH_Num vertices,
    const SCOTCH_Num* starts, const SCOTCH_Num* ends, const SCOTCH_Num* weights,
    const SCOTCH_Num* labels, SCOTCH_Num edges, const SCOTCH_Num* indices,
    const SCOTCH_Num* loads)
  {
    return failure.reject(Operation::GraphBuild)
      ? 1
      : __real_SCOTCH_graphBuild(
          graph, base, vertices, starts, ends, weights, labels, edges, indices, loads);
  }

  int __real_SCOTCH_stratInit(SCOTCH_Strat* strategy);
  int __wrap_SCOTCH_stratInit(SCOTCH_Strat* strategy)
  {
    return failure.reject(Operation::StrategyInit) ? 1
                                                   : __real_SCOTCH_stratInit(strategy);
  }

  int __real_SCOTCH_stratGraphMap(SCOTCH_Strat* strategy, const char* text);
  int __wrap_SCOTCH_stratGraphMap(SCOTCH_Strat* strategy, const char* text)
  {
    ++failure.mappingCalls;
    failure.mapping = text;
    return failure.reject(Operation::StrategyMap)
      ? 1
      : __real_SCOTCH_stratGraphMap(strategy, text);
  }

  int __real_SCOTCH_graphPart(SCOTCH_Graph* graph, SCOTCH_Num count,
    SCOTCH_Strat* strategy, SCOTCH_Num* assignment);
  int __wrap_SCOTCH_graphPart(
    SCOTCH_Graph* graph, SCOTCH_Num count, SCOTCH_Strat* strategy, SCOTCH_Num* assignment)
  {
    return failure.reject(Operation::Partition)
      ? 1
      : __real_SCOTCH_graphPart(graph, count, strategy, assignment);
  }

  int __real_SCOTCH_stratSave(const SCOTCH_Strat* strategy, FILE* stream);
  int __wrap_SCOTCH_stratSave(const SCOTCH_Strat* strategy, FILE* stream)
  {
    return failure.reject(Operation::StrategySave)
      ? 1
      : __real_SCOTCH_stratSave(strategy, stream);
  }

  void __real_SCOTCH_graphExit(SCOTCH_Graph* graph);
  void __wrap_SCOTCH_graphExit(SCOTCH_Graph* graph)
  {
    if (failure.operation)
      ++failure.graphsClosed;
    __real_SCOTCH_graphExit(graph);
  }

  void __real_SCOTCH_stratExit(SCOTCH_Strat* strategy);
  void __wrap_SCOTCH_stratExit(SCOTCH_Strat* strategy)
  {
    if (failure.operation)
      ++failure.strategiesClosed;
    __real_SCOTCH_stratExit(strategy);
  }

  FILE* __real_tmpfile();
  FILE* __wrap_tmpfile()
  {
    if (failure.reject(Operation::TemporaryFile))
      return nullptr;
    auto* stream = __real_tmpfile();
    if (failure.operation && stream)
      failure.stream = stream;
    return stream;
  }

  int __real_fclose(FILE* stream);
  int __wrap_fclose(FILE* stream)
  {
    if (failure.stream && *failure.stream == stream)
    {
      ++failure.filesClosed;
      failure.stream.reset();
    }
    return __real_fclose(stream);
  }

  int __real_fseek(FILE* stream, long offset, int origin);
  int __wrap_fseek(FILE* stream, long offset, int origin)
  {
    return failure.reject(Operation::Rewind) ? -1 : __real_fseek(stream, offset, origin);
  }

  size_t __real_fread(void* output, size_t size, size_t count, FILE* stream);
  size_t __wrap_fread(void* output, size_t size, size_t count, FILE* stream)
  {
    if (failure.operation == Operation::Read && failure.readCalls++ == 0)
      return __real_fread(output, size, std::min(count, size_t{3}), stream);
    return failure.reject(Operation::Read) ? 0
                                           : __real_fread(output, size, count, stream);
  }

  int __real_ferror(FILE* stream);
  int __wrap_ferror(FILE* stream)
  {
    return failure.operation == Operation::Read ? 1 : __real_ferror(stream);
  }
}

namespace Rodin::Tests::Unit
{
  class ScotchFailureTest : public ::testing::Test
  {
    protected:
      void TearDown() override
      {
        failure.operation.reset();
        EXPECT_FALSE(failure.stream);
      }

      template <class Action>
      void expectFailure(
        Operation operation, Action&& action, const char* method, const char* diagnostic)
      {
        failure.arm(operation);
        try
        {
          action();
          ADD_FAILURE() << "The selected backend failure was accepted";
        }
        catch (const Alert::Exception& error)
        {
          const std::string message = error.what();
          EXPECT_NE(
            message.find(std::string("member function ") + method), std::string::npos);
          EXPECT_NE(message.find(diagnostic), std::string::npos);
        }
        failure.operation.reset();
        EXPECT_EQ(failure.failedCalls, 1u);
      }

      std::vector<size_t> assignments(const Scotch::Partitioner& partitioner)
      {
        std::vector<size_t> result;
        for (Index cell = 0; cell < partitioner.getMesh().getCellCount(); ++cell)
          result.push_back(partitioner.getPartition(cell));
        return result;
      }

      void expectValid(const Scotch::Partitioner& partitioner)
      {
        EXPECT_EQ(partitioner.getCount(), 2u);
        std::array<size_t, 2> counts{0, 0};
        for (const auto part : assignments(partitioner))
        {
          ASSERT_LT(part, counts.size());
          ++counts[part];
        }
        EXPECT_GT(counts[0], 0u);
        EXPECT_GT(counts[1], 0u);
      }
  };

  /// @brief Backend graph failures release only initialized resources and preserve prior assignments.
  TEST_F(ScotchFailureTest, PartitionFailuresPreserveAssignmentsAndPermitRetry)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {9});
    mesh.getConnectivity().compute(1, 1);
    Scotch::Partitioner partitioner(mesh);
    SCOTCH_Strat strategy;
    ASSERT_EQ(SCOTCH_stratInit(&strategy), 0);
    std::unique_ptr<SCOTCH_Strat, decltype(&SCOTCH_stratExit)> owner(
      &strategy, &SCOTCH_stratExit);
    ASSERT_EQ(SCOTCH_stratGraphMapBuild(&strategy, SCOTCH_STRATQUALITY, 2, 0.05), 0);
    partitioner.setStrategy(strategy).partition(2);
    for (const auto operation : {Operation::GraphInit, Operation::GraphBuild,
           Operation::StrategyInit, Operation::StrategyMap, Operation::Partition})
    {
      SCOPED_TRACE(static_cast<int>(operation));
      const auto previous = assignments(partitioner);
      const auto mapping = failure.mapping;
      expectFailure(operation, [&] { partitioner.partition(3); }, "partition", "failed");
      EXPECT_EQ(partitioner.getCount(), 2u);
      EXPECT_EQ(assignments(partitioner), previous);
      EXPECT_EQ(failure.graphsClosed, operation == Operation::GraphInit ? 0u : 1u);
      EXPECT_EQ(failure.strategiesClosed,
        operation == Operation::StrategyMap || operation == Operation::Partition ? 1u
                                                                                 : 0u);
      const auto previousCalls = failure.mappingCalls;
      partitioner.partition(2);
      EXPECT_EQ(failure.mappingCalls, previousCalls + 1);
      expectValid(partitioner);
      EXPECT_EQ(failure.mapping, mapping);
    }
  }

  /// @brief File and serialization errors preserve the owned strategy and close temporary streams.
  TEST_F(ScotchFailureTest, StrategyFailuresPreserveConfigurationAndPermitRetry)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {9});
    mesh.getConnectivity().compute(1, 1);
    Scotch::Partitioner partitioner(mesh);
    SCOTCH_Strat strategy;
    ASSERT_EQ(SCOTCH_stratInit(&strategy), 0);
    std::unique_ptr<SCOTCH_Strat, decltype(&SCOTCH_stratExit)> owner(
      &strategy, &SCOTCH_stratExit);
    ASSERT_EQ(SCOTCH_stratGraphMapBuild(&strategy, SCOTCH_STRATQUALITY, 2, 0.05), 0);
    partitioner.setStrategy(strategy).partition(2);
    const auto mapping = failure.mapping;
    SCOTCH_Strat replacement;
    ASSERT_EQ(SCOTCH_stratInit(&replacement), 0);
    std::unique_ptr<SCOTCH_Strat, decltype(&SCOTCH_stratExit)> replacementOwner(
      &replacement, &SCOTCH_stratExit);
    ASSERT_EQ(SCOTCH_stratGraphMapBuild(&replacement, SCOTCH_STRATSPEED, 2, 0.1), 0);
    for (const auto operation : {Operation::TemporaryFile, Operation::StrategySave,
           Operation::Rewind, Operation::Read})
    {
      SCOPED_TRACE(static_cast<int>(operation));
      const auto previous = assignments(partitioner);
      expectFailure(
        operation, [&] { partitioner.setStrategy(replacement); }, "setStrategy",
        operation == Operation::Read ? "strategy read failed"
                                     : "strategy serialization failed");
      EXPECT_EQ(partitioner.getCount(), 2u);
      EXPECT_EQ(assignments(partitioner), previous);
      EXPECT_EQ(failure.filesClosed, operation == Operation::TemporaryFile ? 0u : 1u);
      EXPECT_FALSE(failure.stream);
      const auto previousCalls = failure.mappingCalls;
      partitioner.partition(2);
      EXPECT_EQ(failure.mappingCalls, previousCalls + 1);
      expectValid(partitioner);
      EXPECT_EQ(failure.mapping, mapping);
    }
    partitioner.setStrategy(replacement).partition(2);
    expectValid(partitioner);
    EXPECT_NE(failure.mapping, mapping);
  }
}
