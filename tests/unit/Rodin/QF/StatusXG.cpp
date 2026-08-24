/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
/**
 * @file
 * @brief Development driver reporting node-elimination counts against the
 * published Xiao-Gimbutas tables.
 *
 * Not a test. Degrees are independent and searched in parallel, each flushed
 * as it completes so an interrupted run keeps what it found.
 */
#include <atomic>
#include <cstdio>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "PublishedCounts.h"
#include "QuadratureInvariants.h"
#include "Rodin/QF/NodeElimination.h"

using namespace Rodin;
using namespace Rodin::QF;
using namespace Rodin::Tests::QF;

namespace
{
  std::mutex g_mutex;
}

int main(int argc, char** argv)
{
  const std::string which = (argc > 1) ? argv[1] : "tri";
  const auto g = (which == "tet") ? Geometry::Polytope::Type::Tetrahedron
                                  : Geometry::Polytope::Type::Triangle;
  const size_t maxDegree = (argc > 2) ? std::stoul(argv[2]) : 20;
  const size_t jobs = (argc > 3) ? std::stoul(argv[3]) : 8;
  const size_t d = Geometry::Polytope::Traits(g).getDimension();

  std::printf("deg  ours  XG    bound status oracle\n");
  std::fflush(stdout);

  std::atomic<size_t> next{1};
  std::vector<std::thread> pool;
  const size_t workers = std::min<size_t>(jobs, maxDegree);
  for (size_t w = 0; w < workers; ++w)
    pool.emplace_back([&] {
      for (;;)
      {
        const size_t p = next++;
        if (p > maxDegree)
          return;
        const auto seed = NodeElimination::productSeed(g, p);
        const auto rule = NodeElimination::reduce(g, p, seed);
        const auto report = exactnessSweep(rule, g, p);
        const bool positive = allWeightsPositive(rule);
        const bool inside = allPointsInside(rule, g, 1e-12);
        const bool ok = (report.worstRelativeError < 1e-12) && positive && inside;
        const size_t xg = publishedCount(xiaoGimbutasCounts(g), p);

        const std::lock_guard<std::mutex> lock(g_mutex);
        std::printf("%-4zu %-5zu %-5s %-5zu %-6s %.2e%s%s\n", p, rule.getSize(),
          xg ? std::to_string(xg).c_str() : ".", countingBound(d, p), ok ? "ok" : "FAIL",
          report.worstRelativeError, positive ? "" : " NEGW", inside ? "" : " OUTSIDE");
        std::fflush(stdout);
      }
    });
  for (auto& worker : pool)
    worker.join();
  return 0;
}
