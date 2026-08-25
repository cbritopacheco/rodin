/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
/**
 * @file
 * @brief Development driver writing the XiaoGimbutas coefficient tables.
 *
 * Not a test. Prints one array of XiaoGimbutas.cpp, as pasteable C++, for the
 * named element. Degrees are independent and searched in parallel, each
 * flushed as it completes so an interrupted run keeps what it found.
 *
 * @par Usage
 * @code
 *   RodinXiaoGimbutas <element> [maxDegree] [jobs]
 * @endcode
 * where @c element is tri or tet:
 * @code
 *   RodinXiaoGimbutas tri > triangle.txt      # s_triangle
 *   RodinXiaoGimbutas tet > tetrahedron.txt   # s_tetrahedron
 * @endcode
 *
 * These rules come from node elimination: a larger rule has nodes removed one
 * at a time, which destroys symmetry but reaches lower counts than a symmetric
 * construction can. That is why the family covers only the simplices here, and
 * why its published counts differ from those of WitherdenVincent.cpp for the
 * same element and strength --- the two are compared only against their own
 * column. The symmetric rules, covering every element, are written by
 * RodinWitherdenVincent.
 *
 * Most strengths reached here are below the published count, and several sit
 * exactly on the counting bound, which no rule of that strength can beat.
 */
#include <atomic>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "PublishedCounts.h"
#include "TableEmitter.h"
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
  if (argc < 2)
  {
    std::fprintf(stderr, "usage: RodinXiaoGimbutas <tri|tet> [maxDegree] [jobs]\n");
    return 2;
  }
  const auto g = elementOf(argv[1]);
  if (g != Geometry::Polytope::Type::Triangle &&
    g != Geometry::Polytope::Type::Tetrahedron)
  {
    std::fprintf(stderr,
      "node elimination here covers only the simplices, which is where Xiao"
      " and Gimbutas publish; use RodinWitherdenVincent for '%s'\n",
      argv[1]);
    return 2;
  }
  const auto& published = xiaoGimbutasCounts(g);
  const size_t maxDegree = (argc > 2)
    ? std::stoul(argv[2])
    : (published.empty() ? 15 : published.rbegin()->first);
  const size_t jobs =
    (argc > 3) ? std::stoul(argv[3]) : std::thread::hardware_concurrency();

  std::atomic<size_t> next{1};
  std::vector<std::thread> pool;
  for (size_t w = 0; w < std::max<size_t>(1, std::min(jobs, maxDegree)); ++w)
    pool.emplace_back([&] {
      for (;;)
      {
        const size_t degree = next++;
        if (degree > maxDegree)
          return;
        const auto seed = NodeElimination::productSeed(g, degree);
        const auto rule = NodeElimination::reduce(g, degree, seed);
        const std::lock_guard<std::mutex> lock(g_mutex);
        emitEntry(rule, g, degree);
      }
    });
  for (auto& worker : pool)
    worker.join();
  return 0;
}
