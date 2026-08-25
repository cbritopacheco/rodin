/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */

/**
 * @file
 * @brief Development driver emitting inlinable quadrature coefficients.
 *
 * Not a test. Degrees are independent, so they are searched in parallel, and
 * each is flushed as it completes so an interrupted run keeps what it found.
 */
#include <cstdio>
#include <string>
#include <mutex>
#include <thread>
#include <atomic>
#include <vector>
#include <algorithm>

#include "PublishedCounts.h"
#include "QuadratureInvariants.h"
#include "Rodin/QF/NodeElimination.h"

using namespace Rodin;
using namespace Rodin::QF;

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
  const size_t d = Geometry::Polytope::Traits(g).getDimension();

  std::atomic<size_t> next{1};
  std::vector<std::thread> pool;
  const size_t jobs = (argc > 3) ? std::stoul(argv[3]) : 4;
  const size_t W = std::min<size_t>(jobs, maxDegree);
  for (size_t w = 0; w < W; ++w)
    pool.emplace_back([&] {
      for (;;)
      {
        const size_t p = next++;
        if (p > maxDegree)
          return;
        const auto seed = NodeElimination::productSeed(g, p);
        const auto out = NodeElimination::reduce(g, p, seed);

        // A table is worth having only if what it holds is exact, so each
        // rule is measured against the moment oracle before it is printed.
        const auto report = Rodin::Tests::QF::exactnessSweep(out, g, p);
        const bool positive = Rodin::Tests::QF::allWeightsPositive(out);
        const bool inside = Rodin::Tests::QF::allPointsInside(out, g, 1e-12);
        const std::lock_guard<std::mutex> lock(g_mutex);
        if (report.worstRelativeError >= 1e-12 || !positive || !inside)
        {
          std::printf("    // degree %zu: REJECTED (oracle %.2e%s%s)\n", p,
            report.worstRelativeError, positive ? "" : ", negative weight",
            inside ? "" : ", point outside");
          std::fflush(stdout);
          continue;
        }

        std::printf("    // degree %zu: %zu points, oracle %.2e\n", p,
          out.getSize(), report.worstRelativeError);
        std::printf("      {\n");
        for (size_t q = 0; q < out.getSize(); ++q)
        {
          std::printf("        ");
          for (size_t k = 0; k < d; ++k)
            std::printf("%.17g, ", out.getPoint(q)[static_cast<Eigen::Index>(k)]);
          std::printf("%.17g,\n", out.getWeight(q));
        }
        std::printf("      },\n");
        std::fflush(stdout);
      }
    });
  for (auto& t : pool)
    t.join();
  return 0;
}
