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

#include "Rodin/QF/NodeElimination.h"

using namespace Rodin;
using namespace Rodin::QF;

namespace { std::mutex g_mutex; }

int main(int argc, char** argv)
{
  const std::string which = (argc > 1) ? argv[1] : "tri";
  const auto g = (which == "tet") ? Geometry::Polytope::Type::Tetrahedron
                                  : Geometry::Polytope::Type::Triangle;
  const size_t maxDegree = (argc > 2) ? std::stoul(argv[2]) : 20;
  const size_t d = Geometry::Polytope::Traits(g).getDimension();

  std::atomic<size_t> next{1};
  std::vector<std::thread> pool;
  const size_t W = std::min<size_t>(std::thread::hardware_concurrency(), maxDegree);
  for (size_t w = 0; w < W; ++w)
    pool.emplace_back([&]{
      for (;;)
      {
        const size_t p = next++;
        if (p > maxDegree) return;
        const auto seed = NodeElimination::productSeed(g, p);
        const auto out = NodeElimination::reduce(g, p, seed);
        std::string line = "  // degree " + std::to_string(p) + ": "
          + std::to_string(out.getSize()) + " points\n  {";
        char buf[128];
        for (size_t q = 0; q < out.getSize(); ++q)
        {
          for (size_t k = 0; k < d; ++k)
          {
            std::snprintf(buf, sizeof(buf), "%.17g, ",
              out.getPoint(q)[static_cast<Eigen::Index>(k)]);
            line += buf;
          }
          std::snprintf(buf, sizeof(buf), "%.17g,%s", out.getWeight(q),
            (q + 1 < out.getSize()) ? " " : "");
          line += buf;
        }
        line += "},";
        std::lock_guard<std::mutex> lk(g_mutex);
        std::printf("%s\n", line.c_str());
        std::fflush(stdout);
      }
    });
  for (auto& t : pool) t.join();
  return 0;
}
