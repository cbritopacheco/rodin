/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
/**
 * @file
 * @brief Development driver emitting symmetric rule coefficients.
 */
#include <cstdio>
#include <string>
#include <mutex>
#include <thread>
#include <atomic>
#include <algorithm>
#include "Rodin/QF/SymmetricRuleSolver.h"
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
    : (which == "wedge")          ? Geometry::Polytope::Type::Wedge
                                  : Geometry::Polytope::Type::Triangle;
  const size_t maxDegree = (argc > 2) ? std::stoul(argv[2]) : 10;
  const size_t maxPoints = (argc > 3) ? std::stoul(argv[3]) : 60;
  const size_t d = Geometry::Polytope::Traits(g).getDimension();

  std::atomic<size_t> next{1};
  std::vector<std::thread> pool;
  const size_t W = std::min<size_t>(std::thread::hardware_concurrency(), maxDegree);
  for (size_t w = 0; w < W; ++w)
    pool.emplace_back([&] {
      for (;;)
      {
        const size_t p = next++;
        if (p > maxDegree)
          return;
        SymmetricRuleSolver::Configuration cfg;
        const auto r = SymmetricRuleSolver::search(g, p, maxPoints, 256, &cfg);
        std::lock_guard<std::mutex> lk(g_mutex);
        if (!(r.converged && r.admissible))
        {
          std::printf("  // degree %zu: FAILED res %.2e\n", p, r.residual);
          std::fflush(stdout);
          continue;
        }
        size_t npts = 0;
        for (const auto& o : r.orbits)
          npts += o.getSize();
        std::printf("  // degree %zu: %zu points\n  {", p, npts);
        char buf[128];
        for (const auto& o : r.orbits)
          for (const auto& pt : o.expandPoints(g))
          {
            for (size_t k = 0; k < d; ++k)
            {
              std::snprintf(
                buf, sizeof(buf), "%.17g, ", pt[static_cast<Eigen::Index>(k)]);
              std::printf("%s", buf);
            }
            std::snprintf(buf, sizeof(buf), "%.17g, ", o.getWeight());
            std::printf("%s", buf);
          }
        std::printf("},\n");
        std::fflush(stdout);
      }
    });
  for (auto& t : pool)
    t.join();
  return 0;
}
