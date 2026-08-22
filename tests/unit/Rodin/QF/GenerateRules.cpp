/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Development driver that generates symmetric quadrature rules for
 * inlining.
 *
 * Not a test. Rules are searched offline and their coefficients inlined; this
 * is what produces them.
 *
 * Degrees are mutually independent, so they are searched in parallel, one
 * worker per degree. Results are appended to the output file the moment a
 * degree completes rather than at the end, so a run that is interrupted still
 * yields everything it had found. Determinism is unaffected: each degree is a
 * separate deterministic search and the file is sorted on read, not on write.
 *
 * Usage: RodinGenerateRules <tri|tet|wedge> <maxDegree> [maxPoints] [restarts]
 */

#include <cstdio>
#include <string>
#include <mutex>
#include <thread>
#include <atomic>
#include <vector>
#include <algorithm>
#include <chrono>

#include "Rodin/QF/SymmetricRuleSolver.h"

using namespace Rodin;
using namespace Rodin::QF;

namespace
{
  std::mutex g_mutex;

  std::string format(Geometry::Polytope::Type g, size_t degree,
    const SymmetricRuleSolver::Configuration& config,
    const SymmetricRuleSolver::Result& r, double seconds)
  {
    char buf[256];
    std::string out;
    std::snprintf(buf, sizeof(buf), "d%zu n=%zu res=%.3e t=%.1fs :", degree,
      SymmetricRuleSolver::configurationSize(config), r.residual, seconds);
    out += buf;
    for (const auto& orbit : r.orbits)
    {
      auto b = orbit.getBarycentric();
      std::sort(b.begin(), b.end());
      out += " {";
      for (size_t i = 0; i < b.size(); ++i)
      {
        std::snprintf(buf, sizeof(buf), "%.17g%s", b[i], (i + 1 < b.size()) ? "," : "");
        out += buf;
      }
      for (const auto t : orbit.getTensor())
      {
        std::snprintf(buf, sizeof(buf), "|z=%.17g", t);
        out += buf;
      }
      std::snprintf(buf, sizeof(buf), "|w=%.17g}", orbit.getWeight());
      out += buf;
    }
    (void)g;
    return out;
  }
}

int main(int argc, char** argv)
{
  const std::string which = (argc > 1) ? argv[1] : "tri";
  const auto g = (which == "tet") ? Geometry::Polytope::Type::Tetrahedron
    : (which == "wedge")          ? Geometry::Polytope::Type::Wedge
                                  : Geometry::Polytope::Type::Triangle;
  const size_t maxDegree = (argc > 2) ? std::stoul(argv[2]) : 8;
  const size_t maxPoints = (argc > 3) ? std::stoul(argv[3]) : 96;
  const size_t restarts = (argc > 4) ? std::stoul(argv[4]) : 256;

  const size_t workers = std::min<size_t>(
    std::max<size_t>(std::thread::hardware_concurrency(), 1u), maxDegree);

  std::atomic<size_t> next{1};
  std::vector<std::thread> pool;
  for (size_t w = 0; w < workers; ++w)
  {
    pool.emplace_back([&] {
      for (;;)
      {
        const size_t degree = next++;
        if (degree > maxDegree)
          return;

        const auto t0 = std::chrono::steady_clock::now();
        SymmetricRuleSolver::Configuration config;
        const auto r =
          SymmetricRuleSolver::search(g, degree, maxPoints, restarts, &config);
        const double seconds =
          std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

        std::string line;
        if (r.converged && r.admissible)
        {
          line = format(g, degree, config, r, seconds);
        }
        else
        {
          char buf[192];
          std::snprintf(buf, sizeof(buf), "d%zu FAILED conv=%d adm=%d res=%.3e t=%.1fs",
            degree, static_cast<int>(r.converged), static_cast<int>(r.admissible),
            r.residual, seconds);
          line = buf;
        }

        // Flush per degree: an interrupted run keeps what it has found.
        std::lock_guard<std::mutex> lock(g_mutex);
        std::printf("%s\n", line.c_str());
        std::fflush(stdout);
      }
    });
  }
  for (auto& t : pool)
    t.join();
  return 0;
}
