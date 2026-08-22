#include <gtest/gtest.h>
#include <chrono>
#include <thread>
#include <atomic>
#include <mutex>
#include "Rodin/QF/NodeElimination.h"
#include "QuadratureInvariants.h"
using namespace Rodin; using namespace Rodin::QF;
using namespace Rodin::Geometry; using namespace Rodin::Tests::QF;
namespace { std::mutex mtx; }
TEST(Sweep, ToTwenty)
{
  struct C { const char* n; Polytope::Type g; size_t maxd; };
  for (const auto& c : std::vector<C>{{"tri", Polytope::Type::Triangle, 20},
                                      {"tet", Polytope::Type::Tetrahedron, 12}})
  {
    std::atomic<size_t> next{1};
    std::vector<std::thread> pool;
    const size_t W = std::min<size_t>(std::thread::hardware_concurrency(), c.maxd);
    for (size_t w = 0; w < W; ++w)
      pool.emplace_back([&]{
        for (;;) {
          const size_t p = next++;
          if (p > c.maxd) return;
          const auto seed = NodeElimination::productSeed(c.g, p);
          const auto t0 = std::chrono::steady_clock::now();
          const auto out = NodeElimination::reduce(c.g, p, seed);
          const double s = std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
          const auto rep = exactnessSweep(out, c.g, p);
          const bool ok = allWeightsPositive(out) && allPointsInside(out, c.g);
          std::lock_guard<std::mutex> lk(mtx);
          std::printf("%s d%-2zu %4zu -> %4zu  err %.1e  ok %d  %.0fs\n",
            c.n, p, seed.getSize(), out.getSize(), rep.worstRelativeError, (int)ok, s);
          std::fflush(stdout);
        }
      });
    for (auto& t : pool) t.join();
  }
}
