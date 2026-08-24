#include <gtest/gtest.h>
#include <chrono>
#include "Rodin/QF/SymmetricRuleSolver.h"
#include "QuadratureInvariants.h"
using namespace Rodin;
using namespace Rodin::QF;
using namespace Rodin::Geometry;
using namespace Rodin::Tests::QF;
namespace
{
  struct R
  {
      const std::vector<SymmetricOrbit>& o;
      size_t getSize() const
      {
        return o.size();
      }
      Real getWeight(size_t i) const
      {
        return o[i].getWeight();
      }
      Math::SpatialVector<Real> getPoint(size_t i) const
      {
        Math::SpatialVector<Real> p;
        p.resize(3);
        for (int k = 0; k < 3; ++k)
          p[k] = o[i].getBarycentric()[k];
        return p;
      }
  };
}
TEST(Diag, PyramidRules)
{
  const int wv[] = {1, 5, 6, 10, 15, 24, 31, 44};
  for (size_t p = 1; p <= 4; ++p)
  {
    const auto t0 = std::chrono::steady_clock::now();
    SymmetricRuleSolver::PyramidConfiguration cfg;
    const auto r = SymmetricRuleSolver::searchPyramid(p, 24, 120, &cfg);
    const double s =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (!(r.converged && r.admissible))
    {
      std::printf("pyr d%zu FAILED res %.2e (%.0fs)\n", p, r.residual, s);
      std::fflush(stdout);
      continue;
    }
    const R rule{r.orbits};
    const auto rep = exactnessSweep(rule, Polytope::Type::Pyramid, p);
    std::printf("pyr d%zu -> %2zu pts (WV %d) err %.1e ok %d (%.0fs)\n", p,
      r.orbits.size(), wv[p - 1], rep.worstRelativeError,
      (int)(allWeightsPositive(rule) && allPointsInside(rule, Polytope::Type::Pyramid)),
      s);
    std::fflush(stdout);
  }
}
