#include <cstdio>
#include <string>
#include "PublishedCounts.h"
#include "QuadratureInvariants.h"
#include "Rodin/QF/SymmetricRuleGenerator.h"
using namespace Rodin; using namespace Rodin::QF; using namespace Rodin::Tests::QF;
int main(int argc, char** argv){
  const std::string w = (argc > 1) ? argv[1] : "tri";
  const auto g = (w=="tet") ? Geometry::Polytope::Type::Tetrahedron
    : (w=="wedge") ? Geometry::Polytope::Type::Wedge
    : (w=="pyr") ? Geometry::Polytope::Type::Pyramid
    : (w=="quad") ? Geometry::Polytope::Type::Quadrilateral
    : (w=="hex") ? Geometry::Polytope::Type::Hexahedron
    : Geometry::Polytope::Type::Triangle;
  const size_t maxDeg = (argc>2) ? std::stoul(argv[2]) : 10;
  const size_t maxPts = (argc>3) ? std::stoul(argv[3]) : 120;
  const size_t lo = (argc>4) ? std::stoul(argv[4]) : 1;
  const double budget = (argc>5) ? std::stod(argv[5]) : 0.0;
  const size_t d = Geometry::Polytope::Traits(g).getDimension();
  std::printf("deg  ours  WV    XG    bound cond  status oracle\n");
  for (size_t p = lo; p <= maxDeg; ++p) {
    const auto r = SymmetricRuleGenerator::search(g, p, maxPts, 64, 1e-12, 0, budget);
    const bool ok = r.converged && r.admissible;
    const auto rep = exactnessSweep(r, g, ok ? p : 0);
    const size_t wv = publishedCount(witherdenVincentCounts(g), p);
    const size_t xg = publishedCount(xiaoGimbutasCounts(g), p);
    std::printf("%-4zu %-5s %-5s %-5s %-5zu %-5zu %-6s %.2e%s%s\n", p,
      ok ? std::to_string(r.getSize()).c_str() : "-",
      wv ? std::to_string(wv).c_str() : ".",
      xg ? std::to_string(xg).c_str() : ".", countingBound(d, p),
      SymmetricRuleGenerator::invariantDimension(g, p),
      ok ? "ok" : "FAIL", ok ? rep.worstRelativeError : r.residual,
      (ok && allWeightsPositive(r)) ? "" : " NEGW",
      (ok && allPointsInside(r, g, 1e-12)) ? "" : " OUTSIDE");
    std::fflush(stdout);
  }
  return 0;
}
