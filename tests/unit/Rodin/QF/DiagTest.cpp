#include <cstdio>
#include "Rodin/QF/SymmetricRuleGenerator.h"
using namespace Rodin; using namespace Rodin::QF;
int main(int argc, char** argv){
  const auto g = Geometry::Polytope::Type::Wedge;
  const size_t degree = (argc>1)? std::stoul(argv[1]) : 7;
  const size_t points = (argc>2)? std::stoul(argv[2]) : 35;
  const size_t maxOrbits = (argc>3)? std::stoul(argv[3]) : 7;
  const size_t restarts = (argc>4)? std::stoul(argv[4]) : 128;
  const auto& st = SymmetryGroup::strata(g);
  const size_t cond = SymmetricRuleGenerator::invariantDimension(g, degree);
  auto cands = SymmetricRuleGenerator::decompositions(g, points);
  std::printf("wedge d%zu at %zu pts: %zu candidates, %zu conditions\n",
              degree, points, cands.size(), cond);
  size_t tried = 0;
  for (const auto& c : cands) {
    if (c.size() > maxOrbits) continue;
    size_t unk = c.size(); bool bnd = false;
    for (size_t w : c) { unk += st[w].getFreedom(); bnd |= st[w].boundary; }
    if (bnd) continue;
    const auto r = SymmetricRuleGenerator::solve(g, degree, c, restarts);
    ++tried;
    if (r.converged && r.admissible) {
      std::printf("  SOLVED [");
      for (size_t w : c) std::printf("%zu ", st[w].orbitSize);
      std::printf("] unknowns %zu  residual %.2e\n", unk, r.residual);
      std::fflush(stdout);
    }
  }
  std::printf("  tried %zu interior decompositions of <= %zu orbits\n", tried, maxOrbits);
  return 0;
}
