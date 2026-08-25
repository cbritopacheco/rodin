/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
/**
 * @file
 * @brief Development driver emitting inlinable symmetric rule coefficients.
 *
 * Not a test. Prints the coefficient tables of WitherdenVincent.cpp, ready to
 * paste, for one element over a range of strengths.
 *
 * @par Usage
 * @code
 *   RodinGenerateTables <element> [maxDegree] [seconds] [fromDegree] [jobs]
 * @endcode
 * where @c element is one of tri, quad, tet, wedge, pyr, hex. Only the element
 * is required; the defaults reproduce the shipped tables:
 *
 * @code
 *   RodinGenerateTables tri          # the whole triangle table
 *   RodinGenerateTables wedge 10     # the prism, strengths one to ten
 * @endcode
 *
 * @par Parameters that matter
 * - @c seconds is a per-strength deadline, and is the only parameter that can
 *   change what comes out: it truncates the search, so a strength that would
 *   have been solved can come back larger or not at all. Zero, the default,
 *   means no deadline. Set it only for an unattended sweep that must cover
 *   ground rather than finish a particular strength.
 * - @c jobs defaults to the machine. It changes how long a search takes and
 *   not what it returns, since the ordering of decompositions and the seeding
 *   are fixed. It matters only through the deadline: the hexahedron at
 *   strengths six to nine was thought unreachable when run on two threads
 *   against a fifteen minute deadline, and takes about two minutes given the
 *   machine.
 *
 * @par What to expect
 * The two-dimensional elements are slow per decomposition and fast to
 * enumerate: the triangle at strength twenty spends about ten seconds on each
 * of 756 candidates. The three-dimensional elements are the reverse --- the
 * prism at strength six enumerates 30 806 candidates at about a tenth of a
 * second each, and its published 28-point rule lies well down that list. Both
 * finish; they just fail differently when hurried.
 *
 * Every rule is checked against the moment oracle before it is printed, so a
 * rule that is not exact, positive and interior never reaches a table.
 */
#include <cstdio>
#include <string>
#include <vector>

#include "PublishedCounts.h"
#include "QuadratureInvariants.h"
#include "Rodin/QF/SymmetricRuleGenerator.h"

using namespace Rodin;
using namespace Rodin::QF;
using namespace Rodin::Tests::QF;

namespace
{
  /// @brief Strengths published for each element, which is how far the tables
  /// are worth taking by default.
  size_t defaultMaxDegree(Geometry::Polytope::Type g)
  {
    const auto& published = witherdenVincentCounts(g);
    return published.empty() ? 10 : published.rbegin()->first;
  }
}

int main(int argc, char** argv)
{
  const std::string which = (argc > 1) ? argv[1] : "tri";
  const auto g = (which == "tet") ? Geometry::Polytope::Type::Tetrahedron
    : (which == "wedge")          ? Geometry::Polytope::Type::Wedge
    : (which == "pyr")            ? Geometry::Polytope::Type::Pyramid
    : (which == "quad")           ? Geometry::Polytope::Type::Quadrilateral
    : (which == "hex")            ? Geometry::Polytope::Type::Hexahedron
                                  : Geometry::Polytope::Type::Triangle;
  const size_t maxDegree = (argc > 2) ? std::stoul(argv[2]) : defaultMaxDegree(g);
  const double seconds = (argc > 3) ? std::stod(argv[3]) : 0.0;
  const size_t from = (argc > 4) ? std::stoul(argv[4]) : 1;
  const size_t jobs = (argc > 5) ? std::stoul(argv[5]) : 0;
  const size_t d = Geometry::Polytope::Traits(g).getDimension();

  for (size_t degree = from; degree <= maxDegree; ++degree)
  {
    const auto rule =
      SymmetricRuleGenerator::search(g, degree, 0, 64, 1e-12, jobs, seconds);
    if (!rule.converged || !rule.admissible)
    {
      std::printf("    // degree %zu: NOT FOUND (residual %.2e)%s\n", degree,
        rule.residual, (seconds > 0) ? " -- deadline reached; retry without one" : "");
      std::fflush(stdout);
      continue;
    }

    const auto report = exactnessSweep(rule, g, degree);
    if (report.worstRelativeError >= 1e-12 || !allWeightsPositive(rule) ||
      !allPointsInside(rule, g, 1e-12))
    {
      std::printf(
        "    // degree %zu: REJECTED (oracle %.2e)\n", degree, report.worstRelativeError);
      std::fflush(stdout);
      continue;
    }

    std::printf("    // degree %zu: %zu points, oracle %.2e\n", degree, rule.getSize(),
      report.worstRelativeError);
    std::printf("      {\n");
    for (size_t i = 0; i < rule.getSize(); ++i)
    {
      std::printf("        ");
      for (size_t k = 0; k < d; ++k)
        std::printf("%.17g, ", rule.getPoint(i)[static_cast<Eigen::Index>(k)]);
      std::printf("%.17g,\n", rule.getWeight(i));
    }
    std::printf("      },\n");
    std::fflush(stdout);
  }
  return 0;
}
