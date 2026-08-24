/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
/**
 * @file
 * @brief Development driver emitting inlinable symmetric rule coefficients.
 *
 * Not a test. Each degree is emitted as it completes, so an interrupted run
 * keeps what it found, and every rule is checked against the moment oracle
 * before being printed --- a table is worth having only if what it holds is
 * exact.
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

int main(int argc, char** argv)
{
  const std::string which = (argc > 1) ? argv[1] : "tri";
  const auto g = (which == "tet") ? Geometry::Polytope::Type::Tetrahedron
    : (which == "wedge")          ? Geometry::Polytope::Type::Wedge
    : (which == "pyr")            ? Geometry::Polytope::Type::Pyramid
    : (which == "quad")           ? Geometry::Polytope::Type::Quadrilateral
    : (which == "hex")            ? Geometry::Polytope::Type::Hexahedron
                                  : Geometry::Polytope::Type::Triangle;
  const size_t maxDegree = (argc > 2) ? std::stoul(argv[2]) : 20;
  const size_t maxPoints = (argc > 3) ? std::stoul(argv[3]) : 130;
  const double seconds = (argc > 4) ? std::stod(argv[4]) : 900.0;
  const size_t d = Geometry::Polytope::Traits(g).getDimension();

  for (size_t degree = 1; degree <= maxDegree; ++degree)
  {
    const auto rule =
      SymmetricRuleGenerator::search(g, degree, maxPoints, 64, 1e-12, 0, seconds);
    if (!rule.converged || !rule.admissible)
    {
      std::printf(
        "    // degree %zu: NOT FOUND (residual %.2e)\n", degree, rule.residual);
      std::fflush(stdout);
      continue;
    }

    // A table is only worth having if what it holds is exact.
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
