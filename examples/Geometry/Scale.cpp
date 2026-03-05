/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Alert/Exception.h"
#include "Rodin/Alert/Info.h"
#include "Rodin/Alert/Raise.h"
#include "Rodin/Alert/Warning.h"
#include <cstdlib>
#include <iostream>

#include <Rodin/Alert.h>
#include <Rodin/Geometry.h>

using namespace Rodin;
using namespace Rodin::Geometry;

static void usage(const char* prog)
{
  Alert::Warning()
    << "Invalid parameters.\nUsage: " << prog << " <input.mesh> <output.mesh> <scale>\n"
    << "Example: " << prog << " in.mesh out.mesh 0.5\n" << Alert::Raise;
}

int main(int argc, char** argv)
{
  if (argc != 4)
  {
    usage(argv[0]);
    return EXIT_FAILURE;
  }

  const char* inputFile  = argv[1];
  const char* outputFile = argv[2];

  char* end = nullptr;
  const Real s = std::strtod(argv[3], &end);
  if (end == argv[3] || *end != '\0')
  {
    Alert::Exception() << "Invalid scale: " << argv[3] << "\n" << Alert::Raise;
    usage(argv[0]);
    return EXIT_FAILURE;
  }

  Mesh mesh;
  mesh.load(inputFile, IO::FileFormat::MEDIT);
  mesh.scale(s);
  mesh.save(outputFile, IO::FileFormat::MEDIT);

  return EXIT_SUCCESS;
}
