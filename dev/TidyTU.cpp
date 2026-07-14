/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file TidyTU.cpp
 * @brief Translation unit aggregating the public headers for clang-tidy.
 *
 * The library is header-heavy, so per-file analysis of .cpp sources would
 * miss most of the code. This TU includes the public module headers so a
 * single clang-tidy invocation (with HeaderFilterRegex from .clang-tidy)
 * checks the whole header surface. Used by dev/style_lint and the Style CI
 * workflow; not part of the library build.
 *
 * Optional-dependency headers (PETSc, MMG, MPI) are excluded so the TU
 * analyzes without those packages installed.
 */
#include <Rodin/Types.h>
#include <Rodin/Alert.h>
#include <Rodin/Math.h>
#include <Rodin/Geometry.h>
#include <Rodin/QF.h>
#include <Rodin/FormLanguage.h>
#include <Rodin/Variational.h>
#include <Rodin/Assembly.h>
#include <Rodin/Solver.h>
#include <Rodin/Solid.h>
#include <Rodin/Heart.h>
#include <Rodin/IO.h>
#include <Rodin/Serialization/Export.h>
#include <Rodin/Adaptation/WNGIR.h>

int main()
{
  return 0;
}
