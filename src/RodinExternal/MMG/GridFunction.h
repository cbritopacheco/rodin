/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_RODINEXTERNAL_MMG_GRIDFUNCTION_H
#define RODIN_RODINEXTERNAL_MMG_GRIDFUNCTION_H

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/P1.h"

namespace Rodin::External::MMG
{
  /// GridFunction class for use with the Rodin::External::MMG module.
  template <class Range>
  using GridFunction = Variational::GridFunction<
        Variational::P1<Range, Context::Sequential, Geometry::Mesh<Context::Sequential>>>;

  /// Type alias for a GridFunction with scalar range.
  using ScalarGridFunction = GridFunction<Scalar>;

  /// Type alias for a GridFunction with vector range.
  using VectorGridFunction = GridFunction<Math::Vector>;
}

#endif
