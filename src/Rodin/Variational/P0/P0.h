/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P0_P0_H
#define RODIN_VARIATIONAL_P0_P0_H

#include "Rodin/Types.h"

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Connectivity.h"

#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/FiniteElementSpace.h"

#include "ForwardDecls.h"
#include "P0Element.h"

namespace Rodin::Variational
{
  /**
   * @defgroup P0Specializations P0 Template Specializations
   * @brief Template specializations of the P0 class.
   * @see P0
   */

  template <class Range, class Context, class Mesh = Geometry::Mesh<Context>>
  class P0;
}

#endif

