/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_ASSEMBLY_H
#define RODIN_VARIATIONAL_ASSEMBLY_H

#include <variant>
#include <ostream>

#include <mfem.hpp>

#include "Rodin/Geometry/Element.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
   enum class Backend
   {
      Native
   };

   std::ostream& operator<<(std::ostream& os, const Backend& backend);
}

#endif
