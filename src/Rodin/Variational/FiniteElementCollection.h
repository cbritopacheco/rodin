/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FINITEELEMENTCOLLECTION_H
#define RODIN_VARIATIONAL_FINITEELEMENTCOLLECTION_H

#include <mfem.hpp>

namespace Rodin::Variational
{
   class FiniteElementCollectionBase
   {
      public:
         virtual mfem::FiniteElementCollection& getHandle() = 0;
         virtual const mfem::FiniteElementCollection& getHandle() const = 0;
   };
}

#endif
