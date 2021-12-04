/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FINITEELEMENTSPACE_H
#define RODIN_VARIATIONAL_FINITEELEMENTSPACE_H

#include <mfem.hpp>

#include "Rodin/Mesh.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
   /**
    * @brief Base type for the finite element spaces.
    *
    * @tparam Derived Subclass which derives from FiniteElementSpace.
    */
   template <class Derived>
   class FiniteElementSpace
   {
      public:
         Mesh& getMesh()
         {
            return static_cast<Derived*>(this)->getMesh();
         }

         mfem::FiniteElementCollection& getFEC()
         {
            return static_cast<Derived*>(this)->getFEC();
         }

         mfem::FiniteElementSpace& getFES()
         {
            return static_cast<Derived*>(this)->getFES();
         }

         bool operator=(const FiniteElementSpace& other)
         {
            return this == &other;
         }
   };
}

#endif
