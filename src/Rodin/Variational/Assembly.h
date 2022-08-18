/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_ASSEMBLY_H
#define RODIN_VARIATIONAL_ASSEMBLY_H

#include <mfem.hpp>

namespace Rodin::Variational::Assembly
{
   enum class Type
   {
      Common,
      Device
   };

   struct Common
   {
      const mfem::FiniteElement& fe;
      mfem::ElementTransformation& trans;
      mfem::Vector& vec;
   };

   struct Device
   {
      const mfem::FiniteElementSpace& fes;
      const mfem::Array<int>& markers;
      mfem::Vector& vec;
   };
}

#endif
