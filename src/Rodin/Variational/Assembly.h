/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_ASSEMBLY_H
#define RODIN_VARIATIONAL_ASSEMBLY_H

#include <variant>

#include <mfem.hpp>

#include "Rodin/Geometry/Element.h"

#include "ForwardDecls.h"

namespace Rodin::Variational::Linear::Assembly
{
   enum class Type
   {
      Native, ///< Enumerator corresponding to Linear::Assembly::Native
      Device ///< Enumerator corresponding to Linear::Assembly::Device
   };

   struct Native
   {
      mfem::Vector& vec;
      const Geometry::ElementBase& element;
   };

   struct Device
   {
      const mfem::FiniteElementSpace& fes;
      const mfem::Array<int>& markers;
      mfem::Vector& vec;
   };
}

namespace Rodin::Variational::Bilinear::Assembly
{
   enum class Type
   {
      Native
   };

   struct Native
   {
      mfem::DenseMatrix& matrix;
      const Geometry::ElementBase& element;
   };
}

#endif
