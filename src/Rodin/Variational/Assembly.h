/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_ASSEMBLY_H
#define RODIN_VARIATIONAL_ASSEMBLY_H

#include <mfem.hpp>

#include "ForwardDecls.h"

namespace Rodin::Variational::Linear::Assembly
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

namespace Rodin::Variational::Bilinear::Assembly
{
   enum class Type
   {
      Common
   };

   struct Common
   {
      const mfem::FiniteElement& trial;
      const mfem::FiniteElement& test;
      mfem::ElementTransformation& trans;
      mfem::DenseMatrix& mat;
   };
}

#endif
