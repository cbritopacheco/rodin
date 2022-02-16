/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_MATRIXCOEFFICIENT_H
#define RODIN_VARIATIONAL_MATRIXCOEFFICIENT_H

#include <mfem.hpp>

#include "ForwardDecls.h"
#include "FormLanguage/Base.h"

#include "VectorCoefficient.h"

namespace Rodin::Variational
{
   class MatrixCoefficientBase
      : public FormLanguage::Buildable<mfem::MatrixCoefficient>
   {
      public:
         virtual ~MatrixCoefficientBase() = default;
         virtual Transpose T() const;

         virtual int getRows() const = 0;
         virtual int getColumns() const = 0;
         virtual void build() override = 0;
         virtual mfem::MatrixCoefficient& get() override = 0;
         virtual MatrixCoefficientBase* copy() const noexcept override = 0;
   };
}

#endif
