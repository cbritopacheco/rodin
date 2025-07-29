/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_IO_GRIDFUNCTIONPRINTER_H
#define RODIN_PETSC_IO_GRIDFUNCTIONPRINTER_H

#include "Rodin/IO/ForwardDecls.h"
#include "Rodin/IO/GridFunctionPrinter.h"

#include "Rodin/FormLanguage/Traits.h"

#include "Rodin/PETSc/Math/Vector.h"
#include <petscsystypes.h>

namespace Rodin::IO
{
  template <FileFormat Fmt, class FES>
  class GridFunctionPrinter<Fmt, FES, PETSc::Math::Vector>
    : public GridFunctionPrinterBase<Fmt, FES, PETSc::Math::Vector>
  {
    public:
      using FESType = FES;

      static constexpr FileFormat Format = Fmt;

      using RangeType = typename FormLanguage::Traits<FES>::RangeType;

      using ScalarType = typename FormLanguage::Traits<RangeType>::ScalarType;

      using DataType = PETSc::Math::Vector;

      using Parent = GridFunctionPrinterBase<Format, FES, DataType>;

      using ContextType = typename FormLanguage::Traits<FESType>::ContextType;

      using Parent::Parent;

      void printData(std::ostream& os) override
      {
        const auto& gf = this->getObject();
        const size_t sz = gf.getSize();
        for (size_t i = 0; i < sz; ++i)
          os << gf[i] << "\n";
      }
  };
}
#endif

