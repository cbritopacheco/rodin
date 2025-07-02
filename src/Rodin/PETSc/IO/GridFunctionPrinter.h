/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_IO_GRIDFUNCTIONPRINTER_H
#define RODIN_PETSC_IO_GRIDFUNCTIONPRINTER_H

#include <petsc.h>
#include <petscsystypes.h>
#include <petscvec.h>

#include "Rodin/IO/GridFunctionPrinter.h"
#include "Rodin/FormLanguage/Traits.h"
#include "Rodin/IO/ForwardDecls.h"

namespace Rodin::IO
{
  template <FileFormat Fmt, class FES>
  class GridFunctionPrinter<Fmt, FES, ::Vec>
    : public GridFunctionPrinterBase<Fmt, FES, ::Vec>
  {
    public:
      using FESType = FES;

      static constexpr FileFormat Format = Fmt;

      using RangeType = typename FormLanguage::Traits<FES>::RangeType;

      using ScalarType = typename FormLanguage::Traits<RangeType>::ScalarType;

      using DataType = ::Vec;

      using Parent = GridFunctionPrinterBase<Format, FES, DataType>;

      using ContextType = typename FormLanguage::Traits<FESType>::ContextType;

      using Parent::Parent;

      void printData(std::ostream& os) override
      {
        const auto& gf = this->getObject();
        const auto& data = gf.getData();
        const auto& raw = gf.getRaw();
        for (PetscInt i = 0; i < gf.getSize(); ++i)
          os << raw[i] << "\n";
      }
  };
}
#endif

