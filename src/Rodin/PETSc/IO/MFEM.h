/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_IO_MFEM_H
#define RODIN_PETSC_IO_MFEM_H

#include <petsc.h>
#include <petscsystypes.h>
#include <petscvec.h>

#include "Rodin/Context/Local.h"
#include "Rodin/IO/MFEM.h"

#include "Rodin/MPI/Context.h"

namespace Rodin::IO
{
  template <class FES>
  class GridFunctionPrinter<FileFormat::MFEM, FES, ::Vec>
    : public GridFunctionPrinterBase<FileFormat::MFEM, FES, ::Vec>
  {
    public:
      using FESType = FES;

      static constexpr FileFormat Format = FileFormat::MFEM;

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
        const PetscInt sz = 0;
        VecGetSize(data, &sz);
        const PetscScalar* raw = nullptr;
        VecGetArrayRead(data, &raw);
        for (PetscInt i = 0; i < sz; ++i)
          os << raw[i] << "\n";
        VecRestoreArrayRead(data, &raw);
      }
  };
}
#endif
