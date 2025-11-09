/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_IO_GRIDFUNCTIONPRINTER_H
#define RODIN_PETSC_IO_GRIDFUNCTIONPRINTER_H

#include <petscvec.h>

#include "Rodin/Context/Local.h"

#include "Rodin/IO/ForwardDecls.h"
#include "Rodin/IO/GridFunctionPrinter.h"

#include "Rodin/FormLanguage/Traits.h"

#include "Rodin/MPI/Context/ForwardDecls.h"

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

      using FESMeshType = typename FormLanguage::Traits<FESType>::MeshType;

      using FESMeshContextType = typename FormLanguage::Traits<FESMeshType>::ContextType;

      using Parent::Parent;

      void printData(std::ostream& os) override
      {
        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          const auto& gf = this->getObject();
          const size_t sz = gf.getSize();
          for (size_t i = 0; i < sz; ++i)
            os << gf[i] << '\n';
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          const auto& gf = this->getObject();
          const auto& fes = gf.getFiniteElementSpace();
          const auto& shard = fes.getShard();
          const size_t sz = shard.getSize();
          gf.acquire();
          for (size_t i = 0; i < sz; ++i)
            os << gf[fes.getGlobalIndex(i)] << '\n';
          gf.flush();
        }
        else
        {
          assert(false);
        }
      }
  };
}
#endif

