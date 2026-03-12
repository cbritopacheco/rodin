/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_IO_MEDIT_H
#define RODIN_PETSC_IO_MEDIT_H

#include <petscvec.h>

#include <iomanip>
#include <limits>
#include <type_traits>

#include "Rodin/IO/MEDIT.h"
#include "Rodin/Context/Local.h"
#include "Rodin/MPI/Context/ForwardDecls.h"

namespace Rodin::IO
{
  template <class FES>
  class GridFunctionPrinter<
    FileFormat::MEDIT,
    FES,
    ::Vec> final
    : public GridFunctionPrinterBase<FileFormat::MEDIT, FES, ::Vec>
  {
    public:
      using FESType   = FES;
      using DataType  = ::Vec;
      using Parent    = GridFunctionPrinterBase<FileFormat::MEDIT, FES, DataType>;

      using RangeType  = typename FormLanguage::Traits<FESType>::RangeType;
      using ScalarType = typename FormLanguage::Traits<RangeType>::ScalarType;

      static_assert(std::is_same_v<ScalarType, PetscScalar>,
        "PETSc MEDIT printer: FES scalar type must be PetscScalar.");

      using MeshType        = typename FormLanguage::Traits<FESType>::MeshType;
      using MeshContextType = typename FormLanguage::Traits<MeshType>::ContextType;

      using Parent::Parent;

      void printData(std::ostream& os) override
      {
        os << std::setprecision(std::numeric_limits<PetscReal>::max_digits10);

        const auto& gf  = this->getObject();
        const auto& fes = gf.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();

        if constexpr (std::is_same_v<MeshContextType, Context::MPI>)
        {
          // Ensure ghost values are available for vertex evaluations.
          gf.acquire();
        }

        // MEDIT expects one entry per vertex; for vector-valued solutions MEDIT
        // expects components on the same line separated by spaces.
        //
        // Your base prints:
        //   ... << 1 << " " << (vdim>1 ? Vector : Real) << '\n';
        // which means:
        //   - Real  : print 1 value per vertex line
        //   - Vector: print vdim values per vertex line (same line)
        //
        // So we must print with spaces, NOT one-per-line-per-component.
        const Geometry::Polytope::Traits ts(Geometry::Polytope::Type::Point);
        for (auto it = mesh.getVertex(); !it.end(); ++it)
        {
          const Geometry::Point p(
            *it,
            ts.getVertex(0),
            it->getCoordinates());

          const auto val = gf(p);

          if constexpr (std::is_same_v<RangeType, PetscScalar>)
          {
            os << val << '\n';
          }
          else
          {
            static_assert(
              std::is_same_v<RangeType, Math::Vector<PetscScalar>>,
              "MEDIT PETSc printer expects Range to be PetscScalar or Math::Vector<PetscScalar>."
            );

            const size_t vdim = fes.getVectorDimension();
            for (size_t d = 0; d < vdim; ++d)
            {
              os << val[d];
              if (d + 1 < vdim) os << ' ';
            }
            os << '\n';
          }
        }

        if constexpr (std::is_same_v<MeshContextType, Context::MPI>)
        {
          gf.flush();
        }

        os << '\n';
      }
  };
}

#endif
