/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_IO_GRIDFUNCTIONPRINTER_H
#define RODIN_PETSC_IO_GRIDFUNCTIONPRINTER_H

/**
 * @file
 * @brief Forward declaration of PETSc grid function printer specializations.
 *
 * The concrete specializations are provided in the format-specific headers
 * (HDF5.h, MFEM.h, MEDIT.h).
 *
 * @see Rodin::IO::GridFunctionPrinter
 */

#include <petscvec.h>

#include "Rodin/Context/Local.h"

#include "Rodin/IO/ForwardDecls.h"
#include "Rodin/IO/GridFunctionPrinter.h"

#include "Rodin/FormLanguage/Traits.h"

#include "Rodin/MPI/Context/ForwardDecls.h"

namespace Rodin::IO
{
  /**
   * @brief Partial specialization of @ref Rodin::IO::GridFunctionPrinter
   * for PETSc vector data (`::Vec`).
   *
   * @tparam Fmt  File format enumerator.
   * @tparam FES  Finite element space type.
   */
  template <FileFormat Fmt, class FES>
  class GridFunctionPrinter<Fmt, FES, ::Vec>;
}
#endif
