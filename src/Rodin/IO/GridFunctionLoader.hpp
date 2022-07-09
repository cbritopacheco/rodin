/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GRIDFUNCTIONLOADER_HPP
#define RODIN_VARIATIONAL_GRIDFUNCTIONLOADER_HPP

#include "GridFunctionLoader.h"

#include "Rodin/Variational/GridFunction.h"

namespace Rodin::IO
{
   template <class FEC, class Trait>
   GridFunctionLoaderBase<FEC, Trait>::GridFunctionLoaderBase(
         Variational::GridFunction<FEC, Trait>& gf)
      : m_gf(gf)
   {}

   template <class FEC, class Trait>
   Variational::GridFunction<FEC, Trait>& GridFunctionLoaderBase<FEC, Trait>::getObject()
   {
      return m_gf;
   }

   template <class FEC>
   GridFunctionLoader<GridFunctionFormat::MFEM, FEC, Traits::Serial>
   ::GridFunctionLoader(Variational::GridFunction<FEC, Traits::Serial>& gf)
      : GridFunctionLoaderBase<FEC, Traits::Serial>(gf)
   {}

   template <class FEC>
   IO::Status GridFunctionLoader<GridFunctionFormat::MFEM, FEC, Traits::Serial>
   ::load(std::istream& is)
   {
      auto& gf = GridFunctionLoaderBase<FEC, Traits::Serial>::getObject();
      gf.getHandle() = mfem::GridFunction(
            &gf.getFiniteElementSpace().getMesh().getHandle(), is);
      return {true, {}};
   }

   template <class FEC>
   GridFunctionLoader<GridFunctionFormat::MEDIT, FEC, Traits::Serial>
   ::GridFunctionLoader(Variational::GridFunction<FEC, Traits::Serial>& gf)
      : GridFunctionLoaderBase<FEC, Traits::Serial>(gf)
   {}

   template <class FEC>
   IO::Status GridFunctionLoader<GridFunctionFormat::MEDIT, FEC, Traits::Serial>
   ::load(std::istream& is)
   {
      return {true, {}};
   }
}

#endif
