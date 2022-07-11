/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GRIDFUNCTIONLOADER_H
#define RODIN_VARIATIONAL_GRIDFUNCTIONLOADER_H

#include <map>
#include <optional>
#include <mfem.hpp>

#include "Rodin/Variational/ForwardDecls.h"

#include "ForwardDecls.h"

#include "Loader.h"

namespace Rodin::IO
{
   template <class FEC, class Trait>
   class GridFunctionLoaderBase : public IO::Loader<Variational::GridFunction<FEC, Trait>>
   {
      public:
         GridFunctionLoaderBase(Variational::GridFunction<FEC, Trait>& gf)
            : m_gf(gf)
         {}

      protected:
         Variational::GridFunction<FEC, Trait>& getObject() override
         {
            return m_gf;
         }

      private:
         Variational::GridFunction<FEC, Trait>& m_gf;
   };

   // ---- MFEM Format -------------------------------------------------------
   template <class FEC>
   class GridFunctionLoader<GridFunctionFormat::MFEM, FEC, Traits::Serial>
      : public GridFunctionLoaderBase<FEC, Traits::Serial>
   {
      public:
         GridFunctionLoader(Variational::GridFunction<FEC, Traits::Serial>& gf)
            : GridFunctionLoaderBase<FEC, Traits::Serial>(gf)
         {}

         void load(std::istream& is) override;
   };

   // ---- MEDIT Format ------------------------------------------------------
   template <class FEC>
   class GridFunctionLoader<GridFunctionFormat::MEDIT, FEC, Traits::Serial>
      : public GridFunctionLoaderBase<FEC, Traits::Serial>
   {
      public:
         GridFunctionLoader(Variational::GridFunction<FEC, Traits::Serial>& gf)
            : GridFunctionLoaderBase<FEC, Traits::Serial>(gf)
         {}

         void load(std::istream& is) override;
   };
}

#include "GridFunctionLoader.hpp"

#endif
