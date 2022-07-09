/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GRIDFUNCTION_HPP
#define RODIN_VARIATIONAL_GRIDFUNCTION_HPP

#include "GridFunction.h"

#include "Rodin/IO/GridFunctionLoader.h"

namespace Rodin::Variational
{
   template <class Trait>
   GridFunction<H1, Trait>& GridFunction<H1, Trait>
   ::load(const boost::filesystem::path& filename, IO::GridFunctionFormat fmt)
   {
      IO::Status status;
      mfem::named_ifgzstream input(filename.c_str());
      switch (fmt)
      {
         case IO::GridFunctionFormat::MFEM:
         {
            IO::GridFunctionLoader<IO::GridFunctionFormat::MFEM, H1, Trait> loader(*this);
            status = loader.load(input);
            *this = std::move(loader.getObject());
            break;
         }
         case IO::GridFunctionFormat::MEDIT:
         {
            IO::GridFunctionLoader<IO::GridFunctionFormat::MEDIT, H1, Trait> loader(*this);
            status = loader.load(input);
            *this = std::move(loader.getObject());
            break;
         }
      }

      if (!status.success)
      {
         Alert::Exception() << "Could not load GridFunction from file: " << filename << ". "
                            << (status.error ? status.error->message : "")
                            << Alert::Raise;
      }
      return *this;
   }
}

#endif
