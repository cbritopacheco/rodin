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
#include "Rodin/IO/GridFunctionPrinter.h"

namespace Rodin::Variational
{
   template <class Trait>
   GridFunction<H1, Trait>& GridFunction<H1, Trait>
   ::load(const boost::filesystem::path& filename, IO::GridFunctionFormat fmt)
   {
      mfem::named_ifgzstream input(filename.c_str());
      if (!input)
      {
         Alert::Exception()
            << "Failed to open " << filename << " for reading."
            << Alert::Raise;
      }

      switch (fmt)
      {
         case IO::GridFunctionFormat::MFEM:
         {
            IO::GridFunctionLoader<IO::GridFunctionFormat::MFEM, H1, Trait> loader(*this);
            loader.load(input);
            break;
         }
         case IO::GridFunctionFormat::MEDIT:
         {
            IO::GridFunctionLoader<IO::GridFunctionFormat::MEDIT, H1, Trait> loader(*this);
            loader.load(input);
            break;
         }
      }
      return *this;
   }

   template <class Trait>
   void GridFunction<H1, Trait>
   ::save(const boost::filesystem::path& filename, IO::GridFunctionFormat fmt, int precision)
   const
   {
      std::ofstream output(filename.c_str());
      if (!output)
      {
         Alert::Exception()
            << "Failed to open " << filename << " for writing."
            << Alert::Raise;
      }

      output.precision(precision);
      switch (fmt)
      {
         case IO::GridFunctionFormat::MFEM:
         {
            IO::GridFunctionPrinter<IO::GridFunctionFormat::MFEM, H1, Trait> printer(*this);
            printer.print(output);
            break;
         }
         case IO::GridFunctionFormat::MEDIT:
         {
            IO::GridFunctionPrinter<IO::GridFunctionFormat::MEDIT, H1, Trait> printer(*this);
            printer.print(output);
            break;
         }
      }
   }
}

#endif
