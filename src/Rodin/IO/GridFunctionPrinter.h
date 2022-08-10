/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GRIDFUNCTIONPRINTER_H
#define RODIN_VARIATIONAL_GRIDFUNCTIONPRINTER_H

#include <map>
#include <optional>
#include <mfem.hpp>
#include <boost/filesystem.hpp>

#include "Rodin/Variational/ForwardDecls.h"

#include "ForwardDecls.h"

#include "Printer.h"

namespace Rodin::IO
{
   template <class FES>
   class GridFunctionPrinterBase : public IO::Printer<Variational::GridFunction<FES>>
   {
      public:
         GridFunctionPrinterBase(const Variational::GridFunction<FES>& gf)
            : m_gf(gf)
         {}

      protected:
         const Variational::GridFunction<FES>& getObject() const override
         {
            return m_gf;
         }

      private:
         const Variational::GridFunction<FES>& m_gf;
   };

   template <class FES>
   class GridFunctionPrinter<FileFormat::MFEM, FES>
      : public GridFunctionPrinterBase<FES>
   {
      public:
         GridFunctionPrinter(const Variational::GridFunction<FES>& gf)
            : GridFunctionPrinterBase<FES>(gf)
         {}

         void print(std::ostream& os) override;
   };

   template <class FES>
   class GridFunctionPrinter<FileFormat::MEDIT, FES>
      : public GridFunctionPrinterBase<FES>
   {
      public:
         GridFunctionPrinter(const Variational::GridFunction<FES>& gf)
            : GridFunctionPrinterBase<FES>(gf)
         {}

         void print(std::ostream& os) override;
   };
}

#include "GridFunctionPrinter.hpp"

#endif
