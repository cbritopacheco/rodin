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
   template <class FEC, class Trait>
   class GridFunctionPrinterBase : public IO::Printer<Variational::GridFunction<FEC, Trait>>
   {
      public:
         GridFunctionPrinterBase(const Variational::GridFunction<FEC, Trait>& gf)
            : m_gf(gf)
         {}

      protected:
         const Variational::GridFunction<FEC, Trait>& getObject() const override
         {
            return m_gf;
         }

      private:
         const Variational::GridFunction<FEC, Trait>& m_gf;
   };

   template <class FEC>
   class GridFunctionPrinter<FileFormat::MFEM, FEC, Traits::Serial>
      : public GridFunctionPrinterBase<FEC, Traits::Serial>
   {
      public:
         GridFunctionPrinter(const Variational::GridFunction<FEC, Traits::Serial>& gf)
            : GridFunctionPrinterBase<FEC, Traits::Serial>(gf)
         {}

         void print(std::ostream& os) override;
   };

   template <class FEC>
   class GridFunctionPrinter<FileFormat::MEDIT, FEC, Traits::Serial>
      : public GridFunctionPrinterBase<FEC, Traits::Serial>
   {
      public:
         GridFunctionPrinter(const Variational::GridFunction<FEC, Traits::Serial>& gf)
            : GridFunctionPrinterBase<FEC, Traits::Serial>(gf)
         {}

         void print(std::ostream& os) override;
   };
}

#include "GridFunctionPrinter.hpp"

#endif
