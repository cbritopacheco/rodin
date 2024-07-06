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
#include <boost/filesystem.hpp>

#include "Rodin/Types.h"
#include "Rodin/Context.h"
#include "Rodin/Variational/ForwardDecls.h"

#include "ForwardDecls.h"

#include "Printer.h"

namespace Rodin::IO
{
  template <class FES>
  class GridFunctionPrinterBase : public IO::Printer<Variational::GridFunction<FES>>
  {
    public:
      using FESType = FES;

      using ObjectType = Variational::GridFunction<FESType>;

      using Parent = IO::Printer<ObjectType>;

      GridFunctionPrinterBase(const ObjectType& gf)
        : m_gf(gf)
      {}

      const Variational::GridFunction<FES>& getObject() const override
      {
        return m_gf.get();
      }

    private:
      std::reference_wrapper<const ObjectType> m_gf;
  };
}

#include "GridFunctionPrinter.hpp"

#endif
