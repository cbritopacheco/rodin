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
      GridFunctionPrinterBase(const Variational::GridFunction<FES>& gf)
        : m_gf(gf)
      {}

      const Variational::GridFunction<FES>& getObject() const override
      {
        return m_gf;
      }

    private:
      const Variational::GridFunction<FES>& m_gf;
  };

  template <class Range, class ... Args>
  class GridFunctionPrinter<FileFormat::MFEM, Variational::P1<Range, Context::Serial, Args...>>
    : public GridFunctionPrinterBase<Variational::P1<Range, Context::Serial, Args...>>
  {
    public:
      using FES = Variational::P1<Range, Context::Serial, Args...>;

      GridFunctionPrinter(const Variational::GridFunction<FES>& gf)
        : GridFunctionPrinterBase<FES>(gf)
      {}

      void print(std::ostream& os) override
      {
        const auto& gf = this->getObject();
        const auto& fes = gf.getFiniteElementSpace();
        os << "FiniteElementSpace\n"
           << "FiniteElementCollection: " << "H1_" << fes.getMesh().getDimension() << "D_P1\n"
           << "VDim: " << fes.getVectorDimension() << '\n'
           << "Ordering: 1\n\n";
        const auto& matrix = gf.getData();
        const Scalar* data = matrix.data();
        assert(matrix.size() >= 0);
        for (size_t i = 0; i < static_cast<size_t>(matrix.size()); i++)
          os << data[i] << '\n';
      }
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
