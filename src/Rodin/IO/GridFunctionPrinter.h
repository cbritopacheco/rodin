/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GRIDFUNCTIONPRINTER_H
#define RODIN_VARIATIONAL_GRIDFUNCTIONPRINTER_H

#include <functional>
#include <boost/filesystem.hpp>

#include "Rodin/Alert/Exception.h"
#include "ForwardDecls.h"
#include "Rodin/Variational/ForwardDecls.h"

#include "Printer.h"

namespace Rodin::IO
{
  /**
   * @brief Primary template for GridFunctionPrinter.
   *
   * This primary template provides a default definition so that the compiler
   * does not fail with "implicit instantiation of undefined template" when
   * GridFunctionBase::save() is compiled for FES/Data combinations that do not
   * have a matching partial specialization for a given file format.
   *
   * @note An exception is raised at runtime if this primary template is
   *       actually invoked.
   *
   * @see GridFunctionLoader
   */
  template <FileFormat Fmt, class FES, class Data>
  class GridFunctionPrinter : public IO::Printer<Variational::GridFunction<FES, Data>>
  {
    public:
      using ObjectType = Variational::GridFunction<FES, Data>;

      GridFunctionPrinter(const ObjectType& gf)
        : m_gf(gf)
      {}

      void print(std::ostream&) override
      {
        Alert::Exception()
          << "No GridFunctionPrinter specialization for this format/FES/Data combination."
          << Alert::Raise;
      }

      const ObjectType& getObject() const override { return m_gf.get(); }

    private:
      std::reference_wrapper<const ObjectType> m_gf;
  };

  template <FileFormat Fmt, class FES, class Data>
  class GridFunctionPrinterBase : public IO::Printer<Variational::GridFunction<FES, Data>>
  {
    public:
      using FESType = FES;

      using DataType = Data;

      /**
       * @brief Type of mesh object being printed.
       */
      using ObjectType = Variational::GridFunction<FES, Data>;

      GridFunctionPrinterBase(const ObjectType& gf)
        : m_gf(gf)
      {}

      const ObjectType& getObject() const override
      {
        return m_gf.get();
      }

    private:
      std::reference_wrapper<const ObjectType> m_gf;
  };
}

#endif
