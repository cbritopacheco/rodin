/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IO_GRIDFUNCTIONPRINTER_H
#define RODIN_IO_GRIDFUNCTIONPRINTER_H

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
   * @see <a href="class_rodin_1_1_i_o_1_1_grid_function_loader.html">GridFunctionLoader</a>
   */
  template <FileFormat Fmt, class FES, class Data>
  class GridFunctionPrinter : public IO::Printer<Variational::GridFunction<FES, Data>>
  {
    public:
      /// @brief Grid function type being printed.
      using ObjectType = Variational::GridFunction<FES, Data>;

      /**
       * @brief Constructs a fallback grid function printer.
       * @param[in] gf Grid function to print.
       */
      GridFunctionPrinter(const ObjectType& gf)
        : m_gf(gf)
      {}

      void print(std::ostream& os) override
      {
        Alert::Exception()
          << "No GridFunctionPrinter specialization for this format/FES/Data combination."
          << Alert::Raise;
      }

      const ObjectType& getObject() const override { return m_gf.get(); }

    private:
      std::reference_wrapper<const ObjectType> m_gf;
  };

  /**
   * @brief Base class for concrete grid function printers.
   *
   * Stores the grid function by reference and exposes it to derived printers
   * through @ref getObject().
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref GridFunctionPrinterBase "GridFunctionPrinterBase<FileFormat::MFEM, FES, Math::Vector<Scalar>>" | Base for MFEM grid function printers backed by local vector coefficients. |
   * | @ref GridFunctionPrinterBase "GridFunctionPrinterBase<FileFormat::MFEM, FES, Vec>" | Base for MFEM grid function printers backed by PETSc vectors. |
   * | @ref GridFunctionPrinterBase "GridFunctionPrinterBase<FileFormat::MEDIT, FES, Math::Vector<Scalar>>" | Base for MEDIT grid function printers backed by local vector coefficients. |
   * | @ref GridFunctionPrinterBase "GridFunctionPrinterBase<FileFormat::MEDIT, FES, Vec>" | Base for MEDIT grid function printers backed by PETSc vectors. |
   * | @ref GridFunctionPrinterBase "GridFunctionPrinterBase<FileFormat::HDF5, FES, Math::Vector<Scalar>>" | Base for HDF5 grid function printers backed by local vector coefficients. |
   * | @ref GridFunctionPrinterBase "GridFunctionPrinterBase<FileFormat::HDF5, FES, Vec>" | Base for HDF5 grid function printers backed by PETSc vectors. |
   */
  template <FileFormat Fmt, class FES, class Data>
  class GridFunctionPrinterBase : public IO::Printer<Variational::GridFunction<FES, Data>>
  {
    public:
      /// @brief File format handled by this printer base.
      static constexpr FileFormat Format = Fmt;

      /// @brief Finite element space type.
      using FESType = FES;

      /// @brief Coefficient data storage type.
      using DataType = Data;

      /**
       * @brief Type of mesh object being printed.
       */
      using ObjectType = Variational::GridFunction<FES, Data>;

      /**
       * @brief Constructs a grid function printer base.
       * @param[in] gf Grid function to print.
       */
      GridFunctionPrinterBase(const ObjectType& gf)
        : m_gf(gf)
      {}

      /**
       * @brief Returns the grid function bound to this printer.
       */
      const ObjectType& getObject() const override
      {
        return m_gf.get();
      }

    private:
      std::reference_wrapper<const ObjectType> m_gf;
  };
}

#endif
