/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GRIDFUNCTIONLOADER_H
#define RODIN_VARIATIONAL_GRIDFUNCTIONLOADER_H

#include "Rodin/Variational/ForwardDecls.h"

#include "ForwardDecls.h"

#include "Loader.h"

namespace Rodin::IO
{
  /**
   * @brief Base class for loading grid functions from files or streams.
   *
   * GridFunctionLoaderBase provides the foundation for loading finite element
   * solution data in different file formats. It extends the generic Loader class
   * to handle grid function-specific operations.
   *
   * @tparam FES Finite element space type
   * @tparam Data Data storage type (typically a vector type)
   *
   * Specialized loaders for specific file formats should derive from this class
   * and implement the load(std::istream&) method to parse their respective formats.
   *
   * ## Usage Example
   * ```cpp
   * P1<Real> Vh(mesh);
   * GridFunction<P1<Real>> u(Vh);
   * GridFunctionLoader<FileFormat::MFEM, P1<Real>, Vector<Real>> loader(u);
   * loader.load("solution.gf");
   * ```
   *
   * @see Loader, GridFunctionPrinter
   */
  template <class FES, class Data>
  class GridFunctionLoaderBase : public IO::Loader<Variational::GridFunction<FES, Data>>
  {
    public:
      /**
       * @brief Finite element space type.
       */
      using FESType = FES;

      /**
       * @brief Data storage type.
       */
      using DataType = Data;

      /**
       * @brief Type of grid function object being loaded.
       */
      using ObjectType = Variational::GridFunction<FESType, Data>;

      /**
       * @brief Parent Loader class type.
       */
      using Parent = IO::Loader<ObjectType>;

      /**
       * @brief Constructs a grid function loader for the given grid function.
       * @param[in,out] gf Grid function object to be populated with loaded data
       *
       * The grid function is stored by reference and will be modified during loading.
       */
      GridFunctionLoaderBase(ObjectType& gf)
        : m_gf(gf)
      {}

      /**
       * @brief Gets a reference to the grid function being loaded.
       * @returns Reference to the grid function object
       *
       * Provides access to the grid function for derived classes to populate during loading.
       */
      ObjectType& getObject() override
      {
        return m_gf.get();
      }

    private:
      std::reference_wrapper<ObjectType> m_gf;
  };
}

#endif
