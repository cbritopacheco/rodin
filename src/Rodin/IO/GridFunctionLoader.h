/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GRIDFUNCTIONLOADER_H
#define RODIN_VARIATIONAL_GRIDFUNCTIONLOADER_H

#include <functional>

#include "Rodin/Alert/Exception.h"
#include "Rodin/Variational/ForwardDecls.h"

#include "ForwardDecls.h"

#include "Loader.h"

namespace Rodin::IO
{
  /**
   * @brief Primary template for GridFunctionLoader.
   *
   * This primary template provides a default definition so that the compiler
   * does not fail with "implicit instantiation of undefined template" when
   * GridFunctionBase::load() is compiled for FES/Data combinations that do not
   * have a matching partial specialization for a given file format.
   *
   * @note If you reach this template at runtime it means no specialization
   *       exists for the requested format/FES/Data triple.  An exception is
   *       raised so that the unused switch branches in
   *       GridFunctionBase::load() remain compilable while still producing a
   *       clear error at runtime.
   *
   * @see Loader, GridFunctionPrinter
   */
  template <FileFormat Fmt, class FES, class Data>
  class GridFunctionLoader : public IO::Loader<Variational::GridFunction<FES, Data>>
  {
    public:
      using ObjectType = Variational::GridFunction<FES, Data>;

      GridFunctionLoader(ObjectType& gf)
        : m_gf(gf)
      {}

      void load(std::istream&) override
      {
        Alert::Exception()
          << "No GridFunctionLoader specialization for this format/FES/Data combination."
          << Alert::Raise;
      }

    protected:
      ObjectType& getObject() override { return m_gf.get(); }

    private:
      std::reference_wrapper<ObjectType> m_gf;
  };

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
