/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GRIDFUNCTIONLOADER_H
#define RODIN_VARIATIONAL_GRIDFUNCTIONLOADER_H

#include <map>
#include <optional>

#include "Rodin/Variational/ForwardDecls.h"

#include "ForwardDecls.h"

#include "Loader.h"

namespace Rodin::IO
{
  template <class FES>
  class GridFunctionLoaderBase : public IO::Loader<Variational::GridFunction<FES>>
  {
    public:
      GridFunctionLoaderBase(Variational::GridFunction<FES>& gf)
        : m_gf(gf)
      {}

    protected:
      Variational::GridFunction<FES>& getObject() override
      {
        return m_gf;
      }

    private:
      Variational::GridFunction<FES>& m_gf;
  };

  // ---- MFEM Format -------------------------------------------------------
  template <class FES>
  class GridFunctionLoader<FileFormat::MFEM, FES>
    : public GridFunctionLoaderBase<FES>
  {
    public:
      GridFunctionLoader(Variational::GridFunction<FES>& gf)
        : GridFunctionLoaderBase<FES>(gf)
      {}

      void load(std::istream& is) override;
  };

  // ---- MEDIT Format ------------------------------------------------------
  template <class FES>
  class GridFunctionLoader<FileFormat::MEDIT, FES>
    : public GridFunctionLoaderBase<FES>
  {
    public:
      GridFunctionLoader(Variational::GridFunction<FES>& gf)
        : GridFunctionLoaderBase<FES>(gf)
      {}

      void load(std::istream& is) override;
  };
}

#include "GridFunctionLoader.hpp"

#endif
