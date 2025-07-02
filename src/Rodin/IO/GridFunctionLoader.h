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
  template <class FES, class Data>
  class GridFunctionLoaderBase : public IO::Loader<Variational::GridFunction<FES, Data>>
  {
    public:
      using FESType = FES;

      using DataType = Data;

      using ObjectType = Variational::GridFunction<FESType, Data>;

      using Parent = IO::Loader<ObjectType>;

      GridFunctionLoaderBase(ObjectType& gf)
        : m_gf(gf)
      {}

      ObjectType& getObject() override
      {
        return m_gf.get();
      }

    private:
      std::reference_wrapper<ObjectType> m_gf;
  };
}

#endif
