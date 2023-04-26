/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IO_MFEM_H
#define RODIN_IO_MFEM_H

#include <boost/bimap.hpp>
#include <boost/spirit/home/x3.hpp>

#include "Rodin/Types.h"
#include "Rodin/Alert.h"
#include "Rodin/Context.h"
#include "Rodin/Geometry/Types.h"

#include "MeshLoader.h"
#include "MeshPrinter.h"
#include "ForwardDecls.h"

namespace Rodin::IO::MFEM
{
}

namespace Rodin::IO
{
  template <>
  class MeshLoader<IO::FileFormat::MFEM, Context::Serial>
    : public MeshLoaderBase<Context::Serial>
  {
    public:
      using Object = Rodin::Geometry::Mesh<Context::Serial>;

      MeshLoader(Object& mesh)
        : MeshLoaderBase<Context::Serial>(mesh)
      {}

      void load(std::istream& is) override;
  };

  template <>
  class MeshPrinter<FileFormat::MFEM, Context::Serial>
    : public MeshPrinterBase<Context::Serial>
  {
    public:
      MeshPrinter(const Rodin::Geometry::Mesh<Context::Serial>& mesh)
        : MeshPrinterBase(mesh)
      {}

      void print(std::ostream& os) override;
  };
}
#endif
