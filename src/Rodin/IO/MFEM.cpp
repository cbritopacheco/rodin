/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "MFEM.h"

namespace Rodin::IO::MFEM
{}

namespace Rodin::IO
{
  void MeshLoader<IO::FileFormat::MFEM, Context::Serial>::load(std::istream &is)
  {
  }

  void MeshPrinter<FileFormat::MFEM, Context::Serial>::print(std::ostream &os)
  {
  }
}


