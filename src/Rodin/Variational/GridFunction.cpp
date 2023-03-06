/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "GridFunction.h"

#include "Rodin/IO/GridFunctionLoader.h"

#include "Sum.h"
#include "Mult.h"
#include "Minus.h"
#include "Division.h"
#include "BooleanFunction.h"
#include "FiniteElementSpace.h"


namespace Rodin::Variational
{
  // template <>
  // GridFunction<L2<Context::Serial>>&
  // GridFunction<L2<Context::Serial>>::load(
  //     const boost::filesystem::path& filename, IO::FileFormat fmt)
  // {
  //   mfem::named_ifgzstream input(filename.c_str());
  //   if (!input)
  //   {
  //     Alert::Exception()
  //       << "Failed to open " << filename << " for reading."
  //       << Alert::Raise;
  //   }

  //   switch (fmt)
  //   {
  //     case IO::FileFormat::MFEM:
  //     {
  //       IO::GridFunctionLoader<IO::FileFormat::MFEM, L2<Context::Serial>> loader(*this);
  //       loader.load(input);
  //       break;
  //     }
  //     case IO::FileFormat::MEDIT:
  //     {
  //       IO::GridFunctionLoader<IO::FileFormat::MEDIT, L2<Context::Serial>> loader(*this);
  //       loader.load(input);
  //       break;
  //     }
  //     default:
  //     {
  //       Alert::Exception()
  //         << "Loading from \"" << fmt << "\" format unssuported."
  //         << Alert::Raise;
  //     }
  //   }
  //   return *this;
  // }

  // template <>
  // void GridFunction<L2<Context::Serial>>
  // ::save(const boost::filesystem::path& filename, IO::FileFormat fmt, int precision)
  // const
  // {
  //   std::ofstream output(filename.c_str());
  //   if (!output)
  //   {
  //     Alert::Exception()
  //       << "Failed to open " << filename << " for writing."
  //       << Alert::Raise;
  //   }

  //   output.precision(precision);
  //   switch (fmt)
  //   {
  //     case IO::FileFormat::MFEM:
  //     {
  //       IO::GridFunctionPrinter<IO::FileFormat::MFEM, L2<Context::Serial>> printer(*this);
  //       printer.print(output);
  //       break;
  //     }
  //     case IO::FileFormat::MEDIT:
  //     {
  //       IO::GridFunctionPrinter<IO::FileFormat::MEDIT, L2<Context::Serial>> printer(*this);
  //       printer.print(output);
  //       break;
  //     }
  //     default:
  //     {
  //       Alert::Exception()
  //         << "Saving to \"" << fmt << "\" format unssuported."
  //         << Alert::Raise;
  //     }
  //   }
  // }
}

