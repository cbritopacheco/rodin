/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GRIDFUNCTIONPRINTER_HPP
#define RODIN_VARIATIONAL_GRIDFUNCTIONPRINTER_HPP

#include "GridFunctionPrinter.h"

#include "MEDIT.h"

namespace Rodin::IO
{
   template <class FES>
   void GridFunctionPrinter<FileFormat::MEDIT, FES>::print(std::ostream& os)
   {
      assert(false);
      // int vertexCount = this->getObject().getFiniteElementSpace().getMesh().getHandle().GetNV();
      // int vdim = this->getObject().getFiniteElementSpace().getVectorDimension();
      // int spaceDim = this->getObject().getFiniteElementSpace().getMesh().getSpaceDimension();

      // os << MEDIT::Keyword::MeshVersionFormatted << " " << 2
      //    << '\n'
      //    << MEDIT::Keyword::Dimension << " " << spaceDim
      //    << "\n\n";

      // os << MEDIT::Keyword::SolAtVertices
      //    << '\n'
      //    << vertexCount
      //    << '\n'
      //    << 1 // Only one solution
      //    << " " << ((vdim > 1) ? MEDIT::SolutionType::Vector : MEDIT::SolutionType::Scalar)
      //    << '\n';

      // const auto& mfemGf = this->getObject().getHandle();
      // switch (this->getObject().getFiniteElementSpace().getHandle().GetOrdering())
      // {
      //    case mfem::Ordering::byNODES:
      //    {
      //       for (int i = 0; i < vertexCount; i++)
      //       {
      //          for (int j = 0; j < vdim - 1; j++)
      //             os << mfemGf[i + j * vertexCount] << " ";
      //          os << mfemGf[i + (vdim - 1) * vertexCount] << '\n';
      //       }
      //       break;
      //    }
      //    case mfem::Ordering::byVDIM:
      //    {
      //       for (int i = 0; i < vertexCount; i++)
      //       {
      //          for (int j = 0; j < vdim - 1; j++)
      //             os << mfemGf[j + i * vdim] << " ";
      //          os << mfemGf[(vdim - 1) + i * vdim] << '\n';
      //       }
      //       break;
      //    }
      // }
      // os << '\n' << MEDIT::Keyword::End;
   }
}

#endif
