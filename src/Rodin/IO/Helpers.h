/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IO_HELPERS_H
#define RODIN_IO_HELPERS_H

#include <ostream>
#include <boost/bimap.hpp>

#include "ForwardDecls.h"

namespace Rodin::IO
{
   std::map<std::string, IO::MeshFormat> getMeshFileHeaders();
   static const std::map<std::string, IO::MeshFormat> MeshFileHeaders = getMeshFileHeaders();

   std::optional<IO::MeshFormat> getMeshFormat(std::istream& input);
}

namespace Rodin::IO::Medit
{
   enum class Keyword
   {
      MeshVersionFormatted,
      Dimension,
      Vertices,
      Triangles,
      Tetrahedra,
      Edges,
      SolAtVertices,
      SolAtEdges,
      SolAtTriangles,
      SolAtQuadrilaterals,
      SolAtTetrahedra,
      SolAtPentahedra,
      SolAtHexahedra,
      End
   };

   enum SolutionType
   {
      Scalar = 1,
      Vector = 2,
      Tensor = 3
   };

   boost::bimap<std::string, Keyword> getKeywordMap();
   static const boost::bimap<std::string, Keyword> KeywordMap = getKeywordMap();

   std::ostream& operator<<(std::ostream& os, Keyword kw);
}
#endif
