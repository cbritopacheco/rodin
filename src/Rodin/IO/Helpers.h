/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IO_HELPERS_H
#define RODIN_IO_HELPERS_H

#include <boost/bimap.hpp>

#include "ForwardDecls.h"

namespace Rodin::IO
{
   std::map<std::string, IO::MeshFormat> getMeshFileHeaders();
   static const std::map<std::string, IO::MeshFormat> MeshFileHeaders = getMeshFileHeaders();

   std::optional<IO::MeshFormat> getMeshFormat(std::istream& input, bool seekBeg = true);
}

namespace Rodin::IO::Medit
{
   enum class SolKeyword
   {
      SolAtVertices,
      SolAtEdges,
      SolAtTriangles,
      SolAtQuadrilaterals,
      SolAtTetrahedra,
      SolAtPentahedra,
      SolAtHexahedra
   };
   boost::bimap<std::string, SolKeyword> getSolKeywordMap();
   static const boost::bimap<std::string, SolKeyword> SolKeywordMap = getSolKeywordMap();

   enum class EntityKeyword
   {
      Vertices,
      Triangles,
      Tetrahedra,
      Edges
   };
   boost::bimap<std::string, EntityKeyword> getEntityKeywordMap();
   static const boost::bimap<std::string, EntityKeyword> EntityKeywordMap = getEntityKeywordMap();
}
#endif
