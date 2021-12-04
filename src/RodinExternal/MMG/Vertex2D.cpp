/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <common/mmgcommon.h>

#include "Vertex2D.h"

namespace Rodin::External::MMG
{
   Vertex2D::Vertex2D(MMG5_pPoint point)
      : m_point(point)
   {}

   double& Vertex2D::x()
   {
      return getHandle()->c[0];
   }

   double& Vertex2D::y()
   {
      return getHandle()->c[1];
   }

   const double& Vertex2D::x() const
   {
      return getHandle()->c[0];
   }

   const double& Vertex2D::y() const
   {
      return getHandle()->c[1];
   }

   int& Vertex2D::label()
   {
      return getHandle()->ref;
   }

   bool Vertex2D::isCorner() const
   {
      return getHandle()->tag & MG_CRN;
   }

   bool Vertex2D::isRequired() const
   {
      return getHandle()->tag & MG_REQ;
   }

   MMG5_pPoint& Vertex2D::getHandle()
   {
      return m_point;
   }

   const MMG5_pPoint& Vertex2D::getHandle() const
   {
      return m_point;
   }

   std::ostream& operator<<(std::ostream& os, const Vertex2D& obj)
   {
       os << "(" << obj.x() << ", " << obj.y() << ")";
       return os;
   }
}
