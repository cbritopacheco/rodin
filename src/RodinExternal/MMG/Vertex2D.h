/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_RODININTEGRATION_MMG_VERTEX2D_H
#define RODIN_RODININTEGRATION_MMG_VERTEX2D_H

#include <cassert>
#include <ostream>
#include <mmg/mmg2d/libmmg2d.h>

#include "Vertex.h"

namespace Rodin::External::MMG
{
   class Vertex2D : public Vertex<2>
   {
      public:
         Vertex2D(MMG5_pPoint point);

         /**
          * @returns x-coordinate of the vertex
          */
         double& x();

         /**
          * @returns x-coordinate of the vertex
          */
         const double& x() const;

         /**
          * @returns y-coordinate of the vertex
          */
         double& y();

         /**
          * @returns y-coordinate of the vertex
          */
         const double& y() const;

         int& label();

         bool isCorner() const;

         bool isRequired() const;

         MMG5_pPoint& getHandle() override;
         const MMG5_pPoint& getHandle() const override;

         friend std::ostream& operator<<(std::ostream& os, const Vertex2D& obj);

      private:
         MMG5_pPoint m_point;
   };
}

#endif
