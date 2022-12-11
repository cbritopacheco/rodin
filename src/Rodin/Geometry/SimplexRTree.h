/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_SIMPLEXRTREE_H
#define RODIN_GEOMETRY_SIMPLEXRTREE_H

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>

#include "ForwardDecls.h"

namespace Rodin::Geometry
{
   class SimplexRTreeValue
   {
      public:
         SimplexRTreeValue(Dimension dimension, Region region, Attribute attr, Type type, Index index)
            : m_dimension(dimension), m_region(region), m_attr(attr), m_type(type), m_index(index)
         {}

         constexpr
         Dimension getDimension() const
         {
            return m_dimension;
         }

         constexpr
         Region getRegion() const
         {
            return m_region;
         }

         constexpr
         Attribute getAttribute() const
         {
            return m_attr;
         }

         constexpr
         Type getType() const
         {
            return m_type;
         }

         constexpr
         Index getIndex() const
         {
            return m_index;
         }

      private:
         Dimension m_dimension;
         Region m_region;
         Attribute m_attr;
         Type m_type;
         Index m_index;
   };

   class SimplexRTreeIndexableGetter
   {
      public:
         using Point =
            boost::geometry::model::point<size_t, 5, boost::geometry::cs::cartesian>;
         using result_type = Point;

         Point operator()(const SimplexRTreeValue& v) const
         {
            Point p;
            p.set<0>(v.getDimension());
            p.set<1>(v.getAttribute());
            p.set<2>(static_cast<size_t>(v.getRegion()));
            p.set<3>(v.getIndex());
            p.set<4>(static_cast<size_t>(v.getType()));
            return p;
         }
   };

   using SimplexRTree =
      boost::geometry::index::rtree<
         SimplexRTreeValue, boost::geometry::index::quadratic<16>, SimplexRTreeIndexableGetter>;
}

#endif
