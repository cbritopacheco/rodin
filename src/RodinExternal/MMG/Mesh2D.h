/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_RODININTEGRATION_MMG_MESH2D_H
#define RODIN_RODININTEGRATION_MMG_MESH2D_H

#include <string>
#include <vector>
#include <variant>
#include <optional>

#include <mmg/mmg2d/libmmg2d.h>

#include "Mesh.h"
#include "Vertex2D.h"

namespace Rodin::External::MMG
{
   /**
    * @brief Represents a 2D mesh.
    */
   class Mesh2D : public Mesh<2, Mesh2D>
   {
      public:
         class VertexIterator
         {
            public:
               using iterator_category = std::forward_iterator_tag;
               using difference_type   = std::ptrdiff_t;
               using value_type        = Vertex2D;
               using pointer           = Vertex2D*;

               VertexIterator(Mesh2D& mesh);

               Vertex2D& operator*() const;
               Vertex2D* operator->() const;
               VertexIterator& operator++();
               VertexIterator operator++(int);
               VertexIterator operator+(size_t n);
               VertexIterator& operator+=(size_t n);
               friend bool operator==(const VertexIterator& a, const VertexIterator& b);
               friend bool operator!=(const VertexIterator& a, const VertexIterator& b);

            private:
               std::vector<Vertex2D> m_vertices;
               std::vector<Vertex2D>::iterator m_it;
               size_t m_offset;
         };

         class EdgeIterator;
         class TriangleIterator;

         enum Entity
         {
            Vertex = MMG5_Vertex,
            Edge = MMG5_Edg,
            Triangle = MMG5_Triangle
         };

         static Mesh2D load(const std::string& filename);
         void save(const std::string& filename);

         /**
          * Creates an empty mesh
          */
         Mesh2D();

         Mesh2D(Mesh2D&& other);

         ~Mesh2D();

         template <Entity e = Vertex>
         int count() const;

         MMG5_pMesh& getHandle();

         const MMG5_pMesh& getHandle() const;

      private:
         MMG5_pMesh  m_mesh;
   };

   template <>
   int Mesh2D::count<Mesh2D::Vertex>() const;

   template <>
   int Mesh2D::count<Mesh2D::Edge>() const;

   template <>
   int Mesh2D::count<Mesh2D::Triangle>() const;
}

#endif
