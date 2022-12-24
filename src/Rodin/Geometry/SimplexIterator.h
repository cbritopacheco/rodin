/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_SIMPLEXITERATOR_H
#define RODIN_GEOMETRY_SIMPLEXITERATOR_H

#include <memory>
#include <utility>

#include "ForwardDecls.h"
#include "Element.h"

namespace Rodin::Geometry
{
   class SimplexIterator
   {
      using InternalIterator = std::vector<Index>::iterator;

      public:
         struct Data
         {
            const size_t dimension;
            std::reference_wrapper<const MeshBase> mesh;
            std::vector<Index> indices;
         };

         SimplexIterator(Data params);

         SimplexIterator(const SimplexIterator& other);

         SimplexIterator(SimplexIterator&& other);

         operator bool() const
         {
            return end();
         }

         virtual bool end() const;

         virtual ~SimplexIterator() = default;

         virtual SimplexIterator operator++(int);

         virtual bool operator==(const SimplexIterator& other) const;

         virtual bool operator!=(const SimplexIterator& other) const;

         virtual SimplexIterator& operator++();

         virtual Simplex& operator*() const noexcept;

         virtual Simplex* operator->() const noexcept;

      protected:
         Data& getData()
         {
            return m_data;
         }

         const Data& getData() const
         {
            return m_data;
         }

         InternalIterator& getInternalIterator()
         {
            return m_it;
         }

         const InternalIterator& getInternalIterator() const
         {
            return m_it;
         }

         virtual void update(Simplex* simplex);

         virtual Simplex* generate() const;

      private:
         Data m_data;
         InternalIterator m_it;
         std::unique_ptr<Simplex> m_simplex;
   };

   class ElementIterator : public SimplexIterator
   {
      public:
         ElementIterator(Data data);
         ElementIterator(const ElementIterator& other);
         ElementIterator(ElementIterator&& other);

         virtual Element& operator*() const noexcept override;
         virtual Element* operator->() const noexcept override;

      protected:
         virtual Element* generate() const override;
   };

   class FaceIterator : public SimplexIterator
   {
      public:
         using ParentGeometry = SimplexIterator;

         FaceIterator(Data params);
         FaceIterator(const FaceIterator& other);
         FaceIterator(FaceIterator&& other);

         virtual Face& operator*() const noexcept override;
         virtual Face* operator->() const noexcept override;

      protected:
         virtual Face* generate() const override;
   };

   class VertexIterator final : public SimplexIterator
   {
      public:
         VertexIterator(Data data);
         VertexIterator(const VertexIterator& other);
         VertexIterator(VertexIterator&& other);

         // virtual Vertex& operator*() const noexcept override;
         // virtual Vertex* operator->() const noexcept override;

      protected:
         // virtual Vertex* generate() const override;
   };

   class BoundaryIterator final : public FaceIterator
   {
      public:
         BoundaryIterator(Data params);
         BoundaryIterator(const BoundaryIterator& other);
         BoundaryIterator(BoundaryIterator&& other);

         Boundary& operator*() const noexcept override;
         Boundary* operator->() const noexcept override;

      protected:
         Boundary* generate() const override;
   };

   class InterfaceIterator final : public FaceIterator
   {
      public:
         InterfaceIterator(Data params);
         InterfaceIterator(const InterfaceIterator& other);
         InterfaceIterator(InterfaceIterator&& other);

         Interface& operator*() const noexcept override;
         Interface* operator->() const noexcept override;

      protected:
         Interface* generate() const override;
   };
}

#endif
