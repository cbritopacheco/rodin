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
#include "SimplexRTree.h"
#include "Element.h"

namespace Rodin::Geometry
{
   class SimplexIterator
   {
      public:
         virtual Simplex& operator*() const noexcept = 0;
         virtual Simplex* operator->() const noexcept = 0;
         virtual SimplexIterator& operator++() = 0;
         virtual bool end() const = 0;
   };

   class ElementIterator : public SimplexIterator
   {
      virtual Element& operator*() const noexcept override = 0;
      virtual Element* operator->() const noexcept override = 0;
   };

   class MeshSimplexIterator : public SimplexIterator
   {
      public:
         struct Data
         {
            SimplexRTree::const_iterator it;
            const SimplexRTree::const_iterator end;
            std::reference_wrapper<const MeshBase> mesh;
         };

         MeshSimplexIterator(Data params);

         MeshSimplexIterator(const MeshSimplexIterator& other);

         MeshSimplexIterator(MeshSimplexIterator&& other);

         operator bool() const
         {
            return end();
         }

         virtual bool end() const override;

         virtual ~MeshSimplexIterator() = default;

         virtual MeshSimplexIterator operator++(int);

         virtual bool operator==(const MeshSimplexIterator& other) const;

         virtual bool operator!=(const MeshSimplexIterator& other) const;

         virtual MeshSimplexIterator& operator++() override;

         virtual Simplex& operator*() const noexcept override;

         virtual Simplex* operator->() const noexcept override;

      protected:
         Data& getData()
         {
            return m_data;
         }

         const Data& getData() const
         {
            return m_data;
         }

         virtual void update(Simplex* simplex);

         virtual Simplex* generate() const;

      private:
         Data m_data;
         std::unique_ptr<Simplex> m_simplex;
   };

   class MeshElementIterator : public MeshSimplexIterator
   {
      public:
         MeshElementIterator(Data params);
         MeshElementIterator(const MeshElementIterator& other);
         MeshElementIterator(MeshElementIterator&& other);

         virtual Element& operator*() const noexcept override;
         virtual Element* operator->() const noexcept override;

      protected:
         virtual Element* generate() const override;
   };

   class MeshFaceIterator : public MeshSimplexIterator
   {
      public:
         using ParentGeometry = MeshSimplexIterator;

         MeshFaceIterator(Data params);
         MeshFaceIterator(const MeshFaceIterator& other);
         MeshFaceIterator(MeshFaceIterator&& other);

         virtual Face& operator*() const noexcept override;
         virtual Face* operator->() const noexcept override;

      protected:
         virtual Face* generate() const override;
   };

   class MeshVertexIterator final : public MeshSimplexIterator
   {
      public:
         MeshVertexIterator(Data data);
         MeshVertexIterator(const MeshVertexIterator& other);
         MeshVertexIterator(MeshVertexIterator&& other);

         virtual Vertex& operator*() const noexcept override;
         virtual Vertex* operator->() const noexcept override;

      protected:
         virtual Vertex* generate() const override;
   };

   class MeshBoundaryIterator final : public MeshFaceIterator
   {
      public:
         MeshBoundaryIterator(Data params);
         MeshBoundaryIterator(const MeshBoundaryIterator& other);
         MeshBoundaryIterator(MeshBoundaryIterator&& other);

         Boundary& operator*() const noexcept override;
         Boundary* operator->() const noexcept override;

      protected:
         Boundary* generate() const override;
   };

   class MeshInterfaceIterator final : public MeshFaceIterator
   {
      public:
         MeshInterfaceIterator(Data params);
         MeshInterfaceIterator(const MeshInterfaceIterator& other);
         MeshInterfaceIterator(MeshInterfaceIterator&& other);

         Interface& operator*() const noexcept override;
         Interface* operator->() const noexcept override;

      protected:
         Interface* generate() const override;
   };

   class AdjacentSimplexIterator
   {};

   class AdjacentElementIterator : public AdjacentSimplexIterator
   {};

   class AdjacentFaceIterator : public AdjacentSimplexIterator
   {};

   class SimplexVertexIterator
   {
      public:
         struct Data
         {
            std::reference_wrapper<const MeshBase> mesh;
            std::reference_wrapper<const Simplex> simplex;
         };

         SimplexVertexIterator(Data data);
         SimplexVertexIterator(const SimplexVertexIterator& other);
         SimplexVertexIterator(SimplexVertexIterator&& other);

         Vertex& operator*() const noexcept;
         Vertex& operator->() const noexcept;
         SimplexVertexIterator& operator++();
         SimplexVertexIterator operator++(int);
         bool operator==(const SimplexVertexIterator& other) const;
         bool operator!=(const SimplexVertexIterator& other) const;

      protected:
         Vertex* generate() const;
   };

   class FaceElementIterator : public ElementIterator
   {
      public:
         struct Data
         {
            std::vector<Index> faceIndices;
            std::reference_wrapper<const MeshBase> mesh;
         };

         FaceElementIterator(Data data)
            : m_data(data), m_it(data.faceIndices.begin())
         {}

         FaceElementIterator(const FaceElementIterator& other);

         FaceElementIterator(FaceElementIterator&& other) = default;

         FaceElementIterator operator++(int);
         bool operator==(const FaceElementIterator& other) const;
         bool operator!=(const FaceElementIterator& other) const;

         bool end() const override;
         FaceElementIterator& operator++() override;
         Element& operator*() const noexcept override;
         Element* operator->() const noexcept override;

      protected:
         void update(Element* simplex);
         Element* generate() const;

      private:
         Data m_data;
         std::vector<Index>::iterator m_it;
         std::unique_ptr<Element> m_element;
   };
}

#endif
