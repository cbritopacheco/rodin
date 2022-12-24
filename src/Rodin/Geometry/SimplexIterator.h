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

#include "IndexGenerator.h"
#include "Element.h"
#include "ForwardDecls.h"

namespace Rodin::Geometry
{
   class SimplexIterator
   {
      public:
         SimplexIterator(size_t dimension, const MeshBase& mesh, IndexGeneratorBase&& gen);

         SimplexIterator(const SimplexIterator& other);

         SimplexIterator(SimplexIterator&& other);

         virtual ~SimplexIterator() = default;

         operator bool() const
         {
            return end();
         }

         virtual bool end() const;

         virtual SimplexIterator& operator++();

         virtual Simplex& operator*() const noexcept;

         virtual Simplex* operator->() const noexcept;

      protected:
         size_t getDimension() const
         {
            return m_dimension;
         }

         const MeshBase& getMesh() const
         {
            return m_mesh.get();
         }

         IndexGeneratorBase& getIndexGenerator()
         {
            assert(m_gen);
            return *m_gen;
         }

         const IndexGeneratorBase& getIndexGenerator() const
         {
            assert(m_gen);
            return *m_gen;
         }

         virtual void update(Simplex* simplex);

         virtual Simplex* generate() const;

      private:
         const size_t m_dimension;
         std::reference_wrapper<const MeshBase> m_mesh;
         std::unique_ptr<IndexGeneratorBase> m_gen;
         std::unique_ptr<Simplex> m_simplex;
   };

   class ElementIterator : public SimplexIterator
   {
      public:
         ElementIterator(const MeshBase& mesh, IndexGeneratorBase&& gen);
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

         FaceIterator(const MeshBase& mesh, IndexGeneratorBase&& gen);
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
         VertexIterator(const MeshBase& mesh, IndexGeneratorBase gen);
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
         BoundaryIterator(const MeshBase& mesh, IndexGeneratorBase&& gen);
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
         InterfaceIterator(const MeshBase& mesh, IndexGeneratorBase&& gen);
         InterfaceIterator(const InterfaceIterator& other);
         InterfaceIterator(InterfaceIterator&& other);

         Interface& operator*() const noexcept override;
         Interface* operator->() const noexcept override;

      protected:
         Interface* generate() const override;
   };
}

#endif
