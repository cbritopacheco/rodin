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
   class SimplexIteratorBase
   {
      public:
         virtual ~SimplexIteratorBase() = default;

         virtual bool end() const = 0;

         virtual SimplexIteratorBase& operator++() = 0;

         virtual Simplex& operator*() const noexcept = 0;

         virtual Simplex* operator->() const noexcept = 0;

         virtual size_t getDimension() const = 0;

         virtual const MeshBase& getMesh() const = 0;

         virtual const IndexGeneratorBase& getIndexGenerator() const = 0;
   };

   class SimplexIterator final : public SimplexIteratorBase
   {
      public:
         SimplexIterator(size_t dimension, const MeshBase& mesh, IndexGeneratorBase&& gen);

         SimplexIterator(const SimplexIterator& other);

         SimplexIterator(SimplexIterator&& other);

         bool end() const override;

         SimplexIterator& operator++() override;

         Simplex& operator*() const noexcept override;

         Simplex* operator->() const noexcept override;

         size_t getDimension() const override
         {
            return m_dimension;
         }

         const MeshBase& getMesh() const override
         {
            return m_mesh.get();
         }

         const IndexGeneratorBase& getIndexGenerator() const override
         {
            assert(m_gen);
            return *m_gen;
         }

      private:
         Simplex* generate() const;

         IndexGeneratorBase& getIndexGenerator()
         {
            assert(m_gen);
            return *m_gen;
         }

         const size_t m_dimension;
         std::reference_wrapper<const MeshBase> m_mesh;
         std::unique_ptr<IndexGeneratorBase> m_gen;
         std::unique_ptr<Simplex> m_simplex;
   };

   class ElementIterator final : public SimplexIteratorBase
   {
      public:
         ElementIterator(const MeshBase& mesh, IndexGeneratorBase&& gen);

         ElementIterator(const ElementIterator& other);

         ElementIterator(ElementIterator&& other);

         bool end() const override;

         ElementIterator& operator++() override;

         Element& operator*() const noexcept override;

         Element* operator->() const noexcept override;

         size_t getDimension() const override;

         const MeshBase& getMesh() const override
         {
            return m_mesh.get();
         }

         const IndexGeneratorBase& getIndexGenerator() const override
         {
            assert(m_gen);
            return *m_gen;
         }

      private:
         Element* generate() const;

         IndexGeneratorBase& getIndexGenerator()
         {
            assert(m_gen);
            return *m_gen;
         }

         std::reference_wrapper<const MeshBase> m_mesh;
         std::unique_ptr<IndexGeneratorBase> m_gen;
         std::unique_ptr<Element> m_simplex;
   };

   class FaceIterator : public SimplexIteratorBase
   {
      public:
         FaceIterator(const MeshBase& mesh, IndexGeneratorBase&& gen);

         FaceIterator(const FaceIterator& other);

         FaceIterator(FaceIterator&& other);

         bool end() const override;

         FaceIterator& operator++() override;

         Face& operator*() const noexcept override;

         Face* operator->() const noexcept override;

         size_t getDimension() const override;

         const MeshBase& getMesh() const override
         {
            return m_mesh.get();
         }

         const IndexGeneratorBase& getIndexGenerator() const override
         {
            assert(m_gen);
            return *m_gen;
         }

      private:
         Face* generate() const;

         IndexGeneratorBase& getIndexGenerator()
         {
            assert(m_gen);
            return *m_gen;
         }

         std::reference_wrapper<const MeshBase> m_mesh;
         std::unique_ptr<IndexGeneratorBase> m_gen;
         std::unique_ptr<Face> m_simplex;
   };

   class BoundaryIterator final : public SimplexIteratorBase
   {
      public:
         BoundaryIterator(const MeshBase& mesh, IndexGeneratorBase&& gen);

         BoundaryIterator(const BoundaryIterator& other);

         BoundaryIterator(BoundaryIterator&& other);

         bool end() const override;

         BoundaryIterator& operator++() override;

         Boundary& operator*() const noexcept override;

         Boundary* operator->() const noexcept override;

         size_t getDimension() const override;

         const MeshBase& getMesh() const override
         {
            return m_mesh.get();
         }

         const IndexGeneratorBase& getIndexGenerator() const override
         {
            assert(m_gen);
            return *m_gen;
         }

      private:
         Boundary* generate() const;

         IndexGeneratorBase& getIndexGenerator()
         {
            assert(m_gen);
            return *m_gen;
         }

         std::reference_wrapper<const MeshBase> m_mesh;
         std::unique_ptr<IndexGeneratorBase> m_gen;
         std::unique_ptr<Boundary> m_simplex;
   };

   class InterfaceIterator final : public SimplexIteratorBase
   {
      public:
         InterfaceIterator(const MeshBase& mesh, IndexGeneratorBase&& gen);

         InterfaceIterator(const InterfaceIterator& other);

         InterfaceIterator(InterfaceIterator&& other);

         bool end() const override;

         InterfaceIterator& operator++() override;

         Interface& operator*() const noexcept override;

         Interface* operator->() const noexcept override;

         size_t getDimension() const override;

         const MeshBase& getMesh() const override
         {
            return m_mesh.get();
         }

         const IndexGeneratorBase& getIndexGenerator() const override
         {
            assert(m_gen);
            return *m_gen;
         }

      private:
         Interface* generate() const;

         IndexGeneratorBase& getIndexGenerator()
         {
            assert(m_gen);
            return *m_gen;
         }

         std::reference_wrapper<const MeshBase> m_mesh;
         std::unique_ptr<IndexGeneratorBase> m_gen;
         std::unique_ptr<Interface> m_simplex;
   };

   class VertexIterator : public SimplexIteratorBase
   {};
}

#endif
