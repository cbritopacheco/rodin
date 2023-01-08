/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Mesh.h"
#include "SimplexIterator.h"

namespace Rodin::Geometry::Internal
{
   mfem::Element* makeGeometry(mfem::Geometry::Type t);

   std::unique_ptr<mfem::Element>
   makeElement(const MeshBase& mesh, Index index, Attribute attr, mfem::Geometry::Type type);

   std::unique_ptr<mfem::ElementTransformation>
   makeElementTransformation(const MeshBase& mesh, Index index, Attribute attr);

   std::unique_ptr<mfem::Element>
   makeFace(const MeshBase& mesh, Index index, Attribute attr, mfem::Geometry::Type type);

   std::unique_ptr<mfem::ElementTransformation>
   makeFaceTransformation(const MeshBase& mesh, Index index, Attribute attr);

   mfem::Element* makeGeometry(mfem::Geometry::Type t)
   {
      mfem::Element* element = nullptr;
      switch (t)
      {
         case mfem::Geometry::Type::POINT:
         {
            element = new mfem::Point;
            break;
         }
         case mfem::Geometry::Type::SEGMENT:
         {
            element = new mfem::Segment;
            break;
         }
         case mfem::Geometry::Type::TRIANGLE:
         {
            element = new mfem::Triangle;
            break;
         }
         case mfem::Geometry::Type::SQUARE:
         {
            element = new mfem::Quadrilateral;
            break;
         }
         case mfem::Geometry::Type::TETRAHEDRON:
         {
            element = new mfem::Tetrahedron;
            break;
         }
         case mfem::Geometry::Type::CUBE:
         {
            element = new mfem::Hexahedron;
            break;
         }
         case mfem::Geometry::Type::PRISM:
         {
            element = new mfem::Wedge;
            break;
         }
         case mfem::Geometry::Type::PYRAMID:
         {
            element = new mfem::Pyramid;
            break;
         }
         default:
         {
            assert(false);
            break;
         }
      }
      return element;
   }

   std::unique_ptr<mfem::Element>
   makeElement(const MeshBase& mesh, Index index, Attribute attr, mfem::Geometry::Type type)
   {
      // Build element object
      mfem::Element* element = makeGeometry(type);
      mfem::Array<int> vs;
      mesh.getHandle().GetElementVertices(index, vs);
      element->SetVertices(std::move(vs));
      element->SetAttribute(attr);
      return std::unique_ptr<mfem::Element>(element);
   }

   std::unique_ptr<mfem::Element>
   makeFace(const MeshBase& mesh, Index index, Attribute attr, mfem::Geometry::Type type)
   {
      // Build element object
      mfem::Element* element = makeGeometry(type);
      mfem::Array<int> vs;
      mesh.getHandle().GetFaceVertices(index, vs);
      element->SetVertices(std::move(vs));
      element->SetAttribute(attr);
      return std::unique_ptr<mfem::Element>(element);
   }

   std::unique_ptr<mfem::ElementTransformation> makeElementTransformation(
         const MeshBase& mesh, Index index, Attribute attr)
   {
      mfem::IsoparametricTransformation* trans = new mfem::IsoparametricTransformation;
      trans->Attribute = attr;
      trans->ElementNo = index;
      trans->ElementType = mfem::ElementTransformation::ELEMENT;
      trans->mesh = nullptr;
      trans->Reset();

      const mfem::Mesh& meshHandle = mesh.getHandle();
      const mfem::GridFunction* nodes = meshHandle.GetNodes();

      if (!nodes)
      {
         meshHandle.GetPointMatrix(index, trans->GetPointMat());
         trans->SetFE(
               meshHandle.GetTransformationFEforElementType(
                  meshHandle.GetElementType(index)));
      }
      else
      {
         assert(false);
      }
      return std::unique_ptr<mfem::ElementTransformation>(trans);
   }

   std::unique_ptr<mfem::ElementTransformation> makeFaceTransformation(
         const MeshBase& mesh, Index index, Attribute attr)
   {
      mfem::IsoparametricTransformation* trans = new mfem::IsoparametricTransformation;
      trans->Attribute = attr;
      trans->ElementNo = index;
      trans->ElementType = mfem::ElementTransformation::FACE;
      trans->mesh = nullptr;
      mfem::DenseMatrix& pm = trans->GetPointMat();
      trans->Reset();

      const mfem::Mesh& meshHandle = mesh.getHandle();
      const mfem::GridFunction* nodes = meshHandle.GetNodes();
      const size_t spaceDim = mesh.getSpaceDimension();

      if (!nodes)
      {

         mfem::Array<int> v;
         meshHandle.GetFaceVertices(index, v);
         const int nv = v.Size();
         pm.SetSize(spaceDim, nv);
         for (size_t i = 0; i < spaceDim; i++)
            for (int j = 0; j < nv; j++)
               pm(i, j) = meshHandle.GetVertex(v[j])[i];
         trans->SetFE(
               meshHandle.GetTransformationFEforElementType(
                  meshHandle.GetFaceElementType(index)));
      }
      else
      {
         assert(false);
      }

      return std::unique_ptr<mfem::ElementTransformation>(trans);
   }
}

namespace Rodin::Geometry
{
   // ---- SimplexIterator ---------------------------------------------------
   SimplexIterator::SimplexIterator(size_t dimension, const MeshBase& mesh, IndexGeneratorBase&& gen)
      : m_dimension(dimension), m_mesh(mesh), m_gen(std::move(gen).move())
   {}

   SimplexIterator::SimplexIterator(const SimplexIterator& other)
      : m_dimension(other.m_dimension), m_mesh(other.m_mesh), m_gen(other.m_gen->copy())
   {}

   SimplexIterator::SimplexIterator(SimplexIterator&& other)
      : m_dimension(other.m_dimension), m_mesh(other.m_mesh), m_gen(std::move(other.m_gen)),
        m_simplex(std::move(other.m_simplex))
   {}

   bool SimplexIterator::end() const
   {
      return getIndexGenerator().end();
   }

   SimplexIterator& SimplexIterator::operator++()
   {
      ++getIndexGenerator();
      return *this;
   }

   Simplex& SimplexIterator::operator*() const noexcept
   {
      m_simplex.reset(generate());
      return *m_simplex;
   }

   Simplex* SimplexIterator::operator->() const noexcept
   {
      m_simplex.reset(generate());
      return m_simplex.get();
   }

   Simplex* SimplexIterator::generate() const
   {
      if (end()) return nullptr;
      const auto& gen = getIndexGenerator();
      const auto& index = *gen;
      const auto& dimension = m_dimension;
      const auto& mesh = m_mesh.get();
      if (dimension == mesh.getDimension())
      {
         const auto& attribute = mesh.getAttribute(dimension, index);
         const auto& geometry = mesh.getHandle().GetElementGeometry(index);
         auto element = Internal::makeElement(mesh, index, attribute, geometry);
         auto trans = Internal::makeElementTransformation(mesh, index, attribute);
         return new Element(index, mesh, {std::move(element), std::move(trans)});
      }
      else if (dimension == mesh.getDimension() - 1)
      {
         const auto& attribute = mesh.getAttribute(dimension, index);
         const auto& geometry = mesh.getHandle().GetFaceGeometry(index);
         auto element = Internal::makeFace(mesh, index, attribute, geometry);
         auto trans = Internal::makeFaceTransformation(mesh, index, attribute);
         return new Face(index, mesh, {std::move(element), std::move(trans)});
      }
      else if (dimension == 0)
      {
         assert(false);
         return nullptr;
      }
      else
      {
         assert(false);
         return nullptr;
      }
   }

   // ---- ElementIterator ---------------------------------------------------
   ElementIterator::ElementIterator(const MeshBase& mesh, IndexGeneratorBase&& gen)
      : m_mesh(mesh), m_gen(std::move(gen).move())
   {}

   ElementIterator::ElementIterator(const ElementIterator& other)
      : m_mesh(other.m_mesh), m_gen(other.m_gen->copy())
   {}

   ElementIterator::ElementIterator(ElementIterator&& other)
      : m_mesh(other.m_mesh), m_gen(std::move(other.m_gen))
   {}

   bool ElementIterator::end() const
   {
      return getIndexGenerator().end();
   }

   ElementIterator& ElementIterator::operator++()
   {
      ++getIndexGenerator();
      return *this;
   }

   Element& ElementIterator::operator*() const noexcept
   {
      m_simplex.reset(generate());
      return *m_simplex;
   }

   Element* ElementIterator::operator->() const noexcept
   {
      m_simplex.reset(generate());
      return m_simplex.get();
   }

   size_t ElementIterator::getDimension() const
   {
      return getMesh().getDimension();
   }

   Element* ElementIterator::generate() const
   {
      if (end()) return nullptr;
      const auto& gen = getIndexGenerator();
      const auto& index = *gen;
      const auto& mesh = getMesh();
      const auto& dimension = getDimension();
      const auto& attribute = mesh.getAttribute(dimension, index);
      const auto& geometry = mesh.getHandle().GetElementGeometry(index);
      auto element = Internal::makeElement(mesh, index, attribute, geometry);
      auto trans = Internal::makeElementTransformation(mesh, index, attribute);
      return new Element(index, mesh, {std::move(element), std::move(trans)});
   }

   // ---- FaceIterator ------------------------------------------------------
   FaceIterator::FaceIterator(const MeshBase& mesh, IndexGeneratorBase&& gen)
      : m_mesh(mesh), m_gen(std::move(gen).move())
   {}

   FaceIterator::FaceIterator(const FaceIterator& other)
      : m_mesh(other.m_mesh), m_gen(other.m_gen->copy())
   {}

   FaceIterator::FaceIterator(FaceIterator&& other)
      : m_mesh(other.m_mesh), m_gen(std::move(other.m_gen))
   {}

   bool FaceIterator::end() const
   {
      return getIndexGenerator().end();
   }

   FaceIterator& FaceIterator::operator++()
   {
      ++getIndexGenerator();
      return *this;
   }

   Face& FaceIterator::operator*() const noexcept
   {
      m_simplex.reset(generate());
      return *m_simplex;
   }

   Face* FaceIterator::operator->() const noexcept
   {
      m_simplex.reset(generate());
      return m_simplex.get();
   }

   size_t FaceIterator::getDimension() const
   {
      return getMesh().getDimension() - 1;
   }

   Face* FaceIterator::generate() const
   {
      if (end()) return nullptr;
      const auto& gen = getIndexGenerator();
      const auto& index = *gen;
      const auto& mesh = getMesh();
      const auto& dimension = getDimension();
      const auto& attribute = mesh.getAttribute(dimension, index);
      const auto& geometry = mesh.getHandle().GetFaceGeometry(index);
      auto element = Internal::makeFace(mesh, index, attribute, geometry);
      auto trans = Internal::makeFaceTransformation(mesh, index, attribute);
      return new Face(index, mesh, {std::move(element), std::move(trans)});
   }
}
