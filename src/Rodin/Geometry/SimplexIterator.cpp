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
      trans->Reset();

      const mfem::Mesh& meshHandle = mesh.getHandle();
      const mfem::GridFunction* nodes = meshHandle.GetNodes();

      if (!nodes)
      {
         meshHandle.GetPointMatrix(index, trans->GetPointMat());
         trans->SetFE(
               meshHandle.GetTransformationFEforElementType(
                  meshHandle.GetFaceElementType(index)));
         mfem::DenseMatrix& pm = trans->GetPointMat();
         mfem::Array<int> v;
         meshHandle.GetFaceVertices(index, v);
         const int nv = v.Size();
         const int spaceDim = mesh.getSpaceDimension();
         pm.SetSize(spaceDim, nv);

         for (int i = 0; i < spaceDim; i++)
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
   {
      update(generate());
   }

   SimplexIterator::SimplexIterator(const SimplexIterator& other)
      : m_dimension(other.m_dimension), m_mesh(other.m_mesh), m_gen(other.m_gen->copy())
   {
      update(generate());
   }

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
      update(generate());
      return *this;
   }

   Simplex& SimplexIterator::operator*() const noexcept
   {
      assert(m_simplex);
      return *m_simplex;
   }

   Simplex* SimplexIterator::operator->() const noexcept
   {
      return m_simplex.get();
   }

   void SimplexIterator::update(Simplex* simplex)
   {
      m_simplex.reset(simplex);
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
         return new Element({index, mesh, std::move(element), std::move(trans)});
      }
      else if (dimension == mesh.getDimension() - 1)
      {
         const auto& attribute = mesh.getAttribute(dimension, index);
         const auto& geometry = mesh.getHandle().GetFaceGeometry(index);
         auto element = Internal::makeFace(mesh, index, attribute, geometry);
         auto trans = Internal::makeFaceTransformation(mesh, index, attribute);
         return new Face({index, mesh, std::move(element), std::move(trans)});
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
      : SimplexIterator(mesh.getDimension(), mesh, std::move(gen))
   {}

   ElementIterator::ElementIterator(const ElementIterator& other)
      : SimplexIterator(other)
   {}

   ElementIterator::ElementIterator(ElementIterator&& other)
      : SimplexIterator(std::move(other))
   {}

   Element& ElementIterator::operator*() const noexcept
   {
      assert(dynamic_cast<Element*>(&SimplexIterator::operator*()));
      return static_cast<Element&>(SimplexIterator::operator*());
   }

   Element* ElementIterator::operator->() const noexcept
   {
      assert(dynamic_cast<Element*>(SimplexIterator::operator->()));
      return static_cast<Element*>(SimplexIterator::operator->());
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
      return new Element({index, mesh, std::move(element), std::move(trans)});
   }

   // ---- FaceIterator ------------------------------------------------------
   FaceIterator::FaceIterator(const MeshBase& mesh, IndexGeneratorBase&& gen)
      : SimplexIterator(mesh.getDimension() - 1, mesh, std::move(gen))
   {}

   FaceIterator::FaceIterator(const FaceIterator& other)
      : SimplexIterator(other)
   {}

   FaceIterator::FaceIterator(FaceIterator&& other)
      : SimplexIterator(std::move(other))
   {}

   Face& FaceIterator::operator*() const noexcept
   {
      assert(dynamic_cast<Face*>(&SimplexIterator::operator*()));
      return static_cast<Face&>(SimplexIterator::operator*());
   }

   Face* FaceIterator::operator->() const noexcept
   {
      assert(dynamic_cast<Face*>(SimplexIterator::operator->()));
      return static_cast<Face*>(SimplexIterator::operator->());
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
      return new Face({index, mesh, std::move(element), std::move(trans)});
   }

   // ---- BoundaryIterator --------------------------------------------------
   BoundaryIterator::BoundaryIterator(const MeshBase& mesh, IndexGeneratorBase&& gen)
      : FaceIterator(mesh, std::move(gen))
   {}

   BoundaryIterator::BoundaryIterator(const BoundaryIterator& other)
      : FaceIterator(other)
   {}

   BoundaryIterator::BoundaryIterator(BoundaryIterator&& other)
      : FaceIterator(std::move(other))
   {}

   Boundary& BoundaryIterator::operator*() const noexcept
   {
      assert(dynamic_cast<Boundary*>(&FaceIterator::operator*()));
      return static_cast<Boundary&>(FaceIterator::operator*());
   }

   Boundary* BoundaryIterator::operator->() const noexcept
   {
      assert(dynamic_cast<Boundary*>(SimplexIterator::operator->()));
      return static_cast<Boundary*>(SimplexIterator::operator->());
   }

   Boundary* BoundaryIterator::generate() const
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
      return new Boundary({index, mesh, std::move(element), std::move(trans)});
   }

   // ---- InterfaceIterator -------------------------------------------------
   InterfaceIterator::InterfaceIterator(const MeshBase& mesh, IndexGeneratorBase&& gen)
      : FaceIterator(mesh, std::move(gen))
   {}

   InterfaceIterator::InterfaceIterator(const InterfaceIterator& other)
      : FaceIterator(other)
   {}

   InterfaceIterator::InterfaceIterator(InterfaceIterator&& other)
      : FaceIterator(std::move(other))
   {}

   Interface& InterfaceIterator::operator*() const noexcept
   {
      assert(dynamic_cast<Interface*>(&FaceIterator::operator*()));
      return static_cast<Interface&>(FaceIterator::operator*());
   }

   Interface* InterfaceIterator::operator->() const noexcept
   {
      assert(dynamic_cast<Interface*>(SimplexIterator::operator->()));
      return static_cast<Interface*>(SimplexIterator::operator->());
   }

   Interface* InterfaceIterator::generate() const
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
      return new Interface({index, mesh, std::move(element), std::move(trans)});
   }
}
