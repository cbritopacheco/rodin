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
   SimplexIterator::SimplexIterator(Data data)
      : m_data(std::move(data)),
        m_it(m_data.indices.begin())
   {
      if (m_data.indices.size() > 0)
      {
         update(generate());
      }
      else
      {
         assert(m_data.indices.begin() == m_data.indices.end());
         assert(m_it == m_data.indices.end());
      }
   }

   SimplexIterator::SimplexIterator(const SimplexIterator& other)
      : m_data(other.m_data),
        m_it(other.m_it)
   {
      if (m_data.indices.size() > 0)
      {
         update(generate());
      }
      else
      {
         assert(m_data.indices.begin() == m_data.indices.end());
         assert(m_it == m_data.indices.end());
      }
   }

   SimplexIterator::SimplexIterator(SimplexIterator&& other)
      : m_data(std::move(other.m_data)),
        m_it(std::move(other.m_it)),
        m_simplex(std::move(other.m_simplex))
   {}

   bool SimplexIterator::end() const
   {
      return m_it == m_data.indices.end();
   }

   SimplexIterator& SimplexIterator::operator++()
   {
      m_it++;
      update(generate());
      return *this;
   }

   SimplexIterator SimplexIterator::operator++(int)
   {
      auto r = *this;
      ++(*this);
      return r;
   }

   bool SimplexIterator::operator==(const SimplexIterator& other) const
   {
      assert(m_data.mesh.get() == other.m_data.mesh.get());
      return m_it == other.m_it;
   }

   bool SimplexIterator::operator!=(const SimplexIterator& other) const
   {
      assert(m_data.mesh.get() == other.m_data.mesh.get());
      return m_it != other.m_it;
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
      const auto& index = *m_it;
      const auto& data = getData();
      const auto& dimension = data.dimension;
      const auto& mesh = data.mesh.get();
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
   ElementIterator::ElementIterator(Data data)
      : SimplexIterator(std::move(data))
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
      const auto& it = getInternalIterator();
      assert(it != getData().indices.end());
      const auto& index = *it;
      const auto& data = getData();
      const auto& mesh = data.mesh.get();
      const auto& dimension = data.dimension;
      const auto& attribute = mesh.getAttribute(dimension, index);
      const auto& geometry = mesh.getHandle().GetElementGeometry(index);
      auto element = Internal::makeElement(mesh, index, attribute, geometry);
      auto trans = Internal::makeElementTransformation(mesh, index, attribute);
      return new Element({index, mesh, std::move(element), std::move(trans)});
   }

   // ---- FaceIterator ------------------------------------------------------
   FaceIterator::FaceIterator(Data data)
      : SimplexIterator(std::move(data))
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
      const auto& it = getInternalIterator();
      assert(it != getData().indices.end());
      const auto& index = *it;
      const auto& data = getData();
      const auto& mesh = data.mesh.get();
      const auto& dimension = data.dimension;
      const auto& attribute = mesh.getAttribute(dimension, index);
      const auto& geometry = mesh.getHandle().GetFaceGeometry(index);
      auto element = Internal::makeFace(mesh, index, attribute, geometry);
      auto trans = Internal::makeFaceTransformation(mesh, index, attribute);
      return new Face({index, mesh, std::move(element), std::move(trans)});
   }

   // ---- BoundaryIterator --------------------------------------------------
   BoundaryIterator::BoundaryIterator(Data params)
      : FaceIterator(std::move(params))
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
      const auto& it = getInternalIterator();
      assert(it != getData().indices.end());
      const auto& index = *it;
      const auto& data = getData();
      const auto& mesh = data.mesh.get();
      const auto& dimension = data.dimension;
      const auto& attribute = mesh.getAttribute(dimension, index);
      const auto& geometry = mesh.getHandle().GetFaceGeometry(index);
      auto element = Internal::makeFace(mesh, index, attribute, geometry);
      auto trans = Internal::makeFaceTransformation(mesh, index, attribute);
      return new Boundary({index, mesh, std::move(element), std::move(trans)});
   }

   // ---- InterfaceIterator -------------------------------------------------
   InterfaceIterator::InterfaceIterator(Data params)
      : FaceIterator(std::move(params))
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
      const auto& it = getInternalIterator();
      assert(it != getData().indices.end());
      const auto& index = *it;
      const auto& data = getData();
      const auto& mesh = data.mesh.get();
      const auto& dimension = data.dimension;
      const auto& attribute = mesh.getAttribute(dimension, index);
      const auto& geometry = mesh.getHandle().GetFaceGeometry(index);
      auto element = Internal::makeFace(mesh, index, attribute, geometry);
      auto trans = Internal::makeFaceTransformation(mesh, index, attribute);
      return new Interface({index, mesh, std::move(element), std::move(trans)});
   }
}
