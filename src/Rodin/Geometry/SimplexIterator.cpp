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
   makeMeshElement(const MeshBase& mesh, Index index, Attribute attr, mfem::Geometry::Type type);

   std::unique_ptr<mfem::ElementTransformation>
   makeMeshElementTransformation(const MeshBase& mesh, Index index, Attribute attr);

   std::unique_ptr<mfem::Element>
   makeMeshFace(const MeshBase& mesh, Index index, Attribute attr, mfem::Geometry::Type type);

   std::unique_ptr<mfem::ElementTransformation>
   makeMeshFaceTransformation(const MeshBase& mesh, Index index, Attribute attr);

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
   makeMeshElement(const MeshBase& mesh, Index index, Attribute attr, mfem::Geometry::Type type)
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
   makeMeshFace(const MeshBase& mesh, Index index, Attribute attr, mfem::Geometry::Type type)
   {
      // Build element object
      mfem::Element* element = makeGeometry(type);
      mfem::Array<int> vs;
      mesh.getHandle().GetFaceVertices(index, vs);
      element->SetVertices(std::move(vs));
      element->SetAttribute(attr);
      return std::unique_ptr<mfem::Element>(element);
   }

   std::unique_ptr<mfem::ElementTransformation> makeMeshElementTransformation(
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

   std::unique_ptr<mfem::ElementTransformation> makeMeshFaceTransformation(
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
   // ---- MeshSimplexIterator -----------------------------------------------
   MeshSimplexIterator::MeshSimplexIterator(Data data)
      : m_data(std::move(data))
   {
      update(generate());
   }

   MeshSimplexIterator::MeshSimplexIterator(const MeshSimplexIterator& other)
      : m_data(other.m_data)
   {
      update(generate());
   }

   MeshSimplexIterator::MeshSimplexIterator(MeshSimplexIterator&& other)
      : m_data(std::move(other.m_data)),
        m_simplex(std::move(other.m_simplex))
   {}

   bool MeshSimplexIterator::end() const
   {
      return m_data.it == m_data.end;
   }

   MeshSimplexIterator& MeshSimplexIterator::operator++()
   {
      m_data.it++;
      update(generate());
      return *this;
   }

   MeshSimplexIterator MeshSimplexIterator::operator++(int)
   {
      auto r = *this;
      ++(*this);
      return r;
   }

   bool MeshSimplexIterator::operator==(const MeshSimplexIterator& other) const
   {
      return getData().it == other.getData().it && other.m_data.mesh.get() == other.m_data.mesh.get();
   }

   bool MeshSimplexIterator::operator!=(const MeshSimplexIterator& other) const
   {
      return getData().it != other.getData().it && other.m_data.mesh.get() != other.m_data.mesh.get();
   }

   Simplex& MeshSimplexIterator::operator*() const noexcept
   {
      assert(m_simplex);
      return *m_simplex;
   }

   Simplex* MeshSimplexIterator::operator->() const noexcept
   {
      return m_simplex.get();
   }

   void MeshSimplexIterator::update(Simplex* simplex)
   {
      m_simplex.reset(generate());
   }

   Simplex* MeshSimplexIterator::generate() const
   {
      const auto& data = getData();
      if (getData().it->getDimension() == getData().mesh.get().getDimension())
      {
         auto element =
            Internal::makeMeshElement(
               data.mesh.get(), data.it->getIndex(), data.it->getAttribute(),
               static_cast<mfem::Geometry::Type>(data.it->getType()));
         auto trans =
            Internal::makeMeshElementTransformation(
                  data.mesh.get(), data.it->getIndex(), data.it->getAttribute());
         return new Element(
               {data.it->getIndex(), data.mesh.get(), std::move(element), std::move(trans)});
      }
      else if (getData().it->getDimension() == getData().mesh.get().getDimension() - 1)
      {
         auto element =
            Internal::makeMeshFace(
               data.mesh.get(), data.it->getIndex(), data.it->getAttribute(),
               static_cast<mfem::Geometry::Type>(data.it->getType()));
         auto trans =
            Internal::makeMeshFaceTransformation(
                  data.mesh.get(), data.it->getIndex(), data.it->getAttribute());
         return new Face(
               {data.it->getIndex(), data.mesh.get(), std::move(element), std::move(trans)});
      }
      else if (getData().it->getDimension() == 0)
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

   // ---- MeshElementIterator -----------------------------------------------
   MeshElementIterator::MeshElementIterator(const MeshElementIterator& other)
      : MeshSimplexIterator(other)
   {}

   MeshElementIterator::MeshElementIterator(MeshElementIterator&& other)
      : MeshSimplexIterator(std::move(other))
   {}

   Element& MeshElementIterator::operator*() const noexcept
   {
      assert(dynamic_cast<Element*>(&MeshSimplexIterator::operator*()));
      return static_cast<Element&>(MeshSimplexIterator::operator*());
   }

   Element* MeshElementIterator::operator->() const noexcept
   {
      assert(dynamic_cast<Element*>(MeshSimplexIterator::operator->()));
      return static_cast<Element*>(MeshSimplexIterator::operator->());
   }

   Element* MeshElementIterator::generate() const
   {
      const auto& data = getData();
      auto element =
         Internal::makeMeshElement(
            data.mesh.get(), data.it->getIndex(), data.it->getAttribute(),
            static_cast<mfem::Geometry::Type>(data.it->getType()));
      auto trans =
         Internal::makeMeshElementTransformation(
               data.mesh.get(), data.it->getIndex(), data.it->getAttribute());
      return new Element({data.it->getIndex(), data.mesh.get(), std::move(element), std::move(trans)});
   }

   // ---- FaceIterator ------------------------------------------------------
   MeshFaceIterator::MeshFaceIterator(Data data)
      : MeshSimplexIterator(std::move(data))
   {}

   MeshFaceIterator::MeshFaceIterator(const MeshFaceIterator& other)
      : MeshSimplexIterator(other)
   {}

   MeshFaceIterator::MeshFaceIterator(MeshFaceIterator&& other)
      : MeshSimplexIterator(std::move(other))
   {}

   Face& MeshFaceIterator::operator*() const noexcept
   {
      assert(dynamic_cast<Face*>(&MeshSimplexIterator::operator*()));
      return static_cast<Face&>(MeshSimplexIterator::operator*());
   }

   Face* MeshFaceIterator::operator->() const noexcept
   {
      assert(dynamic_cast<Face*>(MeshSimplexIterator::operator->()));
      return static_cast<Face*>(MeshSimplexIterator::operator->());
   }

   Face* MeshFaceIterator::generate() const
   {
      const auto& data = getData();
      auto element =
         Internal::makeMeshFace(
            data.mesh.get(), data.it->getIndex(), data.it->getAttribute(),
            static_cast<mfem::Geometry::Type>(data.it->getType()));
      auto trans =
         Internal::makeMeshFaceTransformation(
               data.mesh.get(), data.it->getIndex(), data.it->getAttribute());
      return new Face({data.it->getIndex(), data.mesh.get(), std::move(element), std::move(trans)});
   }

   // ---- BoundaryIterator --------------------------------------------------
   MeshBoundaryIterator::MeshBoundaryIterator(Data params)
      : MeshFaceIterator(std::move(params))
   {}

   MeshBoundaryIterator::MeshBoundaryIterator(const MeshBoundaryIterator& other)
      : MeshFaceIterator(other)
   {}

   MeshBoundaryIterator::MeshBoundaryIterator(MeshBoundaryIterator&& other)
      : MeshFaceIterator(std::move(other))
   {}

   Boundary& MeshBoundaryIterator::operator*() const noexcept
   {
      assert(dynamic_cast<Boundary*>(&MeshFaceIterator::operator*()));
      return static_cast<Boundary&>(MeshFaceIterator::operator*());
   }

   Boundary* MeshBoundaryIterator::operator->() const noexcept
   {
      assert(dynamic_cast<Boundary*>(MeshSimplexIterator::operator->()));
      return static_cast<Boundary*>(MeshSimplexIterator::operator->());
   }

   Boundary* MeshBoundaryIterator::generate() const
   {
      const auto& data = getData();
      auto element =
         Internal::makeMeshFace(
            data.mesh.get(), data.it->getIndex(), data.it->getAttribute(),
            static_cast<mfem::Geometry::Type>(data.it->getType()));
      auto trans =
         Internal::makeMeshFaceTransformation(
               data.mesh.get(), data.it->getIndex(), data.it->getAttribute());
      assert(data.mesh.get().isBoundary(data.it->getIndex()));
      return new Boundary({data.it->getIndex(), data.mesh.get(), std::move(element), std::move(trans)});
   }

   // ---- InterfaceIterator --------------------------------------------------
   MeshInterfaceIterator::MeshInterfaceIterator(Data params)
      : MeshFaceIterator(std::move(params))
   {}

   MeshInterfaceIterator::MeshInterfaceIterator(const MeshInterfaceIterator& other)
      : MeshFaceIterator(other)
   {}

   MeshInterfaceIterator::MeshInterfaceIterator(MeshInterfaceIterator&& other)
      : MeshFaceIterator(std::move(other))
   {}

   Interface& MeshInterfaceIterator::operator*() const noexcept
   {
      assert(dynamic_cast<Interface*>(&MeshFaceIterator::operator*()));
      return static_cast<Interface&>(MeshFaceIterator::operator*());
   }

   Interface* MeshInterfaceIterator::operator->() const noexcept
   {
      assert(dynamic_cast<Interface*>(MeshSimplexIterator::operator->()));
      return static_cast<Interface*>(MeshSimplexIterator::operator->());
   }

   Interface* MeshInterfaceIterator::generate() const
   {
      const auto& data = getData();
      auto element =
         Internal::makeMeshFace(
            data.mesh.get(), data.it->getIndex(), data.it->getAttribute(),
            static_cast<mfem::Geometry::Type>(data.it->getType()));
      auto trans =
         Internal::makeMeshFaceTransformation(
               data.mesh.get(), data.it->getIndex(), data.it->getAttribute());
      assert(data.mesh.get().isInterface(data.it->getIndex()));
      return new Interface(
            {data.it->getIndex(), data.mesh.get(), std::move(element), std::move(trans)});
   }

   FaceElementIterator::FaceElementIterator(const FaceElementIterator& other)
      : m_data(other.m_data)
   {
      update(generate());
   }

   void FaceElementIterator::update(Element* simplex)
   {
      m_element.reset(simplex);
   }

   Element* FaceElementIterator::generate() const
   {
      Attribute attr = m_data.mesh.get().getHandle().GetAttribute(*m_it);
      auto element =
         Internal::makeMeshFace(
            m_data.mesh,
            *m_it,
            attr,
            m_data.mesh.get().getHandle().GetElementGeometry(*m_it));
      auto trans =
         Internal::makeMeshFaceTransformation(m_data.mesh.get(), *m_it, attr);
      return new Element({*m_it, m_data.mesh, std::move(element), std::move(trans)});
   }

   bool FaceElementIterator::end() const
   {
      return m_it == m_data.faceIndices.end();
   }

   FaceElementIterator& FaceElementIterator::operator++()
   {
      m_it++;
      return *this;
   }

   FaceElementIterator FaceElementIterator::operator++(int)
   {
      auto r = *this;
      ++(*this);
      return r;
   }

   bool FaceElementIterator::operator==(const FaceElementIterator& other) const
   {
      return m_it == other.m_it;
   }

   bool FaceElementIterator::operator!=(const FaceElementIterator& other) const
   {
      return m_it != other.m_it;
   }

   Element& FaceElementIterator::operator*() const noexcept
   {
      assert(m_element);
      return *m_element;
   }

   Element* FaceElementIterator::operator->() const noexcept
   {
      return m_element.get();
   }
}
