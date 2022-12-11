#ifndef RODIN_MESH_SUBMESH_H
#define RODIN_MESH_SUBMESH_H

#include <map>
#include <optional>
#include <functional>
#include <boost/bimap.hpp>

#include "ForwardDecls.h"
#include "Mesh.h"

namespace Rodin::Geometry
{
   /**
    * @brief A SubMesh object represents a subregion of a Mesh object.
    *
    * A SubMesh object contains a reference to the parent Mesh object. It also
    * contains information regarding the mapping of elements and vertices
    * between the child and parent Mesh.
    *
    * A Mesh which is also a SubMesh may be casted into down to access
    * the SubMesh functionality. For example:
    * @code{.cpp}
    * if (mesh.isSubMesh())
    * {
    *    // Cast is well defined
    *    auto& submesh = static_cast<SubMesh&>(mesh);
    * }
    * @endcode
    *
    */
   template <>
   class SubMesh<Context::Serial> : public Mesh<Context::Serial>
   {
      public:
         SubMesh(const MeshBase& parent);

         SubMesh(const SubMesh& other);

         // SubMesh& initialize(int dim, int sdim, int numVert = 0, int numElem = 0, int numBdrElem = 0)
         // {
         //    getHandle() = mfem::Mesh(dim, numVert, numElem, numBdrElem, sdim);
         //    return *this;
         // }

         /**
          * @brief Adds an element from the parent mesh to the submesh.
          *
          * @param[in] el Element from the parent mesh
          */
         SubMesh& add(const Element& el);

         /**
          * Adds a bounday element from the parent mesh to the submesh.
          *
          * If `this->isSurface() && !getParent().isSurface()` then this will add the
          * boundary element as an element of the submesh. Otherwise, it gets
          * added normally as a boundary element.
          */
         SubMesh& add(const Boundary& el);

         /**
          * @returns Reference to the parent Mesh object
          */
         const MeshBase& getParent() const;

         /**
          * @returns The SubMesh to Mesh vertex map
          */
         const boost::bimap<size_t, size_t>& getVertexMap() const;

         const boost::bimap<size_t, size_t>& getElementMap() const;

         const boost::bimap<size_t, size_t>& getBoundaryElementMap() const;

         bool isSubMesh() const override
         {
            return true;
         }

      private:
         const MeshBase& m_parent;
         boost::bimap<size_t, size_t> m_s2pv;
         boost::bimap<size_t, size_t> m_s2pf;
         boost::bimap<size_t, size_t> m_s2pe;
         boost::bimap<size_t, size_t> m_s2pb;
   };
}

#endif

