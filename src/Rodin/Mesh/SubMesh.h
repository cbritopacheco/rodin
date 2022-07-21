#ifndef RODIN_MESH_SUBMESH_H
#define RODIN_MESH_SUBMESH_H

#include <map>
#include <optional>
#include <functional>
#include <boost/bimap.hpp>

#include "ForwardDecls.h"
#include "Mesh.h"

namespace Rodin
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
   class SubMesh<Traits::Serial> : public Mesh<Traits::Serial>
   {
      public:
         SubMesh(const MeshBase& parent);

         SubMesh(const SubMesh& other);

         /**
          * @brief Sets the SubMesh to Mesh vertex map
          * @param[in] Map going from SubMesh to Mesh vertex indices.
          *
          * Each key-value pair in the map indicates the correspondence between
          * vertex indices in the SubMesh and vertex indices in the parent
          * Mesh.
          */
         SubMesh& setVertexMap(boost::bimap<int, int> s2pv);

         /**
          * @returns Reference to the parent Mesh object
          */
         const MeshBase& getParent() const;

         /**
          * @returns The SubMesh to Mesh vertex map
          */
         const boost::bimap<int, int>& getVertexMap() const;

         bool isSubMesh() const override
         {
            return true;
         }

         SubMesh& initialize(int dim, int sdim, int numVert = 0, int numElem = 0, int numBdrElem = 0)
         {
            getHandle() = mfem::Mesh(dim, numVert, numElem, numBdrElem, sdim);
            return *this;
         }

         void finalize()
         {
            getHandle().Finalize();
         }

         /**
          * @brief Adds an element into the mesh
          *
          * @param[in] el Element from the parent mesh
          */
         SubMesh& add(const ElementBase& el);

      private:
         const MeshBase& m_parent;
         boost::bimap<int, int> m_s2pv;
         boost::bimap<int, int> m_s2pf;
         boost::bimap<int, int> m_s2pe;
   };
}

#endif

