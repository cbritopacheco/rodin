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
         using Mesh<Traits::Serial>::Mesh;

         SubMesh(const SubMesh& other);

         /**
          * @brief Sets the SubMesh to Mesh vertex map
          *
          * Each key-value pair in the map indicates the correspondence between
          * vertex indices in the SubMesh and vertex indices in the parent
          * Mesh.
          */
         SubMesh& setVertexMap(boost::bimap<int, int> s2pv);

         /**
          * @brief Sets the reference to the parent Mesh object
          */
         SubMesh& setParent(const MeshBase& parent);

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

         Mesh<Traits::Serial>& mesh()
         {
            return *this;
         }

      private:
         std::optional<std::reference_wrapper<const MeshBase>> m_parent;
         std::optional<boost::bimap<int, int>> m_s2pv;
   };
}

#endif

