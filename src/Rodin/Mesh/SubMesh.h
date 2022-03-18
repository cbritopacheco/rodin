#ifndef RODIN_MESH_SUBMESH_H
#define RODIN_MESH_SUBMESH_H

#include <map>
#include <optional>
#include <functional>

#include "ForwardDecls.h"
#include "Mesh.h"

namespace Rodin
{
   class SubMesh : public Mesh
   {
      public:
         using Mesh::Mesh;

         SubMesh(const SubMesh& other);

         /**
          * @brief Sets the SubMesh to Mesh vertex map
          *
          * Each key-value pair in the map indicates the correspondence between
          * vertex indices in the SubMesh and vertex indices in the parent
          * Mesh.
          */
         SubMesh& setVertexMap(std::map<int, int> s2pv);

         /**
          * @brief Sets the reference to the parent Mesh object
          */
         SubMesh& setParent(const Mesh& parent);

         /**
          * @returns Reference to the parent Mesh object
          */
         const Mesh& getParent();

         /**
          * @returns The SubMesh to Mesh vertex map
          */
         const std::map<int, int>& getVertexMap();

         bool isSubMesh() const override
         {
            return true;
         }

      private:
         std::optional<std::reference_wrapper<const Mesh>> m_parent;
         std::optional<std::map<int, int>> m_s2pv;
   };
}

#endif

