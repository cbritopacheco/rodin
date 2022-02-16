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

         const std::map<int, int>& getVertexMap();

         SubMesh& setVertexMap(std::map<int, int> s2pv);

         SubMesh& setParent(const Mesh& parent);

         const Mesh& getParent();

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

