/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_SURFACEMESH_H
#define RODIN_EXTERNAL_MMG_SURFACEMESH_H

#include <string>
#include <vector>
#include <variant>
#include <optional>
#include <filesystem>

#include <mmg/mmgs/libmmgs.h>

#include "Mesh.h"

namespace Rodin::External::MMG
{
   /**
    * @brief Triangular surface mesh
    */
   class SurfaceMesh : public MeshBase
   {
       public:
          /**
           * @brief Enumeration describing possible mesh entities.
           */
          enum Entity
          {
             Vertex = MMG5_Vertex,
             Edge = MMG5_Edg,
             Triangle = MMG5_Triangle
          };

          /**
           * @brief Loads the SurfaceMesh object from a text file
           *
           * @param[in] filename Name of file where the mesh will be read from.
           * @returns A Mesh3D object containing the mesh data in the file.
           */
          static SurfaceMesh load(const std::filesystem::path& filename);

          /**
           * @brief Creates an empty mesh.
           */
          SurfaceMesh();

          /**
           * @brief Performs a deep copy of the mesh.
           */
          SurfaceMesh(const SurfaceMesh& other);

          /**
           * @brief Moves constructs the data of another SurfaceMesh object, while
           * invalidating the latter.
           */
          SurfaceMesh(SurfaceMesh&& other);

          /**
           * @brief Destructs the object.
           */
          ~SurfaceMesh();

          /**
           * @brief Gets a count of the specified entities in the mesh.
           * @param[in] e Type of entity to count.
           * @return Number of entities of type SurfaceMesh::Entity.
           */
          int count(SurfaceMesh::Entity e) const;

          void save(const std::filesystem::path& filename) override;

          MMG5_pMesh& getHandle() override;

          const MMG5_pMesh& getHandle() const override;

       private:
          MMG5_pMesh  m_mesh;
   };
}

#endif
