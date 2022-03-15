/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_MESH3D_H
#define RODIN_EXTERNAL_MMG_MESH3D_H

#include <string>
#include <vector>
#include <variant>
#include <optional>
#include <filesystem>

#include <mmg/mmg3d/libmmg3d.h>

#include "Mesh.h"

namespace Rodin::External::MMG
{
   /**
    * @brief Represents a mesh of a 3D domain.
    */
   class Mesh3D : public MeshBase
   {
      public:
         /**
          * @brief Enumeration describing possible mesh entities.
          */
         enum Entity
         {
            Vertex = MMG5_Vertex,
            Edge = MMG5_Edg,
            Triangle = MMG5_Triangle,
            Tetrahedra = MMG5_Tetrahedron
         };

         /**
          * @brief Loads the Mesh3D object from a text file
          *
          * @param[in] filename Name of file where the mesh will be read from.
          * @returns A Mesh3D object containing the mesh data in the file.
          */
         static Mesh3D load(const std::filesystem::path& filename);

         /**
          * @brief Creates an empty mesh.
          */
         Mesh3D();

         /**
          * @brief Performs a deep copy of the mesh.
          */
         Mesh3D(const Mesh3D& other);

         /**
          * @brief Moves constructs the data of another Mesh3D object, while
          * invalidating the latter.
          */
         Mesh3D(Mesh3D&& other);

         /**
          * @brief Destructs the object.
          */
         ~Mesh3D();

         /**
          * @brief Gets a count of the specified entities in the mesh.
          * @param[in] e Type of entity to count.
          * @return Number of entities of type Mesh3D::Entity.
          */
         int count(Mesh3D::Entity e) const;

         void save(const std::filesystem::path& filename) override;

         MMG5_pMesh& getHandle() override;

         const MMG5_pMesh& getHandle() const override;

      private:
         MMG5_pMesh  m_mesh;
   };
}

#endif
