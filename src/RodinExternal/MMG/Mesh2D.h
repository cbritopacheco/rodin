/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_RODININTEGRATION_MMG_MESH2D_H
#define RODIN_RODININTEGRATION_MMG_MESH2D_H

#include <string>
#include <vector>
#include <variant>
#include <optional>
#include <filesystem>

#include <mmg/mmg2d/libmmg2d.h>

#include "Mesh.h"

namespace Rodin::External::MMG
{
   /**
    * @brief Represents a mesh of a 2D domain.
    */
   class Mesh2D : public Mesh<2>
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
          * @brief Loads the Mesh2D object from a text file, assuming it is in
          * the `medit2` file format.
          *
          * @param[in] filename Name of file where the mesh will be read from.
          * @returns A Mesh2D object containing the mesh data in the file.
          */
         static Mesh2D load(const std::filesystem::path& filename);

         /**
          * @brief Writes the Mesh2D object to file using the `medit2` file
          * format.
          *
          * @param[in] filename Name of file where the mesh will be written to.
          */
         void save(const std::filesystem::path& filename);

         /**
          * @brief Creates an empty mesh.
          */
         Mesh2D();

         /**
          * @brief Copies the mesh.
          */
         Mesh2D(const Mesh2D& other);

         /**
          * @brief Moves constructs the data of another Mesh2D object, while
          * invalidating the latter.
          */
         Mesh2D(Mesh2D&& other);

         /**
          * @brief Destructs the object.
          */
         ~Mesh2D();

         /**
          * @brief Gets a count of the specified entities in the mesh.
          * @param[in] e Type of entity to count.
          * @return Number of entities of type Mesh2D::Entity.
          */
         int count(Mesh2D::Entity e) const;

         MMG5_pMesh& getHandle() override;

         const MMG5_pMesh& getHandle() const override;

      private:
         MMG5_pMesh  m_mesh;
   };
}

#endif
