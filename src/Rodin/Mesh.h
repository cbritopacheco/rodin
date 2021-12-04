/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MESH_MESH_H
#define RODIN_MESH_MESH_H

#include <string>
#include <mfem.hpp>

namespace Rodin
{
   class Mesh
   {
      public:
         static Mesh load(const std::string& filename);

         void save(const std::string& filename);

         explicit
         Mesh(const mfem::Mesh& mesh);

         Mesh() = default;

         Mesh(const Mesh& other);

         int getDimension() const;

         mfem::Mesh& getHandle();

         const mfem::Mesh& getHandle() const;

      private:
         mfem::Mesh m_mesh;
   };
}

#endif
