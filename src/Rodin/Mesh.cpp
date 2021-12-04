/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Mesh.h"

namespace Rodin
{
   Mesh Mesh::load(const std::string& filename)
   {
      return Mesh(mfem::Mesh::LoadFromFile(filename.c_str()));
   }

   void Mesh::save(const std::string& filename)
   {
      getHandle().Save(filename.c_str());
   }

   Mesh::Mesh(const mfem::Mesh& mesh)
      : m_mesh(mesh)
   {}

   Mesh::Mesh(const Mesh& other)
      : m_mesh(other.m_mesh)
   {}

   int Mesh::getDimension() const
   {
      return m_mesh.Dimension();
   }

   mfem::Mesh& Mesh::getHandle()
   {
      return m_mesh;
   }

   const mfem::Mesh& Mesh::getHandle() const
   {
      return m_mesh;
   }
}

