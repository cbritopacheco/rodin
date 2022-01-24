/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Variational/GridFunction.h"
#include "Rodin/Variational/FiniteElementSpace.h"

#include "Mesh.h"

namespace Rodin
{
   Mesh Mesh::load(const std::filesystem::path& filename)
   {
      return Mesh(mfem::Mesh::LoadFromFile(filename.c_str()));
   }

   void Mesh::save(const std::filesystem::path& filename)
   {
      getHandle().Save(filename.c_str());
   }

   Mesh::Mesh(const mfem::Mesh& mesh)
      : m_mesh(mesh)
   {}

   Mesh::Mesh(const Mesh& other)
      : m_mesh(other.m_mesh)
   {}

   Mesh& Mesh::displace(const Variational::GridFunctionBase& u)
   {
      assert(u.getFiniteElementSpace().getRangeDimension() == getDimension());
      m_mesh.MoveNodes(u.getHandle());
      return *this;
   }

   double
   Mesh::getMaximumDisplacement(const Variational::GridFunctionBase& u)
   {
      double res;
      m_mesh.CheckDisplacements(u.getHandle(), res);
      return res;
   }

   int Mesh::getDimension() const
   {
      return m_mesh.Dimension();
   }

   double Mesh::getVolume()
   {
      double totalVolume = 0;
      for (int i = 0; i < m_mesh.GetNE(); i++)
         totalVolume += m_mesh.GetElementVolume(i);
      return totalVolume;
   }

   double Mesh::getVolume(int attr)
   {
      double totalVolume = 0;
      for (int i = 0; i < m_mesh.GetNE(); i++)
      {
         if (m_mesh.GetAttribute(i) == attr)
            totalVolume += m_mesh.GetElementVolume(i);
      }
      return totalVolume;
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

