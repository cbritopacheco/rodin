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

#include "Variational/ForwardDecls.h"

namespace Rodin
{
   /**
    * @brief Represents the subdivision of some domain into faces of (possibly)
    * different geometries.
    */
   class Mesh
   {
      public:
         /**
          * @brief Loads a mesh mfem text format from file.
          * @param[in] filename Name of file to read
          */
         static Mesh load(const std::filesystem::path& filename);

         /**
          * @brief Saves the mesh to file in mfem text format.
          * @param[in] filename Name of file to write
          */
         void save(const std::filesystem::path& filename);

         /**
          * @brief Constructs an empty mesh (no elements).
          */
         Mesh() = default;

         Mesh(Mesh&& other) = default;

         Mesh(const Mesh& other);

         Mesh& operator=(Mesh&& other) = default;

         /**
          * @brief Gets the dimension of the elements.
          * @returns Dimension of the mesh.
          */
         int getDimension() const;

         /**
          * @brief Displaces the mesh nodes by the displacement @f$ u @f$.
          * @param[in] u Displacement at each node
          *
          * Given a grid function @f$ u @f$, the method will perform the
          * displacement
          * @f[
          *    x \mapsto x + u(x)
          * @f]
          * at each node @f$ x @f$ of the mesh.
          *
          * @note The range dimension of @f$ u @f$ must be equal to the
          * mesh dimension.
          *
          * @returns Reference to this (for method chaining)
          */
         Mesh& displace(const Variational::GridFunctionBase& u);

         /**
          * @brief Gets the maximum number @f$ t @f$ by which the mesh will
          * remain valid, when displacing by @f$ u @f$.
          * @param[in] u Displacement at each node
          *
          * This function will calculate the maximum number @f$ t @f$ so that
          * the displacement
          * @f[
          *    x \mapsto x + t u(x)
          * @f]
          * gives a valid mesh without actually displacing the mesh.
          *
          * @note The range dimension of @f$ u @f$ must be equal to the
          * mesh dimension.
          *
          * @returns Maximum time so that the mesh remains valid.
          */
         double getMaximumDisplacement(const Variational::GridFunctionBase& u);

         double getVolume();

         mfem::Mesh& getHandle();

         const mfem::Mesh& getHandle() const;

         /**
          * @internal
          */
         explicit
         Mesh(const mfem::Mesh& mesh);
      private:
         mfem::Mesh m_mesh;
   };
}

#endif
