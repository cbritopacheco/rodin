/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MESH_MESH_H
#define RODIN_MESH_MESH_H

#include <set>
#include <string>
#include <filesystem>

#include <mfem.hpp>

#include "Rodin/Variational/ForwardDecls.h"

#include "ForwardDecls.h"

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
          * @brief Loads a mesh from file.
          * @param[in] filename Name of file to read
          */
         static Mesh load(const std::filesystem::path& filename);

         /**
          * @brief Saves the mesh to file.
          * @param[in] filename Name of file to write
          */
         void save(const std::filesystem::path& filename);

         /**
          * @brief Constructs an empty mesh with no elements.
          */
         Mesh() = default;

         /**
          * @brief Move constructs the mesh from another mesh.
          */
         Mesh(Mesh&& other) = default;

         /**
          * @brief Performs a deep copy of another mesh.
          */
         Mesh(const Mesh& other);

         /**
          * @brief Move assigns the mesh from another mesh.
          */
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

         /**
          * @brief Gets the total volume of the mesh.
          */
         double getVolume();

         /**
          * @brief Gets the sum of the volumes of the elements given by the
          * specified attribute.
          * @param[in] attr Attribute of elements
          * @returns Sum of element volumes
          * @note If the element attribute does not exist then this function
          * will return 0 as the volume.
          */
         double getVolume(int attr);

         /**
          * @brief Trims the elements with the given material reference.
          *
          * This function will trim the current mesh and return a Submesh
          * object containing the elements which were not trimmed from the
          * original mesh.
          */
         SubMesh trim(int attr, int bdrLabel);

         /**
          * @brief Trims the elements with the given material references.
          *
          * This function will trim the current mesh and return a Submesh
          * object containing the elements which were not trimmed from the
          * original mesh.
          */
         SubMesh trim(const std::set<int>& attrs, int bdrLabel);

         /**
          * @brief Indicates whether the mesh is a submesh or not.
          * @returns True if mesh is a submesh, false otherwise.
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
         virtual bool isSubMesh() const
         {
            return false;
         }

         /**
          * @brief Performs a uniform refinement over all mesh elements
          */
         void refine()
         {
            getHandle().UniformRefinement();
         }

         /**
          * @internal
          * @brief Gets the underlying handle for the internal mesh.
          * @returns Reference to the underlying mfem::Mesh.
          */
         mfem::Mesh& getHandle();

         /**
          * @internal
          * @brief Gets the underlying handle for the internal mesh.
          * @returns Constant reference to the underlying mfem::Mesh.
          */
         const mfem::Mesh& getHandle() const;

         /**
          * @internal
          * @brief Move constructs a Rodin::Mesh from an mfem::Mesh.
          */
         explicit
         Mesh(mfem::Mesh&& mesh);

      private:
         mfem::Mesh m_mesh;
   };
}

#endif
