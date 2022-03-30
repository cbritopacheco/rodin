/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_RODININTEGRATION_MMG_MESH_H
#define RODIN_RODININTEGRATION_MMG_MESH_H

#include <cassert>
#include <boost/filesystem.hpp>

#include <mmg/libmmg.h>

#include "Configure.h"

namespace Rodin::External::MMG
{
   /**
    * @brief Base class for MMG mesh objects.
    * @tparam Dimension Mesh dimension
    *
    * This class wraps the functionality of the type MMG5_Mesh.
    */
   class MeshBase
   {
      public:
         virtual ~MeshBase() = default;

         /**
          * @brief Gets the maximum memory available to the mesh.
          */
         size_t getMaximumMemory() const
         {
            return getHandle()->memMax;
         }

         /**
          * @brief Gets the current memory usage of the mesh.
          */
         size_t getCurrentMemory() const
         {
            return getHandle()->memCur;
         }

         /**
          * @returns Dimension of the mesh.
          */
         int getDimension() const
         {
            return getHandle()->dim;
         }

         /**
          * @brief Writes the Mesh object to a text file using mmg format
          *
          * @param[in] filename Name of the file the mesh will be written to
          */
         virtual void save(const boost::filesystem::path& filename) = 0;

         /**
          * @internal
          * @returns Reference to underlying mesh handle.
          */
         virtual MMG5_pMesh& getHandle() = 0;

         /**
          * @internal
          * @returns Reference to underlying mesh handle.
          */
         virtual const MMG5_pMesh& getHandle() const = 0;
   };
}

#endif
