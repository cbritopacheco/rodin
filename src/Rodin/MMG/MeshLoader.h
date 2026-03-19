/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file MeshLoader.h
 * @brief MEDIT mesh loader specialization for @ref Rodin::MMG::Mesh.
 */
#ifndef RODIN_EXTERNAL_MMG_MESHLOADER_H
#define RODIN_EXTERNAL_MMG_MESHLOADER_H

#include "Rodin/IO/MEDIT.h"

#include "Mesh.h"

namespace Rodin::MMG
{
  /**
   * @brief Loads an @ref Rodin::MMG::Mesh from the MEDIT format.
   *
   * Extends the base MEDIT loader by additionally parsing MMG-specific sections
   * (corners, ridges, required vertices and required edges) and storing them in
   * the associated @ref MMG::Mesh metadata sets.
   */
  class MeshLoader : public IO::MeshLoader<IO::FileFormat::MEDIT, Context::Local>
  {
   public:
      /// Base MEDIT loader type.
      using Parent = IO::MeshLoader<IO::FileFormat::MEDIT, Context::Local>;

      /**
       * @brief Constructs a loader bound to an MMG mesh instance.
       * @param[in, out] mesh Destination mesh populated by @ref load.
       */
      MeshLoader(MMG::Mesh& mesh)
        : Parent(mesh),
          m_mesh(mesh)
      {}

      /**
       * @brief Parses mesh content from a stream in MEDIT syntax.
       * @param[in, out] is Input stream containing a MEDIT mesh.
       */
      void load(std::istream& is) override;

      /**
       * @brief Gets the destination mesh object.
       * @returns The mesh instance currently bound to this loader.
       */
      MMG::Mesh& getObject() override
      {
        return m_mesh.get();
     }

    private:
      std::reference_wrapper<MMG::Mesh> m_mesh; ///< Destination mesh reference.
  };
}

#endif
