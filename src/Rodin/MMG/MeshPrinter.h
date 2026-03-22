/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file MeshPrinter.h
 * @brief MEDIT mesh printer for @ref Rodin::MMG::Mesh.
 */
#ifndef RODIN_EXTERNAL_MMG_MESHPRINTER_H
#define RODIN_EXTERNAL_MMG_MESHPRINTER_H

#include "Rodin/IO/MEDIT.h"

#include "Mesh.h"

namespace Rodin::MMG
{
  /**
   * @brief Printer for @ref Rodin::MMG::Mesh in MEDIT/MMG-compatible syntax.
   *
   * In addition to standard mesh topology/geometry, this printer serializes
   * MMG-specific sections:
   * - `Corners`
   * - `Ridges`
   * - `RequiredVertices`
   * - `RequiredEdges`
   */
  class MeshPrinter : public IO::Printer<MMG::Mesh>
  {
    public:
      /**
       * @brief Constructs a printer for a given MMG mesh.
       * @param[in] mesh Mesh to print.
       */
      MeshPrinter(const MMG::Mesh& mesh)
        : m_mesh(mesh),
          m_printer(mesh)
      {}

      /**
       * @brief Prints a full MEDIT mesh including MMG-specific sections.
       * @param[out] os Output stream.
       */
      void print(std::ostream& os) override;

      /**
       * @brief Prints the `Corners` section.
       * @param[out] os Output stream.
       */
      void printCorners(std::ostream& os);
      /**
       * @brief Prints the `Ridges` section.
       * @param[out] os Output stream.
       */
      void printRidges(std::ostream& os);
      /**
       * @brief Prints the `RequiredVertices` section.
       * @param[out] os Output stream.
       */
      void printRequiredVertices(std::ostream& os);
      /**
       * @brief Prints the `RequiredEdges` section.
       * @param[out] os Output stream.
       */
      void printRequiredEdges(std::ostream& os);

      /**
       * @brief Gets the mesh bound to this printer.
       * @returns The mesh object being serialized.
       */
      const MMG::Mesh& getObject() const override
      {
        return m_mesh.get();
      }

    private:
      std::reference_wrapper<const MMG::Mesh> m_mesh; ///< Source mesh reference.
      IO::MeshPrinter<IO::FileFormat::MEDIT, Context::Local> m_printer; ///< Base MEDIT printer.
  };
}

#endif
