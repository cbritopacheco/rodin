/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Mesh.h
 * @brief MMG-aware mesh type extending Rodin local meshes.
 */
#ifndef RODIN_EXTERNAL_MMG_MESH_H
#define RODIN_EXTERNAL_MMG_MESH_H

#include "Rodin/Geometry.h"

namespace Rodin::MMG
{
  /**
   * @brief Local mesh enriched with MMG boundary tags and constraints.
   *
   * This class extends @ref Rodin::Geometry::Mesh<Context::Local> with index
   * sets required by MMG workflows:
   * - corners (`MG_CRN`),
   * - ridges (`MG_GEO`),
   * - required vertices (`MG_REQ`),
   * - required edges (`MG_REQ` on edges).
   *
   * These sets are preserved when converting to/from native MMG structures via
   * @ref MMG5::rodinToMesh and @ref MMG5::meshToRodin.
   */
  class Mesh : public Geometry::Mesh<Context::Local>
  {
    public:
      /// Parent class
      using Parent = Geometry::Mesh<Rodin::Context::Local>;

      using Context = typename Parent::Context;

      /// Index set of corner vertices in the mesh.
      using CornerIndex = IndexSet;

      /// Index set of ridge edges in the mesh.
      using RidgeIndex = IndexSet;

      /// Index set of required vertices in the mesh.
      using RequiredVertexIndex = IndexSet;

      /// Index set of required edges in the mesh.
      using RequiredEdgeIndex = IndexSet;

      /**
       * @brief Class used to build MMG::Mesh instances.
       */
      class Builder : public Parent::Builder
      {
        public:
          /**
           * @brief Default constructor.
           */
          Builder() = default;

          /**
           * @brief Deleted copy constructor.
           */
          Builder(const Builder&) = delete;

          /**
           * @brief Move constructor.
           */
          Builder(Builder&& other)
            : Parent::Builder(std::move(other)),
              m_cornerIndex(std::move(other.m_cornerIndex)),
              m_ridgeIndex(std::move(other.m_ridgeIndex)),
              m_requiredVertexIndex(std::move(other.m_requiredVertexIndex)),
              m_requiredEdgeIndex(std::move(other.m_requiredEdgeIndex))
          {}

          /**
           * @brief Move assignment.
           */
          Builder& operator=(Builder&& other);

          /**
           * @brief Marks a vertex as a corner.
           * @param[in] vertexIdx Vertex index in the mesh.
           * @returns Reference to this builder.
           */
          Builder& corner(Index vertexIdx);

          /**
           * @brief Marks an edge as a ridge.
           * @param[in] edgeIdx Edge index in the mesh.
           * @returns Reference to this builder.
           */
          Builder& ridge(Index edgeIdx);

          /**
           * @brief Marks an edge as required.
           * @param[in] edgeIdx Edge index in the mesh.
           * @returns Reference to this builder.
           */
          Builder& requiredEdge(Index edgeIdx);

          /**
           * @brief Marks a vertex as required.
           * @param[in] vertexIdx Vertex index in the mesh.
           * @returns Reference to this builder.
           */
          Builder& requiredVertex(Index vertexIdx);

          /**
           * @brief Finalizes and returns the constructed MMG mesh.
           */
          Mesh finalize();

        private:
          CornerIndex m_cornerIndex;
          RidgeIndex  m_ridgeIndex;
          RequiredVertexIndex m_requiredVertexIndex;
          RequiredEdgeIndex m_requiredEdgeIndex;
      };

      /**
       * @brief Creates a builder for constructing an @ref MMG::Mesh instance.
       */
      static MMG::Mesh::Builder Build()
      {
        return MMG::Mesh::Builder();
      }

      /**
       * @brief Constructs an empty MMG mesh.
       */
      Mesh() = default;

      /**
       * @brief Move-constructs from a base local mesh.
       */
      Mesh(Parent&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Copy constructor.
       */
      Mesh(const Mesh& other)
        : Parent(other),
          m_cornerIndex(other.m_cornerIndex),
          m_requiredVertexIndex(other.m_requiredVertexIndex),
          m_ridgeIndex(other.m_ridgeIndex),
          m_requiredEdgeIndex(other.m_requiredEdgeIndex)
      {}

      /**
       * @brief Move constructor.
       */
      Mesh(Mesh&& other)
        : Parent(std::move(other)),
          m_cornerIndex(std::move(other.m_cornerIndex)),
          m_requiredVertexIndex(std::move(other.m_requiredVertexIndex)),
          m_ridgeIndex(std::move(other.m_ridgeIndex)),
          m_requiredEdgeIndex(std::move(other.m_requiredEdgeIndex))
      {}

      /**
       * @brief Move assignment.
       */
      Mesh& operator=(Mesh&& other)
      {
        Parent::operator=(std::move(other));
        m_cornerIndex = std::move(other.m_cornerIndex);
        m_requiredVertexIndex = std::move(other.m_requiredVertexIndex);
        m_ridgeIndex = std::move(other.m_ridgeIndex);
        m_requiredEdgeIndex = std::move(other.m_requiredEdgeIndex);
        return *this;
      }

      Mesh& operator=(Parent&& other)
      {
        Parent::operator=(std::move(other));
        return *this;
      }

      /**
       * @brief Marks a vertex as a corner.
       * @param[in] vertexIdx Vertex index.
       * @returns Reference to this mesh.
       */
      Mesh& setCorner(Index vertexIdx);

      /**
       * @brief Marks an edge as a ridge.
       * @param[in] edgeIdx Edge index.
       * @returns Reference to this mesh.
       */
      Mesh& setRidge(Index edgeIdx);

      /**
       * @brief Marks an edge as required.
       * @param[in] edgeIdx Edge index.
       * @returns Reference to this mesh.
       */
      Mesh& setRequiredEdge(Index edgeIdx);

      /**
       * @brief Marks a vertex as required.
       * @param[in] vertexIdx Vertex index.
       * @returns Reference to this mesh.
       */
      Mesh& setRequiredVertex(Index vertexIdx);

      /**
       * @brief Gets the corner index set (const).
       */
      const CornerIndex& getCorners() const
      {
        return m_cornerIndex;
      }

      /**
       * @brief Gets the corner index set.
       */
      CornerIndex& getCorners()
      {
        return m_cornerIndex;
      }

      /**
       * @brief Gets the ridge index set (const).
       */
      const RidgeIndex& getRidges() const
      {
        return m_ridgeIndex;
      }

      /**
       * @brief Gets the ridge index set.
       */
      RidgeIndex& getRidges()
      {
        return m_ridgeIndex;
      }

      /**
       * @brief Gets the required-edge index set.
       */
      RequiredEdgeIndex& getRequiredEdges()
      {
        return m_requiredEdgeIndex;
      }

      /**
       * @brief Gets the required-edge index set (const).
       */
      const RequiredEdgeIndex& getRequiredEdges() const
      {
        return m_requiredEdgeIndex;
      }

      /**
       * @brief Gets the required-vertex index set.
       */
      RequiredVertexIndex& getRequiredVertices()
      {
        return m_requiredVertexIndex;
      }

      /**
       * @brief Gets the required-vertex index set (const).
       */
      const RequiredVertexIndex& getRequiredVertices() const
      {
        return m_requiredVertexIndex;
      }

      /**
       * @brief Saves the mesh to disk.
       * @param[in] filename Destination path.
       * @param[in] fmt Explicit file format.
       *
       * For @ref Rodin::IO::FileFormat::MEDIT this writes MMG-specific sections
       * via @ref MeshPrinter. Other formats delegate to the parent mesh
       * implementation.
       */
      void save(
         const boost::filesystem::path& filename,
         IO::FileFormat fmt) const override;

      /**
       * @brief Loads the mesh from disk.
       * @param[in] filename Source path.
       * @param[in] fmt Explicit file format.
       * @returns Reference to this mesh.
       *
       * For @ref Rodin::IO::FileFormat::MEDIT this restores MMG-specific
       * sections via @ref MeshLoader. Other formats delegate to the parent mesh
       * implementation.
       */
      Mesh& load(
         const boost::filesystem::path& filename,
         IO::FileFormat fmt) override;

    private:
      CornerIndex m_cornerIndex;
      RequiredVertexIndex m_requiredVertexIndex;

      RidgeIndex  m_ridgeIndex;
      RequiredEdgeIndex m_requiredEdgeIndex;
  };
}

#endif
