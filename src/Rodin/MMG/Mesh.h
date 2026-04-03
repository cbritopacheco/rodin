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

      /// Mesh context type (always local for @ref MMG::Mesh).
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

          // ----------------------------------------------------------------
          // Forwarding overrides of Parent::Builder methods.
          // These return MMG::Mesh::Builder& so that method chaining can
          // reach MMG-specific setters (corner, ridge, etc.) after calling
          // any inherited builder method.
          // ----------------------------------------------------------------

          /** @copydoc Parent::Builder::reserve */
          Builder& reserve(size_t d, size_t count)
          {
            Parent::Builder::reserve(d, count);
            return *this;
          }

          /** @copydoc Parent::Builder::initialize */
          Builder& initialize(size_t sdim)
          {
            Parent::Builder::initialize(sdim);
            return *this;
          }

          /** @copydoc Parent::Builder::nodes */
          Builder& nodes(size_t n)
          {
            Parent::Builder::nodes(n);
            return *this;
          }

          /** @copydoc Parent::Builder::vertex(std::initializer_list<Real>) */
          Builder& vertex(std::initializer_list<Real> l)
          {
            Parent::Builder::vertex(l);
            return *this;
          }

          /** @copydoc Parent::Builder::vertex(const Real*) */
          Builder& vertex(const Real* data)
          {
            Parent::Builder::vertex(data);
            return *this;
          }

          /** @copydoc Parent::Builder::vertex(const Math::SpatialPoint&) */
          Builder& vertex(const Math::SpatialPoint& x)
          {
            Parent::Builder::vertex(x);
            return *this;
          }

          /** @copydoc Parent::Builder::vertex(const Real (&)[Size]) */
          template <size_t Size>
          Builder& vertex(const Real (&data)[Size])
          {
            Parent::Builder::vertex(data);
            return *this;
          }

          /** @copydoc Parent::Builder::attribute */
          Builder& attribute(
              const std::pair<size_t, Index>& p,
              const Optional<Geometry::Attribute>& attr)
          {
            Parent::Builder::attribute(p, attr);
            return *this;
          }

          /** @copydoc Parent::Builder::polytope(Polytope::Type, std::initializer_list<Index>) */
          Builder& polytope(Geometry::Polytope::Type t, std::initializer_list<Index> vs)
          {
            Parent::Builder::polytope(t, vs);
            return *this;
          }

          /** @copydoc Parent::Builder::polytope(Polytope::Type, const IndexArray&) */
          Builder& polytope(Geometry::Polytope::Type t, const IndexArray& vs)
          {
            Parent::Builder::polytope(t, vs);
            return *this;
          }

          /** @copydoc Parent::Builder::polytope(Polytope::Type, IndexArray&&) */
          Builder& polytope(Geometry::Polytope::Type t, IndexArray&& vs)
          {
            Parent::Builder::polytope(t, std::move(vs));
            return *this;
          }

          // ----------------------------------------------------------------
          // MMG-specific builder methods
          // ----------------------------------------------------------------

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
           * @returns Fully built mesh containing both topology/geometry and
           * MMG metadata (corners, ridges, required entities).
           */
          Mesh finalize();

        private:
          CornerIndex m_cornerIndex; ///< Corner vertex set accumulated during build.
          RidgeIndex  m_ridgeIndex; ///< Ridge edge set accumulated during build.
          RequiredVertexIndex m_requiredVertexIndex; ///< Required-vertex set accumulated during build.
          RequiredEdgeIndex m_requiredEdgeIndex; ///< Required-edge set accumulated during build.
      };

      /**
       * @brief Creates a builder for constructing an @ref MMG::Mesh instance.
       * @returns Fresh builder object.
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

      /**
       * @brief Copy assignment is deleted.
       *
       * Parent @ref Geometry::Mesh<Context::Local> deletes copy assignment.
       */
      Mesh& operator=(const Mesh& other) = delete;

      /**
       * @brief Move-assigns from a parent mesh, clearing MMG metadata.
       */
      Mesh& operator=(Parent&& other)
      {
        Parent::operator=(std::move(other));
        m_cornerIndex.clear();
        m_requiredVertexIndex.clear();
        m_ridgeIndex.clear();
        m_requiredEdgeIndex.clear();
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
       * @returns Immutable set of corner vertex indices.
       */
      const CornerIndex& getCorners() const
      {
        return m_cornerIndex;
      }

      /**
       * @brief Gets the corner index set.
       * @returns Mutable set of corner vertex indices.
       */
      CornerIndex& getCorners()
      {
        return m_cornerIndex;
      }

      /**
       * @brief Gets the ridge index set (const).
       * @returns Immutable set of ridge edge indices.
       */
      const RidgeIndex& getRidges() const
      {
        return m_ridgeIndex;
      }

      /**
       * @brief Gets the ridge index set.
       * @returns Mutable set of ridge edge indices.
       */
      RidgeIndex& getRidges()
      {
        return m_ridgeIndex;
      }

      /**
       * @brief Gets the required-edge index set.
       * @returns Mutable set of required edge indices.
       */
      RequiredEdgeIndex& getRequiredEdges()
      {
        return m_requiredEdgeIndex;
      }

      /**
       * @brief Gets the required-edge index set (const).
       * @returns Immutable set of required edge indices.
       */
      const RequiredEdgeIndex& getRequiredEdges() const
      {
        return m_requiredEdgeIndex;
      }

      /**
       * @brief Gets the required-vertex index set.
       * @returns Mutable set of required vertex indices.
       */
      RequiredVertexIndex& getRequiredVertices()
      {
        return m_requiredVertexIndex;
      }

      /**
       * @brief Gets the required-vertex index set (const).
       * @returns Immutable set of required vertex indices.
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
      CornerIndex m_cornerIndex; ///< Corner vertex indices (`MG_CRN`).
      RequiredVertexIndex m_requiredVertexIndex; ///< Required vertex indices (`MG_REQ` on vertices).

      RidgeIndex  m_ridgeIndex; ///< Ridge edge indices (`MG_GEO`).
      RequiredEdgeIndex m_requiredEdgeIndex; ///< Required edge indices (`MG_REQ` on edges).
  };
}

#endif
