/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_MESH_H
#define RODIN_GEOMETRY_MESH_H

#include <set>
#include <string>
#include <deque>

#include <boost/filesystem.hpp>

#include "Rodin/Math.h"
#include "Rodin/Types.h"
#include "Rodin/IO/ForwardDecls.h"
#include "Rodin/Utility/IsSpecialization.h"
#include "Rodin/Variational/Traits.h"
#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/P1/ForwardDecls.h"

#include "ForwardDecls.h"
#include "Connectivity.h"
#include "Polytope.h"
#include "PolytopeCount.h"
#include "PolytopeIndexed.h"
#include "PolytopeIterator.h"
#include "PolytopeTransformation.h"

/**
 * @ingroup RodinDirectives
 * @brief Requires the precondition that the Mesh object have the specified
 * connectivity computed.
 *
 * Throws an exception if the current mesh instance does not have the
 * connectivity
 * @f[
 *  d \longrightarrow d', \quad 0 \leq d, d' \leq D
 * @f]
 * computed, where @f$ D @f$ is the topological dimension of the mesh.
 */
#define RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, d, dp) \
  if (mesh.getConnectivity().getIncidence(d, dp).size() == 0) \
  { \
    Rodin::Alert::MemberFunctionException(*this, __func__) \
      << Rodin::Alert::Notation::Incidence(d, dp) \
      << " has not been computed and is required to use this function." \
      << Rodin::Alert::Raise; \
  }

/**
 * @ingroup RodinDirectives
 * @brief Requires the precondition that the Mesh object have the specified
 * connectivity computed.
 *
 * Throws an exception if the current mesh instance does not have the
 * connectivity
 * @f[
 *  d \longrightarrow d', \quad 0 \leq d, d' \leq D
 * @f]
 * computed, where @f$ D @f$ is the topological dimension of the mesh.
 */
#define RODIN_GEOMETRY_MESH_REQUIRE_INCIDENCE(d, dp) \
  if (this->getConnectivity().getIncidence(d, dp).size() == 0) \
  { \
    Rodin::Alert::MemberFunctionException(*this, __func__) \
      << Rodin::Alert::Notation::Incidence(d, dp) \
      << " has not been computed and is required to use this function." \
      << Rodin::Alert::Raise; \
  }

/**
 * @ingroup RodinDirectives
 * @brief Requires the precondition that the Mesh object be a SubMesh.
 *
 * Throws an exception if the current Mesh object is not a SubMesh, i.e.
 * `this->isSubMesh()` evaluates to `false`.
 */
#define RODIN_GEOMETRY_MESH_REQUIRE_SUBMESH() \
  if (!this->isSubMesh()) \
  { \
    Rodin::Alert::MemberFunctionException(*this, __func__) \
      << "This instance of Mesh is not a SubMesh " \
      << Rodin::Alert::Notation::Predicate(false, "isSubMesh()") \
      << ". Downcasting to SubMesh is ill-defined." \
      << Rodin::Alert::Raise; \
  }

namespace Rodin::Geometry
{
  class CCL
  {
    public:
      using Component = FlatSet<Index>;

      CCL(std::deque<Component>&& dq)
        : m_components(std::move(dq))
      {}

      inline
      const std::deque<Component>& getComponents() const
      {
        return m_components;
      }

      inline
      auto begin()
      {
        return m_components.begin();
      }

      inline
      auto end()
      {
        return m_components.end();
      }

      inline
      auto begin() const
      {
        return m_components.begin();
      }

      inline
      auto end() const
      {
        return m_components.end();
      }

      inline
      auto cbegin() const
      {
        return m_components.cbegin();
      }

      inline
      auto cend() const
      {
        return m_components.cend();
      }

    private:
      std::deque<Component> m_components;
  };

  /**
   * @defgroup MeshTypes Mesh Types and Template Specializations
   * @brief Different types of mesh and template specializations of the Mesh
   * class.
   * @see Mesh
   */

  /**
  * @brief Abstract base class for Mesh objects.
  */
  class MeshBase
  {
    public:
      virtual ~MeshBase() = default;

      inline
      constexpr
      bool operator==(const MeshBase& other) const
      {
        return this == &other;
      }

      inline
      constexpr
      bool operator!=(const MeshBase& other) const
      {
        return this != &other;
      }

      inline
      CCL ccl(std::function<Boolean(const Polytope&, const Polytope&)> p) const
      {
        return ccl(p, getDimension());
      }

      inline
      CCL ccl(std::function<Boolean(const Polytope&, const Polytope&)> p, size_t d) const
      {
        return ccl(p, d, FlatSet<Attribute>{});
      }

      inline
      CCL ccl(std::function<Boolean(const Polytope&, const Polytope&)> p,
          size_t d, Attribute attr) const
      {
        return ccl(p, d, FlatSet<Attribute>{ attr });
      }

      virtual CCL ccl(std::function<Boolean(const Polytope&, const Polytope&)> p,
          size_t d,
          const FlatSet<Attribute>& attrs) const;

      virtual MeshBase& scale(Scalar c) = 0;

      virtual MeshBase& load(
        const boost::filesystem::path& filename,
        IO::FileFormat fmt = IO::FileFormat::MFEM) = 0;

      virtual void save(
        const boost::filesystem::path& filename,
        IO::FileFormat fmt = IO::FileFormat::MFEM, size_t precison = 16) const = 0;

      virtual void flush() = 0;

      /**
       * @brief Indicates if the mesh is empty or not.
       *
       * An empty mesh is defined as a mesh with no vertices.
       */
      inline
      bool isEmpty() const
      {
        return getVertexCount() == 0;
      }

      /**
       * @brief Indicates whether the mesh is a surface or not.
       *
       * A mesh is considered a surface mesh if has codimension of 1, meaning
       * the difference between its space dimension and its topological
       * dimension is 1.
       *
       * @returns True if mesh is a surface, false otherwise.
       */
      bool isSurface() const;

      /**
       * @brief Indicates whether the mesh is a submesh or not.
       * @returns True if mesh is a submesh, false otherwise.
       *
       * A Mesh which is also a SubMesh may be casted into down to access
       * the SubMesh functionality. For example:
       * @code{.cpp}
       * if (mesh.isSubMesh())
       * {
       *   // Cast is well defined
       *   auto& submesh = static_cast<SubMesh&>(mesh);
       * }
       * @endcode
       *
       */
      virtual bool isSubMesh() const = 0;

      virtual bool isInterface(Index faceIdx) const = 0;

      /**
       * @brief Determines whether a face of the mesh is on the boundary.
       */
      virtual bool isBoundary(Index faceIdx) const = 0;

      /**
       * @brief Gets the dimension of the cells.
       * @returns Dimension of the cells.
       * @see getSpaceDimension() const
       */
      virtual size_t getDimension() const = 0;

      /**
       * @brief Gets the dimension of the ambient space
       * @returns Dimension of the space which the mesh is embedded in
       * @see getDimension() const
       */
      virtual size_t getSpaceDimension() const = 0;

      /**
       * @brief Gets the number of vertices in the mesh.
       */
      inline
      size_t getVertexCount() const
      {
        return getCount(0);
      }

      /**
       * @brief Gets the number of faces in the mesh.
       */
      inline
      size_t getFaceCount() const
      {
        return getCount(getDimension() - 1);
      }

      /**
       * @brief Gets the number of cells in the mesh.
       */
      inline
      size_t getCellCount() const
      {
        return getCount(getDimension());
      }

      inline
      Attribute getFaceAttribute(Index index) const
      {
        return getAttribute(getDimension() - 1, index);
      }

      inline
      Attribute getCellAttribute(Index index) const
      {
        return getAttribute(getDimension(), index);
      }

      /**
       * @brief Gets the total volume of the mesh.
       * @returns Sum of all cell volumes.
       */
      Scalar getVolume();

      /**
       * @brief Gets the sum of the volumes of the cells given by the
       * specified attribute.
       * @param[in] attr Attribute of cells
       * @returns Sum of element volumes with given attribute
       * @note If the element attribute does not exist then this function
       * will return 0 as the volume.
       */
      Scalar getVolume(Attribute attr);

      /**
       * @brief Gets the total perimeter of the mesh.
       * @returns Sum of all element perimeters.
       */
      Scalar getPerimeter();

      /**
       * @brief Gets the sum of the perimeters of the cells given by the
       * specified attribute.
       * @param[in] attr Attribute of cells
       * @returns Sum of element perimeters with given attribute
       * @note If the element attribute does not exist then this function
       * will return 0 as the perimeter.
       */
      Scalar getPerimeter(Attribute attr);

      /**
       * @brief Gets the labels of the domain cells in the mesh.
       * @returns Set of all the attributes in the mesh object.
       * @see getBoundaryAttributes() const
       */
      virtual const FlatSet<Attribute>& getAttributes(size_t d) const = 0;

      /**
       * @brief Gets a FaceIterator for the boundary faces.
       */
      virtual FaceIterator getBoundary() const = 0;

      /**
       * @brief Gets a FaceIterator for the interface faces.
       */
      virtual FaceIterator getInterface() const = 0;

      /**
       * @brief Gets the count of polytope of the given dimension.
       * @param[in] dimension Polytope dimension
       */
      virtual size_t getCount(size_t dimension) const = 0;

      /**
       * @brief Gets the count of polytope of the given type.
       * @param[in] dim Polytope type
       */
      virtual size_t getCount(Polytope::Type g) const = 0;

      /**
       * @brief Gets an CellIterator to the cells of the mesh.
       */
      virtual CellIterator getCell(Index idx = 0) const = 0;

      /**
       * @brief Gets a FaceIterator to the faces of the mesh.
       */
      virtual FaceIterator getFace(Index idx = 0) const = 0;

      /**
       * @brief Gets a VertexIterator to the vertices of the mesh.
       */
      virtual VertexIterator getVertex(Index idx = 0) const = 0;

      /**
       * @brief Gets a PolytopeIterator to the polytopes of the given dimension
       * of the mesh.
       * @param[in] dimension Polytope dimension
       */
      virtual PolytopeIterator getPolytope(size_t dimension, Index idx = 0) const = 0;

      /**
       * @brief Gets the PolytopeTransformation associated to the @f$ (d, i)
       * @f$-polytope.
       * @param[in] d Polytope dimension
       * @param[in] idx Polytope index
       */
      virtual const PolytopeTransformation& getPolytopeTransformation(
          size_t dimension, Index idx) const = 0;

      /**
       * Gets the geometry type of the @f$ (d, i) @f$-polytope.
       * @param[in] d Polytope dimension
       * @param[in] idx Polytope index
       */
      virtual Polytope::Type getGeometry(size_t dimension, Index idx) const = 0;

      /**
       * Gets the attribute of the @f$ (d, i) @f$-polytope.
       * @param[in] d Polytope dimension
       * @param[in] idx Polytope index
       */
      virtual Attribute getAttribute(size_t dimension, Index index) const = 0;

      /**
       * Sets the attribute of the @f$ (d, i) @f$-polytope.
       * @param[in] p Pair indicating polytope dimension and index
       * @param[in] attr Attribute of polytope
       */
      virtual MeshBase& setAttribute(const std::pair<size_t, Index>& p, Attribute attr) = 0;

      /**
       * @brief Gets a reference to the mesh connectivity.
       */
      virtual MeshConnectivity& getConnectivity() = 0;

      /**
       * @brief Gets a constant reference to the mesh connectivity.
       */
      virtual const MeshConnectivity& getConnectivity() const = 0;

      /**
       * @brief Gets the space coordinates of the vertex at the given index.
       * @param[in] idx Vertex index
       */
      virtual Eigen::Map<const Math::SpatialVector> getVertexCoordinates(Index idx) const = 0;

      /**
       * @brief Sets the space coordinate of the vertex at the given index for
       * the given coordinate index.
       * @param[in] idx Vertex index
       * @param[in] s New coordinate
       * @param[in] i Coordinate index
       *
       * For example, the following code sets the coordinates of the 0-vertex
       * to @f$ (0, 5, 10) @f$ in a mesh embedded in three dimensional space.
       *
       * @code{.cpp}
       *   Mesh mesh;
       *   // Add vertices...
       *   mesh.setVertexCoordinates(0, 0.0, 0);
       *   mesh.setVertexCoordinates(0, 5.0, 1);
       *   mesh.setVertexCoordinates(0, 10.0, 2);
       * @endcode
       */
      virtual MeshBase& setVertexCoordinates(Index idx, Scalar s, size_t i) = 0;

      /**
       * @brief Sets the space coordinate of the vertex at the given index for
       * the given coordinate index.
       * @param[in] idx Vertex index
       * @param[in] coords New coordinates
       */
      virtual MeshBase& setVertexCoordinates(Index idx, const Math::SpatialVector& coords) = 0;

      virtual SubMeshBase& asSubMesh() = 0;

      virtual const SubMeshBase& asSubMesh() const = 0;
  };

  /// Index containing the indices of boundary cells.
  using BoundaryIndex = IndexSet;

  /// Index containing the attribute numbers of the polytopes.
  using AttributeIndex = PolytopeIndexed<Geometry::Attribute>;

  /// Index containing the transformations of the polytopes.
  using TransformationIndex = PolytopeIndexed<std::unique_ptr<PolytopeTransformation>>;

  /**
   *
   * @ingroup MeshTypes
   * @brief Represents the subdivision of some domain into faces of (possibly)
   * different geometries.
   */
  template <>
  class Mesh<Context::Serial> : public MeshBase
  {
    public:
      /**
       * @brief Class used to build Mesh<Context::Serial> instances.
       */
      class Builder
      {
        public:
          /**
           * @brief Default constructor.
           */
          Builder() = default;

          virtual ~Builder() = default;

          /**
           * @brief Deleted copy constructor.
           */
          Builder(const Builder&) = delete;

          /**
           * @brief Move constructor.
           */
          Builder(Builder&&) = default;

          /**
           * @brief Deleted copy assignment.
           */
          Builder& operator=(const Builder&) = delete;

          /**
           * @brief Move assignment.
           */
          Builder& operator=(Builder&&) = default;

          /**
           * @brief Reserves memory accross the data structure for the polytopes of
           * the given dimension.
           */
          Builder& reserve(size_t d, size_t count);

          /**
           * @brief Initializes construction of the mesh object.
           */
          Builder& initialize(size_t sdim);

          /**
           * @brief Sets the number of nodes in the mesh.
           */
          Builder& nodes(size_t n);

          /**
           * @brief Adds vertex with coordinates given by the fixed size array.
           *
           * @note This method requires nodes(size_t) to be called beforehand.
           */
          template <size_t Size>
          inline
          Builder& vertex(const Scalar (&data)[Size])
          {
            assert(Size == m_sdim);
            m_vertices.col(m_nodes++) = data;
            return *this;
          }

          /**
           * @brief Adds vertex with coordinates given by the initializer list.
           *
           * @note This method requires nodes(size_t) to be called beforehand.
           */
          Builder& vertex(std::initializer_list<Scalar> l);

          /**
           * @brief Adds vertex with coordinates given by the array pointer.
           *
           * @note This method requires nodes(size_t) to be called beforehand.
           */
          Builder& vertex(const Scalar* data);

          /**
           * @brief Adds vertex with coordinates given by the mapped memory.
           *
           * @note This method requires nodes(size_t) to be called beforehand.
           */
          Builder& vertex(const Eigen::Map<const Math::Vector>& x);

          /**
           * @brief Adds vertex with coordinates given by the vector.
           *
           * @note This method requires nodes(size_t) to be called beforehand.
           */
          Builder& vertex(Math::Vector&& x);

          /**
           * @brief Adds vertex with coordinates given by the vector.
           *
           * This method requires nodes(size_t) to be called beforehand.
           */
          Builder& vertex(const Math::Vector& x);

          /**
           * @brief Sets the attribute of the given polytope.
           */
          Builder& attribute(const std::pair<size_t, Index>& p, Attribute attr);

          /**
           * @brief Adds polytope defined by the given vertices.
           */
          Builder& polytope(Polytope::Type t, std::initializer_list<Index> vs)
          {
            return polytope(t, IndexArray({ vs }));
          }

          /**
           * @brief Adds polytope defined by the given vertices.
           */
          Builder& polytope(Polytope::Type t, const IndexArray& vs);

          /**
           * @brief Adds polytope defined by the given vertices.
           */
          Builder& polytope(Polytope::Type t, IndexArray&& vs);

          /**
           * @brief Finalizes construction of the Mesh<Context::Serial> object.
           */
          Mesh finalize();

          Builder& setVertices(Math::Matrix&& connectivity);

          Builder& setConnectivity(MeshConnectivity&& connectivity);

          Builder& setAttributeIndex(AttributeIndex&& connectivity);

          Builder& setTransformationIndex(TransformationIndex&& connectivity);

          inline
          MeshConnectivity& getConnectivity()
          {
            return m_connectivity;
          }

          inline
          const MeshConnectivity& getConnectivity() const
          {
            return m_connectivity;
          }

        private:
          size_t m_sdim;
          size_t m_nodes;

          Math::PointMatrix m_vertices;
          MeshConnectivity m_connectivity;

          AttributeIndex m_attributeIndex;
          TransformationIndex m_transformationIndex;

          std::vector<FlatSet<Attribute>> m_attributes;
      };

      /**
       * @brief Generates a Builder instance to build a Mesh object.
       */
      inline
      static Builder Build()
      {
        return Builder();
      }

      /**
       * @brief Generates a uniform grid for a given geometry.
       */
      static Mesh UniformGrid(Polytope::Type g, size_t h, size_t w);

      /**
      * @brief Constructs an empty mesh with no cells.
      */
      Mesh()
        : m_sdim(0)
      {}

      Mesh(const boost::filesystem::path& filename, IO::FileFormat fmt = IO::FileFormat::MFEM)
      {
        load(filename, fmt);
      }

      /**
      * @brief Performs a copy of another mesh.
      */
      Mesh(const Mesh& other)
        : m_sdim(other.m_sdim),
          m_vertices(other.m_vertices),
          m_connectivity(other.m_connectivity),
          m_attributeIndex(other.m_attributeIndex),
          m_attributes(other.m_attributes)
      {}

      /**
      * @brief Move constructs the mesh from another mesh.
      */
      Mesh(Mesh&& other) = default;

      Mesh& operator=(const Mesh& other) = delete;

      /**
      * @brief Move assigns the mesh from another mesh.
      */
      Mesh& operator=(Mesh&&) = default;

      /**
       * @brief Displaces the mesh nodes by the displacement @f$ u @f$.
       * @param[in] u Displacement at each node
       *
       * Given a vector valued function @f$ \vec{u} @f$, the method will perform the
       * displacement
       * @f[
       *   x \mapsto x + \vec{u}(x)
       * @f]
       * at each node @f$ x @f$ of the mesh.
       *
       * @note The vector dimension of @f$ \vec{u} @f$ must be equal to the
       * space dimension.
       *
       * @returns Reference to this (for method chaining)
       */
      template <class FunctionDerived>
      Mesh& displace(const Variational::FunctionBase<FunctionDerived>& u)
      {
        for (auto it = getVertex(); !it.end(); ++it)
        {
          const Geometry::Point p(*it, it->getTransformation(),
              Polytope::getVertices(Polytope::Type::Point).col(0), it->getCoordinates());
          m_vertices.col(it->getIndex()) += u(p);
        }
        return *this;
      }

      const PolytopeIndexed<Attribute>& getAttributeIndex() const
      {
        return m_attributeIndex;
      }

      const PolytopeIndexed<std::unique_ptr<PolytopeTransformation>>& getTransformationIndex() const
      {
        return m_transformationIndex;
      }

      inline
      const Math::PointMatrix& getVertices() const
      {
        return m_vertices;
      }

      /**
      * @brief Skins the mesh to obtain its boundary mesh
      * @returns SubMesh object to the boundary region of the mesh
      *
      * This function "skins" the mesh to return the boundary as a new SubMesh
      * object. The resulting mesh will be embedded in the original space
      * dimension.
      *
      * The lower dimensional polytopes of dimension @f$ 1 \leq d \leq D - 2
      * @f$ are included if the connectivity @f$ (D - 1) \longrightarrow d @f$
      * is already computed in the mesh.
      */
      virtual SubMesh<Context::Serial> skin() const;

      /**
      * @brief Trims the cells with the given attribute.
      * @param[in] attr Attribute to trim
      * @returns SubMesh of the remaining region mesh
      *
      * Convenience function to call trim(const std::FlatSet<Attribute>&) with
      * only one attribute.
      */
      virtual SubMesh<Context::Serial> trim(Attribute attr) const;

      /**
      * @brief Trims the cells with the given attribute.
      * @param[in] attrs Attributes to trim
      * @returns SubMesh object to the remaining region of the mesh
      *
      * This function will trim discard all the cells that have an attribute
      * in the given set of attributes.
      *
      * The lower dimensional polytopes of dimension @f$ 1 \leq d \leq D - 1
      * @f$ are included if the connectivity @f$ D \longrightarrow d @f$ is
      * already computed in the mesh.
      *
      * @returns A SubMesh object consisting of cells that have attributes
      * not in the given set.
      */
      virtual SubMesh<Context::Serial> trim(const FlatSet<Attribute>& attrs) const;

      /**
      * @brief Keeps the cells with the given attribute.
      * @param[in] attr Attribute to keep
      * @returns SubMesh of the remaining region mesh
      *
      * Convenience function to call keep(const std::FlatSet<Attribute>&) with
      * only one attribute.
      */
      virtual SubMesh<Context::Serial> keep(Attribute attr) const;

      /**
      * @brief Trims the cells with the given attributes.
      * @param[in] attrs Attributes to trim
      * @returns SubMesh object to the remaining region of the mesh
      *
      * This function will trim keep only the cells that have an attribute
      * in the given set of attributes.
      *
      * The lower dimensional polytopes of dimension @f$ 1 \leq d \leq D - 1
      * @f$ are included if the connectivity @f$ D \longrightarrow d @f$ is
      * already computed in the mesh.
      *
      * @returns A SubMesh object consisting of cells that have attributes
      * not in the given set.
      */
      virtual SubMesh<Context::Serial> keep(const FlatSet<Attribute>& attrs) const;

      virtual Mesh& trace(const Map<std::pair<Attribute, Attribute>, Attribute>& tmap);

      SubMeshBase& asSubMesh() override;

      const SubMeshBase& asSubMesh() const override;

      /**
      * @brief Loads a mesh from file in the given format.
      * @param[in] filename Name of file to read
      * @param[in] fmt Mesh file format
      * @returns Reference to this (for method chaining)
      */
      virtual Mesh& load(
        const boost::filesystem::path& filename,
        IO::FileFormat fmt = IO::FileFormat::MFEM) override;

      /**
      * @brief Saves a mesh to file in the given format.
      * @param[in] filename Name of file to write
      * @param[in] fmt Mesh file format
      * @returns Reference to this (for method chaining)
      */
      virtual void save(
        const boost::filesystem::path& filename,
        IO::FileFormat fmt = IO::FileFormat::MFEM, size_t precison = 16) const override;

      virtual Mesh& scale(Scalar c) override;

      virtual Mesh& setAttribute(const std::pair<size_t, Index>&, Attribute attr) override;

      virtual size_t getCount(size_t dim) const override;

      virtual size_t getCount(Polytope::Type g) const override;

      virtual FaceIterator getBoundary() const override;

      virtual FaceIterator getInterface() const override;

      virtual CellIterator getCell(Index idx = 0) const override;

      virtual FaceIterator getFace(Index idx = 0) const override;

      virtual VertexIterator getVertex(Index idx = 0) const override;

      virtual PolytopeIterator getPolytope(size_t dimension, Index idx = 0) const override;

      virtual bool isSubMesh() const override
      {
        return false;
      }

      virtual bool isInterface(Index faceIdx) const override;

      virtual bool isBoundary(Index faceIdx) const override;

      virtual size_t getDimension() const override;

      virtual size_t getSpaceDimension() const override;

      virtual const PolytopeTransformation& getPolytopeTransformation(
          size_t dimension, Index idx) const override;

      virtual Polytope::Type getGeometry(size_t dimension, Index idx) const override;

      virtual Attribute getAttribute(size_t dimension, Index index) const override;

      virtual MeshConnectivity& getConnectivity() override
      {
        return m_connectivity;
      }

      virtual const MeshConnectivity& getConnectivity() const override
      {
        return m_connectivity;
      }

      virtual void flush() override
      {
        m_transformationIndex.clear();
      }

      virtual Eigen::Map<const Math::SpatialVector> getVertexCoordinates(Index idx) const override;

      virtual const FlatSet<Attribute>& getAttributes(size_t d) const override;

      virtual Mesh& setVertexCoordinates(Index idx, Scalar xi, size_t i) override;

      virtual Mesh& setVertexCoordinates(Index idx, const Math::SpatialVector& coords) override;

    private:
      static const GeometryIndexed<Math::Matrix> s_vertices;

      size_t m_sdim;

      Math::PointMatrix m_vertices;
      MeshConnectivity m_connectivity;

      PolytopeIndexed<Geometry::Attribute> m_attributeIndex;
      mutable PolytopeIndexed<std::unique_ptr<PolytopeTransformation>> m_transformationIndex;

      std::vector<FlatSet<Attribute>> m_attributes;
  };

  /// Type alias for Mesh<Context::Serial>
  using SerialMesh = Mesh<Context::Serial>;
}

#include "Mesh.hpp"

#endif
