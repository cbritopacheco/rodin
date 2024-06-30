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
#include "Rodin/Configure.h"
#include "Rodin/Threads/Mutable.h"
#include "Rodin/IO/ForwardDecls.h"
#include "Rodin/Context/Sequential.h"
#include "Rodin/Utility/IsSpecialization.h"
#include "Rodin/Variational/Traits.h"
#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/P1/ForwardDecls.h"

#include "ForwardDecls.h"
#include "Connectivity.h"
#include "Point.h"
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

namespace Rodin::FormLanguage
{
  template <class Context>
  struct Traits<Geometry::Mesh<Context>>
  {
    using ContextType = Context;
  };
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

      inline
      size_t getCount() const
      {
        return m_components.size();
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

      virtual std::optional<Point> inclusion(const Point& p) const = 0;

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

      virtual SubMeshBase& asSubMesh() = 0;

      virtual const SubMeshBase& asSubMesh() const = 0;

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
        return getPolytopeCount(0);
      }

      /**
       * @brief Gets the number of faces in the mesh.
       */
      inline
      size_t getFaceCount() const
      {
        return getPolytopeCount(getDimension() - 1);
      }

      /**
       * @brief Gets the number of cells in the mesh.
       */
      inline
      size_t getCellCount() const
      {
        return getPolytopeCount(getDimension());
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
      Scalar getVolume() const;

      /**
       * @brief Gets the sum of the volumes of the cells given by the
       * specified attribute.
       * @param[in] attr Attribute of cells
       * @returns Sum of element volumes with given attribute
       * @note If the element attribute does not exist then this function
       * will return 0 as the volume.
       */
      Scalar getVolume(Attribute attr) const;

      Scalar getVolume(const FlatSet<Attribute>& attr) const;

      /**
       * @brief Gets the total perimeter of the mesh.
       * @returns Sum of all element perimeters.
       */
      Scalar getPerimeter() const;

      /**
       * @brief Gets the sum of the perimeters of the cells given by the
       * specified attribute.
       * @param[in] attr Attribute of cells
       * @returns Sum of element perimeters with given attribute
       * @note If the element attribute does not exist then this function
       * will return 0 as the perimeter.
       */
      Scalar getPerimeter(Attribute attr) const;

      Scalar getPerimeter(const FlatSet<Attribute>& attr) const;

      Scalar getArea() const;

      Scalar getArea(Attribute attr) const;

      Scalar getArea(const FlatSet<Attribute>& attr) const;

      Scalar getMeasure(size_t d) const;

      Scalar getMeasure(size_t d, Attribute attr) const;

      Scalar getMeasure(size_t d, const FlatSet<Attribute>& attr) const;

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
      virtual size_t getPolytopeCount(size_t dimension) const = 0;

      /**
       * @brief Gets the count of polytope of the given type.
       * @param[in] dim Polytope type
       */
      virtual size_t getPolytopeCount(Polytope::Type g) const = 0;

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
      virtual const PolytopeTransformation& getPolytopeTransformation(size_t dimension, Index idx) const = 0;

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
      virtual ConnectivityBase& getConnectivity() = 0;

      /**
       * @brief Gets a constant reference to the mesh connectivity.
       */
      virtual const ConnectivityBase& getConnectivity() const = 0;

      /**
       * @brief Gets the space coordinates of the vertex at the given index.
       * @param[in] idx Vertex index
       */
      virtual Eigen::Map<const Math::SpatialVector<Scalar>> getVertexCoordinates(Index idx) const = 0;

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
      virtual MeshBase& setVertexCoordinates(Index idx, const Math::SpatialVector<Scalar>& coords) = 0;

      virtual MeshBase& setPolytopeTransformation(
          const std::pair<size_t, Index> p, PolytopeTransformation* trans) = 0;

      virtual const Context::Base& getContext() const = 0;
  };

  /// Type alias for Mesh<Context::Sequential>
  using SequentialMesh = Mesh<Context::Sequential>;

  /// Index containing the indices of boundary cells.
  using BoundaryIndex = IndexSet;

  /// Index containing the attribute numbers of the polytopes.
  using AttributeIndex = PolytopeIndexed<Geometry::Attribute>;

  /// Index containing the transformations of the polytopes.
  using TransformationIndex =
    std::vector<Threads::Mutable<std::vector<PolytopeTransformation*>>>;

  /**
   *
   * @ingroup MeshTypes
   * @brief Represents the subdivision of some domain into faces of (possibly)
   * different geometries.
   */
  template <>
  class Mesh<Context::Sequential> : public MeshBase
  {
    public:
      using Parent = MeshBase;
      using Context = Context::Sequential;

      /**
       * @brief Class used to build Mesh<Context::Sequential> instances.
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
          Builder& vertex(const Eigen::Map<const Math::Vector<Scalar>>& x);

          /**
           * @brief Adds vertex with coordinates given by the vector.
           *
           * @note This method requires nodes(size_t) to be called beforehand.
           */
          Builder& vertex(Math::Vector<Scalar>&& x);

          /**
           * @brief Adds vertex with coordinates given by the vector.
           *
           * This method requires nodes(size_t) to be called beforehand.
           */
          Builder& vertex(const Math::Vector<Scalar>& x);

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
           * @brief Finalizes construction of the Mesh<Context::Sequential> object.
           */
          Mesh finalize();

          Builder& setVertices(Math::Matrix<Scalar>&& connectivity);

          Builder& setConnectivity(Connectivity<Context>&& connectivity);

          Builder& setAttributeIndex(AttributeIndex&& connectivity);

          Builder& setTransformationIndex(TransformationIndex&& connectivity);

          inline
          Connectivity<Context>& getConnectivity()
          {
            return m_connectivity;
          }

          inline
          const Connectivity<Context>& getConnectivity() const
          {
            return m_connectivity;
          }

        private:
          size_t m_sdim;
          size_t m_nodes;

          Math::PointMatrix m_vertices;
          Connectivity<Context> m_connectivity;

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

      inline
      static Mesh UniformGrid(Polytope::Type g, std::initializer_list<size_t> l)
      {
        Array<size_t> shape(l.size());
        std::copy(l.begin(), l.end(), shape.begin());
        return UniformGrid(g, shape);
      }

      /**
       * @brief Generates a uniform grid for a given geometry.
       */
      static Mesh UniformGrid(Polytope::Type g, const Array<size_t>& shape);

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
      Mesh(const Mesh& other);

      /**
      * @brief Move constructs the mesh from another mesh.
      */
      Mesh(Mesh&& other);

      virtual ~Mesh();

      Mesh& operator=(const Mesh& other) = delete;

      /**
      * @brief Move assigns the mesh from another mesh.
      */
      Mesh& operator=(Mesh&&);

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


      virtual void flush() override
      {
        for (auto& mt : m_transformationIndex)
          mt.write([](auto& obj) { obj.clear(); });
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
      virtual SubMesh<Context> skin() const;

      /**
      * @brief Trims the cells with the given attribute.
      * @param[in] attr Attribute to trim
      * @returns SubMesh of the remaining region mesh
      *
      * Convenience function to call trim(const std::FlatSet<Attribute>&) with
      * only one attribute.
      */
      virtual SubMesh<Context> trim(Attribute attr) const;

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
      virtual SubMesh<Context> trim(const FlatSet<Attribute>& attrs) const;

      /**
      * @brief Keeps the cells with the given attribute.
      * @param[in] attr Attribute to keep
      * @returns SubMesh of the remaining region mesh
      *
      * Convenience function to call keep(const std::FlatSet<Attribute>&) with
      * only one attribute.
      */
      virtual SubMesh<Context> keep(Attribute attr) const;

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
      virtual SubMesh<Context> keep(const FlatSet<Attribute>& attrs) const;

      inline
      Mesh& trace(const Map<std::pair<Attribute, Attribute>, Attribute>& tmap)
      {
        return trace(tmap, FlatSet<Attribute>{});
      }

      inline
      Mesh& trace(const Map<std::pair<Attribute, Attribute>, Attribute>& tmap, Attribute attr)
      {
        return trace(tmap, FlatSet<Attribute>{ attr });
      }

      virtual Mesh& trace(const Map<std::pair<Attribute, Attribute>, Attribute>& tmap, const FlatSet<Attribute>& attrs);

      virtual std::optional<Point> inclusion(const Point& p) const override;

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

      const AttributeIndex& getAttributeIndex() const
      {
        return m_attributeIndex;
      }

      const TransformationIndex& getTransformationIndex() const
      {
        return m_transformationIndex;
      }

      inline
      const Math::PointMatrix& getVertices() const
      {
        return m_vertices;
      }

      inline
      const Context& getContext() const override
      {
        return m_context;
      }

      virtual size_t getPolytopeCount(size_t dim) const override;

      virtual size_t getPolytopeCount(Polytope::Type g) const override;

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

      virtual Polytope::Type getGeometry(size_t dimension, Index idx) const override;

      virtual Attribute getAttribute(size_t dimension, Index index) const override;

      virtual Connectivity<Context>& getConnectivity() override
      {
        return m_connectivity;
      }

      virtual const Connectivity<Context>& getConnectivity() const override
      {
        return m_connectivity;
      }

      virtual Eigen::Map<const Math::SpatialVector<Scalar>> getVertexCoordinates(Index idx) const override;

      virtual const FlatSet<Attribute>& getAttributes(size_t d) const override;

      virtual Mesh& setAttribute(const std::pair<size_t, Index>&, Attribute attr) override;

      virtual Mesh& setVertexCoordinates(Index idx, Scalar xi, size_t i) override;

      virtual Mesh& setVertexCoordinates(Index idx, const Math::SpatialVector<Scalar>& coords) override;

      virtual Mesh& setPolytopeTransformation(
          const std::pair<size_t, Index> p, PolytopeTransformation* trans) override;

      virtual const PolytopeTransformation& getPolytopeTransformation(
          size_t dimension, Index idx) const override;

    protected:
      PolytopeTransformation* getDefaultPolytopeTransformation(size_t d, Index i) const;

    private:
      size_t m_sdim;

      Math::PointMatrix m_vertices;
      Connectivity<Context> m_connectivity;

      AttributeIndex m_attributeIndex;
      mutable TransformationIndex m_transformationIndex;

      std::vector<FlatSet<Attribute>> m_attributes;

      Context m_context;
  };
}

#include "Mesh.hpp"

#endif
