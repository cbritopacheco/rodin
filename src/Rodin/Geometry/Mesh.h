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

#include <mfem.hpp>

#include "Rodin/Configure.h"

#include <boost/filesystem.hpp>

#include "Rodin/IO/ForwardDecls.h"
#include "Rodin/Variational/ForwardDecls.h"

#include "ForwardDecls.h"
#include "SimplexIterator.h"

#ifdef RODIN_USE_MPI
#include <boost/mpi.hpp>
#endif

namespace Rodin
{
   /**
    * @brief Namespace containing utilities for handling serial and parallel contexts.
    */
   namespace Context
   {
      /**
       * @brief Serial context metadata.
       */
      struct Serial   {};

      /**
       * @brief Parallel context metadata.
       */
      struct Parallel {};
   }
}

namespace Rodin::Geometry
{
   /**
    * @brief Abstract base class for Mesh objects.
    */
   class MeshBase
   {
      public:
         virtual ~MeshBase() = default;

         virtual BoundaryIterator getBoundary() const = 0;

         virtual InterfaceIterator getInterface() const = 0;

         virtual size_t getElementCount() const = 0;

         virtual size_t getCount(size_t dim) const = 0;

         virtual ElementIterator getElement(size_t idx = 0) const = 0;

         virtual FaceIterator getFace(size_t idx = 0) const = 0;

         virtual VertexIterator getVertex(size_t idx = 0) const = 0;

         virtual SimplexIterator getSimplex(size_t dimension, Index idx) const = 0;

         virtual bool isInterface(Index faceIdx) const = 0;

         virtual bool isBoundary(Index faceIdx) const = 0;

         virtual Attribute getAttribute(size_t dimension, Index index) const = 0;

         MeshBase& update();

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
          * @note The vector dimension of @f$ u @f$ must be equal to the
          * space dimension.
          *
          * @returns Reference to this (for method chaining)
          */
         MeshBase& displace(const Variational::GridFunctionBase& u);

         /**
          * @brief Edits all elements in the mesh via the given function.
          * @param[in] f Function which takes an ElementView to modify each
          * element.
          */
         // MeshBase& edit(std::function<void(const MeshElementIterator)> f);

         // /**
         //  * @brief Performs connected-component labelling.
         //  * @param[in] p Function which returns true if two adjacent elements
         //  * belong to the same component, false otherwise.
         //  * @returns List of sets, each set representing a component containing
         //  * the indices of its elements.
         //  *
         //  * @note Both elements passed to the function will always be adjacent
         //  * to each other, i.e. it is not necessary to verify this is the case.
         //  */
         // std::deque<std::set<int>> ccl(
         //       std::function<bool(const Element&, const Element&)> p) const;

         /**
          * @brief Edits the specified elements in the mesh via the given function.
          * @param[in] f Function which takes an ElementView to modify each
          * element.
          * @param[in] elements Set of indices corresponding to the elements
          * which will be modified.
          */
         // MeshBase& edit(std::function<void(ElementView)> f, const std::set<int>& elements);

         /**
          * @brief Obtains a set of elements satisfying the given condition.
          * @param[in] condition Function which returns true if the element
          * satisfies the condition.
          * @returns Set containing the element indices satisfying the
          * condition.
          */
         // std::set<int> where(std::function<bool(const Element&)> condition) const;

         // std::set<int> where(std::function<bool(const Point&)> condition) const;

         /**
          * @brief Indicates whether the mesh is a surface or not.
          * @returns True if mesh is a surface, false otherwise.
          *
          * A surface mesh is a mesh of codimension 1. That is, the difference
          * between its space dimension and dimension is 1.
          */
         bool isSurface() const;

         /**
          * @brief Gets the total volume of the mesh.
          * @returns Sum of all element volumes.
          */
         double getVolume();

         /**
          * @brief Gets the sum of the volumes of the elements given by the
          * specified attribute.
          * @param[in] attr Attribute of elements
          * @returns Sum of element volumes with given attribute
          * @note If the element attribute does not exist then this function
          * will return 0 as the volume.
          */
         double getVolume(Attribute attr);

         /**
          * @brief Gets the total perimeter of the mesh.
          * @returns Sum of all element perimeters.
          */
         double getPerimeter();

         /**
          * @brief Gets the sum of the perimeters of the elements given by the
          * specified attribute.
          * @param[in] attr Attribute of elements
          * @returns Sum of element perimeters with given attribute
          * @note If the element attribute does not exist then this function
          * will return 0 as the perimeter.
          */
         double getPerimeter(Attribute attr);

         /**
          * @brief Gets the labels of the domain elements in the mesh.
          * @returns Set of all the attributes in the mesh object.
          * @see getBoundaryAttributes() const
          */
         std::set<int> getAttributes() const;

         /**
          * @brief Gets the labels of the boundary elements in the mesh.
          * @returns Set of all the boundary attributes in the mesh object.
          * @see getAttributes() const
          */
         std::set<int> getBoundaryAttributes() const;

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
          * @note The vector dimension of @f$ u @f$ must be equal to the
          * space dimension.
          *
          * @returns Maximum time so that the mesh remains valid.
          */
         double getMaximumDisplacement(const Variational::GridFunctionBase& u);

         bool operator==(const MeshBase& other) const
         {
            return this == &other;
         }

         bool operator!=(const MeshBase& other) const
         {
            return this != &other;
         }

         /**
          * @brief Gets the dimension of the elements.
          * @returns Dimension of the elements.
          * @see getSpaceDimension() const
          */
         virtual size_t getDimension() const;

         /**
          * @brief Gets the dimension of the ambient space
          * @returns Dimension of the space which the mesh is embedded in
          * @see getDimension() const
          */
         virtual size_t getSpaceDimension() const;

         /**
          * @brief Indicates if the MeshBase object has been parallelized.
          * @returns True if the MeshBase has been parallelized, false
          * otherwise.
          * @see Mesh<Traits::Parallel>
          * @see Mesh<Traits::Serial>::parallelize()
          */
         virtual bool isParallel() const = 0;

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
         virtual bool isSubMesh() const = 0;

         /**
          * @internal
          * @brief Gets the underlying handle for the internal mesh.
          * @returns Reference to the underlying mfem::Mesh.
          */
         virtual mfem::Mesh& getHandle() = 0;

         /**
          * @internal
          * @brief Gets the underlying handle for the internal mesh.
          * @returns Constant reference to the underlying mfem::Mesh.
          */
         virtual const mfem::Mesh& getHandle() const = 0;
   };

   using SerialMesh = Mesh<Context::Serial>;

   /**
    * @brief Represents the subdivision of some domain into faces of (possibly)
    * different geometries.
    */
   template <>
   class Mesh<Context::Serial> : public MeshBase
   {
      public:
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
         Mesh(const Mesh& other) = default;

         /**
          * @brief Move assigns the mesh from another mesh.
          */
         Mesh& operator=(Mesh&& other) = default;

         virtual Mesh& initialize(size_t dim, size_t sdim);

         virtual Mesh& vertex(const std::vector<double>& x);

         virtual Mesh& face(
               Type geom,
               const std::vector<int>& vs,
               Attribute attr = RODIN_DEFAULT_SIMPLEX_ATTRIBUTE);

         virtual Mesh& element(
               Type geom,
               const std::vector<int>& vs,
               Attribute attr = RODIN_DEFAULT_SIMPLEX_ATTRIBUTE);

         virtual Mesh& finalize();

         /**
          * @brief Loads a mesh from file in the given format.
          * @param[in] filename Name of file to read
          * @param[in] fmt Mesh file format
          * @returns Reference to this (for method chaining)
          */
         virtual Mesh& load(
               const boost::filesystem::path& filename,
               IO::FileFormat fmt = IO::FileFormat::MFEM);

         /**
          * @brief Saves a mesh to file in the given format.
          * @param[in] filename Name of file to write
          * @param[in] fmt Mesh file format
          * @returns Reference to this (for method chaining)
          */
         virtual void save(
               const boost::filesystem::path& filename,
               IO::FileFormat fmt = IO::FileFormat::MFEM, size_t precison = 16) const;

         /**
          * @brief Skins the mesh to obtain its boundary mesh
          * @returns SubMesh object to the boundary region of the mesh
          *
          * This function will "skin" the mesh to return the mesh boundary as a
          * new SubMesh object. The new mesh will be embedded in the original
          * space dimension.
          */
         virtual SubMesh<Context::Serial> skin() const;

         /**
          * @brief Trims the elements with the given material reference.
          * @param[in] attr Attribute to trim
          * @returns SubMesh object to the remaining region of the mesh
          *
          * Convenience function to call trim(const std::set<int>&, int) with
          * only one attribute.
          */
         virtual SubMesh<Context::Serial> trim(int attr);

         /**
          * @brief Trims the elements with the given material references.
          * @param[in] attrs Attributes to trim
          * @returns SubMesh object to the remaining region of the mesh
          *
          * This function will trim the current mesh and return a Submesh
          * object containing the elements which were not trimmed from the
          * original mesh.
          */
         virtual SubMesh<Context::Serial> trim(const std::set<int>& attrs);

         virtual SubMesh<Context::Serial> keep(int attr);

         virtual SubMesh<Context::Serial> keep(const std::set<int>& attrs);

         /**
          * @internal
          * @brief Move constructs a Rodin::Mesh from an mfem::Mesh.
          */
         explicit Mesh(mfem::Mesh&& mesh);

#ifdef RODIN_USE_MPI
         /**
          * @brief Parallelizes the mesh by setting the MPI communicator
          * object.
          * @param[in] comm MPI communicator
          *
          * This method parallelizes the serial mesh object and returns a new
          * parallelized mesh. The parallelized mesh is independent of the
          * serial mesh instance.
          */
         Mesh<Context::Parallel> parallelize(boost::mpi::communicator comm);
#endif

         size_t getElementCount() const override;

         size_t getCount(size_t dim) const override;

         BoundaryIterator getBoundary() const override;

         InterfaceIterator getInterface() const override;

         ElementIterator getElement(size_t idx = 0) const override;

         FaceIterator getFace(size_t idx = 0) const override;

         VertexIterator getVertex(size_t idx = 0) const override;

         SimplexIterator getSimplex(size_t dimension, size_t idx) const override;

         bool isInterface(Index faceIdx) const override;

         bool isBoundary(Index faceIdx) const override;

         bool isParallel() const override
         {
            return false;
         }

         mfem::Mesh& getHandle() override;

         const mfem::Mesh& getHandle() const override;

         virtual bool isSubMesh() const override
         {
            return false;
         }

         Attribute getAttribute(size_t dimension, Index index) const override;

      private:
         mfem::Mesh m_mesh;
         std::map<Index, Index> m_f2b;
   };

#ifdef RODIN_USE_MPI
   using ParallelMesh = Mesh<Context::Parallel>;

   template <>
   class Mesh<Context::Parallel> : public MeshBase
   {
      public:
         /**
          * Constructs a parallel mesh object from a serial mesh instance.
          * @param[in] comm MPI communicator
          * @param[in] serialMesh Serial mesh instance
          *
          * This method parallelizes the serial mesh object and constructs a
          * new parallelized mesh. The parallelized mesh is independent of the
          * serial mesh instance.
          */
         Mesh(boost::mpi::communicator comm, Mesh<Context::Serial>& serialMesh)
            :  m_comm(comm),
               m_mesh(mfem::ParMesh(m_comm, serialMesh.getHandle()))
         {}

         /**
          * @brief Deleted copy constructor.
          */
         Mesh(const Mesh&) = delete;

         /**
          * @brief Move constructor.
          */
         Mesh(Mesh&&) = default;

         /**
          * @brief Gets a constant reference to the MPI communicator object.
          * @returns Constant reference to the MPI communicator object.
          */
         const boost::mpi::communicator& getMPIComm() const
         {
            return m_comm;
         }

         virtual bool isSubMesh() const override
         {
            return false;
         }

         bool isParallel() const override
         {
            return true;
         }

         mfem::ParMesh& getHandle() override
         {
            return m_mesh;
         }

         const mfem::ParMesh& getHandle() const override
         {
            return m_mesh;
         }
      private:
         boost::mpi::communicator m_comm;
         mfem::ParMesh m_mesh;
   };
#endif
}

#include "Mesh.hpp"

#endif
