/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_PK_PK_HPP
#define RODIN_VARIATIONAL_PK_PK_HPP

#include "Pk.h"

#include "Rodin/Geometry/Connectivity.h"

namespace Rodin::Variational
{
  // ========================================================================
  // Scalar Pk Space Implementation
  // ========================================================================

  template <size_t K, class Scalar>
  Pk<K, Scalar, Geometry::Mesh<Context::Local>>::Pk(const MeshType& mesh)
    : m_mesh(mesh), m_totalDOFs(0)
  {
    build();
  }

  template <size_t K, class Scalar>
  Pk<K, Scalar, Geometry::Mesh<Context::Local>>::Pk(const Pk& other)
    : Parent(other),
      m_mesh(other.m_mesh),
      m_dofs(other.m_dofs),
      m_totalDOFs(other.m_totalDOFs)
  {}

  template <size_t K, class Scalar>
  Pk<K, Scalar, Geometry::Mesh<Context::Local>>::Pk(Pk&& other)
    : Parent(std::move(other)),
      m_mesh(std::move(other.m_mesh)),
      m_dofs(std::move(other.m_dofs)),
      m_totalDOFs(std::move(other.m_totalDOFs))
  {}

  template <size_t K, class Scalar>
  Pk<K, Scalar, Geometry::Mesh<Context::Local>>& 
  Pk<K, Scalar, Geometry::Mesh<Context::Local>>::operator=(Pk&& other)
  {
    if (this != &other)
    {
      Parent::operator=(std::move(other));
      m_mesh = std::move(other.m_mesh);
      m_dofs = std::move(other.m_dofs);
      m_totalDOFs = std::move(other.m_totalDOFs);
    }
    return *this;
  }

  template <size_t K, class Scalar>
  Pk<K, Scalar, Geometry::Mesh<Context::Local>>& 
  Pk<K, Scalar, Geometry::Mesh<Context::Local>>::operator=(const Pk& other)
  {
    if (this != &other)
    {
      Parent::operator=(other);
      m_mesh = other.m_mesh;
      m_dofs = other.m_dofs;
      m_totalDOFs = other.m_totalDOFs;
    }
    return *this;
  }

  template <size_t K, class Scalar>
  const typename Pk<K, Scalar, Geometry::Mesh<Context::Local>>::ElementType& 
  Pk<K, Scalar, Geometry::Mesh<Context::Local>>::getFiniteElement(size_t d, Index i) const
  {
    const auto g = getMesh().getGeometry(d, i);
    switch (g)
    {
      case Geometry::Polytope::Type::Point:
      {
        static thread_local const ElementType s_element(Geometry::Polytope::Type::Point);
        return s_element;
      }
      case Geometry::Polytope::Type::Segment:
      {
        static thread_local const ElementType s_element(Geometry::Polytope::Type::Segment);
        return s_element;
      }
      case Geometry::Polytope::Type::Triangle:
      {
        static thread_local const ElementType s_element(Geometry::Polytope::Type::Triangle);
        return s_element;
      }
      case Geometry::Polytope::Type::Quadrilateral:
      {
        static thread_local const ElementType s_element(Geometry::Polytope::Type::Quadrilateral);
        return s_element;
      }
      case Geometry::Polytope::Type::Tetrahedron:
      {
        static thread_local const ElementType s_element(Geometry::Polytope::Type::Tetrahedron);
        return s_element;
      }
      case Geometry::Polytope::Type::Wedge:
      {
        static thread_local const ElementType s_element(Geometry::Polytope::Type::Wedge);
        return s_element;
      }
    }
    assert(false);
    static thread_local const ElementType s_null;
    return s_null;
  }

  template <size_t K, class Scalar>
  size_t Pk<K, Scalar, Geometry::Mesh<Context::Local>>::getSize() const
  {
    return m_totalDOFs;
  }

  template <size_t K, class Scalar>
  size_t Pk<K, Scalar, Geometry::Mesh<Context::Local>>::getVectorDimension() const
  {
    return 1;
  }

  template <size_t K, class Scalar>
  const typename Pk<K, Scalar, Geometry::Mesh<Context::Local>>::MeshType& 
  Pk<K, Scalar, Geometry::Mesh<Context::Local>>::getMesh() const
  {
    return m_mesh.get();
  }

  template <size_t K, class Scalar>
  const IndexArray& Pk<K, Scalar, Geometry::Mesh<Context::Local>>::getDOFs(size_t d, Index i) const
  {
    return m_dofs[d][i];
  }

  template <size_t K, class Scalar>
  Index Pk<K, Scalar, Geometry::Mesh<Context::Local>>::getGlobalIndex(
      const std::pair<size_t, Index>& idx, Index local) const
  {
    const auto& [d, i] = idx;
    const auto& dofs = getDOFs(d, i);
    assert(local < static_cast<Index>(dofs.size()));
    return dofs(local);
  }

  template <size_t K, class Scalar>
  void Pk<K, Scalar, Geometry::Mesh<Context::Local>>::build()
  {
    const auto& mesh = getMesh();
    const size_t meshDim = mesh.getDimension();

    // Initialize DOF arrays for each dimension
    m_dofs.resize(meshDim + 1);

    if constexpr (K == 0)
    {
      // P0: one DOF per cell
      const size_t numCells = mesh.getCellCount();
      m_dofs[meshDim].reserve(numCells);
      for (size_t i = 0; i < numCells; ++i)
        m_dofs[meshDim].push_back(IndexArray{{static_cast<Index>(i)}});
      m_totalDOFs = numCells;
    }
    else if constexpr (K == 1)
    {
      // P1: one DOF per vertex
      m_totalDOFs = mesh.getVertexCount();

      // Build DOF arrays for each dimension
      for (size_t d = 0; d <= meshDim; ++d)
      {
        const size_t n = mesh.getConnectivity().getCount(d);
        m_dofs[d].reserve(n);
        for (size_t i = 0; i < n; ++i)
        {
          const auto& polytope = mesh.getConnectivity().getPolytope(d, i);
          m_dofs[d].push_back(polytope);
        }
      }
    }
    else
    {
      // For K>=2, we need to build a more complex DOF map
      // DOFs are located at:
      // - vertices (1 DOF per vertex)
      // - interior points of edges (K-1 DOFs per edge for K>=2)
      // - interior points of faces (depends on face geometry and K)
      // - interior points of cells (depends on cell geometry and K)

      Index currentDOF = 0;

      // Step 1: Number vertex DOFs (dimension 0)
      const size_t numVertices = mesh.getConnectivity().getCount(0);
      m_dofs[0].reserve(numVertices);
      for (size_t i = 0; i < numVertices; ++i)
      {
        m_dofs[0].push_back(IndexArray{{currentDOF++}});
      }

      // Step 2: Number edge interior DOFs (dimension 1)
      if (meshDim >= 1)
      {
        const size_t numEdges = mesh.getConnectivity().getCount(1);
        m_dofs[1].reserve(numEdges);
        
        // For edges, we need K-1 interior DOFs for K>=2
        const size_t edgeInteriorDOFs = (K > 1) ? (K - 1) : 0;
        
        for (size_t i = 0; i < numEdges; ++i)
        {
          const auto& edgeVertices = mesh.getConnectivity().getPolytope(1, i);
          IndexArray edgeDOFs(edgeVertices.size() + edgeInteriorDOFs);
          
          // First, add the vertex DOFs
          for (size_t j = 0; j < edgeVertices.size(); ++j)
          {
            edgeDOFs(j) = edgeVertices(j);
          }
          
          // Then add the interior DOFs
          for (size_t j = 0; j < edgeInteriorDOFs; ++j)
          {
            edgeDOFs(edgeVertices.size() + j) = currentDOF++;
          }
          
          m_dofs[1].push_back(edgeDOFs);
        }
      }

      // Step 3: Number face interior DOFs (dimension 2)
      if (meshDim >= 2)
      {
        const size_t numFaces = mesh.getConnectivity().getCount(2);
        m_dofs[2].reserve(numFaces);
        
        for (size_t i = 0; i < numFaces; ++i)
        {
          const auto& faceVertices = mesh.getConnectivity().getPolytope(2, i);
          const auto faceGeometry = mesh.getGeometry(2, i);
          
          // Get a reference element to determine number of interior DOFs
          PkElement<K, Scalar> faceElement(faceGeometry);
          const size_t faceTotalDOFs = faceElement.getCount();
          
          // Face interior DOFs = total DOFs - vertex DOFs - edge interior DOFs
          const size_t numFaceEdges = faceVertices.size(); // For simplicial faces
          const size_t edgeInteriorDOFs = (K > 1) ? (K - 1) : 0;
          const size_t faceInteriorDOFs = faceTotalDOFs - faceVertices.size() - numFaceEdges * edgeInteriorDOFs;
          
          IndexArray faceDOFs(faceTotalDOFs);
          
          // Add vertex DOFs
          for (size_t j = 0; j < faceVertices.size(); ++j)
          {
            faceDOFs(j) = faceVertices(j);
          }
          
          // Add edge interior DOFs (we need to get the edge DOFs from connectivity)
          // This is simplified - assumes ordering matches reference element
          size_t offset = faceVertices.size();
          for (size_t j = 0; j < numFaceEdges * edgeInteriorDOFs; ++j)
          {
            // For now, allocate new DOFs - proper implementation would map to edge DOFs
            faceDOFs(offset + j) = currentDOF++;
          }
          
          // Add face interior DOFs
          offset += numFaceEdges * edgeInteriorDOFs;
          for (size_t j = 0; j < faceInteriorDOFs; ++j)
          {
            faceDOFs(offset + j) = currentDOF++;
          }
          
          m_dofs[2].push_back(faceDOFs);
        }
      }

      // Step 4: Number cell DOFs (dimension meshDim)
      const size_t numCells = mesh.getCellCount();
      m_dofs[meshDim].reserve(numCells);

      for (size_t i = 0; i < numCells; ++i)
      {
        const auto& cellVertices = mesh.getConnectivity().getPolytope(meshDim, i);
        const auto cellGeometry = mesh.getGeometry(meshDim, i);
        
        PkElement<K, Scalar> cellElement(cellGeometry);
        const size_t cellTotalDOFs = cellElement.getCount();
        
        IndexArray cellDOFs(cellTotalDOFs);
        
        // Add vertex DOFs
        size_t localIdx = 0;
        for (size_t j = 0; j < cellVertices.size(); ++j)
        {
          cellDOFs(localIdx++) = cellVertices(j);
        }
        
        // Add edge, face, and cell interior DOFs
        // For simplicity, allocate the remaining DOFs
        // Proper implementation would map to actual edge/face DOFs
        for (size_t j = localIdx; j < cellTotalDOFs; ++j)
        {
          cellDOFs(j) = currentDOF++;
        }
        
        m_dofs[meshDim].push_back(cellDOFs);
      }

      m_totalDOFs = currentDOF;
    }
  }

  // ========================================================================
  // Vector Pk Space Implementation
  // ========================================================================

  template <size_t K, class Scalar>
  Pk<K, Math::Vector<Scalar>, Geometry::Mesh<Context::Local>>::Pk(
      const Geometry::Mesh<ContextType>& mesh, size_t vdim)
    : m_mesh(mesh), m_vdim(vdim), m_totalDOFs(0)
  {
    build();
  }

  template <size_t K, class Scalar>
  Pk<K, Math::Vector<Scalar>, Geometry::Mesh<Context::Local>>::Pk(const Pk& other)
    : Parent(other),
      m_mesh(other.m_mesh),
      m_vdim(other.m_vdim),
      m_dofs(other.m_dofs),
      m_totalDOFs(other.m_totalDOFs)
  {}

  template <size_t K, class Scalar>
  Pk<K, Math::Vector<Scalar>, Geometry::Mesh<Context::Local>>::Pk(Pk&& other)
    : Parent(std::move(other)),
      m_mesh(std::move(other.m_mesh)),
      m_vdim(std::move(other.m_vdim)),
      m_dofs(std::move(other.m_dofs)),
      m_totalDOFs(std::move(other.m_totalDOFs))
  {}

  template <size_t K, class Scalar>
  Pk<K, Math::Vector<Scalar>, Geometry::Mesh<Context::Local>>& 
  Pk<K, Math::Vector<Scalar>, Geometry::Mesh<Context::Local>>::operator=(Pk&& other)
  {
    if (this != &other)
    {
      Parent::operator=(std::move(other));
      m_mesh = std::move(other.m_mesh);
      m_vdim = std::move(other.m_vdim);
      m_dofs = std::move(other.m_dofs);
      m_totalDOFs = std::move(other.m_totalDOFs);
    }
    return *this;
  }

  template <size_t K, class Scalar>
  Pk<K, Math::Vector<Scalar>, Geometry::Mesh<Context::Local>>& 
  Pk<K, Math::Vector<Scalar>, Geometry::Mesh<Context::Local>>::operator=(const Pk& other)
  {
    if (this != &other)
    {
      Parent::operator=(other);
      m_mesh = other.m_mesh;
      m_vdim = other.m_vdim;
      m_dofs = other.m_dofs;
      m_totalDOFs = other.m_totalDOFs;
    }
    return *this;
  }

  template <size_t K, class Scalar>
  const typename Pk<K, Math::Vector<Scalar>, Geometry::Mesh<Context::Local>>::ElementType& 
  Pk<K, Math::Vector<Scalar>, Geometry::Mesh<Context::Local>>::getFiniteElement(size_t d, Index i) const
  {
    const auto& g = getMesh().getGeometry(d, i);
    switch (g)
    {
      case Geometry::Polytope::Type::Point:
      {
        static thread_local std::array<ElementType, RODIN_MAXIMAL_SPACE_DIMENSION + 1> s_elements =
        {
          ElementType(Geometry::Polytope::Type::Point, 0),
          ElementType(Geometry::Polytope::Type::Point, 1),
          ElementType(Geometry::Polytope::Type::Point, 2),
          ElementType(Geometry::Polytope::Type::Point, 3)
        };
        return s_elements[m_vdim];
      }
      case Geometry::Polytope::Type::Segment:
      {
        static thread_local std::array<ElementType, RODIN_MAXIMAL_SPACE_DIMENSION + 1> s_elements =
        {
          ElementType(Geometry::Polytope::Type::Segment, 0),
          ElementType(Geometry::Polytope::Type::Segment, 1),
          ElementType(Geometry::Polytope::Type::Segment, 2),
          ElementType(Geometry::Polytope::Type::Segment, 3)
        };
        return s_elements[m_vdim];
      }
      case Geometry::Polytope::Type::Triangle:
      {
        static thread_local std::array<ElementType, RODIN_MAXIMAL_SPACE_DIMENSION + 1> s_elements =
        {
          ElementType(Geometry::Polytope::Type::Triangle, 0),
          ElementType(Geometry::Polytope::Type::Triangle, 1),
          ElementType(Geometry::Polytope::Type::Triangle, 2),
          ElementType(Geometry::Polytope::Type::Triangle, 3)
        };
        return s_elements[m_vdim];
      }
      case Geometry::Polytope::Type::Quadrilateral:
      {
        static thread_local std::array<ElementType, RODIN_MAXIMAL_SPACE_DIMENSION + 1> s_elements =
        {
          ElementType(Geometry::Polytope::Type::Quadrilateral, 0),
          ElementType(Geometry::Polytope::Type::Quadrilateral, 1),
          ElementType(Geometry::Polytope::Type::Quadrilateral, 2),
          ElementType(Geometry::Polytope::Type::Quadrilateral, 3)
        };
        return s_elements[m_vdim];
      }
      case Geometry::Polytope::Type::Tetrahedron:
      {
        static thread_local std::array<ElementType, RODIN_MAXIMAL_SPACE_DIMENSION + 1> s_elements =
        {
          ElementType(Geometry::Polytope::Type::Tetrahedron, 0),
          ElementType(Geometry::Polytope::Type::Tetrahedron, 1),
          ElementType(Geometry::Polytope::Type::Tetrahedron, 2),
          ElementType(Geometry::Polytope::Type::Tetrahedron, 3)
        };
        return s_elements[m_vdim];
      }
      case Geometry::Polytope::Type::Wedge:
      {
        static thread_local std::array<ElementType, RODIN_MAXIMAL_SPACE_DIMENSION + 1> s_elements =
        {
          ElementType(Geometry::Polytope::Type::Wedge, 0),
          ElementType(Geometry::Polytope::Type::Wedge, 1),
          ElementType(Geometry::Polytope::Type::Wedge, 2),
          ElementType(Geometry::Polytope::Type::Wedge, 3)
        };
        return s_elements[m_vdim];
      }
    }
    assert(false);
    static thread_local ElementType s_null(Geometry::Polytope::Type::Point, 0);
    return s_null;
  }

  template <size_t K, class Scalar>
  size_t Pk<K, Math::Vector<Scalar>, Geometry::Mesh<Context::Local>>::getSize() const
  {
    return m_totalDOFs;
  }

  template <size_t K, class Scalar>
  size_t Pk<K, Math::Vector<Scalar>, Geometry::Mesh<Context::Local>>::getVectorDimension() const
  {
    return m_vdim;
  }

  template <size_t K, class Scalar>
  const typename Pk<K, Math::Vector<Scalar>, Geometry::Mesh<Context::Local>>::MeshType& 
  Pk<K, Math::Vector<Scalar>, Geometry::Mesh<Context::Local>>::getMesh() const
  {
    return m_mesh.get();
  }

  template <size_t K, class Scalar>
  const IndexArray& Pk<K, Math::Vector<Scalar>, Geometry::Mesh<Context::Local>>::getDOFs(
      size_t d, Index i) const
  {
    return m_dofs[d][i];
  }

  template <size_t K, class Scalar>
  Index Pk<K, Math::Vector<Scalar>, Geometry::Mesh<Context::Local>>::getGlobalIndex(
      const std::pair<size_t, Index>& idx, Index local) const
  {
    const auto& [d, i] = idx;
    const auto& dofs = getDOFs(d, i);
    assert(local < static_cast<Index>(dofs.size()));
    return dofs(local);
  }

  template <size_t K, class Scalar>
  void Pk<K, Math::Vector<Scalar>, Geometry::Mesh<Context::Local>>::build()
  {
    const auto& mesh = getMesh();
    const size_t meshDim = mesh.getDimension();

    // Initialize DOF arrays for each dimension
    m_dofs.resize(meshDim + 1);

    if constexpr (K == 0)
    {
      // P0: vdim DOFs per cell
      const size_t numCells = mesh.getCellCount();
      m_dofs[meshDim].reserve(numCells);

      Index currentDOF = 0;
      for (size_t i = 0; i < numCells; ++i)
      {
        IndexArray cellDOFs(m_vdim);
        for (size_t c = 0; c < m_vdim; ++c)
          cellDOFs(c) = currentDOF++;
        m_dofs[meshDim].push_back(cellDOFs);
      }
      m_totalDOFs = currentDOF;

      // Lower dimensional entities have no DOFs
      for (size_t d = 0; d < meshDim; ++d)
      {
        const size_t n = mesh.getConnectivity().getCount(d);
        m_dofs[d].resize(n);
      }
    }
    else if constexpr (K == 1)
    {
      // P1: vdim DOFs per vertex
      const size_t vn = mesh.getVertexCount();
      m_totalDOFs = vn * m_vdim;

      // Build DOF arrays for each dimension
      for (size_t d = 0; d <= meshDim; ++d)
      {
        const size_t n = mesh.getConnectivity().getCount(d);
        m_dofs[d].reserve(n);
        for (size_t i = 0; i < n; ++i)
        {
          const auto& polytope = mesh.getConnectivity().getPolytope(d, i);
          const size_t count = polytope.size();
          IndexArray dofs(count * m_vdim);
          for (size_t local = 0; local < count * m_vdim; ++local)
          {
            const size_t q = local / m_vdim;
            const size_t r = local % m_vdim;
            assert(q < count);
            dofs.coeffRef(local) = polytope(q) + r * vn;
          }
          m_dofs[d].push_back(dofs);
        }
      }
    }
    else
    {
      // For K>=2, build DOF map based on scalar Pk structure and replicate for each vector component
      // We'll use the same global numbering strategy as the scalar version
      
      Index currentDOF = 0;

      // Step 1: Number vertex DOFs (vdim DOFs per vertex)
      const size_t numVertices = mesh.getConnectivity().getCount(0);
      m_dofs[0].reserve(numVertices);
      for (size_t i = 0; i < numVertices; ++i)
      {
        IndexArray vertexDOFs(m_vdim);
        for (size_t c = 0; c < m_vdim; ++c)
        {
          vertexDOFs(c) = currentDOF++;
        }
        m_dofs[0].push_back(vertexDOFs);
      }

      // Step 2: Number edge interior DOFs (vdim * (K-1) DOFs per edge)
      if (meshDim >= 1)
      {
        const size_t numEdges = mesh.getConnectivity().getCount(1);
        m_dofs[1].reserve(numEdges);
        
        const size_t edgeInteriorDOFs = (K > 1) ? (K - 1) : 0;
        
        for (size_t i = 0; i < numEdges; ++i)
        {
          const auto& edgeVertices = mesh.getConnectivity().getPolytope(1, i);
          IndexArray edgeDOFs((edgeVertices.size() + edgeInteriorDOFs) * m_vdim);
          
          // Add vertex DOFs for each component
          for (size_t j = 0; j < edgeVertices.size(); ++j)
          {
            for (size_t c = 0; c < m_vdim; ++c)
            {
              edgeDOFs(j * m_vdim + c) = m_dofs[0][edgeVertices(j)](c);
            }
          }
          
          // Add edge interior DOFs
          size_t offset = edgeVertices.size() * m_vdim;
          for (size_t j = 0; j < edgeInteriorDOFs; ++j)
          {
            for (size_t c = 0; c < m_vdim; ++c)
            {
              edgeDOFs(offset + j * m_vdim + c) = currentDOF++;
            }
          }
          
          m_dofs[1].push_back(edgeDOFs);
        }
      }

      // Step 3: Number face interior DOFs
      if (meshDim >= 2)
      {
        const size_t numFaces = mesh.getConnectivity().getCount(2);
        m_dofs[2].reserve(numFaces);
        
        for (size_t i = 0; i < numFaces; ++i)
        {
          const auto& faceVertices = mesh.getConnectivity().getPolytope(2, i);
          const auto faceGeometry = mesh.getGeometry(2, i);
          
          PkElement<K, Scalar> faceElement(faceGeometry);
          const size_t faceTotalScalarDOFs = faceElement.getCount();
          const size_t faceTotalDOFs = faceTotalScalarDOFs * m_vdim;
          
          const size_t numFaceEdges = faceVertices.size();
          const size_t edgeInteriorDOFs = (K > 1) ? (K - 1) : 0;
          const size_t faceInteriorScalarDOFs = faceTotalScalarDOFs - faceVertices.size() - numFaceEdges * edgeInteriorDOFs;
          
          IndexArray faceDOFs(faceTotalDOFs);
          
          // Add vertex DOFs
          for (size_t j = 0; j < faceVertices.size(); ++j)
          {
            for (size_t c = 0; c < m_vdim; ++c)
            {
              faceDOFs(j * m_vdim + c) = m_dofs[0][faceVertices(j)](c);
            }
          }
          
          // Add edge and face interior DOFs (simplified allocation)
          size_t offset = faceVertices.size() * m_vdim;
          size_t remainingDOFs = faceTotalDOFs - offset;
          for (size_t j = 0; j < remainingDOFs; ++j)
          {
            faceDOFs(offset + j) = currentDOF++;
          }
          
          m_dofs[2].push_back(faceDOFs);
        }
      }

      // Step 4: Number cell DOFs
      const size_t numCells = mesh.getCellCount();
      m_dofs[meshDim].reserve(numCells);

      for (size_t i = 0; i < numCells; ++i)
      {
        const auto& cellVertices = mesh.getConnectivity().getPolytope(meshDim, i);
        const auto cellGeometry = mesh.getGeometry(meshDim, i);
        
        PkElement<K, Scalar> cellElement(cellGeometry);
        const size_t cellTotalScalarDOFs = cellElement.getCount();
        const size_t cellTotalDOFs = cellTotalScalarDOFs * m_vdim;
        
        IndexArray cellDOFs(cellTotalDOFs);
        
        // Add vertex DOFs
        for (size_t j = 0; j < cellVertices.size(); ++j)
        {
          for (size_t c = 0; c < m_vdim; ++c)
          {
            cellDOFs(j * m_vdim + c) = m_dofs[0][cellVertices(j)](c);
          }
        }
        
        // Add remaining DOFs (edges, faces, cell interior)
        size_t offset = cellVertices.size() * m_vdim;
        for (size_t j = offset; j < cellTotalDOFs; ++j)
        {
          cellDOFs(j) = currentDOF++;
        }
        
        m_dofs[meshDim].push_back(cellDOFs);
      }

      m_totalDOFs = currentDOF;
    }
  }
}

#endif
