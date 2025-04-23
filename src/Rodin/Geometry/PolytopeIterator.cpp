/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Mesh.h"
#include "PolytopeIterator.h"

namespace Rodin::Geometry
{
  // ---- PolytopeIterator -------------------------------------------------------
  PolytopeIterator::PolytopeIterator(
      size_t dimension, const MeshBase& mesh, IndexGeneratorBase&& gen)
    : Parent(dimension, mesh, std::move(gen))
  {}

  PolytopeIterator& PolytopeIterator::operator=(CellIterator it)
  {
    PolytopeIterator tmp(it.getDimension(), it.getMesh(), std::move(it.getIndexGenerator()));
    *this = std::move(tmp);
    return *this;
  }

  PolytopeIterator& PolytopeIterator::operator=(FaceIterator it)
  {
    PolytopeIterator tmp(it.getDimension(), it.getMesh(), std::move(it.getIndexGenerator()));
    *this = std::move(tmp);
    return *this;
  }

  PolytopeIterator& PolytopeIterator::operator=(VertexIterator it)
  {
    PolytopeIterator tmp(it.getDimension(), it.getMesh(), std::move(it.getIndexGenerator()));
    *this = std::move(tmp);
    return *this;
  }

  Polytope* PolytopeIterator::generate() const
  {
    if (this->end()) return nullptr;
    const auto& gen = getIndexGenerator();
    const auto& index = *gen;
    return new Polytope(getDimension(), index, getMesh());
  }

  // ---- CellIterator -------------------------------------------------------
  CellIterator::CellIterator(const MeshBase& mesh, IndexGeneratorBase&& gen)
    : Parent(mesh.getDimension(), mesh, std::move(gen))
  {}

  Cell* CellIterator::generate() const
  {
    if (this->end()) return nullptr;
    const auto& gen = getIndexGenerator();
    const auto& index = *gen;
    return new Cell(index, getMesh());
  }

  // ---- FaceIterator -------------------------------------------------------
  FaceIterator::FaceIterator(const MeshBase& mesh, IndexGeneratorBase&& gen)
    : Parent(mesh.getDimension() - 1, mesh, std::move(gen))
  {}

  Face* FaceIterator::generate() const
  {
    if (this->end()) return nullptr;
    const auto& gen = getIndexGenerator();
    const auto& index = *gen;
    return new Face(index, getMesh());
  }

  // ---- VertexIterator -----------------------------------------------------
  VertexIterator::VertexIterator(const MeshBase& mesh, IndexGeneratorBase&& gen)
    : Parent(0, mesh, std::move(gen))
  {}

  Vertex* VertexIterator::generate() const
  {
    if (this->end()) return nullptr;
    const auto& gen = getIndexGenerator();
    const auto& index = *gen;
    return new Vertex(index, getMesh());
  }
}
