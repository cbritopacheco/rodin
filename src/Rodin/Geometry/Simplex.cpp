#include "Mesh.h"

#include "Simplex.h"
#include "SimplexTransformation.h"
#include "../Variational/MFEM.h"

namespace Rodin::Geometry
{
  bool operator<(const Simplex& lhs, const Simplex& rhs)
  {
    return lhs.getIndex() < rhs.getIndex();
  }

  // ---- Simplex -----------------------------------------------------------
  Simplex::Simplex(
      size_t dimension,
      Index index,
      const MeshBase& mesh,
      const std::vector<Index>& vertices,
      Attribute attr)
    :  m_dimension(dimension), m_index(index), m_mesh(mesh),
      m_vertices(vertices), m_attr(attr)
  {
    if (m_dimension == mesh.getDimension())
    {
      m_type = static_cast<Geometry::Type>(mesh.getHandle().GetElementGeometry(index));
    }
    else if (m_dimension == mesh.getDimension() - 1)
    {
      m_type = static_cast<Geometry::Type>(mesh.getHandle().GetFaceGeometry(index));
    }
    else if (m_dimension == 0)
    {
      m_type = Geometry::Type::Point;
    }
    else
    {
      assert(false);
    }
  }

  // VertexIterator Simplex::getVertices() const
  // {
  //   assert(false);
  // }

  SimplexIterator Simplex::getAdjacent() const
  {
    assert(false);
    return SimplexIterator(0, getMesh(), EmptyIndexGenerator());
  }

  SimplexIterator Simplex::getIncident() const
  {
    assert(false);
    return SimplexIterator(0, getMesh(), EmptyIndexGenerator());
  }

  const SimplexTransformation& Simplex::getTransformation() const
  {
    return m_mesh.get().getSimplexTransformation(m_dimension, m_index);
  }

  Scalar Simplex::getVolume() const
  {
    assert(false);
    return 0;
    // mfem::ElementTransformation& trans = getTransformation();
    // const mfem::IntegrationRule& ir =
    //   mfem::IntRules.Get(static_cast<int>(getGeometry()), trans.OrderJ());
    // double volume = 0.0;
    // for (int j = 0; j < ir.GetNPoints(); j++)
    // {
    //   const mfem::IntegrationPoint &ip = ir.IntPoint(j);
    //   trans.SetIntPoint(&ip);
    //   volume += ip.weight * trans.Weight();
    // }
    // return volume;
  }

  // VertexIterator Simplex::getVertices() const
  // {
  //   // mfem::Array<int> vs;
  //   // m_data.element->GetVertices(vs);
  //   // return std::vector<int>(vs.begin(), vs.end());
  //   assert(false);
  // }

  // ---- Element -----------------------------------------------------------
  Element::Element(
      Index index,
      const MeshBase& mesh, const std::vector<Index>& vertices, Attribute attr)
    : Simplex(mesh.getDimension(), index, mesh, vertices, attr)
  {}

  // ---- Face --------------------------------------------------------------
  Face::Face(
      Index index,
      const MeshBase& mesh, const std::vector<Index>& vertices, Attribute attr)
    : Simplex(mesh.getDimension() - 1, index, mesh, vertices, attr)
  {}

  bool Face::isBoundary() const
  {
    return getMesh().isBoundary(getIndex());
  }

  bool Face::isInterface() const
  {
    return getMesh().isInterface(getIndex());
  }

  // ---- Vertex -------------------------------------------------------------
  Vertex::Vertex(
      Index index,
      const MeshBase& mesh, const Math::Vector& coordinates, Attribute attr)
    : Simplex(0, index, mesh, {index}, attr), m_coordinates(coordinates)
  {}

  Scalar Vertex::operator()(size_t i) const
  {
    assert(i < static_cast<size_t>(m_coordinates.size()));
    return m_coordinates(i);
  }

  // ---- Point --------------------------------------------------------------
  Point::Point(const Simplex& simplex, const SimplexTransformation& trans, const Math::Vector& rc)
    : m_simplex(simplex), m_trans(trans), m_rc(rc)
  {}

  const Math::Vector& Point::getPhysical() const
  {
    if (!m_pc.has_value())
    {
      m_pc.emplace(m_trans.get().transform(m_rc));
      assert(m_pc->size() == static_cast<int>(m_simplex.get().getMesh().getSpaceDimension()));
    }
    assert(m_pc.has_value());
    return m_pc.value();
  }

  const Math::Matrix& Point::getJacobian() const
  {
    if (!m_jacobian.has_value())
    {
      m_jacobian.emplace(m_trans.get().jacobian(getReference()));
    }
    assert(m_jacobian.has_value());
    return m_jacobian.value();
  }

  const Math::Matrix& Point::getInverseJacobian() const
  {
    if (!m_inverseJacobian.has_value())
    {
      mfem::IntegrationPoint ip = Variational::Internal::vec2ip(m_rc);
      m_trans.get().getHandle().SetIntPoint(&ip);
      mfem::DenseMatrix tmp(m_trans.get().getHandle().InverseJacobian());
      m_inverseJacobian.emplace(Eigen::Map<Math::Matrix>(tmp.Data(), tmp.NumRows(), tmp.NumCols()));
    }
    assert(m_inverseJacobian.has_value());
    return m_inverseJacobian.value();
  }

  size_t Point::getDimension(Coordinates coords) const
  {
    switch (coords)
    {
      case Coordinates::Physical:
      {
        return m_simplex.get().getMesh().getSpaceDimension();
      }
      case Coordinates::Reference:
      {
        return m_simplex.get().getMesh().getDimension();
      }
      default:
      {
        assert(false);
        return 0;
      }
    }
  }
}
