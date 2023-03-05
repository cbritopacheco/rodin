#include "Rodin/Variational/MFEM.h"
#include "Rodin/Variational/QuadratureRule.h"

#include "Mesh.h"
#include "SimplexTransformation.h"

#include "Simplex.h"

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
    mfem::ElementTransformation& trans = getTransformation().getHandle();
    const Variational::QuadratureRule& qr =
      Variational::QuadratureRule::get(getGeometry(), trans.OrderJ());
    Scalar volume = 0.0;
    for (size_t i = 0; i < qr.size(); i++)
      volume += qr.getWeight(i) * trans.Weight();
    return volume;
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
    : m_simplex(simplex), m_trans(trans), m_rc(rc), m_ip(Variational::Internal::vec2ip(m_rc))
  {
    m_trans.get().getHandle().SetIntPoint(&m_ip);
  }

  const Math::Vector& Point::getCoordinates(Coordinates coords) const
  {
    switch (coords)
    {
      case Coordinates::Physical:
      {
        if (!m_pc.has_value())
        {
          Math::Vector pc(getSimplex().getMesh().getSpaceDimension());
          mfem::Vector tmp(pc.data(), pc.size());
          m_trans.get().getHandle().Transform(m_ip, tmp);
          m_pc.emplace(std::move(pc));
        }
        assert(m_pc.has_value());
        return m_pc.value();
      }
      case Coordinates::Reference:
      {
        return m_rc.get();
      }
    }
  }

  const Math::Matrix& Point::getJacobian() const
  {
    if (!m_jacobian.has_value())
    {
      const size_t rdim = getSimplex().getDimension();
      const size_t sdim = getSimplex().getMesh().getSpaceDimension();
      Math::Matrix jacobian(sdim, rdim);
      mfem::DenseMatrix tmp(jacobian.data(), jacobian.rows(), jacobian.cols());
      assert(&m_trans.get().getHandle().GetIntPoint() == &m_ip);
      tmp = m_trans.get().getHandle().Jacobian();
      m_jacobian.emplace(std::move(jacobian));
    }
    assert(m_jacobian.has_value());
    return m_jacobian.value();
  }

  const Math::Matrix& Point::getJacobianInverse() const
  {
    if (!m_inverseJacobian.has_value())
    {
      const size_t rdim = getSimplex().getDimension();
      const size_t sdim = getSimplex().getMesh().getSpaceDimension();
      Math::Matrix inv(rdim, sdim);
      mfem::DenseMatrix tmp(inv.data(), inv.rows(), inv.cols());
      assert(&m_trans.get().getHandle().GetIntPoint() == &m_ip);
      tmp = m_trans.get().getHandle().InverseJacobian();
      m_inverseJacobian.emplace(std::move(inv));
    }
    assert(m_inverseJacobian.has_value());
    return m_inverseJacobian.value();
  }

  Scalar Point::getDistortion() const
  {
    if (!m_distortion.has_value())
    {
      m_distortion.emplace(m_trans.get().getHandle().Weight());
    }
    assert(m_distortion.has_value());
    return m_distortion.value();
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
