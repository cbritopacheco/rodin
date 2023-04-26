#include "Rodin/Variational/MFEM.h"
#include "Rodin/Variational/QuadratureRule.h"

#include "Mesh.h"
#include "SimplexTransformation.h"

#include "Simplex.h"

namespace Rodin::Geometry
{
  bool operator<(const Polytope& lhs, const Polytope& rhs)
  {
    return lhs.getIndex() < rhs.getIndex();
  }

  // ---- Simplex -----------------------------------------------------------
  Polytope::Polytope(size_t dimension, Index index, const MeshBase& mesh)
    : m_dimension(dimension), m_index(index), m_mesh(mesh)
  {}

  Attribute Polytope::getAttribute() const
  {
    return getMesh().getAttribute(getDimension(), getIndex());
  }

  Polytope::Geometry Polytope::getGeometry() const
  {
    return getMesh().getGeometry(getDimension(), getIndex());
  }

  VertexIterator Polytope::getVertex() const
  {
    const auto& vertices = getVertices();
    return VertexIterator(
        getMesh(), IteratorIndexGenerator(vertices.begin(), vertices.end()));
  }

  const Array<Index>& Polytope::getVertices() const
  {
    return m_mesh.get().getConnectivity().getPolytope(getDimension(), getIndex());
  }

  PolytopeIterator Polytope::getAdjacent() const
  {
    assert(false);
    return PolytopeIterator(0, getMesh(), EmptyIndexGenerator());
  }

  PolytopeIterator Polytope::getIncident() const
  {
    assert(false);
    return PolytopeIterator(0, getMesh(), EmptyIndexGenerator());
  }

  const PolytopeTransformation& Polytope::getTransformation() const
  {
    return m_mesh.get().getPolytopeTransformation(m_dimension, m_index);
  }

  Scalar Polytope::getVolume() const
  {
    mfem::ElementTransformation& trans = getTransformation().getHandle();
    const Variational::QuadratureRule& qr =
      Variational::QuadratureRule::get(*this, trans.OrderJ());
    Scalar volume = 0.0;
    for (size_t i = 0; i < qr.size(); i++)
      volume += qr.getWeight(i) * trans.Weight();
    return volume;
  }

  // ---- Element -----------------------------------------------------------
  Element::Element(Index index, const MeshBase& mesh)
    : Polytope(mesh.getDimension(), index, mesh)
  {}

  // ---- Face --------------------------------------------------------------
  Face::Face(Index index, const MeshBase& mesh)
    : Polytope(mesh.getDimension() - 1, index, mesh)
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
  Vertex::Vertex(Index index, const MeshBase& mesh, const Math::Vector& coordinates)
    : Polytope(0, index, mesh), m_coordinates(coordinates)
  {}

  Scalar Vertex::operator()(size_t i) const
  {
    assert(i < static_cast<size_t>(m_coordinates.size()));
    return m_coordinates(i);
  }

  // ---- Point --------------------------------------------------------------
  Point::Point(const Polytope& simplex, const PolytopeTransformation& trans, const Math::Vector& rc)
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

    return m_pc.value(); // Some compilers complain, so return any value
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
