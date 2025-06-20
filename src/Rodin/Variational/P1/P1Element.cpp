#include "P1Element.h"

namespace Rodin::Variational
{
  const Geometry::GeometryIndexed<Math::PointMatrix> RealP1Element::s_nodes =
  {
    { Geometry::Polytope::Type::Point,
      Math::PointMatrix{{0}} },
    { Geometry::Polytope::Type::Segment,
      Math::PointMatrix{{0, 1}} },
    { Geometry::Polytope::Type::Triangle,
      Math::PointMatrix{{0, 1, 0},
                        {0, 0, 1}} },
    { Geometry::Polytope::Type::Quadrilateral,
      Math::PointMatrix{{0, 1, 0, 1},
                        {0, 0, 1, 1}} },
    { Geometry::Polytope::Type::Tetrahedron,
      Math::PointMatrix{{0, 1, 0, 0},
                        {0, 0, 1, 0},
                        {0, 0, 0, 1}} },
    { Geometry::Polytope::Type::Wedge,
      Math::PointMatrix{{0, 1, 0, 0, 1, 0},
                        {0, 0, 1, 0, 0, 1},
                        {0, 0, 0, 1, 1, 1}} },
  };

  // ComplexP1Element --------------------------------------------------------

  const Geometry::GeometryIndexed<Math::PointMatrix> ComplexP1Element::s_nodes =
  {
    { Geometry::Polytope::Type::Point,
      Math::PointMatrix{{0}} },
    { Geometry::Polytope::Type::Segment,
      Math::PointMatrix{{0, 1}} },
    { Geometry::Polytope::Type::Triangle,
      Math::PointMatrix{{0, 1, 0},
                        {0, 0, 1}} },
    { Geometry::Polytope::Type::Quadrilateral,
      Math::PointMatrix{{0, 1, 0, 1},
                        {0, 0, 1, 1}} },
    { Geometry::Polytope::Type::Tetrahedron,
      Math::PointMatrix{{0, 1, 0, 0},
                        {0, 0, 1, 0},
                        {0, 0, 0, 1}} },
    { Geometry::Polytope::Type::Wedge,
      Math::PointMatrix{{0, 1, 0, 0, 1, 0},
                        {0, 0, 1, 0, 0, 1},
                        {0, 0, 0, 1, 1, 1}} },
  };

  const std::array<Geometry::GeometryIndexed<Math::PointMatrix>, RODIN_P1_MAX_VECTOR_DIMENSION>
  VectorP1Element::s_nodes = Internal::initVectorP1Nodes();

  const std::array<Geometry::GeometryIndexed<std::vector<VectorP1Element::LinearForm>>, RODIN_P1_MAX_VECTOR_DIMENSION>
  VectorP1Element::s_ls = Internal::initVectorP1LinearForms();

  const std::array<Geometry::GeometryIndexed<std::vector<VectorP1Element::BasisFunction>>, RODIN_P1_MAX_VECTOR_DIMENSION>
  VectorP1Element::s_basis = Internal::initVectorP1Basis();

  const std::array<Geometry::GeometryIndexed<std::vector<VectorP1Element::JacobianFunction>>, RODIN_P1_MAX_VECTOR_DIMENSION>
  VectorP1Element::s_jacobian = Internal::initVectorP1Jacobian();

  namespace Internal
  {
    std::array<Geometry::GeometryIndexed<Math::PointMatrix>, RODIN_P1_MAX_VECTOR_DIMENSION>
    initVectorP1Nodes()
    {
      std::array<Geometry::GeometryIndexed<Math::PointMatrix>, RODIN_P1_MAX_VECTOR_DIMENSION> res;
      for (size_t vdim = 1; vdim < RODIN_P1_MAX_VECTOR_DIMENSION; vdim++)
      {
        for (auto g : Geometry::Polytope::Types)
        {
          switch (g)
          {
            case Geometry::Polytope::Type::Point:
            {
              res[vdim][g] = Math::PointMatrix::Zero(1, vdim);
              break;
            }
            default:
            {
              auto sfe = RealP1Element(g);
              const size_t n = Geometry::Polytope::getVertexCount(g);
              const size_t d = Geometry::Polytope::getGeometryDimension(g);
              res[vdim][g].resize(d, n * vdim);
              for (size_t i = 0; i < n * vdim; i++)
                res[vdim][g].col(i) = sfe.getNode(i / vdim);
              break;
            }
          }
        }
      }
      return res;
    }

    std::array<Geometry::GeometryIndexed<std::vector<VectorP1Element::LinearForm>>, RODIN_P1_MAX_VECTOR_DIMENSION>
    initVectorP1LinearForms()
    {
      std::array<Geometry::GeometryIndexed<std::vector<VectorP1Element::LinearForm>>, RODIN_P1_MAX_VECTOR_DIMENSION> res;
      for (size_t vdim = 1; vdim < RODIN_P1_MAX_VECTOR_DIMENSION; vdim++)
      {
        for (auto g : Geometry::Polytope::Types)
        {
          const size_t n = Geometry::Polytope::getVertexCount(g);
          res[vdim][g].resize(n * vdim);
          for (size_t i = 0; i < n * vdim; i++)
            res[vdim][g][i] = VectorP1Element::LinearForm(vdim, i, g);
        }
      }
      return res;
    }

    std::array<Geometry::GeometryIndexed<std::vector<VectorP1Element::BasisFunction>>, RODIN_P1_MAX_VECTOR_DIMENSION>
    initVectorP1Basis()
    {
      std::array<Geometry::GeometryIndexed<std::vector<VectorP1Element::BasisFunction>>, RODIN_P1_MAX_VECTOR_DIMENSION> res;
      for (size_t vdim = 1; vdim < RODIN_P1_MAX_VECTOR_DIMENSION; vdim++)
      {
        for (auto g : Geometry::Polytope::Types)
        {
          const size_t n = Geometry::Polytope::getVertexCount(g);
          res[vdim][g].resize(n * vdim);
          for (size_t i = 0; i < n * vdim; i++)
            res[vdim][g][i] = VectorP1Element::BasisFunction(vdim, i, g);
        }
      }
      return res;
    }

    std::array<Geometry::GeometryIndexed<std::vector<VectorP1Element::JacobianFunction>>, RODIN_P1_MAX_VECTOR_DIMENSION>
    initVectorP1Jacobian()
    {
      std::array<Geometry::GeometryIndexed<std::vector<VectorP1Element::JacobianFunction>>, RODIN_P1_MAX_VECTOR_DIMENSION> res;
      for (size_t vdim = 1; vdim < RODIN_P1_MAX_VECTOR_DIMENSION; vdim++)
      {
        for (auto g : Geometry::Polytope::Types)
        {
          const size_t n = Geometry::Polytope::getVertexCount(g);
          res[vdim][g].resize(n * vdim);
          for (size_t i = 0; i < n * vdim; i++)
            res[vdim][g][i] = VectorP1Element::JacobianFunction(vdim, i, g);
        }
      }
      return res;
    }
  }
}


namespace Rodin::Variational
{
  void VectorP1Element::BasisFunction::operator()(Math::Vector<Real>& out, const Math::SpatialVector<Real>& r) const
  {
    out = Math::Vector<Real>::Zero(m_vdim);
    const size_t i = m_i % m_vdim;
    const size_t k = m_i / m_vdim;
    assert(k < Geometry::Polytope::getVertexCount(m_g));
    switch (m_g)
    {
      case Geometry::Polytope::Type::Point:
      {
        out.coeffRef(i) = 1;
        return;
      }
      case Geometry::Polytope::Type::Segment:
      {
        switch (k)
        {
          case 0:
          {
            out.coeffRef(i) = 1 - r.x();
            return;
          }
          case 1:
          {
            out.coeffRef(i) = r.x();
            return;
          }
          default:
          {
            assert(false);
            out.setConstant(NAN);
            return;
          }
        }
      }
      case Geometry::Polytope::Type::Triangle:
      {
        switch (k)
        {
          case 0:
          {
            out.coeffRef(i) = -r.x() - r.y() + 1;
            return;
          }
          case 1:
          {
            out.coeffRef(i) = r.x();
            return;
          }
          case 2:
          {
            out.coeffRef(i) = r.y();
            return;
          }
          default:
          {
            assert(false);
            out.setConstant(Math::nan<Real>());
            return;
          }
        }
      }
      case Geometry::Polytope::Type::Quadrilateral:
      {
        switch (k)
        {
          case 0:
          {
            out.coeffRef(i) = r.x() * r.y() - r.x() - r.y() + 1;
            return;
          }
          case 1:
          {
            out.coeffRef(i) = r.x() * (1 - r.y());
            return;
          }
          case 2:
          {
            out.coeffRef(i) = r.y() * (1 - r.x());
            return;
          }
          case 3:
          {
            out.coeffRef(i) = r.x() * r.y();
            return;
          }
          default:
          {
            assert(false);
            out.setConstant(Math::nan<Real>());
            return;
          }
        }
      }
      case Geometry::Polytope::Type::Tetrahedron:
      {
        switch (k)
        {
          case 0:
          {
            out.coeffRef(i) = -r.x() - r.y() - r.z() + 1;
            return;
          }
          case 1:
          {
            out.coeffRef(i) = r.x();
            return;
          }
          case 2:
          {
            out.coeffRef(i) = r.y();
            return;
          }
          case 3:
          {
            out.coeffRef(i) = r.z();
            return;
          }
          default:
          {
            assert(false);
            out.setConstant(Math::nan<Real>());
            return;
          }
        }
      }
      case Geometry::Polytope::Type::Wedge:
      {
        switch (k)
        {
          case 0:
          {
            out.coeffRef(i) = r.x() * r.z() - r.x() + r.y() * r.z() - r.y() - r.z() + 1;
            return;
          }
          case 1:
          {
            out.coeffRef(i) = r.x() * (1 - r.z());
            return;
          }
          case 2:
          {
            out.coeffRef(i) = r.y() * (1 - r.z());
            return;
          }
          case 3:
          {
            out.coeffRef(i) = r.z() * (1 - r.x() - r.y());
            return;
          }
          case 4:
          {
            out.coeffRef(i) = r.x() * r.z();
            return;
          }
          case 5:
          {
            out.coeffRef(i) = r.y() * r.z();
            return;
          }
          default:
          {
            assert(false);
            out.setConstant(Math::nan<Real>());
            return;
          }
        }
      }
    }
    assert(false);
    out.setConstant(NAN);
  }

  void VectorP1Element::JacobianFunction::operator()(Math::SpatialMatrix<Real>& out, const Math::SpatialVector<Real>& rc) const
  {
    out = Math::SpatialMatrix<Real>::Zero(m_vdim, Geometry::Polytope::getGeometryDimension(m_g));
    const size_t i = m_i % m_vdim;
    const size_t k = m_i / m_vdim;
    assert(k < Geometry::Polytope::getVertexCount(m_g));
    switch (m_g)
    {
      case Geometry::Polytope::Type::Point:
      {
        return;
      }
      case Geometry::Polytope::Type::Segment:
      {
        switch (k)
        {
          case 0:
          {
            out.row(i) << -1;
            return;
          }
          case 1:
          {
            out.row(i) << 1;
            return;
          }
          default:
          {
            assert(false);
            out.setConstant(Math::nan<Real>());
            return;
          }
        }
      }
      case Geometry::Polytope::Type::Triangle:
      {
        switch (k)
        {
          case 0:
          {
            out.row(i) << -1, -1;
            return;
          }
          case 1:
          {
            out.row(i) << 1, 0;
            return;
          }
          case 2:
          {
            out.row(i) << 0, 1;
            return;
          }
          default:
          {
            assert(false);
            out.setConstant(Math::nan<Real>());
            return;
          }
        }
      }
      case Geometry::Polytope::Type::Quadrilateral:
      {
        switch (k)
        {
          case 0:
          {
            out.row(i) << rc.y() - 1, rc.x() - 1;
            return;
          }
          case 1:
          {
            out.row(i) << 1 - rc.y(), -rc.x();
            return;
          }
          case 2:
          {
            out.row(i) << -rc.y(), 1 - rc.x();
            return;
          }
          case 3:
          {
            out.row(i) << rc.y(), rc.x();
            return;
          }
          default:
          {
            assert(false);
            out.setConstant(Math::nan<Real>());
            return;
          }
        }
      }
      case Geometry::Polytope::Type::Tetrahedron:
      {
        switch (k)
        {
          case 0:
          {
            out.row(i) << -1, -1, -1;
            return;
          }
          case 1:
          {
            out.row(i) << 1, 0, 0;
            return;
          }
          case 2:
          {
            out.row(i) << 0, 1, 0;
            return;
          }
          case 3:
          {
            out.row(i) << 0, 0, 1;
            return;
          }
          default:
          {
            assert(false);
            out.setConstant(Math::nan<Real>());
            return;
          }
        }
      }
      case Geometry::Polytope::Type::Wedge:
      {
        switch (k)
        {
          case 0:
          {
            out.row(i) << rc.z() - 1, rc.z() - 1, rc.x() + rc.y() - 1;
            return;
          }
          case 1:
          {
            out.row(i) << 1 - rc.z(), 0, -rc.x();
            return;
          }
          case 2:
          {
            out.row(i) << 0, 1 - rc.z(), -rc.y();
            return;
          }
          case 3:
          {
            out.row(i) << -rc.z(), -rc.z(), 1 - rc.x() - rc.y();
            return;
          }
          case 4:
          {
            out.row(i) << rc.z(), 0, rc.x();
            return;
          }
          case 5:
          {
            out.row(i) << 0, rc.z(), rc.y();
            return;
          }
          default:
          {
            assert(false);
            out.setConstant(Math::nan<Real>());
            return;
          }
        }
      }
    }
    assert(false);
    out.setConstant(Math::nan<Real>());
  }
}
