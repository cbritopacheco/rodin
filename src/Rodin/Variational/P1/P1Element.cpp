#include "P1Element.h"

namespace Rodin::Variational
{
  const Geometry::GeometryIndexed<std::vector<ScalarP1Element::LinearForm>>
  ScalarP1Element::s_ls =
  {
    { Geometry::Polytope::Type::Point,
      {
        { 0, Geometry::Polytope::Type::Point }
      }
    },
    { Geometry::Polytope::Type::Segment,
      {
        { 0, Geometry::Polytope::Type::Segment },
        { 1, Geometry::Polytope::Type::Segment }
      }
    },
    { Geometry::Polytope::Type::Triangle,
      {
        { 0, Geometry::Polytope::Type::Triangle },
        { 1, Geometry::Polytope::Type::Triangle },
        { 2, Geometry::Polytope::Type::Triangle }
      }
    },
    { Geometry::Polytope::Type::Quadrilateral,
      {
        { 0, Geometry::Polytope::Type::Quadrilateral },
        { 1, Geometry::Polytope::Type::Quadrilateral },
        { 2, Geometry::Polytope::Type::Quadrilateral },
        { 3, Geometry::Polytope::Type::Quadrilateral }
      }
    },
    { Geometry::Polytope::Type::Tetrahedron,
      {
        { 0, Geometry::Polytope::Type::Tetrahedron },
        { 1, Geometry::Polytope::Type::Tetrahedron },
        { 2, Geometry::Polytope::Type::Tetrahedron },
        { 3, Geometry::Polytope::Type::Tetrahedron }
      }
    }
  };

  const Geometry::GeometryIndexed<std::vector<ScalarP1Element::BasisFunction>>
  ScalarP1Element::s_basis =
  {
    { Geometry::Polytope::Type::Point,
      {
        { 0, Geometry::Polytope::Type::Point }
      }
    },
    { Geometry::Polytope::Type::Segment,
      {
        { 0, Geometry::Polytope::Type::Segment },
        { 1, Geometry::Polytope::Type::Segment }
      }
    },
    { Geometry::Polytope::Type::Triangle,
      {
        { 0, Geometry::Polytope::Type::Triangle },
        { 1, Geometry::Polytope::Type::Triangle },
        { 2, Geometry::Polytope::Type::Triangle }
      }
    },
    { Geometry::Polytope::Type::Quadrilateral,
      {
        { 0, Geometry::Polytope::Type::Quadrilateral },
        { 1, Geometry::Polytope::Type::Quadrilateral },
        { 2, Geometry::Polytope::Type::Quadrilateral },
        { 3, Geometry::Polytope::Type::Quadrilateral }
      }
    },
    { Geometry::Polytope::Type::Tetrahedron,
      {
        { 0, Geometry::Polytope::Type::Tetrahedron },
        { 1, Geometry::Polytope::Type::Tetrahedron },
        { 2, Geometry::Polytope::Type::Tetrahedron },
        { 3, Geometry::Polytope::Type::Tetrahedron }
      }
    }
  };

  const Geometry::GeometryIndexed<std::vector<ScalarP1Element::GradientFunction>>
  ScalarP1Element::s_gradient =
  {
    { Geometry::Polytope::Type::Point,
      {
        { 0, Geometry::Polytope::Type::Point }
      }
    },
    { Geometry::Polytope::Type::Segment,
      {
        { 0, Geometry::Polytope::Type::Segment },
        { 1, Geometry::Polytope::Type::Segment }
      }
    },
    { Geometry::Polytope::Type::Triangle,
      {
        { 0, Geometry::Polytope::Type::Triangle },
        { 1, Geometry::Polytope::Type::Triangle },
        { 2, Geometry::Polytope::Type::Triangle }
      }
    },
    { Geometry::Polytope::Type::Quadrilateral,
      {
        { 0, Geometry::Polytope::Type::Quadrilateral },
        { 1, Geometry::Polytope::Type::Quadrilateral },
        { 2, Geometry::Polytope::Type::Quadrilateral },
        { 3, Geometry::Polytope::Type::Quadrilateral }
      }
    },
    { Geometry::Polytope::Type::Tetrahedron,
      {
        { 0, Geometry::Polytope::Type::Tetrahedron },
        { 1, Geometry::Polytope::Type::Tetrahedron },
        { 2, Geometry::Polytope::Type::Tetrahedron },
        { 3, Geometry::Polytope::Type::Tetrahedron }
      }
    }
  };

  const Geometry::GeometryIndexed<Math::PointMatrix> ScalarP1Element::s_nodes =
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
              auto sfe = ScalarP1Element(g);
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
  Scalar ScalarP1Element::BasisFunction::operator()(const Math::SpatialVector& r) const
  {
    switch (m_g)
    {
      case Geometry::Polytope::Type::Point:
      {
        return 1;
      }
      case Geometry::Polytope::Type::Segment:
      {
        switch (m_i)
        {
          case 0:
          {
            return 1 - r.x();
          }
          case 1:
          {
            return r.x();
          }
          default:
          {
            assert(false);
            return NAN;
          }
        }
      }
      case Geometry::Polytope::Type::Triangle:
      {
        switch (m_i)
        {
          case 0:
          {
            return -r.x() - r.y() + 1;
          }
          case 1:
          {
            return r.x();
          }
          case 2:
          {
            return r.y();
          }
          default:
          {
            assert(false);
            return NAN;
          }
        }
      }
      case Geometry::Polytope::Type::Quadrilateral:
      {
        switch (m_i)
        {
          case 0:
          {
            auto x = r.x();
            auto y = r.y();
            return x * y - x - y + 1;
          }
          case 1:
          {
            return r.x() * (1 - r.y());
          }
          case 2:
          {
            return r.y() * (1 - r.x());
          }
          case 3:
          {
            return r.x() * r.y();
          }
          default:
          {
            assert(false);
            return NAN;
          }
        }
      }
      case Geometry::Polytope::Type::Tetrahedron:
      {
        switch (m_i)
        {
          case 0:
          {
            return - r.x() - r.y() - r.z() + 1;
          }
          case 1:
          {
            return r.x();
          }
          case 2:
          {
            return r.y();
          }
          case 3:
          {
            return r.z();
          }
          default:
          {
            assert(false);
            return NAN;
          }
        }
      }
    }
    assert(false);
    return NAN;
  }

  void
  ScalarP1Element::GradientFunction::operator()(Math::SpatialVector& out, const Math::SpatialVector& r) const
  {
    switch (m_g)
    {
      case Geometry::Polytope::Type::Point:
      {
        out.resize(1);
        out.coeffRef(0) = 0;
        return;
      }
      case Geometry::Polytope::Type::Segment:
      {
        out.resize(1);
        switch (m_i)
        {
          case 0:
          {
            out.coeffRef(0) = -1;
            return;
          }
          case 1:
          {
            out.coeffRef(0) = 1;
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
        out.resize(2);
        switch (m_i)
        {
          case 0:
          {
            out.setConstant(-1);
            return;
          }
          case 1:
          {
            out.coeffRef(0) = 1;
            out.coeffRef(1) = 0;
            return;
          }
          case 2:
          {
            out.coeffRef(0) = 0;
            out.coeffRef(1) = 1;
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
      case Geometry::Polytope::Type::Quadrilateral:
      {
        out.resize(2);
        switch (m_i)
        {
          case 0:
          {
            out.coeffRef(0) = r.y() - 1;
            out.coeffRef(1) = r.x() - 1;
            return;
          }
          case 1:
          {
            out.coeffRef(0) = 1 - r.y();
            out.coeffRef(1) = -r.y();
            return;
          }
          case 2:
          {
            out.coeffRef(0) = 1 - r.x();
            out.coeffRef(1) = -r.x();
            return;
          }
          case 3:
          {
            out.coeffRef(0) = r.y();
            out.coeffRef(1) = r.x();
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
      case Geometry::Polytope::Type::Tetrahedron:
      {
        out.resize(3);
        switch (m_i)
        {
          case 0:
          {
            out.setConstant(-1);
            return;
          }
          case 1:
          {
            out.coeffRef(0) = 1;
            out.coeffRef(1) = 0;
            out.coeffRef(2) = 0;
            return;
          }
          case 2:
          {
            out.coeffRef(0) = 0;
            out.coeffRef(1) = 1;
            out.coeffRef(2) = 0;
            return;
          }
          case 3:
          {
            out.coeffRef(0) = 0;
            out.coeffRef(1) = 0;
            out.coeffRef(2) = 1;
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
    }
    assert(false);
    out.setConstant(NAN);
  }

  void VectorP1Element::BasisFunction::operator()(Math::Vector& out, const Math::SpatialVector& r) const
  {
    out = Math::Vector::Zero(m_vdim);
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
            out.setConstant(NAN);
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
            out.setConstant(NAN);
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
            out.setConstant(NAN);
            return;
          }
        }
      }
    }
    assert(false);
    out.setConstant(NAN);
  }

  void VectorP1Element::JacobianFunction::operator()(Math::SpatialMatrix& out, const Math::SpatialVector& rc) const
  {
    out = Math::SpatialMatrix::Zero(m_vdim, Geometry::Polytope::getGeometryDimension(m_g));
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
            out.setConstant(NAN);
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
            out.setConstant(NAN);
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
            out.setConstant(NAN);
            return;
          }
        }
      }
    }
    assert(false);
    out.setConstant(NAN);
  }
}
