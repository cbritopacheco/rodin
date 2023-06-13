#include "P1Element.h"

namespace Rodin::Variational
{
  const Geometry::GeometryIndexed<std::vector<ScalarP1Element::LinearForm>>
  ScalarP1Element::s_ls =
  {
    { Geometry::Polytope::Geometry::Point,
      {
        { 0, Geometry::Polytope::Geometry::Point }
      }
    },
    { Geometry::Polytope::Geometry::Segment,
      {
        { 0, Geometry::Polytope::Geometry::Segment },
        { 1, Geometry::Polytope::Geometry::Segment }
      }
    },
    { Geometry::Polytope::Geometry::Triangle,
      {
        { 0, Geometry::Polytope::Geometry::Triangle },
        { 1, Geometry::Polytope::Geometry::Triangle },
        { 2, Geometry::Polytope::Geometry::Triangle }
      }
    },
    { Geometry::Polytope::Geometry::Quadrilateral,
      {
        { 0, Geometry::Polytope::Geometry::Quadrilateral },
        { 1, Geometry::Polytope::Geometry::Quadrilateral },
        { 2, Geometry::Polytope::Geometry::Quadrilateral },
        { 3, Geometry::Polytope::Geometry::Quadrilateral }
      }
    },
    { Geometry::Polytope::Geometry::Tetrahedron,
      {
        { 0, Geometry::Polytope::Geometry::Tetrahedron },
        { 1, Geometry::Polytope::Geometry::Tetrahedron },
        { 2, Geometry::Polytope::Geometry::Tetrahedron },
        { 3, Geometry::Polytope::Geometry::Tetrahedron }
      }
    }
  };

  const Geometry::GeometryIndexed<std::vector<ScalarP1Element::BasisFunction>>
  ScalarP1Element::s_basis =
  {
    { Geometry::Polytope::Geometry::Point,
      {
        { 0, Geometry::Polytope::Geometry::Point }
      }
    },
    { Geometry::Polytope::Geometry::Segment,
      {
        { 0, Geometry::Polytope::Geometry::Segment },
        { 1, Geometry::Polytope::Geometry::Segment }
      }
    },
    { Geometry::Polytope::Geometry::Triangle,
      {
        { 0, Geometry::Polytope::Geometry::Triangle },
        { 1, Geometry::Polytope::Geometry::Triangle },
        { 2, Geometry::Polytope::Geometry::Triangle }
      }
    },
    { Geometry::Polytope::Geometry::Quadrilateral,
      {
        { 0, Geometry::Polytope::Geometry::Quadrilateral },
        { 1, Geometry::Polytope::Geometry::Quadrilateral },
        { 2, Geometry::Polytope::Geometry::Quadrilateral },
        { 3, Geometry::Polytope::Geometry::Quadrilateral }
      }
    },
    { Geometry::Polytope::Geometry::Tetrahedron,
      {
        { 0, Geometry::Polytope::Geometry::Tetrahedron },
        { 1, Geometry::Polytope::Geometry::Tetrahedron },
        { 2, Geometry::Polytope::Geometry::Tetrahedron },
        { 3, Geometry::Polytope::Geometry::Tetrahedron }
      }
    }
  };

  const Geometry::GeometryIndexed<std::vector<ScalarP1Element::GradientFunction>>
  ScalarP1Element::s_gradient =
  {
    { Geometry::Polytope::Geometry::Point,
      {
        { 0, Geometry::Polytope::Geometry::Point }
      }
    },
    { Geometry::Polytope::Geometry::Segment,
      {
        { 0, Geometry::Polytope::Geometry::Segment },
        { 1, Geometry::Polytope::Geometry::Segment }
      }
    },
    { Geometry::Polytope::Geometry::Triangle,
      {
        { 0, Geometry::Polytope::Geometry::Triangle },
        { 1, Geometry::Polytope::Geometry::Triangle },
        { 2, Geometry::Polytope::Geometry::Triangle }
      }
    },
    { Geometry::Polytope::Geometry::Quadrilateral,
      {
        { 0, Geometry::Polytope::Geometry::Quadrilateral },
        { 1, Geometry::Polytope::Geometry::Quadrilateral },
        { 2, Geometry::Polytope::Geometry::Quadrilateral },
        { 3, Geometry::Polytope::Geometry::Quadrilateral }
      }
    },
    { Geometry::Polytope::Geometry::Tetrahedron,
      {
        { 0, Geometry::Polytope::Geometry::Tetrahedron },
        { 1, Geometry::Polytope::Geometry::Tetrahedron },
        { 2, Geometry::Polytope::Geometry::Tetrahedron },
        { 3, Geometry::Polytope::Geometry::Tetrahedron }
      }
    }
  };

  const Geometry::GeometryIndexed<Math::Matrix> ScalarP1Element::s_nodes =
  {
    { Geometry::Polytope::Geometry::Point,
      Math::Matrix{{0}} },
    { Geometry::Polytope::Geometry::Segment,
      Math::Matrix{{0, 1}} },
    { Geometry::Polytope::Geometry::Triangle,
      Math::Matrix{{0, 1, 0},
                   {0, 0, 1}} },
    { Geometry::Polytope::Geometry::Quadrilateral,
      Math::Matrix{{0, 1, 0, 1},
                   {0, 0, 1, 1}} },
    { Geometry::Polytope::Geometry::Tetrahedron,
      Math::Matrix{{0, 1, 0, 0},
                   {0, 0, 1, 0},
                   {0, 0, 0, 1}} },
  };

  const std::array<Geometry::GeometryIndexed<Math::Matrix>, RODIN_P1_MAX_VECTOR_DIMENSION>
  VectorP1Element::s_nodes = Internal::initVectorP1Nodes();

  const std::array<Geometry::GeometryIndexed<std::vector<VectorP1Element::LinearForm>>, RODIN_P1_MAX_VECTOR_DIMENSION>
  VectorP1Element::s_ls = Internal::initVectorP1LinearForms();

  const std::array<Geometry::GeometryIndexed<std::vector<VectorP1Element::BasisFunction>>, RODIN_P1_MAX_VECTOR_DIMENSION>
  VectorP1Element::s_basis = Internal::initVectorP1Basis();

  const std::array<Geometry::GeometryIndexed<std::vector<VectorP1Element::JacobianFunction>>, RODIN_P1_MAX_VECTOR_DIMENSION>
  VectorP1Element::s_jacobian = Internal::initVectorP1Jacobian();

  namespace Internal
  {
    std::array<Geometry::GeometryIndexed<Math::Matrix>, RODIN_P1_MAX_VECTOR_DIMENSION>
    initVectorP1Nodes()
    {
      std::array<Geometry::GeometryIndexed<Math::Matrix>, RODIN_P1_MAX_VECTOR_DIMENSION> res;
      for (size_t vdim = 1; vdim < RODIN_P1_MAX_VECTOR_DIMENSION; vdim++)
      {
        for (auto g : Geometry::Polytope::Geometries)
        {
          switch (g)
          {
            case Geometry::Polytope::Geometry::Point:
            {
              res[vdim][g] = Math::Matrix::Zero(1, vdim);
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
        for (auto g : Geometry::Polytope::Geometries)
        {
          auto sfe = ScalarP1Element(g);
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
        for (auto g : Geometry::Polytope::Geometries)
        {
          auto sfe = ScalarP1Element(g);
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
        for (auto g : Geometry::Polytope::Geometries)
        {
          auto sfe = ScalarP1Element(g);
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

