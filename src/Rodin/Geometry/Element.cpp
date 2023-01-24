#include "Element.h"
#include "Mesh.h"

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

  mfem::ElementTransformation& Simplex::getTransformation() const
  {
    if (!m_trans)
    {
      const auto& mesh = getMesh();
      const auto index = getIndex();
      const auto dimension = getDimension();
      const auto attribute = getAttribute();
      if (dimension == getMesh().getDimension())
      {
        mfem::IsoparametricTransformation* trans = new mfem::IsoparametricTransformation;
        trans->Attribute = attribute;
        trans->ElementNo = index;
        trans->ElementType = mfem::ElementTransformation::ELEMENT;
        trans->mesh = nullptr;
        trans->Reset();
        const mfem::Mesh& meshHandle = mesh.getHandle();
        const mfem::GridFunction* nodes = meshHandle.GetNodes();
        if (!nodes)
        {
          meshHandle.GetPointMatrix(getIndex(), trans->GetPointMat());
          trans->SetFE(
              meshHandle.GetTransformationFEforElementType(
                meshHandle.GetElementType(getIndex())));
        }
        else
        {
          assert(false);
        }
        m_trans.reset(trans);
      }
      else if (dimension == getMesh().getDimension() - 1)
      {
        mfem::IsoparametricTransformation* trans = new mfem::IsoparametricTransformation;
        trans->Attribute = attribute;
        trans->ElementNo = index;
        trans->ElementType = mfem::ElementTransformation::FACE;
        trans->mesh = nullptr;
        mfem::DenseMatrix& pm = trans->GetPointMat();
        trans->Reset();
        const mfem::Mesh& meshHandle = mesh.getHandle();
        const mfem::GridFunction* nodes = meshHandle.GetNodes();
        const size_t spaceDim = mesh.getSpaceDimension();
        if (!nodes)
        {

          mfem::Array<int> v;
          meshHandle.GetFaceVertices(index, v);
          const int nv = v.Size();
          pm.SetSize(spaceDim, nv);
          for (size_t i = 0; i < spaceDim; i++)
            for (int j = 0; j < nv; j++)
              pm(i, j) = meshHandle.GetVertex(v[j])[i];
          trans->SetFE(
              meshHandle.GetTransformationFEforElementType(
                meshHandle.GetFaceElementType(index)));
        }
        else
        {
          assert(false);
        }
        m_trans.reset(trans);
      }
      else if (dimension == 0)
      {
        assert(false);
      }
      else
      {
        assert(false);
      }
    }
    return *m_trans;
  }

  std::vector<Geometry::Point> Simplex::getIntegrationRule(int order) const
  {
    const mfem::IntegrationRule* ir =
      &mfem::IntRules.Get(getTransformation().GetGeometryType(), order);
    std::vector<Geometry::Point> res;
    res.reserve(ir->GetNPoints());
    for (int i = 0; i < ir->GetNPoints(); i++)
      res.emplace_back(*this, ir->IntPoint(i));
    return res;
  }

  double Simplex::getVolume() const
  {
    mfem::ElementTransformation& trans = getTransformation();
    const mfem::IntegrationRule& ir =
      mfem::IntRules.Get(static_cast<int>(getGeometry()), trans.OrderJ());
    double volume = 0.0;
    for (int j = 0; j < ir.GetNPoints(); j++)
    {
      const mfem::IntegrationPoint &ip = ir.IntPoint(j);
      trans.SetIntPoint(&ip);
      volume += ip.weight * trans.Weight();
    }
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

  // ---- Point -------------------------------------------------------------
  Point::Point(const Simplex& element, const mfem::IntegrationPoint& ip)
    : m_element(element), m_ip(ip)
  {
    m_element.get().getTransformation().SetIntPoint(&m_ip);
    m_element.get().getTransformation().Transform(m_ip, m_physical);
    assert(static_cast<size_t>(m_physical.Size()) == element.getMesh().getSpaceDimension());
  }

  Point::Point(const Simplex& element, mfem::IntegrationPoint&& ip)
    : m_element(element), m_ip(std::move(ip))
  {}

  size_t Point::getDimension(Coordinates coords) const
  {
    switch (coords)
    {
      case Coordinates::Physical:
      {
        return m_element.get().getMesh().getSpaceDimension();
      }
      case Coordinates::Reference:
      {
        return m_element.get().getMesh().getDimension();
      }
      default:
      {
        assert(false);
        return 0;
      }
    }
  }
}
