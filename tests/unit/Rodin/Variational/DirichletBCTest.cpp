#include <cmath>
#include <memory>

#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  namespace
  {
    constexpr Real BoundaryTolerance = 1e-12;
    constexpr Attribute LeftAttribute = 1;
    constexpr Attribute RightAttribute = 2;
    constexpr Attribute BottomAttribute = 3;
    constexpr Attribute TopAttribute = 4;

    Geometry::Mesh<Context::Local> makeUnitSquareMesh(size_t n)
    {
      assert(n >= 2);
      Geometry::Mesh<Context::Local> mesh =
        LocalMesh::UniformGrid(Polytope::Type::Triangle, { n, n });
      mesh.scale(1.0 / (n - 1));
      mesh.getConnectivity().compute(mesh.getDimension() - 1, mesh.getDimension());
      return mesh;
    }

    void labelBoundaryAttributes(
      Geometry::Mesh<Context::Local>& mesh,
      Attribute left = LeftAttribute,
      Attribute right = RightAttribute,
      Attribute bottom = BottomAttribute,
      Attribute top = TopAttribute)
    {
      for (auto it = mesh.getBoundary(); !it.end(); ++it)
      {
        Real x = 0;
        Real y = 0;
        size_t count = 0;
        for (const auto vertex : it->getVertices())
        {
          const auto coords = mesh.getVertexCoordinates(vertex);
          x += coords(0);
          y += coords(1);
          count++;
        }

        x /= count;
        y /= count;

        if (std::abs(x) < BoundaryTolerance)
          mesh.setAttribute(it.key(), left);
        else if (std::abs(x - 1.0) < BoundaryTolerance)
          mesh.setAttribute(it.key(), right);
        else if (std::abs(y) < BoundaryTolerance)
          mesh.setAttribute(it.key(), bottom);
        else if (std::abs(y - 1.0) < BoundaryTolerance)
          mesh.setAttribute(it.key(), top);
      }
    }

    FlatSet<Index> getBoundaryVertices(
      const Geometry::Mesh<Context::Local>& mesh,
      const FlatSet<Attribute>& attrs)
    {
      FlatSet<Index> res;
      for (auto it = mesh.getBoundary(); !it.end(); ++it)
      {
        const auto attr = it->getAttribute();
        if (!attr || attrs.find(*attr) == attrs.end())
          continue;

        for (const auto vertex : it->getVertices())
          res.insert(vertex);
      }
      return res;
    }
  }

  TEST(Rodin_Variational_Real_P1_SanityTest, TriangularUniformGrid2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    const size_t D = mesh.getDimension();

    mesh.getConnectivity().compute(D - 1, D);

    P1 fes(mesh);

    TrialFunction u(fes);

    RealFunction c = 1;

    DirichletBC dbc(u, c);
    dbc.assemble();

    EXPECT_EQ(std::get<IndexMap<Real>>(dbc.getDOFs()).size(), 4);
  }

  TEST(Rodin_Variational_Real_P1_SanityTest, TriangularUniformGrid16)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 16, 16 });
    const size_t D = mesh.getDimension();
    const Attribute attr = 1;

    mesh.getConnectivity().compute(D - 1, D);

    for (auto it = mesh.getBoundary(); !it.end(); ++it)
    {
      mesh.setAttribute(it.key(), attr);
    }

    P1 fes(mesh);

    TrialFunction u(fes);

    RealFunction c = 1;

    DirichletBC dbc(u, c);
    dbc.on(attr);
    dbc.assemble();

    EXPECT_EQ(std::get<IndexMap<Real>>(dbc.getDOFs()).size(), 60);
  }

  TEST(Rodin_Variational_Real_P1_SanityTest, DirichletBCShapeFunctionIdentity)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    const size_t D = mesh.getDimension();
    const Attribute attr = 1;

    mesh.getConnectivity().compute(D - 1, D);

    for (auto it = mesh.getBoundary(); !it.end(); ++it)
    {
      mesh.setAttribute(it.key(), attr);
    }

    P1 fes(mesh);

    TrialFunction u(fes);
    TrialFunction v(fes);

    DirichletBC dbc(u, v);
    dbc.on(attr);
    dbc.assemble();

    using IdentifiedDOFs = DirichletBCBase<Real>::IdentifiedDOFs;
    ASSERT_TRUE(std::holds_alternative<IdentifiedDOFs>(dbc.getDOFs()));

    const auto& dofs = std::get<IdentifiedDOFs>(dbc.getDOFs());
    EXPECT_EQ(dofs.size(), 4);

    for (const auto& [slave, pair] : dofs)
    {
      const auto& masters = pair.first;
      const auto& coeffs = pair.second;

      ASSERT_EQ(masters.size(), 1);
      ASSERT_EQ(coeffs.size(), 1);
      EXPECT_EQ(masters.coeff(0), slave);
      EXPECT_DOUBLE_EQ(coeffs.coeff(0), 1.0);
    }
  }

  TEST(Rodin_Variational_Real_P1_SanityTest, DirichletBCShapeFunctionPointEvaluation)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    const size_t D = mesh.getDimension();
    const Attribute attr = 1;

    mesh.getConnectivity().compute(D - 1, D);

    for (auto it = mesh.getBoundary(); !it.end(); ++it)
      mesh.setAttribute(it.key(), attr);

    P1 fes(mesh);

    TrialFunction u(fes);
    TrialFunction v(fes);

    DirichletBC dbc(u, (RealFunction(1.0) + F::x) * v);
    dbc.on(attr);
    dbc.assemble();

    using IdentifiedDOFs = DirichletBCBase<Real>::IdentifiedDOFs;
    ASSERT_TRUE(std::holds_alternative<IdentifiedDOFs>(dbc.getDOFs()));

    const auto& dofs = std::get<IdentifiedDOFs>(dbc.getDOFs());
    EXPECT_EQ(dofs.size(), 4);

    for (const auto& [slave, pair] : dofs)
    {
      const auto& masters = pair.first;
      const auto& coeffs = pair.second;

      ASSERT_EQ(masters.size(), 1);
      ASSERT_EQ(coeffs.size(), 1);
      EXPECT_EQ(masters.coeff(0), slave);

      const auto coords = mesh.getVertexCoordinates(slave);
      EXPECT_DOUBLE_EQ(coeffs.coeff(0), 1.0 + coords(0));
    }
  }

  TEST(Rodin_Variational_ShapeFunctionSetPoint, H1MatchesBasis)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D - 1, 0);

    H1<2, Real> fes(std::integral_constant<size_t, 2>{}, mesh);
    TestFunction v(fes);

    Math::SpatialVector<Real> rc(2);
    rc[0] = 0.2;
    rc[1] = 0.3;

    auto it = mesh.getPolytope(D, 0);
    const Geometry::Point p(*it, rc);

    v.setIntegrationPoint(IntegrationPoint(p));

    const auto& fe = fes.getFiniteElement(D, 0);
    for (size_t local = 0; local < fe.getCount(); local++)
      EXPECT_NEAR(v.getBasis(local), fe.getBasis(local)(rc), 1e-14);
  }

  TEST(Rodin_Variational_ShapeFunctionSetPoint, H1GradMatchesBasisGradient)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D - 1, 0);

    H1<2, Real> fes(std::integral_constant<size_t, 2>{}, mesh);
    TestFunction v(fes);
    auto gv = Grad(v);

    Math::SpatialVector<Real> rc(2);
    rc[0] = 0.2;
    rc[1] = 0.3;

    auto it = mesh.getPolytope(D, 0);
    const Geometry::Point p(*it, rc);

    gv.setIntegrationPoint(IntegrationPoint(p));

    const auto& fe = fes.getFiniteElement(D, 0);
    const auto JinvT = p.getJacobianInverse().transpose();
    for (size_t local = 0; local < fe.getCount(); local++)
    {
      Math::SpatialVector<Real> ref(D);
      for (size_t d = 0; d < D; d++)
        ref(d) = fe.getBasis(local).template getDerivative<1>(d)(rc);

      const auto expected = JinvT * ref;
      const auto actual = gv.getBasis(local);
      for (size_t d = 0; d < D; d++)
        EXPECT_NEAR(actual(d), expected(d), 1e-14);
    }
  }

  TEST(Rodin_Variational_DirichletBC, ValueAPIAndBoundarySelection)
  {
    Mesh mesh = makeUnitSquareMesh(4);
    labelBoundaryAttributes(mesh);

    P1 fes(mesh);
    TrialFunction u(fes);

    auto g = F::x + 2 * F::y;
    DirichletBC dbc(u, g);
    dbc.on(LeftAttribute, RightAttribute);

    const FlatSet<Attribute> attrs{ LeftAttribute, RightAttribute };
    EXPECT_EQ(dbc.getAttributes(), attrs);
    EXPECT_FALSE(dbc.isComponent());
    EXPECT_FALSE(dbc.getValueUUID().has_value());
    EXPECT_EQ(dbc.getOperand().getUUID(), u.getUUID());

    dbc.assemble();

    ASSERT_TRUE(std::holds_alternative<IndexMap<Real>>(dbc.getDOFs()));
    const auto& dofs = std::get<IndexMap<Real>>(dbc.getDOFs());
    const auto expectedVertices = getBoundaryVertices(mesh, attrs);
    ASSERT_EQ(dofs.size(), expectedVertices.size());

    for (const auto vertex : expectedVertices)
    {
      const auto it = dofs.find(vertex);
      ASSERT_NE(it, dofs.end());

      const auto coords = mesh.getVertexCoordinates(vertex);
      EXPECT_DOUBLE_EQ(it->second, coords(0) + 2 * coords(1));
    }

    std::unique_ptr<DirichletBCBase<Real>> copy(dbc.copy());
    ASSERT_NE(copy, nullptr);
    ASSERT_TRUE(std::holds_alternative<IndexMap<Real>>(copy->getDOFs()));
    EXPECT_EQ(copy->getOperand().getUUID(), u.getUUID());
    EXPECT_FALSE(copy->getValueUUID().has_value());
    EXPECT_EQ(std::get<IndexMap<Real>>(copy->getDOFs()).size(), dofs.size());
  }

  TEST(Rodin_Variational_DirichletBC, ValueFlatSetSelection)
  {
    Mesh mesh = makeUnitSquareMesh(4);
    labelBoundaryAttributes(mesh);

    P1 fes(mesh);
    TrialFunction u(fes);

    const FlatSet<Attribute> attrs{ BottomAttribute, TopAttribute };
    DirichletBC dbc(u, RealFunction(3.0));
    dbc.on(attrs);
    dbc.assemble();

    ASSERT_TRUE(std::holds_alternative<IndexMap<Real>>(dbc.getDOFs()));
    const auto& dofs = std::get<IndexMap<Real>>(dbc.getDOFs());
    const auto expectedVertices = getBoundaryVertices(mesh, attrs);
    ASSERT_EQ(dofs.size(), expectedVertices.size());

    for (const auto vertex : expectedVertices)
    {
      const auto it = dofs.find(vertex);
      ASSERT_NE(it, dofs.end());
      EXPECT_DOUBLE_EQ(it->second, 3.0);
    }
  }

  TEST(Rodin_Variational_DirichletBC, IdentificationAPIAndScaledRows)
  {
    Mesh mesh = makeUnitSquareMesh(4);
    labelBoundaryAttributes(mesh);

    P1 fes(mesh);
    TrialFunction u(fes);
    TrialFunction v(fes);

    const FlatSet<Attribute> attrs{ LeftAttribute, TopAttribute };
    DirichletBC dbc(u, (RealFunction(1.0) + F::x) * v);
    dbc.on(attrs);

    EXPECT_EQ(dbc.getAttributes(), attrs);
    EXPECT_FALSE(dbc.isComponent());
    ASSERT_TRUE(dbc.getValueUUID().has_value());
    EXPECT_EQ(*dbc.getValueUUID(), v.getUUID());
    EXPECT_EQ(dbc.getOperand().getUUID(), u.getUUID());

    dbc.assemble();

    using IdentifiedDOFs = DirichletBCBase<Real>::IdentifiedDOFs;
    ASSERT_TRUE(std::holds_alternative<IdentifiedDOFs>(dbc.getDOFs()));
    const auto& dofs = std::get<IdentifiedDOFs>(dbc.getDOFs());
    const auto expectedVertices = getBoundaryVertices(mesh, attrs);
    ASSERT_EQ(dofs.size(), expectedVertices.size());

    for (const auto vertex : expectedVertices)
    {
      const auto it = dofs.find(vertex);
      ASSERT_NE(it, dofs.end());

      const auto& masters = it->second.first;
      const auto& coeffs = it->second.second;
      ASSERT_EQ(masters.size(), 1);
      ASSERT_EQ(coeffs.size(), 1);
      EXPECT_EQ(masters.coeff(0), vertex);

      const auto coords = mesh.getVertexCoordinates(vertex);
      EXPECT_DOUBLE_EQ(coeffs.coeff(0), 1.0 + coords(0));
    }

    std::unique_ptr<DirichletBCBase<Real>> copy(dbc.copy());
    ASSERT_NE(copy, nullptr);
    ASSERT_TRUE(copy->getValueUUID().has_value());
    EXPECT_EQ(*copy->getValueUUID(), v.getUUID());
    ASSERT_TRUE(std::holds_alternative<IdentifiedDOFs>(copy->getDOFs()));
    EXPECT_EQ(std::get<IdentifiedDOFs>(copy->getDOFs()).size(), dofs.size());
  }

  TEST(Rodin_Variational_DirichletBC, IdentificationComponentSelection)
  {
    Mesh mesh = makeUnitSquareMesh(4);
    labelBoundaryAttributes(mesh);

    P1 scalarFes(mesh);
    P1 vectorFes(mesh, mesh.getSpaceDimension());
    TrialFunction u(scalarFes);
    TrialFunction v(vectorFes);

    const FlatSet<Attribute> attrs{
      LeftAttribute, RightAttribute, BottomAttribute, TopAttribute
    };
    const auto expectedVertices = getBoundaryVertices(mesh, attrs);
    const auto vertexCount = mesh.getVertexCount();

    DirichletBC xbc(u, v.x());
    xbc.on(attrs);
    xbc.assemble();

    using IdentifiedDOFs = DirichletBCBase<Real>::IdentifiedDOFs;
    ASSERT_TRUE(std::holds_alternative<IdentifiedDOFs>(xbc.getDOFs()));
    const auto& xdofs = std::get<IdentifiedDOFs>(xbc.getDOFs());
    ASSERT_EQ(xdofs.size(), expectedVertices.size());
    for (const auto vertex : expectedVertices)
    {
      const auto it = xdofs.find(vertex);
      ASSERT_NE(it, xdofs.end());
      ASSERT_EQ(it->second.first.size(), 1);
      ASSERT_EQ(it->second.second.size(), 1);
      EXPECT_EQ(it->second.first.coeff(0), vertex);
      EXPECT_DOUBLE_EQ(it->second.second.coeff(0), 1.0);
    }

    DirichletBC ybc(u, v.y());
    ybc.on(attrs);
    ybc.assemble();

    ASSERT_TRUE(std::holds_alternative<IdentifiedDOFs>(ybc.getDOFs()));
    const auto& ydofs = std::get<IdentifiedDOFs>(ybc.getDOFs());
    ASSERT_EQ(ydofs.size(), expectedVertices.size());
    for (const auto vertex : expectedVertices)
    {
      const auto it = ydofs.find(vertex);
      ASSERT_NE(it, ydofs.end());
      ASSERT_EQ(it->second.first.size(), 1);
      ASSERT_EQ(it->second.second.size(), 1);
      EXPECT_EQ(it->second.first.coeff(0), vertex + vertexCount);
      EXPECT_DOUBLE_EQ(it->second.second.coeff(0), 1.0);
    }
  }
}
