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

  /// @brief Verifies triangular uniform grid 2 for variational real P1 sanity test by checking exact expected values, form assembly.
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

  /// @brief Verifies triangular uniform grid 16 for variational real P1 sanity test by checking exact expected values, form assembly.
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

  /// @brief Verifies dirichlet BC shape function identity for variational real P1 sanity test by checking tolerance-based numerical results, exact expected values, true predicates.
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

  /// @brief Verifies dirichlet BC shape function point evaluation for variational real P1 sanity test by checking tolerance-based numerical results, exact expected values, true predicates.
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

  /// @brief Verifies H1 matches basis for variational point integration point by checking tolerance-based numerical results.
  TEST(Rodin_Variational_PointIntegrationPoint, H1MatchesBasis)
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

    const IntegrationPoint ip(p);
    v.setIntegrationPoint(ip);

    const auto& fe = fes.getFiniteElement(D, 0);
    for (size_t local = 0; local < fe.getCount(); local++)
      EXPECT_NEAR(v.getBasis(local), fe.getBasis(local)(rc), 1e-14);
  }

  /// @brief Verifies H1 grad matches basis gradient for variational point integration point by checking tolerance-based numerical results.
  TEST(Rodin_Variational_PointIntegrationPoint, H1GradMatchesBasisGradient)
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

    const IntegrationPoint ip(p);
    gv.setIntegrationPoint(ip);

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

  /// @brief Verifies value API and boundary selection for variational dirichlet BC by checking tolerance-based numerical results, exact expected values, true predicates.
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

  /// @brief Verifies value flat set selection for variational dirichlet BC by checking tolerance-based numerical results, exact expected values, true predicates.
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

  /// @brief Verifies identification API and scaled rows for variational dirichlet BC by checking tolerance-based numerical results, exact expected values, true predicates.
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

  /// @brief Verifies identification component selection for variational dirichlet BC by checking tolerance-based numerical results, exact expected values, true predicates.
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

  /**
   * @brief Verifies that `DirichletBC(u, -v)` produces coefficient -1.
   *
   * Mathematical content:
   * @f[
   *   C_{sj} = \ell_s^u(-\varphi_j^v) = -\delta_{sj}
   *   \;\Longrightarrow\; \text{coeff} = -1.
   * @f]
   */
  TEST(Rodin_Variational_DirichletBC, IdentificationUnaryMinus)
  {
    Mesh mesh = makeUnitSquareMesh(4);
    labelBoundaryAttributes(mesh);

    P1 fes(mesh);
    TrialFunction u(fes);
    TrialFunction v(fes);

    const FlatSet<Attribute> attrs{ LeftAttribute, RightAttribute };
    DirichletBC dbc(u, -v);
    dbc.on(attrs);
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
      ASSERT_EQ(it->second.first.size(), 1);
      ASSERT_EQ(it->second.second.size(), 1);
      EXPECT_EQ(it->second.first.coeff(0), vertex);
      EXPECT_DOUBLE_EQ(it->second.second.coeff(0), -1.0);
    }
  }

  /**
   * @brief `DirichletBC(u, c*v)` with constant scalar `c` produces
   *        coefficient `c` for the matching DOF.
   */
  TEST(Rodin_Variational_DirichletBC, IdentificationScalarMultiplication)
  {
    Mesh mesh = makeUnitSquareMesh(8);
    labelBoundaryAttributes(mesh);

    P1 fes(mesh);
    TrialFunction u(fes);
    TrialFunction v(fes);

    const Real c = 3.5;
    const FlatSet<Attribute> attrs{ LeftAttribute, RightAttribute };
    DirichletBC dbc(u, RealFunction(c) * v);
    dbc.on(attrs);
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
      ASSERT_EQ(it->second.first.size(), 1);
      ASSERT_EQ(it->second.second.size(), 1);
      EXPECT_EQ(it->second.first.coeff(0), vertex);
      EXPECT_DOUBLE_EQ(it->second.second.coeff(0), c);
    }
  }

  /**
   * @brief Self-identification `DirichletBC(u, u).on(Γ)` produces, for each
   * slave DOF on Γ, the same DOF as master with coefficient 1 — a row
   * @f$ u_s - u_s = 0 @f$ which is algebraically vacuous (singular slave
   * row) but assembled correctly.
   *
   * This test only checks the assembled `IdentifiedDOFs` map; using the BC
   * in a Problem solve would yield a singular system. For periodicity use
   * `PeriodicBC` with an explicit DOF adjacency map.
   */
  TEST(Rodin_Variational_DirichletBC, IdentificationSelfDegenerate)
  {
    Mesh mesh = makeUnitSquareMesh(4);
    labelBoundaryAttributes(mesh);

    P1 fes(mesh);
    TrialFunction u(fes);

    const FlatSet<Attribute> attrs{ LeftAttribute };
    DirichletBC dbc(u, u);
    dbc.on(attrs);

    EXPECT_EQ(dbc.getOperand().getUUID(), u.getUUID());
    ASSERT_TRUE(dbc.getValueUUID().has_value());
    EXPECT_EQ(*dbc.getValueUUID(), u.getUUID());

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
      ASSERT_EQ(it->second.first.size(), 1);
      ASSERT_EQ(it->second.second.size(), 1);
      EXPECT_EQ(it->second.first.coeff(0), vertex);
      EXPECT_DOUBLE_EQ(it->second.second.coeff(0), 1.0);
    }
  }

  /**
   * @brief Vector-to-vector identification: `DirichletBC(u_vec, v_vec)`.
   *
   * Both u and v are vector P1 trial functions on the same vector FES.
   * The slave/master DOFs share the same global indices (same FES) and
   * each master appears exactly once with coefficient 1.
   */
  TEST(Rodin_Variational_DirichletBC, IdentificationVectorVector)
  {
    Mesh mesh = makeUnitSquareMesh(4);
    labelBoundaryAttributes(mesh);

    const size_t d = mesh.getSpaceDimension();
    P1 vfes(mesh, d);
    TrialFunction u(vfes);
    TrialFunction v(vfes);

    const FlatSet<Attribute> attrs{ LeftAttribute };
    DirichletBC dbc(u, v);
    dbc.on(attrs);
    dbc.assemble();

    using IdentifiedDOFs = DirichletBCBase<Real>::IdentifiedDOFs;
    ASSERT_TRUE(std::holds_alternative<IdentifiedDOFs>(dbc.getDOFs()));
    const auto& dofs = std::get<IdentifiedDOFs>(dbc.getDOFs());

    // Each boundary vertex contributes d slave DOFs (one per component).
    const auto expectedVertices = getBoundaryVertices(mesh, attrs);
    ASSERT_EQ(dofs.size(), expectedVertices.size() * d);

    // For each slave we expect exactly one master DOF (the same global
    // index, since slave and master FES coincide) with coefficient 1.
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

  /**
   * @brief When `.on()` is never called, the identification BC defaults to
   * ALL boundary faces.
   */
  TEST(Rodin_Variational_DirichletBC, IdentificationDefaultBoundaryAllFaces)
  {
    Mesh mesh = makeUnitSquareMesh(4);
    labelBoundaryAttributes(mesh);

    P1 fes(mesh);
    TrialFunction u(fes);
    TrialFunction v(fes);

    DirichletBC dbc(u, v);
    EXPECT_TRUE(dbc.getAttributes().empty());

    dbc.assemble();

    using IdentifiedDOFs = DirichletBCBase<Real>::IdentifiedDOFs;
    ASSERT_TRUE(std::holds_alternative<IdentifiedDOFs>(dbc.getDOFs()));
    const auto& dofs = std::get<IdentifiedDOFs>(dbc.getDOFs());

    const FlatSet<Attribute> all{
      LeftAttribute, RightAttribute, BottomAttribute, TopAttribute
    };
    const auto expectedVertices = getBoundaryVertices(mesh, all);
    EXPECT_EQ(dofs.size(), expectedVertices.size());
  }

  /// @brief Verifies identification on tagged interior face for variational dirichlet BC by checking tolerance-based numerical results, exact expected values, true predicates.
  TEST(Rodin_Variational_DirichletBC, IdentificationOnTaggedInteriorFace)
  {
    constexpr Attribute InterfaceAttribute = 9;

    Mesh mesh = makeUnitSquareMesh(4);
    const size_t faceDim = mesh.getDimension() - 1;

    FlatSet<Index> interfaceVertices;
    bool tagged = false;
    for (auto it = mesh.getFace(); !it.end(); ++it)
    {
      const Index face = it->getIndex();
      if (mesh.isBoundary(face))
        continue;

      mesh.setAttribute({ faceDim, face }, InterfaceAttribute);
      for (const auto vertex : it->getVertices())
        interfaceVertices.insert(vertex);
      tagged = true;
      break;
    }
    ASSERT_TRUE(tagged);
    ASSERT_FALSE(interfaceVertices.empty());

    P1 fes(mesh);
    TrialFunction u(fes);
    TrialFunction v(fes);

    DirichletBC dbc(u, RealFunction(2.0) * v);
    dbc.on(InterfaceAttribute);
    dbc.assemble();

    using IdentifiedDOFs = DirichletBCBase<Real>::IdentifiedDOFs;
    ASSERT_TRUE(std::holds_alternative<IdentifiedDOFs>(dbc.getDOFs()));
    const auto& dofs = std::get<IdentifiedDOFs>(dbc.getDOFs());
    ASSERT_EQ(dofs.size(), interfaceVertices.size());

    for (const auto vertex : interfaceVertices)
    {
      const auto it = dofs.find(vertex);
      ASSERT_NE(it, dofs.end());
      ASSERT_EQ(it->second.first.size(), 1);
      ASSERT_EQ(it->second.second.size(), 1);
      EXPECT_EQ(it->second.first.coeff(0), vertex);
      EXPECT_DOUBLE_EQ(it->second.second.coeff(0), 2.0);
    }
  }

  /**
   * @brief getValue() returns a reference whose getUUID() matches v's UUID
   * when the value is a plain shape function.
   */
  TEST(Rodin_Variational_DirichletBC, IdentificationGetValueReturnsShapeFunction)
  {
    Mesh mesh = makeUnitSquareMesh(4);
    labelBoundaryAttributes(mesh);

    P1 fes(mesh);
    TrialFunction u(fes);
    TrialFunction v(fes);

    DirichletBC dbc(u, v);
    dbc.on(LeftAttribute);

    const auto& valueRef = dbc.getValue();
    EXPECT_EQ(valueRef.getUUID(), v.getUUID());
  }

  /**
   * @brief Verify isComponent() returns false for both BC kinds — the
   * `Component<...>` case is handled inside `getBasis`,
   * not via the legacy `isComponent()` flag.
   */
  TEST(Rodin_Variational_DirichletBC, IsComponentAlwaysFalse)
  {
    Mesh mesh = makeUnitSquareMesh(4);
    labelBoundaryAttributes(mesh);

    P1 sfes(mesh);
    P1 vfes(mesh, mesh.getSpaceDimension());
    TrialFunction us(sfes);
    TrialFunction vv(vfes);

    DirichletBC valueBC(us, F::x);
    EXPECT_FALSE(valueBC.isComponent());

    DirichletBC identBC(us, vv.x());
    EXPECT_FALSE(identBC.isComponent());
  }

  /**
   * @brief `DirichletBC(u, A(v) + B(v))`: the assembler should propagate
   * pointwise basis evaluation through `Sum` and produce the sum of the two
   * coefficient rows. For Lagrange same-FES with both summands being v itself,
   * this collapses to coefficient 2 at the matching DOF.
   */
  TEST(Rodin_Variational_DirichletBC, IdentificationSumOfShapeFunctions)
  {
    Mesh mesh = makeUnitSquareMesh(4);
    labelBoundaryAttributes(mesh);

    P1 fes(mesh);
    TrialFunction u(fes);
    TrialFunction v(fes);

    const FlatSet<Attribute> attrs{ LeftAttribute };
    // A(v) = 1*v + 1*v = 2*v at every basis. (A SUM of two ShapeFunction
    // expressions sharing the same underlying TrialFunction.)
    DirichletBC dbc(u, RealFunction(1.0) * v + RealFunction(1.0) * v);
    dbc.on(attrs);
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
      ASSERT_EQ(it->second.first.size(), 1);
      ASSERT_EQ(it->second.second.size(), 1);
      EXPECT_EQ(it->second.first.coeff(0), vertex);
      EXPECT_DOUBLE_EQ(it->second.second.coeff(0), 2.0);
    }
  }

  /**
   * @brief `DirichletBC(u, A(v), d)` keeps the generic shape-function
   * identification row for `A(v)` and stores the known defect through the
   * slave FE's own DOF functional.
   */
  TEST(Rodin_Variational_DirichletBC, AffineIdentificationStoresDefectValues)
  {
    Mesh mesh = makeUnitSquareMesh(4);
    labelBoundaryAttributes(mesh);

    P1 fes(mesh);
    TrialFunction u(fes);
    TrialFunction v(fes);

    constexpr Real gamma = 2.5;
    RealFunction defect = [](const Point& p)
    {
      return 3.0 + 2.0 * p.x();
    };

    const FlatSet<Attribute> attrs{ TopAttribute };
    DirichletBC dbc(u, RealFunction(gamma) * v, defect);
    dbc.on(attrs);
    dbc.assemble();

    using IdentifiedDOFs = DirichletBCBase<Real>::IdentifiedDOFs;
    ASSERT_TRUE(std::holds_alternative<IdentifiedDOFs>(dbc.getDOFs()));
    const auto& dofs = std::get<IdentifiedDOFs>(dbc.getDOFs());
    const auto& values = dbc.getIdentificationValues();
    const auto expectedVertices = getBoundaryVertices(mesh, attrs);
    ASSERT_EQ(dofs.size(), expectedVertices.size());
    ASSERT_EQ(values.size(), expectedVertices.size());

    for (const auto vertex : expectedVertices)
    {
      const auto it = dofs.find(vertex);
      ASSERT_NE(it, dofs.end());
      ASSERT_EQ(it->second.first.size(), 1);
      ASSERT_EQ(it->second.second.size(), 1);
      EXPECT_EQ(it->second.first.coeff(0), vertex);
      EXPECT_DOUBLE_EQ(it->second.second.coeff(0), gamma);

      const auto valueIt = values.find(vertex);
      ASSERT_NE(valueIt, values.end());
      const auto coords = mesh.getVertexCoordinates(vertex);
      EXPECT_NEAR(valueIt->second, 3.0 + 2.0 * coords(0), 1e-14);
    }
  }

  /**
   * @brief Re-assembling after changing the boundary attribute set
   * regenerates the IdentifiedDOFs map for the new region.
   */
  TEST(Rodin_Variational_DirichletBC, IdentificationReassemblesAfterOn)
  {
    Mesh mesh = makeUnitSquareMesh(4);
    labelBoundaryAttributes(mesh);

    P1 fes(mesh);
    TrialFunction u(fes);
    TrialFunction v(fes);

    DirichletBC dbc(u, v);
    using IdentifiedDOFs = DirichletBCBase<Real>::IdentifiedDOFs;

    dbc.on(LeftAttribute);
    dbc.assemble();
    const size_t leftSize = std::get<IdentifiedDOFs>(dbc.getDOFs()).size();
    EXPECT_EQ(leftSize, getBoundaryVertices(mesh, { LeftAttribute }).size());

    dbc.on(LeftAttribute, RightAttribute);
    dbc.assemble();
    const size_t bothSize = std::get<IdentifiedDOFs>(dbc.getDOFs()).size();
    EXPECT_EQ(bothSize,
        getBoundaryVertices(mesh, { LeftAttribute, RightAttribute }).size());
    EXPECT_GT(bothSize, leftSize);
  }

  /// @brief Verifies H1 matches basis for variational shape function set point by checking tolerance-based numerical results.
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

    v.setPoint(p);

    const auto& fe = fes.getFiniteElement(D, 0);
    for (size_t local = 0; local < fe.getCount(); local++)
      EXPECT_NEAR(v.getBasis(local), fe.getBasis(local)(rc), 1e-14);
  }

  /// @brief Verifies H1 grad matches basis gradient for variational shape function set point by checking tolerance-based numerical results.
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

    gv.setPoint(p);

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
}
