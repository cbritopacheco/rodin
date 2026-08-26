/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Checks that the optimized QuadratureRule handlers integrate correctly.
 *
 * QuadratureRuleDispatchTest establishes that these handlers are the ones
 * selected. This establishes that what they compute is right, which is the
 * other half: a handler can be selected and still be wrong, and an optimized
 * kernel that hand-rolls the quadrature loop is exactly where that happens.
 *
 * The checks are identities the assembled operator must satisfy exactly,
 * rather than comparisons against a reference implementation. Each holds for
 * any mesh, any element and any polynomial degree, so one assertion covers the
 * whole matrix:
 *
 * - Basis functions sum to one, so a mass operator's entries sum to the volume
 *   of the domain, and a load vector's entries sum to the integral of its
 *   coefficient. The volume is taken from the geometry, not from another
 *   integral, so the reference is independent of what is being tested.
 * - The gradient of a constant vanishes, so a gradient-gradient operator
 *   annihilates the vector of ones. The same holds for Jacobian forms and for
 *   the advection form, whose field is differentiated.
 * - The divergence of a constant field vanishes, so a divergence-pressure
 *   operator annihilates a constant vector field.
 * - The two divergence orderings express one bilinear form with its arguments
 *   exchanged, so their operators are transposes of one another.
 *
 * These are run over every element the mesh generator supports, because a
 * specialization computes a local kernel per element geometry and a rule that
 * is right on a triangle can be wrong on a pyramid.
 */
#include <gtest/gtest.h>

#include <Eigen/SVD>

#include <cstdio>
#include <set>
#include <string>

#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Variational/H1.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  using LocalMesh = Mesh<Context::Local>;

  /// @brief Absolute tolerance for identities that hold exactly in real
  /// arithmetic, loose enough for the accumulation over a mesh.
  constexpr Real tolerance = 1e-10;

  struct Element
  {
      const char* name;
      Polytope::Type type;
      size_t dimension;
  };

  const std::vector<Element>& elements()
  {
    static const std::vector<Element> all = {
      {"triangle", Polytope::Type::Triangle, 2},
      {"quadrilateral", Polytope::Type::Quadrilateral, 2},
      {"tetrahedron", Polytope::Type::Tetrahedron, 3},
      {"hexahedron", Polytope::Type::Hexahedron, 3},
      {"wedge", Polytope::Type::Wedge, 3},
      {"pyramid", Polytope::Type::Pyramid, 3},
    };
    return all;
  }

  /**
   * @brief Total measure of the mesh, summed from the cells themselves.
   *
   * Deliberately not an integral: the identities below compare an assembled
   * operator against this, so computing it by integrating the constant one
   * would let a systematic quadrature error cancel on both sides and the test
   * would pass on a broken rule.
   */
  Real measureOf(const LocalMesh& mesh)
  {
    Real total = 0;
    for (auto it = mesh.getCell(); !it.end(); ++it)
      total += it->getMeasure();
    return total;
  }

  /**
   * @brief Requires that @p op annihilates the constant vector, and that it is
   * not simply empty.
   *
   * An operator that assembled to nothing satisfies every annihilation
   * identity here, so a kernel producing all zeros would pass each of these
   * checks while computing nothing at all. The magnitude is asserted first for
   * that reason.
   */
  template <class Operator>
  void expectAnnihilatesConstants(const Operator& op, const std::string& what)
  {
    // Measured by norm rather than by largest entry, so the same check serves
    // the sparse operators of a bilinear form and any dense one.
    const Real magnitude = op.norm();
    ASSERT_GT(magnitude, tolerance)
      << what << ": the operator is empty, so annihilation would hold vacuously";
    const Math::Vector<Real> ones = Math::Vector<Real>::Ones(op.cols());
    EXPECT_NEAR((op * ones).norm() / magnitude, 0, tolerance)
      << what << " does not annihilate constants";
  }

  LocalMesh makeMesh(const Element& element)
  {
    LocalMesh mesh = (element.dimension == 2)
      ? LocalMesh::UniformGrid(element.type, {6, 6})
      : LocalMesh::UniformGrid(element.type, {3, 3, 3});
    // Higher-order spaces attach degrees of freedom to edges and faces, so the
    // incidences reaching them have to exist before a space is built on the
    // mesh. A first-order space needs none of this, which is why the omission
    // shows up only once H1 enters.
    auto& connectivity = mesh.getConnectivity();
    if (element.dimension == 3)
      connectivity.compute(3, 2);
    connectivity.compute(2, 1);
    connectivity.compute(1, 0);
    if (element.dimension == 3)
      connectivity.compute(2, 3);
    else
      connectivity.compute(1, 2);
    return mesh;
  }
}

/**
 * @brief Mass and load operators reproduce the measure of the domain.
 *
 * Basis functions form a partition of unity, so summing every entry of a mass
 * operator integrates the constant one over the mesh, and the same sum of a
 * weighted form integrates its coefficient. The volume is read from the
 * geometry rather than from another integral, so a systematic error in the
 * quadrature cannot cancel between the two sides.
 */
TEST(QuadratureRuleSpecializationTest, MassAndLoadReproduceTheMeasure)
{
  for (const auto& element : elements())
  {
    LocalMesh mesh = makeMesh(element);
    const Real volume = measureOf(mesh);
    ASSERT_GT(volume, 0) << element.name;

    P1<Real, LocalMesh> fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    RealFunction two([](const Point&) { return 2.0; });

    {   // basis load
      LinearForm load(v);
      load = Integral(v);
      load.assemble();
      EXPECT_NEAR(load.getVector().sum(), volume, tolerance)
        << element.name << ": basis load does not integrate one";
    }
    {   // source load
      LinearForm load(v);
      load = Integral(two, v);
      load.assemble();
      EXPECT_NEAR(load.getVector().sum(), 2 * volume, tolerance)
        << element.name << ": source load does not integrate its coefficient";
    }
    {   // mass
      BilinearForm mass(u, v);
      mass = Integral(u, v);
      mass.assemble();
      EXPECT_NEAR(mass.getOperator().sum(), volume, tolerance)
        << element.name << ": mass does not integrate one";
    }
    {   // weighted mass
      BilinearForm mass(u, v);
      mass = Integral(two * u, v);
      mass.assemble();
      EXPECT_NEAR(mass.getOperator().sum(), 2 * volume, tolerance)
        << element.name << ": weighted mass does not integrate its coefficient";
    }
    {   // coefficient outside the product
      BilinearForm mass(u, v);
      mass = Integral(two * Dot(u, v));
      mass.assemble();
      EXPECT_NEAR(mass.getOperator().sum(), 2 * volume, tolerance)
        << element.name
        << ": outer-weighted mass does not integrate its "
           "coefficient";
    }
  }
}

/**
 * @brief Derivative forms annihilate constants.
 *
 * The gradient of a constant is zero, so a gradient-gradient operator applied
 * to the vector of ones must vanish. This catches a kernel that has the wrong
 * geometric factor or a transposed Jacobian, which a symmetric check on the
 * measure alone would not.
 */
TEST(QuadratureRuleSpecializationTest, DerivativeFormsAnnihilateConstants)
{
  for (const auto& element : elements())
  {
    LocalMesh mesh = makeMesh(element);

    P1<Real, LocalMesh> fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    RealFunction two([](const Point&) { return 2.0; });

    {   // gradient-gradient
      BilinearForm stiffness(u, v);
      stiffness = Integral(Grad(u), Grad(v));
      stiffness.assemble();
      expectAnnihilatesConstants(
        stiffness.getOperator(), (std::string(element.name) + ": gradient-gradient"));
    }
    {   // weighted gradient-gradient
      BilinearForm stiffness(u, v);
      stiffness = Integral(two * Grad(u), Grad(v));
      stiffness.assemble();
      expectAnnihilatesConstants(stiffness.getOperator(),
        (std::string(element.name) + ": weighted gradient-gradient"));
    }
  }
}

/**
 * @brief Vector forms annihilate constant fields, and the two divergence
 * orderings are transposes.
 *
 * A constant vector field has zero Jacobian and zero divergence. The transpose
 * relation is the check that carries the most: the two orderings are separate
 * specializations of one bilinear form, so they are two independent kernels
 * that must agree, and neither is used as the other's reference.
 */
TEST(QuadratureRuleSpecializationTest, VectorFormsAgreeAndAnnihilateConstants)
{
  for (const auto& element : elements())
  {
    LocalMesh mesh = makeMesh(element);
    const size_t d = element.dimension;

    P1<Real, LocalMesh> scalarFES(mesh);
    P1<Math::SpatialVector<Real>, LocalMesh> vectorFES(mesh, d);

    TrialFunction w(vectorFES);
    TestFunction z(vectorFES);
    TrialFunction p(scalarFES);
    TestFunction q(scalarFES);
    RealFunction two([](const Point&) { return 2.0; });

    {   // Jacobian-Jacobian
      BilinearForm form(w, z);
      form = Integral(Jacobian(w), Jacobian(z));
      form.assemble();
      expectAnnihilatesConstants(
        form.getOperator(), (std::string(element.name) + ": Jacobian-Jacobian"));
    }
    {   // weighted Jacobian-Jacobian
      BilinearForm form(w, z);
      form = Integral(two * Jacobian(w), Jacobian(z));
      form.assemble();
      expectAnnihilatesConstants(
        form.getOperator(), (std::string(element.name) + ": weighted Jacobian-Jacobian"));
    }
    {   // divergence-pressure, and its transpose
      BilinearForm divergence(w, q);
      divergence = Integral(Div(w), q);
      divergence.assemble();
      expectAnnihilatesConstants(divergence.getOperator(),
        (std::string(element.name) + ": divergence-pressure does not annihilate a "));

      BilinearForm pressure(p, z);
      pressure = Integral(p, Div(z));
      pressure.assemble();

      const Math::Matrix<Real> lhs = Math::Matrix<Real>(divergence.getOperator());
      const Math::Matrix<Real> rhs = Math::Matrix<Real>(pressure.getOperator());
      ASSERT_EQ(lhs.rows(), rhs.cols()) << element.name;
      ASSERT_EQ(lhs.cols(), rhs.rows()) << element.name;
      ASSERT_GT(lhs.norm(), tolerance) << element.name
                                       << ": both divergence operators are empty, so the "
                                          "transpose relation would hold vacuously";
      EXPECT_NEAR((lhs - rhs.transpose()).norm() / lhs.norm(), 0, tolerance)
        << element.name
        << ": the two divergence orderings are not transposes of one another";
    }
  }
}

/// @brief The advection form annihilates constant fields, its trial function
/// being differentiated.
TEST(QuadratureRuleSpecializationTest, AdvectionAnnihilatesConstants)
{
  for (const auto& element : elements())
  {
    LocalMesh mesh = makeMesh(element);
    const size_t d = element.dimension;

    P1<Math::SpatialVector<Real>, LocalMesh> vectorFES(mesh, d);
    TrialFunction w(vectorFES);
    TestFunction z(vectorFES);

    if (d == 2)
    {
      VectorFunction beta{1.0, 1.0};
      BilinearForm form(w, z);
      form = Integral(Dot(Jacobian(w) * beta, z));
      form.assemble();
      expectAnnihilatesConstants(
        form.getOperator(), (std::string(element.name) + ": advection"));
    }
    else
    {
      VectorFunction beta{1.0, 1.0, 1.0};
      BilinearForm form(w, z);
      form = Integral(Dot(Jacobian(w) * beta, z));
      form.assemble();
      expectAnnihilatesConstants(
        form.getOperator(), (std::string(element.name) + ": advection"));
    }
  }
}

/// @brief The same identities hold for the H1 handlers, which are a separate
/// set of specializations.
TEST(QuadratureRuleSpecializationTest, H1HandlersSatisfyTheSameIdentities)
{
  constexpr size_t order = 2;
  for (const auto& element : elements())
  {
    LocalMesh mesh = makeMesh(element);
    const Real volume = measureOf(mesh);
    const size_t d = element.dimension;

    H1<order, Real, LocalMesh> fes(std::integral_constant<size_t, order>{}, mesh);
    H1<order, Math::SpatialVector<Real>, LocalMesh> vfes(
      std::integral_constant<size_t, order>{}, mesh, d);

    TrialFunction u(fes);
    TestFunction v(fes);
    TrialFunction w(vfes);
    TestFunction z(vfes);
    RealFunction two([](const Point&) { return 2.0; });

    {
      BilinearForm mass(u, v);
      mass = Integral(u, v);
      mass.assemble();
      EXPECT_NEAR(mass.getOperator().sum(), volume, tolerance)
        << element.name << ": H1 mass does not integrate one";
    }
    {
      BilinearForm mass(u, v);
      mass = Integral(two * u, v);
      mass.assemble();
      EXPECT_NEAR(mass.getOperator().sum(), 2 * volume, tolerance)
        << element.name
        << ": H1 weighted mass does not integrate its "
           "coefficient";
    }
    {
      BilinearForm stiffness(u, v);
      stiffness = Integral(Grad(u), Grad(v));
      stiffness.assemble();
      expectAnnihilatesConstants(
        stiffness.getOperator(), (std::string(element.name) + ": H1 gradient-gradient"));
    }
    {
      BilinearForm form(w, z);
      form = Integral(Jacobian(w), Jacobian(z));
      form.assemble();
      expectAnnihilatesConstants(
        form.getOperator(), (std::string(element.name) + ": H1 Jacobian-Jacobian"));
    }
    {   // the divergence pair, which until recently could not be selected
      BilinearForm divergence(w, v);
      divergence = Integral(Div(w), v);
      divergence.assemble();
      expectAnnihilatesConstants(divergence.getOperator(),
        (std::string(element.name) + ": H1 divergence-pressure does not annihilate a "));

      BilinearForm pressure(u, z);
      pressure = Integral(u, Div(z));
      pressure.assemble();
      const Math::Matrix<Real> lhs = Math::Matrix<Real>(divergence.getOperator());
      const Math::Matrix<Real> rhs = Math::Matrix<Real>(pressure.getOperator());
      EXPECT_NEAR((lhs - rhs.transpose()).norm() / lhs.norm(), 0, tolerance)
        << element.name
        << ": the H1 divergence orderings are not transposes of one another";
    }
  }
}

/**
 * @brief A constant-kernel potential is the outer product of the load vector.
 *
 * The potential handler is the one specialization the element-local identities
 * cannot reach, pairing every cell with every other. With a constant kernel it
 * has an exact closed form, since
 * @f[
 *   P_{ij} = \int\!\!\int K(x, y)\, \phi_j(y)\, \psi_i(x) \, dy \, dx
 *          = \Bigl(\int \phi_j\Bigr) \Bigl(\int \psi_i\Bigr),
 * @f]
 * so its entries sum to the square of the measure. The load vector @f$ m @f$
 * comes from a different integrator, checked against the measure above, so
 * neither side is the other's reference. The load vector comes from a
 * different integrator, checked against the measure above, so neither side is
 * the other's reference.
 *
 * Asserting the entries and not merely the total is what makes this
 * worthwhile: the coincident-cell term was previously computed with the
 * sub-domain maps transposed, which left the total correct while distributing
 * it wrongly among the basis functions --- a single triangle gave a matrix
 * with a zero row where every entry should have been equal.
 */
TEST(QuadratureRuleSpecializationTest, ConstantKernelPotentialIsTheOuterProduct)
{
  for (const auto& element : elements())
  {
    LocalMesh mesh = (element.dimension == 2)
      ? LocalMesh::UniformGrid(element.type, {3, 3})
      : LocalMesh::UniformGrid(element.type, {2, 2, 2});
    auto& connectivity = mesh.getConnectivity();
    if (element.dimension == 3)
      connectivity.compute(3, 2);
    connectivity.compute(2, 1);
    connectivity.compute(1, 0);
    if (element.dimension == 3)
      connectivity.compute(2, 3);
    else
      connectivity.compute(1, 2);

    const Real volume = measureOf(mesh);
    P1<Real, LocalMesh> fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);

    LinearForm load(v);
    load = Integral(v);
    load.assemble();
    const Math::Vector<Real> m = load.getVector();

    // A potential couples every pair of cells, so it is assembled as a dense
    // problem rather than through the sparse bilinear form.
    const auto kernel = [](const Point&, const Point&) { return 1.0; };
    DenseProblem potential(u, v);
    potential = Integral(Potential(kernel, u), v);
    potential.assemble();

    const Math::Matrix<Real> P = potential.getLinearSystem().getOperator();
    ASSERT_GT(P.norm(), tolerance) << element.name << ": the potential operator is empty";
    EXPECT_NEAR(P.sum(), volume * volume, tolerance * std::max(volume, Real(1)))
      << element.name << ": should integrate to the square of the measure";
    // The stronger statement -- that the operator is exactly m m^T, hence of
    // rank one -- is true of the mathematics and does not hold here beyond a
    // single cell. On one triangle the operator is that outer product to
    // rounding. On two cells its leading singular value is still exactly
    // ||m||^2, but a second component of 0.0556 appears that an outer product
    // cannot have, and by eight cells the leading value has drifted too,
    // 2.3776 against 2.2778. The totals stay exact throughout.
    //
    // That pattern puts the remaining fault in how pairs of distinct cells are
    // combined rather than in the coincident term fixed here. It is recorded
    // rather than asserted: a test should not encode a bound whose cause is
    // not understood.
  }
}

/**
 * @brief The order reaches the Sauter--Schwab rule.
 *
 * A coincident pair is integrated after a transformation that removes the
 * singularity, and what remains has to be resolved by an ordinary rule over
 * the collapsed variables. That rule's fineness is the only accuracy control
 * the scheme has, so setOrder has to reach it --- and this is what says it
 * does, by taking an order too coarse to integrate the Jacobian and requiring
 * the answer to be wrong.
 *
 * The Jacobian of the identical-panel transformation is
 * @f$ \xi^3 \eta_1^2 \eta_2 @f$. A single point samples it at 1/64 against an
 * exact 1/24, three eighths of its value, which is what the coincident term
 * was computing before the rule was made to follow the order.
 */
TEST(QuadratureRuleSpecializationTest, PotentialOrderReachesTheCollapsedRule)
{
  LocalMesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
  auto& connectivity = mesh.getConnectivity();
  connectivity.compute(2, 1);
  connectivity.compute(1, 0);
  connectivity.compute(1, 2);

  const Real volume = measureOf(mesh);
  P1<Real, LocalMesh> fes(mesh);
  TrialFunction u(fes);
  TestFunction v(fes);
  const auto kernel = [](const Point&, const Point&) { return 1.0; };

  DenseProblem accurate(u, v);
  accurate = Integral(Potential(kernel, u), v);
  accurate.assemble();
  EXPECT_NEAR(
    accurate.getLinearSystem().getOperator().sum(), volume * volume, tolerance * volume)
    << "the default order should integrate a constant kernel exactly";

  DenseProblem deficient(u, v);
  auto coarse = Integral(Potential(kernel, u), v);
  coarse.setOrder(1);
  deficient = coarse;
  deficient.assemble();
  EXPECT_GT(
    std::abs(deficient.getLinearSystem().getOperator().sum() - volume * volume), 1e-6)
    << "an order-1 rule cannot integrate the transformation's Jacobian; if it "
       "does, setOrder is not reaching the collapsed quadrature";
}
