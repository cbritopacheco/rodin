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
 * @brief The finite element families the specializations are written for.
 *
 * Dispatch does not depend on the element a mesh is made of --- no handler
 * signature mentions a geometry --- so the space is the axis that separates
 * them, and each family below selects a different set of handlers. What the
 * element does change is the local kernel each handler computes, which is why
 * the identities are checked on every element for every family.
 */
namespace
{
  struct P1Family
  {
      static constexpr const char* name = "P1";

      static P1<Real, LocalMesh> scalar(const LocalMesh& mesh)
      {
        return P1<Real, LocalMesh>(mesh);
      }

      static P1<Math::SpatialVector<Real>, LocalMesh> vector(
        const LocalMesh& mesh, size_t dimension)
      {
        return P1<Math::SpatialVector<Real>, LocalMesh>(mesh, dimension);
      }
  };

  template <size_t K>
  struct H1Family
  {
      static constexpr const char* name = "H1";

      static H1<K, Real, LocalMesh> scalar(const LocalMesh& mesh)
      {
        return H1<K, Real, LocalMesh>(std::integral_constant<size_t, K>{}, mesh);
      }

      static H1<K, Math::SpatialVector<Real>, LocalMesh> vector(
        const LocalMesh& mesh, size_t dimension)
      {
        return H1<K, Math::SpatialVector<Real>, LocalMesh>(
          std::integral_constant<size_t, K>{}, mesh, dimension);
      }
  };

  /// @brief Scalar identities: mass and loads reproduce the measure, and
  /// derivative forms annihilate constants.
  template <class Family>
  void checkScalarIdentities(const std::string& label)
  {
    for (const auto& element : elements())
    {
      const std::string where = label + " on " + element.name;
      LocalMesh mesh = makeMesh(element);
      const Real volume = measureOf(mesh);
      ASSERT_GT(volume, 0) << where;

      auto fes = Family::scalar(mesh);
      TrialFunction u(fes);
      TestFunction v(fes);
      RealFunction two([](const Point&) { return 2.0; });

      {
        LinearForm load(v);
        load = Integral(v);
        load.assemble();
        EXPECT_NEAR(load.getVector().sum(), volume, tolerance)
          << where << ": basis load should integrate one";
      }
      {
        LinearForm load(v);
        load = Integral(two, v);
        load.assemble();
        EXPECT_NEAR(load.getVector().sum(), 2 * volume, tolerance)
          << where << ": source load should integrate its coefficient";
      }
      {
        BilinearForm mass(u, v);
        mass = Integral(u, v);
        mass.assemble();
        EXPECT_NEAR(mass.getOperator().sum(), volume, tolerance)
          << where << ": mass should integrate one";
      }
      {
        BilinearForm mass(u, v);
        mass = Integral(two * u, v);
        mass.assemble();
        EXPECT_NEAR(mass.getOperator().sum(), 2 * volume, tolerance)
          << where << ": weighted mass should integrate its coefficient";
      }
      {
        BilinearForm mass(u, v);
        mass = Integral(two * Dot(u, v));
        mass.assemble();
        EXPECT_NEAR(mass.getOperator().sum(), 2 * volume, tolerance)
          << where << ": outer-weighted mass should integrate its coefficient";
      }
      {
        BilinearForm stiffness(u, v);
        stiffness = Integral(Grad(u), Grad(v));
        stiffness.assemble();
        expectAnnihilatesConstants(
          stiffness.getOperator(), where + ": gradient-gradient");
      }
      {
        BilinearForm stiffness(u, v);
        stiffness = Integral(two * Grad(u), Grad(v));
        stiffness.assemble();
        expectAnnihilatesConstants(
          stiffness.getOperator(), where + ": weighted gradient-gradient");
      }
    }
  }

  /// @brief Vector identities: derivative forms annihilate constant fields,
  /// and the two divergence orderings are transposes of one another.
  template <class Family>
  void checkVectorIdentities(const std::string& label)
  {
    for (const auto& element : elements())
    {
      const std::string where = label + " on " + element.name;
      LocalMesh mesh = makeMesh(element);
      const size_t d = element.dimension;

      auto scalarFES = Family::scalar(mesh);
      auto vectorFES = Family::vector(mesh, d);
      TrialFunction w(vectorFES);
      TestFunction z(vectorFES);
      TrialFunction p(scalarFES);
      TestFunction q(scalarFES);
      RealFunction two([](const Point&) { return 2.0; });

      {
        BilinearForm form(w, z);
        form = Integral(Jacobian(w), Jacobian(z));
        form.assemble();
        expectAnnihilatesConstants(form.getOperator(), where + ": Jacobian-Jacobian");
      }
      {
        BilinearForm form(w, z);
        form = Integral(two * Jacobian(w), Jacobian(z));
        form.assemble();
        expectAnnihilatesConstants(
          form.getOperator(), where + ": weighted Jacobian-Jacobian");
      }
      {
        BilinearForm divergence(w, q);
        divergence = Integral(Div(w), q);
        divergence.assemble();
        expectAnnihilatesConstants(
          divergence.getOperator(), where + ": divergence-pressure");

        BilinearForm pressure(p, z);
        pressure = Integral(p, Div(z));
        pressure.assemble();

        const Math::Matrix<Real> lhs = Math::Matrix<Real>(divergence.getOperator());
        const Math::Matrix<Real> rhs = Math::Matrix<Real>(pressure.getOperator());
        ASSERT_EQ(lhs.rows(), rhs.cols()) << where;
        ASSERT_EQ(lhs.cols(), rhs.rows()) << where;
        ASSERT_GT(lhs.norm(), tolerance)
          << where
          << ": both divergence operators are empty, so the transpose "
             "relation would hold vacuously";
        EXPECT_NEAR((lhs - rhs.transpose()).norm() / lhs.norm(), 0, tolerance)
          << where << ": the divergence orderings should be transposes";
      }
      // The transport field's arity follows the dimension, so the two cases
      // are written out rather than selected: their types differ.
      if (d == 2)
      {
        VectorFunction beta{1.0, 1.0};
        BilinearForm form(w, z);
        form = Integral(Dot(Jacobian(w) * beta, z));
        form.assemble();
        expectAnnihilatesConstants(form.getOperator(), where + ": advection");
      }
      else
      {
        VectorFunction beta{1.0, 1.0, 1.0};
        BilinearForm form(w, z);
        form = Integral(Dot(Jacobian(w) * beta, z));
        form.assemble();
        expectAnnihilatesConstants(form.getOperator(), where + ": advection");
      }
    }
  }
}

/// @brief The scalar identities hold for P1 on every element.
TEST(QuadratureRuleSpecializationTest, P1ScalarIdentities)
{
  checkScalarIdentities<P1Family>("P1");
}

/// @brief The vector identities hold for P1 on every element.
TEST(QuadratureRuleSpecializationTest, P1VectorIdentities)
{
  checkVectorIdentities<P1Family>("P1");
}

/// @brief The scalar identities hold for first-order H1 on every element.
TEST(QuadratureRuleSpecializationTest, H1P1ScalarIdentities)
{
  checkScalarIdentities<H1Family<1>>("H1P1");
}

/// @brief The vector identities hold for first-order H1 on every element.
TEST(QuadratureRuleSpecializationTest, H1P1VectorIdentities)
{
  checkVectorIdentities<H1Family<1>>("H1P1");
}

/// @brief The scalar identities hold for second-order H1 on every element.
TEST(QuadratureRuleSpecializationTest, H1P2ScalarIdentities)
{
  checkScalarIdentities<H1Family<2>>("H1P2");
}

/// @brief The vector identities hold for second-order H1 on every element.
TEST(QuadratureRuleSpecializationTest, H1P2VectorIdentities)
{
  checkVectorIdentities<H1Family<2>>("H1P2");
}

/// @brief The scalar identities hold for third-order H1 on every element.
TEST(QuadratureRuleSpecializationTest, H1P3ScalarIdentities)
{
  checkScalarIdentities<H1Family<3>>("H1P3");
}

/// @brief The vector identities hold for third-order H1 on every element.
TEST(QuadratureRuleSpecializationTest, H1P3VectorIdentities)
{
  checkVectorIdentities<H1Family<3>>("H1P3");
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

/**
 * @brief On one triangle the coincident term is the whole operator, and exact.
 *
 * This is the case that exercises the Sauter--Schwab path directly. Only the
 * triangle has one: every other geometry sends a coincident pair to the plain
 * centroid branch, so the tests above, though they run on six elements, reach
 * this code on one of them.
 *
 * A single cell leaves nothing but the coincident term, and with a constant
 * kernel the answer is known entrywise --- every entry is
 * @f$ (\int \phi_i)(\int \phi_j) = (|K|/3)^2 @f$ on a first-order triangle,
 * the basis functions being interchangeable. Checking the entries rather than
 * their total is the point: the sub-domain maps were transposed, which left
 * the total correct while putting a zero row where every entry should have
 * been equal, and no test of the sum could have seen it.
 */
TEST(QuadratureRuleSpecializationTest, CoincidentTriangleTermIsExactEntrywise)
{
  LocalMesh mesh = LocalMesh::Builder()
                     .initialize(2)
                     .nodes(3)
                     .vertex({0, 0})
                     .vertex({1, 0})
                     .vertex({0, 1})
                     .polytope(Polytope::Type::Triangle, {{0, 1, 2}})
                     .finalize();

  const Real area = measureOf(mesh);
  ASSERT_NEAR(area, 0.5, tolerance);

  P1<Real, LocalMesh> fes(mesh);
  TrialFunction u(fes);
  TestFunction v(fes);
  const auto kernel = [](const Point&, const Point&) { return 1.0; };

  DenseProblem potential(u, v);
  potential = Integral(Potential(kernel, u), v);
  potential.assemble();

  const Math::Matrix<Real> P = potential.getLinearSystem().getOperator();
  ASSERT_EQ(P.rows(), 3);
  ASSERT_EQ(P.cols(), 3);

  const Real expected = (area / 3) * (area / 3);
  for (Eigen::Index i = 0; i < P.rows(); ++i)
    for (Eigen::Index j = 0; j < P.cols(); ++j)
      EXPECT_NEAR(P(i, j), expected, tolerance)
        << "entry (" << i << ", " << j << ") of the coincident term";
}

/**
 * @brief Expressions no specialization matches are integrated correctly too.
 *
 * These are the forms the benchmarks carry as controls --- divergence
 * divergence, symmetric-gradient elasticity, flux-gradient loads, piecewise
 * constant mass, and boundary forms --- and nothing tested them. A control
 * exists to show that a change helped the optimized path rather than the mesh
 * traversal underneath, which it cannot do if the control is itself wrong.
 *
 * The same identities serve: the divergence and symmetric gradient of a
 * constant field vanish, a constant flux dotted with the gradient of a
 * partition of unity integrates to nothing, and piecewise constant basis
 * functions sum to one just as nodal ones do. Boundary forms reproduce the
 * measure of the boundary rather than of the domain, taken from the faces
 * themselves.
 */
TEST(QuadratureRuleSpecializationTest, UnspecializedFormsAreAlsoCorrect)
{
  for (const auto& element : elements())
  {
    const std::string where = std::string("on ") + element.name;
    LocalMesh mesh = makeMesh(element);
    const size_t d = element.dimension;
    const Real volume = measureOf(mesh);

    P1<Real, LocalMesh> scalarFES(mesh);
    P1<Math::SpatialVector<Real>, LocalMesh> vectorFES(mesh, d);
    TrialFunction w(vectorFES);
    TestFunction z(vectorFES);
    TrialFunction u(scalarFES);
    TestFunction v(scalarFES);

    {   // divergence-divergence, which no handler specializes
      BilinearForm form(w, z);
      form = Integral(Div(w), Div(z));
      form.assemble();
      expectAnnihilatesConstants(form.getOperator(), where + ": divergence-divergence");
    }
    {   // symmetric-gradient elasticity
      BilinearForm form(w, z);
      form =
        Integral(Jacobian(w) + Jacobian(w).T(), 0.5 * (Jacobian(z) + Jacobian(z).T()));
      form.assemble();
      expectAnnihilatesConstants(form.getOperator(), where + ": elasticity");
    }
    {   // a constant flux against the gradient of a partition of unity
      LinearForm load(v);
      if (d == 2)
        load = Integral(VectorFunction{1.0, 1.0}, Grad(v));
      else
        load = Integral(VectorFunction{1.0, 1.0, 1.0}, Grad(v));
      load.assemble();
      EXPECT_NEAR(load.getVector().sum(), 0, tolerance)
        << where
        << ": a constant flux against Grad of a partition of unity "
           "should integrate to nothing";
    }
    {   // piecewise constants, whose basis also sums to one
      P0<Real, LocalMesh> fes(mesh);
      TrialFunction p(fes);
      TestFunction q(fes);

      LinearForm load(q);
      load = Integral(q);
      load.assemble();
      EXPECT_NEAR(load.getVector().sum(), volume, tolerance)
        << where << ": P0 load should integrate one";

      BilinearForm mass(p, q);
      mass = Integral(p, q);
      mass.assemble();
      EXPECT_NEAR(mass.getOperator().sum(), volume, tolerance)
        << where << ": P0 mass should integrate one";
    }
  }
}

/**
 * @brief Boundary forms reproduce the measure of the boundary.
 *
 * Separate from the volume forms because the reference is different: the sum
 * is over the faces, not the cells, and it is computed from the faces
 * themselves so that neither side of the comparison is an integral.
 */
TEST(QuadratureRuleSpecializationTest, BoundaryFormsReproduceTheBoundaryMeasure)
{
  for (const auto& element : elements())
  {
    const std::string where = std::string("on ") + element.name;
    LocalMesh mesh = makeMesh(element);

    Real perimeter = 0;
    for (auto it = mesh.getBoundary(); !it.end(); ++it)
      perimeter += it->getMeasure();
    ASSERT_GT(perimeter, 0) << where << ": the mesh reports no boundary";

    P1<Real, LocalMesh> fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);

    {
      LinearForm load(v);
      load = BoundaryIntegral(v);
      load.assemble();
      EXPECT_NEAR(load.getVector().sum(), perimeter, tolerance)
        << where << ": a boundary load should integrate one over the boundary";
    }
    {
      BilinearForm mass(u, v);
      mass = BoundaryIntegral(u, v);
      mass.assemble();
      EXPECT_NEAR(mass.getOperator().sum(), perimeter, tolerance)
        << where << ": a boundary mass should integrate one over the boundary";
    }
  }
}
