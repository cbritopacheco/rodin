/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <functional>
#include "Rodin/Variational.h"
#include "../../../QuadratureReference.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  /**
   * @brief Qualifies coefficient linearity and ambient-gradient convection mapping.
   *
   * ## Mathematical foundation
   *
   * For real-valued scalar nodal bases, component-diagonal entries are
   * @f$A_{bc,ac}=\int_K(\nabla\phi_a)^T b\,\phi_b\,dx@f$.
   * Consequently @f$A(ib)=iA(b)@f$: conjugation belongs to the test slot,
   * not the coefficient. This identity is exact independently of quadrature
   * exactness. The entry-loop oracle independently evaluates the expression.
   *
   * ## Scope
   *
   * Complex-field checks use degrees one through six, mixed pairs (1,2), (2,1),
   * (2,3), (3,2), and component counts zero through three. Real-field entry-loop
   * checks use degrees one through three and two components. All eight cell
   * families are included; complex embedded segment/triangle cells and an
   * ambient-three Point also exercise rectangular mappings.
   * One-point tests qualify the local operator, not integration exactness or
   * a PDE convergence rate. The generic entry-loop oracle covers small
   * operators; the higher-degree complex cases use the exact phase identity.
   * A separate component-block identity checks real and complex degrees one,
   * two, three and six, including the four mixed pairs above, without a tolerance.
   */
  class ConvectionQuadratureTest : public ::testing::Test
  {
    protected:
      /**
       * @brief Checks coefficient linearity and an independent entry loop.
       * @param trial Complex trial space with real-valued scalar basis.
       * @param test Complex test space with the same component count.
       * @param order Requested quadrature order.
       * @param reference Whether to evaluate the complete entry-loop oracle.
       */
      void check(auto& trial, auto& test, size_t order, bool reference)
      {
        const auto& mesh = trial.getMesh();
        const auto cell = *mesh.getCell();
        const size_t ambient = mesh.getSpaceDimension();
        const VectorFunction realDirection(ambient, [ambient](const Point&) {
          Math::SpatialVector<Real> value(ambient);
          for (size_t j = 0; j < ambient; ++j)
            value[j] = Real(j + 1);
          return value;
        });
        const VectorFunction imaginaryDirection(ambient, [ambient](const Point&) {
          Math::SpatialVector<Complex> value(ambient);
          for (size_t j = 0; j < ambient; ++j)
            value[j] = Complex(0, Real(j + 1));
          return value;
        });
        TrialFunction u(trial);
        TestFunction v(test);
        auto realForm = Integral(Jacobian(u) * realDirection, v);
        const auto integrand = Dot(Jacobian(u) * imaginaryDirection, v);
        Integral imaginaryForm(integrand);
        static_assert(decltype(realForm)::Specialized);
        static_assert(decltype(imaginaryForm)::Specialized);
        realForm.setOrder(order);
        imaginaryForm.setOrder(order);
        realForm.setPolytope(cell);
        imaginaryForm.setPolytope(cell);
        const size_t nr = u.getDOFs(cell);
        const size_t nt = v.getDOFs(cell);
        Math::Matrix<Complex> actual(nt, nr);
        bool nonzero = false;
        for (size_t te = 0; te < nt; ++te)
        {
          for (size_t tr = 0; tr < nr; ++tr)
          {
            const Complex expected = Complex(0, 1) * realForm.integrate(tr, te);
            actual(te, tr) = imaginaryForm.integrate(tr, te);
            ASSERT_EQ(actual(te, tr), expected) << "entry " << te << ',' << tr;
            nonzero = nonzero || expected != Complex(0);
          }
        }
        if (cell.getDimension() == 0 || trial.getVectorDimension() == 0)
        {
          EXPECT_FALSE(nonzero);
        }
        else
        {
          EXPECT_TRUE(nonzero);
        }
        if (reference)
        {
          QuadratureReference<decltype(imaginaryForm)> oracle(integrand);
          oracle.assemble(cell, order);
          // Expression and specialized loops associate products differently.
          constexpr Real RelativeEntryLoopBudget = 1e-12;
          const auto& expected = oracle.getOperator();
          const Real scale = std::max(Real(1), expected.norm());
          EXPECT_LE((actual - expected).norm(), RelativeEntryLoopBudget * scale);
        }
      }

      /**
       * @brief Exercises H1 orders and mixed-degree directions on one mesh.
       * @param mesh Native, embedded or zero-dimensional cell mesh.
       */
      void checkOrders(LocalMesh& mesh)
      {
        for (size_t d = 1; d <= mesh.getDimension(); ++d)
        {
          for (size_t lower = 0; lower < d; ++lower)
            mesh.getConnectivity().compute(d, lower);
        }
        for (const size_t components : {size_t(0), size_t(1), size_t(2), size_t(3)})
        {
          SCOPED_TRACE(components);
          P1<Math::SpatialVector<Complex>> linear(mesh, components);
          check(linear, linear, 0, true);
          const auto degrees = [&]<size_t TrialDegree, size_t TestDegree>() {
            SCOPED_TRACE(TrialDegree);
            SCOPED_TRACE(TestDegree);
            H1<TrialDegree, Math::SpatialVector<Complex>> trial(
              std::integral_constant<size_t, TrialDegree>{}, mesh, components);
            H1<TestDegree, Math::SpatialVector<Complex>> test(
              std::integral_constant<size_t, TestDegree>{}, mesh, components);
            check(trial, test, 0, TrialDegree <= 2 && TestDegree <= 2);
            if constexpr (TrialDegree <= 2 && TestDegree <= 2)
              check(trial, test, 4, true);
          };
          degrees.template operator()<1, 1>();
          degrees.template operator()<2, 2>();
          degrees.template operator()<3, 3>();
          degrees.template operator()<4, 4>();
          degrees.template operator()<5, 5>();
          degrees.template operator()<6, 6>();
          degrees.template operator()<1, 2>();
          degrees.template operator()<2, 1>();
          degrees.template operator()<2, 3>();
          degrees.template operator()<3, 2>();
        }
      }
  };

  TEST_F(ConvectionQuadratureTest, ComplexCoefficientLinearity)
  {
    for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(int(geometry));
      const Array<size_t> shape =
        Array<size_t>::Constant(Polytope::Traits(geometry).getDimension(), 2);
      auto mesh = LocalMesh::UniformGrid(geometry, shape);
      checkOrders(mesh);
    }
  }

  /**
   * @brief Checks same-type re-entrant coefficient evaluation across ambient dimensions.
   *
   * A coefficient binds a different integrator of the exact same C++ type on
   * a different mesh. Its thread-local working storage is consequently shared.
   * Both ambient-dimension directions must preserve the outer operator exactly.
   * The scalar basis degree and component count remain fixed during each bind;
   * this is not recursive mutation of the same integrator or finite element space.
   */
  TEST_F(ConvectionQuadratureTest, ReentrantCoefficientScratchLifetime)
  {
    auto segment =
      LocalMesh::UniformGrid(Polytope::Type::Segment, Array<size_t>::Constant(1, 2));
    auto hex =
      LocalMesh::UniformGrid(Polytope::Type::Hexahedron, Array<size_t>::Constant(3, 2));
    for (auto* mesh : {&segment, &hex})
    {
      for (size_t d = 1; d <= mesh->getDimension(); ++d)
      {
        for (size_t lower = 0; lower < d; ++lower)
          mesh->getConnectivity().compute(d, lower);
      }
    }
    const auto checkScalar = [&]<class Scalar>() {
      using Value = Math::SpatialVector<Scalar>;
      using Callable = std::function<Value(const Point&)>;
      const Callable constant = [](const Point& point) {
        Value value(point.getDimension());
        for (size_t j = 0; j < point.getDimension(); ++j)
        {
          if constexpr (std::is_same_v<Scalar, Complex>)
            value[j] = Complex(Real(j + 1), Real(j + 2));
          else
            value[j] = Real(j + 1);
        }
        return value;
      };
      for (const bool reverse : {false, true})
      {
        SCOPED_TRACE(reverse);
        auto& outerMesh = reverse ? segment : hex;
        auto& innerMesh = reverse ? hex : segment;
        H1<2, Value> outerSpace(std::integral_constant<size_t, 2>{}, outerMesh, 2);
        H1<2, Value> innerSpace(std::integral_constant<size_t, 2>{}, innerMesh, 2);
        TrialFunction outerU(outerSpace), innerU(innerSpace);
        TestFunction outerV(outerSpace), innerV(innerSpace);
        const VectorFunction innerDirection(innerMesh.getSpaceDimension(), constant);
        auto inner = Integral(Jacobian(innerU) * innerDirection, innerV);
        inner.setOrder(0);
        const auto innerCell = *innerMesh.getCell();
        const Callable reentrant = [&](const Point& point) {
          inner.setPolytope(innerCell);
          return constant(point);
        };
        const VectorFunction outerDirection(outerMesh.getSpaceDimension(), reentrant);
        const VectorFunction plainDirection(outerMesh.getSpaceDimension(), constant);
        auto actual = Integral(Jacobian(outerU) * outerDirection, outerV);
        auto expected = Integral(Jacobian(outerU) * plainDirection, outerV);
        static_assert(std::is_same_v<decltype(actual), decltype(inner)>);
        actual.setOrder(4);
        expected.setOrder(4);
        const auto cell = *outerMesh.getCell();
        expected.setPolytope(cell);
        const size_t n = outerU.getDOFs(cell);
        Math::Matrix<Scalar> initial(n, n);
        for (size_t te = 0; te < n; ++te)
        {
          for (size_t tr = 0; tr < n; ++tr)
            initial(te, tr) = expected.integrate(tr, te);
        }
        ASSERT_NO_THROW(actual.setPolytope(cell));
        for (size_t te = 0; te < n; ++te)
        {
          for (size_t tr = 0; tr < n; ++tr)
            EXPECT_EQ(actual.integrate(tr, te), initial(te, tr));
        }
      }
    };
    checkScalar.template operator()<Real>();
    checkScalar.template operator()<Complex>();
  }

  TEST_F(ConvectionQuadratureTest, EmbeddedAndPointMappings)
  {
    auto point = LocalMesh::Builder().initialize(3).nodes(1).vertex({1, 2, 3}).finalize();
    checkOrders(point);
    auto segment = LocalMesh::Builder()
                     .initialize(3)
                     .nodes(2)
                     .vertex({0, 0, 0})
                     .vertex({1, 1, 1})
                     .polytope(Polytope::Type::Segment, {{0, 1}})
                     .finalize();
    checkOrders(segment);
    auto triangle = LocalMesh::Builder()
                      .initialize(3)
                      .nodes(3)
                      .vertex({0, 0, 0})
                      .vertex({1, 0, 1})
                      .vertex({0, 1, 1})
                      .polytope(Polytope::Type::Triangle, {{0, 1, 2}})
                      .finalize();
    checkOrders(triangle);
  }

  /**
   * @brief Checks real fields and directions through the independent expression loop.
   *
   * The form-language Dot constraint currently rejects a complex direction
   * paired with real vector trial/test ranges. This test does not widen that
   * separate expression contract or claim coverage of the rejected pairing.
   */
  TEST_F(ConvectionQuadratureTest, RealFieldsAndDirections)
  {
    for (const auto geometry :
      {Polytope::Type::Point, Polytope::Type::Segment, Polytope::Type::Triangle,
        Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
        Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(int(geometry));
      const size_t dimension = Polytope::Traits(geometry).getDimension();
      auto mesh = [&]() {
        if (dimension == 0)
          return LocalMesh::Builder().initialize(3).nodes(1).vertex({0, 0, 0}).finalize();
        return LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 2));
      }();
      for (size_t d = 1; d <= dimension; ++d)
      {
        for (size_t lower = 0; lower < d; ++lower)
          mesh.getConnectivity().compute(d, lower);
      }
      const auto checkReal = [&](auto& space) {
        TrialFunction u(space);
        TestFunction v(space);
        const size_t ambient = mesh.getSpaceDimension();
        const VectorFunction direction(ambient, [ambient](const Point&) {
          Math::SpatialVector<Real> value(ambient);
          for (size_t j = 0; j < ambient; ++j)
            value[j] = Real(j + 1);
          return value;
        });
        const auto integrand = Dot(Jacobian(u) * direction, v);
        Integral integral(integrand);
        static_assert(decltype(integral)::Specialized);
        const auto cell = *mesh.getCell();
        const size_t n = u.getDOFs(cell);
        for (const size_t order : {size_t(0), size_t(4)})
        {
          integral.setOrder(order);
          integral.setPolytope(cell);
          Math::Matrix<Real> actual(n, n);
          for (size_t te = 0; te < n; ++te)
          {
            for (size_t tr = 0; tr < n; ++tr)
              actual(te, tr) = integral.integrate(tr, te);
          }
          QuadratureReference<decltype(integral)> reference(integrand);
          reference.assemble(cell, order);
          constexpr Real RelativeEntryLoopBudget = 1e-12;
          const auto& expected = reference.getOperator();
          EXPECT_LE((actual - expected).norm(),
            RelativeEntryLoopBudget * std::max(Real(1), expected.norm()));
          integral.setPolytope(cell);
          for (size_t te = 0; te < n; ++te)
          {
            for (size_t tr = 0; tr < n; ++tr)
              EXPECT_EQ(integral.integrate(tr, te), actual(te, tr));
          }
        }
      };
      P1<Math::SpatialVector<Real>> linear(mesh, 2);
      checkReal(linear);
      const auto degree = [&]<size_t Degree>() {
        SCOPED_TRACE(Degree);
        H1<Degree, Math::SpatialVector<Real>> space(
          std::integral_constant<size_t, Degree>{}, mesh, 2);
        checkReal(space);
      };
      degree.template operator()<1>();
      degree.template operator()<2>();
      degree.template operator()<3>();
    }
  }

  /**
   * @brief Checks every H1 component block against its one-component operator.
   *
   * H1 vector bases use the component-expanded local order. For equal component
   * counts and possibly different scalar degrees,
   * @f$A_{bm+c,am+e}=\delta_{ce}S_{ba}@f$. Each matching block must reproduce
   * the one-component operator exactly; all cross-component entries vanish.
   * This checks rectangular local operators and matrix addressing independently
   * of coefficient-phase identities. A zero component count gives an empty
   * operator. Point cells have zero derivatives rather than a refinement limit.
   */
  TEST_F(ConvectionQuadratureTest, ComponentDiagonalEntries)
  {
    for (const auto geometry :
      {Polytope::Type::Point, Polytope::Type::Segment, Polytope::Type::Triangle,
        Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
        Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(int(geometry));
      const size_t dimension = Polytope::Traits(geometry).getDimension();
      auto mesh = [&]() {
        if (dimension == 0)
          return LocalMesh::Builder().initialize(3).nodes(1).vertex({0, 0, 0}).finalize();
        return LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 2));
      }();
      for (size_t d = 1; d <= dimension; ++d)
      {
        for (size_t lower = 0; lower < d; ++lower)
          mesh.getConnectivity().compute(d, lower);
      }
      const auto cell = *mesh.getCell();
      const auto check = [&]<class Scalar, size_t TrialDegree, size_t TestDegree>() {
        SCOPED_TRACE(TrialDegree);
        SCOPED_TRACE(TestDegree);
        using TrialSpace = H1<TrialDegree, Math::SpatialVector<Scalar>>;
        using TestSpace = H1<TestDegree, Math::SpatialVector<Scalar>>;
        TrialSpace scalarTrial(std::integral_constant<size_t, TrialDegree>{}, mesh, 1);
        TestSpace scalarTest(std::integral_constant<size_t, TestDegree>{}, mesh, 1);
        TrialFunction scalarU(scalarTrial);
        TestFunction scalarV(scalarTest);
        const size_t ambient = mesh.getSpaceDimension();
        const VectorFunction direction(ambient, [ambient](const Point&) {
          Math::SpatialVector<Scalar> value(ambient);
          for (size_t axis = 0; axis < ambient; ++axis)
          {
            if constexpr (std::is_same_v<Scalar, Complex>)
              value[axis] = Complex(Real(axis + 1), Real(axis + 2));
            else
              value[axis] = Real(axis + 1);
          }
          return value;
        });
        auto scalar = Integral(Jacobian(scalarU) * direction, scalarV);
        scalar.setOrder(0);
        scalar.setPolytope(cell);
        const size_t trialCount = scalarU.getDOFs(cell);
        const size_t testCount = scalarV.getDOFs(cell);
        Math::Matrix<Scalar> expected(testCount, trialCount);
        for (size_t te = 0; te < testCount; ++te)
        {
          for (size_t tr = 0; tr < trialCount; ++tr)
            expected(te, tr) = scalar.integrate(tr, te);
        }
        for (const size_t components : {size_t(0), size_t(1), size_t(2), size_t(3)})
        {
          SCOPED_TRACE(components);
          TrialSpace trial(
            std::integral_constant<size_t, TrialDegree>{}, mesh, components);
          TestSpace test(std::integral_constant<size_t, TestDegree>{}, mesh, components);
          TrialFunction u(trial);
          TestFunction v(test);
          auto integral = Integral(Jacobian(u) * direction, v);
          integral.setOrder(0);
          integral.setPolytope(cell);
          ASSERT_EQ(u.getDOFs(cell), components * trialCount);
          ASSERT_EQ(v.getDOFs(cell), components * testCount);
          for (size_t te = 0; te < components * testCount; ++te)
          {
            for (size_t tr = 0; tr < components * trialCount; ++tr)
            {
              const Scalar value = te % components == tr % components
                ? expected(te / components, tr / components)
                : Scalar(0);
              ASSERT_EQ(integral.integrate(tr, te), value) << te << ',' << tr;
            }
          }
        }
      };
      const auto fields = [&]<class Scalar>() {
        check.template operator()<Scalar, 1, 1>();
        check.template operator()<Scalar, 2, 2>();
        check.template operator()<Scalar, 3, 3>();
        check.template operator()<Scalar, 6, 6>();
        check.template operator()<Scalar, 1, 2>();
        check.template operator()<Scalar, 2, 1>();
        check.template operator()<Scalar, 2, 3>();
        check.template operator()<Scalar, 3, 2>();
      };
      fields.template operator()<Real>();
      fields.template operator()<Complex>();
    }
  }

  /**
   * @brief Checks independent-expression reference copies before and after binding.
   *
   * The specialized directional integrator need not expose its expression.
   * Reference copies retain the explicitly supplied expression and order policy,
   * but bind their own mapped samples. The comparison is exact because all
   * objects use the same reference entry algorithm and accumulation order.
   */
  TEST_F(ConvectionQuadratureTest, ExplicitReferenceExpressionCopies)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    mesh.getConnectivity().compute(2, 0);
    mesh.getConnectivity().compute(2, 1);
    H1<2, Math::SpatialVector<Complex>> space(
      std::integral_constant<size_t, 2>{}, mesh, 3);
    TrialFunction u(space);
    TestFunction v(space);
    const VectorFunction direction(size_t(2), [](const Point&) {
      return Math::SpatialVector<Complex>{Complex(1, 2), Complex(3, 4)};
    });
    const auto integrand = Dot(Jacobian(u) * direction, v);
    Integral integral(integrand);
    integral.setOrder(4);
    ReferenceIntegral reference(integral, integrand);
    const auto before = std::unique_ptr<decltype(reference)>(reference.copy());
    const auto cell = *mesh.getCell();
    reference.setPolytope(cell);
    const auto after = std::unique_ptr<decltype(reference)>(reference.copy());
    ASSERT_GT(mesh.getPolytopeCount(2), size_t(1));
    const auto other = mesh.getPolytope(2, 1);
    ReferenceIntegral expected(integral, integrand);
    expected.setPolytope(*other);
    for (auto* copy : {before.get(), after.get()})
    {
      ASSERT_TRUE(copy->getOrder(*other));
      EXPECT_EQ(*copy->getOrder(*other), size_t(4));
      copy->setPolytope(*other);
    }
    reference.setPolytope(*other);
    reference.setPolytope(cell);
    for (auto* copy : {before.get(), after.get()})
    {
      EXPECT_EQ(copy->getPolytope().getIndex(), other->getIndex());
      bool nonzero = false;
      for (size_t te = 0; te < v.getDOFs(*other); ++te)
      {
        for (size_t tr = 0; tr < u.getDOFs(*other); ++tr)
        {
          const Complex value = expected.integrate(tr, te);
          nonzero = nonzero || value != Complex(0);
          EXPECT_EQ(copy->integrate(tr, te), value);
        }
      }
      EXPECT_TRUE(nonzero);
    }
  }

  /**
   * @brief Checks the analytic P1 unit-segment operator for a complex direction.
   *
   * For @f$b=2+3i@f$, the two trial derivatives are @f$-1,+1@f$ and
   * each test basis integrates to @f$1/2@f$. Thus each matching-component
   * entry is exactly @f$\pm(1+3i/2)@f$; off-component entries vanish.
   */
  TEST_F(ConvectionQuadratureTest, AnalyticP1Segment)
  {
    auto mesh = LocalMesh::Builder()
                  .initialize(1)
                  .nodes(2)
                  .vertex({0})
                  .vertex({1})
                  .polytope(Polytope::Type::Segment, {{0, 1}})
                  .finalize();
    mesh.getConnectivity().compute(1, 0);
    P1<Math::SpatialVector<Complex>> space(mesh, 3);
    TrialFunction u(space);
    TestFunction v(space);
    const VectorFunction direction(size_t(1), [](const Point&) {
      Math::SpatialVector<Complex> value(1);
      value[0] = Complex(2, 3);
      return value;
    });
    auto integral = Integral(Jacobian(u) * direction, v);
    integral.setOrder(0);
    integral.setPolytope(*mesh.getCell());
    BilinearForm assembled(u, v);
    assembled = integral;
    assembled.assemble();
    const Math::Matrix<Complex> global(assembled.getOperator());
    ASSERT_EQ(global.rows(), 6);
    ASSERT_EQ(global.cols(), 6);
    const auto cell = *mesh.getCell();
    const std::pair<size_t, Index> owner{cell.getDimension(), cell.getIndex()};
    std::array<bool, 6> visited{};
    for (size_t local = 0; local < 6; ++local)
    {
      const Index index = space.getGlobalIndex(owner, local);
      ASSERT_LT(index, visited.size());
      EXPECT_FALSE(visited[index]);
      visited[index] = true;
    }
    for (size_t te = 0; te < 6; ++te)
    {
      for (size_t tr = 0; tr < 6; ++tr)
      {
        const Complex expected = tr % 3 == te % 3
          ? (tr / 3 == 0 ? -Complex(1, 1.5) : Complex(1, 1.5))
          : Complex(0);
        EXPECT_EQ(integral.integrate(tr, te), expected);
        EXPECT_EQ(
          global(space.getGlobalIndex(owner, te), space.getGlobalIndex(owner, tr)),
          expected);
      }
    }
  }
}
