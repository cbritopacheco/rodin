/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <functional>
#include "Rodin/Geometry/ParametricTransformation.h"
#include "Rodin/Variational.h"
#include "../../../QuadratureReference.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  struct WeightedJacobianGeometry
  {
    const char* name;
    Polytope::Type type;
    bool embedded;
    bool curved;
  };

  /**
   * @brief Qualifies weighted vector Jacobians independently of kernel symmetry.
   *
   * ## Architecture
   *
   * The generic expression-entry oracle computes @f$(C J(u),J(v))@f$ without
   * production gradient reuse or triangular mirroring. Scalar real/complex
   * coefficients and a nonsymmetric component matrix exercise distinct
   * coefficient contracts. Every entry is checked after formula replacement
   * and rebinding, for equal and mixed degrees, real/complex fields, all volume
   * geometries, embedded cells and the empty reference derivative on points.
   */
  class WeightedH1JacobianTest : public ::testing::TestWithParam<WeightedJacobianGeometry>
  {
    protected:
      template <size_t TrialDegree, size_t TestDegree, class Scalar>
      static void check(const LocalMesh& mesh)
      {
        for (const size_t components : {size_t(2), size_t(3)})
        {
          SCOPED_TRACE(components);
          SCOPED_TRACE(TrialDegree);
          SCOPED_TRACE(TestDegree);
          H1<TrialDegree, Math::SpatialVector<Scalar>> trial(
            std::integral_constant<size_t, TrialDegree>{}, mesh, components);
          H1<TestDegree, Math::SpatialVector<Scalar>> test(
            std::integral_constant<size_t, TestDegree>{}, mesh, components);
          TrialFunction u(trial);
          TestFunction v(test);
          const auto checkCoefficient = [&](const auto& coefficient) {
            const auto integrand = Dot(coefficient * Jacobian(u), Jacobian(v));
            Integral integral(integrand);
            QuadratureReference<decltype(integral)> reference(integrand);
            using Value = typename decltype(reference)::Scalar;
            // The two entry algorithms differ only by floating-point association.
            constexpr Real RelativeEntryTolerance = 1e-12;
            for (const size_t order : {size_t(4), size_t(2), size_t(6)})
            {
              integral.setOrder(order);
              for (auto cell = mesh.getCell(); cell; ++cell)
              {
                integral.setPolytope(*cell);
                reference.assemble(*cell, order);
                const auto& expected = reference.getOperator();
                Math::Matrix<Value> actual(expected.rows(), expected.cols());
                for (Eigen::Index te = 0; te < actual.rows(); ++te)
                {
                  for (Eigen::Index tr = 0; tr < actual.cols(); ++tr)
                    actual(te, tr) = integral.integrate(tr, te);
                }
                const Real scale = std::max(Real(1), expected.norm());
                ASSERT_LE((actual - expected).norm(), RelativeEntryTolerance * scale);
                auto wrong = actual;
                wrong(0, 0) += scale;
                ASSERT_GT((wrong - expected).norm(), RelativeEntryTolerance * scale);
              }
            }
          };
          checkCoefficient(RealFunction(Real(2)));
          // A nonconstant positive viscosity qualifies quadrature accumulation,
          // independently of the component-factorization control below.
          checkCoefficient(
            RealFunction([](const Point& point) { return Real(2) + point.x(); }));
          if constexpr (std::is_same_v<Scalar, Complex>)
            checkCoefficient(ComplexFunction(Complex(Real(0.25), Real(0.5))));
          Math::SpatialMatrix<Scalar> matrix(components, components);
          matrix.setZero();
          for (size_t c = 0; c < components; ++c)
            matrix(c, c) = Scalar(1);
          matrix(0, 1) = Scalar(Real(0.5));
          matrix(1, 0) = Scalar(Real(-0.25));
          if constexpr (std::is_same_v<Scalar, Complex>)
            matrix(0, 1) += Complex(0, Real(0.25));
          checkCoefficient(MatrixFunction(matrix));
        }
      }

      /**
       * @brief Checks scalar-weighted component factorization exactly.
       *
       * The one-component operator supplies @f$S_{ba}@f$ independently of
       * component replication. Every entry must satisfy
       * @f$A_{bm+e,am+c}=\delta_{ec}S_{ba}@f$, including mixed degrees,
       * empty ranges, complex weights and changing coefficient/formula state.
       * For multisample rules the weight at the final logical index vanishes:
       * @f$\nu_{N_q-1}=0@f$. Completed entries must still be replicated even
       * when their final contribution is zero. No coordinate matching is used.
       *
       * @tparam TrialDegree Trial polynomial degree.
       * @tparam TestDegree Test polynomial degree.
       * @tparam Scalar Field and coefficient scalar type.
       * @param mesh Native, embedded or curved geometric fixture.
       */
      template <size_t TrialDegree, size_t TestDegree, class Scalar>
      static void checkBlocks(const LocalMesh& mesh)
      {
        using Range = Math::SpatialVector<Scalar>;
        H1<TrialDegree, Range> scalarTrial(
          std::integral_constant<size_t, TrialDegree>{}, mesh, 1);
        H1<TestDegree, Range> scalarTest(
          std::integral_constant<size_t, TestDegree>{}, mesh, 1);
        TrialFunction scalarU(scalarTrial);
        TestFunction scalarV(scalarTest);
        Real scale = 1;
        size_t evaluations = 0;
        size_t zeroFinalSamples = 0;
        const auto coefficient = [&]() {
          const auto valueAt = [&](const Point& point) {
            const Real x = point.getDimension() == 0 ? Real(0) : point.x();
            if constexpr (std::is_same_v<Scalar, Complex>)
              return scale * Complex(Real(2) + x, Real(1));
            else
              return scale * (Real(2) + x);
          };
          const auto callable = [&, valueAt](const auto& sample) {
            ++evaluations;
            using Sample = std::remove_cvref_t<decltype(sample)>;
            if constexpr (std::is_same_v<Sample, IntegrationPoint>)
            {
              const auto* formula = sample.getQuadratureFormula();
              if (formula && formula->getSize() > 1 &&
                sample.getIndex() + 1 == formula->getSize())
              {
                ++zeroFinalSamples;
                return Scalar(0);
              }
              return valueAt(sample.getPoint());
            }
            else
            {
              return valueAt(sample);
            }
          };
          if constexpr (std::is_same_v<Scalar, Complex>)
            return ComplexFunction(callable);
          else
            return RealFunction(callable);
        }();
        auto scalar = Integral(coefficient * Jacobian(scalarU), Jacobian(scalarV));
        for (const size_t components : {size_t(0), size_t(1), size_t(2), size_t(3)})
        {
          SCOPED_TRACE(TrialDegree);
          SCOPED_TRACE(TestDegree);
          SCOPED_TRACE(components);
          H1<TrialDegree, Range> trial(
            std::integral_constant<size_t, TrialDegree>{}, mesh, components);
          H1<TestDegree, Range> test(
            std::integral_constant<size_t, TestDegree>{}, mesh, components);
          TrialFunction u(trial);
          TestFunction v(test);
          auto actual = Integral(coefficient * Jacobian(u), Jacobian(v));
          for (size_t step = 0; step < 3; ++step)
          {
            scale = Real(step + 1);
            const size_t order = step == 1 ? 4 : 0;
            scalar.setOrder(order);
            actual.setOrder(order);
            for (auto cell = mesh.getCell(); cell; ++cell)
            {
              scalar.setPolytope(*cell);
              const size_t previous = evaluations;
              const size_t previousZero = zeroFinalSamples;
              actual.setPolytope(*cell);
              if (components == 0)
              {
                ASSERT_EQ(evaluations, previous);
              }
              else if (order == 4 && cell->getDimension() != 0)
              {
                ASSERT_EQ(zeroFinalSamples, previousZero + 1);
              }
              const size_t nr = scalarU.getDOFs(*cell), nt = scalarV.getDOFs(*cell);
              ASSERT_EQ(u.getDOFs(*cell), components * nr);
              ASSERT_EQ(v.getDOFs(*cell), components * nt);
              bool nonzero = false;
              for (size_t b = 0; b < nt; ++b)
              {
                for (size_t a = 0; a < nr; ++a)
                {
                  const Scalar expected = scalar.integrate(a, b);
                  nonzero = nonzero || expected != Scalar(0);
                  for (size_t e = 0; e < components; ++e)
                  {
                    for (size_t c = 0; c < components; ++c)
                    {
                      ASSERT_EQ(actual.integrate(a * components + c, b * components + e),
                        c == e ? expected : Scalar(0));
                    }
                  }
                }
              }
              if (cell->getDimension() != 0)
              {
                ASSERT_TRUE(nonzero);
              }
            }
          }
        }
      }

      /**
       * @brief Constructs the parameter's native or embedded geometric fixture.
       * @returns Mesh with explicit lower-dimensional incidence and optional curvature.
       */
      LocalMesh makeMesh() const
      {
        const auto parameter = GetParam();
        auto mesh = [&]() {
          if (parameter.type == Polytope::Type::Point)
          {
            auto builder = LocalMesh::Builder();
            builder.initialize(3).nodes(1).vertex({0, 0, 0});
            return builder.finalize();
          }
          if (parameter.embedded)
          {
            auto builder = LocalMesh::Builder();
            builder.initialize(3).nodes(3);
            builder.vertex({0, 0, 0}).vertex({1, 0, 1});
            if (parameter.type == Polytope::Type::Segment)
              return builder.polytope(parameter.type, {0, 1}).finalize();
            builder.vertex({0, 1, 1});
            return builder.polytope(parameter.type, {0, 1, 2}).finalize();
          }
          const size_t dimension = Polytope::Traits(parameter.type).getDimension();
          Array<size_t> shape(dimension);
          std::fill(shape.begin(), shape.end(), 2);
          return LocalMesh::UniformGrid(parameter.type, shape);
        }();
        const size_t dimension = mesh.getDimension();
        for (size_t d = 1; d <= dimension; ++d)
        {
          for (size_t lower = 0; lower < d; ++lower)
            mesh.getConnectivity().compute(d, lower);
        }
        if (parameter.curved)
        {
          const RealH1Element<2> element(parameter.type);
          PointCloud nodes(3, element.getCount());
          for (size_t a = 0; a < element.getCount(); ++a)
          {
            const auto& r = element.getNode(a);
            nodes(0, a) = r[0];
            nodes(1, a) = r[1];
            nodes(2, a) = r[0] + r[1] + Real(0.25) * r[0] * r[1];
          }
          mesh.setPolytopeTransformation(
            {2, 0}, new ParametricTransformation<RealH1Element<2>>(nodes, element));
        }
        return mesh;
      }
  };

  TEST_P(WeightedH1JacobianTest, MatchesOriginalEntryLoop)
  {
    auto mesh = makeMesh();
    check<1, 1, Real>(mesh);
    check<2, 2, Real>(mesh);
    check<1, 2, Real>(mesh);
    check<2, 1, Real>(mesh);
    check<2, 3, Real>(mesh);
    check<3, 2, Real>(mesh);
    check<1, 1, Complex>(mesh);
    check<2, 2, Complex>(mesh);
    check<1, 2, Complex>(mesh);
    check<2, 1, Complex>(mesh);
    check<2, 3, Complex>(mesh);
    check<3, 2, Complex>(mesh);
  }

  /**
   * @brief Scalar weights preserve all component blocks and empty-range semantics.
   *
   * Every supported degree is checked in both scalar fields on all parameter
   * geometries. Mixed pairs also exercise rectangular scalar blocks.
   */
  TEST_P(WeightedH1JacobianTest, ScalarCoefficientComponentBlocks)
  {
    auto mesh = makeMesh();
    const auto degrees = [&]<class Scalar>() {
      checkBlocks<1, 1, Scalar>(mesh);
      checkBlocks<2, 2, Scalar>(mesh);
      checkBlocks<3, 3, Scalar>(mesh);
      checkBlocks<4, 4, Scalar>(mesh);
      checkBlocks<5, 5, Scalar>(mesh);
      checkBlocks<6, 6, Scalar>(mesh);
      checkBlocks<1, 2, Scalar>(mesh);
      checkBlocks<2, 1, Scalar>(mesh);
      checkBlocks<2, 3, Scalar>(mesh);
      checkBlocks<3, 2, Scalar>(mesh);
    };
    degrees.template operator()<Real>();
    degrees.template operator()<Complex>();
  }

  /**
   * @brief Nested coefficient evaluation preserves invocation-local gradient shapes.
   *
   * Different integrator owners of the same C++ type share thread-local scratch.
   * A segment/hexahedron callback changes its active physical dimension and
   * scalar basis count. The outer operator must equal its non-nested control
   * entry for entry, in both dimension directions and after repeated binding.
   * Thus @f$A_{\mathrm{nested}}=A_{\mathrm{plain}}@f$ is checked exactly,
   * independently of the selected quadrature rule's integration accuracy.
   * Equal and mixed degrees cover shared and separate trial/test storage;
   * scalar and nonsymmetric matrix coefficients cover both contraction paths.
   * Real scalar-weighted cases cover degrees one through six; complex and
   * matrix-weighted cases cover degrees one through three and both mixed pairs.
   * This does not permit recursive mutation of the same integrator owner.
   */
  TEST(WeightedH1JacobianRegression, ReentrantCoefficientScratchLifetime)
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
    const auto check = [&]<size_t TrialDegree, size_t TestDegree, class Scalar,
                         class CoefficientScalar, bool Matrix>() {
      using Value = std::conditional_t<Matrix, Math::SpatialMatrix<CoefficientScalar>,
        CoefficientScalar>;
      using Callable = std::function<Value(const Point&)>;
      const Callable constant = [](const Point&) {
        if constexpr (Matrix)
        {
          Value value(2, 2);
          value(0, 0) = CoefficientScalar(2);
          value(1, 1) = CoefficientScalar(3);
          value(0, 1) = CoefficientScalar(1);
          value(1, 0) = CoefficientScalar(-1);
          if constexpr (std::is_same_v<CoefficientScalar, Complex>)
            value(0, 1) += Complex(0, 1);
          return value;
        }
        else if constexpr (std::is_same_v<CoefficientScalar, Complex>)
          return Complex(2, 1);
        else
          return Real(2);
      };
      const auto coefficient = [](const Callable& callable) {
        if constexpr (Matrix)
          return MatrixFunction(2, 2, callable);
        else if constexpr (std::is_same_v<CoefficientScalar, Complex>)
          return ComplexFunction(callable);
        else
          return RealFunction(callable);
      };
      for (const bool reverse : {false, true})
      {
        SCOPED_TRACE(reverse);
        SCOPED_TRACE(TrialDegree);
        SCOPED_TRACE(TestDegree);
        SCOPED_TRACE(Matrix);
        auto& outerMesh = reverse ? segment : hex;
        auto& innerMesh = reverse ? hex : segment;
        using Range = Math::SpatialVector<Scalar>;
        H1<TrialDegree, Range> outerTrial(
          std::integral_constant<size_t, TrialDegree>{}, outerMesh, 2);
        H1<TestDegree, Range> outerTest(
          std::integral_constant<size_t, TestDegree>{}, outerMesh, 2);
        H1<TrialDegree, Range> innerTrial(
          std::integral_constant<size_t, TrialDegree>{}, innerMesh, 2);
        H1<TestDegree, Range> innerTest(
          std::integral_constant<size_t, TestDegree>{}, innerMesh, 2);
        TrialFunction outerU(outerTrial), innerU(innerTrial);
        TestFunction outerV(outerTest), innerV(innerTest);
        auto inner = Integral(coefficient(constant) * Jacobian(innerU), Jacobian(innerV));
        inner.setOrder(0);
        const auto innerCell = *innerMesh.getCell();
        size_t callbacks = 0;
        const Callable nested = [&](const Point& point) {
          ++callbacks;
          inner.setPolytope(innerCell);
          return constant(point);
        };
        auto actual = Integral(coefficient(nested) * Jacobian(outerU), Jacobian(outerV));
        auto expected =
          Integral(coefficient(constant) * Jacobian(outerU), Jacobian(outerV));
        static_assert(std::is_same_v<decltype(actual), decltype(inner)>);
        static_assert(decltype(actual)::Specialized);
        const auto cell = *outerMesh.getCell();
        for (const size_t order : {size_t(0), size_t(4), size_t(0)})
        {
          actual.setOrder(order);
          expected.setOrder(order);
          expected.setPolytope(cell);
          const size_t nr = outerU.getDOFs(cell), nt = outerV.getDOFs(cell);
          Math::Matrix<Scalar> initial(nt, nr);
          bool nonzero = false;
          for (size_t te = 0; te < nt; ++te)
          {
            for (size_t tr = 0; tr < nr; ++tr)
            {
              initial(te, tr) = expected.integrate(tr, te);
              nonzero = nonzero || initial(te, tr) != Scalar(0);
            }
          }
          ASSERT_TRUE(nonzero);
          const size_t previous = callbacks;
          ASSERT_NO_THROW(actual.setPolytope(cell));
          ASSERT_GT(callbacks, previous);
          for (size_t te = 0; te < nt; ++te)
          {
            for (size_t tr = 0; tr < nr; ++tr)
              ASSERT_EQ(actual.integrate(tr, te), initial(te, tr));
          }
        }
      }
    };
    const auto degrees = [&]<class Scalar, class CoefficientScalar, bool Matrix>() {
      check.template operator()<1, 1, Scalar, CoefficientScalar, Matrix>();
      check.template operator()<2, 2, Scalar, CoefficientScalar, Matrix>();
      check.template operator()<2, 3, Scalar, CoefficientScalar, Matrix>();
      check.template operator()<3, 2, Scalar, CoefficientScalar, Matrix>();
    };
    degrees.template operator()<Real, Real, false>();
    degrees.template operator()<Complex, Real, false>();
    degrees.template operator()<Complex, Complex, false>();
    degrees.template operator()<Real, Real, true>();
    degrees.template operator()<Complex, Complex, true>();
    check.template operator()<4, 4, Real, Real, false>();
    check.template operator()<5, 5, Real, Real, false>();
    check.template operator()<6, 6, Real, Real, false>();
  }

  using G = Polytope::Type;
  INSTANTIATE_TEST_SUITE_P(AllGeometries, WeightedH1JacobianTest,
    ::testing::Values(WeightedJacobianGeometry{"Point", G::Point, false, false},
      WeightedJacobianGeometry{"Segment", G::Segment, false, false},
      WeightedJacobianGeometry{"Triangle", G::Triangle, false, false},
      WeightedJacobianGeometry{"Quadrilateral", G::Quadrilateral, false, false},
      WeightedJacobianGeometry{"Tetrahedron", G::Tetrahedron, false, false},
      WeightedJacobianGeometry{"Pyramid", G::Pyramid, false, false},
      WeightedJacobianGeometry{"Hexahedron", G::Hexahedron, false, false},
      WeightedJacobianGeometry{"Wedge", G::Wedge, false, false},
      WeightedJacobianGeometry{"EmbeddedSegment", G::Segment, true, false},
      WeightedJacobianGeometry{"EmbeddedTriangle", G::Triangle, true, false},
      WeightedJacobianGeometry{"CurvedEmbeddedTriangle", G::Triangle, true, true}),
    [](const auto& info) { return std::string(info.param.name); });
}
