/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
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
  };

  TEST_P(WeightedH1JacobianTest, MatchesOriginalEntryLoop)
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
