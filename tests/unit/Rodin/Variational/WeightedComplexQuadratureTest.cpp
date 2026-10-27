/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Scalar coefficient promotion on complex P1/H1 kernels.
 * The identity Integral(c A(u), A(v)) = c Integral(A(u), A(v))
 * is checked entrywise for real and genuinely complex constants.
 * A nonzero reference trace prevents an omitted operator from passing.
 */

#include <gtest/gtest.h>
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  struct WeightedGeometry
  {
      const char* name;
      Polytope::Type type;
      size_t dim;
  };

  template <class Trial, class Test>
  void checkWeighted(const Trial& trial, const Test& test, const Polytope& cell)
  {
    auto reference = Integral(trial, test);
    auto realWeighted = Integral(0.25 * trial, test);
    auto complexWeighted = Integral(Complex(0.25, 0.5) * trial, test);
    reference.setOrder(8);
    realWeighted.setOrder(8);
    complexWeighted.setOrder(8);
    reference.setPolytope(cell);
    realWeighted.setPolytope(cell);
    complexWeighted.setPolytope(cell);
    const size_t count = trial.getDOFs(cell);
    Real trace = 0;
    for (size_t i = 0; i < count; ++i)
    {
      trace += std::real(reference.integrate(i, i));
      for (size_t j = 0; j < count; ++j)
      {
        const auto entry = reference.integrate(j, i);
        EXPECT_LT(std::abs(realWeighted.integrate(j, i) - 0.25 * entry), 1e-12);
        EXPECT_LT(
          std::abs(complexWeighted.integrate(j, i) - Complex(0.25, 0.5) * entry), 1e-12);
      }
    }
    EXPECT_GT(trace, 0);
  }

  template <class FES>
  void checkSpace(const FES& space, const Polytope& cell)
  {
    TrialFunction u(space);
    TestFunction v(space);
    checkWeighted(u, v, cell);
    if constexpr (FormLanguage::IsVectorRange<
                    typename FormLanguage::Traits<FES>::RangeType>::Value)
      checkWeighted(Jacobian(u), Jacobian(v), cell);
    else
      checkWeighted(Grad(u), Grad(v), cell);
  }

  template <size_t K>
  void checkH1(const LocalMesh& mesh, const Polytope& cell)
  {
    SCOPED_TRACE(::testing::Message() << "H1 degree=" << K);
    const H1<K, Complex> scalar(std::integral_constant<size_t, K>{}, mesh);
    const H1<K, Math::SpatialVector<Complex>> vector(
      std::integral_constant<size_t, K>{}, mesh, mesh.getDimension());
    checkSpace(scalar, cell);
    checkSpace(vector, cell);
  }

  class WeightedComplexQuadratureTest : public ::testing::TestWithParam<WeightedGeometry>
  {};

  TEST_P(WeightedComplexQuadratureTest, ScalarAndVectorP1H1PromoteCoefficients)
  {
    const auto geometry = GetParam();
    Array<size_t> shape(geometry.dim);
    std::fill(shape.begin(), shape.end(), 2);
    auto mesh = LocalMesh::UniformGrid(geometry.type, shape);
    for (size_t d = 1; d <= geometry.dim; ++d)
    {
      mesh.getConnectivity().compute(d, 0);
      mesh.getConnectivity().compute(geometry.dim, d);
      for (size_t lower = 1; lower < d; ++lower)
        mesh.getConnectivity().compute(d, lower);
    }
    auto cell = mesh.getCell();
    const P1<Complex> scalar(mesh);
    const P1<Math::SpatialVector<Complex>> vector(mesh, geometry.dim);
    checkSpace(scalar, *cell);
    checkSpace(vector, *cell);
    checkH1<1>(mesh, *cell);
    checkH1<2>(mesh, *cell);
    checkH1<3>(mesh, *cell);
    checkH1<4>(mesh, *cell);
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, WeightedComplexQuadratureTest,
    ::testing::Values(WeightedGeometry{"Segment", Polytope::Type::Segment, 1},
      WeightedGeometry{"Triangle", Polytope::Type::Triangle, 2},
      WeightedGeometry{"Quadrilateral", Polytope::Type::Quadrilateral, 2},
      WeightedGeometry{"Tetrahedron", Polytope::Type::Tetrahedron, 3},
      WeightedGeometry{"Pyramid", Polytope::Type::Pyramid, 3},
      WeightedGeometry{"Hexahedron", Polytope::Type::Hexahedron, 3},
      WeightedGeometry{"Wedge", Polytope::Type::Wedge, 3}),
    [](const auto& info) { return std::string(info.param.name); });
}
