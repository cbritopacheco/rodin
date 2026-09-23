/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief L2 projection convergence of real and complex scalar/vector P0.
 *
 * A smooth affine field is projected through the mass problem
 * @f$(u_h,v_h)=(f,v_h)@f$. Its cell means converge in L2 at order one.
 * Since this uses the assembled mass matrix, it also validates the complex
 * sesquilinear and vector component paths of the discontinuous space.
 */

#include <cmath>
#include <cstdint>
#include <functional>
#include <string>
#include <type_traits>

#include <gtest/gtest.h>

#include "Convergence.h"
#include "Rodin/Assembly.h"
#include "Rodin/Solver/CG.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Solver;
using namespace Rodin::Variational;

namespace Rodin::Tests::Convergence::H::P0Projection
{
  template <class Scalar>
  Scalar coefficient(Real real, Real imaginary)
  {
    if constexpr (std::is_same_v<Scalar, Real>)
      return real;
    else
      return Scalar(real, imaginary);
  }

  template <class FES, class Exact>
  Real project(const LocalMesh& mesh, FES& space, const Exact& exact)
  {
    TrialFunction u(space);
    TestFunction v(space);
    auto mass = Integral(u, v);
    auto load = Integral(exact, v);
    mass.setOrder(10);
    load.setOrder(10);
    Problem projection(u, v);
    projection = mass - load;
    CG solver(projection);
    solver.setTolerance(1e-13).setMaxIterations(20000).solve();
    EXPECT_TRUE(solver.success());
    return ErrorNorm::computeL2(mesh, u.getSolution(), exact, 12);
  }

  template <class Scalar>
  Real scalarError(const LocalMesh& mesh)
  {
    P0<Scalar, LocalMesh> space(mesh);
    const auto field = [dim = mesh.getDimension()](const Point& point)
    {
      return coefficient<Scalar>(1.25, 0.5)
        + coefficient<Scalar>(1, -0.75) * point(0)
        + coefficient<Scalar>(0.3, 0.2) * point(dim - 1);
    };
    if constexpr (std::is_same_v<Scalar, Real>)
      return project(mesh, space, RealFunction(field));
    else
      return project(mesh, space, ComplexFunction(field));
  }

  template <class Scalar>
  Real vectorError(const LocalMesh& mesh)
  {
    P0<Math::SpatialVector<Scalar>, LocalMesh> space(mesh, 2);
    const auto field = [dim = mesh.getDimension()](const Point& point)
    {
      Math::SpatialVector<Scalar> value(2);
      value(0) = coefficient<Scalar>(1.25, 0.5)
        + coefficient<Scalar>(1, -0.75) * point(0);
      value(1) = coefficient<Scalar>(-0.75, 0.25)
        + coefficient<Scalar>(0.8, 0.2)
        * point(dim - 1);
      return value;
    };
    return project(mesh, space, VectorFunction<decltype(field)>(2, field));
  }

  template <class Error>
  void checkRates(Polytope::Type geometry, Error error)
  {
    UniformGridHierarchy hierarchy(geometry, {3, 5, 9});
    NormHistory history;
    for (size_t level : hierarchy.getLevels())
    {
      auto mesh = hierarchy.makeMesh(level);
      history.append(hierarchy.getMeshSize(level), error(mesh));
    }
    ASSERT_EQ(history.getSize(), 3);
    for (size_t i = 1; i < history.getSize(); ++i)
    {
      const Real coarse = history.getSample(i - 1).error;
      const Real fine = history.getSample(i).error;
      ASSERT_TRUE(std::isfinite(coarse));
      ASSERT_TRUE(std::isfinite(fine));
      ASSERT_GT(coarse, fine);
      const Real rate = history.getAlgebraicRate(i);
      SCOPED_TRACE(::testing::Message()
        << "L2 " << coarse << " -> " << fine << ", rate " << rate);
      EXPECT_GT(rate, 0.8);
      EXPECT_LT(rate, 1.2);
    }
  }

  class P0ProjectionTest : public ::testing::TestWithParam<Polytope::Type> {};

  TEST_P(P0ProjectionTest, RealScalar)
  {
    checkRates(GetParam(), scalarError<Real>);
  }

  TEST_P(P0ProjectionTest, ComplexScalar)
  {
    checkRates(GetParam(), scalarError<Complex>);
  }

  TEST_P(P0ProjectionTest, RealVector)
  {
    checkRates(GetParam(), vectorError<Real>);
  }

  TEST_P(P0ProjectionTest, ComplexVector)
  {
    checkRates(GetParam(), vectorError<Complex>);
  }

  /** P0g contains constants exactly, so its assembled projection is exact. */
  class P0gProjectionTest : public ::testing::TestWithParam<Polytope::Type> {};

  template <class Scalar>
  Real globalScalarError(const LocalMesh& mesh)
  {
    P0g<Scalar, LocalMesh> space(mesh);
    const auto value = coefficient<Scalar>(1.25, -0.75);
    const auto field = [value](const Point&) { return value; };
    if constexpr (std::is_same_v<Scalar, Real>)
      return project(mesh, space, RealFunction(field));
    else
      return project(mesh, space, ComplexFunction(field));
  }

  template <class Scalar>
  Real globalVectorError(const LocalMesh& mesh)
  {
    P0g<Math::SpatialVector<Scalar>, LocalMesh> space(mesh, 2);
    const auto field = [](const Point&)
    {
      Math::SpatialVector<Scalar> value(2);
      value(0) = coefficient<Scalar>(1.25, -0.75);
      value(1) = coefficient<Scalar>(-0.5, 0.25);
      return value;
    };
    return project(mesh, space, VectorFunction<decltype(field)>(2, field));
  }

  TEST_P(P0gProjectionTest, ExactGlobalConstants)
  {
    UniformGrid grid(GetParam());
    const auto mesh = grid.makeMesh(4);
    EXPECT_LT(globalScalarError<Real>(mesh), 1e-10);
    EXPECT_LT(globalScalarError<Complex>(mesh), 1e-10);
    EXPECT_LT(globalVectorError<Real>(mesh), 1e-10);
    EXPECT_LT(globalVectorError<Complex>(mesh), 1e-10);
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, P0ProjectionTest,
    ::testing::Values(
      Polytope::Type::Segment,
      Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron,
      Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron,
      Polytope::Type::Wedge),
    [](const ::testing::TestParamInfo<Polytope::Type>& info)
    {
      return std::string(UniformGrid::getGeometryName(info.param));
    });

  INSTANTIATE_TEST_SUITE_P(AllGeometries, P0gProjectionTest,
    ::testing::Values(
      Polytope::Type::Segment,
      Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron,
      Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron,
      Polytope::Type::Wedge),
    [](const ::testing::TestParamInfo<Polytope::Type>& info)
    {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
}
