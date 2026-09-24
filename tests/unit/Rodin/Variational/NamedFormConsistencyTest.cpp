/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file NamedFormConsistencyTest.cpp
 * @brief Checks every named form against the Integral or integrator path that
 * computes the same operator, over cells, attributes and boundary regions.
 */
#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <type_traits>

#include "Rodin/Assembly.h"
#include "Rodin/Solid.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  namespace
  {
    constexpr Real tolerance = 1e-12;

    /// @brief Attribute of the cells left of x = 1/2.
    constexpr Attribute Left = 1;

    /// @brief Attribute of the cells right of x = 1/2.
    constexpr Attribute Right = 2;

    /// @brief Attribute of the boundary faces on y = 0.
    constexpr Attribute Bottom = 3;

    /// @brief Attribute of the other boundary faces.
    constexpr Attribute Rest = 4;

    /// @brief Uniform grid of one geometry on the unit hypercube, with the
    /// connectivity the forms and the boundary iteration need.
    LocalMesh geometryMesh(Polytope::Type geometry)
    {
      LocalMesh mesh;
      switch (geometry)
      {
        case Polytope::Type::Segment:
          mesh = LocalMesh::UniformGrid(geometry, {3});
          break;
        case Polytope::Type::Triangle:
        case Polytope::Type::Quadrilateral:
          mesh = LocalMesh::UniformGrid(geometry, {3, 3});
          break;
        case Polytope::Type::Tetrahedron:
        case Polytope::Type::Hexahedron:
        case Polytope::Type::Pyramid:
        case Polytope::Type::Wedge:
          mesh = LocalMesh::UniformGrid(geometry, {2, 2, 2});
          break;
        default:
          assert(false);
          return {};
      }
      const size_t dimension = mesh.getDimension();
      for (size_t d = dimension; d > 0; --d)
        mesh.getConnectivity().compute(d, d - 1);
      if (dimension > 1)
        mesh.getConnectivity().compute(dimension - 1, dimension);
      return mesh;
    }

    /// @brief Mean of the vertices of a polytope.
    template <class PolytopeType>
    Math::SpatialPoint centroid(const PolytopeType& polytope)
    {
      Math::SpatialPoint sum = Math::SpatialPoint::Zero(2);
      size_t count = 0;
      for (auto it = polytope.getVertex(); it; ++it, ++count)
        sum += it->getCoordinates();
      return sum / static_cast<Real>(count);
    }

    /// @brief Triangulated unit square with its cells split into Left and
    /// Right and its boundary into Bottom and Rest.
    LocalMesh attributedSquare()
    {
      LocalMesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {5, 5});
      mesh.scale(Real(1) / 4);
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      mesh.getConnectivity().compute(1, 2);
      for (auto it = mesh.getCell(); it; ++it)
      {
        mesh.setAttribute({it->getDimension(), it->getIndex()},
          centroid(*it)(0) < Real(0.5) ? Left : Right);
      }
      for (auto it = mesh.getBoundary(); it; ++it)
      {
        mesh.setAttribute({it->getDimension(), it->getIndex()},
          std::abs(centroid(*it)(1)) < Real(1e-12) ? Bottom : Rest);
      }
      return mesh;
    }

    /// @brief A coefficient whose order is unknown to the form language, so
    /// that both paths fall back to their default quadrature order, and
    /// which no quadrature of that order integrates exactly.
    RealFunction<std::function<Real(const Point&)>> smooth(Real shift = 1)
    {
      return RealFunction<std::function<Real(const Point&)>>(
        [shift](const Point& p) { return shift + std::exp(p.x()); });
    }

    /// @brief Expects two operators to agree entry by entry.
    template <class Actual, class Expected>
    void expectNear(
      const Actual& actual, const Expected& expected, const std::string& where)
    {
      const Math::Matrix<Real> a = actual;
      const Math::Matrix<Real> e = expected;
      ASSERT_EQ(a.rows(), e.rows()) << where;
      ASSERT_EQ(a.cols(), e.cols()) << where;
      Real largest = 0;
      for (Eigen::Index i = 0; i < a.rows(); ++i)
        for (Eigen::Index j = 0; j < a.cols(); ++j)
          largest = std::max(largest, std::abs(a(i, j) - e(i, j)));
      EXPECT_LE(largest, tolerance)
        << where << ": largest entry-wise difference " << largest;
      EXPECT_GT(e.cwiseAbs().maxCoeff(), 0)
        << where << ": the reference operator is zero";
    }

    /// @brief Name of a geometry for failure messages.
    std::string name(Polytope::Type geometry)
    {
      return "geometry " + std::to_string(static_cast<int>(geometry));
    }
  }

  /// @brief Parameterized fixture running the coefficient forms over one cell geometry.
  class Rodin_Variational_NamedFormConsistency
    : public ::testing::TestWithParam<Polytope::Type>
  {};

  /// @brief Verifies the weighted mass form matches Integral(c * u, v) for constant, callable and grid function coefficients.
  TEST_P(Rodin_Variational_NamedFormConsistency, WeightedMass)
  {
    auto mesh = geometryMesh(GetParam());
    const auto where = name(GetParam());
    P1 scalar(mesh);
    GridFunction gf(scalar);
    gf.project(smooth());
    const auto f = smooth();

    auto check = [&](auto& fes, const std::string& space) {
      TrialFunction u(fes);
      TestFunction v(fes);
      {
        MassForm actual(Real(2.5), u, v);
        BilinearForm expected(u, v);
        expected = Integral(Real(2.5) * u, v);
        expected.assemble();
        expectNear(actual.getOperator(), expected.getOperator(),
          where + " " + space + " constant");
      }
      {
        MassForm actual(f, u, v);
        BilinearForm expected(u, v);
        expected = Integral(f * u, v);
        expected.assemble();
        expectNear(actual.getOperator(), expected.getOperator(),
          where + " " + space + " callable");
      }
      {
        MassForm actual(gf, u, v);
        BilinearForm expected(u, v);
        expected = Integral(gf * u, v);
        expected.assemble();
        expectNear(actual.getOperator(), expected.getOperator(),
          where + " " + space + " grid function");
      }
    };

    check(scalar, "P1");
    H1 quadratic(std::integral_constant<size_t, 2>{}, mesh);
    check(quadratic, "H1P2");
  }

  /// @brief Verifies the unweighted and weighted diffusion forms match Integral(Grad(u), Grad(v)) and Integral(a * Grad(u), Grad(v)).
  TEST_P(Rodin_Variational_NamedFormConsistency, Diffusion)
  {
    auto mesh = geometryMesh(GetParam());
    const auto where = name(GetParam());
    P1 scalar(mesh);
    GridFunction gf(scalar);
    gf.project(smooth());
    const auto f = smooth();

    auto check = [&](auto& fes, const std::string& space) {
      TrialFunction u(fes);
      TestFunction v(fes);
      {
        DiffusionForm actual(u, v);
        BilinearForm expected(u, v);
        expected = Integral(Grad(u), Grad(v));
        expected.assemble();
        expectNear(actual.getOperator(), expected.getOperator(),
          where + " " + space + " unweighted");
      }
      {
        DiffusionForm actual(Real(0.75), u, v);
        BilinearForm expected(u, v);
        expected = Integral(Real(0.75) * Grad(u), Grad(v));
        expected.assemble();
        expectNear(actual.getOperator(), expected.getOperator(),
          where + " " + space + " constant");
      }
      {
        DiffusionForm actual(f, u, v);
        BilinearForm expected(u, v);
        expected = Integral(f * Grad(u), Grad(v));
        expected.assemble();
        expectNear(actual.getOperator(), expected.getOperator(),
          where + " " + space + " callable");
      }
      {
        DiffusionForm actual(gf, u, v);
        BilinearForm expected(u, v);
        expected = Integral(gf * Grad(u), Grad(v));
        expected.assemble();
        expectNear(actual.getOperator(), expected.getOperator(),
          where + " " + space + " grid function");
      }
    };

    check(scalar, "P1");
    H1 quadratic(std::integral_constant<size_t, 2>{}, mesh);
    check(quadratic, "H1P2");
  }

  /// @brief Verifies the Helmholtz form matches Integral(a * Grad(u), Grad(v)) + Integral(c * u, v) for mixed kinds of coefficients.
  TEST_P(Rodin_Variational_NamedFormConsistency, Helmholtz)
  {
    auto mesh = geometryMesh(GetParam());
    const auto where = name(GetParam());
    P1 scalar(mesh);
    GridFunction gf(scalar);
    gf.project(smooth());
    const auto f = smooth();

    auto check = [&](auto& fes, const std::string& space) {
      TrialFunction u(fes);
      TestFunction v(fes);
      {
        HelmholtzForm actual(Real(0.75), Real(-4), u, v);
        BilinearForm expected(u, v);
        expected = Integral(Real(0.75) * Grad(u), Grad(v)) + Integral(Real(-4) * u, v);
        expected.assemble();
        expectNear(actual.getOperator(), expected.getOperator(),
          where + " " + space + " constants");
      }
      {
        HelmholtzForm actual(gf, f, u, v);
        BilinearForm expected(u, v);
        expected = Integral(gf * Grad(u), Grad(v)) + Integral(f * u, v);
        expected.assemble();
        expectNear(actual.getOperator(), expected.getOperator(),
          where + " " + space + " functions");
      }
      {
        HelmholtzForm actual(f, Real(2), u, v);
        BilinearForm expected(u, v);
        expected = Integral(f * Grad(u), Grad(v)) + Integral(Real(2) * u, v);
        expected.assemble();
        expectNear(
          actual.getOperator(), expected.getOperator(), where + " " + space + " mixed");
      }
    };

    check(scalar, "P1");
    H1 quadratic(std::integral_constant<size_t, 2>{}, mesh);
    check(quadratic, "H1P2");
  }

  /// @brief Verifies the linear elasticity form matches LinearElasticityIntegral(u, v)(lambda, mu).
  TEST_P(Rodin_Variational_NamedFormConsistency, LinearElasticity)
  {
    if (GetParam() == Polytope::Type::Segment)
      GTEST_SKIP() << "linear elasticity needs a mesh of dimension two or three";

    auto mesh = geometryMesh(GetParam());
    const auto where = name(GetParam());
    P1 scalar(mesh);
    GridFunction gf(scalar);
    gf.project(smooth());
    const auto f = smooth();

    P1 fes(mesh, mesh.getSpaceDimension());
    TrialFunction u(fes);
    TestFunction v(fes);
    {
      LinearElasticityForm actual(Real(1.5), Real(0.5), u, v);
      BilinearForm expected(u, v);
      expected = LinearElasticityIntegral(u, v)(Real(1.5), Real(0.5));
      expected.assemble();
      expectNear(actual.getOperator(), expected.getOperator(), where + " constants");
    }
    {
      GridFunction gfMu(scalar);
      gfMu.project(smooth(3));
      LinearElasticityForm actual(gf, gfMu, u, v);
      BilinearForm expected(u, v);
      expected = LinearElasticityIntegral(u, v)(gf, gfMu);
      expected.assemble();
      expectNear(actual.getOperator(), expected.getOperator(), where + " grid functions");
    }
    {
      const auto mu = smooth(3);
      LinearElasticityForm actual(f, mu, u, v);
      BilinearForm expected(u, v);
      expected = LinearElasticityIntegral(u, v)(f, mu);
      expected.assemble();
      expectNear(actual.getOperator(), expected.getOperator(), where + " callables");
    }
    // Lambda and mu of different types, and a function next to a constant:
    // each deduces its own coefficient slot of the integrator.
    {
      LinearElasticityForm actual(gf, f, u, v);
      BilinearForm expected(u, v);
      expected = LinearElasticityIntegral(u, v)(gf, f);
      expected.assemble();
      expectNear(actual.getOperator(), expected.getOperator(),
        where + " grid function and callable");
    }
    {
      LinearElasticityForm actual(f, gf, u, v);
      BilinearForm expected(u, v);
      expected = LinearElasticityIntegral(u, v)(f, gf);
      expected.assemble();
      expectNear(actual.getOperator(), expected.getOperator(),
        where + " callable and grid function");
    }
    {
      LinearElasticityForm actual(gf, Real(0.5), u, v);
      BilinearForm expected(u, v);
      expected = LinearElasticityIntegral(u, v)(gf, Real(0.5));
      expected.assemble();
      expectNear(actual.getOperator(), expected.getOperator(),
        where + " grid function and constant");
    }
    {
      LinearElasticityForm actual(Real(1.5), f, u, v);
      BilinearForm expected(u, v);
      expected = LinearElasticityIntegral(u, v)(Real(1.5), f);
      expected.assemble();
      expectNear(
        actual.getOperator(), expected.getOperator(), where + " constant and callable");
    }
  }

  /// @brief Instantiates the coefficient consistency checks over every cell geometry.
  INSTANTIATE_TEST_SUITE_P(AllCellGeometries, Rodin_Variational_NamedFormConsistency,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge));

  /// @brief Verifies the weighted mass form matches Integral(c * u, v) on a vector-valued space.
  TEST(Rodin_Variational_NamedFormConsistency, WeightedVectorMass)
  {
    auto mesh = attributedSquare();
    P1 scalar(mesh);
    GridFunction gf(scalar);
    gf.project(smooth());
    P1 fes(mesh, 2);
    TrialFunction u(fes);
    TestFunction v(fes);
    MassForm actual(gf, u, v);
    BilinearForm expected(u, v);
    expected = Integral(gf * u, v);
    expected.assemble();
    expectNear(actual.getOperator(), expected.getOperator(), "vector P1");
  }

  /// @brief Verifies the coefficient forms match their integrals when trial and test spaces differ in order.
  TEST(Rodin_Variational_NamedFormConsistency, MixedOrders)
  {
    auto mesh = attributedSquare();
    const auto f = smooth();
    H1 trialFES(std::integral_constant<size_t, 2>{}, mesh);
    H1 testFES(std::integral_constant<size_t, 1>{}, mesh);
    TrialFunction u(trialFES);
    TestFunction v(testFES);
    {
      MassForm actual(f, u, v);
      BilinearForm expected(u, v);
      expected = Integral(f * u, v);
      expected.assemble();
      expectNear(actual.getOperator(), expected.getOperator(), "mass");
    }
    {
      DiffusionForm actual(f, u, v);
      BilinearForm expected(u, v);
      expected = Integral(f * Grad(u), Grad(v));
      expected.assemble();
      expectNear(actual.getOperator(), expected.getOperator(), "diffusion");
    }
    {
      HelmholtzForm actual(f, f, u, v);
      BilinearForm expected(u, v);
      expected = Integral(f * Grad(u), Grad(v)) + Integral(f * u, v);
      expected.assemble();
      expectNear(actual.getOperator(), expected.getOperator(), "helmholtz");
    }
  }

  /// @brief Verifies every form restricted to a cell attribute matches its integral restricted to the same attribute, and follows a change of attribute.
  TEST(Rodin_Variational_NamedFormConsistency, CellAttributes)
  {
    auto mesh = attributedSquare();
    const auto f = smooth();
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);

    {
      MassForm actual(f, u, v);
      for (Attribute attribute : {Left, Right})
      {
        actual.over(attribute);
        BilinearForm expected(u, v);
        expected = Integral(f * u, v).over(attribute);
        expected.assemble();
        expectNear(actual.getOperator(), expected.getOperator(),
          "mass over " + std::to_string(attribute));
      }
    }
    {
      DiffusionForm actual(f, u, v);
      for (Attribute attribute : {Left, Right})
      {
        actual.over(attribute);
        BilinearForm expected(u, v);
        expected = Integral(f * Grad(u), Grad(v)).over(attribute);
        expected.assemble();
        expectNear(actual.getOperator(), expected.getOperator(),
          "diffusion over " + std::to_string(attribute));
      }
    }
    {
      HelmholtzForm actual(f, Real(3), u, v);
      for (Attribute attribute : {Left, Right})
      {
        actual.over(attribute);
        BilinearForm expected(u, v);
        expected = Integral(f * Grad(u), Grad(v)).over(attribute) +
          Integral(Real(3) * u, v).over(attribute);
        expected.assemble();
        expectNear(actual.getOperator(), expected.getOperator(),
          "helmholtz over " + std::to_string(attribute));
      }
    }
    {
      P1 vfes(mesh, 2);
      TrialFunction w(vfes);
      TestFunction z(vfes);
      const auto mu = smooth(3);
      LinearElasticityForm actual(f, mu, w, z);
      for (Attribute attribute : {Left, Right})
      {
        actual.over(attribute);
        BilinearForm expected(w, z);
        expected = LinearElasticityIntegral(w, z)(f, mu).over(attribute);
        expected.assemble();
        expectNear(actual.getOperator(), expected.getOperator(),
          "elasticity over " + std::to_string(attribute));
      }
    }
    {
      MassForm actual(f, u, v);
      actual.over(Left, Right);
      MassForm unrestricted(f, u, v);
      expectNear(
        actual.getOperator(), unrestricted.getOperator(), "mass over both halves");
    }
  }

  /// @brief Verifies a scatter map over the boundary region matches BoundaryIntegral, with and without attributes, across reassemblies and both backends.
  TEST(Rodin_Variational_NamedFormConsistency, BoundaryRegion)
  {
    auto mesh = attributedSquare();
    const auto f = smooth();
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);

    // The kernels read the dimension off the polytope, so a cell form's
    // kernel integrates a face just as well; the region is the iteration's.
    const MassForm unweighted(u, v);
    const MassForm weighted(f, u, v);
    const FlatSet<Attribute> all;
    const FlatSet<Attribute> bottom{Bottom};

    BilinearForm expectedAll(u, v);
    expectedAll = BoundaryIntegral(u, v);
    expectedAll.assemble();

    BilinearForm expectedBottom(u, v);
    expectedBottom = BoundaryIntegral(f * u, v).over(Bottom);
    expectedBottom.assemble();

    Assembly::SequentialIteration seq(mesh, Region::Boundary);
    {
      Assembly::ScatterMap<Real> map;
      Math::SparseMatrix<Real> op;
      for (size_t pass = 0; pass < 2; ++pass)
      {
        map.assemble(op, unweighted.getKernel(), fes, fes, seq, all);
        expectNear(op, expectedAll.getOperator(), "sequential, whole boundary");
      }
      // Narrowing the selection and widening it back must rebuild both ways.
      for (size_t pass = 0; pass < 2; ++pass)
      {
        map.assemble(op, weighted.getKernel(), fes, fes, seq, bottom);
        expectNear(op, expectedBottom.getOperator(), "sequential, bottom");
      }
      map.assemble(op, unweighted.getKernel(), fes, fes, seq, all);
      expectNear(op, expectedAll.getOperator(), "sequential, whole boundary again");
    }

#ifdef RODIN_USE_OPENMP
    Assembly::OpenMPIteration omp(mesh, Region::Boundary);
    {
      Assembly::ScatterMap<Real> map;
      Math::SparseMatrix<Real> op;
      for (size_t pass = 0; pass < 2; ++pass)
      {
        map.assemble(op, unweighted.getKernel(), fes, fes, omp, all, 2);
        expectNear(op, expectedAll.getOperator(), "OpenMP, whole boundary");
      }
      for (size_t pass = 0; pass < 2; ++pass)
      {
        map.assemble(op, weighted.getKernel(), fes, fes, omp, bottom, 2);
        expectNear(op, expectedBottom.getOperator(), "OpenMP, bottom");
      }
      map.assemble(op, unweighted.getKernel(), fes, fes, omp, all, 2);
      expectNear(op, expectedAll.getOperator(), "OpenMP, whole boundary again");
    }
#endif
  }

  /// @brief Verifies the named forms agree with their integrals on reassembly through both backends.
  TEST(Rodin_Variational_NamedFormConsistency, ReassemblyAcrossBackends)
  {
    auto mesh = attributedSquare();
    const auto f = smooth();
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    HelmholtzForm form(f, Real(-2), u, v);
    form.over(Right);

    BilinearForm expected(u, v);
    expected =
      Integral(f * Grad(u), Grad(v)).over(Right) + Integral(Real(-2) * u, v).over(Right);
    expected.assemble();

    using Form = decltype(form);
    typename Form::OperatorType sequentialOperator;
    Assembly::Sequential<typename Form::OperatorType, Form> sequential;
    for (size_t pass = 0; pass < 3; ++pass)
    {
      sequential.execute(sequentialOperator, form);
      expectNear(sequentialOperator, expected.getOperator(), "sequential");
    }

#ifdef RODIN_USE_OPENMP
    typename Form::OperatorType openMPOperator;
    Assembly::OpenMP<typename Form::OperatorType, Form> openMP;
    for (size_t pass = 0; pass < 3; ++pass)
    {
      openMP.execute(openMPOperator, form);
      expectNear(openMPOperator, expected.getOperator(), "OpenMP");
    }
#endif
  }
}
