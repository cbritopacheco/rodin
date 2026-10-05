/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file @brief Generic bilinear quadrature value-reuse regression tests.
 * The test-only reference retains the original quadrature/trial/test loop.
 * Every local entry must agree exactly, with unchanged summation order and
 * trial-first, test-second conjugation. Assembled comparisons use DOF maps,
 * never assumed coefficient layouts. Callback counts are integer assertions.
 */
#include <gtest/gtest.h>
#include "Rodin/Assembly.h"
#include "Rodin/Variational.h"
#include "Rodin/Variational/H1.h"
#include "Rodin/Variational/P0g.h"
#include "../../../convergence/Convergence.h"
#include "../../../QuadratureReference.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  namespace
  {
    /** @brief Original entry evaluation, independent of the production rule. */
    template <class I>
    auto reference(const I& integral, const Polytope& cell, size_t order)
    {
      QuadratureReference<I> baseline(integral);
      baseline.assemble(cell, order);
      return baseline.getOperator();
    }

    template <class I>
    void checkLocal(I& integral, const Polytope& cell, size_t order = 6)
    {
      integral.setOrder(order);
      const auto expected = reference(integral, cell, order);
      integral.setPolytope(cell);
      for (Eigen::Index te = 0; te < expected.rows(); ++te)
        for (Eigen::Index tr = 0; tr < expected.cols(); ++tr)
          EXPECT_EQ(integral.integrate(tr, te), expected(te, tr))
            << "cell " << cell.getIndex() << ", entry " << te << ',' << tr;
    }

    template <class F, class U, class V, class I>
    void checkAssembly(F& space, U& u, V& v, I& integral)
    {
      using Scalar = typename FormLanguage::Traits<typename I::IntegrandType>::ScalarType;
      Math::Matrix<Scalar> expected(space.getSize(), space.getSize());
      expected.setZero();
      const auto& mesh = space.getMesh();
      const size_t d = mesh.getDimension();
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        const auto local = reference(integral, *cell, 6);
        const auto& dofs = space.getDOFs(d, cell->getIndex());
        for (Eigen::Index te = 0; te < local.rows(); ++te)
          for (Eigen::Index tr = 0; tr < local.cols(); ++tr)
            expected(dofs(te), dofs(tr)) += local(te, tr);
      }
      BilinearForm form(u, v);
      integral.setOrder(6);
      form = integral;
      for (size_t repeat = 0; repeat < 2; ++repeat)
      {
        form.assemble();
        const Math::Matrix<Scalar> actual = form.getOperator();
        EXPECT_LE((actual - expected).norm(), 1e-12 * std::max(Real(1), expected.norm()));
      }
    }

    template <class F>
    void checkSpace(F& space)
    {
      TrialFunction u(space);
      TestFunction v(space);
      RealFunction coefficient([](const Point& p) {
        Real value = 1;
        for (size_t j = 0; j < p.getDimension(); ++j)
          value += p(j);
        return value;
      });
      auto mass = Integral(coefficient * (u + u), v - 0.25 * v);
      for (auto cell = space.getMesh().getCell(); cell; ++cell)
      {
        checkLocal(mass, *cell);
        if constexpr (std::is_same_v<typename FormLanguage::Traits<F>::ScalarType,
                        Complex>)
        {
          const Complex a(1, 2), b(2, -1);
          auto phased =
            Integral((u + u) * ComplexFunction(a), (v + v) * ComplexFunction(b));
          checkLocal(phased, *cell);
          auto unscaled = Integral(u + u, v + v);
          unscaled.setOrder(6);
          unscaled.setPolytope(*cell);
          const size_t count =
            space.getDOFs(cell->getDimension(), cell->getIndex()).size();
          for (size_t te = 0; te < count; ++te)
            for (size_t tr = 0; tr < count; ++tr)
              EXPECT_LE(std::abs(phased.integrate(tr, te) -
                          a * std::conj(b) * unscaled.integrate(tr, te)),
                1e-12);
        }
        if constexpr (!FormLanguage::IsVectorRange<
                        typename FormLanguage::Traits<F>::RangeType>::Value)
        {
          auto conjugated = Integral(u + u, Conjugate(v + v));
          checkLocal(conjugated, *cell);
        }
        auto copy = mass;
        checkLocal(copy, *cell, 3);
        auto moved = std::move(copy);
        checkLocal(moved, *cell);
      }
      checkAssembly(space, u, v, mass);
    }

    template <size_t K, class Scalar>
    void checkH1(Polytope::Type geometry)
    {
      auto mesh = Convergence::UniformGrid(geometry).makeMesh(2);
      const size_t d = mesh.getDimension();
      H1<K, Scalar, LocalMesh> scalar(std::integral_constant<size_t, K>{}, mesh);
      H1<K, Math::SpatialVector<Scalar>, LocalMesh> vector(
        std::integral_constant<size_t, K>{}, mesh, d);
      checkSpace(scalar);
      checkSpace(vector);
      TrialFunction u(vector);
      TestFunction v(vector);
      auto shear = Integral(
        0.5 * (Jacobian(u) + Jacobian(u).T()), 0.5 * (Jacobian(v) + Jacobian(v).T()));
      auto divergence = Integral(Div(u), Div(v));
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        checkLocal(shear, *cell);
        checkLocal(divergence, *cell);
      }
      checkAssembly(vector, u, v, shear);
    }

    class QuadratureExpression : public ::testing::TestWithParam<Polytope::Type>
    {};

    TEST_P(QuadratureExpression, H1RealOrdersOneTwoThree)
    {
      checkH1<1, Real>(GetParam());
      checkH1<2, Real>(GetParam());
      checkH1<3, Real>(GetParam());
    }

    TEST_P(QuadratureExpression, H1ComplexOrdersOneTwoThree)
    {
      checkH1<1, Complex>(GetParam());
      checkH1<2, Complex>(GetParam());
      checkH1<3, Complex>(GetParam());
    }

    template <class Scalar>
    void checkOtherSpaces(Polytope::Type geometry)
    {
      auto mesh = Convergence::UniformGrid(geometry).makeMesh(2);
      const size_t d = mesh.getDimension();
      P1<Scalar, LocalMesh> p1(mesh);
      P1<Math::SpatialVector<Scalar>, LocalMesh> vp1(mesh, d);
      P0<Scalar, LocalMesh> p0(mesh);
      P0<Math::SpatialVector<Scalar>, LocalMesh> vp0(mesh, d);
      P0g<Scalar, LocalMesh> p0g(mesh);
      P0g<Math::SpatialVector<Scalar>, LocalMesh> vp0g(mesh, d);
      checkSpace(p1);
      checkSpace(vp1);
      checkSpace(p0);
      checkSpace(vp0);
      checkSpace(p0g);
      checkSpace(vp0g);
    }

    TEST_P(QuadratureExpression, ConstantAndP1RealComplexScalarVector)
    {
      checkOtherSpaces<Real>(GetParam());
      checkOtherSpaces<Complex>(GetParam());
    }

    INSTANTIATE_TEST_SUITE_P(AllGeometries, QuadratureExpression,
      ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
        Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
        Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge));

    /** @brief Reuse is established by integer counts, not timing thresholds. */
    TEST(QuadratureExpressionRegression, CountsAndCoefficientRebinds)
    {
      auto mesh = Convergence::UniformGrid(Polytope::Type::Triangle).makeMesh(3);
      H1<2, Math::SpatialVector<Real>, LocalMesh> trialSpace(
        std::integral_constant<size_t, 2>{}, mesh, 2);
      H1<1, Math::SpatialVector<Real>, LocalMesh> testSpace(
        std::integral_constant<size_t, 1>{}, mesh, 2);
      TrialFunction u(trialSpace);
      TestFunction v(testSpace);
      size_t calls = 0;
      Real scale = 1;
      RealFunction coefficient([&](const Point& p) {
        ++calls;
        return scale * (1 + p(0));
      });
      auto integral = Integral(
        Jacobian(u) + Jacobian(u).T(), coefficient * (Jacobian(v) + Jacobian(v).T()));
      integral.setOrder(6);
      for (Real value : {Real(1), Real(2), Real(-0.5)})
      {
        scale = value;
        for (auto cell = mesh.getCell(); cell; ++cell)
        {
          calls = 0;
          integral.setPolytope(*cell);
          const size_t nt =
            testSpace.getDOFs(mesh.getDimension(), cell->getIndex()).size();
          const size_t nq =
            QF::PolytopeQuadratureFormula::get(6, cell->getGeometry()).getSize();
          EXPECT_EQ(calls, nq * nt);
          const auto expected = reference(integral, *cell, 6);
          for (Eigen::Index te = 0; te < expected.rows(); ++te)
            for (Eigen::Index tr = 0; tr < expected.cols(); ++tr)
              EXPECT_EQ(integral.integrate(tr, te), expected(te, tr));
        }
      }
    }

    TEST(QuadratureExpressionRegression, CurvedCellAndFacetValues)
    {
      auto mesh = Convergence::UniformGrid(Polytope::Type::Triangle).makeMesh(2);
      const auto cell = mesh.getCell(0);
      RealH1Element<2> geometryFE(Polytope::Type::Triangle);
      PointCloud nodes(2, geometryFE.getCount());
      for (size_t local = 0; local < geometryFE.getCount(); ++local)
      {
        const auto& rc = geometryFE.getNode(local);
        Math::SpatialPoint x;
        cell->getTransformation().transform(x, rc);
        x(1) += 0.05 * rc(0) * (1 - rc(0) - rc(1));
        nodes(0, local) = x(0);
        nodes(1, local) = x(1);
      }
      mesh.setPolytopeTransformation({2, cell->getIndex()},
        new ParametricTransformation<RealH1Element<2>>(std::move(nodes), geometryFE));
      H1<2, Complex, LocalMesh> space(std::integral_constant<size_t, 2>{}, mesh);
      checkSpace(space);
      TrialFunction u(space);
      TestFunction v(space);
      auto boundary = BoundaryIntegral(u + u, v + v);
      auto faces = FaceIntegral(u + u, v + v);
      auto interfaces = InterfaceIntegral(u + u, v + v);
      for (auto face = mesh.getPolytope(1); face; ++face)
      {
        checkLocal(faces, *face);
        if (mesh.isBoundary(face->getIndex()))
          checkLocal(boundary, *face);
        else
          checkLocal(interfaces, *face);
      }
    }
  }
}
