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
    /** Observes retention requests without altering mesh-cache lifetimes. */
    class ObservedMesh : public LocalMesh
    {
      public:
        explicit ObservedMesh(LocalMesh&& mesh)
          : LocalMesh(std::move(mesh)),
            requests(0)
        {}

        const PolytopeQuadrature& getQuadrature(
          size_t d, Index i, const QF::QuadratureFormulaBase& qf) const override
        {
          ++requests;
          return LocalMesh::getQuadrature(d, i, qf);
        }

        mutable size_t requests;
    };

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
      const auto& mesh = static_cast<const ObservedMesh&>(cell.getMesh());
      const size_t requests = mesh.requests;
      integral.setPolytope(cell);
      EXPECT_EQ(mesh.requests, requests)
        << "Local assembly must not retain mapped points for every mesh cell";
      const auto& formula = QF::PolytopeQuadratureFormula::get(order, cell.getGeometry());
      for (const auto* ip : {&integral.getIntegrand().getLHS().getIntegrationPoint(),
             &integral.getIntegrand().getRHS().getIntegrationPoint()})
      {
        // The binding remains readable after setPolytope returns. ASan must
        // reject a shape expression retaining a stack-local integration point.
        EXPECT_EQ(ip->getQuadratureFormula(), &formula);
        EXPECT_EQ(ip->getIndex(), formula.getSize() - 1);
        EXPECT_EQ(ip->getPoint().getPolytope().getDimension(), cell.getDimension());
        EXPECT_EQ(ip->getPoint().getPolytope().getIndex(), cell.getIndex());
        EXPECT_EQ(&ip->getPoint().getPolytope().getMesh(), &cell.getMesh());
      }
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
        const auto& observed = static_cast<const ObservedMesh&>(mesh);
        const size_t requests = observed.requests;
        form.assemble();
        EXPECT_EQ(observed.requests, requests);
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
          // The sum selects generic quadrature: specialised bare H1/P1
          // kernels tabulate directly and do not bind the expression tree.
          auto linear = Integral(v + v);
          linear.setOrder(3);
          linear.setPolytope(*cell);
          const auto* bound = &linear.getIntegrand().getIntegrationPoint();
          const auto& formula =
            QF::PolytopeQuadratureFormula::get(3, cell->getGeometry());
          EXPECT_EQ(bound->getQuadratureFormula(), &formula);
          EXPECT_EQ(bound->getIndex(), formula.getSize() - 1);
          EXPECT_EQ(bound->getPoint().getPolytope().getDimension(), cell->getDimension());
          EXPECT_EQ(bound->getPoint().getPolytope().getIndex(), cell->getIndex());
          EXPECT_EQ(&bound->getPoint().getPolytope().getMesh(), &cell->getMesh());
          auto movedLinear = std::move(linear);
          EXPECT_EQ(&movedLinear.getIntegrand().getIntegrationPoint(), bound);
          auto copiedLinear = movedLinear;
          copiedLinear.setPolytope(*cell);
          EXPECT_NE(&copiedLinear.getIntegrand().getIntegrationPoint(), bound);
          EXPECT_EQ(copiedLinear.getIntegrand().getIntegrationPoint().getIndex(),
            formula.getSize() - 1);
        }
        auto copy = mass;
        checkLocal(copy, *cell, 3);
        const auto* boundTrial = &copy.getIntegrand().getLHS().getIntegrationPoint();
        const auto* boundTest = &copy.getIntegrand().getRHS().getIntegrationPoint();
        auto moved = std::move(copy);
        EXPECT_EQ(&moved.getIntegrand().getLHS().getIntegrationPoint(), boundTrial);
        EXPECT_EQ(&moved.getIntegrand().getRHS().getIntegrationPoint(), boundTest);
        checkLocal(moved, *cell);
      }
      checkAssembly(space, u, v, mass);
    }

    template <size_t K, class Scalar>
    void checkH1(Polytope::Type geometry)
    {
      ObservedMesh mesh(Convergence::UniformGrid(geometry).makeMesh(2));
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

    template <class F>
    void checkFunctionLifetime(const FunctionBase<F>& function, const ObservedMesh& mesh)
    {
      using Rule = QuadratureRule<FunctionBase<F>>;
      using Scalar = typename Rule::ScalarType;
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        const auto& qf = QF::PolytopeQuadratureFormula::get(1, cell->getGeometry());
        const auto& borrowed = cell->getQuadrature(qf);
        Scalar expected = 0;
        for (size_t qp = 0; qp < borrowed.getSize(); ++qp)
        {
          const auto& p = borrowed.getPoint(qp);
          expected += qf.getWeight(qp) * p.getDistortion() * function(p);
        }

        const size_t requests = mesh.requests;
        Rule rule(function);
        rule.setPolytope(*cell);
        EXPECT_EQ(rule.compute(), expected);
        Rule copy(rule);
        copy.setPolytope(*cell);
        EXPECT_EQ(copy.compute(), expected);
        Rule moved(std::move(rule));
        EXPECT_EQ(moved.compute(), expected);
        EXPECT_EQ(mesh.requests, requests);
        // Internal ownership must not evict a publicly borrowed quadrature.
        EXPECT_EQ(&cell->getQuadrature(qf), &borrowed);
      }
    }

    TEST_P(QuadratureExpression, FunctionOwnershipCopyMoveAndBorrowedCacheLifetime)
    {
      ObservedMesh mesh(Convergence::UniformGrid(GetParam()).makeMesh(3));
      const RealFunction real([](const Point& p) { return 1 + p(0); });
      const ComplexFunction complex(
        [](const Point& p) { return Complex(1 + p(0), 2 - p(0)); });
      checkFunctionLifetime(real, mesh);
      checkFunctionLifetime(complex, mesh);
    }

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
      ObservedMesh mesh(Convergence::UniformGrid(geometry).makeMesh(2));
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

    TEST_P(QuadratureExpression, LinearBindingAllRangesAndH1Orders)
    {
      ObservedMesh mesh(Convergence::UniformGrid(GetParam()).makeMesh(2));
      const auto check = [&]<class F>(const F& space) {
        TestFunction v(space);
        // A component is scalar-valued, but still binds the complete parent
        // shape expression. This avoids a quadratic-size high-order matrix.
        auto integrand = [&] {
          using Range = typename FormLanguage::Traits<F>::RangeType;
          if constexpr (FormLanguage::IsMatrixRange<Range>::Value)
            return Component(v + v, 1, 2);
          else if constexpr (FormLanguage::IsVectorRange<Range>::Value)
            return Component(v + v, 2);
          else
            return v + v;
        }();
        auto integral = Integral(integrand);
        const IntegrationPoint* retained = nullptr;
        for (auto cell = mesh.getCell(); cell; ++cell)
        {
          integral.setOrder(1);
          integral.setPolytope(*cell);
          const auto* bound = &integral.getIntegrand().getIntegrationPoint();
          if (retained)
          {
            EXPECT_EQ(bound, retained);
          }
          retained = bound;
          const auto verify = [&](const auto& rule, size_t order) {
            const auto& ip = rule.getIntegrand().getIntegrationPoint();
            const auto& formula =
              QF::PolytopeQuadratureFormula::get(order, cell->getGeometry());
            EXPECT_EQ(ip.getQuadratureFormula(), &formula);
            EXPECT_EQ(ip.getIndex(), formula.getSize() - 1);
            EXPECT_EQ(ip.getPoint().getPolytope().getDimension(), cell->getDimension());
            EXPECT_EQ(ip.getPoint().getPolytope().getIndex(), cell->getIndex());
            EXPECT_EQ(&ip.getPoint().getPolytope().getMesh(), &mesh);
          };
          verify(integral, 1);
          auto copy = integral;
          copy.setOrder(2);
          copy.setPolytope(*cell);
          const auto* independent = &copy.getIntegrand().getIntegrationPoint();
          EXPECT_NE(independent, bound);
          verify(copy, 2);
          // Rebinding the copy cannot mutate the original borrowed context.
          verify(integral, 1);
          auto moved = std::move(copy);
          EXPECT_EQ(&moved.getIntegrand().getIntegrationPoint(), independent);
          verify(moved, 2);
          moved.setOrder(1);
          moved.setPolytope(*cell);
          EXPECT_EQ(&moved.getIntegrand().getIntegrationPoint(), independent);
          verify(moved, 1);
        }
      };
      const auto ranges = [&]<class Scalar>() {
        using Vector = Math::SpatialVector<Scalar>;
        using Matrix = Math::SpatialMatrix<Scalar>;
        check(P0<Matrix, LocalMesh>(mesh, 2, 3));
        check(P0g<Matrix, LocalMesh>(mesh, 2, 3));
        check(P1<Matrix, LocalMesh>(mesh, 2, 3));
        const auto order = [&]<size_t K>(std::integral_constant<size_t, K> degree) {
          SCOPED_TRACE(::testing::Message() << "degree=" << K);
          check(H1<K, Scalar, LocalMesh>(degree, mesh));
          check(H1<K, Vector, LocalMesh>(degree, mesh, 3));
          check(H1<K, Matrix, LocalMesh>(degree, mesh, 2, 3));
        };
        order(std::integral_constant<size_t, 1>{});
        order(std::integral_constant<size_t, 2>{});
        order(std::integral_constant<size_t, 3>{});
        order(std::integral_constant<size_t, 4>{});
        order(std::integral_constant<size_t, 5>{});
        order(std::integral_constant<size_t, 6>{});
      };
      ranges.template operator()<Real>();
      ranges.template operator()<Complex>();
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
          const size_t nq =
            QF::PolytopeQuadratureFormula::get(6, cell->getGeometry()).getSize();
          // The function cache binds once per quadrature point, independently
          // of both trial and test counts. Rebinding still refreshes its value.
          EXPECT_EQ(calls, nq);
          const auto expected = reference(integral, *cell, 6);
          for (Eigen::Index te = 0; te < expected.rows(); ++te)
            for (Eigen::Index tr = 0; tr < expected.cols(); ++tr)
              EXPECT_EQ(integral.integrate(tr, te), expected(te, tr));
        }
      }
    }

    TEST(QuadratureExpressionRegression, CurvedCellAndFacetValues)
    {
      ObservedMesh mesh(Convergence::UniformGrid(Polytope::Type::Triangle).makeMesh(2));
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
