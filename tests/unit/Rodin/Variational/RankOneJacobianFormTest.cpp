#include <gtest/gtest.h>

#include <array>
#include <type_traits>

#include "Rodin/Variational.h"
#include "Rodin/Variational/H1.h"
#include "Rodin/Variational/P1.h"
#include "Rodin/Assembly/Default.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  namespace
  {
    class ScaledIdentityMatrix final
      : public MatrixFunctionBase<Real, ScaledIdentityMatrix>
    {
      public:
        ScaledIdentityMatrix(size_t dimension, Real scale)
          : m_dimension(dimension),
            m_scale(scale)
        {}

        size_t getRows() const
        {
          return m_dimension;
        }

        size_t getColumns() const
        {
          return m_dimension;
        }

        Math::SpatialMatrix<Real> getValue(const Point&) const
        {
          return m_scale * Math::SpatialMatrix<Real>::Identity(m_dimension, m_dimension);
        }

        Optional<size_t> getOrder(const Polytope&) const noexcept
        {
          return 0;
        }

        ScaledIdentityMatrix* copy() const noexcept override
        {
          return new ScaledIdentityMatrix(*this);
        }

      private:
        size_t m_dimension;
        Real m_scale;
    };

    LocalMesh unitSquare(size_t n)
    {
      LocalMesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {n, n});
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      return mesh;
    }

    LocalMesh unitCube(size_t n)
    {
      LocalMesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {n, n, n});
      mesh.getConnectivity().compute(3, 2);
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      return mesh;
    }

    LocalMesh uniformGrid(Polytope::Type geometry, size_t n)
    {
      const size_t dimension = Polytope::Traits(geometry).getDimension();
      Array<size_t> shape(dimension);
      for (size_t i = 0; i < dimension; ++i)
        shape(i) = n;
      LocalMesh mesh = LocalMesh::UniformGrid(geometry, shape);
      for (size_t d = dimension; d > 0; --d)
        mesh.getConnectivity().compute(d, d - 1);
      return mesh;
    }

    template <class FES>
    void expectIdentityCoefficientEqualsDivergence(FES& fes, size_t dimension)
    {
      TrialFunction u(fes);
      TestFunction v(fes);
      const IdentityMatrix A(dimension);

      BilinearForm rankOne(u, v);
      rankOne = Integral(Dot(Dot(A, Jacobian(u)), Dot(A, Jacobian(v))));
      rankOne.assemble();

      BilinearForm divergence(u, v);
      divergence = Integral(Div(u), Div(v));
      divergence.assemble();

      EXPECT_EQ(rankOne.getOperator().rows(), divergence.getOperator().rows());
      const Math::Matrix<Real> actual = rankOne.getOperator();
      const Math::Matrix<Real> expected = divergence.getOperator();
      for (Eigen::Index i = 0; i < actual.rows(); ++i)
      {
        for (Eigen::Index j = 0; j < actual.cols(); ++j)
        {
          EXPECT_NEAR(actual(i, j), expected(i, j), 1e-12)
            << "entry (" << i << ", " << j << ")";
        }
      }
    }

    template <class FES>
    void expectDistinctCoefficientsUseBothSides(FES& fes, size_t dimension)
    {
      TrialFunction u(fes);
      TestFunction v(fes);
      const ScaledIdentityMatrix A(dimension, 1);
      const ScaledIdentityMatrix B(dimension, 2);

      BilinearForm rankOne(u, v);
      rankOne = Integral(Dot(Dot(A, Jacobian(u)), Dot(B, Jacobian(v))));
      rankOne.assemble();

      BilinearForm divergence(u, v);
      divergence = Integral(Div(u), Div(v));
      divergence.assemble();

      const Math::Matrix<Real> actual = rankOne.getOperator();
      const Math::Matrix<Real> expected = divergence.getOperator();
      ASSERT_EQ(actual.rows(), expected.rows());
      ASSERT_EQ(actual.cols(), expected.cols());
      for (Eigen::Index i = 0; i < actual.rows(); ++i)
      {
        for (Eigen::Index j = 0; j < actual.cols(); ++j)
        {
          EXPECT_NEAR(actual(i, j), Real(2) * expected(i, j), 1e-12)
            << "entry (" << i << ", " << j << ")";
        }
      }
    }
  }

  namespace
  {
    /**
     * @brief Checks all ambient columns of embedded rank-one Jacobian forms.
     *
     * Basis derivatives and the rectangular inverse Jacobian independently
     * form the physical Frobenius contractions. Axis-aligned and oblique
     * segments and a surface in three dimensions expose omitted columns.
     * Shared and distinct coefficients exercise both contraction branches.
     *
     * @tparam Scalar Real or complex field scalar.
     * @tparam HighOrder Selects degree-two H1 rather than P1.
     */
    template <class Scalar, bool HighOrder>
    void checkEmbeddedRankOne()
    {
      constexpr size_t Components = 3;
      constexpr size_t Order = 6;
      const std::array<Math::SpatialPoint, 4> endpoints = {Math::SpatialPoint{1, 0, 0},
        Math::SpatialPoint{0, 1, 0}, Math::SpatialPoint{0, 0, 1},
        Math::SpatialPoint{1, 2, 3}};
      for (size_t fixture = 0; fixture <= endpoints.size(); ++fixture)
      {
        SCOPED_TRACE(fixture);
        const bool surface = fixture == endpoints.size();
        auto builder = LocalMesh::Builder();
        builder.initialize(3).nodes(surface ? 3 : 2).vertex({0, 0, 0});
        if (surface)
          builder.vertex({0, 1, 0}).vertex({0, 0, 1});
        else
          builder.vertex(endpoints[fixture]);
        const auto geometry =
          surface ? Polytope::Type::Triangle : Polytope::Type::Segment;
        auto mesh = surface ? builder.polytope(geometry, {0, 1, 2}).finalize()
                            : builder.polytope(geometry, {0, 1}).finalize();
        if (surface)
        {
          mesh.getConnectivity().compute(2, 1);
          mesh.getConnectivity().compute(1, 0);
        }
        const auto cell = *mesh.getCell();
        auto space = [&] {
          if constexpr (HighOrder)
          {
            return H1<2, Math::SpatialVector<Scalar>>(
              std::integral_constant<size_t, 2>{}, mesh, Components);
          }
          else
            return P1<Math::SpatialVector<Scalar>>(mesh, Components);
        }();
        TrialFunction u(space);
        TestFunction v(space);
        Math::SpatialMatrix<Scalar> first(Components, 3), second(Components, 3);
        for (size_t component = 0; component < Components; ++component)
        {
          for (size_t axis = 0; axis < 3; ++axis)
          {
            first(component, axis) = Scalar(1 + component + axis);
            second(component, axis) = Scalar(2 + 2 * component - Real(axis));
            if constexpr (std::is_same_v<Scalar, Complex>)
            {
              first(component, axis) += Scalar(0, Real(component + 1) / 2);
              second(component, axis) -= Scalar(0, Real(axis + 1) / 4);
            }
          }
        }
        const auto callable = [](const Math::SpatialMatrix<Scalar>& value) {
          return [value](const Geometry::Point&) { return value; };
        };
        const auto A = MatrixFunction(Components, 3, callable(first));
        const auto B = MatrixFunction(Components, 3, callable(second));
        const auto& formula = QF::QuadratureFormula::get(Order, geometry);
        const Geometry::PolytopeQuadrature quadrature(cell, formula);
        const auto& element =
          space.getFiniteElement(cell.getDimension(), cell.getIndex());
        const auto check = [&](auto& integral, bool shared) {
          integral.setOrder(Order);
          integral.setPolytope(cell);
          Real expectedSquaredNorm = 0;
          for (size_t test = 0; test < element.getCount(); ++test)
          {
            for (size_t trial = 0; trial < element.getCount(); ++trial)
            {
              Scalar expected = 0;
              for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
              {
                const auto& point = quadrature.getPoint(qp);
                const auto& reference = point.getReferenceCoordinates();
                const auto trialReference =
                  element.getBasis(trial).getJacobian()(reference);
                const auto testReference =
                  element.getBasis(test).getJacobian()(reference);
                const auto inverse = point.getJacobianInverse();
                Scalar trialAction = 0, testAction = 0;
                for (size_t component = 0; component < Components; ++component)
                {
                  for (size_t axis = 0; axis < 3; ++axis)
                  {
                    for (size_t r = 0; r < cell.getDimension(); ++r)
                    {
                      trialAction += first(component, axis) *
                        trialReference(component, r) * inverse(r, axis);
                      testAction += (shared ? first : second)(component, axis) *
                        testReference(component, r) * inverse(r, axis);
                    }
                  }
                }
                expected += formula.getWeight(qp) * point.getDistortion() * trialAction *
                  Math::conj(testAction);
              }
              expectedSquaredNorm += std::norm(expected);
              EXPECT_NEAR(std::abs(integral.integrate(trial, test) - expected), 0,
                Real(1e-11) * (Real(1) + std::abs(expected)));
            }
          }
          EXPECT_GT(expectedSquaredNorm, 0);
        };
        auto shared = Integral(Dot(Dot(A, Jacobian(u)), Dot(A, Jacobian(v))));
        check(shared, true);
        const auto distinctIntegrand = Dot(Dot(A, Jacobian(u)), Dot(B, Jacobian(v)));
        using Rule = QuadratureRule<std::decay_t<decltype(distinctIntegrand)>>;
        static_assert(Rule::IsRankOneJacobianForm);
        auto distinct = Integral(distinctIntegrand);
        check(distinct, false);
      }
    }
  }

  /// @brief P1 rank-one forms retain ambient directions on embedded cells.
  TEST(Rodin_Variational_RankOneJacobianForm, EmbeddedP1Real)
  {
    checkEmbeddedRankOne<Real, false>();
  }

  /// @brief Complex P1 rank-one forms retain ambient directions and test conjugation.
  TEST(Rodin_Variational_RankOneJacobianForm, EmbeddedP1Complex)
  {
    checkEmbeddedRankOne<Complex, false>();
  }

  /// @brief H1 rank-one forms retain ambient directions on embedded cells.
  TEST(Rodin_Variational_RankOneJacobianForm, EmbeddedH1Real)
  {
    checkEmbeddedRankOne<Real, true>();
  }

  /// @brief Complex H1 rank-one forms retain ambient directions and test conjugation.
  TEST(Rodin_Variational_RankOneJacobianForm, EmbeddedH1Complex)
  {
    checkEmbeddedRankOne<Complex, true>();
  }

  /**
   * @brief The rank-one Jacobian form selects its own quadrature rule.
   *
   * A specialization that fails to match is not an error: the expression falls
   * back to the generic rule, which evaluates the coefficient once per basis
   * pair instead of once per quadrature point. That degrades silently, so the
   * selection is asserted at compile time rather than inferred from timing.
   */
  TEST(Rodin_Variational_RankOneJacobianForm, SpecializationIsSelected)
  {
    auto mesh = unitSquare(2);
    P1 fes(mesh, 2);
    TrialFunction u(fes);
    TestFunction v(fes);
    const IdentityMatrix A(2);

    auto integrand = Dot(Dot(A, Jacobian(u)), Dot(A, Jacobian(v)));
    using IntegrandType = std::decay_t<decltype(integrand)>;
    using Rule = QuadratureRule<IntegrandType>;

    static_assert(Rule::IsRankOneJacobianForm,
      "the rank-one Jacobian specialization was not selected; the expression "
      "would fall back to the generic rule and lose coefficient hoisting");
    SUCCEED();
  }

  /// @brief The same expression selects the H1 specialization for k greater than one.
  TEST(Rodin_Variational_RankOneJacobianForm, H1SpecializationIsSelected)
  {
    auto mesh = unitSquare(2);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, 2);
    TrialFunction u(fes);
    TestFunction v(fes);
    const IdentityMatrix A(2);

    auto integrand = Dot(Dot(A, Jacobian(u)), Dot(A, Jacobian(v)));
    using IntegrandType = std::decay_t<decltype(integrand)>;
    using Rule = QuadratureRule<IntegrandType>;

    static_assert(Rule::IsRankOneJacobianForm,
      "the rank-one H1 Jacobian specialization was not selected");
    SUCCEED();
  }

  /**
   * @brief With @f$A=I@f$ the rank-one form reduces to the squared trace.
   *
   * @f$I:\nabla u = \operatorname{tr}\nabla u = \operatorname{div} u@f$, so the
   * form is @f$\int (\operatorname{div}u)(\operatorname{div}v)@f$, which is
   * assembled independently by the divergence form and must agree.
   */
  TEST(Rodin_Variational_RankOneJacobianForm, IdentityCoefficientIsDivergenceForm)
  {
    auto mesh = unitSquare(3);
    P1 fes(mesh, 2);
    expectIdentityCoefficientEqualsDivergence(fes, 2);
  }

  /// @brief Distinct same-typed P1 coefficients are evaluated on their respective sides.
  TEST(Rodin_Variational_RankOneJacobianForm, P1DistinctCoefficientsUseBothSides)
  {
    auto mesh = unitSquare(3);
    P1 fes(mesh, 2);
    expectDistinctCoefficientsUseBothSides(fes, 2);
  }

  /// @brief Distinct same-typed H1 coefficients are evaluated on their respective sides.
  TEST(Rodin_Variational_RankOneJacobianForm, H1DistinctCoefficientsUseBothSides)
  {
    auto mesh = unitSquare(2);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, 2);
    expectDistinctCoefficientsUseBothSides(fes, 2);
  }

  /// @brief H1<2> tabulated gradients produce the divergence form in 2D.
  TEST(Rodin_Variational_RankOneJacobianForm, H1Quadratic2DIsDivergenceForm)
  {
    auto mesh = unitSquare(2);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, 2);
    expectIdentityCoefficientEqualsDivergence(fes, 2);
  }

  /// @brief H1<2> tabulated gradients produce the divergence form in 3D.
  TEST(Rodin_Variational_RankOneJacobianForm, H1Quadratic3DIsDivergenceForm)
  {
    auto mesh = unitCube(2);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, 3);
    expectIdentityCoefficientEqualsDivergence(fes, 3);
  }

  /// @brief H1<2> tabulation is correct on every supported cell geometry.
  TEST(Rodin_Variational_RankOneJacobianForm, H1QuadraticAcrossGeometries)
  {
    for (const auto geometry :
      {Polytope::Type::Segment, Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
        Polytope::Type::Tetrahedron, Polytope::Type::Wedge, Polytope::Type::Pyramid})
    {
      auto mesh = uniformGrid(geometry, 2);
      const size_t dimension = mesh.getDimension();
      H1 fes(std::integral_constant<size_t, 2>{}, mesh, dimension);
      expectIdentityCoefficientEqualsDivergence(fes, dimension);
    }
  }

  /// @brief Different H1 trial and test orders use their respective tabulations.
  TEST(Rodin_Variational_RankOneJacobianForm, H1MixedOrderIsDivergenceForm)
  {
    auto mesh = unitSquare(2);
    H1 trialFES(std::integral_constant<size_t, 2>{}, mesh, 2);
    H1 testFES(std::integral_constant<size_t, 3>{}, mesh, 2);
    TrialFunction u(trialFES);
    TestFunction v(testFES);
    const IdentityMatrix A(2);

    BilinearForm rankOne(u, v);
    rankOne = Integral(Dot(Dot(A, Jacobian(u)), Dot(A, Jacobian(v))));
    rankOne.assemble();

    BilinearForm divergence(u, v);
    divergence = Integral(Div(u), Div(v));
    divergence.assemble();

    const Math::Matrix<Real> actual = rankOne.getOperator();
    const Math::Matrix<Real> expected = divergence.getOperator();
    ASSERT_EQ(actual.rows(), expected.rows());
    ASSERT_EQ(actual.cols(), expected.cols());
    for (Eigen::Index i = 0; i < actual.rows(); ++i)
    {
      for (Eigen::Index j = 0; j < actual.cols(); ++j)
      {
        EXPECT_NEAR(actual(i, j), expected(i, j), 1e-12)
          << "entry (" << i << ", " << j << ")";
      }
    }
  }

  /// @brief The assembled form is symmetric positive semi-definite.
  TEST(Rodin_Variational_RankOneJacobianForm, IsSymmetricPositiveSemiDefinite)
  {
    auto mesh = unitSquare(3);
    P1 fes(mesh, 2);
    TrialFunction u(fes);
    TestFunction v(fes);
    const IdentityMatrix A(2);

    BilinearForm bf(u, v);
    bf = Integral(Dot(Dot(A, Jacobian(u)), Dot(A, Jacobian(v))));
    bf.assemble();

    const Math::Matrix<Real> m = bf.getOperator();
    for (Eigen::Index i = 0; i < m.rows(); ++i)
    {
      EXPECT_GE(m(i, i), -1e-14) << "diagonal " << i;
      for (Eigen::Index j = 0; j < i; ++j)
        EXPECT_NEAR(m(i, j), m(j, i), 1e-12);
    }

    // x^T M x = sum over quadrature of (A : J u_x)^2 >= 0.
    Math::Vector<Real> x(m.rows());
    for (Eigen::Index i = 0; i < m.rows(); ++i)
      x(i) = std::sin(Real(i) * Real(0.7)) + Real(0.3);
    EXPECT_GE(x.dot(m * x), -1e-12);
  }

  /// @brief An explicit order is honoured by the specialization.
  TEST(Rodin_Variational_RankOneJacobianForm, HonoursExplicitOrder)
  {
    auto mesh = unitSquare(2);
    P1 fes(mesh, 2);
    TrialFunction u(fes);
    TestFunction v(fes);
    const IdentityMatrix A(2);

    auto integ = Integral(Dot(Dot(A, Jacobian(u)), Dot(A, Jacobian(v))));
    const auto it = mesh.getPolytope(2, 0);
    EXPECT_FALSE(integ.getOrder(*it).has_value());
    integ.setOrder(6);
    EXPECT_EQ(integ.getOrder(*it).value_or(0), 6u);
  }
}
