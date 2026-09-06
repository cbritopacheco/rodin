#include <gtest/gtest.h>

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
        for (Eigen::Index j = 0; j < actual.cols(); ++j)
          EXPECT_NEAR(actual(i, j), expected(i, j), 1e-12)
            << "entry (" << i << ", " << j << ")";
    }
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
      for (Eigen::Index j = 0; j < actual.cols(); ++j)
        EXPECT_NEAR(actual(i, j), expected(i, j), 1e-12)
          << "entry (" << i << ", " << j << ")";
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
