/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Complex named-form deduction and implemented-rule equivalence.
 */
#include <gtest/gtest.h>
#include <type_traits>
#include "Rodin/Variational.h"
#include "../../../convergence/Convergence.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  namespace
  {
    /// @brief Existing named-form relative Frobenius budget on dimensionless meshes.
    constexpr Real OperatorTolerance = 1e-12;

    /**
     * @brief Checks named operators against explicitly matched Integral rules.
     *
     * Nonzero imaginary coefficients must remain complex rather than being
     * offered to a numeric real-lifting overload. Mixed Helmholtz arguments
     * exercise both function/value directions. Constant-space diffusion is
     * identically zero. The Point case certifies algebraic values, not a
     * zero-dimensional PDE refinement rate. Trial/test functions and all
     * coefficient captures remain alive throughout assembly.
     *
     * @param space Scalar or vector complex space.
     */
    void checkForms(auto& space)
    {
      using Space = std::decay_t<decltype(space)>;
      using Range = typename FormLanguage::Traits<Space>::RangeType;
      auto& mesh = space.getMesh();
      const auto cell = *mesh.getCell();
      const size_t degree =
        space.getFiniteElement(cell.getDimension(), cell.getIndex()).getOrder();
      const size_t massOrder = 2 * degree;
      const size_t diffusionOrder = degree == 0 ? 0 : 2 * degree - 2;
      const ComplexFunction diffusion([](const Point& point) {
        Real sum = 0;
        for (size_t j = 0; j < point.getDimension(); ++j)
          sum += point(j);
        return Complex(2 + sum, 0.5 + Real(0.25) * sum);
      });
      const ComplexFunction reaction([](const Point& point) {
        Real sum = 0;
        for (size_t j = 0; j < point.getDimension(); ++j)
          sum += point(j) * point(j);
        return Complex(3 + sum, -0.25 - Real(0.125) * sum);
      });
      TrialFunction trial(space);
      TestFunction test(space);
      const auto compare = [&](auto&& named, const auto& reference) {
        BilinearForm expected(trial, test);
        expected = reference;
        expected.assemble();
        const Math::Matrix<Complex> actual = named.getOperator();
        const Math::Matrix<Complex> oracle = expected.getOperator();
        ASSERT_EQ(actual.rows(), oracle.rows());
        ASSERT_EQ(actual.cols(), oracle.cols());
        const Real error = (actual - oracle).norm() / std::max(Real(1), oracle.norm());
        EXPECT_LE(error, OperatorTolerance)
          << "first named entry " << actual(0, 0) << ", Integral entry " << oracle(0, 0);
        named.assemble();
        const Math::Matrix<Complex> reassembled = named.getOperator();
        EXPECT_LE((reassembled - oracle).norm() / std::max(Real(1), oracle.norm()),
          OperatorTolerance);
      };
      MassForm mass(reaction, trial, test);
      const auto massIntegral = [&](const auto& coefficient) {
        auto integral = Integral(coefficient * trial, test);
        integral.setOrder(massOrder);
        return integral;
      };
      compare(mass, massIntegral(reaction));
      EXPECT_GT(Math::Matrix<Complex>(mass.getOperator()).imag().norm(), 0);
      if (mesh.getDimension() == 0)
      {
        // Independent algebraic oracle: the single Point basis equals one.
        EXPECT_EQ(Math::Matrix<Complex>(mass.getOperator())(0, 0), Complex(17, -2));
      }
      // Numeric lifting remains available alongside the complex function path.
      compare(MassForm(Real(3), trial, test), massIntegral(Real(3)));
      if constexpr (std::is_same_v<Range, Complex>)
      {
        const auto diffusionIntegral = [&](const auto& coefficient) {
          auto integral = Integral(coefficient * Grad(trial), Grad(test));
          integral.setOrder(diffusionOrder);
          return integral;
        };
        DiffusionForm gradient(diffusion, trial, test);
        compare(gradient, diffusionIntegral(diffusion));
        compare(DiffusionForm(Real(2), trial, test), diffusionIntegral(Real(2)));
        const auto reference = [&](const auto& first, const auto& zeroth) {
          return diffusionIntegral(first) + massIntegral(zeroth);
        };
        compare(HelmholtzForm(diffusion, reaction, trial, test),
          reference(diffusion, reaction));
        compare(
          HelmholtzForm(diffusion, Real(3), trial, test), reference(diffusion, Real(3)));
        compare(
          HelmholtzForm(Real(2), reaction, trial, test), reference(Real(2), reaction));
        compare(
          HelmholtzForm(Real(2), Real(3), trial, test), reference(Real(2), Real(3)));
      }
    }
  }

  class ComplexNamedFormTest : public ::testing::TestWithParam<Polytope::Type>
  {
    protected:
    LocalMesh makeMesh() const
    {
      if (GetParam() == Polytope::Type::Point)
        return LocalMesh::Builder().initialize(3).nodes(1).vertex({1, 2, 3}).finalize();
      return Convergence::UniformGrid(GetParam()).makeMesh(2);
    }
  };

  TEST_P(ComplexNamedFormTest, ScalarComplexFunctions)
  {
    auto mesh = makeMesh();
    P0<Complex> constant(mesh);
    P0g<Complex, LocalMesh> globalConstant(mesh);
    P1<Complex> linear(mesh);
    H1<1, Complex> orderOne(std::integral_constant<size_t, 1>{}, mesh);
    H1<2, Complex> quadratic(std::integral_constant<size_t, 2>{}, mesh);
    checkForms(constant);
    checkForms(globalConstant);
    checkForms(linear);
    checkForms(orderOne);
    checkForms(quadratic);
  }

  TEST_P(ComplexNamedFormTest, VectorComplexMass)
  {
    auto mesh = makeMesh();
    using Vector = Math::SpatialVector<Complex>;
    const size_t components = mesh.getSpaceDimension();
    P0<Vector> constant(mesh, components);
    P0g<Vector, LocalMesh> globalConstant(mesh, components);
    P1<Vector> linear(mesh, components);
    H1<1, Vector> orderOne(std::integral_constant<size_t, 1>{}, mesh, components);
    H1<2, Vector> quadratic(std::integral_constant<size_t, 2>{}, mesh, components);
    checkForms(constant);
    checkForms(globalConstant);
    checkForms(linear);
    checkForms(orderOne);
    checkForms(quadratic);
  }

  /**
   * @brief Constant bases have exactly zero broken diffusion entries.
   *
   * Real and complex cellwise/global constants exercise the named kernel on
   * Point and every cell geometry. This algebraic identity requires neither a
   * numerical tolerance nor a positive-definiteness or PDE-solve assertion.
   */
  TEST_P(ComplexNamedFormTest, ConstantNamedDiffusionIsZero)
  {
    auto mesh = makeMesh();
    const auto check = [&](auto&& space) {
      using Scalar =
        typename FormLanguage::Traits<std::decay_t<decltype(space)>>::ScalarType;
      TrialFunction trial(space);
      TestFunction test(space);
      DiffusionForm form(Real(2), trial, test);
      for (size_t repeat = 0; repeat < 2; ++repeat)
      {
        form.assemble();
        const Math::Matrix<Scalar> matrix = form.getOperator();
        EXPECT_GT(matrix.rows(), 0);
        for (Eigen::Index row = 0; row < matrix.rows(); ++row)
        {
          for (Eigen::Index col = 0; col < matrix.cols(); ++col)
          {
            EXPECT_EQ(matrix(row, col), Scalar(0));
          }
        }
      }
    };
    check(P0<Real>(mesh));
    check(P0g<Real, LocalMesh>(mesh));
    check(P0<Complex>(mesh));
    check(P0g<Complex, LocalMesh>(mesh));
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, ComplexNamedFormTest,
    ::testing::Values(Polytope::Type::Point, Polytope::Type::Segment,
      Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron, Polytope::Type::Pyramid, Polytope::Type::Hexahedron,
      Polytope::Type::Wedge));
}
