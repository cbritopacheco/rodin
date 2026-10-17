/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Local quadrature ownership checked against analytic field contractions. */
#include <gtest/gtest.h>
#include <functional>
#include <memory>
#include <vector>
#include <Rodin/Variational.h>
#include <Rodin/Geometry.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  template <class Scalar, bool HigherOrder>
  void checkLocalQuadratureLifecycle(size_t dimension)
  {
    auto mesh = dimension == 1 ? LocalMesh::UniformGrid(Polytope::Type::Segment, {3})
                               : LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    mesh.getConnectivity().compute(1, 0);
    if (dimension == 2)
    {
      mesh.getConnectivity().compute(1, 2);
      mesh.getConnectivity().compute(2, 1);
    }
    using ScalarSpace = std::conditional_t<HigherOrder, H1<2, Scalar>, P1<Scalar>>;
    using VectorSpace = std::conditional_t<HigherOrder,
      H1<2, Math::SpatialVector<Scalar>>, P1<Math::SpatialVector<Scalar>>>;
    auto scalarSpace = [&] {
      if constexpr (HigherOrder)
        return ScalarSpace(std::integral_constant<size_t, 2>{}, mesh);
      else
        return ScalarSpace(mesh);
    }();
    auto vectorSpace = [&] {
      if constexpr (HigherOrder)
        return VectorSpace(std::integral_constant<size_t, 2>{}, mesh, dimension);
      else
        return VectorSpace(mesh, dimension);
    }();
    auto otherVectorSpace = [&] {
      if constexpr (HigherOrder)
        return VectorSpace(std::integral_constant<size_t, 2>{}, mesh, dimension);
      else
        return VectorSpace(mesh, dimension);
    }();
    const Scalar factor = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(1, -2);
      else
        return Scalar(1);
    }();
    const auto scalarFunction = [&] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return ComplexFunction([factor, dimension](const Point& p) {
          return factor * (1 + p.x() + (dimension == 2 ? 2 * p.y() : 0));
        });
      else
        return RealFunction([factor, dimension](const Point& p) {
          return factor * (1 + p.x() + (dimension == 2 ? 2 * p.y() : 0));
        });
    }();
    const auto vectorFunction =
      VectorFunction(dimension, [factor, dimension](const Point& p) {
        Math::SpatialVector<Scalar> value(dimension);
        value[0] = factor * (1 + 3 * p.x() - (dimension == 2 ? p.y() : 0));
        if (dimension == 2)
          value[1] = factor * (-2 + 2 * p.x() + 4 * p.y());
        return value;
      });
    GridFunction scalarField(scalarSpace), vectorField(vectorSpace),
      otherVectorField(otherVectorSpace);
    scalarField = scalarFunction;
    vectorField = vectorFunction;
    otherVectorField = vectorFunction;
    TrialFunction u(scalarSpace), w(vectorSpace);
    TestFunction v(scalarSpace), z(vectorSpace), otherZ(otherVectorSpace);
    RealFunction coefficient([](const Point& p) { return 1 + p.x(); });
    using BilinearBase = LocalBilinearFormIntegratorBase<Scalar>;
    using LinearBase = LinearFormIntegratorBase<Scalar>;
    using BilinearCheck = std::function<void(BilinearBase&, const Polytope&)>;
    using LinearCheck = std::function<void(LinearBase&, const Polytope&)>;
    std::vector<std::pair<std::unique_ptr<BilinearBase>, BilinearCheck>> bilinearRules;
    std::vector<std::pair<std::unique_ptr<LinearBase>, LinearCheck>> linearRules;
    auto lifecycle = [&](auto rule, auto action, auto density) {
      using Rule = std::decay_t<decltype(rule)>;
      constexpr size_t order = 8;
      rule.setOrder(order);
      auto verify = [action, density, order](auto& integrator, const Polytope& cell) {
        EXPECT_EQ(&integrator.getPolytope(), &cell);
        EXPECT_EQ(integrator.getOrder(cell), order);
        const auto& formula =
          QF::PolytopeQuadratureFormula::get(order, cell.getGeometry());
        Scalar expected = 0;
        for (size_t sample = 0; sample < formula.getSize(); ++sample)
        {
          const Point point(cell, formula.getPoint(sample));
          expected += formula.getWeight(sample) * point.getDistortion() * density(point);
        }
        EXPECT_LT(std::abs(action(integrator, cell) - expected), 1e-10);
      };
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        rule.setPolytope(*cell);
        verify(rule, *cell);
        Rule copied(rule);
        copied.setPolytope(*cell);
        Rule moved(std::move(copied));
        verify(moved, *cell);
        std::unique_ptr<Rule> clone(rule.copy());
        clone->setPolytope(*cell);
        verify(*clone, *cell);
      }
      if constexpr (std::is_base_of_v<BilinearBase, Rule>)
        bilinearRules.emplace_back(std::unique_ptr<BilinearBase>(rule.copy()), verify);
      else
        linearRules.emplace_back(std::unique_ptr<LinearBase>(rule.copy()), verify);
    };
    auto bilinear = [&](auto rule, const auto& trial, const auto& test, auto density) {
      auto action = [&](auto& integrator, const Polytope& cell) {
        const auto& columns =
          trial.getFiniteElementSpace().getDOFs(dimension, cell.getIndex());
        const auto& rows =
          test.getFiniteElementSpace().getDOFs(dimension, cell.getIndex());
        Scalar value = 0;
        for (Index row = 0; row < rows.size(); ++row)
        {
          for (Index column = 0; column < columns.size(); ++column)
          {
            const Scalar testCoefficient = test.getData()[rows[row]];
            const Scalar conjugate = [&] {
              if constexpr (std::is_same_v<Scalar, Complex>)
                return std::conj(testCoefficient);
              else
                return testCoefficient;
            }();
            value += conjugate * integrator.integrate(column, row) *
              trial.getData()[columns[column]];
          }
        }
        return value;
      };
      lifecycle(std::move(rule), action, density);
    };
    auto linear = [&](auto rule, const auto& field, auto density) {
      auto action = [&](auto& integrator, const Polytope& cell) {
        const auto& dofs =
          field.getFiniteElementSpace().getDOFs(dimension, cell.getIndex());
        Scalar value = 0;
        for (Index local = 0; local < dofs.size(); ++local)
        {
          const Scalar testCoefficient = field.getData()[dofs[local]];
          if constexpr (std::is_same_v<Scalar, Complex>)
            value += std::conj(testCoefficient) * integrator.integrate(local);
          else
            value += testCoefficient * integrator.integrate(local);
        }
        return value;
      };
      lifecycle(std::move(rule), action, density);
    };
    const Real amplitude = std::norm(Complex(factor));
    const Real gradientNorm = dimension == 1 ? 1 : 5;
    const Real jacobianNorm = dimension == 1 ? 9 : 30;
    const Real divergence = dimension == 1 ? 3 : 7;
    auto scalarMass = [amplitude, dimension](const Point& p) {
      return Scalar(
        amplitude * std::pow(1 + p.x() + (dimension == 2 ? 2 * p.y() : 0), 2));
    };
    auto vectorMass = [amplitude, dimension](const Point& p) {
      const Real first = 1 + 3 * p.x() - (dimension == 2 ? p.y() : 0);
      const Real second = dimension == 2 ? -2 + 2 * p.x() + 4 * p.y() : 0;
      return Scalar(amplitude * (first * first + second * second));
    };
    bilinear(Integral(u, v), scalarField, scalarField, scalarMass);
    bilinear(Integral(w, z), vectorField, vectorField, vectorMass);
    bilinear(Integral(Grad(u), Grad(v)), scalarField, scalarField,
      [=](const Point&) { return Scalar(gradientNorm * amplitude); });
    bilinear(Integral(coefficient * Grad(u), Grad(v)), scalarField, scalarField,
      [=](const Point& p) { return Scalar(gradientNorm * amplitude * (1 + p.x())); });
    bilinear(Integral(Jacobian(w), Jacobian(otherZ)), vectorField, otherVectorField,
      [=](const Point&) { return Scalar(jacobianNorm * amplitude); });
    bilinear(Integral(coefficient * Jacobian(w), Jacobian(otherZ)), vectorField,
      otherVectorField,
      [=](const Point& p) { return Scalar(jacobianNorm * amplitude * (1 + p.x())); });
    bilinear(Integral(Div(w), Div(z)), vectorField, vectorField,
      [=](const Point&) { return Scalar(divergence * divergence * amplitude); });
    bilinear(Integral(coefficient * Div(w), Div(z)), vectorField, vectorField,
      [=](const Point& p) {
        return Scalar(divergence * divergence * amplitude * (1 + p.x()));
      });
    bilinear(Integral(Div(w), v), vectorField, scalarField, [=](const Point& p) {
      return Scalar(
        divergence * amplitude * (1 + p.x() + (dimension == 2 ? 2 * p.y() : 0)));
    });
    bilinear(Integral(u, Div(z)), scalarField, vectorField, [=](const Point& p) {
      return Scalar(
        divergence * amplitude * (1 + p.x() + (dimension == 2 ? 2 * p.y() : 0)));
    });
    const auto rankOneCoefficient =
      MatrixFunction(dimension, dimension, [dimension](const Point& p) {
        Math::SpatialMatrix<Real> value(dimension, dimension);
        value(0, 0) = 1 + p.x();
        if (dimension == 2)
        {
          value(0, 1) = 2;
          value(1, 0) = 3;
          value(1, 1) = 4;
        }
        return value;
      });
    if constexpr (std::is_same_v<Scalar, Real>)
    {
      auto density = [=](const Point& p) {
        const Real action = dimension == 1 ? 3 * (1 + p.x()) : 23 + 3 * p.x();
        return Scalar(amplitude * action * action);
      };
      bilinear(Integral(Dot(rankOneCoefficient, Jacobian(w)),
                 Dot(rankOneCoefficient, Jacobian(z))),
        vectorField, vectorField, density);
      bilinear(Integral(Dot(rankOneCoefficient, Jacobian(w)),
                 Dot(rankOneCoefficient, Jacobian(otherZ))),
        vectorField, otherVectorField, density);
    }
    if constexpr (std::is_same_v<Scalar, Complex>)
    {
      const auto complexCoefficient =
        MatrixFunction(dimension, dimension, [dimension](const Point& p) {
          Math::SpatialMatrix<Scalar> value(dimension, dimension);
          value(0, 0) = Scalar(0.5, 0.25) * (1 + p.x());
          if (dimension == 2)
          {
            value(0, 1) = Scalar(-0.5, 0.25);
            value(1, 0) = Scalar(0.375, -0.75);
            value(1, 1) = Scalar(1, -2);
          }
          return value;
        });
      auto density = [=](const Point& p) {
        const Scalar action = dimension == 1
          ? Scalar(1.5, 0.75) * (1 + p.x())
          : Scalar(6.75 + 1.5 * p.x(), -9 + 0.75 * p.x());
        return Scalar(amplitude * std::norm(action));
      };
      bilinear(Integral(Dot(complexCoefficient, Jacobian(w)),
                 Dot(complexCoefficient, Jacobian(z))),
        vectorField, vectorField, density);
      bilinear(Integral(Dot(complexCoefficient, Jacobian(w)),
                 Dot(complexCoefficient, Jacobian(otherZ))),
        vectorField, otherVectorField, density);
    }
    linear(Integral(scalarFunction, v), scalarField, scalarMass);
    linear(Integral(vectorFunction, z), vectorField, vectorMass);
    auto checkRuntimeRules = [&](auto& rules) {
      for (auto& [rule, verify] : rules)
      {
        for (auto cell = mesh.getCell(); cell; ++cell)
        {
          rule->setPolytope(*cell);
          verify(*rule, *cell);
          using Base = typename std::remove_reference_t<decltype(rule)>::element_type;
          std::unique_ptr<Base> copy(rule->copy());
          copy->setPolytope(*cell);
          verify(*copy, *cell);
        }
      }
    };
    checkRuntimeRules(bilinearRules);
    checkRuntimeRules(linearRules);
  }
}

TEST(QuadratureLifecycle, RealP1)
{
  for (size_t dimension : {size_t(1), size_t(2)})
    checkLocalQuadratureLifecycle<Real, false>(dimension);
}
TEST(QuadratureLifecycle, ComplexP1)
{
  for (size_t dimension : {size_t(1), size_t(2)})
    checkLocalQuadratureLifecycle<Complex, false>(dimension);
}
TEST(QuadratureLifecycle, RealH1)
{
  for (size_t dimension : {size_t(1), size_t(2)})
    checkLocalQuadratureLifecycle<Real, true>(dimension);
}
TEST(QuadratureLifecycle, ComplexH1)
{
  for (size_t dimension : {size_t(1), size_t(2)})
    checkLocalQuadratureLifecycle<Complex, true>(dimension);
}
