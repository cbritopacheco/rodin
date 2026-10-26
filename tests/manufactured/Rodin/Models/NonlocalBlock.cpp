/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Manufactured mixed nonlocal solves for real/complex dense/sparse systems.
 */
#include <gtest/gtest.h>
#include "Rodin/Variational.h"
#include "Rodin/Solver/SparseLU.h"
#include "Rodin/Solver/PartialPivLU.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Manufactured
{
  /**
   * @brief Two exact constant fields solve a mass-plus-global-coupling system.
   * @tparam Scalar Real or complex coefficient type.
   * @tparam Operator Dense or sparse matrix storage.
   * @param preassembledConstrained Use preassembled local forms and exact boundary data.
   */
  template <class Scalar, class Operator>
  void checkManufacturedNonlocalBlock(bool preassembledConstrained = false)
  {
    constexpr Real Tolerance = Real(1e-11);
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    mesh.getConnectivity().compute(0, 1);
    P1<Scalar> fes(mesh);
    TrialFunction u(fes), p(fes);
    TestFunction v(fes), q(fes);
    using System = Math::LinearSystem<Operator, Math::Vector<Scalar>>;
    Problem<System, decltype(u), decltype(v), decltype(p), decltype(q)> problem(
      u, v, p, q);
    const Scalar factor = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(1, -2);
      else
        return Scalar(1);
    }();
    const Scalar exactU = Scalar(2) * factor;
    const Scalar exactP = Scalar(-3) / factor;
    auto constant = [](Scalar value) {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return ComplexFunction(value);
      else
        return RealFunction(value);
    };
    auto kernel = [](Scalar value) {
      return [value](const Point&, const Point&) { return value; };
    };
    const Scalar c00 = Scalar(2) * factor, c01 = Scalar(-3) * factor;
    const Scalar c10 = Scalar(0.5) * factor, c11 = Scalar(4) * factor;
    // The interval has length two: integrating either exact field gives twice its value.
    const Scalar loadU = exactU + Scalar(2) * (c00 * exactU + c01 * exactP);
    const Scalar loadP = exactP + Scalar(2) * (c10 * exactU + c11 * exactP);
    typename decltype(problem)::ProblemBodyType body;
    if (preassembledConstrained)
    {
      using Form = BilinearForm<typename decltype(u)::SolutionType, decltype(fes),
        decltype(fes), Operator>;
      Form massU(u, v), massP(p, q), coupling(p, v);
      massU = Integral(u, v);
      massP = Integral(p, q);
      coupling = Integral(Scalar(2) * p, v);
      massU.assemble();
      massP.assemble();
      coupling.assemble();
      LinearForm loadFirst(v), loadSecond(q);
      loadFirst = Integral(constant(loadU + Scalar(2) * exactP), v);
      loadSecond = Integral(constant(loadP), q);
      loadFirst.assemble();
      loadSecond.assemble();
      body.getBFs().add(massU);
      body.getLFs().add(loadFirst);
      body.getLFs().add(loadSecond);
      body.getDBCs().add(DirichletBC(u, constant(exactU)));
      body.getDBCs().add(DirichletBC(p, constant(exactP)));
      body.getBFs().add(massP);
      body.getBFs().add(coupling);
    }
    else
    {
      body = Integral(u, v) + Integral(p, q) - Integral(constant(loadU), v) -
        Integral(constant(loadP), q);
    }
    // Potential borrows its kernel callable; retain all four until the solve completes.
    const auto k00 = kernel(c00), k01 = kernel(c01), k10 = kernel(c10), k11 = kernel(c11);
    body.getGlobalBFIs().add(Integral(Potential(k00, u), v));
    body.getGlobalBFIs().add(Integral(Potential(k01, p), v));
    body.getGlobalBFIs().add(Integral(Potential(k10, u), q));
    body.getGlobalBFIs().add(Integral(Potential(k11, p), q));
    problem = body;
    if constexpr (std::is_base_of_v<Eigen::SparseMatrixBase<Operator>, Operator>)
    {
      Solver::SparseLU solver(problem);
      solver.solve();
      ASSERT_TRUE(solver.success());
    }
    else
    {
      Solver::PartialPivLU solver(problem);
      solver.solve();
    }
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      const Point point(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
      EXPECT_LT(std::abs(u.getSolution()(point) - exactU), Tolerance);
      EXPECT_LT(std::abs(p.getSolution()(point) - exactP), Tolerance);
    }
  }

  /// @brief Real sparse storage recovers both independently prescribed constant fields.
  TEST(Rodin_Manufactured_NonlocalBlock, SparseReal)
  {
    checkManufacturedNonlocalBlock<Real, Math::SparseMatrix<Real>>();
  }

  /// @brief Complex sparse storage preserves nonreal kernels and solution components.
  TEST(Rodin_Manufactured_NonlocalBlock, SparseComplex)
  {
    checkManufacturedNonlocalBlock<Complex, Math::SparseMatrix<Complex>>();
  }

  /// @brief Real dense storage has the same manufactured physical solution.
  TEST(Rodin_Manufactured_NonlocalBlock, DenseReal)
  {
    checkManufacturedNonlocalBlock<Real, Math::Matrix<Real>>();
  }

  /// @brief Complex dense storage has the same manufactured physical solution.
  TEST(Rodin_Manufactured_NonlocalBlock, DenseComplex)
  {
    checkManufacturedNonlocalBlock<Complex, Math::Matrix<Complex>>();
  }

  /// @brief Dense real preassembled blocks and loads retain exact boundary elimination.
  TEST(Rodin_Manufactured_NonlocalBlock, DenseRealPreassembledBoundary)
  {
    checkManufacturedNonlocalBlock<Real, Math::Matrix<Real>>(true);
  }

  /// @brief Dense complex preassembled blocks preserve nonreal boundary values and coupling.
  TEST(Rodin_Manufactured_NonlocalBlock, DenseComplexPreassembledBoundary)
  {
    checkManufacturedNonlocalBlock<Complex, Math::Matrix<Complex>>(true);
  }
}

namespace
{
  /// @brief Subtracting a constant nonlocal kernel retains an exact constant solution.
  template <class Scalar>
  void checkSubtractedNonlocalMass()
  {
    const Scalar coefficient = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(0.125, 0.25);
      else
        return Scalar(0.125);
    }();
    const Scalar exact = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(2, -3);
      else
        return Scalar(2);
    }();
    const auto kernel = [coefficient](const Point&, const Point&) { return coefficient; };
    const auto load = [exact, coefficient] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return ComplexFunction(exact * (Scalar(1) - coefficient));
      else
        return RealFunction(exact * (Scalar(1) - coefficient));
    }();
    for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      const size_t dimension = Polytope::Traits(geometry).getDimension();
      for (const size_t nodes : {size_t(2), size_t(3)})
      {
        auto mesh =
          LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, nodes));
        mesh.scale(Real(1) / (nodes - 1));
        P1<Scalar> space(mesh);
        TrialFunction u(space);
        TestFunction v(space);
        BilinearForm form(u, v);
        form = Integral(u, v);
        form -= Integral(Potential(kernel, u), v);
        form.assemble();
        Problem problem(u, v);
        typename decltype(problem)::ProblemBodyType body;
        body.getBFs().add(form);
        body.getLFIs().add(-Integral(load, v));
        problem = body;
        Solver::SparseLU solver(problem);
        solver.solve();
        ASSERT_TRUE(solver.success());
        const auto& solution = u.getSolution();
        for (size_t i = 0; i < solution.getSize(); ++i)
          EXPECT_LT(std::abs(solution[i] - exact), 1e-10);
      }
    }
  }
}

TEST(Rodin_Manufactured_NonlocalSubtraction, RealConstantSolution)
{
  checkSubtractedNonlocalMass<Real>();
}

TEST(Rodin_Manufactured_NonlocalSubtraction, ComplexConstantSolution)
{
  checkSubtractedNonlocalMass<Complex>();
}
