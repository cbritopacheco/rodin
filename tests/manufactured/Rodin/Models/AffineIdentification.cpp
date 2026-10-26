/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Exact mixed mass solutions with affine identification and fixed master data.
 */
#include <gtest/gtest.h>
#include "Rodin/Variational.h"
#include "Rodin/Solver/PartialPivLU.h"
#include "Rodin/Solver/SparseLU.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  /// @brief Mixed constraints recover independently prescribed constant fields.
  template <class Scalar, class Operator>
  void checkAffineIdentification()
  {
    const Scalar exactP = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(-1, 0.5);
      else
        return Scalar(-1);
    }();
    const Scalar exactU = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(17.0 / 8, 2.5);
      else
        return Scalar(2);
    }();
    const Scalar gamma = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(2, -0.25);
      else
        return Scalar(2);
    }();
    const Scalar defect = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(4, 1.25);
      else
        return Scalar(4);
    }();
    auto constant = [](Scalar value) {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return ComplexFunction(value);
      else
        return RealFunction(value);
    };
    for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      const size_t dimension = Polytope::Traits(geometry).getDimension();
      for (const size_t nodes : {size_t(2), size_t(3)})
      {
        SCOPED_TRACE(nodes);
        auto mesh =
          LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, nodes));
        mesh.scale(Real(1) / (nodes - 1));
        mesh.getConnectivity().compute(dimension - 1, dimension);
        for (auto face = mesh.getBoundary(); face; ++face)
        {
          const Point center(*face, Polytope::Traits(face->getGeometry()).getCentroid());
          mesh.setAttribute(
            {dimension - 1, face->getIndex()}, std::abs(center(0)) < 1e-12 ? 31 : 37);
        }
        P1<Scalar> space(mesh);
        TrialFunction u(space), p(space);
        TestFunction v(space), q(space);
        using System = Math::LinearSystem<Operator, Math::Vector<Scalar>>;
        Problem<System, decltype(u), decltype(v), decltype(p), decltype(q)> problem(
          u, v, p, q);
        auto identification = DirichletBC(u, constant(gamma) * p, constant(defect));
        identification.on({31, 37});
        auto fixed = DirichletBC(p, constant(exactP));
        fixed.on(31);
        problem = Integral(u, v) + Integral(p, q) - Integral(constant(exactU), v) -
          Integral(constant(exactP), q) + identification + fixed;
        auto check = [&](auto& instance) {
          if constexpr (std::is_base_of_v<Eigen::SparseMatrixBase<Operator>, Operator>)
          {
            Solver::SparseLU solver(instance);
            solver.solve();
            ASSERT_TRUE(solver.success());
          }
          else
          {
            Solver::PartialPivLU solver(instance);
            solver.solve();
          }
          for (size_t i = 0; i < space.getSize(); ++i)
          {
            EXPECT_LT(std::abs(u.getSolution()[i] - exactU), 1e-10);
            EXPECT_LT(std::abs(p.getSolution()[i] - exactP), 1e-10);
          }
        };
        check(problem);
        check(problem);
        auto copied = problem;
        check(copied);
        auto moved = std::move(copied);
        check(moved);
      }
    }
  }
}

TEST(Rodin_Manufactured_AffineIdentification, SparseReal)
{
  checkAffineIdentification<Real, Math::SparseMatrix<Real>>();
}

TEST(Rodin_Manufactured_AffineIdentification, DenseReal)
{
  checkAffineIdentification<Real, Math::Matrix<Real>>();
}

TEST(Rodin_Manufactured_AffineIdentification, SparseComplex)
{
  checkAffineIdentification<Complex, Math::SparseMatrix<Complex>>();
}

TEST(Rodin_Manufactured_AffineIdentification, DenseComplex)
{
  checkAffineIdentification<Complex, Math::Matrix<Complex>>();
}
