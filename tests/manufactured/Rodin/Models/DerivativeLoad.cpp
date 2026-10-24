/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Manufactured mass solves with coordinate-derivative loads.
 */
#include <gtest/gtest.h>
#include "Rodin/Variational.h"
#include "Rodin/Solver/SparseLU.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  /// @brief An affine coefficient produces an independently known constant solution.
  template <class Scalar>
  void checkDerivativeLoad()
  {
    const Scalar factor = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(2, -3);
      else
        return Scalar(2);
    }();
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
        P1<Scalar> space(mesh);
        GridFunction coefficient(space);
        coefficient = [factor, dimension](const Point& point) {
          Scalar value = 1;
          for (size_t axis = 0; axis < dimension; ++axis)
            value += factor * Scalar(axis + 1) * point(axis);
          return value;
        };
        TrialFunction u(space);
        TestFunction v(space);
        Problem problem(u, v);
        for (size_t axis = 0; axis < dimension; ++axis)
        {
          SCOPED_TRACE(axis);
          problem = Integral(u, v) - Integral(Derivative(axis, coefficient), v);
          Solver::SparseLU solver(problem);
          solver.solve();
          ASSERT_TRUE(solver.success());
          const auto& solution = u.getSolution();
          for (size_t i = 0; i < solution.getSize(); ++i)
            EXPECT_LT(std::abs(solution[i] - factor * Scalar(axis + 1)), 1e-10);
        }
      }
    }
  }
}

TEST(Rodin_Manufactured_DerivativeLoad, RealAffineCoefficient)
{
  checkDerivativeLoad<Real>();
}

TEST(Rodin_Manufactured_DerivativeLoad, ComplexAffineCoefficient)
{
  checkDerivativeLoad<Complex>();
}
