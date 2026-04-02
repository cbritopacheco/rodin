/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2024.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Variational_Trace, IdentityMatrix_2D)
  {
    // Trace of the 2x2 identity matrix should be 2
    auto I = IdentityMatrix(2);
    auto tr = Trace(I);

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const Math::Vector<Real> rc{{ 0.5, 0.5 }};
    Point p(*it, rc);

    EXPECT_NEAR(tr.getValue(p), 2.0, 1e-10);
  }

  TEST(Rodin_Variational_Trace, IdentityMatrix_3D)
  {
    // Trace of the 3x3 identity matrix should be 3
    auto I = IdentityMatrix(3);
    auto tr = Trace(I);

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const Math::Vector<Real> rc{{ 0.3, 0.3 }};
    Point p(*it, rc);

    EXPECT_NEAR(tr.getValue(p), 3.0, 1e-10);
  }

  TEST(Rodin_Variational_Trace, ScaledIdentityMatrix)
  {
    // Trace of alpha * I should be alpha * n
    RealFunction alpha(5.0);
    auto I = IdentityMatrix(2);
    auto scaled = alpha * I;
    auto tr = Trace(scaled);

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const Math::Vector<Real> rc{{ 0.4, 0.4 }};
    Point p(*it, rc);

    EXPECT_NEAR(tr.getValue(p), 10.0, 1e-10);
  }

  TEST(Rodin_Variational_Trace, CopyConstruction)
  {
    auto I = IdentityMatrix(2);
    auto tr = Trace(I);
    auto tr_copy = tr;

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const Math::Vector<Real> rc{{ 0.5, 0.5 }};
    Point p(*it, rc);

    EXPECT_NEAR(tr_copy.getValue(p), 2.0, 1e-10);
  }

  TEST(Rodin_Variational_Trace, MoveConstruction)
  {
    auto I = IdentityMatrix(2);
    auto tr = Trace(I);
    auto tr_moved = std::move(tr);

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const Math::Vector<Real> rc{{ 0.5, 0.5 }};
    Point p(*it, rc);

    EXPECT_NEAR(tr_moved.getValue(p), 2.0, 1e-10);
  }

  TEST(Rodin_Variational_Trace, PolymorphicCopy)
  {
    auto I = IdentityMatrix(3);
    auto tr = Trace(I);
    std::unique_ptr<decltype(tr)> copy(tr.copy());

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const Math::Vector<Real> rc{{ 0.2, 0.2 }};
    Point p(*it, rc);

    EXPECT_NEAR(copy->getValue(p), 3.0, 1e-10);
  }
}
