/**
 * @file ConjugateTest.cpp
 * @brief Tests for the Conjugate operator.
 */
#include <gtest/gtest.h>
#include <cmath>
#include <complex>

#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

TEST(Rodin_Variational_Conjugate, RealValueUnchanged)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  ComplexFunction f(Complex(3.0, 0.0));
  auto conj = Conjugate(f);

  Complex val = conj.getValue(p);
  EXPECT_NEAR(val.real(), 3.0, 1e-10);
  EXPECT_NEAR(val.imag(), 0.0, 1e-10);
}

TEST(Rodin_Variational_Conjugate, ImaginaryPart)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  ComplexFunction f(Complex(2.0, 5.0));
  auto conj = Conjugate(f);

  Complex val = conj.getValue(p);
  EXPECT_NEAR(val.real(), 2.0, 1e-10);
  EXPECT_NEAR(val.imag(), -5.0, 1e-10);
}

TEST(Rodin_Variational_Conjugate, PureImaginary)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  ComplexFunction f(Complex(0.0, 4.0));
  auto conj = Conjugate(f);

  Complex val = conj.getValue(p);
  EXPECT_NEAR(val.real(), 0.0, 1e-10);
  EXPECT_NEAR(val.imag(), -4.0, 1e-10);
}

TEST(Rodin_Variational_Conjugate, Copy)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  ComplexFunction f(Complex(1.0, 2.0));
  auto conj = Conjugate(f);
  auto copy = conj;

  Complex val = copy.getValue(p);
  EXPECT_NEAR(val.real(), 1.0, 1e-10);
  EXPECT_NEAR(val.imag(), -2.0, 1e-10);
}

TEST(Rodin_Variational_Conjugate, GetOrder)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);

  ComplexFunction f(Complex(3.0, -7.0));
  auto conj = Conjugate(f);
  auto order = conj.getOrder(*it);
  EXPECT_TRUE(order.has_value());
  EXPECT_EQ(*order, 0u);
}
