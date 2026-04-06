/**
 * @file ReImTest.cpp
 * @brief Tests for the Re and Im complex part extractors.
 */
#include <gtest/gtest.h>
#include <cmath>
#include <complex>

#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

// --- Re tests ---

TEST(Rodin_Variational_Re, RealPart)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  ComplexFunction f(Complex(3.0, 4.0));
  auto re = Re(f);

  EXPECT_NEAR(re.getValue(p), 3.0, 1e-10);
}

TEST(Rodin_Variational_Re, PureImaginary)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  ComplexFunction f(Complex(0.0, 5.0));
  auto re = Re(f);

  EXPECT_NEAR(re.getValue(p), 0.0, 1e-10);
}

TEST(Rodin_Variational_Re, NegativeValues)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  ComplexFunction f(Complex(-2.5, 1.0));
  auto re = Re(f);

  EXPECT_NEAR(re.getValue(p), -2.5, 1e-10);
}

TEST(Rodin_Variational_Re, Copy)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  ComplexFunction f(Complex(7.0, 8.0));
  auto re = Re(f);
  auto copy = re;

  EXPECT_NEAR(copy.getValue(p), 7.0, 1e-10);
}

// --- Im tests ---

TEST(Rodin_Variational_Im, ImaginaryPart)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  ComplexFunction f(Complex(3.0, 4.0));
  auto im = Im(f);

  EXPECT_NEAR(im.getValue(p), 4.0, 1e-10);
}

TEST(Rodin_Variational_Im, PureReal)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  ComplexFunction f(Complex(5.0, 0.0));
  auto im = Im(f);

  EXPECT_NEAR(im.getValue(p), 0.0, 1e-10);
}

TEST(Rodin_Variational_Im, NegativeValues)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  ComplexFunction f(Complex(1.0, -6.0));
  auto im = Im(f);

  EXPECT_NEAR(im.getValue(p), -6.0, 1e-10);
}

TEST(Rodin_Variational_Im, Copy)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  ComplexFunction f(Complex(7.0, 8.0));
  auto im = Im(f);
  auto copy = im;

  EXPECT_NEAR(copy.getValue(p), 8.0, 1e-10);
}

// --- Re + Im identity test ---

TEST(Rodin_Variational_ReIm, ReImReconstruct)
{
  // Test: Re(f)^2 + Im(f)^2 = |f|^2
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  ComplexFunction f(Complex(3.0, 4.0));
  auto re = Re(f);
  auto im = Im(f);

  Real reVal = re.getValue(p);
  Real imVal = im.getValue(p);

  // |f|^2 = 3^2 + 4^2 = 25
  EXPECT_NEAR(reVal * reVal + imVal * imVal, 25.0, 1e-10);
}
