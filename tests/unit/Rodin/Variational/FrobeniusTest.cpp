/**
 * @file FrobeniusTest.cpp
 * @brief Tests for the Frobenius norm operator.
 */
#include <gtest/gtest.h>
#include <cmath>

#include "Rodin/Variational.h"
#include "Rodin/Variational/ComplexFunction.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

/// @brief Verifies identity 2 x 2 for variational frobenius by checking tolerance-based numerical results.
TEST(Rodin_Variational_Frobenius, Identity2x2)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  IdentityMatrix I(2);
  auto frob = Frobenius(I);

  // ||I_2||_F = sqrt(1² + 1²) = sqrt(2)
  EXPECT_NEAR(frob.getValue(p), std::sqrt(2.0), 1e-10);
}

/// @brief Verifies identity 3 x 3 for variational frobenius by checking tolerance-based numerical results.
TEST(Rodin_Variational_Frobenius, Identity3x3)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  IdentityMatrix I(3);
  auto frob = Frobenius(I);

  // ||I_3||_F = sqrt(3)
  EXPECT_NEAR(frob.getValue(p), std::sqrt(3.0), 1e-10);
}

/// @brief Verifies copy for variational frobenius by checking tolerance-based numerical results, copy semantics.
TEST(Rodin_Variational_Frobenius, Copy)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  IdentityMatrix I(2);
  auto frob = Frobenius(I);
  auto copy = frob;

  EXPECT_NEAR(copy.getValue(p), std::sqrt(2.0), 1e-10);
}

/**
 * @brief Complex spatial-vector norms remain real on every cell geometry.
 *
 * The constant field @f$z=(1,i)@f$ has
 * @f$\|z\|_2=\sqrt{2}@f$, independently of the evaluation point or cell map.
 * Its bilinear self-contraction vanishes, so it detects omission of complex
 * conjugation in the underlying norm implementation.
 */
TEST(Rodin_Variational_Frobenius, ComplexSpatialVectorAllGeometries)
{
  for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
         Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
         Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    const auto mesh = dimension == 1 ? LocalMesh::UniformGrid(geometry, {2})
      : dimension == 2               ? LocalMesh::UniformGrid(geometry, {2, 2})
                                     : LocalMesh::UniformGrid(geometry, {2, 2, 2});
    const auto cell = mesh.getPolytope(dimension, 0);
    const Point point(*cell, Polytope::Traits(geometry).getCentroid());
    const VectorFunction field(size_t(2), [](const Point&) {
      return Math::SpatialVector<Complex>({Complex(1, 0), Complex(0, 1)});
    });
    const auto norm = Frobenius(field);
    EXPECT_TRUE((std::is_same_v<decltype(norm.getValue(point)), Real>));
    EXPECT_DOUBLE_EQ(std::real(norm.getValue(point)), std::sqrt(Real(2)));
    EXPECT_DOUBLE_EQ(std::imag(norm.getValue(point)), 0);
  }
}

/**
 * @brief Scalar Frobenius norms are real magnitudes on cells and vertices.
 *
 * The scalar specialization must evaluate @f$|3+4i|=|-3-4i|=5@f$.
 * Constant and callable complex fields are checked alongside a real field,
 * including copied expressions. Cell geometry and entity dimension do not
 * alter this pointwise identity.
 */
TEST(Rodin_Variational_Frobenius, ComplexScalarAllGeometries)
{
  for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
         Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
         Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    const auto mesh = dimension == 1 ? LocalMesh::UniformGrid(geometry, {2})
      : dimension == 2               ? LocalMesh::UniformGrid(geometry, {2, 2})
                                     : LocalMesh::UniformGrid(geometry, {2, 2, 2});
    const ComplexFunction constant(Complex(3, 4));
    const ComplexFunction callable([](const Point&) { return Complex(-3, -4); });
    const RealFunction real(Real(-5));
    const auto check = [](const auto& field, const Point& point) {
      const auto norm = Frobenius(field);
      const auto copy = norm;
      EXPECT_TRUE((std::is_same_v<typename decltype(norm)::ScalarType, Real>));
      EXPECT_TRUE((std::is_same_v<decltype(norm.getValue(point)), Real>));
      EXPECT_DOUBLE_EQ(std::real(norm.getValue(point)), 5);
      EXPECT_DOUBLE_EQ(std::imag(norm.getValue(point)), 0);
      EXPECT_DOUBLE_EQ(std::real(copy(point)), 5);
      EXPECT_DOUBLE_EQ(std::imag(copy(point)), 0);
    };
    const auto cell = mesh.getPolytope(dimension, 0);
    const auto vertex = mesh.getPolytope(0, 0);
    for (const Point& point : {Point(*cell, Polytope::Traits(geometry).getCentroid()),
           Point(*vertex, Polytope::Traits(Polytope::Type::Point).getCentroid())})
    {
      check(constant, point);
      check(callable, point);
      check(real, point);
    }
  }
}

/**
 * @brief Complex matrix Frobenius norms have real values and scalar traits.
 *
 * The constant matrix @f$A=\begin{pmatrix}1&i\\3+4i&0\end{pmatrix}@f$ has squared
 * Frobenius norm @f$1+1+25=27@f$. Constant and callable coefficients,
 * including copied expressions, must retain that magnitude on cells and
 * vertices. This also checks the matrix branch of the real-valued norm
 * expression independently of the scalar and vector branches.
 */
TEST(Rodin_Variational_Frobenius, ComplexMatrixAllGeometries)
{
  Math::SpatialMatrix<Complex> value(2, 2);
  value(0, 0) = Complex(1, 0);
  value(0, 1) = Complex(0, 1);
  value(1, 0) = Complex(3, 4);
  value(1, 1) = Complex(0, 0);
  const MatrixFunction constant(value);
  const MatrixFunction callable(
    size_t(2), size_t(2), [value](const Point&) { return value; });
  for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
         Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
         Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    const auto mesh = dimension == 1 ? LocalMesh::UniformGrid(geometry, {2})
      : dimension == 2               ? LocalMesh::UniformGrid(geometry, {2, 2})
                                     : LocalMesh::UniformGrid(geometry, {2, 2, 2});
    const auto check = [](const auto& field, const Point& point) {
      const auto norm = Frobenius(field);
      const auto copy = norm;
      EXPECT_TRUE((std::is_same_v<typename decltype(norm)::ScalarType, Real>));
      EXPECT_TRUE((std::is_same_v<decltype(norm.getValue(point)), Real>));
      EXPECT_DOUBLE_EQ(std::real(norm.getValue(point)), std::sqrt(Real(27)));
      EXPECT_DOUBLE_EQ(std::imag(norm.getValue(point)), 0);
      EXPECT_DOUBLE_EQ(std::real(copy(point)), std::sqrt(Real(27)));
    };
    const auto cell = mesh.getPolytope(dimension, 0);
    const auto vertex = mesh.getPolytope(0, 0);
    for (const Point& point : {Point(*cell, Polytope::Traits(geometry).getCentroid()),
           Point(*vertex, Polytope::Traits(Polytope::Type::Point).getCentroid())})
    {
      check(constant, point);
      check(callable, point);
    }
  }
}
