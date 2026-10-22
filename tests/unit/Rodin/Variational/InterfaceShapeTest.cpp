/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include <gtest/gtest.h>
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"
#include "Rodin/Variational/H1.h"
#include "Rodin/Assembly/Sequential.h"
#ifdef RODIN_USE_OPENMP
#include "Rodin/Assembly/OpenMP.h"
#endif

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

/// @brief A global P1 basis uses the appropriate local index on each incident cell.
TEST(InterfaceShape, ContinuousValues)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  mesh.getConnectivity().compute(1, 2);
  P1<Real> fes(mesh);
  TrialFunction u(fes);
  auto average = Average(u);
  auto jump = Jump(u);
  for (auto face = mesh.getFace(); face; ++face)
  {
    if (mesh.getConnectivity().getIncidence({1, 2}, face->getIndex()).size() != 2)
      continue;
    const Point point(*face, Polytope::Traits(face->getGeometry()).getCentroid());
    const IntegrationPoint ip(point);
    average.setPolytope(*face).setIntegrationPoint(ip);
    jump.setPolytope(*face).setIntegrationPoint(ip);
    const auto& indices = *average.getDOFIndices();
    Real sum = 0;
    for (size_t local = 0; local < static_cast<size_t>(indices.size()); ++local)
    {
      EXPECT_NEAR(jump.getBasis(local), 0, 1e-13);
      sum += average.getBasis(local);
    }
    EXPECT_NEAR(sum, 1, 1e-13);
  }
}

/// @brief Discontinuous constants retain both cell DOFs with opposite jump signs.
TEST(InterfaceShape, DiscontinuousConstants)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  mesh.getConnectivity().compute(1, 2);
  P0<Real> fes(mesh);
  TrialFunction u(fes);
  auto average = Average(u);
  auto jump = Jump(u);
  for (auto face = mesh.getFace(); face; ++face)
  {
    if (mesh.getConnectivity().getIncidence({1, 2}, face->getIndex()).size() != 2)
      continue;
    const Point point(*face, Polytope::Traits(face->getGeometry()).getCentroid());
    const IntegrationPoint ip(point);
    average.setPolytope(*face).setIntegrationPoint(ip);
    jump.setPolytope(*face).setIntegrationPoint(ip);
    ASSERT_EQ(average.getDOFIndices()->size(), 2);
    EXPECT_EQ(average.getBasis(0), 0.5);
    EXPECT_EQ(average.getBasis(1), 0.5);
    EXPECT_EQ(jump.getBasis(0) + jump.getBasis(1), 0);
    EXPECT_EQ(std::abs(jump.getBasis(0)), 1);
  }
}

/// @brief Interface assembly scatters both discontinuous cell supports and wrapper maps.
TEST(InterfaceShape, AssembledConstants)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  mesh.getConnectivity().compute(1, 2);
  P0<Real> fes(mesh);
  TrialFunction u(fes);
  TestFunction v(fes);
  Math::Matrix<Real> expected(fes.getSize(), fes.getSize());
  Math::Vector<Real> load(fes.getSize());
  expected.setZero();
  load.setZero();
  for (auto face = mesh.getFace(); face; ++face)
  {
    const auto& incident = mesh.getConnectivity().getIncidence({1, 2}, face->getIndex());
    if (incident.size() != 2)
      continue;
    const Index first = *incident.begin();
    const Index second = *std::next(incident.begin());
    const Real measure = face->getMeasure();
    expected(first, first) += measure;
    expected(second, second) += measure;
    expected(first, second) -= measure;
    expected(second, first) -= measure;
    load(first) += measure / 2;
    load(second) += measure / 2;
  }
  BilinearForm form(u, v);
  form = -(-InterfaceIntegral(Jump(u), Jump(v)));
  form.assemble();
  const Math::Matrix<Real> actual = form.getOperator();
  EXPECT_LT((actual - expected).norm(), 1e-12);
  const auto checkBackend = [&]<template <class, class> class Backend>() {
    Math::SparseMatrix<Real> matrix;
    Backend<Math::SparseMatrix<Real>, decltype(form)> backend;
    if constexpr (requires { backend.setThreadCount(size_t(2)); })
      backend.setThreadCount(2);
    backend.execute(
      matrix, {fes, fes, form.getLocalIntegrators(), form.getGlobalIntegrators()});
    EXPECT_LT((Math::Matrix<Real>(matrix) - expected).norm(), 1e-12);
  };
  checkBackend.template operator()<Assembly::Sequential>();
#ifdef RODIN_USE_OPENMP
  checkBackend.template operator()<Assembly::OpenMP>();
#endif
  LinearForm linear(v);
  linear = -(-2 * InterfaceIntegral(Average(v) / RealFunction(2)));
  linear.assemble();
  EXPECT_LT((linear.getVector() - load).norm(), 1e-12);
}

/// @brief Composed averages and jumps preserve the union ordering after copies and moves.
TEST(InterfaceShape, ComposedValues)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  mesh.getConnectivity().compute(1, 2);
  P1<Real> fes(mesh);
  TrialFunction u(fes);
  auto expression = 2 * Average(u) + Jump(u);
  auto copy = expression;
  auto moved = std::move(copy);
  for (auto face = mesh.getFace(); face; ++face)
  {
    if (mesh.getConnectivity().getIncidence({1, 2}, face->getIndex()).size() != 2)
      continue;
    const Point point(*face, Polytope::Traits(face->getGeometry()).getCentroid());
    const IntegrationPoint ip(point);
    moved.setPolytope(*face).setIntegrationPoint(ip);
    Real sum = 0;
    for (size_t local = 0; local < static_cast<size_t>(moved.getDOFIndices()->size());
         ++local)
      sum += moved.getBasis(local);
    EXPECT_NEAR(sum, 2, 1e-13);
  }
}

namespace
{
  /**
   * @brief Compares two-cell gradients with individual finite-element derivatives.
   * @tparam Degree Polynomial degree; zero selects native P1.
   * @param geometry Mesh cell geometry.
   */
  template <size_t Degree>
  void checkInterfaceGradients(Polytope::Type geometry)
  {
    auto mesh = LocalMesh::UniformGrid(geometry, {2, 2});
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    mesh.getConnectivity().compute(1, 2);
    auto fes = [&] {
      if constexpr (Degree == 0)
        return P1<Real>(mesh);
      else
        return H1<Degree, Real>(std::integral_constant<size_t, Degree>{}, mesh);
    }();
    TrialFunction u(fes);
    auto average = Average(Grad(u));
    auto jump = Jump(Grad(u));
    for (auto face = mesh.getFace(); face; ++face)
    {
      const auto& incident =
        mesh.getConnectivity().getIncidence({1, 2}, face->getIndex());
      if (incident.size() != 2)
        continue;
      const Point point(*face, Polytope::Traits(face->getGeometry()).getCentroid());
      const IntegrationPoint ip(point);
      average.setPolytope(*face).setIntegrationPoint(ip);
      jump.setPolytope(*face).setIntegrationPoint(ip);
      const IndexArray indices = *average.getDOFIndices();
      std::vector<Math::SpatialVector<Real>> averages(
        indices.size(), Math::SpatialVector<Real>::Zero(2));
      auto jumps = averages;
      size_t side = 0;
      for (const Index index : incident)
      {
        const auto cell = mesh.getPolytope(2, index);
        Math::SpatialPoint reference;
        cell->getTransformation().inverse(reference, point.getPhysicalCoordinates());
        const Point cellPoint(*cell, reference);
        const auto& element = fes.getFiniteElement(2, index);
        const IndexArray dofs = fes.getDOFs(2, index);
        for (size_t local = 0; local < element.getCount(); ++local)
        {
          Math::SpatialVector<Real> derivative(2);
          for (size_t coordinate = 0; coordinate < 2; ++coordinate)
          {
            derivative[coordinate] =
              element.getBasis(local).template getDerivative<1>(coordinate)(reference);
          }
          const auto gradient = cellPoint.getJacobianInverse().transpose() * derivative;
          const auto found = std::find(indices.begin(), indices.end(), dofs[local]);
          ASSERT_NE(found, indices.end());
          const size_t position = std::distance(indices.begin(), found);
          averages[position] += gradient * 0.5;
          jumps[position] += gradient * (side == 0 ? 1 : -1);
        }
        ++side;
      }
      for (size_t local = 0; local < static_cast<size_t>(indices.size()); ++local)
      {
        EXPECT_LT((average.getBasis(local) - averages[local]).norm(), 1e-11);
        EXPECT_LT((jump.getBasis(local) - jumps[local]).norm(), 1e-11);
      }
    }
  }
}

/// @brief Native P1 gradients retain the opposite vertex contribution on a face.
TEST(InterfaceShape, P1Gradients)
{
  checkInterfaceGradients<0>(Polytope::Type::Triangle);
}

/// @brief Higher-order face gradients include cell modes with zero value trace.
TEST(InterfaceShape, CubicGradients)
{
  checkInterfaceGradients<3>(Polytope::Type::Triangle);
}

/// @brief Tensor-product cell bases use independent side-local coordinates.
TEST(InterfaceShape, QuadrilateralGradients)
{
  checkInterfaceGradients<2>(Polytope::Type::Quadrilateral);
}

/// @brief Binding rejects a cell, a boundary face, and a face from a foreign mesh.
TEST(InterfaceShape, InvalidBindings)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  mesh.getConnectivity().compute(1, 2);
  auto other = mesh;
  P1<Real> fes(mesh);
  TrialFunction u(fes);
  auto average = Average(u);
  const Point cellPoint(*mesh.getCell(), Math::SpatialPoint{0.2, 0.3});
  const IntegrationPoint cellIP(cellPoint);
  EXPECT_THROW(average.setIntegrationPoint(cellIP), Alert::Exception);
  EXPECT_THROW(average.setPolytope(*mesh.getCell()), Alert::Exception);
  EXPECT_THROW(average.setPolytope(*other.getFace()), Alert::Exception);
  for (auto face = mesh.getFace(); face; ++face)
  {
    if (mesh.getConnectivity().getIncidence({1, 2}, face->getIndex()).size() == 1)
    {
      EXPECT_THROW(average.setPolytope(*face), Alert::Exception);
      break;
    }
  }
}

/// @brief Matrix-valued complex traces preserve component shape and continuity.
TEST(InterfaceShape, ComplexMatrixValues)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);
  mesh.getConnectivity().compute(1, 2);
  H1<2, Math::SpatialMatrix<Complex>> fes(
    std::integral_constant<size_t, 2>{}, mesh, 2, 3);
  TrialFunction u(fes);
  auto average = Average(u);
  auto jump = Jump(u);
  for (auto face = mesh.getFace(); face; ++face)
  {
    if (mesh.getConnectivity().getIncidence({1, 2}, face->getIndex()).size() != 2)
      continue;
    const Point point(*face, Polytope::Traits(face->getGeometry()).getCentroid());
    const IntegrationPoint ip(point);
    average.setPolytope(*face).setIntegrationPoint(ip);
    jump.setPolytope(*face).setIntegrationPoint(ip);
    for (size_t local = 0; local < static_cast<size_t>(average.getDOFIndices()->size());
         ++local)
    {
      const auto value = average.getBasis(local);
      EXPECT_EQ(value.rows(), 2);
      EXPECT_EQ(value.cols(), 3);
      EXPECT_LT(jump.getBasis(local).norm(), 1e-12);
    }
  }
}

/// @brief A triangle/quad interface combines different cell supports with a face trace.
TEST(InterfaceShape, MixedGeometryAndFaceSupport)
{
  auto mesh = LocalMesh::Builder()
                .initialize(2)
                .nodes(5)
                .vertex({0, 0})
                .vertex({1, 0})
                .vertex({1, 1})
                .vertex({0, 1})
                .vertex({2, 0})
                .polytope(Polytope::Type::Quadrilateral, {0, 1, 2, 3})
                .polytope(Polytope::Type::Triangle, {1, 4, 2})
                .finalize();
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 2);
  P1<Real> fes(mesh);
  TrialFunction u(fes);
  auto average = Average(u);
  auto jump = Jump(u);
  auto combined = Average(u) + u;
  for (auto face = mesh.getFace(); face; ++face)
  {
    if (mesh.getConnectivity().getIncidence({1, 2}, face->getIndex()).size() != 2)
      continue;
    const Point point(*face, Math::SpatialPoint{0.3});
    const IntegrationPoint ip(point);
    average.setPolytope(*face).setIntegrationPoint(ip);
    jump.setPolytope(*face).setIntegrationPoint(ip);
    combined.setPolytope(*face).setIntegrationPoint(ip);
    ASSERT_EQ(average.getDOFIndices()->size(), 5);
    Real total = 0;
    for (size_t local = 0; local < 5; ++local)
    {
      EXPECT_NEAR(jump.getBasis(local), 0, 1e-13);
      total += combined.getBasis(local);
    }
    EXPECT_NEAR(total, 2, 1e-13);
  }
}

/// @brief Union numbering accepts empty supports, rejects duplicates, and can be rebound.
TEST(InterfaceShape, UnionNumbering)
{
  DOFMapping<> mapping;
  mapping.setLHS(IndexArray{});
  mapping.setRHS(IndexArray{});
  EXPECT_EQ(mapping.getIndices().size(), 0);
  mapping.setLHS(IndexArray{{1, 1}});
  EXPECT_THROW(mapping.setRHS(IndexArray{{2}}), Alert::Exception);
  mapping.setLHS(IndexArray{{3, 1}});
  mapping.setRHS(IndexArray{{2, 3}});
  ASSERT_EQ(mapping.getIndices().size(), 3);
  EXPECT_EQ(mapping.getIndices()[0], 1);
  EXPECT_EQ(mapping.getIndices()[1], 2);
  EXPECT_EQ(mapping.getIndices()[2], 3);
  EXPECT_EQ(mapping.getLHS(0), Optional<size_t>(1));
  EXPECT_FALSE(mapping.getRHS(0));
  EXPECT_FALSE(mapping.getLHS(1));
  EXPECT_EQ(mapping.getRHS(1), Optional<size_t>(0));
  EXPECT_EQ(mapping.getLHS(2), Optional<size_t>(0));
  EXPECT_EQ(mapping.getRHS(2), Optional<size_t>(1));
  mapping.setLHS(IndexArray{{1}});
  EXPECT_THROW(mapping.setRHS(IndexArray{{2, 2}}), Alert::Exception);
}

/// @brief Complex scaling of real integrators clones their public leaf and attribute contracts.
TEST(InterfaceShape, PromotedIntegratorCopies)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  mesh.getConnectivity().compute(1, 2);
  P0<Real> fes(mesh);
  TrialFunction u(fes);
  TestFunction v(fes);
  auto matrix = Complex(0, 2) * InterfaceIntegral(Jump(u), Jump(v));
  auto vector = Complex(0, 2) * InterfaceIntegral(Average(v));
  matrix.over(7);
  vector.over(7);
  auto matrixCopy = matrix;
  auto vectorCopy = vector;
  EXPECT_EQ(matrixCopy.getAttributes().count(7), 1);
  EXPECT_EQ(vectorCopy.getAttributes().count(7), 1);
  for (auto face = mesh.getFace(); face; ++face)
  {
    if (mesh.getConnectivity().getIncidence({1, 2}, face->getIndex()).size() != 2)
      continue;
    matrixCopy.setPolytope(*face);
    vectorCopy.setPolytope(*face);
    EXPECT_EQ(matrixCopy.getTrialDOFIndices()->size(), 2);
    EXPECT_EQ(vectorCopy.getTestDOFIndices()->size(), 2);
    EXPECT_NEAR(std::abs(matrixCopy.integrate(0, 0) - Complex(0, 2 * face->getMeasure())),
      0, 1e-13);
    EXPECT_NEAR(
      std::abs(vectorCopy.integrate(0) - Complex(0, face->getMeasure())), 0, 1e-13);
  }
}
