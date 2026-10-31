/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TESTS_UNIT_KOKKOS_VARIATIONAL_COMMON_H
#define RODIN_TESTS_UNIT_KOKKOS_VARIATIONAL_COMMON_H

#include <cmath>

#include <gtest/gtest.h>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Kokkos.h>

#include "../TestContext.h"

namespace Rodin::Tests::Unit
{
  namespace DeviceGeometry = Rodin::Kokkos::Geometry;
  using HostMesh = Geometry::Mesh<Context::Local>;
  using DeviceMesh = Geometry::Mesh<TestContext>;
  using Polytope = Geometry::Polytope;

  /**
   * @brief Restricts an integral to attributes, keeping its type.
   * @param integral Integral to restrict.
   * @param attributes Attributes of the integrated polytopes.
   * @return Restricted integral.
   */
  template <class Form, class... Attributes>
  Form restricted(Form integral, Attributes... attributes)
  {
    integral.over(Geometry::Attribute(attributes)...);
    return integral;
  }

  /**
   * @brief Compares a device assembled operator with its host reference.
   * @param host Operator assembled on the local mesh.
   * @param device Operator assembled on the resident mesh.
   */
  inline void expectSameOperator(
    const Math::SparseMatrix<Real>& host, const Math::SparseMatrix<Real>& device)
  {
    ASSERT_EQ(host.rows(), device.rows());
    ASSERT_EQ(host.cols(), device.cols());
    ASSERT_GT(host.norm(), Real(0));
    const Math::SparseMatrix<Real> difference = host - device;
    EXPECT_LE(difference.norm(), Real(1e-12) * host.norm());
  }

  /**
   * @brief Assembles a host and a device form and compares their operators.
   * @param hu Host trial function.
   * @param hv Host test function.
   * @param h Host integral.
   * @param du Device trial function.
   * @param dv Device test function.
   * @param d Device integral of the same expression.
   */
  template <class HostTrial, class HostTest, class HostIntegral, class DeviceTrial,
    class DeviceTest, class DeviceIntegral>
  void expectSameForm(const HostTrial& hu, const HostTest& hv, const HostIntegral& h,
    const DeviceTrial& du, const DeviceTest& dv, const DeviceIntegral& d)
  {
    static_assert(
      requires(const DeviceIntegral& i) { i.assemble().values; },
      "The device expression must lower to a Kokkos integral.");
    Variational::BilinearForm hf(hu, hv);
    hf = h;
    hf.assemble();
    Variational::BilinearForm df(du, dv);
    df = d;
    df.assemble();
    expectSameOperator(hf.getOperator(), df.getOperator());
  }

  /**
   * @brief Compares a device assembled vector with its host reference.
   * @param host Vector assembled on the local mesh.
   * @param device Vector assembled on the resident mesh.
   */
  inline void expectSameVector(
    const Math::Vector<Real>& host, const Math::Vector<Real>& device)
  {
    ASSERT_EQ(host.size(), device.size());
    ASSERT_GT(host.norm(), Real(0));
    EXPECT_LE((host - device).norm(), Real(1e-12) * host.norm());
  }

  /**
   * @brief Assembles a host and a device linear form and compares their vectors.
   * @param hv Host test function.
   * @param h Host integral.
   * @param dv Device test function.
   * @param d Device integral of the same expression.
   */
  template <class HostTest, class HostIntegral, class DeviceTest, class DeviceIntegral>
  void expectSameLinearForm(const HostTest& hv, const HostIntegral& h,
    const DeviceTest& dv, const DeviceIntegral& d)
  {
    static_assert(
      requires(const DeviceIntegral& i) { i.assemble().values; },
      "The device expression must lower to a Kokkos integral.");
    Variational::LinearForm hf(hv);
    hf = h;
    hf.assemble();
    Variational::LinearForm df(dv);
    df = d;
    df.assemble();
    expectSameVector(hf.getVector(), df.getVector());
  }

  /**
   * @brief Sets smooth vertex values of a P1 grid function.
   * @param gf Grid function, component r of vertex i stored at i + r * n.
   * @param mesh Host mesh of the vertex coordinates.
   * @param vdim Vector dimension of the grid function.
   */
  template <class GF>
  void fill(GF& gf, const HostMesh& mesh, size_t vdim)
  {
    auto& data = gf.getData();
    const size_t n = mesh.getVertexCount();
    for (size_t i = 0; i < n; ++i)
    {
      const auto x = mesh.getVertexCoordinates(i);
      for (size_t r = 0; r < vdim; ++r)
      {
        data.coeffRef(i + r * n) =
          1 + Real(0.3) * r + (r + 1) * x(0) - Real(0.5) * x(1) * x(1) + x(0) * x(1);
      }
    }
  }

  /**
   * @brief Two quadrilaterals and two triangles with tagged cells.
   * @return Mixed mesh with cell attributes 1 and 2.
   */
  inline HostMesh makeMixedMesh()
  {
    auto mesh = HostMesh::Builder()
                  .initialize(2)
                  .nodes(8)
                  .vertex({0, 0})
                  .vertex({1, 0})
                  .vertex({0, 1})
                  .vertex({1.2, 1.1})
                  .vertex({2, 0})
                  .vertex({2.1, 1})
                  .vertex({0.4, 2})
                  .vertex({1.5, 2.2})
                  .polytope(Polytope::Type::Quadrilateral, {0, 1, 3, 2})
                  .polytope(Polytope::Type::Triangle, {1, 4, 3})
                  .polytope(Polytope::Type::Triangle, {4, 5, 3})
                  .polytope(Polytope::Type::Quadrilateral, {2, 3, 7, 6})
                  .finalize();
    mesh.getConnectivity().compute(2, 0);
    for (Index i = 0; i < mesh.getCellCount(); ++i)
      mesh.setAttribute({2, i}, Geometry::Attribute(1 + i % 2));
    return mesh;
  }

  /**
   * @brief Distorted triangle grid.
   * @return Uniform triangle grid with displaced interior vertices.
   */
  inline HostMesh makeTriangleMesh()
  {
    auto mesh = HostMesh::UniformGrid(Polytope::Type::Triangle, {6, 5});
    mesh.scale(Real(1) / 5);
    for (Index i = 0; i < mesh.getVertexCount(); ++i)
    {
      const auto x = mesh.getVertexCoordinates(i);
      Math::SpatialPoint y(2);
      y(0) = x(0) + Real(0.03) * std::sin(Real(3) * x(1));
      y(1) = x(1) + Real(0.02) * std::cos(Real(2) * x(0));
      mesh.setVertexCoordinates(i, y);
    }
    return mesh;
  }

  /**
   * @brief Distorted tetrahedral grid.
   * @return Uniform tetrahedral grid with displaced vertices.
   */
  inline HostMesh makeTetrahedronMesh()
  {
    auto mesh = HostMesh::UniformGrid(Polytope::Type::Tetrahedron, {3, 4, 3});
    mesh.scale(Real(1) / 3);
    for (Index i = 0; i < mesh.getVertexCount(); ++i)
    {
      const auto x = mesh.getVertexCoordinates(i);
      Math::SpatialPoint y(3);
      y(0) = x(0) + Real(0.02) * std::sin(Real(2) * x(1));
      y(1) = x(1) + Real(0.03) * std::cos(Real(3) * x(2));
      y(2) = x(2) + Real(0.02) * std::sin(Real(2) * x(0));
      mesh.setVertexCoordinates(i, y);
    }
    return mesh;
  }

  /**
   * @brief Distorted quadrilateral grid.
   * @return Uniform quadrilateral grid with displaced vertices, so the maps are bilinear.
   */
  inline HostMesh makeQuadrilateralMesh()
  {
    auto mesh = HostMesh::UniformGrid(Polytope::Type::Quadrilateral, {6, 5});
    mesh.scale(Real(1) / 5);
    for (Index i = 0; i < mesh.getVertexCount(); ++i)
    {
      const auto x = mesh.getVertexCoordinates(i);
      Math::SpatialPoint y(2);
      y(0) = x(0) + Real(0.03) * std::sin(Real(3) * x(1));
      y(1) = x(1) + Real(0.02) * std::cos(Real(2) * x(0));
      mesh.setVertexCoordinates(i, y);
    }
    return mesh;
  }

  /**
   * @brief Distorted hexahedral grid.
   * @return Uniform hexahedral grid with displaced vertices, so the maps are trilinear.
   */
  inline HostMesh makeHexahedronMesh()
  {
    auto mesh = HostMesh::UniformGrid(Polytope::Type::Hexahedron, {3, 4, 3});
    mesh.scale(Real(1) / 3);
    for (Index i = 0; i < mesh.getVertexCount(); ++i)
    {
      const auto x = mesh.getVertexCoordinates(i);
      Math::SpatialPoint y(3);
      y(0) = x(0) + Real(0.02) * std::sin(Real(2) * x(1));
      y(1) = x(1) + Real(0.03) * std::cos(Real(3) * x(2));
      y(2) = x(2) + Real(0.02) * std::sin(Real(2) * x(0));
      mesh.setVertexCoordinates(i, y);
    }
    return mesh;
  }
}

#endif
