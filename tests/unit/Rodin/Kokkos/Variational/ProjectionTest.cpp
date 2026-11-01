/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Common.h"

namespace Rodin::Tests::Unit
{
  using namespace Rodin::Variational;

  /// @brief Degree of a device H1 space.
  template <size_t K>
  constexpr std::integral_constant<size_t, K> Degree{};

  /**
   * @brief Projects on a host and a device grid function and compares them.
   * @param hfes Host space.
   * @param dfes Device space built on a copy of the host mesh.
   * @param host Projection applied to the host grid function.
   * @param device Same projection applied to the device grid function.
   *
   * Both grid functions start from the same coefficients, so that the
   * degrees of freedom outside the selection are compared too.
   */
  template <class HostFES, class DeviceFES, class HostProjection, class DeviceProjection>
  void expectSameProjection(const HostFES& hfes, const DeviceFES& dfes,
    const HostProjection& host, const DeviceProjection& device)
  {
    GridFunction hg(hfes);
    fillCoefficients(hg.getData());
    const Math::Vector<Real> initial = hg.getData();
    GridFunction dg(dfes);
    dg.getData() = initial;
    host(hg);
    device(dg);
    EXPECT_GT((hg.getData() - initial).norm(), Real(0));
    expectSameVector(hg.getData(), dg.getData());
  }

  /**
   * @brief Scalar projections match the host.
   * @param source Mesh with its closure, see withClosure().
   * @param make Builds the space on a host or a device mesh.
   */
  template <class Make>
  void checkScalarProjections(const HostMesh& source, const Make& make)
  {
    const DeviceMesh mesh(TestContext{}, source);
    const auto hfes = make(source);
    const auto dfes = make(mesh);
    P1 hsh(source);
    P1 dsh(mesh);
    GridFunction hs(hsh);
    fill(hs, source, 1);
    GridFunction ds(dsh);
    ds.getData() = hs.getData();

    const auto hf = [](const Geometry::Point& p) {
      return 1 + p.x() * p.y() - Real(0.3) * p.x() * p.y() * p.y();
    };
    const auto df = KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real
    {
      return 1 + p.x() * p.y() - Real(0.3) * p.x() * p.y() * p.y();
    };
    const RealFunction hfn(hf);
    const RealFunction dfn(df);

    expectSameProjection(
      hfes, dfes, [&](auto& g) { g = hfn; }, [&](auto& g) { g = dfn; });
    expectSameProjection(
      hfes, dfes, [](auto& g) { g = Real(2.5); }, [](auto& g) { g = Real(2.5); });
    expectSameProjection(
      hfes, dfes, [&](auto& g) { g = hfn * hs; }, [&](auto& g) { g = dfn * ds; });
    expectSameProjection(
      hfes, dfes, [&](auto& g) { g.project(Geometry::Region::Boundary, hf); },
      [&](auto& g) { g.project(Geometry::Region::Boundary, df); });
    expectSameProjection(
      hfes, dfes,
      [&](auto& g) { g.project(Geometry::Region::Faces, hfn, Geometry::Attribute(2)); },
      [&](auto& g) { g.project(Geometry::Region::Faces, dfn, Geometry::Attribute(2)); });
    expectSameProjection(
      hfes, dfes,
      [&](auto& g) {
        g.project(Geometry::Region::Interface, hs, FlatSet<Geometry::Attribute>{1, 3});
      },
      [&](auto& g) {
        g.project(Geometry::Region::Interface, ds, FlatSet<Geometry::Attribute>{1, 3});
      });
  }

  /**
   * @brief Scalar projections on tagged cells match the host.
   * @param make Builds the space on a host or a device mesh.
   */
  template <class Make>
  void checkCellAttributeProjections(const Make& make)
  {
    const auto source = withClosure(makeMixedMesh());
    const DeviceMesh mesh(TestContext{}, source);
    const auto hfes = make(source);
    const auto dfes = make(mesh);
    const RealFunction hfn(
      [](const Geometry::Point& p) { return 2 - p.x() * p.x() * p.y(); });
    const RealFunction dfn(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 2 - p.x() * p.x() * p.y();
    });
    expectSameProjection(
      hfes, dfes,
      [&](auto& g) { g.project(Geometry::Region::Cells, hfn, Geometry::Attribute(1)); },
      [&](auto& g) { g.project(Geometry::Region::Cells, dfn, Geometry::Attribute(1)); });
  }

  /**
   * @brief Vector projections match the host.
   * @param source Mesh with its closure, see withClosure().
   * @param make Builds the vector space on a host or a device mesh.
   */
  template <class Make>
  void checkVectorProjections(const HostMesh& source, const Make& make)
  {
    const DeviceMesh mesh(TestContext{}, source);
    const size_t d = source.getSpaceDimension();
    const auto hfes = make(source, d);
    const auto dfes = make(mesh, d);
    P1 hvh(source, d);
    P1 dvh(mesh, d);
    GridFunction hs(hvh);
    fill(hs, source, d);
    GridFunction ds(dvh);
    ds.getData() = hs.getData();

    const auto hf0 = [](const Geometry::Point& p) { return 1 + p.x() * p.y() * p.y(); };
    const auto hf1 = [](const Geometry::Point& p) { return p.x() - Real(0.5) * p.y(); };
    const auto df0 = KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real
    {
      return 1 + p.x() * p.y() * p.y();
    };
    const auto df1 = KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real
    {
      return p.x() - Real(0.5) * p.y();
    };
    const auto dv =
      KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->DeviceGeometry::Vector
    {
      DeviceGeometry::Vector v(d);
      for (Index r = 0; r < d; ++r)
        v(r) = (r + 1) * p.x() * p.x() - p.y();
      return v;
    };
    const auto hv = [d](const Geometry::Point& p) {
      Math::SpatialVector<Real> v(d);
      for (Index r = 0; r < d; ++r)
        v(r) = (r + 1) * p.x() * p.x() - p.y();
      return v;
    };
    Math::SpatialVector<Real> c(d);
    for (Index r = 0; r < d; ++r)
      c(r) = Real(0.5) - r;

    const auto compare = [&](const auto& hw, const auto& dw) {
      expectSameProjection(
        hfes, dfes, [&](auto& g) { g = hw; }, [&](auto& g) { g = dw; });
    };
    if (d == 2)
      compare(VectorFunction{hf0, hf1}, VectorFunction{df0, df1});
    else
      compare(VectorFunction{hf0, hf1, Real(3)}, VectorFunction{df0, df1, Real(3)});
    expectSameProjection(hfes, dfes, [&](auto& g) { g = c; }, [&](auto& g) { g = c; });
    expectSameProjection(hfes, dfes, [&](auto& g) { g = hs; }, [&](auto& g) { g = ds; });
    expectSameProjection(
      hfes, dfes, [&](auto& g) { g.project(Geometry::Region::Boundary, hv); },
      [&](auto& g) { g.project(Geometry::Region::Boundary, dv); });
  }

  /// @brief Builds a scalar P1 space, named for the device lambdas of nvcc.
  struct ScalarP1
  {
    /// @param mesh Host or device mesh. @return Space.
    template <class Mesh>
    auto operator()(const Mesh& mesh) const
    {
      return P1(mesh);
    }
  };

  /// @brief Builds a vector P1 space.
  struct VectorP1
  {
    /// @param mesh Host or device mesh. @param d Vector dimension. @return Space.
    template <class Mesh>
    auto operator()(const Mesh& mesh, size_t d) const
    {
      return P1(mesh, d);
    }
  };

  /// @brief Builds a scalar H1 space of degree K.
  template <size_t K>
  struct ScalarH1
  {
    /// @param mesh Host or device mesh. @return Space.
    template <class Mesh>
    auto operator()(const Mesh& mesh) const
    {
      return H1(Degree<K>, mesh);
    }
  };

  /// @brief Builds a vector H1 space of degree K.
  template <size_t K>
  struct VectorH1
  {
    /// @param mesh Host or device mesh. @param d Vector dimension. @return Space.
    template <class Mesh>
    auto operator()(const Mesh& mesh, size_t d) const
    {
      return H1(Degree<K>, mesh, d);
    }
  };

  TEST(Kokkos_Variational_Projection, P1OnTriangles)
  {
    checkScalarProjections(withClosure(makeTriangleMesh()), ScalarP1{});
    checkVectorProjections(withClosure(makeTriangleMesh()), VectorP1{});
  }

  TEST(Kokkos_Variational_Projection, P1OnTetrahedra)
  {
    checkScalarProjections(withClosure(makeTetrahedronMesh()), ScalarP1{});
    checkVectorProjections(withClosure(makeTetrahedronMesh()), VectorP1{});
  }

  TEST(Kokkos_Variational_Projection, P1OnMixedMesh)
  {
    checkScalarProjections(withClosure(makeMixedMesh()), ScalarP1{});
    checkVectorProjections(withClosure(makeMixedMesh()), VectorP1{});
    checkCellAttributeProjections(ScalarP1{});
  }

  TEST(Kokkos_Variational_Projection, P1OnHexahedra)
  {
    checkScalarProjections(withClosure(makeHexahedronMesh()), ScalarP1{});
    checkVectorProjections(withClosure(makeHexahedronMesh()), VectorP1{});
  }

  TEST(Kokkos_Variational_Projection, QuadraticOnTriangles)
  {
    checkScalarProjections(withClosure(makeTriangleMesh()), ScalarH1<2>{});
    checkVectorProjections(withClosure(makeTriangleMesh()), VectorH1<2>{});
  }

  TEST(Kokkos_Variational_Projection, QuadraticOnTetrahedra)
  {
    checkScalarProjections(withClosure(makeTetrahedronMesh()), ScalarH1<2>{});
    checkVectorProjections(withClosure(makeTetrahedronMesh()), VectorH1<2>{});
  }

  TEST(Kokkos_Variational_Projection, QuadraticOnMixedMesh)
  {
    checkScalarProjections(withClosure(makeMixedMesh()), ScalarH1<2>{});
    checkCellAttributeProjections(ScalarH1<2>{});
  }

  TEST(Kokkos_Variational_Projection, QuadraticOnHexahedra)
  {
    checkScalarProjections(withClosure(makeHexahedronMesh()), ScalarH1<2>{});
  }

  TEST(Kokkos_Variational_Projection, CubicOnTetrahedra)
  {
    checkScalarProjections(withClosure(makeTetrahedronMesh()), ScalarH1<3>{});
    checkVectorProjections(withClosure(makeTetrahedronMesh()), VectorH1<3>{});
  }

  TEST(Kokkos_Variational_Projection, CubicOnQuadrilaterals)
  {
    checkScalarProjections(withClosure(makeQuadrilateralMesh()), ScalarH1<3>{});
  }
}
