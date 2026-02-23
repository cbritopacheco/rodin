/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Variational.h"
#include "Rodin/Eikonal/FMM.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Manufactured::Eikonal
{
  template <size_t NX, size_t NY, size_t NZ = 1>
  class EikonalManufacturedTest : public ::testing::TestWithParam<Polytope::Type>
  {
  protected:
    Mesh<Context::Local> buildMesh() const
    {
      const auto geom = GetParam();
      if constexpr (NZ == 1)
      {
        Mesh<Context::Local> mesh = Mesh<Context::Local>::UniformGrid(geom, { NX, NY });
        mesh.scale(1.0 / (NX - 1));
        mesh.getConnectivity().compute(1, 2);
        mesh.getConnectivity().compute(0, 1);
        mesh.getConnectivity().compute(0, 0);
        return mesh;
      }
      else
      {
        Mesh<Context::Local> mesh = Mesh<Context::Local>::UniformGrid(geom, { NX, NY, NZ });
        mesh.scale(1.0 / (NX - 1));
        mesh.getConnectivity().compute(2, 3);
        mesh.getConnectivity().compute(0, 1);
        mesh.getConnectivity().compute(0, 0);
        return mesh;
      }
    }

    template <class GF>
    static void checkDistance(const Mesh<Context::Local>& mesh, const GF& u, Real tol)
    {
      const size_t dim = mesh.getDimension();
      Math::SpatialVector<Real> center(dim == 2 ? Math::SpatialVector<Real>{{0.5, 0.5}}
                                                : Math::SpatialVector<Real>{{0.5, 0.5, 0.5}});

      for (auto it = mesh.getVertex(); !it.end(); ++it)
      {
        const auto coord = mesh.getVertexCoordinates(it->getIndex());
        const Real r = (coord - center).norm();
        const Real err = std::abs(u[it->getIndex()] - r);
        EXPECT_LT(err, tol);
      }
    }
  };

  using Eikonal2DTest = EikonalManufacturedTest<33, 33, 1>;
  using Eikonal3DTest = EikonalManufacturedTest<33, 33, 33>;

  TEST_P(Eikonal2DTest, ConstantSpeed_PointSource)
  {
    auto mesh = buildMesh();
    P1 vh(mesh);
    GridFunction u(vh);

    auto speed = [](const Geometry::Point&) -> Real { return 1.0; };

    // seed at center
    std::vector<Index> interface;
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      if ((coord - Math::SpatialVector<Real>{{0.5, 0.5}}).norm() < 1e-6)
      {
        interface.push_back(it->getIndex());
        break;
      }
    }

    ASSERT_FALSE(interface.empty());

    Rodin::Eikonal::FMM fmm(u, speed);
    fmm.seed(interface).solve();

    const Real tol = 5e-2;
    checkDistance(mesh, u, tol);
  }

  TEST_P(Eikonal3DTest, ConstantSpeed_PointSource)
  {
    auto mesh = buildMesh();
    P1 vh(mesh);

    GridFunction u(vh);

    auto speed = [](const Geometry::Point&) -> Real { return 1.0; };

    std::vector<Index> interface;
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto coord = mesh.getVertexCoordinates(it->getIndex());
      if ((coord - Math::SpatialVector<Real>{{0.5, 0.5, 0.5}}).norm() < 1e-5)
      {
        interface.push_back(it->getIndex());
        break;
      }
    }
    ASSERT_FALSE(interface.empty());

    Rodin::Eikonal::FMM fmm(u, speed);

    fmm.seed(interface);
    fmm.solve();

    const Real tol = 9e-2;
    checkDistance(mesh, u, tol);
  }

  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage2D,
    Eikonal2DTest,
    ::testing::Values(
      Polytope::Type::Triangle)
  );

  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage3D,
    Eikonal3DTest,
    ::testing::Values(
      Polytope::Type::Tetrahedron)
  );
}
