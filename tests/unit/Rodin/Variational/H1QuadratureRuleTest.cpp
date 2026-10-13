#include <gtest/gtest.h>

#include "Rodin/Variational.h"
#include "Rodin/Assembly/Default.h"
#include "Rodin/Geometry/Mesh.h"
#include "../../../QuadratureReference.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  /**
   * @brief Distinguishes shared and independent physical-gradient preparation.
   *
   * Every local entry equals the original quadrature entry loop on affine
   * volume cells. Equal-degree rules share reference tables, while both mixed
   * degree directions must retain distinct gradient buffers. Rebinding across
   * cells and changing orders must not retain the previous physical gradients.
   */
  TEST(Rodin_Variational_H1QuadratureRule, SharedAndDistinctGradientTables)
  {
    const auto check = []<size_t TrialDegree, size_t TestDegree, class Scalar>() {
      for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
             Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
             Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
      {
        SCOPED_TRACE(static_cast<int>(geometry));
        SCOPED_TRACE(TrialDegree);
        SCOPED_TRACE(TestDegree);
        const size_t dimension = Polytope::Traits(geometry).getDimension();
        Array<size_t> shape(dimension);
        std::fill(shape.begin(), shape.end(), 2);
        auto mesh = LocalMesh::UniformGrid(geometry, shape);
        for (size_t d = 1; d <= dimension; ++d)
        {
          for (size_t lower = 0; lower < d; ++lower)
            mesh.getConnectivity().compute(d, lower);
        }
        H1<TrialDegree, Scalar> trial(
          std::integral_constant<size_t, TrialDegree>{}, mesh);
        H1<TestDegree, Scalar> test(std::integral_constant<size_t, TestDegree>{}, mesh);
        TrialFunction u(trial);
        TestFunction v(test);
        auto integral = Integral(Grad(u), Grad(v));
        QuadratureReference reference(integral);
        for (const size_t order : {size_t(8), size_t(2), size_t(4)})
        {
          SCOPED_TRACE(order);
          integral.setOrder(order);
          for (auto cell = mesh.getCell(); cell; ++cell)
          {
            integral.setPolytope(*cell);
            reference.assemble(*cell, order);
            const auto& expected = reference.getOperator();
            for (Eigen::Index te = 0; te < expected.rows(); ++te)
            {
              for (Eigen::Index tr = 0; tr < expected.cols(); ++tr)
                ASSERT_EQ(integral.integrate(tr, te), expected(te, tr));
            }
          }
        }
      }
    };
    check.template operator()<1, 1, Real>();
    check.template operator()<2, 2, Real>();
    check.template operator()<1, 2, Real>();
    check.template operator()<2, 1, Real>();
    check.template operator()<1, 1, Complex>();
    check.template operator()<2, 2, Complex>();
    check.template operator()<1, 2, Complex>();
    check.template operator()<2, 1, Complex>();
  }

  /**
   * @brief Preserves every component of scalar-weighted H1 mass on rebinding.
   *
   * A dyadic coefficient permits exact comparisons on the directly accumulated
   * lower triangle. The unweighted kernel mirrors its upper triangle, whereas
   * the weighted rule accumulates it independently; roundoff equality between
   * these different upper-triangle algorithms is not required. Every diagonal
   * must be positive, so dropping
   * the same component in both kernels cannot satisfy the scaling identity.
   * Component counts deliberately differ from the geometric dimension.
   * Orders change while the same rules traverse all cells, exercising retained
   * storage and formula replacement for scalar, vector and matrix fields.
   */
  TEST(Rodin_Variational_H1QuadratureRule, MassComponentContractAndRebinding)
  {
    const auto checkDegree = []<size_t K>() {
      for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
             Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
             Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
      {
        SCOPED_TRACE(static_cast<int>(geometry));
        SCOPED_TRACE(K);
        const size_t dimension = Polytope::Traits(geometry).getDimension();
        Array<size_t> shape(dimension);
        std::fill(shape.begin(), shape.end(), 2);
        auto mesh = LocalMesh::UniformGrid(geometry, shape);
        for (size_t d = 1; d <= dimension; ++d)
        {
          for (size_t lower = 0; lower < d; ++lower)
            mesh.getConnectivity().compute(d, lower);
        }
        const auto checkRange = [&]<class Range>() {
          const auto space = [&]() {
            if constexpr (FormLanguage::IsMatrixRange<Range>::Value)
              return H1<K, Range>(std::integral_constant<size_t, K>{}, mesh, 2, 3);
            else if constexpr (FormLanguage::IsVectorRange<Range>::Value)
              return H1<K, Range>(std::integral_constant<size_t, K>{}, mesh, 3);
            else
              return H1<K, Range>(std::integral_constant<size_t, K>{}, mesh);
          }();
          TrialFunction u(space);
          TestFunction v(space);
          auto mass = Integral(u, v);
          auto weighted = Integral(Real(0.25) * u, v);
          for (const size_t order : {size_t(8), size_t(2), size_t(4)})
          {
            SCOPED_TRACE(order);
            mass.setOrder(order);
            weighted.setOrder(order);
            for (auto cell = mesh.getCell(); cell; ++cell)
            {
              mass.setPolytope(*cell);
              weighted.setPolytope(*cell);
              const size_t count = u.getDOFs(*cell);
              for (size_t te = 0; te < count; ++te)
              {
                ASSERT_GT(std::real(mass.integrate(te, te)), 0);
                ASSERT_GT(std::real(weighted.integrate(te, te)), 0);
                for (size_t tr = 0; tr <= te; ++tr)
                {
                  EXPECT_EQ(
                    weighted.integrate(tr, te), Real(0.25) * mass.integrate(tr, te));
                }
              }
            }
          }
        };
        checkRange.template operator()<Real>();
        checkRange.template operator()<Complex>();
        checkRange.template operator()<Math::SpatialVector<Real>>();
        checkRange.template operator()<Math::SpatialVector<Complex>>();
        checkRange.template operator()<Math::SpatialMatrix<Real>>();
        checkRange.template operator()<Math::SpatialMatrix<Complex>>();
      }
    };
    checkDegree.template operator()<1>();
    checkDegree.template operator()<2>();
  }

  /// @brief Verifies mixed order grad grad assembles for variational H1 quadrature rule by checking exact expected values, form assembly.
  TEST(Rodin_Variational_H1QuadratureRule, MixedOrder_GradGrad_Assembles)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    mesh.getConnectivity().compute(1, 0);
    H1<2, Real> fesTr(std::integral_constant<size_t, 2>{}, mesh);
    H1<1, Real> fesTe(std::integral_constant<size_t, 1>{}, mesh);

    TrialFunction u(fesTr);
    TestFunction v(fesTe);

    BilinearForm bf(u, v);
    bf = Integral(Grad(u), Grad(v));
    bf.assemble();

    const auto& mat = bf.getOperator();
    EXPECT_EQ(mat.rows(), static_cast<Eigen::Index>(fesTe.getSize()));
    EXPECT_EQ(mat.cols(), static_cast<Eigen::Index>(fesTr.getSize()));
    EXPECT_GT(mat.norm(), 0.0);
  }

  /// @brief Verifies mixed order jacobian assembles for variational H1 quadrature rule by checking exact expected values, form assembly.
  TEST(Rodin_Variational_H1QuadratureRule, MixedOrder_Jacobian_Assembles)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    mesh.getConnectivity().compute(1, 0);
    const size_t vdim = 2;
    H1<2, Math::SpatialVector<Real>> fesTr(std::integral_constant<size_t, 2>{}, mesh, vdim);
    H1<1, Math::SpatialVector<Real>> fesTe(std::integral_constant<size_t, 1>{}, mesh, vdim);

    TrialFunction u(fesTr);
    TestFunction v(fesTe);

    BilinearForm bf(u, v);
    bf = Integral(Jacobian(u), Jacobian(v));
    bf.assemble();

    const auto& mat = bf.getOperator();
    EXPECT_EQ(mat.rows(), static_cast<Eigen::Index>(fesTe.getSize()));
    EXPECT_EQ(mat.cols(), static_cast<Eigen::Index>(fesTr.getSize()));
    EXPECT_GT(mat.norm(), 0.0);
  }

  /// @brief Verifies vector mass matches grid function linear form for variational H1 quadrature rule by checking tolerance-based numerical results, exact expected values, form assembly.
  TEST(Rodin_Variational_H1QuadratureRule, VectorMass_MatchesGridFunctionLinearForm)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    mesh.getConnectivity().compute(1, 0);

    const size_t vdim = 2;
    H1<2, Math::SpatialVector<Real>> fes(std::integral_constant<size_t, 2>{}, mesh, vdim);
    TrialFunction u(fes);
    TestFunction v(fes);

    GridFunction gf(fes);
    auto& x = gf.getData();
    ASSERT_EQ(x.size(), static_cast<Eigen::Index>(fes.getSize()));
    for (Eigen::Index i = 0; i < x.size(); ++i)
      x(i) = static_cast<Real>(i + 1);

    BilinearForm mass(u, v);
    mass = Integral(Dot(u, v));
    mass.assemble();

    LinearForm load(v);
    load = Integral(Dot(gf, v));
    load.assemble();

    const auto residual = mass.getOperator() * x - load.getVector();
    EXPECT_NEAR(residual.norm(), 0.0, 1e-13);
  }
}
