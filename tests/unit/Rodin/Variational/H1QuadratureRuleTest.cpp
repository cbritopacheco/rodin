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
   * @brief Qualifies vector-Jacobian gradient preparation across table lifetimes.
   *
   * Every entry is compared with the independent expression-entry quadrature
   * loop. Equal degrees permit shared physical gradients; both mixed-degree
   * directions require distinct buffers. Rebinding, quadrature replacement,
   * zero-dimensional cells and curved embedded mappings exercise the physical
   * gradient lifetime. Component counts are independent of mesh dimension.
   */
  TEST(Rodin_Variational_H1QuadratureRule, JacobianGradientTablesAndMappings)
  {
    const auto checkMesh = [](const LocalMesh& mesh) {
      const auto check = [&]<size_t TrialDegree, size_t TestDegree, class Scalar>() {
        for (const size_t components : {size_t(2), size_t(3)})
        {
          SCOPED_TRACE(TrialDegree);
          SCOPED_TRACE(TestDegree);
          SCOPED_TRACE(components);
          H1<TrialDegree, Math::SpatialVector<Scalar>> trial(
            std::integral_constant<size_t, TrialDegree>{}, mesh, components);
          H1<TestDegree, Math::SpatialVector<Scalar>> test(
            std::integral_constant<size_t, TestDegree>{}, mesh, components);
          TrialFunction u(trial);
          TestFunction v(test);
          const auto integrand = Dot(Jacobian(u), Jacobian(v));
          Integral integral(integrand);
          QuadratureReference<decltype(integral)> reference(integrand);
          // Generic and specialized loops associate dot products differently.
          constexpr Real RelativeEntryLoopTolerance = 1e-12;
          for (const size_t order : {size_t(8), size_t(2), size_t(4)})
          {
            SCOPED_TRACE(order);
            integral.setOrder(order);
            for (auto cell = mesh.getCell(); cell; ++cell)
            {
              SCOPED_TRACE(cell->getIndex());
              integral.setPolytope(*cell);
              reference.assemble(*cell, order);
              const auto& expected = reference.getOperator();
              Math::Matrix<Scalar> actual(expected.rows(), expected.cols());
              for (Eigen::Index te = 0; te < actual.rows(); ++te)
              {
                for (Eigen::Index tr = 0; tr < actual.cols(); ++tr)
                  actual(te, tr) = integral.integrate(tr, te);
              }
              const Real budget =
                RelativeEntryLoopTolerance * std::max(Real(1), expected.norm());
              ASSERT_LE((actual - expected).norm(), budget);
              auto wrong = actual;
              wrong(0, 0) += std::max(Real(1), expected.norm());
              ASSERT_GT((wrong - expected).norm(), budget);
            }
          }
        }
      };
      check.template operator()<1, 1, Real>();
      check.template operator()<2, 2, Real>();
      check.template operator()<3, 3, Real>();
      check.template operator()<1, 2, Real>();
      check.template operator()<2, 1, Real>();
      check.template operator()<2, 3, Real>();
      check.template operator()<3, 2, Real>();
      check.template operator()<1, 1, Complex>();
      check.template operator()<2, 2, Complex>();
      check.template operator()<3, 3, Complex>();
      check.template operator()<1, 2, Complex>();
      check.template operator()<2, 1, Complex>();
      check.template operator()<2, 3, Complex>();
      check.template operator()<3, 2, Complex>();
    };
    for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      const size_t dim = Polytope::Traits(geometry).getDimension();
      Array<size_t> shape(dim);
      std::fill(shape.begin(), shape.end(), 2);
      auto mesh = LocalMesh::UniformGrid(geometry, shape);
      for (size_t d = 1; d <= dim; ++d)
      {
        for (size_t lower = 0; lower < d; ++lower)
          mesh.getConnectivity().compute(d, lower);
      }
      checkMesh(mesh);
    }
    auto builder = LocalMesh::Builder();
    builder.initialize(3).nodes(3);
    builder.vertex({0, 0, 0}).vertex({1, 0, 1}).vertex({0, 1, 1});
    auto embedded = builder.polytope(Polytope::Type::Triangle, {0, 1, 2}).finalize();
    embedded.getConnectivity().compute(2, 1);
    checkMesh(embedded);
    RealH1Element<2> geometryElement(Polytope::Type::Triangle);
    PointCloud nodes(3, geometryElement.getCount());
    for (size_t a = 0; a < geometryElement.getCount(); ++a)
    {
      const auto& r = geometryElement.getNode(a);
      nodes(0, a) = r[0];
      nodes(1, a) = r[1];
      nodes(2, a) = r[0] + r[1] + Real(0.25) * r[0] * r[1];
    }
    embedded.setPolytopeTransformation(
      {2, 0}, new ParametricTransformation<RealH1Element<2>>(nodes, geometryElement));
    checkMesh(embedded);
    auto pointBuilder = LocalMesh::Builder();
    pointBuilder.initialize(3).nodes(1).vertex({0, 0, 0});
    checkMesh(pointBuilder.finalize());
  }

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
