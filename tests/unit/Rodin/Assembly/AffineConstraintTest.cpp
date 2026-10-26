/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Independent adjoint-projection oracles for complex affine identification. */
#include <gtest/gtest.h>
#include <algorithm>
#include <array>
#include <vector>
#include "Rodin/Variational.h"
#include "Rodin/Assembly/Sequential.h"
#ifdef RODIN_USE_OPENMP
#include "Rodin/Assembly/OpenMP.h"
#endif

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  template <class Scalar>
  class AffineConstraints : public ::testing::Test
  {};
  using AffineScalars = ::testing::Types<Real, Complex>;
  TYPED_TEST_SUITE(AffineConstraints, AffineScalars);

  /// @brief Native dense/sparse backends obey the independent sesquilinear projection for one and two fields.
  template <class Scalar, class Operator, template <class, class> class Backend>
  void checkAffineConstraints()
  {
    const Scalar phase = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(1, -0.5);
      else
        return Scalar(1);
    }();
    const Scalar coefficient = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(0.25, 0.5);
      else
        return Scalar(0.25);
    }();
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    mesh.scale(0.5);
    mesh.getConnectivity().compute(0, 1);
    std::vector<std::pair<Real, Index>> vertices;
    for (auto vertex = mesh.getVertex(); vertex; ++vertex)
    {
      const Real x = mesh.getVertexCoordinates(vertex->getIndex()).x();
      vertices.emplace_back(x, vertex->getIndex());
      if (x == 0)
        mesh.setAttribute(vertex.key(), 11);
      if (x == 1)
        mesh.setAttribute(vertex.key(), 12);
    }
    std::sort(vertices.begin(), vertices.end());
    ASSERT_EQ(vertices.size(), 3u);
    using System = Math::LinearSystem<Operator, Math::Vector<Scalar>>;
    using Mapping = std::array<std::array<Index, 3>, 2>;
    auto verify = [&](auto& assembler, System& system, const auto& input,
                    const Mapping& indices) {
      if constexpr (requires { assembler.setThreadCount(size_t(2)); })
        assembler.setThreadCount(2);
      Math::Matrix<Scalar> raw = Math::Matrix<Scalar>::Zero(6, 6);
      Math::Vector<Scalar> load(6), exact(6);
      const Real mass[3][3] = {{Real(1) / 6, Real(1) / 12, 0},
        {Real(1) / 12, Real(1) / 3, Real(1) / 12}, {0, Real(1) / 12, Real(1) / 6}};
      const Real moments[3] = {Real(1) / 4, Real(1) / 2, Real(1) / 4};
      for (size_t component = 0; component < 2; ++component)
      {
        for (size_t node = 0; node < 3; ++node)
        {
          exact[indices[component][node]] = Scalar(component + 1) * phase;
          load[indices[component][node]] = Scalar(component + 1) * phase * moments[node];
          for (size_t other = 0; other < 3; ++other)
            raw(indices[component][node], indices[component][other]) = mass[node][other];
        }
      }
      ASSERT_LT((raw * exact - load).norm(), 1e-13);
      const Index slave = indices[0][0], master = indices[1][0];
      const std::array<Index, 2> fixed{indices[0][2], indices[1][2]};
      Math::Matrix<Scalar> projection = Math::Matrix<Scalar>::Identity(6, 6);
      projection.row(slave).setZero();
      projection(slave, master) = coefficient;
      Math::Vector<Scalar> lift = Math::Vector<Scalar>::Zero(6);
      lift[slave] = (Scalar(1) - Scalar(2) * coefficient) * phase;
      for (const auto dof : fixed)
      {
        projection.row(dof).setZero();
        lift[dof] = exact[dof];
      }
      Math::Matrix<Scalar> expected = projection.adjoint() * raw * projection;
      Math::Vector<Scalar> expectedLoad = projection.adjoint() * (load - raw * lift);
      expected(slave, slave) = 1;
      expected(slave, master) = -coefficient;
      expectedLoad[slave] = lift[slave];
      for (const auto dof : fixed)
      {
        expected(dof, dof) = 1;
        expectedLoad[dof] = exact[dof];
      }
      for (size_t repeat = 0; repeat < 2; ++repeat)
      {
        assembler.execute(system, input);
        const Math::Matrix<Scalar> actual = system.getOperator();
        EXPECT_LT((actual - expected).norm(), 1e-12);
        EXPECT_LT((system.getVector() - expectedLoad).norm(), 1e-12);
        const Math::Vector<Scalar> recovered =
          actual.fullPivLu().solve(system.getVector());
        EXPECT_LT((recovered - exact).norm(), 1e-11);
      }
    };
    {
      SCOPED_TRACE("single vector field");
      P1<Math::SpatialVector<Scalar>> space(mesh, 2);
      TrialFunction u(space);
      TestFunction v(space);
      using ProblemType = Problem<System, decltype(u), decltype(v)>;
      using Body = typename ProblemType::ProblemBodyType;
      Math::SpatialMatrix<Scalar> transform(2, 2);
      transform(0, 0) = 0;
      transform(0, 1) = coefficient;
      transform(1, 0) = 0;
      transform(1, 1) = 1;
      const auto forcing = VectorFunction(size_t(2), [phase](const Point&) {
        return Math::SpatialVector<Scalar>{phase, Scalar(2) * phase};
      });
      const auto offset = VectorFunction(size_t(2), [phase, coefficient](const Point&) {
        return Math::SpatialVector<Scalar>{
          (Scalar(1) - Scalar(2) * coefficient) * phase, Scalar(0)};
      });
      Body body = Integral(u, v) - Integral(forcing, v) +
        DirichletBC(u, MatrixFunction(transform) * u, offset).on(11) +
        DirichletBC(u, forcing).on(12);
      Assembly::ProblemAssemblyInput<Body, decltype(u), decltype(v)> input(body, u, v);
      Backend<System, ProblemType> assembler;
      System system;
      Mapping indices;
      for (size_t node = 0; node < 3; ++node)
      {
        for (size_t component = 0; component < 2; ++component)
          indices[component][node] = space.getDOFs(0, vertices[node].second)[component];
      }
      verify(assembler, system, input, indices);
    }
    {
      SCOPED_TRACE("two scalar fields");
      P1<Scalar> space(mesh);
      TrialFunction u(space), p(space);
      TestFunction v(space), q(space);
      using ProblemType =
        Problem<System, decltype(u), decltype(v), decltype(p), decltype(q)>;
      using Body = typename ProblemType::ProblemBodyType;
      const auto known = [](Scalar value) {
        if constexpr (std::is_same_v<Scalar, Complex>)
          return ComplexFunction(value);
        else
          return RealFunction(value);
      };
      Body body = Integral(u, v) + Integral(p, q) - Integral(known(phase), v) -
        Integral(known(Scalar(2) * phase), q) +
        DirichletBC(
          u, coefficient * p, known((Scalar(1) - Scalar(2) * coefficient) * phase))
          .on(11) +
        DirichletBC(u, known(phase)).on(12) +
        DirichletBC(p, known(Scalar(2) * phase)).on(12);
      Tuple<std::reference_wrapper<decltype(u)>, std::reference_wrapper<decltype(p)>>
        trials(std::ref(u), std::ref(p));
      Tuple<std::reference_wrapper<decltype(v)>, std::reference_wrapper<decltype(q)>>
        tests(std::ref(v), std::ref(q));
      std::array<size_t, 2> offsets{0, 3};
      boost::bimap<FormLanguage::Base::UUID, size_t> trialMap, testMap;
      trialMap.insert({u.getUUID(), 0});
      trialMap.insert({p.getUUID(), 1});
      testMap.insert({v.getUUID(), 0});
      testMap.insert({q.getUUID(), 1});
      Assembly::ProblemAssemblyInput<Body, decltype(u), decltype(v), decltype(p),
        decltype(q)>
        input(body, trials, tests, offsets, offsets, trialMap, testMap, 6, 6);
      Backend<System, ProblemType> assembler;
      System system;
      Mapping indices;
      for (size_t node = 0; node < 3; ++node)
      {
        const Index dof = space.getDOFs(0, vertices[node].second)[0];
        indices[0][node] = dof;
        indices[1][node] = 3 + dof;
      }
      verify(assembler, system, input, indices);
    }
  }

  TYPED_TEST(AffineConstraints, SequentialDense)
  {
    checkAffineConstraints<TypeParam, Math::Matrix<TypeParam>, Assembly::Sequential>();
  }
  TYPED_TEST(AffineConstraints, SequentialSparse)
  {
    checkAffineConstraints<TypeParam, Math::SparseMatrix<TypeParam>,
      Assembly::Sequential>();
  }
#ifdef RODIN_USE_OPENMP
  TYPED_TEST(AffineConstraints, OpenMPDense)
  {
    checkAffineConstraints<TypeParam, Math::Matrix<TypeParam>, Assembly::OpenMP>();
  }
  TYPED_TEST(AffineConstraints, OpenMPSparse)
  {
    checkAffineConstraints<TypeParam, Math::SparseMatrix<TypeParam>, Assembly::OpenMP>();
  }
#endif
}
